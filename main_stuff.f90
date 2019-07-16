!-------------------------------------------------------------------------------------
! Copyright (c) 2016 Daniel C. Elton
!
! This software is licensed under The MIT License (MIT)
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!-------------------------------------------------------------------------------------
module main_stuff
use consts
use system_mod
use mpi
Implicit none
!----------------------------------------------------------------------------------
!--"global" variables used by PIMD.f90 and others --------------------------------
!----------------------------------------------------------------------------------
!-old variables from TTMF code -----------------------------
double precision, dimension(:,:), allocatable :: RR, dRR
double precision, dimension(:), allocatable :: chg
double precision, dimension(3,3) :: virt
double precision, dimension(3) :: dip_mom
double precision, dimension(:,:), allocatable :: dip_momI, dip_momE
real, dimension(:,:), allocatable ::  dip_mom_all_times
real, dimension(:,:,:), allocatable :: Hvelocities
integer :: i, iw, iat,  iO, ih1, ih2, narg, ia, read_method
character(len=2) :: ch2
character(len=250) :: finp, fconfig, fsave
!-new variables:--------------------------------------------
double precision, dimension(:,:), allocatable :: VV, dRRold, dRRnew
double precision, dimension(3) :: summom, sumvel, sum_dip, avg_box, avg_box2, sum_box, sum_box2
double precision :: delt, delt2,  uk,  imassO, imassH
double precision :: temp, sum_temp, sum_press, sys_temp, avg_vel, init_energy, sum_RMSenergy
double precision :: tot_energy, sum_tot_energy,  sum_tot_energy2, sum_energy2, sum_simple_energy, sum_simple_energy2, simple_energy
double precision :: specific_heat, avg_temp, init_temp, sum_dip2, sum_simple_press, simple_sys_press, sum_pot_en_per_mol
double precision :: dielectric_constant, diel_prefac, dielectric_error, sys_press, init_time=0d0, sum_dip_mag
double precision, dimension(1000)  :: dielectric_running
integer  :: dielectric_index = 1
double precision ::  isotherm_compress
integer, dimension(:), allocatable :: seed
character(len=125) :: dip_file
character(len=11)  :: bead_thermostat_type
integer :: checkpoint_freq, hist_out_freq
integer :: run_timesteps, t, t_freq, tp_freq, td_freq, ti_freq, m, clock, eq_timesteps, ttt, tr
logical :: dip_out, coord_out, TD_out, momenta_out, TP_out, Edip_out
logical :: BAROSTAT, PEQUIL, BOXSIZEOUT, THERMOSTAT, GENVEL, INPCONFIGURATION, IMAGEDIPOLESOUT
logical :: DIELECTRICOUT, WRITECHECKPOINTS, images_out
logical :: CALCIRSPECTRA, CALCDIFFUSION, CALCDOS, RESTART, ENERGYOUT, HISTOUT, SIESTA_MON_CALC
logical ::  BEADTHERMOSTAT, CENTROIDTHERMOSTAT, CALC_RADIUS_GYRATION, CHARGESOUT, EXISTS, CALCGEOMETRY

!logical i/o units
integer :: lun, lunXYZ, lunBOXSIZEOUT, lunIMAGES, lunDIELECTRIC, lunCHARGES, lunHISTOUT
integer :: lundip_out, luncoord_out, lunCHARGESOUT, lunmomenta_out,  lunTP_out
integer :: lunEdip_out, lunIMAGEDIPOLESOUT, lunTD_out, lunimages_out, lunENERGYOUT

!N-H variables
double precision, save :: tau, tau_centroid, s, sbead
integer, save :: global_chain_length, bead_chain_length

!Nose-Hoover chain velocities for all beads are stored here
double precision, dimension(:,:,:,:), allocatable :: vxi_beads

!Nose-Hoover global chain velocities are stored here
double precision, dimension(:), allocatable     :: vxi_global

!-Berendsen thermostat variables
double precision :: tau_P, ref_P, press, CompFac, scale_factor

!- Variables for the paralleziation / PIMD
double precision, dimension(:,:,:), allocatable :: RRt, PPt, dip_momIt, dip_momEt, dRRt
double precision, dimension(:,:), allocatable :: RRc, PPc
double precision, dimension(:), allocatable :: Upott, Virialt, virialct
double precision ::  Upot,  virial, virialc, omegan, kTN, iNbeads, setNMfreq
double precision :: radiusH, radiusO, sum_radiusO, sum_radiusH
integer :: Nnodes, pid, j, k, ierr=0, Nbatches, counti, bat
integer :: status2(MPI_STATUS_SIZE)

! time variables
character(8)  :: date
character(10) :: time
double precision :: seconds

!variables for multiple timestep / contraction
integer :: intra_timesteps, num_SIESTA_nodes
double precision :: deltfast, delt2fast
character(len=200) :: siesta_name, PIMD_type, mon_siesta_name

 contains

!----------------------------------------------------------------------------------
!---------- Periodic Boundary conditions for the beads ---------------------------
!----------------------------------------------------------------------------------
!The box is assume to span [-L/2, L/2]
!The PBCs follow the bead centroid - if the bead centroid crosses the edge of the box,
!then all the beads move with it. This means that at any time, some beads may lie
!outside the box. The potential(RR, ...) subroutine must be able to handle situations
!where beads are outside the box
subroutine PBCs(RRt, RRc)
 Implicit none
 double precision, dimension(3, Natoms,Nbeads), intent(inout) :: RRt
 double precision, dimension(3, Natoms), intent(inout) :: RRc
 double precision, dimension(3, Natoms) :: shifts
 integer :: i

	!store the shifts here
	shifts = box(1)*anint(RRc*boxi(1)) !assume square box for now

 	!Correct the centroids first
	RRc = RRc - shifts

	!move the beads
	do i = 1,Nbeads
		RRt(:,:,i) = RRt(:,:,i) - shifts
	enddo

end subroutine PBCs



!----------------------------------------------------------------------------------!
!-------------- calling bead thermostat -------------------------------------------
!----------------------------------------------------------------------------------!
subroutine bead_thermostat
 use NormalModes
 use Langevin
 Implicit None
 double precision, dimension(3,Nbeads) :: PPtr
 double precision :: uk_bead, tau, imass
 Integer :: i, j, k
!- Calling Langevin thermostat -----------------------------------------------------
 if (bead_thermostat_type .eq. 'Langevin') then
	call Langevin_NM(PPt)
 endif

!- Nose-Hoover coupling in normal mode space ---------------------------------------
 if (bead_thermostat_type .eq. 'Nose-Hoover') then
   do i = 1,Natoms
	if (mod(i+2,3) .eq. 0) then
		imass = imassO
	else
		imass = imassH
	endif
	PPtr = real2NM(PPt(:,i,:),Nbeads)
	if (CENTROIDTHERMOSTAT) then
		tau = tau_centroid
		do k = 1, 3
			uk_bead =  .5d0*imass*PPtr(k,1)**2
			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,1,k), tau, delt2, 1, Nbeads*temp)
			PPtr(k,1) = PPtr(k,1)*sbead
		enddo
	endif
	do j = 2, Nbeads
		tau = 1d0/omegalist(j)
		do k = 1, 3
			uk_bead =  .5d0*imass*PPtr(k,j)**2/MassScaleFactor(j)
			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,j,k), tau, delt2, 1, Nbeads*temp)
			PPtr(k,j) = PPtr(k,j)*sbead
		enddo
	enddo
	PPt(:,i,:) = NM2real(PPtr,Nbeads)
   enddo
 endif


!----------------- Old Nose-Hoover scheme coupling in real space-------------------
! do i = 1, Natoms
!	do j = 1, Nbeads
!		do k = 1, 3
!			if (mod(i+2,3) .eq. 0) then
!				uk_bead = ( 1d0/ (2d0*massO)  )*PPt(k,i,j)**2
!			else
!				uk_bead = ( 1d0/ (2d0*massH)  )*PPt(k,i,j)**2
!			endif
!			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,j,k), tau_centroid, delt2, 1, Nbeads*temp)
!			PPt(k,i,j) = PPt(k,i,j)*sbead
!		enddo
!	enddo
!enddo
end subroutine bead_thermostat



!----------------------------------------------------------------------------------!
!-------------- calculate average radius of gyration (for testing/debugging) ------
!----------------------------------------------------------------------------------!
subroutine calc_radius_of_gyration(RRt, RRc)
 Implicit None
 double precision, dimension(3,Natoms,Nbeads),intent(in)  :: RRt
 double precision, dimension(3,Natoms),intent(in)         :: RRc
 integer    	    :: i, j, iH1, iH2, iO

 radiusH = 0d0
 radiusO = 0d0

 do i = 1, Nwaters
	iO = 3*i-2;  iH1 = 3*i-1;  iH2=3*i
	do j = 1, Nbeads
	radiusO = radiusO + sqrt( sum( (RRt(:,iO ,j) - RRc(:,iO ))**2, 1)  )
	radiusH = radiusH + sqrt( sum( (RRt(:,iH1,j) - RRc(:,iH1))**2, 1)  )
	radiusH = radiusH + sqrt( sum( (RRt(:,iH2,j) - RRc(:,iH2))**2, 1)  )
	enddo
 enddo

 radiusO = radiusO/dble(Nbeads*Nwaters)
 radiusH = radiusH/dble(2*Nbeads*Nwaters)

 sum_radiusO = sum_radiusO + radiusO
 sum_radiusH = sum_radiusH + radiusH

end subroutine calc_radius_of_gyration



!----------------------------------------------------------------------------------!
!--------------  Berendson Pressure Coupling (J. Chem. Phys. 81, 3684, 1984)-------
!----------------------------------------------------------------------------------!
subroutine Pcouple
	scale_factor = ( 1 + CompFac*(simple_sys_press - press)  )**.333333 !sum_press/tr

	RRt = RRt*scale_factor
	RRc = RRc*scale_factor
	box = box*scale_factor
	volume = box(1)*box(2)*box(3)
	boxi = 1d0/box
end subroutine Pcouple






!----------------------------------------------------------------------------------!
!-------------- Initialize bead positions -----------------------------------------
!----------------------------------------------------------------------------------!
subroutine initialize_beads
use NormalModes
double precision, dimension(:,:), allocatable :: RRtemp
double precision :: avgrO, avgrH
allocate(RRtemp(3,Nbeads))

!predict average radius of the ring polymer
if (Nbeads .gt. 1) then
	avgrO  = 2*PI*hbar/(PI*Sqrt(24*massO*KB_amuA2ps2perK*temp))
	avgrH  = 2*PI*hbar/(PI*Sqrt(24*massH*KB_amuA2ps2perK*temp))
	write(lunTP_out,'(a,f10.5,a4)') "Predicted average radius of (N_beads -> infinity) Oxygen ring polymer is   ", avgrO, " Ang"
	write(lunTP_out,'(a,f10.5,a4)') "Predicted average radius of (N_beads -> infinity) Hydrogen ring polymer is ", avgrH, " Ang"
endif

do i=1, Nwaters
		Call gen_rand_ring(RRtemp,massO,temp,Nbeads)
		do j = 1,3
			RRt(j,3*i-2,:)  =  RRc(j,3*i-2) + RRtemp(j,:)
		enddo
		Call gen_rand_ring(RRtemp,massH,temp,Nbeads)
		do j = 1,3
			RRt(j,3*i-1,:)  =  RRc(j,3*i-1) + RRtemp(j,:)
		enddo
		Call gen_rand_ring(RRtemp,massH,temp,Nbeads)
		do j = 1,3
			RRt(j,3*i-0,:)  =  RRc(j,3*i-0) + RRtemp(j,:)
		enddo
enddo

RRc = sum(RRt,3)/Nbeads !adjust RRc after bead generation

call calc_radius_of_gyration(RRt,RRc)

deallocate(RRtemp)

end subroutine initialize_beads


!---------------------------------------------------------------------------------
!---------------- Generate initial bead *momentum* ------------------------------
!---------------------------------------------------------------------------------
subroutine initialize_velocities
use NormalModes
use math
summom = 0
if (INPCONFIGURATION .and. .not. GENVEL) then
 write(lunTP_out,*) "NOTE: using velocities from configuration file"
endif
if ( (INPCONFIGURATION) .and. (GENVEL) ) then
	write(lunTP_out,*) 'WARNING: You selected to input a configuration file and generate velocites.'
	write(lunTP_out,*) '		the velocities in the configuration file will be overwritten.'
endif

if (GENVEL) then

!The program will generate Maxwell-Boltzmann velocities . With a small number of molecules,
!the initial temperature will never be exactly what is inputted
!(variances of +/- 10 K for 128 molecules) For this reason, the program regenerates
!the velocities until the temperature is within 1 K of what was specified.
   do
	summom = 0

	do i=1, Nwaters
		!generate the bead momenta from the canonical distribution for a free ring polymer
		Call gen_rand_ring_momenta(PPt(:,3*i-2,:),massO,temp*Nbeads,Nbeads)
		Call gen_rand_ring_momenta(PPt(:,3*i-1,:),massH,temp*Nbeads,Nbeads)
		Call gen_rand_ring_momenta(PPt(:,3*i-0,:),massH,temp*Nbeads,Nbeads)
		!sum up the bead momenta
		do j = 1, Nbeads
			summom = summom + PPt(:,3*i-2,j) + PPt(:,3*i-1,j) + PPt(:,3*i-0,j)
		enddo
	enddo

	!remove center of momentum from system
	do i = 1, Natoms
		do j = 1, Nbeads
			PPt(:,i,j) = PPt(:,i,j) - summom/(Natoms*Nbeads)
	 	enddo
	enddo

	!calculate centroid velocities
	PPc = sum(PPt, 3)/Nbeads !centroid momenta

	!update kinetic energy
	call calc_uk_centroid

	sys_temp = TEMPFACTOR*uk/(Natoms)

	!write(*,*) sys_temp

	if ( abs(sys_temp - temp) .lt. 1 ) then
		exit
	endif
 enddo

else if ( (.not. GENVEL) .and. (.not. INPCONFIGURATION)) then
  PPt = 0
  PPc = 0
  write(lunTP_out,*) "NOTE: Initial velocities set to zero."
endif

end subroutine initialize_velocities


!---------------------------------------------------------------------------------
!------------------------------ calc total kinetic energy -----------------------
!---------------------------------------------------------------------------------
subroutine calc_uk
 use NormalModes
 double precision, dimension(3,Nbeads) :: PPtr


if ((setNMfreq .eq. 0) .or. (Nbeads .eq. 1) ) then

 uk = 0
 do j = 1, Nbeads
	do i = 1,Nwaters
		uk = uk + imassO*sum( PPt(:,3*i-2,j)**2 )
		uk = uk + imassH*sum( PPt(:,3*i-1,j)**2 )
		uk = uk + imassH*sum( PPt(:,3*i-0,j)**2 )
	enddo
 enddo
 uk = .5d0*uk


else
!if using CMD then transform into normal mode space to calculate kinetic energy..
 uk = 0
 do i = 1,Nwaters
	PPtr =  real2NM(PPt(:,3*i-2,:),Nbeads)
	do j = 1, Nbeads
		uk = uk + imassO*sum( PPtr(:,j)**2 )/MassScaleFactor(j)
	enddo
	PPtr =  real2NM(PPt(:,3*i-1,:),Nbeads)
	do j = 1, Nbeads
		uk = uk + imassH*sum( PPtr(:,j)**2 )/MassScaleFactor(j)
	enddo
	PPtr =  real2NM(PPt(:,3*i-0,:),Nbeads)
	do j = 1, Nbeads
		uk = uk + imassH*sum( PPtr(:,j)**2 )/MassScaleFactor(j)
	enddo
  enddo
  uk = .5d0*uk

endif

end subroutine calc_uk

!---------------------------------------------------------------------------------
!------------------------------ calc bead kinetic energy -----------------------
!---------------------------------------------------------------------------------
subroutine calc_uk_centroid
 uk = 0
 do i = 1,Nwaters
 		uk = uk + imassO*sum( PPc(:,3*i-2)**2 )
 		uk = uk + imassH*sum( PPc(:,3*i-1)**2 )
 		uk = uk + imassH*sum( PPc(:,3*i-0)**2 )
 enddo
 uk = .5*uk

end subroutine calc_uk_centroid



end module main_stuff
