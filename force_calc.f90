!-----------------------------------------------------------------------------------
! Copyright (c) 2015-2016 Daniel C. Elton
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
module force_calc
 use dans_timer
 contains

!--------------------------------------------------------------------------------------
!---------------------------Full PIMD MPI force calculation --------------------------
!--------------------------------------------------------------------------------------
subroutine full_bead_forces
 use main_stuff
 use math, only:str
 Implicit None

 do bat = 0, Nbatches - 1 !batch index

    !we want to send a (3 x Natoms) array to each processor
    if (pid .eq. 0) then

 	do i = 1, Nnodes - 1 !node index
		k = bat*Nnodes + i  !bead index
		!masternode send image coordinates
		Call MPI_Send(RRt(:,:,k), counti, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, ierr)
		!masternode send centroid coordinates to calculate centroid virial
		Call MPI_Send(RRc, 3*Natoms, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, ierr)
	enddo
	!masternode force & virial calculation
	k = bat*Nnodes + Nnodes

	call potential(RRt(:,:,k), RRc, Upott(k), dRRt(:,:,k), virt, virialct(k), &
			dip_momIt(:,:,k), dip_momEt(:,:,k), chg, t, BAROSTAT, trim(sys_label)//trim(str(pid)))

	virialt(k) = virt(1,1) + virt(2,2) + virt(3,3)

	!recieve stuff from nodes
#ifdef parallel
	do i = 1, Nnodes - 1
		k = bat*Nnodes + i
		!masternode receive derivatives
		call MPI_Recv(dRRt(:,:,k), counti, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		!masternode recieve energies
		call MPI_Recv(Upott(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		!masternode receive virials
		call MPI_Recv(virialt(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		!masternode receive centroid virials
		call MPI_Recv(virialct(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		!masternode recieve dipole moments
		if (dip_out .or. TD_out) call MPI_Recv(dip_momIt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		if (Edip_out) call MPI_Recv(dip_momEt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
	enddo
#endif
    else
#ifdef parallel
    !slavenode recieve coords from master node
		call MPI_Recv(RR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
		!slavenode recieve centroid coords from master node
		call MPI_Recv(RRc, 3*Natoms, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
		!slavenode force & virial calculation
		call potential(RR, RRc, Upot, dRR, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT,trim(sys_label)//trim(str(pid)))
		virial = virt(1,1) + virt(2,2) + virt(3,3)
		!slavenode send back derivative
		call MPI_Send(dRR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		!slavenode send back energy
		call MPI_Send(Upot, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		!slavenode send back virial
		call MPI_Send(virial, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		!slavenode send back centroid virial
		call MPI_Send(virialc, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		!slavenode send back dipole moments
		if (dip_out .or. TD_out) call MPI_Send(dip_momI, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		if (Edip_out) call MPI_Send(dip_momE, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
#endif
    endif
    Call MPI_Barrier(MPI_COMM_WORLD, ierr)

 enddo! j = 0, Nbatches - 1 !batch index

if (pid .eq. 0) then
 Upot    = sum(Upott)    	!potential energy for the ENTIRE system
 virial  = sum(virialt)  	!virial for the ENTIRE system (kcal/mol)
 virialc = sum(virialct)/Nbeads !centroid virial for the ENTIRE system (kcal/mol)

 !computation of dipole moments for the SIESTA case
 if (pot_model .eq. 6) then
	call calc_monomer_dip_moments(dip_momIt, RRt)
 endif

endif

end subroutine full_bead_forces


!---------------------------------------------------------------------
!-----------------Call correct potential ----------------------------
!---------------------------------------------------------------------
subroutine potential(RR, RRc, Upot, dRR, virt, virialc, dip_momI, Edip_mom, chg, t, BAROSTAT, sys_lab)
 use consts
 use system_mod
 use pot_mod
 use fsiesta
 use dans_timer
 implicit none
 double precision, dimension(3, Natoms), intent(in) :: RR, RRc
 double precision, intent(out) :: Upot, virialc
 double precision, dimension(3, Natoms), intent(out) :: dRR
 double precision, dimension(3, 3), intent(out) :: virt
 double precision, dimension(Natoms), intent(out)  :: chg
 double precision, dimension(3, NWaters), intent(out)  ::  dip_momI, Edip_mom
 character(len=*), intent(in) :: sys_lab
 double precision, dimension(3) :: dip_mom
 integer, intent(in) :: t
 logical, intent(in) :: BAROSTAT

 siesta_box = 0.0
 siesta_box(1,1) = box(1)
 siesta_box(2,2) = box(2)
 siesta_box(3,3) = box(3)

 if (BAROSTAT) then
	!All the stuff that depends on volume needs to be rescaled if using the barostat
	p4V = FOURPI/volume

	!scale VdW long range correction due to box size change and add correction
	!multiplying by a correction factor is slightly more efficient than recalculating the entire Uvdw_lrc term each timestep
	Uvdw_lrc = Uvdw_lrc0*(volume_init/volume)

	!If the volume is changing than the Ewald k-vectors have to be reset every timestep
	call ewald_set(.false.)
 endif!(BAROSTAT)

if (pot_model==2 .or. pot_model==3) then

	call pot_ttm(RR, RRc, Upot, dRR, virt, virialc, dip_momI, Edip_mom, chg, t)

 else if (pot_model==4 .or. pot_model==5) then

	call pot_spc(RR, Upot, dRR, virt, dip_momI, chg)

 else if (pot_model==6) then

	call start_timer("SIESTA")
    call siesta_forces( trim(sys_lab), Natoms, RR, cell=siesta_box, energy=Upot, fa=dRR)
    call stop_timer("SIESTA")

    !write(*,*) "energy (eV):    ", Upot
    !write(*,*) "foces (ev/Ang): ", -1.0*dRR

    Upot = Upot*EVTOKCALPERMOLE
    !negative sign is to be consistent with rest of code: dRR is potenial derivative, not force
    dRR = -1d0*dRR*EVTOKCALPERMOLE

 endif

end subroutine potential


!---------------------------------------------------------------------
!-----------------SIESTA monomer force calculation ------------------
!---------------------------------------------------------------------
subroutine siesta_monomer(r1, dr1, e1)
 use fsiesta
 use system_mod
 use dans_timer
 use consts
 implicit none
 double precision, dimension(3, 3), intent(in) :: r1
 double precision, dimension(3, 3), intent(out) :: dr1
 double precision, intent(out) :: e1

 siesta_box = 0.0
 siesta_box(1,1) = box(1)
 siesta_box(2,2) = box(2)
 siesta_box(3,3) = box(3)

 call start_timer("SIESTAmonomer")
 call siesta_forces( "monomer", 3, r1, cell=siesta_box, energy=e1, fa=dr1)
 call stop_timer("SIESTAmonomer")

 e1 = e1*EVTOKCALPERMOLE
 dr1 = -1d0*dr1*EVTOKCALPERMOLE

end subroutine siesta_monomer


!---------------------------------------------------------------------
!- Calculate dipole moments using the TIP4P/2005 charges and m-site
!- for the coordinates obtained from a SIESTA calculation
!---------------------------------------------------------------------
subroutine calc_monomer_dip_moments(dip_momIt, RRt)
 use consts
 use pot_mod
 use system_mod !source of Nbeads, box, boxi
 implicit none
 double precision, dimension(3, Nwaters, Nbeads), intent(out) :: dip_momIt
 double precision, dimension(3, 3*Nwaters, Nbeads), intent(in) :: RRt
 double precision, dimension(3,3)   :: r1
 double precision, dimension(3,3,3) :: dq3
 double precision, dimension(3) :: roh1, roh2, r3, summ, q3, Rm, rh1m, rh2m
 double precision :: e1, chgH1, chgH2, tmp
 integer :: iw, j, io, iH1, iH2
 double precision, parameter :: rOM = .1546
 double precision, parameter :: qH_TIP4P = .5564
 double precision, parameter :: qO_TIP4P = -1.1128

 tmp = 0.5d0*gammaM/(1.d0-gammaM)

 do j = 1, Nbeads
	do iw = 1, Nwaters
		iO=3*iw-2; iH1 = 3*iw-1; iH2=3*iw-0

  		r1(1:3, 1:3) = RRt(1:3, (/iO, iH1, iH2/), j)

  		!get charges
		call dms_nasa(r1, q3, dq3,box,boxi)

		q3 = q3*CHARGECON
		chgH1 = q3(2) + tmp*(q3(2)+q3(3))
		chgH2 = q3(3) + tmp*(q3(2)+q3(3))

		roh1 = RRt(:, iH1, j) - RRt(:, iO, j)
		roh1 = roh1 - box*anint(roh1*boxi)!PBC
		roh2 = RRt(:, iH2, j) - RRt(:, iO, j)
  		roh2 = roh2 - box*anint(roh2*boxi)!PBC

		dip_momIt(:,iw,j) = chgH1*roh1 + chgH2*roh2

     	!get m site position if necessary
		if ((pot_model .eq. 3) .or. (pot_model .eq. 4)) then
			Rm = 0.5d0*gammaM*( roh1(:) + roh2(:) ) + RRt(:, iO, j)

			rh1m = RRt(:, iH1, j) - Rm
			rh1m = rh1m - box*anint(rh1m*boxi)!PBC
			rh2m = RRt(:, iH2, j) - Rm
			rh2m = rh2m - box*anint(rh2m*boxi)!PBC

			dip_momIt(:,iw,j) = chgH1*rh1m + chgH2*rh2m
		else
			dip_momIt(:,iw,j) = chgH1*roh1 + chgH2*roh2
		endif

		!!! TIP4P/2005f dipole moment calculation
		!summ = roh1 + roh2!find vector to M-site
		!r3=(summ/sqrt(dot_product(summ,summ)))*rOM
		!dip_momIt(:,iw,j) = (qH_TIP4P*roh1 + qH_TIP4P*roh2 + qO_TIP4P*r3)*a2m*e2coul/3.33564e-30!conv. to Debye
	enddo
 enddo

end subroutine calc_monomer_dip_moments


end module force_calc
