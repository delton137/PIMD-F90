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
! TH SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!-------------------------------------------------------------------------------------
module InputOutput
 use consts
 use main_stuff
 use lun_management
 use spectral_properties
Implicit none

 contains
!----------------------------------------------------------------------------------
!---------------- Read input file ------------------------------------------------
!----------------------------------------------------------------------------------
subroutine read_input_file
 use input_file_reader
 implicit none

 call getarg(1, finp)
 narg = command_argument_count()
 if (narg .eq. 0) then
    write(*,*) "ERROR: no input file argument."
	call abort()
 endif

 call io_assign(lun)
 open(lun, file=finp, status='old')

 call read_inputfile()

 select case (trim(PIMD_type))
	case("full", "fullPIMD")
		MONOMERPIMD = .false.
		CONTRACTION = .false.
	case("contracted", "conPIMD", "contractedPIMD")
		MONOMERPIMD = .false.
		CONTRACTION = .true.
	case("monomerPIMD", "monPIMD")
		MONOMERPIMD = .true.
		CONTRACTION = .false.
	case default
		write(*,*) "InputOutput: ERROR: invalid PIMD type. Must be 'full', 'contracted', or 'monomerPIMD'"
 end select


end subroutine read_input_file


!----------------------------------------------------------------------------------
!---------------- Initialize some variables for all nodes ------------------------
!----------------------------------------------------------------------------------
subroutine read_and_initialize_all_nodes
 integer :: NbeadsIn, slen

 !read the number of atoms, dimension of box and atomic coordinates
 call io_assign(lunXYZ)
 open(lunXYZ,file=fconfig,status='old', action="read")

 slen = len(trim(fconfig))

 !Read file suffix and determine if its an .xyz or image (.img) file
 if (fconfig(slen-3:slen) .eq. ".img") then
	INPCONFIGURATION = .true.
 else if (fconfig(slen-3:slen) .eq. ".xyz") then
	INPCONFIGURATION = .false.
 else
	write(*,*) "ERROR: could not determine input file format. File must be *.xyz or *.img"
	call abort()
 endif

 if (INPCONFIGURATION) then

	read(lunXYZ,'(2i10)') Natoms, NbeadsIn

	if (.not.(NbeadsIn .eq. Nbeads)) then
		write(*,*) "ERROR: the number of beads in configuration image not equal to that specifed in input file"
		write(*,*) " are you sure you have the right file?"
		call abort()
	endif

	read(lunXYZ,'(f12.6,2x,f12.6,3(1x,f12.6))') init_time, box

 else
	!usually the box size is in the first line of a raw .xyz
	!but it might be in the second line. Try a few different combos here.
	read(lunXYZ,*,IOSTAT=ierr) box(1:3)
	read(lunXYZ,*) Natoms
	if (ierr .ne. 0) then
		rewind(lunXYZ)
		read(lunXYZ,*) Natoms
		read(lunXYZ,*,IOSTAT=ierr) box(1:3)
		if (ierr .ne. 0) then
			rewind(lunXYZ)
			read(lunXYZ,*,IOSTAT=ierr) Natoms, box(1:3)
			if (ierr .ne. 0) then
				write(*,*) "ERROR: could not read box size from input file. The best I could come up with is box = ", box
				call abort()
			endif
		endif
	endif
 endif

!Automatic adujustment of cutoff if necessary
 if (Rc .eq. -1) then
 	write(lunTP_out,*) 'you did not specify a cuttoff, so I am settint it to Lmin/2 =', minval(box)/2.0
	Rc = minval(box)/2.0
 endif
 if (Rc1 .eq. -1) then
	Rc1 = 0.8*Rc
 endif

!These parameters are used later on
Rc2 = Rc * Rc
Nwaters = Natoms/3
volume = box(1)*box(2)*box(3)
volume_init = volume
delt = delt/1000d0 !***CONVERT fs - > ps ***
delt2 = delt/2d0
boxi = 1.d0 / box
!inverse masses
imassO = DBLE(1/massO)
imassH = DBLE(1/massH)
iNbeads = 1d0/DBLE(Nbeads)
counti = 3*Natoms
omegan = KB_amuA2ps2perK*temp*Nbeads/hbar
kTN = KB_amuA2ps2perK*temp*Nbeads
s = 1
sbead = 1

Nbatches = Nbeads/Nnodes

!--- Slave-node-only allocations ----
if (pid .ne. 0) then
	allocate(RR(3, Natoms))
	allocate(VV(3, Natoms))
	allocate(dRR(3, Natoms))
endif

!--- All-node allocations ------
allocate(dip_momI(3, Nwaters))
allocate(dip_momE(3, Nwaters))
allocate(chg(Natoms))
!allocate(tx_dip(3,4*Nwaters, 4))
allocate(RRc(3, Natoms))

end subroutine read_and_initialize_all_nodes


!----------------------------------------------------------------------------------
!-------------------------  Initialize SIESTA ------------------------------------
! note: call siesta_units( "ang", 'kcal/mol' ) doesn't work- Siesta cant do kcal/mol
!----------------------------------------------------------------------------------
subroutine init_siesta
 use fsiesta
 use math, only: str
 Implicit None
 character(len=400) :: sys_command
 character(len=100) :: sys_label_i
 character(len=300) :: pipe_name
 integer :: nodes_per_process

 if (pid .eq. 0) then

 if ((SIESTA_MON_CALC) .and. (PIMD_type .eq. "full")) then
	write(*,*) "ERROR : SIESTA_MON_CALC does not work with full PIMD"
	call abort()
 endif

 if (SIESTA_MON_CALC) then
	if (num_SIESTA_nodes .eq. 1) then
		call siesta_launch( trim(mon_siesta_name), "monomer") !launch serial SIESTA process
	elseif (num_SIESTA_nodes .gt. 1) then
		call siesta_launch( trim(mon_siesta_name), "monomer", nnodes=num_SIESTA_nodes ) !launch parallel SIESTA process
	else
		write(*,*) "InputOuput: ERROR: invalid number of SIESTA nodes!!"
		call abort()
	endif
 endif

 pot_model = 3
 guess_initdip = .false. !disable this, because it could cause problems
 call init_pot !initialize ttm for dipoles calculations
 pot_model = 6

 call system("export GFORTRAN_UNBUFFERED_ALL=y")
 endif !pid .eq. 0

 if ( PIMD_type .eq. "full") then
	call sleep(pid) !stagger the system calls from different MPI processes a bit - important to do this!

	!setup Nnodes SIESTA proccesses sharing num_SIESTA_nodes
	nodes_per_process = floor(real(num_SIESTA_nodes)/real(Nnodes))

	sys_label_i = trim(sys_label)//trim(str(pid))

	!make copies of .fdf
	sys_command = "cp "//trim(sys_label)//".fdf "//trim(sys_label_i)//".fdf"

	call system(trim(sys_command))

	sys_command = "sed -i -- 's/"//trim(sys_label)//"/"//trim(sys_label_i)//"/g' "//trim(sys_label_i)//".fdf"

	call system(trim(sys_command))

	call siesta_launch(trim(siesta_name), trim(sys_label_i), nnodes=nodes_per_process ) !launch parallel SIESTA process

 elseif (Nnodes .eq. 1) then
	if (num_SIESTA_nodes .eq. 1) then
		call siesta_launch( trim(siesta_name), trim(sys_label)) !launch serial SIESTA process
	elseif (num_SIESTA_nodes .gt. 1) then
		call siesta_launch(trim(siesta_name), trim(sys_label), nnodes=num_SIESTA_nodes ) !launch parallel SIESTA process
	else
		write(*,*) "InputOuput: ERROR: invalid number of SIESTA nodes!!"
		call abort()
	endif
 else
	write(*,*) "You specified running with ", PIMD_type, " and Nnodes = ", Nnodes
	write(*,*) "This configuration is not supported. With contractiom PIMD must be on a single node. "
	call abort()
 endif

end subroutine init_siesta

!----------------------------------------------------------------------------------
!---------- Error handling  / master node allocations ----------------------------
!----------------------------------------------------------------------------------
subroutine master_node_init
	use Langevin
	use NormalModes
    Implicit None

if (Nbeads .lt. 1) then
	write(*,*) "ERROR : invalid number of beads!! "
endif

if ( .not.( (bead_thermostat_type .eq. 'Langevin') .or. (bead_thermostat_type .eq. 'Nose-Hoover') &
	.or. (bead_thermostat_type .eq. 'none') ) ) then
		write(*,*) "ERROR: Invalid bead thermostat selection. Possible options: "
		write(*,*) "		'Nose-Hoover', 'Langevin', or 'none'  "
	call abort()
endif

if ((CENTROIDTHERMOSTAT .or. BEADTHERMOSTAT) .and. (Nbeads .eq. 1)) then
	write(*,*) "InputOutput: WARNING: Bead/centroid thermostating does not make much sense with 1 bead."
	write(*,*) "InputOutput: WARNING: The dynamics will be unphysical since every atomic DOF will be thermostated. "
	write(*,*) "InputOutput: WARNING: --- The bead thermostat is being turned off !! ----"
	if (CENTROIDTHERMOSTAT) CENTROIDTHERMOSTAT = .false.
	if (BEADTHERMOSTAT) BEADTHERMOSTAT = .false.
endif

if (BEADTHERMOSTAT .and. .not. (THERMOSTAT)) then
	write(*,*) "WARNING: running bead thermostating without a global thermostat is not recommended."
	write(*,*) "You may observe abnormally large temperature fluctuations in the system."
endif

if (CENTROIDTHERMOSTAT.and. .not. (BEADTHERMOSTAT)) then
	write(*,*) "WARNING: You are thermostating the centroid but not thermostating the other modes."
	write(*,*) "There is not really any good reason for doing this. Consider a different scheme."
	call abort()
endif

if ( Rc .gt. minval(box)/2 ) then
	write(*,*) 'ERROR: cutoff radius is greater than half the smallest box dimension (', minval(box), ')'
	write(*,*) 'suggest changing to', minval(box)/2.0
	call abort()
endif
if (rc1 .lt. 0) then
 	write(*,*) "ERROR: start of shifted cutoff cannot be less than zero!!"
	call abort()
endif
if (Rc .lt. 0) then
 	write(*,*) "ERROR: Coloumb cutoff cannot be less than zero!!"
	call abort()
endif
if (rc1 .gt. Rc) then
 	write(*,*) "ERROR: start of shifted cutoff cannot be greater than Coloumb cuttoff!!"
	call abort()
endif
if ( (massH .lt. 0) .or. (massO .lt. 0)) then
 	write(*,*) "Invalid mass!!"
	call abort()
endif

if ( Nnodes .gt. Nbeads) then
	write(*,*) "ERROR : The number of processors is greater than the number of beads! Assuming this is an error! "
	call abort()
endif

if ( (pot_model .gt. 6) .or. (pot_model .lt. 1) ) then
	write(*,*) "ERROR: Invalid potential model !"
	call abort()
endif

if (.not. (CONTRACTION) ) then
	if (.not. (mod(Nbeads,Nnodes) .eq. 0)) then
		write(*,*) "ERROR: the number of beads must be a multiple of the number of nodes."
     	   write(*,'(a,i4,a,i4,a)') "To run on ", Nnodes, " nodes I suggest using ", Nbeads - mod(Nbeads,Nnodes), " beads"
	call abort()
	endif
else
	if (Nnodes .gt. 1) then
		write(*,*) "ERROR: When running with ring polymer contraction, a maximum of 1 nodes can be used."
		call abort()
	endif
endif

if ((CONTRACTION) .and. (Nbeads .eq. 1)) then
	write(*,*) "NOTE: You have specified contraction with one bead. "
endif

if (INPCONFIGURATION .and. GENVEL) then
	write(*,*) "NOTE: Generating velocities not necessary when using inputted configuration image"
	write(*,*) " using velocities in image file"
	GENVEL = .false.
endif
if (INPCONFIGURATION .and. (.not.(eq_timesteps .eq. 0))) then
	write(*,*) "NOTE: when inputing image, equilibration not necessary. Setting equilibration timesteps to 0"
	eq_timesteps = 0
endif


	CompFac = ((4.477d-5)*delt)/(3*tau_P) !Barostat var. (contains compressibility of H2O)
	sum_temp = 0
	sum_press = 0
	sum_pot_en_per_mol = 0
	sum_tot_energy = 0
	sum_simple_energy = 0
	sum_simple_energy2 = 0
	sum_simple_press = 0
	sum_tot_energy2 = 0
	sum_RMSenergy = 0
	sum_box = 0
	sum_box2 = 0
	tr = 0
	ttt = 0
	diel_prefac =  debyeSI**2 / (3 * kbSI * vac_permSI * a2m**3)

	!initialize random number generator
 	CALL RANDOM_SEED(size = m) !get size of seed for the system
 	ALLOCATE(seed(m))
	call system_clock(count=clock)
	seed = clock + 357 * (/ (i - 1, i = 1, m) /)
 	call random_seed(put = seed)  !put in the seed

	!write(lunTP_out,'(a,i4,a,i4,a)') "Running with ", Nbeads, " beads on ", Nnodes, " nodes"

	!Master node allocations
	!only the master node (pid = 0) stores a fully copy of the
	! coords / vel for all beads and the centroid
	allocate(RRt(3, Natoms,Nbeads))
	allocate(PPt(3, Natoms,Nbeads))
	allocate(dRRt(3, Natoms,Nbeads))
	allocate(dip_momIt(3, Nwaters,Nbeads))
	allocate(dip_momEt(3, Nwaters,Nbeads))
	allocate(Upott(Nbeads))
	allocate(Virialt(Nbeads))
	allocate(virialct(Nbeads))
	allocate(PPc(3, Natoms))
    if (CALCIRSPECTRA) allocate(dip_mom_all_times(3, run_timesteps))
    if (CALCDOS)  	   allocate(Hvelocities(3, run_timesteps, 2*Nwaters))

	dRRt = 0

 	if (CONTRACTION) deltfast = delt/intra_timesteps
 	if (CONTRACTION) delt2fast = deltfast/2d0

	if (CONTRACTION) then
		call InitNormalModes(Nbeads, omegan, deltfast, setNMfreq, lunTP_out)
	else
		call InitNormalModes(Nbeads, omegan, delt, setNMfreq, lunTP_out)
	endif

    if (.not.(bead_chain_length==2)) then
        write(*,*) "NOTE:  bead Nose Hoover chain lenght must =2. Setting to 2"
        bead_chain_length = 2
    endif
    if (.not.(global_chain_length==2)) then
        write(*,*) "NOTE: global  Nose Hoover chain lenght must =2. Setting to 2"
        global_chain_length = 2
    endif

	if (THERMOSTAT)  then
		allocate(vxi_global(global_chain_length))
		vxi_global = 1 !set chain velocities to zero initially
	endif
	if (BEADTHERMOSTAT)  then
		allocate(vxi_beads(bead_chain_length,natoms,Nbeads,3))
		vxi_beads = 0 !set chain velocities to zero initially
	endif
	if (bead_thermostat_type .eq. 'Langevin') call Init_Langevin_NM(delt2, CENTROIDTHERMOSTAT, tau_centroid, Nbeads, Nbeads*temp)

	!write Temp/Press file header
     call date_and_time(DATE=date,TIME=time)
     write(lunTP_out,'(a)') " started on  "//date(1:4)//"-"//date(5:6)//"-"//date(7:8)//&
		" at "//time(1:2)//":"//time(3:4)//":"//time(5:10)

	call print_basic_run_info

end subroutine master_node_init




!----------------------------------------------------------------------------------
!---------------- read in and intialize RRc from ------------------------
!----------------------------------------------------------------------------------
Subroutine read_coords_and_init
 double precision, dimension(3) :: avgRR

 avgRR = 0
if (INPCONFIGURATION) then
	call load_checkpoint(lunXYZ)
else
	if (read_method .eq. 0) then
		do i=1, Nwaters
			iO = 3*i-2
			read(lunXYZ,*)ch2, RRc(1:3, iO)
		enddo
		do i=1, Nwaters
			ih1 = 3*i-1; ih2=3*i
			read(lunXYZ,*)ch2, RRc(1:3, ih1)
			read(lunXYZ,*)ch2, RRc(1:3, ih2)
		enddo
	else if (read_method==1) then
		do i=1, Natoms
			read(lunXYZ,*)ch2, RRc(1:3, i)
			avgRR = avgRR    + RRc(:, i)
		enddo
	else
		write(*,*) "ERROR: Invalid read method!!"
		call abort()
	endif

	avgRR = avgRR/Natoms
	do i=1, Natoms
		RRc(:,i) = RRc(:,i) - avgRR
	enddo

	call initialize_beads      ! Initialize RRt
	call initialize_velocities ! Initialize PPt
endif
 call io_close(lunXYZ)


EndSubroutine read_coords_and_init


!-----------------------------------------------------------------------------------------
!--------------Open write out files -----------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine open_files
 Implicit none
 logical :: EXISTS

 if (TP_out) then
	call io_open(lunTP_out,'out_'//TRIM(fsave)//'_TempPress.dat',APPEND=RESTART)
 else
	lunTP_out = io_get_stdout()
 endif

 if (coord_out)      call io_open(luncoord_out,'out_'//TRIM(fsave)//'_coord.xyz',APPEND=RESTART)
 if (momenta_out)    call io_open(lunmomenta_out,'out_'//TRIM(fsave)//'_momenta.dat',APPEND=RESTART)
 if (images_out)     call io_open(lunimages_out,'out_'//TRIM(fsave)//'_images_coord.xyz',APPEND=RESTART)
 if (dip_out)        call io_open(lundip_out,'out_'//TRIM(fsave)//'_dip.dat',APPEND=RESTART)
 if (TD_out)         call io_open(lunTD_out,'out_'//TRIM(fsave)//'_tot_dip.dat',APPEND=RESTART)
 if  (Edip_out)      call io_open(lunEdip_out,'out_'//TRIM(fsave)//'_Edip.dat',APPEND=RESTART)
 if (BOXSIZEOUT)     call io_open(lunBOXSIZEOUT,'out_'//TRIM(fsave)//'_box.dat',APPEND=RESTART)
 if (CHARGESOUT)     call io_open(lunCHARGESOUT,'out_'//TRIM(fsave)//'_chgs.dat',APPEND=RESTART)
 if(IMAGEDIPOLESOUT) call io_open(lunIMAGEDIPOLESOUT,'out_'//TRIM(fsave)//'_images_dip.dat',APPEND=RESTART)
 if (ENERGYOUT)      call io_open(lunENERGYOUT,'out_'//TRIM(fsave)//'_energy.dat',APPEND=RESTART)
 if (HISTOUT)        call io_open(lunHISTOUT,'out_'//TRIM(fsave)//'_histogram.dat',APPEND=RESTART)

end subroutine open_files


!----------------------------------------------------------------------------------
!-------------- Calculate thermodynamic info and write out to file(s) ------------
!----------------------------------------------------------------------------------
subroutine write_out
 use geometry_calculator
 use diffusion_calculator
 use dans_timer
 use estimators
 implicit none

 !for accuracy, pressure is computed at every timestep, since P fluctuations are large
 !total energy and temperature is also calculated every timestep
 !in some instances, slight speedups were sacrificed for code readability

 !reset averaging after equilbration ends
 if (t .eq. eq_timesteps + 1) then
 	write(lunTP_out,*) "#----end of equilibration---restarting averages----------------------------------"
	!store average temp during equil and final energy after equil
	init_energy = tot_energy
	init_temp = sum_temp/tr

	!start new averaging
	tr     = 0
	ttt    = 0
	sum_temp = 0
	sum_press = 0
	sum_dip_mag = 0
	sum_dip = 0
	sum_dip2 = 0
	sum_tot_energy = 0
	sum_simple_energy = 0
	sum_simple_energy2 = 0
	sum_simple_press = 0
	sum_tot_energy2 = 0
	sum_RMSenergy = 0
	sum_radiusH = 0
	sum_radiusO = 0
	sum_pot_en_per_mol = 0
	if (CALCGEOMETRY) call write_out_geometry(lunTP_out,Nbeads)
	if (CALCGEOMETRY) call reset_geometry
 endif

 !---- write columns header ------------------------------------------------------------------
 if ((t .eq. 1) .or. (t .eq. eq_timesteps + 1)) then
	write(lunTP_out,'(a)',advance='no') " |time (ps)"
	write(lunTP_out,'(a)',advance='no') "|temp (K)"
	!write(lunTP_out,'(a)',advance='no') "|press.(bar)|"
	write(lunTP_out,'(a)',advance='no') "|AvgTemp|"
	!write(lunTP_out,'(a)',advance='no') "|avgP(vir)"
	!write(lunTP_out,'(a)',advance='no') "|P(simp)"
	write(lunTP_out,'(a)',advance='no') "|avgP(simp)"
	!write(lunTP_out,'(a)',advance='no') "|Tot E(vir)"
	!write(lunTP_out,'(a)',advance='no') "|avgTotE(vir)"
	write(lunTP_out,'(a)',advance='no') "|Tot E (simple)"
	write(lunTP_out,'(a)',advance='no') "|avgTotE(simp)"
	write(lunTP_out,'(a)',advance='no') "|avgPotE/mol"
!	if (CALC_RADIUS_GYRATION)    write(lunTP_out,'(a,a)',advance='no') "|r_O", "|r_H "
	if (CALC_RADIUS_GYRATION)    write(lunTP_out,'(a,a)',advance='no') "|avg r_O", "|avg r_H "
	if (DIELECTRICOUT)           write(lunTP_out,'(a)',advance='no') "|eps(0)"
    write(lunTP_out,'(a)',advance="no") "|avgDipMom(D)"
	!write(lunTP_out,'(a)',advance=no") " density kg/m^3 "
	write(lunTP_out,'(a)',advance='no') "|[ energies in kcal/mol ]"
	write(lunTP_out,'(a)',advance='yes') ""
 endif

 tr = tr + 1

 sys_temp = TEMPFACTOR*uk/(Natoms*Nbeads*Nbeads)

! call calc_uk_centroid
! write(lunTP_out,*) "centroid temp =", TEMPFACTOR*uk/(Natoms)
!uk = 0
!do j = 1, Nbeads
!	do i = 1,Nwaters
!		uk = uk + imassO*sum( PPt(:,3*i-2,j)**2 )
!		uk = uk + imassH*sum( PPt(:,3*i-1,j)**2 )
!		uk = uk + imassH*sum( PPt(:,3*i-0,j)**2 )
!	enddo
!enddo
!uk = .5d0*uk
! write(lunTP_out,*) "naive bead temp =", TEMPFACTOR*uk/(Natoms)
!!- pressure / total energy calculation : old classical case -
!sys_press =  PRESSCON*(1/(3*volume))*( 2*uk -	 MASSCON*( virt(1,1)+virt(2,2)+virt(3,3) )  )

 call quantum_virial_estimators(RRt, virial, virialc, tot_energy, sys_press, sys_temp, Upot)
 call simple_quantum_estimators(RRt, virial, simple_energy, simple_sys_press, sys_temp, Upot)

 sum_simple_energy = sum_simple_energy + simple_energy
 sum_simple_energy2 = sum_simple_energy2 + simple_energy**2
 sum_simple_press  = sum_simple_press + simple_sys_press
 sum_press         = sum_press + sys_press
 sum_tot_energy    = sum_tot_energy + tot_energy
 sum_temp          = sum_temp + sys_temp
 sum_tot_energy2   = sum_tot_energy2 + tot_energy**2
 sum_RMSenergy     = sum_RMSenergy + (tot_energy - sum_tot_energy/tr)**2
 sum_pot_en_per_mol = sum_pot_en_per_mol + sum(Upott)/dble(Nwaters)/dble(Nbeads)

 !!debugging view
 !write(*,*) "Upot   "  , Upot
 !write(*,*) "virial "  , virial
 !write(*,*) "virialc"  , virialc
 !write(*,*) "simple P" , simple_sys_press
 !write(*,*) "virial P" , sys_press
 !write(*,*) "simple E" , simple_energy
 !write(*,*) "virial E" , tot_energy

 !first, convert dipoles in all images into Debye
 dip_momIt = dip_momIt*DEBYE/CHARGECON
 dip_momEt = dip_momEt*DEBYE/CHARGECON

 !calculate dipole moment by averaging over all beads
 do iw=1,Nwaters
	do j = 1, 3
		dip_momI(j,iw) = sum(dip_momIt(j,iw,:))/Nbeads
	enddo
 enddo

 !caculate total dipole moment (in Debye)
 dip_mom(:) = sum(dip_momI(:,:), dim=2)
 sum_dip  = sum_dip  + dip_mom

 do iw=1,Nwaters
	sum_dip_mag = sum_dip_mag + dsqrt(dot_product(dip_momI(:,iw), dip_momI(:, iw)))
 enddo

 !update quantities for dielectric constant
 !it really isn't necessary to do this every timestep, so we do it every 10 steps
 if (DIELECTRICOUT .and. ( mod(t,10) .eq. 0 )  ) then
	sum_dip2 = sum_dip2 + sum(dip_mom**2)
	ttt = ttt + 1
 endif

 if (CALCGEOMETRY .or. HISTOUT) then
	call start_timer("calc_geometry")
	call calc_geometry(RRc, RRt)
	call stop_timer("calc_geometry")
 endif

 if (CALC_RADIUS_GYRATION) call calc_radius_of_gyration(RRt,RRc)

 !-------------------------------------------------------------
 !print out temperature, pressure, average press, energies & dielectric constant
 !-------------------------------------------------------------
  if (mod(t,tp_freq) == 0) then
    write(lunTP_out,'(1a)',advance='no')    "@"
    write(lunTP_out,'(1f10.3)',advance='no')  tr*delt + init_time
	write(lunTP_out,'(1f10.2)',advance='no')  sys_temp
	!write(lunTP_out,'(1f10.2)',advance='no') sys_press
	write(lunTP_out,'(1f10.2)',advance='no') sum_temp/tr
	!write(lunTP_out,'(1f11.2)',advance='no') sum_press/tr
    !write(lunTP_out,'(1f10.2)',advance='no') simple_sys_press
	write(lunTP_out,'(1f10.2)',advance='no')  sum_simple_press/tr
	!write(lunTP_out,'(2f10.2)',advance='no') tot_energy
	!write(lunTP_out,'(1f10.2)',advance='no') sum_tot_energy/tr
    write(lunTP_out,'(1f12.3)',advance='no')  simple_energy
	write(lunTP_out,'(1f12.3)',advance='no')  sum_simple_energy/tr
	write(lunTP_out,'(f12.3)',advance='no') sum_pot_en_per_mol/tr

	if (CALC_RADIUS_GYRATION) then
		!write(lunTP_out,'(1x,f6.4,1x,f6.4)',advance='no')  radiusO, radiusH
		write(lunTP_out,'(1x,f6.4,1x,f6.4)',advance='no')  sum_radiusO/tr, sum_radiusH/tr
	endif

	!calculate dielectric constant using current volume and average temperature of the run
	if (DIELECTRICOUT) then
		dielectric_constant = diel_prefac*(  sum_dip2/ttt - sum( (sum_dip/ttt)**2 )  )/volume/(sum_temp/tr)
		write(lunTP_out,'(1x,f6.2)',advance='no') dielectric_constant
		!if (mod(t,run_timesteps/1000) .eq. 0) then
		!	dielectric_running(dielectric_index) = dielectric_constant
		!	dielectric_index = dielectric_index + 1
		!endif
	endif

	write(lunTP_out,'(1x,f8.4)',advance='no') sum_dip_mag/tr/Nwaters !avg dipole moment

	!! feature to output the current density (for debugging the barostat)
	!write(lunTP_out,'(1x,f10.6)',advance='no') Nwaters*(massO+2*massH)*amu2grams/(box(1)*box(2)*box(3)*(a2m*100)**3)

	!advance to next line
	write(lunTP_out,'(a)') ""

	call pflush(lunTP_out) !flush I/O buffer
 endif

 !-------------------------------------------------------------
 !write out / store data during run phase
 !-------------------------------------------------------------
 if (t .gt. eq_timesteps) then

	if (CALCIRSPECTRA) dip_mom_all_times(:, tr) = real(dip_mom(:))

    if (CALCDOS) then
        !calculate hydrogen centroid momenta
		PPc = sum(PPt,3)/Nbeads
		do i = 1, Nwaters
			Hvelocities(:, tr, 2*i-0) = real(PPc(:, 3*i-1)*imassH)
			Hvelocities(:, tr, 2*i-1) = real(PPc(:, 3*i-2)*imassH)
		enddo
    endif!(CALCDOS)

    if (CALCDIFFUSION) then
        call start_timer("calc_diffusion")
		call calc_diff_RMSD(RRc, run_timesteps)
		call stop_timer("calc_diffusion")
    endif

    if (mod(t,t_freq)  == 0 ) then
		!coordinate output
		if (coord_out) then
			call save_XYZ(luncoord_out, RRc, Upot, read_method, t, delt)
			call pflush(luncoord_out)
		endif

		!velocity output
		if (momenta_out) then
			call save_XYZ(lunmomenta_out, PPc, Upot, read_method, t, delt)
		endif

		!dipoles output
		if (dip_out) then
			do iw=1,Nwaters
				write(lundip_out,'(4(1x,f12.4))') dip_momI(:,iw) , &
				dsqrt(dot_product(dip_momI(:,iw), dip_momI(:, iw)))
            enddo
			call pflush(lundip_out)
		endif!(dip_out)

		!electronic dipoles output
		if (Edip_out) then
			do iw=1,Nwaters
				do j = 1, 3
                    dip_momE(j,iw) = sum(dip_momEt(j,iw,:))/Nbeads
                enddo
             write(lunEdip_out,'(4(1x,f12.4))') dip_momE(:,iw) , &
					dsqrt(dot_product(dip_momE(:,iw), dip_momE(:, iw)))
			enddo
			call pflush(lunEdip_out)
		endif!(Edip_out)

		!charges out
		if (CHARGESOUT) then
			do iw = 1, 3*Nwaters
				write(lunCHARGESOUT,*) chg(iw)*0.20819434d0*DEBYE/CHARGECON
			enddo
		endif!(CHARGESOUT)

	endif!mod(t,t_out)==0

	!images output
	if (mod(t,ti_freq)  == 0) then

		if (images_out) then
			do i = 1, Nbeads
				call save_XYZ(lunimages_out, RRt(:,:,i), Upot, read_method, t, delt)
			enddo
			call pflush(lunimages_out)
		endif

		if (IMAGEDIPOLESOUT) then
			do i = 1, Nbeads
				do iw = 1,Nwaters
					write(lunIMAGEDIPOLESOUT,'(4(1x,f12.4))') dip_momIt(:,iw,i), &
				 	dsqrt(dot_product(dip_momIt(:,iw,i), dip_momIt(:,iw,i)))
				enddo
			enddo
		endif
	endif

	!total dipole moment output
	if (mod(t,td_freq)  == 0 .and. TD_out ) then
		write(lunTD_out,'(3f12.4)') dip_mom
		call pflush(lunTD_out)
  	endif

	if ((mod(t,checkpoint_freq) .eq. 0).and.(WRITECHECKPOINTS)) then
		call io_open(lun,'out_'//TRIM(fsave)//'_checkpoint_image.img',REPLACE=.true.)
		call save_checkpoint(lun, RRt, PPt, Upot, t,delt)
		call io_close(lun)
    endif

	if ((mod(t,hist_out_freq) .eq. 0).and.(HISTOUT).and.(t.gt.200)) then
		if (CALCGEOMETRY) call write_out_geometry(lunTP_out,Nbeads)
		call write_out_histogram(lunHISTOUT, Nbeads)
		if (CALCDIFFUSION .and. (t.gt.10)) call write_out_diffusion(lunTP_out, delt, fsave)
        if (CALCDOS .and. (t .gt. 10))    call calc_DOS(Hvelocities(:,1:tr,:),box,delt,fsave,avg_temp/Nbeads)
	endif

    !energy output to file
    if ( (mod(t,10) .eq. 0) .and. (ENERGYOUT) ) then
		write(lunENERGYOUT,'(1f16.4)') simple_energy
	endif

 endif!t .gt. eq_timesteps

!box size output
if (BAROSTAT) then
	!for accuracy, box size computed at every timestep
    sum_box  = sum_box + box
    sum_box2 = sum_box2 + box**2
    if (BOXSIZEOUT .and. (mod(t,t_freq) .eq. 0) ) write(lunBOXSIZEOUT,*) sum_box/t
endif

end subroutine write_out

!----------------------------------------------------------------------------------
!----------Print information about the run ---------------------------------------
!----------------------------------------------------------------------------------
subroutine print_run
  use dans_timer
  use geometry_calculator
  use diffusion_calculator
 implicit none

 !write Temp/Press file header
 call date_and_time(DATE=date,TIME=time)
 write(lunTP_out,'(a)') " call abort()ped on  "//date(1:4)//"-"//date(5:6)//"-"//date(7:8)//&
           " at "//time(1:2)//":"//time(3:4)//":"//time(5:10)
 call print_basic_run_info
 call print_timing_report(lunTP_out)

 call get_time("Total time", seconds)
 write(lunTP_out,'(a50, f10.2)') "ps/hour = ", (  (run_timesteps + eq_timesteps)*delt/seconds  )*3600
 write(lunTP_out,'(a50, f10.2)') "ps/day = ",  (  (run_timesteps + eq_timesteps)*delt/seconds  )*3600*24

 write(lunTP_out,*) "#------- thermodynamics: -------------------------------"
 avg_temp =  sum_temp/tr

 write(lunTP_out,'(a50, 3f10.2)') "Average temperature during run (K) = ", avg_temp/Nbeads
 write(lunTP_out,'(a50, 3f10.2)') "Average pressure during run (bar) = ", sum_press/tr
 write(lunTP_out,'(a50, 3f10.2)') "Average total energy during run (kcal/mole) = ", sum_tot_energy/tr
 write(lunTP_out,'(a50, 3f10.2)') "Estimated energy drift (kcal/mole/ps) = ",(tot_energy - init_energy)/(run_timesteps*delt)
 write(lunTP_out,'(a50, 3f10.2)') "Temp drift (K/ps) = ", (avg_temp - init_temp)/(run_timesteps*delt)
 write(lunTP_out,'(a50, 3f10.2)') "RMS energy fluctuation  (kcal/mol) = ", dsqrt( sum_tot_energy2/tr - (sum_tot_energy/tr)**2 )

 specific_heat = dsqrt(  sum_tot_energy2/tr - (sum_tot_energy/tr)**2  ) /( kb*avg_temp )
 write(lunTP_out,'(a50, f10.2)') "Specific heat C_V (only valid in NVT) (cal/g) = ", specific_heat/(1000*(massO+2*massH))

 if (Nbeads .gt. 1) then
	specific_heat = dsqrt(  sum_simple_energy2/tr - (sum_simple_energy/tr)**2  ) /( kb*avg_temp )
	write(lunTP_out,'(a50, f10.2)') "Quantum specific heat C_V (only valid in NVT) (cal/g) = ", specific_heat/(1000*(massO+2*massH))
 endif

 if (BAROSTAT) then
    avg_box2 = sum_box2/t
    avg_box  = sum_box/t

	write(lunTP_out,'(a50, 3f10.6)') "average box size (over entire run) (Ang) = ", avg_box
	write(lunTP_out,'(a50, f10.2)') "average density (g/cm^3) = ", &
		Nwaters*(massO+2*massH)*amu2grams/(avg_box(1)*avg_box(2)*avg_box(3)*(a2m*100)**3)

	isotherm_compress = &
		(avg_box2(1)**3 - (avg_box(1)*avg_box(2)*avg_box(3))**2 )*(10d-7)/(1.38d0*avg_temp*avg_box(1)*avg_box(2)*avg_box(3))

    write(lunTP_out,'(a50, f10.2)') "Isothermal compressibility =", isotherm_compress
    write(lunTP_out,'(a)') "Note : isothermal compressibility only valid in NPT, assumes cubic box, not validated"
 else
	write(lunTP_out,'(a50, f10.2)') "density (fixed) (g/cm^3) = ", Nwaters*(massO+2*massH)*amu2grams/(volume*(a2m*100)**3)
 endif

 if (CALCGEOMETRY)  call write_out_geometry(lunTP_out,Nbeads)
 if (CALCDIFFUSION .and. (run_timesteps.gt.10)) call write_out_diffusion(lunTP_out, delt, fsave)

 if (DIELECTRICOUT) then
    write(lunTP_out,'(a50 )')         "#---------- dielectric constant data  ---------------"
	write(lunTP_out,'(a50, f10.2)') " dielectric constant ", dielectric_constant
	if (run_timesteps .gt. 1000) then
		dielectric_error = sum( (dielectric_running(500:1000) - dielectric_constant)**2 ) /500
	  	write(lunTP_out,'(a50, f16.2)') "estimated error = +/-", dielectric_error
	endif
	write(lunTP_out,'(a50, f16.2)') " sum M^2 (Debye^2) ", sum_dip2
	write(lunTP_out,'(a50, f16.2)') " sum M (Debye^2) ", sum_dip
	write(lunTP_out,'(a50, f16.2)') " average <M^2> (Debye^2) ", sum_dip2/ttt
	write(lunTP_out,'(a50, f16.2)') " average <M>^2 (Debye^2) ", sum(sum_dip**2)/ttt
	write(lunTP_out,'(a50, i5)') " points used to compute dielectric constant: ", ttt
 endif

 if (CALCIRSPECTRA .and. (run_timesteps .gt. 10)) call calc_infrared_spectrum(dip_mom_all_times,box,delt,fsave,avg_temp/Nbeads)
 if (CALCDOS .and. (run_timesteps .gt. 10))       call calc_DOS(Hvelocities,box,delt,fsave,avg_temp/Nbeads)

 if (WRITECHECKPOINTS) then
	call io_open(lun,'out_'//TRIM(fsave)//'_finalchk.img')
	call save_checkpoint(lun, RRt, PPt, Upot, t, delt)
	call io_close(lun)
 endif

end subroutine print_run

!----------------------------------------------------------------------------------
!----------Print basic information about the run ---------------------------------
!----------------------------------------------------------------------------------
subroutine print_basic_run_info
 if (pot_model .eq. 2) write(lunTP_out,'(a50,a)') "Model = ", "TTM2F"
 if (pot_model .eq. 3) write(lunTP_out,'(a50,a)') "Model = ", "TTM3F"
 if (pot_model .eq. 4) write(lunTP_out,'(a50,a)') "Model = ", "qSPCfw"
 if (pot_model .eq. 5) write(lunTP_out,'(a50,a)') "Model = ", "SPCf"
 if (pot_model .eq. 6) write(lunTP_out,'(a50,a,a,a,a,i5)') "Model = ", trim(siesta_name), " sys_label = ", &
									trim(sys_label), " tot num SIESTA nodes = ", num_SIESTA_nodes
 write(lunTP_out,'(a50,i4,a,i4,a)') "Running with ", Nbeads, " beads on ", Nnodes, " nodes"
 write(lunTP_out,'(a50, i6)') "Num Molecules = ", Nwaters
 write(lunTP_out,'(a50, f10.3,a3)') "timestep = ", delt*1000, " fs"
 write(lunTP_out,'(a50, f8.4,a21,f8.4)') "mass of hydrogen = ", massH, "  mass of oxygen = ", massO
 write(lunTP_out,'(a50,a)') "type of run = ", trim(PIMD_type)
 if (THERMOSTAT) write(lunTP_out,'(a50, a )')  " type of global thermostat = ", "Nose-Hoover"
 if (.not. THERMOSTAT) write(lunTP_out,'(a50, a )')  " type of global thermostat = ", "none"
 if (.not. THERMOSTAT) write(lunTP_out,'(a50, a3)') "global thermostat tau = ", "n/a"
 if (THERMOSTAT) write(lunTP_out,'(a50, f10.3,a3)') "global thermostat tau = ", tau, " ps"
 if (BEADTHERMOSTAT) write(lunTP_out,'(a50, a )')  " type of bead thermostat = ", bead_thermostat_type
 if (CENTROIDTHERMOSTAT) write(lunTP_out,'(a50, a )')  " centroid thermostating = ", "yes"
 if (.not. CENTROIDTHERMOSTAT) write(lunTP_out,'(a50, a )')  " centroid thermostating = ", "no"
 if (.not. BEADTHERMOSTAT) write(lunTP_out,'(a50, a )')  " type of bead thermostat = ", "none"
 if (BEADTHERMOSTAT) write(lunTP_out,'(a50, f10.3,a3)')  " centroid thermostat tau = ", tau_centroid, " ps"
 if (.not. BEADTHERMOSTAT) write(lunTP_out,'(a50, a3)') " centroid thermostat tau = ", "n/a"
 if (BAROSTAT) write(lunTP_out,'(a50, f10.3,a3)') "Barostat tau = ", tau_p, " ps"
 if (.not. BAROSTAT) write(lunTP_out,'(a50, a3)') "Barostat tau = ", "n/a"
end subroutine print_basic_run_info

!----------------------------------------------------------------------------------
!----------Print information about the potential----------------------------------
!----------------------------------------------------------------------------------
subroutine print_pot(RR, Upot, dRR, virt, dip_momI, chg, lunTP_out)
use system_mod
use consts
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR
double precision, intent(in) :: Upot
double precision, dimension(3, Natoms), intent(in) :: dRR
double precision, dimension(3,3), intent(in) :: virt
double precision, dimension(3, Nwaters), intent(in) :: dip_momI
double precision, dimension(3 ) :: dip_mom
double precision, dimension( Natoms ), intent(in) :: chg
integer, intent(in) :: lunTP_out
integer :: iw

dip_mom(:) = sum(dip_momI(:, 1:Nwaters), dim=2)

write( lunTP_out,'(a50,f14.6)')"Potential Energy (kcal/mol) = ", Upot
write( lunTP_out,'(a50,f14.6)')"Estimated Enthalpy of Vaporization &
(kJ/mol) = ", Upot*4.184 + 1000*8.3144
write( lunTP_out,'(a50,f14.6)')"monomer Energy = ", Umon
write( lunTP_out,'(a50,f14.6)')"vdw Energy = ", Uvdw
write( lunTP_out,'(a50,f14.6)')"Long range vdw Energy = ", Uvdw_lrc
write( lunTP_out,'(a50,f14.6)')"electrostatic  Energy = ", Uelec
write( lunTP_out,'(a50,f14.6)')"induced Energy = ", Uind
write( lunTP_out,'(a50,3(1x,f12.3))')"Dipole moment [Debye] : ",dip_mom*DEBYE/CHARGECON
write( lunTP_out,'(a/3(10x,f14.5,1x))')"Virial tensor",virt(:,:)
!write(*,'(/a/)')"DERIVATIVES"
!do iw=1, Nwaters
!   write(*,'(a2,3x,3(f12.6,2x))')"O ",dRR(:, 3*iw-2)
!   write(*,'(a2,3x,3(f12.6,2x))')"H ",dRR(:, 3*iw-1)
!   write(*,'(a2,3x,3(f12.6,2x))')"H ",dRR(:, 3*iw-0)
!enddo
end subroutine print_pot


!-----------------------------------------------------------------------------------------
!--------------Write out coords----------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine save_XYZ(iun, RR, Upot, read_method, t, delt)
use system_mod
use consts
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR
double precision ::  delt, Upot
integer ::  i, iO, ih1, ih2, t
integer :: iun, read_method

write(iun,'(i10)') Natoms !, angle
write(iun,*) ""
!write(iun,'(f12.6,2x,f12.6,3(1x,f12.6))') t*delt, box
!if (read_method==0) then
!    do i=1, Nwaters
!       iO = 3*i-2
!       write(iun,'(a2,3(1x,f12.6))')'O ',RR(:, iO)
!    enddo
!    do i=1, Nwaters
!       ih1 = 3*i-1
!       ih2 = 3*i
!       write(iun,'(a2,3(1x,f12.6))')'H ',RR(:, ih1)
!       write(iun,'(a2,3(1x,f12.6))')'H ',RR(:, ih2)
!    enddo
! else if (read_method==1) then
   do i=1, Nwaters
      iO = 3*i-2
      ih1 = 3*i-1
      ih2 = 3*i
      write(iun,'(a2,3(1x,f12.6))')'O ',RR(:, iO)
      write(iun,'(a2,3(1x,f12.6))')'H ',RR(:, ih1)
      write(iun,'(a2,3(1x,f12.6))')'H ',RR(:, ih2)
   enddo
!endif
end subroutine save_XYZ


!-----------------------------------------------------------------------------------------
!--------------Write out configuration "image" ------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine save_checkpoint(iun, RRt, PPt, Upot, t, delt)
use system_mod
use consts
implicit none
double precision, dimension(3, Natoms,Nbeads), intent(in) :: RRt, PPt
double precision, intent(in) :: delt, Upot
integer ::  i, j, k, iO, ih1, ih2, t
integer, intent(in) :: iun

write(iun,'(2i10)') Natoms, Nbeads
write(iun,'(f12.6,2x,f12.6,3(1x,f12.6))') t*delt + init_time,  box
   do i=1, Nwaters
      iO = 3*i-2
      ih1 = 3*i-1
      ih2 = 3*i
	do j = 1, Nbeads
	      write(iun,'(a2,6(1x,f12.6))')'O ',RRt(:, iO,  j), PPt(:, iO,  j)*imassO
	enddo
	do j = 1, Nbeads
	 	write(iun,'(a2,6(1x,f12.6))')'H ',RRt(:, ih1, j), PPt(:, ih1, j)*imassH
	enddo
	do j = 1, Nbeads
	      write(iun,'(a2,6(1x,f12.6))')'H ',RRt(:, ih2, j), PPt(:, ih2, j)*imassH
	enddo
   enddo
 if (allocated(vxi_beads)) then
   do i = 1,size(vxi_beads,2)
		do j = 1, size(vxi_beads,3)
			do k = 1, size(vxi_beads,4)
				write(iun,*) vxi_beads(:,i,j,k)
			enddo
		enddo
	enddo
 endif
 if(allocated(vxi_global)) then
	write(iun,*) vxi_global
 endif
end subroutine save_checkpoint


!-----------------------------------------------------------------------------------------
!--------------load configuration "image" (header is loaded earlier) --------------------
!-----------------------------------------------------------------------------------------
subroutine load_checkpoint(iun)
 use system_mod
 use consts
 implicit none
 integer ::  i, j, k, iO, ih1, ih2
 integer, intent(in) :: iun

 do i=1, Nwaters
      iO = 3*i-2
      ih1 = 3*i-1
      ih2 = 3*i
	do j = 1, Nbeads
	      read(iun,'(a2,6(1x,f12.6))')  ch2, RRt(1:3, iO,  j), PPt(1:3, iO,  j)
	      PPt(1:3, iO, j) =  PPt(1:3, iO,  j)*massO
	enddo
	do j = 1, Nbeads
		read(iun,'(a2,6(1x,f12.6))') ch2 ,RRt(1:3, ih1, j), PPt(1:3, ih1, j)
		PPt(1:3, ih1, j) = PPt(1:3, ih1, j)*massH
	enddo
	do j = 1, Nbeads
	        read(iun,'(a2,6(1x,f12.6))') ch2 ,RRt(1:3, ih2, j), PPt(1:3, ih2, j)
		PPt(1:3, ih2, j) =  PPt(1:3, ih2, j)*massH
	enddo
 enddo
 if ((BEADTHERMOSTAT) .and. (bead_thermostat_type .eq. "Nose-Hoover")) then
	do i = 1,size(vxi_beads,2)
		do j = 1, size(vxi_beads,3)
			do k = 1, size(vxi_beads,4)
				read(iun,*) vxi_beads(:,i,j,k)
			enddo
		enddo
	enddo
 endif
 if (THERMOSTAT) then
	read(iun,*) vxi_global
 endif

 !calculate centroid positions
 RRc = sum(RRt,3)/Nbeads
 !calculate centroid momenta
 PPc = sum(PPt,3)/Nbeads
 !check PBCs
 call PBCs(RRt, RRc)

end subroutine load_checkpoint


!----------------------------------------------------------------------
!--------------Shutdown program --------------------------------------
!----------------------------------------------------------------------
subroutine shutdown
 use system_mod
 use fsiesta

 if (pid .eq. 0)  then
	call print_run
 endif

 if (pot_model == 6) call siesta_quit('all')

 call io_close_all() !close ALL open file units

 call MPI_Barrier(MPI_COMM_WORLD, ierr)
 call MPI_Finalize(ierr)

end subroutine shutdown





end module InputOutput
