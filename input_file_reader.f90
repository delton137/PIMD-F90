!-----------------------------------------------------------------------------------
! variable format ('split line') input file reader 
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
module input_file_reader
 implicit none
 integer :: MAX_LINE_LENGTH = 500

 contains 

!---------------------------------------------------------------
!---- Read input file 
!---------------------------------------------------------------
subroutine read_inputfile()
 use main_stuff
 implicit none 
 character(len=MAX_LINE_LENGTH) :: line, buffer
 character(len=100) :: key
 integer :: pos, lrerr
 logical :: is_comment
 
 is_comment(i) = (i .eq. 33) .or. (i .eq. 35) .or. (i .eq. 59)

 !call io_assign(lun) 
 !open(lun, file=finp, status='old')
 
 !--------------  default values for variables ----------------
 coord_out   = .false.     ! centroid coordinates output 
 momenta_out = .false.     ! centroid velocities output  
 dip_out     = .false.     ! all dipoles output 
 Edip_out    = .false.     ! electronic (polarization) dipole output
 TD_out      = .false.     ! total dipole ouput 
 images_out  = .false.    ! coordinates for all images output
 IMAGEDIPOLESOUT = .false. ! output dipoles for all PIMD images 
 BOXSIZEOUT  = .false.     ! box size running average output 
 TP_out      = .false.     ! Temp/Press (.true. = to file, .false. = to terminal)
 ENERGYOUT     = .true.   ! output total energy to seperate file
 HISTOUT       = .true.   ! output OH histograms to seperate file 
 hist_out_freq = 1000 
 
 CALCGEOMETRY  = .true. !computes averge geometry of h2o molecules and outputs at end
 CALCDIFFUSION = .true. !computes diffusion constant of oxygen atoms
 read_method   = 1 !read_method(=0,1) (0 for OOOO....HHHHH and 1 for OHHOHHOHH...)
 CALCIRSPECTRA = .false. !store dipole moments and calculate IR spectra at end of run
 CALCDOS       = .true. 
 
 CALC_RADIUS_GYRATION = .true. ! output avg radius of gyration (as column in Temp/Press file)
 DIELECTRICOUT = .true.   ! output dielectric constant from running averages <M^2> - <M>^2 (as column in Temp/Press file)
 CHARGESOUT    = .false.   ! write out charges on atoms in seperate file?
 WRITECHECKPOINTS = .true.! configurations of all beads to seperate file 
 td_freq = 1		   ! Total dipole output frequency 
 tp_freq = 10 		   ! Temp/Press output frequency
 ti_freq = 2000		   ! all images output frequency 
 t_freq  = 10          ! Output frequency for everything else
 checkpoint_freq = 2000! checkpoint output frequency
 
 siesta_name = 'siesta' !name of siesta executable 
 num_SIESTA_nodes = 1 !  # SIESTA processors 
 
 Rc = -1 		   !realspace cutoff (Ang)
 rc1 = -1	       !start of VdW force switch (Ang)
 eps_ewald = 1.d-6 !eps for aewald

 polar_maxiter = 15    !max polarization dipole iterations
 polar_sor     = 0.7   !TTM pol factor
 polar_eps     = 1.d-5 !accuracy to converge dipoles to 
 guess_initdip = .true.   !guess initial dipoles (this is normal for ttm calculations, but may cause errors in other cases) 
 print_dipiters = .false.  !print info about polarization dipole convergence for debugging 
 
 GENVEL = .true.	       !generate velocities 
 THERMOSTAT = .false.      ! Global Nose-Hoover thermostat? 
 BEADTHERMOSTAT = .false.  ! Bead thermostat? 
 CENTROIDTHERMOSTAT = .false. !thermostat bead centroid? 
 bead_thermostat_type = 'none' !type of bead thermostat - 'Nose-Hoover', 'Langevin' or 'none'
 tau = .1               !  tau for global thermostat (ps)
 tau_centroid = .1      !  tau for centroid thermostat (ps)
 global_chain_length = 2  ! global Nose-Hoover chain length
 bead_chain_length = 2    ! bead Nose-Hoover chain length
 temp = 300           ! Temperature (Kelvin) 
 BAROSTAT = .false.       ! Berendson barostat ?
 tau_P = .2           ! tau for barostat (ps)
 press = 1            ! reference pressure (bar) 
 PEQUIL = .false.		  ! pressure coupling during equilibration only?

 setNMfreq = 0        ! frequency (inverse cm) to scale normal modes to (0 for none/RPMD)
 PIMD_type = 'full'   ! type of PIMD run (full", "contracted", "monomerPIMD") 
 intra_timesteps = 10 ! ratio of slow timestep / fast timestep for contraction scheme
 massO = 15.994		  ! mass of Oxygen (au)
 massH = 1.008		  ! mass of Hyrdrogen
 RESTART = .false.    ! restart? (append to previous files) 

 SIESTA_MON_CALC = .false. !SIESTA monomer calculation 
 
 Nbeads       = 1      ! Number of beads 


 !---- Required variables ------------------------------------- 
 fconfig      = '-1'    ! input filename - either a centroid .xyz or a full bead configuration
 fsave        = '-1'    ! run name - this will be appended to all the output files
 eq_timesteps = -1      ! number of steps to equilibrate		 
 run_timesteps = -1      ! number of steps to run 
 delt      = -1 		! timestep in fs
 pot_model = -1 		! Model: 2=ttm21f 3=ttm3f 4= qspcfw 5=spcf 6=SIESTA
 
 !Required for SIESTA runs
 sys_label    = '-1'    !SIESTA label of .fdf file
 mon_siesta_name = '-1'

 !------------ read input file line by line --------------------
 do while (ierr .eq. 0)
	
	read(lun, '(A)', iostat=ierr) line  !read line
	line = adjustl(line) !remove leading spaces
 
    if (is_comment(ichar(line(1:1)))) goto 1000 !check if line is commented
	
	if (trim(line) .eq. '') goto 1000 !check if blank line
	
	!call make_lower(line) !make line lower case
	
	pos = scan(line, ' ')  !find position of first space
	key = line(1:pos)      !get key 
	buffer = line(pos+1:)  !store rest of line here

	select case (key)
	
        case('checkpoint_freq')
		READ(buffer, * ) checkpoint_freq
	
		case('CALC_RADIUS_GYRATION')
		READ(buffer, * ) CALC_RADIUS_GYRATION
		
		case('HISTOUT')
		READ(buffer, * )  HISTOUT
		
		case('ENERGYOUT')
		READ(buffer, * )  ENERGYOUT
		
		case('TP_out')
		READ(buffer, * ) TP_out
		
		case('BOXSIZEOUT')
		READ(buffer, * ) BOXSIZEOUT
		
		case('IMAGEDIPOLESOUT')
		READ(buffer, * ) IMAGEDIPOLESOUT
		
		case('images_out')
		READ(buffer, * ) images_out
		
		case('TD_out')
		READ(buffer, * ) TD_out
		
		case('Edip_out')
		READ(buffer, * ) Edip_out
		
		case('dip_out')
		READ(buffer, * ) dip_out
		
		case('momenta_out')
		READ(buffer, * ) momenta_out
		
		case('coord_out')
		READ(buffer, * ) coord_out

		case('DIELECTRICOUT')
		READ(buffer, * ) 

		case('CHARGESOUT')
		READ(buffer, * ) CHARGESOUT
		
		case('WRITECHECKPOINTS')
		READ(buffer, * ) WRITECHECKPOINTS

		case('td_freq')
		READ(buffer, * ) td_freq
		
		case('tp_freq')
		READ(buffer, * ) tp_freq
		
		case('ti_freq')
		READ(buffer, * ) ti_freq
		
		case('t_freq')
		READ(buffer, * )  t_freq

		case('siesta_name')
		READ(buffer, * ) siesta_name
		
		case('num_SIESTA_nodes')
		READ(buffer, * ) num_SIESTA_nodes
 
		case('Rc')
		READ(buffer, * ) Rc

		case('rc1')
		READ(buffer, * ) rc1
		
		case('eps_ewald')
		READ(buffer, * ) eps_ewald
		
		case('polar_maxiter')
		READ(buffer, * ) polar_maxiter
		
		case('polar_sor')
		READ(buffer, * ) polar_sor
		
		case('polar_eps')
		READ(buffer, * ) polar_eps

		case('guess_initdip')
		READ(buffer, * ) guess_initdip

		case('print_dipiters')
		READ(buffer, * ) print_dipiters

		case('GENVEL')
		READ(buffer, * ) GENVEL
		
		case('THERMOSTAT')
		READ(buffer, * ) THERMOSTAT
		
		case('BEADTHERMOSTAT')
		READ(buffer, * ) BEADTHERMOSTAT
		
		case('CENTROIDTHERMOSTAT')
		READ(buffer, * ) CENTROIDTHERMOSTAT
		
		case('bead_thermostat_type')
		READ(buffer, * ) bead_thermostat_type 

		case('tau')
		READ(buffer, *)  tau

		case('tau_centroid')
		READ(buffer, * )   tau_centroid
		
		case('setNMfreq')
		READ(buffer, * )  setNMfreq
		
		case('PIMD_type')
		READ(buffer, * )  PIMD_type
		
		case('intra_timesteps')
		READ(buffer, * )  intra_timesteps
	
		case('massO')
		READ(buffer, * )  massO
		
		case(' massH')
		READ(buffer, * )   massH
		
		case('press')
		READ(buffer, * ) press 
		
		case('tau_P')
		READ(buffer, * )  tau_P
		
		case('BAROSTAT ')
		READ(buffer, * )  BAROSTAT 
		
		case('temp')
		READ(buffer, * ) temp  
		
		case('bead_chain_length')
		READ(buffer, * )  bead_chain_length
		
		case('global_chain_length')
		READ(buffer, * )  global_chain_length

		case('RESTART')
		READ(buffer, * ) RESTART
		
		case('CALCGEOMETRY ')
		READ(buffer, * ) CALCGEOMETRY 
		
		case('CALCDIFFUSION')
		READ(buffer, * ) CALCDIFFUSION
		
		case('read_method')
		READ(buffer, * ) read_method
		
		case('CALCIRSPECTRA')
		READ(buffer, * ) CALCIRSPECTRA
		
		case('CALCDOS')
		READ(buffer, * ) CALCDOS
 
		case('SIESTA_MON_CALC')
		READ(buffer, * ) SIESTA_MON_CALC
		
		case('fconfig')
		READ(buffer, * ) fconfig
		
		case('Nbeads')
		READ(buffer, * ) Nbeads
		
		case('fsave')
		READ(buffer, * ) fsave
		
		case('eq_timesteps')
		READ(buffer, * ) eq_timesteps
		
		case('run_timesteps')
		READ(buffer, * ) run_timesteps
		
		case('delt')
		READ(buffer, * ) delt 
		
		case('pot_model')
		READ(buffer, * ) pot_model
		
		case('sys_label')
		READ(buffer, * )  sys_label
		
		case('mon_siesta_name')
		READ(buffer, * )  mon_siesta_name
		
		case('hist_out_freq')
		READ(buffer, * )  hist_out_freq
		
		case('old_style','OLD_STYLE','OLD','old') 
		rewind(lun) 
		read(lun,*)
		call read_old_style_input_file()
		checkpoint_freq = t_freq
		
	    case default	
		write(*,'(a)') "ERROR: could not understand the following line: ", buffer

	end select
 
1000  continue

 enddo

 !------- Check to make sure required variables were present  
 if (fconfig    .eq. '-1')  call io_abort("fconfig")       ! input filename - either a centroid .xyz or a full bead configuration
 if (Nbeads       .eq. -1)  call io_abort("Nbeads")        ! Number of beads 
 if (fsave        .eq. '-1')call io_abort("fsave")         ! run name - this will be appended to all the output files
 if (eq_timesteps .eq. -1)  call io_abort("eq_timesteps")  ! number of steps to equilibrate		 
 if (run_timesteps .eq. -1)  call io_abort("run_timesteps") ! number of steps to run 
 if (delt         .eq. -1)  call io_abort("delt")   	   ! timestep in fs
 if (pot_model    .eq. -1)  call io_abort("pot_model") 	   ! Model: 2.eq.ttm21f 3.eq.ttm3f 4.eq. qspcfw 5.eq.spcf 6.eq.SIESTA
 
 !Required for SIESTA runs
 if ( (pot_model .eq. 6) .and. (sys_label .eq. '-1')) call io_abort("sys_lable") !SIESTA label of .fdf file 
 if ( (pot_model .eq. 6) .and. (SIESTA_MON_CALC) .and. (mon_siesta_name .eq. '-1')) call io_abort("mon_siesta_name") 
 
end subroutine read_inputfile


!-----------------------------------------------------------------
!-- Read the old style of input file 
!-----------------------------------------------------------------
subroutine read_old_style_input_file()
 use main_stuff
 implicit none 
 read(lun,*) 
 read(lun,*)
 read(lun,*)
 read(lun,*)fconfig
 read(lun,*)fsave
 read(lun,*)Nbeads
 read(lun,*)eq_timesteps
 read(lun,*)run_timesteps
 read(lun,*)delt
 read(lun,*) 
 read(lun,*)
 read(lun,*)coord_out
 read(lun,*)momenta_out
 read(lun,*)dip_out
 read(lun,*)Edip_out
 read(lun,*)TD_out
 read(lun,*)images_out
 read(lun,*)IMAGEDIPOLESOUT
 read(lun,*)BOXSIZEOUT
 read(lun,*)TP_out
 read(lun,*)CALC_RADIUS_GYRATION
 read(lun,*)DIELECTRICOUT
 read(lun,*)CHARGESOUT
 read(lun,*)WRITECHECKPOINTS
 read(lun,*)td_freq
 read(lun,*)tp_freq
 read(lun,*)ti_freq
 read(lun,*)t_freq
 read(lun,*) 
 read(lun,*)
 read(lun,*)pot_model 
 read(lun,*)siesta_name, sys_label, num_SIESTA_nodes
 read(lun,*)Rc, rc1, eps_ewald
 read(lun,*)polar_maxiter, polar_sor, polar_eps, guess_initdip, print_dipiters
 read(lun,*)GENVEL
 read(lun,*)THERMOSTAT
 read(lun,*)BEADTHERMOSTAT
 read(lun,*)CENTROIDTHERMOSTAT
 read(lun,*)bead_thermostat_type
 read(lun,*)tau
 read(lun,*)tau_centroid
 read(lun,*)global_chain_length
 read(lun,*)bead_chain_length
 read(lun,*)temp
 read(lun,*)BAROSTAT
 read(lun,*)tau_P
 read(lun,*)press 
 read(lun,*)PEQUIL
 read(lun,*) 
 read(lun,*)
 read(lun,*)setNMfreq
 read(lun,*)PIMD_type
 read(lun,*)intra_timesteps
 read(lun,*)massO
 read(lun,*)massH
 read(lun,*)RESTART
end subroutine read_old_style_input_file




!-----------------------------------------------------------------
!--  Make a string lower case
!-----------------------------------------------------------------
subroutine io_abort(message)
	implicit none
	character(len=*), intent(in) :: message
	write(*,*) "ERROR: You are missing ", message, " in your input file!"
	stop 
end subroutine 

 
!-----------------------------------------------------------------
!--  Make a string lower case
!-----------------------------------------------------------------
subroutine make_lower(astr)  
 Implicit None
 character :: astr(*)
 integer   :: i, j
 logical   :: is_upper

 is_upper(i) = (i .ge. 65) .and. (i .le. 90)

 do i = 1, len(astr)
	j = ichar(astr(i)) 
	if (is_upper(j)) then
		astr(i) = achar(j + 32)
	endif
 enddo
 
EndSubroutine make_lower

!! Other one line character functions  
!is_digit(i) = (i .ge. 48) .and. (i .le. 57)
!is_upper(i) = (i .ge. 65) .and. (i .le. 90)
!is_lower(i) = (i .ge. 97) .and. (i .le. 122)
!is_alpha(i) = is_upper(i) .or. is_lower(i)
!is_alnum(i) = is_digit(i) .or. is_alpha(i)
!is_delstr(i)  = (i .eq. 34) .or. (i .eq. 39) .or. (i .eq. 96) !String delimiters: "  '  `
  
End Module input_file_reader