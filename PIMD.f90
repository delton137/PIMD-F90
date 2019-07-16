!----------------------------------------------------------------------------------
!----------------- PIMD / RPMD F90 program ---------------------------------------
!----------------------------------------------------------------------------------
!Changelog: 
! TTM3F code edited and expanded by Dan Elton 2/13
! Added PIMD parallelization 12/13
! Switched from velocity to momentum variables 3/14
! Addition of SIESTA force calculation 10/15
! new variable format input file 2/16
!
!To run the program use the following : 
! serial   : ./PIMD.x inputfilename > error_output_file
! parallel :  mpirun -np _ ./PIMD.x inputfilename > error_output_file
!
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
program PIMD
 use system_mod
 use consts
 use main_stuff
 use InputOutput
 use mpi
 use MD
 use dans_timer

 implicit none
!main variables can be found in main_stuff.f90

!---------------------- Read in input file ----------------------------------------
 call read_input_file

!---------------------- Start MPI and find number of nodes and pid --------------
 call MPI_Init(ierr)
 if (.not. (ierr .eq. 0)) write(*,*) "WARNING: MPI did not intialize correctly."
 call MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)

!---------------- Initializations ---------------------------------------------
 call read_and_initialize_all_nodes
 if (pot_model .eq. 6) then 
	call init_siesta
 else 
	call init_pot
 endif
 
!------------------Master node stuff --------------------------------------------
 if (pid .eq. 0) then
	call start_timer("Total time") 
	call open_files
	call master_node_init
	call read_coords_and_init  ! Read in coord data to RRc 
	call calc_uk 		       ! Initalize kinetic energy
 endif!(pid .eq. 0)

 !----------- MD Part of the program ---------------------------------------------
 call start_timer("MD")
 
 do t = 1, run_timesteps + eq_timesteps
	
	if (CONTRACTION) then
		call contracted_MD
	else if (MONOMERPIMD) then 
		call monomer_PIMD
	else
		call PIMD_VV
	endif
	
 enddo!t = 1, run_timesteps + eq_timesteps
 
 call stop_timer("MD")

 call shutdown 

end program PIMD