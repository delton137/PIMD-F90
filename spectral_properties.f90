!---------------------- calc infrared spectrum --------------------------------
! this is a module for calculating the infrared spectrum and 
! velocity-velocity autocorrelation function at the end of a run
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
module spectral_properties
 use consts
 implicit none
 Integer 			 		 :: NumPointsOut = 350 !maximum number of points to use 
 double precision, parameter :: MaxFreqOut = 7000  !maximum frequency to plot (cm^-1)
 
 contains 
 
!-------------------------------------------------------------
!-- compute infrared spectrum and write out -----------------
!-------------------------------------------------------------
subroutine calc_infrared_spectrum(dip_moms, box, timestep, fsave, temp)
 use lun_management 
 implicit none
 real, dimension(:,:), intent(in) :: dip_moms 
 double precision, dimension(3), intent(in) :: box ! box size in Ang
 double precision, intent(in) :: timestep          ! times step in PS
 double precision, intent(in) :: temp	           ! temp in Kelvin
 character(len=*), intent(in) :: fsave             !file label 
 integer :: trun, tcor, ix, i, t, n, tread, lun
 double precision :: omega, IR, magn, avgMag, MinFreqOut, vol, PointsAvailable
 complex, dimension(:), allocatable :: aux1,vcross,ACF
 !double precision :: pi, A0,A1,A2,A3, window
 integer :: NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 tread = size(dip_moms,2)

 !find closest power of 2 less than tread
 trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  	 )

 allocate(aux1(0:2*trun-1)) !Auxiliary array, for fourier transform
 allocate(Vcross(0:2*trun-1)) ! Stores the cross correlation in reciprocal space
 allocate(ACF(0:2*trun-1))  ! stores  the autocorrelation Fucntion (real space)

 ACF=0 

 !find correlation function 
 do ix=1,3
	! Call to the direct FFT
    aux1(0:trun-1) = cmplx(dip_moms(ix,1:trun))
    n=size(aux1)
    call four1(aux1,n,-1)
    Vcross=aux1*conjg(aux1)          
    call four1(Vcross,n,1)
    ACF = ACF + Vcross
 enddo
 
 ! Save the correlation function to file
 call io_open(lun, "out_"//trim(fsave)//"_dip_corr_function.dat")
 do i = 1, tread
 	write(lun,*) i*timestep, real(ACF(i))!/real(ACF(1))
 enddo	
 call io_close(lun)

 !! Fourier transform the ACF
 aux1=cmplx(ACF)
 call four1(aux1,n,-1)
 n=size(aux1) !need this 
  
 !Save the IR spectrum in the file
 call io_open(lun, "out_"//trim(fsave)//"_IR_spectra.dat")

 MinFreqOut =  1d0/(timestep*ps2s*Cspeed*n)
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)
 
 if (PointsAvailable .lt. NumPointsOut) then 
	NumPointsOut = PointsAvailable
	NumAvgOver = 1
 else 
	NumAvgOver = floor(PointsAvailable/NumPointsOut)
 endif 
	
 !write(*,*) "#-------------infrared spectrum ------------------------"
 !write(*,*) "trun = ", trun
 !write(*,*) "n, tread = ",  n, tread
 !write(*,*) "MinFreqOut = ",  MinFreqOut
 !write(*,*) "MaxFreqOut = ",  MaxFreqOut
 !write(*,*) "points available = ", PointsAvailable
 !write(*,*) "Averaging over", NumAvgOver


do t = 0, NumPointsOut-1

  avgMag = 0 
  do i = 1, NumAvgOver
	omega=( t*NumAvgOver+i )/(timestep*n*ps2s) !get freq in 1/s (Hz
	
	avgMag = avgMag + (omega**2)*real(aux1(t*NumAvgOver+i)) 
	
! 	!avgMag = avgMag + omega*tanh(hbar*omega/(Kb*2.0d0*Temp))*real(aux1(t*NumAvgOver+i)) 
	!avgMag = avgMag + omega*(1d0 - dexp(-(hbar*omega/(Kb*Temp))))*real(aux1(t*NumAvgOver+i))
  enddo
  avgMag = avgMag/real(NumAvgOver)

  
  IR = (2d0*(3.14159d0)*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8*kb*Temp*eps0) 

  ! Use the prefactor with harmonic quantum correction (See Ramirez paper)
  !IR = (2d0*(3.14159d0)*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8*hbar) 
 
  !IR = 3.14159d0*(Debye2SI**2)*avgMag/(3d0*vol*hbar*2.99d8)
  
  IR = IR*ps2s ! convert units of dt from Fourier transform into seconds  
  
  IR = IR*.01 !convert 1/m to 1/cm

  IR = IR*2   !fudge factor -D. Elton ) 

  omega=( floor((t+.5)*numAvgOver) )/(timestep*n*ps2s) !get central freq in 1/s (Hz
  omega = omega/Cspeed  	! convert frequency to cm-1

  write(lun,*)  omega, IR
enddo

 call io_close(lun)

EndSubroutine calc_infrared_spectrum


!-------------------------------------------------------------------------
!- vel-vel ACF spectrum (aka "density of states") and write out 
!-------------------------------------------------------------------------
subroutine calc_DOS(Hvelocities, box, timestep, fsave, temp)
 use lun_management 
 use math
 implicit none
 real, dimension(:,:,:), intent(in) :: Hvelocities    ! Hvelocities stored in (3,Ntimesteps,Nhyd) array
 double precision, dimension(3), intent(in) :: box ! box size in Ang
 double precision, intent(in) :: timestep          ! times step in PS
 double precision, intent(in) :: temp	           ! temp in Kelvin
 character(len=*), intent(in) :: fsave             !file label 

 integer :: trun, tcor, ix, i, iH, t, n,  tread, lun, Nhyd, PointsAvailable	
 double precision :: omega, IR, magn, avgMag, MinFreqOut, vol
 double precision, dimension(:), allocatable :: ACF, output, DFT, allfreqs
 double precision, dimension(:), allocatable :: spectra_smoothed, allfreqs_smoothed
 complex, dimension(:), allocatable :: aux1
 integer :: NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 tread = size(Hvelocities,2)
 Nhyd  = size(Hvelocities,3)

 allocate(output(tread))! Stores the cross correlation in reciprocal space
 allocate(ACF(tread))   ! stores the autocorrelation Function (real space)
 allocate(allfreqs(tread))   
 allocate(DFT(tread))   ! stores the autocorrelation Function (real space)

 ACF=0  
 !find correlation function 
 do iH = 1, Nhyd
	do ix = 1, 3
		call calc_corr_function(Hvelocities(ix, 1:tread, iH), output)
		ACF = ACF + output
	enddo
 enddo
 ACF = ACF/Nhyd/ACF(1)
 
 !! Save the correlation function to file--------------------
 !call io_open(lun, "out_"//trim(fsave)//"_vel_vel_corr_function.dat")
 !do i = 1, tread
 ! 	write(lun,*) i*timestep, real(ACF(i))
 !enddo	
 !call io_close(lun)
 !-----------------------------------------------------------

 call calc_DFT(ACF, DFT, allfreqs, timestep, size(ACF))   

 allfreqs = allfreqs/(ps2s*Cspeed) !convert to cm^-1
 
 MinFreqOut =  allfreqs(2)
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)
  
 if (PointsAvailable .lt. NumPointsOut) NumPointsOut = PointsAvailable
  
 allocate(allfreqs_smoothed(NumPointsOut))
 allocate(spectra_smoothed(NumPointsOut)) 
  
 allfreqs_smoothed = block_average(allfreqs(1:PointsAvailable),  NumPointsOut)
 spectra_smoothed  = block_average(DFT(1:PointsAvailable), NumPointsOut)
 
 !! Save the DOS to file--------------------
 call io_open(lun, "out_"//trim(fsave)//"_DOS.dat")
 do i = 1, NumPointsOut
	write(lun,'(2g12.4)')  allfreqs_smoothed(i), spectra_smoothed(i)
 enddo
 call io_close(lun)
 !!----------------------------------------

EndSubroutine calc_DOS


!-------------------------------------------------------------------------
!- alternative infrared spectrum calculation 
!-------------------------------------------------------------------------
subroutine calc_infrared_spectrum2(dip_moms, box, timestep, fsave, temp)
 use lun_management 
 use math
 implicit none
 real, dimension(:,:), intent(in) :: dip_moms 
 double precision, dimension(3), intent(in) :: box ! box size in Ang
 double precision, intent(in) :: timestep          ! times step in PS
 double precision, intent(in) :: temp	           ! temp in Kelvin
 character(len=*), intent(in) :: fsave             !file label 

 integer :: trun, tcor, ix, i, iH, t, n,  tread, lun, Nhyd, PointsAvailable	
 double precision :: omega, IR, magn, avgMag, MinFreqOut, vol
 double precision, dimension(:), allocatable :: ACF, output, DFT, allfreqs
 double precision, dimension(:), allocatable :: spectra_smoothed, allfreqs_smoothed
 complex, dimension(:), allocatable :: aux1
 integer :: NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 tread = size(dip_moms,2)

 allocate(output(tread))! Stores the cross correlation in reciprocal space
 allocate(ACF(tread))   ! stores the autocorrelation Function (real space)
 allocate(allfreqs(tread))   
 allocate(DFT(tread))   ! stores the autocorrelation Function (real space)

 ACF=0  
 !find correlation function 
 do ix = 1, 3
	call calc_corr_function(dip_moms(ix, 1:tread), output)
	ACF = ACF + output
 enddo
 ACF = ACF/Nhyd/ACF(1)
 
 call calc_DFT(ACF, DFT, allfreqs, timestep, size(ACF))   

 allfreqs = allfreqs/ps2s !convert to Hz
 
 !apply quantum harmonic correction 
 do i = 1, PointsAvailable
	DFT(i) = allfreqs(i)*tanh(hbarSI*allfreqs(i)/(Kb*2.0d0*Temp))*DFT(i) 
 enddo
 !Use the prefactor with harmonic quantum correction (See Ramirez paper)
 DFT = (2d0*3.14159d0*(Debye2SI**2)*DFT)/(3d0*vol*2.99d8*hbarSI*Cspeed) 
  
 allfreqs = allfreqs/Cspeed !convert to cm^-1
 
 if (size(allfreqs) .gt. 1) MinFreqOut =  allfreqs(2)
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)
  
 if (PointsAvailable .lt. NumPointsOut) NumPointsOut = PointsAvailable
  
 allocate(allfreqs_smoothed(NumPointsOut))
 allocate(spectra_smoothed(NumPointsOut)) 
  
 allfreqs_smoothed = block_average(allfreqs(1:PointsAvailable),  NumPointsOut)
 spectra_smoothed  = block_average(DFT(1:PointsAvailable), NumPointsOut)
 
 !! Save the DOS to file--------------------
 call io_open(lun, "out_"//trim(fsave)//"_IR.dat")
 do i = 1, NumPointsOut
	write(lun,'(2g12.4)')  allfreqs_smoothed(i), spectra_smoothed(i)
 enddo
 call io_close(lun)
 !!----------------------------------------

EndSubroutine calc_infrared_spectrum2

 
EndModule spectral_properties