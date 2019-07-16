!-----------------------------------------------------------------------------------
!-----------------Normal mode calculation module -----------------------------------
!-----------------This module is self contained ------------------------------------
!-----------------------------------------------------------------------------------
! Copyright (c) 2014 Daniel C. Elton 
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

module NormalModes
 use consts 
 Implicit none
 double precision, dimension(:), allocatable, save   :: A11, A12, A21, A22
 double precision, dimension(:), allocatable, save   :: omegalist, MassScaleFactor
 double precision, dimension(:,:), allocatable, save :: C 

 contains
!---------------------------------------------------------------------
!----------- Evolve one ring a time period delta --------------------
!---------------------------------------------------------------------
subroutine EvolveRing(RR, PP, Nbeads, mass)
 Implicit none
 double precision, dimension(3,Nbeads), intent(inout)  :: RR, PP
 double precision, dimension(3,Nbeads)  :: RRtr, PPtr 
 double precision, dimension(3) :: temp2
 double precision, intent(in)   :: mass
 integer, intent(in) :: Nbeads
 integer              :: i, j, k

!Transform into normal mode representation 
RRtr =  real2NM(RR,Nbeads)
PPtr =  real2NM(PP,Nbeads) 

!Propagate Normal Modes
!j = the normal mode index. The number of normal modes always equals the number of beads
do j = 1, Nbeads
	!The normal mode evolution is simply the evolution of a harmonic oscillator
	!The matrix elements of the propagation matrix are stored in Aij(j) (minus the mass terms)
	!this may not be the best way to do this , but ... I wanted to avoid matrix multiplication	
	temp2     =  A11(j)*RRtr(:,j)     		     + A12(j)*PPtr(:,j)/(mass*MassScaleFactor(j))
	PPtr(:,j) =  A21(j)*mass*MassScaleFactor(j)*RRtr(:,j) + A22(j)*PPtr(:,j)  
 	RRtr(:,j) = temp2
enddo

!Transform normal modes back into real space

RR =  NM2real(RRtr,Nbeads)
PP =  NM2real(PPtr,Nbeads) 

end subroutine EvolveRing

!--------------------------------------------------------------------------
!------------ Initialize Normal Mode matrices ----------------------------
!--------------------------------------------------------------------------
subroutine InitNormalModes(Nbeads,omegan,delta,setNMfreq, lunTP_out)
 use consts
 Implicit None
 double precision, intent(in) :: omegan, delta, setNMfreq
 double precision              :: omega
 integer, intent(in) 	 :: Nbeads 
 integer 			 :: i, j, k, l, lunTP_out

 allocate(A11(Nbeads)) 
 allocate(A12(Nbeads)) 
 allocate(A21(Nbeads)) 
 allocate(A22(Nbeads)) 
 allocate(omegalist(Nbeads))
 allocate(MassScaleFactor(Nbeads))
 allocate(C(Nbeads,Nbeads))

 if (mod(Nbeads,2) .eq. 0) then
 	l = (Nbeads-2)/2  !even Nbeads case
 else
 	l = (Nbeads-1)/2  !odd Nbeads case
 endif


!---------------- construct matrix -------------------------------------------- 
!The modes are stored in the following order: 
!--zero mode
!--negative modes 
!--positive modes
!--extra mode for when Nbeads is even
! This is maybe a bit confusing to index when generating the matrix, 
! etc, but is a nice ordering for later on.

 do i = 1, Nbeads  !i = column index
	C(1,i) = Sqrt(1d0/Nbeads) !zero mode

 	do k = 1, l
		C(k+1  ,i) = Sqrt(2d0/Nbeads)*Cos(k*i*2d0*PI/Nbeads) !-k modes 
		C(k+l+1,i) = Sqrt(2d0/Nbeads)*Sin(k*i*2d0*PI/Nbeads) !+k modes
	enddo

	!extra mode for even Nbeads case only
	if (mod(Nbeads,2) .eq. 0) then
		C(2*l+2,i) = Sqrt(1d0/Nbeads)*((-1)**i)
	endif 
 enddo


!-------------- RPMD frequencies (setNMfreq = 0) --------------------------------------
!-------------- This case is always done for reference --------------------------------
!correct omega=0 mode for the omega-> 0 limit
	omegalist(1) = 0d0
	A11(1) = 1d0
	A12(1) = delta
	A21(1) = 0d0
	A22(1) = 1d0

!negative k modes
 do k = 1,l
	omega = 2*omegan*Sin(PI*( abs(k)/real(Nbeads)) )  
    	omegalist(k+1) = omega
	A11(k+1) = Cos(omega*delta)
	A12(k+1) = (1/omega)*Sin(omega*delta)
	A21(k+1) = -omega*Sin(omega*delta)
	A22(k+1) = Cos(omega*delta)	
 enddo

!positive k modes
 do k = 1,l
	omega = 2*omegan*Sin(PI*( abs(k)/real(Nbeads)) )  
  	omegalist(k+l+1) = omega
	A11(k+l+1) = Cos(omega*delta)
	A12(k+l+1) = (1/omega)*Sin(omega*delta)
	A21(k+l+1) = -omega*Sin(omega*delta)
	A22(k+l+1) = Cos(omega*delta)	
 enddo

!extra mode for even Nbeads case only
 if (mod(Nbeads,2) .eq. 0) then 
	k = l+1
	omega = 2*omegan*Sin(PI*( abs(k)/real(Nbeads)) )  
    	omegalist(2*l + 2) = omega
	A11(2*l+2) = Cos(omega*delta)
	A12(2*l+2) = (1/omega)*Sin(omega*delta)
	A21(2*l+2) = -omega*Sin(omega*delta)
	A22(2*l+2) = Cos(omega*delta)	
 endif 

 MassScaleFactor = 1

if (setNMfreq .eq. 0) then
 	write(lunTP_out,*) "Running usuing RPMD. All beads have physical mass."
 	write(lunTP_out,*) "Normal mode frequencies: (cm^-1)"
 	do k = 1,Nbeads
		omega = omegalist(k) 
		write(lunTP_out,'(f10.2)')  33.33333d0*omega/(2*PI) 
	enddo
endif 


!------------------------- ACMD / CMD / PIMD case (adiabaticity > 1) ------------------------------------
if (.not. (setNMfreq .eq. 0)) then
!correct omega=0 mode for the omega-> 0 limit
	omegalist(1) = 0d0
	A11(1) = 1d0
	A12(1) = delta
	A21(1) = 0d0
	A22(1) = 1d0

	omega = (2*PI)*setNMfreq/33.33333d0 !conv. cm-1 -> 1/ps
	
	do i = 1,Nbeads
		if (omegalist(i) == 0) then 
			MassScaleFactor(i) = 1d0
		else
			MassScaleFactor(i) = (omegalist(i)/omega)**2
		endif
		!write(lunTP_out,*) "mass scale factor", i, " = ", massScaleFactor(i)
	enddo

	write(lunTP_out,*) "Adiabaticity = ", omegan/omega 
 	write(lunTP_out,'(a,f8.3,a8,f10.2,a6)') "All normal modes scaled to ", omega/(2*PI), " 1/ps = ", 33.33333d0*omega/(2*PI), " cm^-1"
	write(lunTP_out,'(a,f8.2,a)') "Timestep should probably not be larger than ", 	(((2*PI)/omega)/4d0)*1000, " fs"
	
  	omegalist(2:2*l+1) = omega
	A11(2:2*l+1) = Cos(omega*delta)
	A12(2:2*l+1) = (1/omega)*Sin(omega*delta)
	A21(2:2*l+1) = -omega*Sin(omega*delta)
	A22(2:2*l+1) = Cos(omega*delta)

!extra mode for even Nbeads case only
 if (mod(Nbeads,2) .eq. 0) then 
    	omegalist(2*l + 2) = omega
	A11(2*l+2) = Cos(omega*delta)
	A12(2*l+2) = (1/omega)*Sin(omega*delta)
	A21(2*l+2) = -omega*Sin(omega*delta)
	A22(2*l+2) = Cos(omega*delta)	
 endif

endif


!check for resonances 
 do i = 1, Nbeads
	if ( abs(1/delta - omegalist(i)/TWOPI)/(1/delta)  .lt. .01 ) then
		write(lunTP_out,'(a,i2,a)') "WARNING: normal mode frequency", i, " differs from 1/timestep by less than 1 percent"
		write(lunTP_out,'(a,f10.4,a,f10.4)') "WARNING: Normal mode frequency ", i, " = ", omegalist(i)/TWOPI, " 1/timestep = ", 1/delta
		write(lunTP_out,*) "WARNING: You may encounter resonances which cause the simulation to fail"
	endif
 enddo

end subroutine InitNormalModes

!-----------------------------------------------------------------------------------------
!---------------- generate one ring polymer sampled from the-----------------------------
!---------------- free ring distribution at temperature T -------------------------------
!-----------------------------------------------------------------------------------------
subroutine gen_rand_ring(RR,mass,temp,Nbeads)
 use math
 Implicit None
 double precision, parameter :: hbar=6.35077993041d0 !hbar in amu*ang^2/ps
 double precision, parameter :: KB_amuA2ps2perK = .831447148d0
 double precision, dimension(3,Nbeads), intent(out) :: RR
 double precision, dimension(3,Nbeads)	         :: RRtr
 double precision, intent(in) :: mass, temp
 double precision 		 :: std_dev, omega, omegan
 integer 			 		 :: j, k, l
 integer, intent(in) 		 :: Nbeads

 omegan = KB_amuA2ps2perK*temp*Nbeads/hbar

 do j = 1,Nbeads
		omega = omegalist(j)
		if (omega .eq. 0d0) then 
			RRtr(:,j) = 0d0
		else 
  			std_dev = Sqrt((omegan*hbar)/(mass*MassScaleFactor(j)*omega**2))  
			RRtr(:,j) = (/ rand_norm(std_dev), rand_norm(std_dev), rand_norm(std_dev) /)	
		endif	
 enddo

 RR = NM2real(RRtr, Nbeads) 

end subroutine gen_rand_ring


!-----------------------------------------------------------------------------------------
!---------------- generage one ring polymer momenta from the-----------------------------
!---------------- free ring distribution at temperature T -------------------------------
!-----------------------------------------------------------------------------------------
subroutine gen_rand_ring_momenta(PP,mass,temp,Nbeads)
 use math
 Implicit None
 double precision, parameter :: KB_amuA2ps2perK = .831447148d0
 double precision, dimension(3,Nbeads), intent(out) :: PP
 double precision, dimension(3,Nbeads)	         :: PPtr
 double precision, intent(in) :: mass, temp
 integer 			 :: j, k
 integer, intent(in) 	 :: Nbeads

 do j = 1,Nbeads
	do k = 1, 3
		PPtr(k,j) = rand_norm(Sqrt(KB_amuA2ps2perK*mass*MassScaleFactor(j)*temp))
	enddo
 enddo

 PP = NM2real(PPtr, Nbeads) 

end subroutine gen_rand_ring_momenta


	


!---------------------------------------------------------------------
!------------ generate initial coordinates for a pure normal mode ---
!------------ (extra code for testing / visualization) --------------
!---------------------------------------------------------------------
!
!subroutine MakeNormalMode(RR, i, Nbeads) 
! Implicit none
! double precision, dimension(3,Nbeads), intent(out) :: RR
! double precision, dimension(3,Nbeads)		 :: RRtr
! integer, intent(in)	::  i, Nbeads
! integer 		::  j,k
!
! RRtr = 0
! RRtr(:,i) = 6d0
!
! do j = 1,Nbeads
!	do k = 1, Nbeads
!		RR(1,j) = RR(1,j) + C(k,j)*RRtr(1,k)
!		RR(2,j) = RR(2,j) + C(k,j)*RRtr(2,k)
!		RR(3,j) = RR(3,j) + C(k,j)*RRtr(3,k)
!	enddo
! enddo
!
!end subroutine MakeNormalMode
!---------------------------------------------------------------------
!-------------------- convert real to normal mode coordinates -------- 
!---------------------------------------------------------------------
function real2NM(AA,Nbeads) 
 Implicit none
 double precision, dimension(3,Nbeads), intent(in)  ::AA
 double precision, dimension(3,Nbeads)  :: real2NM
 integer, intent(in) :: Nbeads
 integer :: j, k 
	real2NM = 0 
	do j = 1,Nbeads
		do k = 1, Nbeads
			real2NM(1,j) = real2NM(1,j) + C(j,k)*AA(1,k)
			real2NM(2,j) = real2NM(2,j) + C(j,k)*AA(2,k)
			real2NM(3,j) = real2NM(3,j) + C(j,k)*AA(3,k)
		enddo
	enddo
 end function real2NM

!---------------------------------------------------------------------
!-------------------- normal mode to real coordinates ----------------
!---------------------------------------------------------------------
function NM2real(AAtr, Nbeads) 
 Implicit none
 double precision, dimension(3,Nbeads), intent(in)  :: AAtr
 double precision, dimension(3,Nbeads)  :: NM2real
 integer, intent(in) :: Nbeads
 integer :: j, k 
	NM2real = 0 
	do j = 1,Nbeads
		do k = 1, Nbeads
			NM2real(1,j) = NM2real(1,j) + C(k,j)*AAtr(1,k)
			NM2real(2,j) = NM2real(2,j) + C(k,j)*AAtr(2,k)
			NM2real(3,j) = NM2real(3,j) + C(k,j)*AAtr(3,k)
		enddo
	enddo
end function NM2real

end module
