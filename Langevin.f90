!-----------------------------------------------------------------------------------
!----------------- Langevin Thermostat for the beads -------------------------------
!----------------- This module is self contained -----------------------------------
!------------------ Ref: JCP 133 124104 (2010) -------------------------------------
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
module Langevin 
 Implicit none
 double precision, dimension(:), allocatable, save :: c1, c2
 
contains

!--------------- Langevin subroutine ---------------------------------
subroutine Langevin_NM(PPt)
 use NormalModes
 use math
 use consts
 use system_mod !source of Nbeads
 Implicit none
 double precision, dimension(3,Natoms,Nbeads), intent(inout) :: PPt 
 double precision, dimension(3,Nbeads)  :: PPtr
 double precision :: sqrtmass
 integer :: i2, j, k
 
 do i2 = 1, Natoms

	!convert into normal mode space
	PPtr = real2NM(PPt(:,i2,:),Nbeads) 

	do j = 1, Nbeads 
		if (mod(i2+2,3) .eq. 0) then
			sqrtmass = sqrt(massO*MassScaleFactor(j)) 
		else 
			sqrtmass = sqrt(massH*MassScaleFactor(j))
		endif
		do k = 1, 3
			PPtr(k, j) = c1(j)*PPtr(k, j) + c2(j)*rand_norm(sqrtmass)
 		enddo 
	enddo
	
	!convert back into real space
	PPt(:,i2,:) = NM2real(PPtr,Nbeads) 

 enddo



end subroutine Langevin_NM


!--------------- Initialization --------------------------------------
subroutine Init_Langevin_NM(delt2, CENTROIDTHERMOSTAT, tau_bead, Nbeads, temp)
 use NormalModes
 use consts
 Implicit none
 double precision, intent(in) :: delt2, tau_bead, temp
 double precision :: gammak
 integer, intent(in) :: Nbeads
 integer :: i 
 logical, intent(in) :: CENTROIDTHERMOSTAT

 allocate(c1(Nbeads))
 allocate(c2(Nbeads))

 do i = 1, Nbeads
	if (i .eq. 1) then 
		if (CENTROIDTHERMOSTAT) then
			gammak = 1/tau_bead !zero mode case
 		else
			gammak = 0 !no centroid thermostat
		endif
 	else 
		gammak = 2*omegalist(i)! other modes
	endif

	c1(i) = dexp(-delt2*gammak)
	c2(i) = dsqrt(1 - c1(i)**2) 
	c2(i) = dsqrt(KB_amuA2ps2perK*temp)*c2(i) 

 enddo
	
end subroutine Init_Langevin_NM

end module
