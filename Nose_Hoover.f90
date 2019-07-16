!---------------------------------------------------------------------------!
!------- General Nose-Hoover subroutine to propagate one chain  ------------
!------- 2014 Dan Elton. This subroutine is self-contained.     ------------
!------- It propagates the chain a half time step, so it should ------------
!------- be called twice per timestep.                          ------------
!------- It only works with chain_length=2 !!! 
!--------The part of the code for chain_length>2 has a bug which leads to temperature 
!------- oscillations and has been commented out 
!------- This code could be modified to support chains of length 1      
!------- but this does not reproduce the Canoncial ensembel & is thus prohibited      
!------- Ref:  "Understanding Molecular Simulation" by D. Frankel, Appendix E.2
!------- Ref for extension to arbitrary chain length:            ----------
!-------    Martyna, et. al. Mol. Phys. Vol. 87 5 1117 (1996)    -----------
!---------------------------------------------------------------------------!
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
subroutine Nose_Hoover(s, uk, chain_length, vxi, tau, delt2, L, temp)
 Implicit none 
 double precision, intent(out) 			    :: s   	      !variable to couple to
 double precision, intent(inout) 		    :: uk 	      !kinetic energy
 integer                 			        :: chain_length   !chain length
 double precision, dimension(chain_length), intent(inout) :: vxi !chain velocities
 double precision, intent(in) 		        :: tau 	      !kinetic energy
 double precision, intent(in) 		        :: delt2      !half timestep
 integer, intent(in) 			        :: L	      !# DOF being coupled to 
 double precision, intent(in) 		        :: temp	      !target temperature !DOUBLE PRECISION
 double precision 		 :: kT, delt4, delt8, M1, M2
 double precision, parameter  :: KB_amuA2ps2perK = .831447148d0
 double precision, parameter  :: hbar=6.35077993041d0 !hbar in amu*ang^2/ps
 integer 			  :: i

 !if (chain_length .lt. 2) then 
 !	 write(*,*) "ERROR: Nose Hoover chain length cannot be less than two."
 !	 stop 
 !endif
 
 
 !code only has been validated for chain length = 2
 !the part for higher chain lenght has been commented out
 !longer chains may need a different mass (M3, M4, etc) 
 
 chain_length = 2

 delt4 = delt2/2d0
 delt8 = delt2/4d0

 kT  = KB_amuA2ps2perK*temp 
 M1  = L*kT*tau**2 !tau = 1/omega
!! M2  = /(kT*(32/hbar)**2) ! = kt/omegan**2 massive thermostating scheme recommended by Tuckerman book assuming 32 beads here. (pg 476)
 M2 = M1!/L !Traditional thermostating 

 !if (chain_length .gt. 2) then
 !	i = chain_length
 !	vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4

 !	do i = chain_length - 1, 3
!		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
!		vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4
!		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
 !	enddo
 !	vxi(2)  = vxi(2)*dExp(-vxi(3)*delt8)
 !endif

 vxi(2) = vxi(2) + (M1*vxi(1)**2 - kT)*delt4/M2
! if (chain_length .gt. 2) vxi(2) = vxi(2)*dExp(-vxi(3)*delt8)
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)
 vxi(1) = vxi(1) + (2*uk - L*kT)*delt4/M1
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)

 !we don't bother calculating chain positions, as they are not needed
 
 s = dExp(-vxi(1)*delt2)

 !VV(i) = VV(i)*s
 uk    = uk*s*s
 
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)
 vxi(1) = vxi(1) + (2*uk - L*kT)*delt4/M1
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)

! if (chain_length .gt. 2) vxi(2) = vxi(2)*dExp(-vxi(3)*delt8)
 vxi(2) =  vxi(2) + (M1*vxi(1)**2 - kT)*delt4/M2

! if (chain_length .gt. 2) then
!	vxi(2)  = vxi(2)*dExp(-vxi(3)*delt8)
!
 !	do i = 3, chain_length - 1
!		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
!		vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4
!		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
!	enddo
!	 i = chain_length
 !	vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4
! endif 

end subroutine Nose_Hoover