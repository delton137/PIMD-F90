!-------------------------------------------------------------------------------------
! Copyright (c) 2014-2016 Daniel C. Elton 
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
module consts
!double precision, parameter :: DENSCON = 1.66053872801494672163d0
!double precision, parameter :: CHARGECON = 18.22261720426243437986d0
 double precision, parameter :: CHARGECON = 18.22234397655801030455d0
 double precision, parameter :: DEBYE  = 4.8033324d0
 double precision, parameter :: MASSCON=418.4d0!at.weig->(kcal/mol)*(ps/Ang)^2 also, !kcal/(mol*ang) -> amu*ang/(particle*ps^2)
 double precision, parameter :: MASSCONi=1d0/418.4d0
 double precision, parameter :: PRESSCON2=69477d0   ! kcal/mole/Ang**3 to bar
 double precision, parameter :: PRESSCON =166.05387280d0  ! amu/(ps^2*ang) to bar
 double precision, parameter :: EVTOKCALPERMOLE = 23.0600d0 !ev to kcal/mol 
!double precision, parameter :: KCALKLVN=503.21659237841666799453d0
!double precision, parameter :: HBAR =0.0151704885d0 ! kcal/mol * ps
 double precision, parameter :: TEMPFACTOR =0.80181484527d0 ! factor to convert kinetic energy amu*Ang^2/ps^2 to K (using T = (2/3) 1/(k_B) * <uk>)
 double precision, parameter :: kb = 0.0019872041d0 !Boltzmann constant in kcal/mole/K  
 double precision, parameter :: KB_amuA2ps2perK = .831447148d0
 double precision, parameter :: hbar=6.35077993041d0 !hbar in amu*ang^2/ps
 double precision, parameter :: eps0 = 8.854187d-12 !SI units

 
 double precision, parameter ::     PI= 3.14159265358979323844d0
 double precision, parameter ::  TWOPI= 6.28318530717958647688d0
 double precision, parameter :: FOURPI=12.56637061435917295376d0
 double precision, parameter ::    TSP= 1.12837916709551257390d0 !  2/sqrt(pi)
 double precision, parameter ::   SQPI= 1.77245385090551602729d0 !  sqrt(pi)
 double precision, parameter :: RADIAN= 57.29577951308232087721d0

 double precision, parameter :: kbSI = 1.3806488d-23 ! J/Kelvin
 double precision, parameter :: vac_permSI = 8.854187817620d-12 ! F/m = C^2/J/m
 double precision, parameter :: a2m=1e-10 ! ang -> meters
 double precision, parameter :: amu2grams = 1.66053892d-24 ! amu -> grams
 double precision, parameter :: debyeSI=3.33564d-30  ! Debye -> SI
 double precision, parameter :: e2coul=1.60217646d-19 ! electron charge -> C
 double precision, parameter :: Cspeed=3.00d10 ! cm/s
 double precision, parameter :: hbarSI=1.0545718*10d-34  ! J*s
 double precision, parameter :: ps2s=1d-12 ! 1fs in s
 double precision, parameter :: Debye2SI = 3.33564d-30

end module consts

