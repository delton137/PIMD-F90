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
module system_mod
   !storage space for global variables

   logical :: guess_initdip
   logical :: print_dipiters

   integer :: pot_model !  =2 for ttm21f and =3 for ttm3f

   double precision, dimension(3)   :: box, boxi
   double precision, dimension(3,3) :: siesta_box

   double precision :: Umon, Uvdw, Uelec, Uind, Uvdw_lrc0, Uvdw_lrc

   double precision :: Rc, rc1, Rc2, volume, volume_init, eps_ewald
   double precision :: polar_sor, polar_eps
   integer :: polar_maxiter !changed to integer by D. Elton
   integer :: Natoms, Nwaters, Nbeads
   integer, dimension(:,:), allocatable :: neigh_list
   integer :: predict_step
   double precision :: const_ts1, const_ts2, const_ts3

   logical :: CONTRACTION, MONOMERPIMD
   character(len=200) :: sys_label 
   
   double precision :: massO = 15.994d0
   double precision :: massH = 1.008d0

end module system_mod
