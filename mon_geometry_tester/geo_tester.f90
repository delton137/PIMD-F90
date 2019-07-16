! monomer geometry tester - prints out energy and coords for monomers for potential energy surface fitting
!
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
program geo_tester
 use fsiesta
 implicit none
 integer, parameter :: num_SIESTA_nodes = 24
 character(100), parameter :: mon_siesta_name = "siesta3.2"
 character(100), parameter :: file_tag = "PBE_DZP"
! double precision, parameter :: minr=.2, rstep = .12 
! integer, parameter          :: mint=4,  maxt = 180, tstep = 4
 integer, parameter          :: nr = 13
 integer, parameter          :: nt = 6
 double precision, dimension(nr) :: r_values
 double precision, dimension(nt) :: t_values
 double precision, dimension(3,3) :: RR, dr1
 integer, parameter :: lunc = 22, lune = 23
 double precision :: r1, r2, r2x, r2y, e1, theta
 double precision, parameter :: Deg2Rad = 3.141592654d0/180d0 
 integer :: ri1, ri2, thi
 
 !launch SIESTA
 if (num_SIESTA_nodes .eq. 1) then  
    call siesta_launch( trim(mon_siesta_name), "monomer") !launch serial SIESTA process
 elseif (num_SIESTA_nodes .gt. 1) then  
    call siesta_launch( trim(mon_siesta_name), "monomer", nnodes=num_SIESTA_nodes ) !launch parallel SIESTA process  
 endif
 
 r_values = (/ .65, 0.75, 0.85, 0.95, .975, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5, 1.6, 1.7/) 
 t_values = (/ 85, 95, 100, 105, 110, 115/) 

 !Oxygen coord
 RR(:,1) = (/ 0.0, 0.0, 0.0 /) 

 !open output files
 open(lunc, file="out_"//trim(file_tag)//".ANI", status='unknown') 
 open(lune, file="out_"//trim(file_tag)//".energies", status='unknown') 
 
 write(*,*) "performing", (nr**2)*nt, " monomer calculations"
 
 !main loop 
 do ri1 = 1, nr
	r1 = r_values(ri1)
	do ri2 = 1, nr 
		r2 = r_values(ri2)
		do thi = 1, nt
            theta = t_values(thi)
            
			!set up coords
			RR(:,2) = (/ dble(r1), 0.0d0, 0.0d0 /) !1st Hydrogen on x axis
			r2x = r2*dcos(dble(theta)*Deg2Rad)
			r2y = r2*dsin(dble(theta)*Deg2Rad)
			RR(:,3) = (/ dble(r2x), dble(r2y), 0.0d0 /)
 
			!write(*,*) r1, r2, theta
			
		        call siesta_monomer(RR, dr1, e1)

			!write out
			write(lunc,'(i10)') 3  
			write(lunc,*) " "
			write(lunc,'(a2,3(1x,f12.6))')'O ', RR(:, 1)
			write(lunc,'(a2,3(1x,f12.6))')'H ', RR(:, 2)
			write(lunc,'(a2,3(1x,f12.6))')'H ', RR(:, 3)
			
			write(lune,'(f15.8)')  e1
		enddo
	enddo
 enddo

 call siesta_quit('all')
 
end program geo_tester


!---------------------------------------------------------------------
!-----------------SIESTA monomer force calculation ------------------
!---------------------------------------------------------------------
subroutine siesta_monomer(r1, dr1, e1)
 use fsiesta
 implicit none 
 double precision, dimension(3, 3), intent(in) :: r1
 double precision, dimension(3, 3), intent(out) :: dr1
 double precision, dimension(3, 3) :: siesta_box
 double precision, intent(out) :: e1
 double precision, parameter :: EVTOKCALPERMOLE = 23.0600d0 !ev to kcal/mol 

 siesta_box = 0.0
 siesta_box(1,1) = 10.0
 siesta_box(2,2) = 10.0
 siesta_box(3,3) = 10.0
  
 call siesta_forces("monomer", 3, r1, cell=siesta_box, energy=e1, fa=dr1)

 e1 = e1*EVTOKCALPERMOLE
 dr1 = -1d0*dr1*EVTOKCALPERMOLE    
  
end subroutine siesta_monomer

