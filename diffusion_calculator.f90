!---------------------- calc diffusion module ----------------------------------------
! this is a self contained module for calculating the diffusion constant 
! it does this in a very naive way, by making a plot of RMSD of oxygen atoms 
! and then fitting a straight line to it, beyond the 'balistic' regime 
! RMSD
!-------------------------------------------------------------------------------------
! Copyright (c) 2015 Daniel C. Elton 
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
module diffusion_calculator
 use system_mod !source of Nwaters, box, and boxi
 implicit none
 double precision, dimension(:), allocatable, save :: aMS
 integer, save :: ntimes
 integer, parameter :: t_out = 20 !how often to save RMS (for plotting)

 contains 
 
!-------------------------------------------------------------
!-- compute RMSD - to be called once per timestep
!-------------------------------------------------------------
subroutine calc_diff_RMSD(RRc, num_timesteps)
 implicit none
 integer, intent(in) :: num_timesteps 
 integer :: i, ia, ti
 double precision, dimension(3,Natoms), intent(in) :: RRc
 double precision, dimension(3,Nwaters) :: disps, Oxy2
 double precision, dimension(:,:), allocatable, save :: netdisp
 double precision, dimension(:,:), allocatable, save :: Oxy
 integer, save :: t=0 

 ntimes = floor(real(num_timesteps)/real(t_out))

 do ia = 1, Nwaters
	Oxy2(:,ia) = RRc(:,3*ia) 
 enddo

 t = t  + 1

 if (.not. allocated(aMS)) then 
    t=0
	allocate(netdisp(3,Nwaters))
	allocate(Oxy(3,Nwaters))
	allocate(aMS(ntimes))
	Oxy = Oxy2
 endif 
 
 !calculate all displacements
 disps = Oxy2 - Oxy 

 !Periodic boundary conditions 
 disps = disps - box(1)*anint(disps*boxi(1))

 !calculate net disps 
 netdisp = netdisp + disps 

 !calclate average RMS displacement
  if (mod(t,t_out) == 0) then
	ti  = int( t/t_out ) + 1
	Do i = 1, Nwaters
		aMS(ti) = aMS(ti) + sum( netdisp(:,i)**2 )
	enddo
	aMS(ti) = aMS(ti)/Nwaters
 endif
 
 !Copy current Oxygen coordinates to old 
 Oxy = Oxy2
 
endsubroutine calc_diff_RMSD

!-------------------------------------------------------------
!-- write out diffusion stuff -------------------------------
!-------------------------------------------------------------
subroutine write_out_diffusion(lun_out, timestep, fsave)
 use lun_management 
 implicit none 
 character(len=*) :: fsave
 double precision :: xbar,ybar,sumSQx,sumSQy, crossprod, b, a, errb, erra, D, errD, timestep
 real, dimension(:), allocatable :: times
 integer :: ntimesfit, ft, lun_out, iun_dif, i
 real, parameter :: frac=0.5!fraction of data to use
 
 allocate(times(ntimes))
 
 do i = 1, ntimes
	times(i) = i*t_out*timestep
 enddo

 !----------- Peform linear regression to last 50% of aMS(t) data ------------
 ft = floor(frac*ntimes)+1
 ntimesfit = ntimes - ft 

 xbar = sum(  times(ft:ntimes)/(  ntimesfit - 1) ) 
 ybar = sum(   aMS(ft:ntimes)/(  ntimesfit - 1) )
 
 sumSQx = sum( (times(ft:ntimes) - xbar)**2 ) 
 sumSQy = sum( (aMS(ft:ntimes) - ybar)**2  )
 crossprod = sum( (aMS(ft:ntimes) - ybar) * (times(ft:ntimes) - xbar) ) 

 b = crossprod/sumSQx
 a = ybar - b*xbar 

 errb = ( sumSQy/sumSQx - b**2 )/( ntimesfit - 2 ) 
 errb = sqrt(errb)

 D = 10*b/6 
 errD = errb/6 
 
 write(lun_out,*) "#---------- diffusion report ----------------------------" 
 write(lun_out,*) "Diffusion constant = ", D ," +/- ", errD, "10^(-5) cm^2/s"

 call io_assign(iun_dif)
 open(iun_dif,file ="out_"//trim(fsave)//"_diffusion.xvg",status="unknown")

 !-------------------Write out data + fit to file --------------------------
 !format for xmgrace
 !to plot use xmgrace -nosafe -nevice-nxy diffusion.xvg command
 write(iun_dif,*) '@ xaxis label char size 1.5'
 write(iun_dif,*) '@ yaxis label char size 1.5'
 write(iun_dif,*) '@ title "Mean square disps vs time" '
 write(iun_dif,*) '@ xaxis  label "t (ps)"  '
 write(iun_dif,*) '@ yaxis  label "<d^2>(t)"  '
 write(iun_dif,'(A16,f10.5,A5,f10.5,A17)') '@ subtitle "D = ', D ,' +/- ', errD, ' 10^(-5) cm^2/s "'
 do i  = 1,ntimes 
	write(iun_dif,'(3f16.4)') times(i), aMS(i), b*times(i) + a
 enddo

#ifdef FC_HAVE_FLUSH
 call flush(iun_dif) 
#endif

 call io_close(iun_dif)

EndSubroutine write_out_diffusion


endmodule diffusion_calculator