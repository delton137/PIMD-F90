!---------------------- calc geometry module --------------------------------------
! this is a self contained module for calculating the avg. geometry during a run
! it is designed to be callable by any MD program 
!
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

module geometry_calculator
implicit none
integer, parameter          :: NUM_BINS = 80
double precision, parameter :: MIN_BIN  = .6  !Ang
double precision, parameter :: MAX_BIN  = 1.6 !Ang
double precision, parameter :: delta = (MAX_BIN-MIN_BIN)/NUM_BINS !Ang
double precision, save :: sum_HH=0,   sum_OH=0,  sum_HOH=0    !centroid-centroid distances
double precision, save :: sum_HH2=0,  sum_OH2=0, sum_HOH2=0   !centroid-centroid distances squared
double precision, save :: bsum_HH=0,  bsum_OH=0,  bsum_HOH=0  !bead-bead distances
double precision, save :: bsum_HH2=0, bsum_OH2=0, bsum_HOH2=0 !bead-bead distances squared
double precision, save :: max_OH=0, bmax_OH=0
integer, save :: num_total=0
integer, save :: num_warn=0
integer, parameter :: MAX_NUM_WARNINGS=20
double precision, dimension(NUM_BINS),save :: dOHhist=0, bdOHhist=0

 contains 
 
!-------------------------------------------------------------
!-- compute geometry ----------------------------------------
!-------------------------------------------------------------
subroutine calc_geometry(RR, RRt)
 use system_mod !source of Natoms, Nbeads, and boxi
 implicit none
 double precision, dimension(3,Natoms,Nbeads), intent(in) :: RRt ! 3xNatomsxNbeads
 double precision, dimension(3,Natoms), intent(in) :: RR ! 3xNatoms
 double precision, dimension(3) :: ROH1, ROH2, RHH
 double precision :: d_HH, d_OH1, d_OH2, d_HOH
 integer ::  ih1, ih2, io, i, j
 double precision, parameter :: dOHmax = 1.5d0
 
 do i = 1, Nwaters
	ih1 = 3*i
	ih2 = 3*i-1
	io = 3*i-2
	
	ROH1 = RR(:,io)  - RR(:,ih1)
	ROH1 = ROH1 - box*anint(ROH1*boxi)!PBC 

	ROH2 = RR(:,io)  - RR(:,ih2)
	ROH2 = ROH2 - box*anint(ROH2*boxi)!PBC
	
	RHH  = RR(:,ih2) - RR(:,ih1)
	RHH  = RHH - box*anint(RHH*boxi) !PBC
	
	d_HH  = dsqrt( dot_product(RHH,RHH) ) 
	d_OH1 = dsqrt( sum(ROH1**2) ) 
	d_OH2 = dsqrt( sum(ROH2**2) )  

	if (d_OH1 .gt. dOHmax) then 
		if (num_warn .lt. MAX_NUM_WARNINGS) write (*,*) "SEVERE WARNING centroid dOH > ", dOHmax, " = ", d_OH1
		num_warn = num_warn + 1
	endif
	if (d_OH2 .gt. dOHmax) then 
		if (num_warn .lt. MAX_NUM_WARNINGS) write (*,*) "SEVERE WARNING centroid dOH > ", dOHmax, " = ", d_OH2
		num_warn = num_warn + 1
	endif
	if (d_OH1 .gt. max_OH) max_OH = d_OH1
	if (d_OH2 .gt. max_OH) max_OH = d_OH1
	call binit(dOHhist, d_OH1)
	call binit(dOHhist, d_OH2)
	
	d_HOH = dacos( (d_OH1**2 + d_OH2**2 - d_HH**2)/(2d0*d_OH1*d_OH2) )
	
	sum_HH  = sum_HH  + d_HH    	
	sum_OH  = sum_OH  + d_OH1 + d_OH2
	sum_HOH = sum_HOH + d_HOH
	
	sum_HH2  = sum_HH2  + d_HH**2   	
	sum_OH2  = sum_OH2  + d_OH1**2 + d_OH2**2
	sum_HOH2 = sum_HOH2 + d_HOH**2

 enddo 
 
 !repeat calculations for all bead images 
  do j = 1, Nbeads 
	do i = 1, Nwaters
		ih1 = 3*i
		ih2 = 3*i-1
		io  = 3*i-2
	
		ROH1 = RRt(:,io,j)  - RRt(:,ih1,j)
		ROH1 = ROH1 - box*anint(ROH1*boxi)!PBC 

		ROH2 = RRt(:,io,j)  - RRt(:,ih2,j)
		ROH2 = ROH2 - box*anint(ROH2*boxi)!PBC

		RHH  = RRt(:,ih2,j) - RRt(:,ih1,j)
		RHH  = RHH - box*anint(RHH*boxi) !PBC
		
		d_HH  = dsqrt( sum(RHH**2) ) 
		d_OH1 = dsqrt( sum(ROH1**2) ) 
		d_OH2 = dsqrt( sum(ROH2**2) )  
		d_HOH = dacos( (d_OH1**2 + d_OH2**2 - d_HH**2)/(2.0*d_OH1*d_OH2) )
		
		if (d_OH1 .gt. dOHmax) then 
		if (num_warn .lt. MAX_NUM_WARNINGS) write (*,*) "WARNING bead centroid dOH > ", dOHmax, " = ", d_OH1
			num_warn = num_warn + 1
		endif
		if (d_OH2 .gt. dOHmax) then 
			if (num_warn .lt. MAX_NUM_WARNINGS) write (*,*) "WARNING bead centroid dOH > ", dOHmax, " = ", d_OH2
			num_warn = num_warn + 1
		endif
		
		if (d_OH1 .gt. bmax_OH) bmax_OH = d_OH1
		if (d_OH2 .gt. bmax_OH) bmax_OH = d_OH2
		call binit(bdOHhist, d_OH1)
		call binit(bdOHhist, d_OH2)
	
		bsum_HH  = bsum_HH  + d_HH    	
		bsum_OH  = bsum_OH  + d_OH1 + d_OH2
		
		if(abs(d_HOH) .gt. 0) bsum_HOH = bsum_HOH + d_HOH
	
		bsum_HH2  = bsum_HH2  + d_HH**2  	
		bsum_OH2  = bsum_OH2  + d_OH1**2 + d_OH2**2
		bsum_HOH2 = bsum_HOH2 + d_HOH**2

	enddo 
 enddo
 
 if (num_warn .eq. MAX_NUM_WARNINGS) write (*,*) "Max dOH warnings reached : supressing further warnings."
 
 num_total = num_total + Nwaters

endsubroutine calc_geometry



!-------------------------------------------------------------
!-- write out geometry --------------------------------------
!-------------------------------------------------------------
subroutine write_out_geometry(iun, Nbeads)
 implicit none 
 integer, intent(in) :: iun, Nbeads
 double precision :: bnum_total
 double precision :: avg_OH, avg_HH, avg_HOH
 double precision :: avg_OH2, avg_HH2, avg_HOH2, RMS_OH, RMS_HH, RMS_HOH
 double precision :: bavg_OH, bavg_HH, bavg_HOH
 double precision :: bavg_OH2,bavg_HH2,bavg_HOH2, bRMS_OH, bRMS_HH, bRMS_HOH
 integer :: i
 
 bnum_total = num_total*Nbeads
 
 avg_HH   = sum_HH/num_total 
 avg_OH   = sum_OH/(2.0*num_total)
 avg_HOH  = (sum_HOH/num_total)*180d0/3.141592d0
 avg_HH2  = sum_HH2/num_total
 avg_OH2  = sum_OH2/(2.0*num_total)
 avg_HOH2 = (sum_HOH2/num_total)*(180d0/3.141592d0)**2

 RMS_HH  = dsqrt(avg_HH2 - avg_HH**2) 
 RMS_OH  = dsqrt(avg_OH2 - avg_OH**2) 
 RMS_HOH = dsqrt(avg_HOH2 - avg_HOH**2) 
 
 
 bavg_HH   = bsum_HH/bnum_total 
 bavg_OH   = bsum_OH/(2.0*bnum_total)
 bavg_HOH  = (bsum_HOH/bnum_total)*180d0/3.141592d0
 bavg_HH2  = bsum_HH2/bnum_total
 bavg_OH2  = bsum_OH2/(2.0*bnum_total)
 bavg_HOH2 = (bsum_HOH2/bnum_total)*(180d0/3.141592d0)**2

 bRMS_HH  = dsqrt(bavg_HH2  - bavg_HH**2) 
 bRMS_OH  = dsqrt(bavg_OH2  - bavg_OH**2) 
 bRMS_HOH = dsqrt(bavg_HOH2 - bavg_HOH**2) 

 write(iun,'(a)')         "#---------- h2o geometry report -------------------------"
 write(iun,'(a)')  "                                 centroid-centroid                | bead - bead  (Ang,deg) " 
 write(iun,'(a40,f12.3,a4,2f12.3,a4,f12.3)') "average HH distance : ", avg_HH, " +/-", RMS_HH, bavg_HH, " +/-", bRMS_HH  
 write(iun,'(a40,f12.3,a4,2f12.3,a4,f12.3)') "average OH distance : ", avg_OH, " +/-", RMS_OH, bavg_OH, " +/-", bRMS_OH  
 write(iun,'(a40,f12.3,a4,2f12.3,a4,f12.3)') "  average HOH angle : ", avg_HOH, " +/-", RMS_HOH, bavg_HOH, " +/-", bRMS_HOH  
 write(iun,'(a40,f12.3,a11,f12.3)')           " max OH distance    : ", max_OH, " max bead OH = ", bmax_OH
 !write(iun,'(a)')         "#---------- final histograms -------------------------"
 !call write_out_histogram(iun, Nbeads)
 
endsubroutine write_out_geometry

!-------------------------------------------------------------
!-- write out histogram to seperate file --------------------
!-------------------------------------------------------------
subroutine write_out_histogram(iun, Nbeads)
 implicit none
 integer, intent(in) :: iun, Nbeads
 integer :: i 

!write(iun,'(a40,f12.6,a3)') "histogram bin size = ", delta, " Ang"
 write(iun,'(a80)') "dist (Ang)    dOH histograms (centroid-centroid , bead-bead)"
 do i = 1, NUM_BINS
	write(iun,'(3f16.10)') MIN_BIN + delta*i, dOHhist(i)/(2*num_total), bdOHhist(i)/(2*num_total*Nbeads)
 enddo
	
 call flush(iun) 

end subroutine write_out_histogram

!-------------------------------------------------------------
!-- reset geometry calculation ------------------------------
!-------------------------------------------------------------
subroutine reset_geometry
 implicit none 
 sum_HH=0;    sum_OH=0;   sum_HOH=0  
 sum_HH2=0;   sum_OH2=0;  sum_HOH2=0 
 bsum_HH=0;   bsum_OH=0;  bsum_HOH=0  
 bsum_HH2=0;  bsum_OH2=0; bsum_HOH2=0; num_total=0
 max_OH=0;    bmax_OH=0;  bdOHhist=0;  dOHhist=0 
endsubroutine

!-------------------------------------------------------------
!-- histogram binning ---------------------------------------
!-------------------------------------------------------------
subroutine binit(hist, thing_to_bin)
 implicit none 
 double precision, dimension(:), intent(inout) :: hist
 double precision, intent(in) :: thing_to_bin
 double precision :: lever, remain
 integer :: bin
  
 if ((thing_to_bin .gt. MIN_BIN) .and. (thing_to_bin .lt. MAX_BIN)) then

	bin = anint( (thing_to_bin-MIN_BIN)/delta)
	!remain = mod( (thing_to_bin-MIN_BIN), delta)

	if (bin .eq. 0) bin = 1
	if (bin .gt. NUM_BINS) bin = NUM_BINS
	
	!!--- bining with splitting between two closest bins ----
	!lever = remain/delta
	!hist(bin)  		        		 = hist(bin) 	+ 1 - lever
	!hist(bin+int(sign(1.d0,lever)))  = hist(bin+1) + remain/delta

	!!--- normal bining 
	hist(bin) = hist(bin) + 1
	
 endif 
 
endsubroutine binit

endmodule geometry_calculator