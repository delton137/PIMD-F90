!-------------------------  m_RDF module --------------------------------------
!--- UNDER CONSTRUCTION / NOT COMPLETE 
!-- this is a self contained module for calculating RDFs during an MD run
!-- it is designed to be callable by any MD program 
!-- it currently contains subroutines for OO, HH, OH RDFs and average geometry 
!-----------------------------------------------------------------------------

module RDF_mod
Implicit None

real,dimension(3), save  :: length, halflength
Integer,save  :: Ndiv
real,save  :: delta, qO, qH, rOM
integer :: nsteps
integer,save  :: Na, Nmol
real,parameter :: e2coul=1.60217646e-19
real,parameter :: ang2m=1e-10
character(120) :: fileheader

real,dimension(:),allocatable :: gKr, gKr2, histOO, histOH, histHH,

 contains 
 
 
subroutine calc_geometry
 
 Nmol = 

 do i = 1, Nmol
	sum_HH  = 
	sum_OH  = 
	sum_OO  = 
	sum_HOH = 
 enddo 

endsubroutine 
 
 
!-----------------------------------------------------------------------
!------------------------ O-O RDF -------------------------------------
!-----------------------------------------------------------------------
  subroutine RDF_OO(xa)
    implicit none
    integer,intent(in) :: Nmol
    real,dimension(3,Nmol),intent(in) :: xa
    real,dimension(Ndiv),intent(out) :: hist
    
    ! INTERNAL VARS
    real :: rho, vol, distance, tmp, histgas
    integer :: ia, ja, ix, i, Ntot ! loop flags
    
    ! ARRAYS
    real,dimension(Ndiv) :: histliquid ! pre-normalization histogram

    histliquid = 0.0d0

    rho = Nmol / (length(1)*length(2)*length(3))

! ADJUST COORDINATES FOR PERIODIC BOUNDING
    do ia=1,Nmol-1                 ! do for every molecule
       do ja=ia+1,Nmol             ! to every other molecule
          distance=0.0d0            
          do ix=1,3
             if(abs(xa(ix,ia)-xa(ix,ja)) >= minval(halflength)) then
                if(xa(ix,ia)-xa(ix,ja) > 0) then        
                   tmp=xa(ix,ja)+minval(halflength)*2.0d0       
                else                                    
                   tmp=xa(ix,ja)-minval(halflength)*2.0d0       
                end if
                distance=distance+(xa(ix,ia)-tmp)**2  
             else                                       
                tmp=xa(ix,ja)                           
                distance=distance+(xa(ix,ia)-tmp)**2    
             endif
          enddo
          distance=sqrt(distance)                      
 ! ACTUALIZE HISTOGRAM
          if(distance .lt. minval(halflength)) then
             call updatehist(distance,histliquid)
          endif

       enddo
    enddo
! NORMALIZE HISTOGRAM
    do i=1,Ndiv    
       vol=(4.0d0/3.0d0)*pi*((i*delta)**3-((i-1)*delta)**3)            
       histgas=rho*vol 
       hist(i)= hist(i) + histliquid(i)/(histgas*Nmol)
    enddo
  end subroutine RDF_OO


!-------------------------------------------------------------------------------
!-------------------- OH RDF --------------------------------------------------
!-------------------------------------------------------------------------------
subroutine RDF_OH(Nmol, Oxy, Hydro, hist)
    implicit none

    ! DUMMY ARGUMENTS
    integer,intent(in) :: Nmol
    real,dimension(3,Nmol),intent(in) :: Oxy
    real,dimension(3,2*Nmol),intent(in) :: Hydro
    real,dimension(Ndiv),intent(inout) :: hist
    
    ! INTERNAL VARS
    real :: rho, vol, distance, tmp, histgas
    integer :: ia, ja, ix, i ! loop flags
    ! ARRAYS
    real,dimension(Ndiv) :: histliquid ! pre-normalization histogram

    ! INITIALIZE VARIABLES
    histliquid = 0.0d0
    rho = 2*Nmol / (length(1)*length(2)*length(3))

    ! ADJUST FOR PERIODIC BOUNDING AND CALCULATE DISTANCE
    do ia=1,Nmol
       do ja=1,2*Nmol 
          distance=0.0d0
          do ix=1,3
             if(abs(Oxy(ix,ia)-Hydro(ix,ja)) >= minval(halflength)) then 
                if(Oxy(ix,ia)-Hydro(ix,ja) > 0) then  
                   tmp=Hydro(ix,ja)+minval(halflength)*2.0d0  
                else                                  
                   tmp=Hydro(ix,ja)-minval(halflength)*2.0d0       
                end if
                distance=distance+(Oxy(ix,ia)-tmp)**2 
             else                                       
                tmp=Hydro(ix,ja)                            
                distance=distance+(Oxy(ix,ia)-tmp)**2   
             endif
          enddo
          distance=sqrt(distance)                     

 ! ACTUALIZE HISTOGRAM
          if(distance .lt. minval(halflength)) then
             call updatehist(distance,histliquid)     
          endif
       enddo
    enddo

! NORMALIZE HISTOGRAM
    do i=1,Ndiv   
       vol=(4.0d0/3.0d0)*pi*((i*delta)**3-((i-1)*delta)**3) 
       histgas=rho*vol
       hist(i)= hist(i) + histliquid(i)/(histgas*2*Nmol)
    enddo
    
end subroutine RDF_OH



!-----------------------------------------------------------------------
!------------- determine how histogram is binned ----------------------
!-----------------------------------------------------------------------
subroutine updatehist(x, hist)
    implicit none
    real :: x,lever
    real, dimension(Ndiv) :: hist
    integer :: i,j,count

     count=floor(x/delta)     
    !lever=mod(x,delta)/delta
    !there is a factor of two because each distance is for 2 atoms
    if (count .gt. 0) then
		histliquid(count) = hist(count) + 2
		!histliquid(count)=hist(count)+2.0d0*(1.0d0-lever)
		!histliquid(count+1)=hist(count+1)+2.0d0*lever
	endif
end subroutine updatehist

end module 
