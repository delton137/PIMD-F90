!-------------------------------------------------------------------------------------
! portable flush()
!  
! a working flush() is critical when using UNIX pipes to communicate between programs
! 
! by default, ouput to files/pipes is buffered (often in chunks of 8192 bytes) 
! flush() forces i/o buffer to be written to disk
!
! the most common Fortran compiles (gfotran, ifort) have a working flush via "call flush()"
! in Fortran 2003, the intrinsic flush() is recommended 
!
! We also try to use the POSIX fsync function 
! if this causes an error, simply comment out the call to fsync at the bottom of the subroutine. 
!
! NAG (Numerical Algorithms Group) requires a special module to be called 
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
subroutine pflush(unit)
#ifdef __NAG__
      use f90_unix_io, only : flush
#endif
 implicit none 
 !interface for POSIX fsync function
 interface
	function fsync (fd) bind(c,name="fsync")
		use iso_c_binding, only: c_int
		integer(c_int), value :: fd
		integer(c_int) :: fsync
	end function fsync
 end interface 
 integer, intent(in) :: unit
 integer :: ret
 
#if defined(F2003) || defined(__F2003__)
      flush(unit)
#elif defined(GFORTRAN) || defined(__GFORTRAN__)
      call flush(unit)
#elif defined(XLF)
      if (unit.eq.6 .or. unit.eq.0) then
        call flush_(unit)
      else
        flush(unit)
      endif
#elif defined (FC_HAVE_FLUSH)
      call flush(unit)
#else
	!try anyway
	flush(unit)
    call flush(unit)
#endif

ret = fsync(fnum(unit))
if (ret /= 0) stop "see pflush.f90: Error calling POSIX FSYNC"

end subroutine pflush