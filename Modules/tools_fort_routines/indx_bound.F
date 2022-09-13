      subroutine indx_bound (x,n, x0,i)
      implicit none
      integer n, i                  ! bounded by
      real(kind=8) x(n), x0         !         x(i) <= x0 < x(i+1)

cf2py intent(in) x,n,x0
cf2py intent(out) i

      if (x0 < x(1)) then
        i=0                         ! if x0 is outside the full range
      elseif (x0 > x(n)) then       ! of x(1) ... x(n), then return
        i=n                         ! i=0 or i=n.
      else
        i=int( ( x(n)-x0 +n*(x0-x(1)) )/(x(n)-x(1)) )
        if (x(i+1)<x0) then
          do while (x(i+1) < x0)    ! This algorithm computes "i" as
            i=i+1                   ! linear interpolation between x(1)
          enddo                     ! and x(n) which should yield the
        elseif (x(i) > x0) then     ! correct value for "i" right a way
          do while (x(i) > x0)      ! because array elements x(i) are
            i=i-1                   ! equidistantly spaced.  The while
          enddo                     ! loops are here merely to address
        endif                       ! possible roundoff errors.

        if (x(i+1)-x0 < 0 .or. x0-x(i) < 0) then
          write(*,'(1x,A,5F12.6)') '### ERROR: indx_bound :: ',
     &                   x(i), x0, x(i+1), x0-x(i), x(i+1)-x0
          stop
        endif
      endif
      
      return
      end 

