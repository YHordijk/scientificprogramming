      program test_reorder
!     Test the sorting and reordering of a vector, ensuring that all
!     functionality is exercised
      integer n
      parameter (n=10)
      real*8 a(n),tmp(1)
      integer indx(n),next(n), indx_result(n)
      data a/1.0d0, 3.0d0, 6.0d0, 2.0d0, 9.0d0,
     &       7.0d0, 0.0d0, 5.0d0, 8.0d0, 4.0d0/
      data indx_result/7,1,4,2,10,8,3,6,9,5/
      integer i

      call sort_vector(n,a,next,indx,1)
      do i = 1,n
        if (indx(i) .ne. indx_result(i)) call exit(1)
      end do
      call reorder_array(1,n,a,1,indx,tmp)
      do i = 2,n
        if (a(i) .lt. a(i-1)) call exit(1)
      end do
      stop
      end
      subroutine dswap(n,a,ia,b,ib)
      integer i,n,ia,ib
      real*8 a(n),b(n)
      real*8 x
      do i = 1,n
        x = a(i)
        a(i) = b(i)
        b(i) = x
      end do
      return
      end

