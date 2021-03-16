!      Copyright (c) 2021 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

C***********************************************************************      
      subroutine reorder_array(n,m,a,lda,indx,tmp)
      implicit none

!     Reorder the columns of an array in place using an index array.
!     indx(j) is the index of the column in the array (the source column)
!     that goes in column j (the destination column).
!
!     As reordering is a permutation, each index occurs only once. This
!     means that the permutations exist in one or more loops, where the
!     columns are moved along to the next index in the loop. For example,
!     indx(1)=6, indx(6)=3, indx(3)=1. The first destination column (1)
!     is stored, then the source column (6) is copied to the destination.
!     The new destination column is now column 6, as it has already been
!     copied. The next source column (3) is copied into column 6, and 
!     column 3 is the new destination. For this column, the source column
!     is 1, which is in temporary storage, so it is copied in and the loop
!     is finished.
!     Tracking is done by negating the index, to indicate it's been
!     used, i.e. the column in the array has been copied into.
!
!     No move is done for a loop of 1 (for obvious reasons), and for a
!     loop of 2, a simple swap can be done efficiently. 
!
!     The number of moves is equal to m + l - s where m is the number
!     of columns, l is the number of loops of length > 2, and s is the
!     number of loops of length 1
!
!     Requires subroutine dswap from blas library
!
!     Note: this function could be modified to eliminate the typing of
!     a and tmp, so that arrays of any type could be reordered. To do
!     this, remove the first 'else if' branch below (which calls dswap),
!     remove the 'implicit none', and remove the type declarations for
!     a and tmp. This implies that a and tmp inherit their type from
!     the calling function. This adds a move for each loop of length 2.
!
!     Written by Kenneth G. Dyall. Last updated March 7, 2021

! arguments      
      integer n                 ! number of elements in column (rows)
      integer m                 ! number of columns
      integer lda               ! leading dimension of array a
      real*8 a(lda,m)           ! array to be reordered
      integer indx(m)           ! index array
      real*8 tmp(n)             ! temporary storage for column
! local variables
      integer i, j, k           ! column indices

      do k = 1,m
        j = indx(k)
        if (j .lt. 0) cycle
        if (j .eq. k) then      ! index points to the same column
          indx(k) = -indx(k)
        else if (indx(j) .eq. k) then ! Loop of 2, swap columns
          call dswap(n,a(1,j),1,a(1,k),1)
          indx(k) = -indx(k)
          indx(j) = -indx(j)
        else if (indx(k) .gt. 0) then ! Loop > 2, traverse loop         
          i = k
          tmp(1:n) = a(1:n,i)   ! Copy first column in loop to temp.
          do while (indx(i) .gt. 0)
            j = indx(i)
            if (j .eq. k) then
              a(1:n,i) = tmp(1:n) ! Loop ends; copy from temp into column
            else
              a(1:n,i) = a(1:n,j) ! Move column to new location
            end if
            indx(i) = -indx(i)
            i = j
          end do
        end if
      end do
      indx(1:m) = -indx(1:m)    ! Restore index array

      return
      end

C***********************************************************************            
      subroutine sort_vector(n,a,next,indx,direction)
      implicit none

!     Sort a vector in ascending or descending order, and return an index
!     array giving the location of each element.
!     A linked list is constructed in array next, in which the value of
!     next(i) is the index of the next element in the list.
!     The linked list is constructed in ascending order of values of a.
!     After constructing the list, it is traversed to generate an index
!     array, for reordering in ascending or descending order.
!
!     The number of comparisons required is at most n(n+1)/2, if the
!     vector is already in ascending order. However, it is n if the 
!     vector is in descending order (this is the minimum number). 
!
!     This subroutine does not reorder the elements of the vector, which
!     can be done with reorder_array.
!
!     Written by Kenneth G. Dyall. Last updated March 1, 2021

! argumnts
      integer n                 ! length of vector a
      real*8 a(n)               ! vector to sort
      integer next(n)           ! pointer to next element in sorted list
      integer indx(n)           ! index array constructed from list
      integer direction         ! +1 for ascending sort, -1 for descending
! local variables
      integer ifirst            ! first element of linked list
      integer i, l, m           ! loop or array indices
      logical larger            ! result of comparison
      integer mstart, mend      ! first and last indices of loop for index
      
      ifirst = 1
      do i = 1,n
        next(i) = 0
      end do
      do i = 2,n
        l = 0
        m = ifirst
        larger = a(i) .ge. a(m)
        ! Walk through ordering list to find location of current element
        do while (larger .and. next(m) .ne. 0)
          l = m
          m = next(l)
          larger = a(i) .ge. a(m)
        end do
        if (larger .and. next(m) .eq. 0) then
          ! Loop terminated at end of list. Append to list
          next(m) = i
        else
          ! Loop terminated on element comparison. Insert into list
          next(i) = m
          if (m .eq. ifirst) then
            ifirst = i
          else
            next(l) = i
          end if
        end if
      end do

!     Now construct index array for ascending or descending order

      if (direction .gt. 0) then
        mstart = 1
        mend = n
      else
        mstart = n
        mend = 1
      end if
      
      i = ifirst
      do m = mstart, mend, direction
        indx(m) = i
        i = next(i)
      end do

      return
      end
