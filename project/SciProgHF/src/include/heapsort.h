      interface heapsort
!----------------------------------------------------------------------
! for integers
!----------------------------------------------------------------------
       subroutine heapsorti (iarr,ndim,linvrt,imask,nmask)
        implicit none
        integer,intent(in) :: ndim
        integer,intent(inout) :: iarr(ndim)
        logical,intent(in),optional :: linvrt
        integer,intent(in),optional :: nmask
        integer,intent(inout),optional :: imask(*)
       end subroutine heapsorti

!----------------------------------------------------------------------
! for double precisision real
!----------------------------------------------------------------------
       subroutine heapsortd (darr,ndim,linvrt,imask,nmask)
        integer,intent(in) :: ndim
        real(8),intent(inout) :: darr(ndim)
        logical,intent(in),optional :: linvrt
        integer,intent(in),optional :: nmask
        integer,intent(inout),optional :: imask(*)
       end subroutine heapsortd

      end interface heapsort
