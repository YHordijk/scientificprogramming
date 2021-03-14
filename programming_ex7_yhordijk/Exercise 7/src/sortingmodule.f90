!implemented shellsort
!adapted from 1D version: https://rosettacode.org/wiki/Sorting_algorithms/Shell_sort#Fortran

module sortingmodule

    public :: sort, insertionSort, Shell_Sort
    private 

    interface sort 
        module procedure shellSort1D
        module procedure shellSort2D
    end interface

    contains 


    subroutine shellSort2D(array, col)
        !shell sort subroutine
        real*8, intent(inout) :: array(:,:)
        integer, intent(in), optional :: col
        
        integer :: i, j, h, n(2), c
        real*8, allocatable :: temp(:)

        !set c to 1 as default
        if (present(col)) then
            c = col
        else
            c = 1
        endif
        
        n = shape(array)

        !temp row vector
        allocate(temp(n(2)))

        h = n(1) / 2
        do while (h > 0)
            do i = h+1, n(1)
               j = i
               !must assign temp subarray since the array may be editted during the next 
               !while loop and array(i,:) may not give desired result
               temp = array(i,:)
               !here we compare values from the c-th column
               do while (j >= h+1 .and. array(j-h,c) > temp(c))
                !move the whole row:
                  array(j,:) = array(j-h,:)
                  j = j - h
               enddo
               array(j,:) = temp
            enddo
            if (h == 2) then
               h = 1
            else
               h = h * 5 / 11
            endif      
        enddo
    end subroutine shellSort2D


    subroutine shellSort1D(array)
        !shell sort subroutine
        real*8, intent(inout) :: array(:)

        integer :: i, j, h
        real*8 :: temp

        h = size(array) / 2
        do while (h > 0)
            do i = h+1, size(array)
               j = i
               temp = array(i)
               do while (j >= h+1 .and. array(j-h) > temp)
                  array(j) = array(j-h)
                  j = j - h
               enddo
               array(j) = temp
            enddo
            if (h == 2) then
               h = 1
            else
               h = h * 5 / 11
            endif      
        enddo
    end subroutine shellSort1D


    subroutine insertionSort(array, col)
        !insertion sort subroutine
        real(8), intent(inout) :: array(:,:)
        integer, intent(in), optional :: col

        real(8), allocatable :: temp(:)
        integer :: i, j, n(2), c

        if (present(col)) then
            c = col
        else
            c = 1
        endif

        n = shape(array)

        allocate(temp(n(2)))

        i = 2
        do while(i < n(1) + 1)
            j = i
            do while(j > 1 .and. array(j-1,c) > array(j,c))
                temp = array(j,:)
                array(j,:) = array(j-1,:)
                array(j-1,:) = temp
                j = j - 1
            enddo
            i = i + 1
        enddo
    end subroutine
end module 