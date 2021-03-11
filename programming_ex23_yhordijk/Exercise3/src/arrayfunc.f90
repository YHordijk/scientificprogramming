!module containing some array functions. Includes printing and getting a column from a 2d array
module arrayfunc

    implicit none

    private

    public :: arrayprint, vectorprint, getcolumn

    interface arrayprint
        module procedure arrayprint
    end interface
    interface vectorprint
        module procedure vectorprint
    end interface

    
    contains

    subroutine vectorprint(a)
        integer :: i
        real*8, intent(in) :: a(:)
    
        write(6, fmt="(a)", advance="no") '|'
        do i = 1, size(a, 1)
            write(6, fmt="(f8.3A)", advance="no") a(i), ' '
        end do 
        write(6,*) '|'
    end subroutine vectorprint

    subroutine arrayprint(a)
        integer :: i, j
        real*8, intent(in) :: a(:,:)

        do i = 1, size(a, 1)
            write(6, fmt="(a)", advance="no") '|'
            do j = 1, size(a, 2)
                write(6, fmt="(f8.3A)", advance="no") a(i,j), ' '
            end do 
            write(6,*) '|'
        end do
    end subroutine arrayprint

    function getcolumn(a, i)
        integer :: i, j
        real*8 :: a(:,:)
        real*8 :: getcolumn(size(a,1))
    
        do j = 1, size(a, 1)
            getcolumn(j) = a(j,i)
        end do
    end function getcolumn

end module arrayfunc

