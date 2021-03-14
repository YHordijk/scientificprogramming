!module for testing the sorting algorithm

!         data representation:

! The program takes two-dimensional real arrays. The data is given in files where in the first
! row the number of rows (n(1)) and columns (n(2)) is given, along with the parameter c indicating which column
! to sort. Following the header are 2*n(1) lines containing n(2) elements. The first n(1) rows are 
! the input, and the last n(1) rows are the truth as implemented in the numpy.argsort function for python (stable timsort).

!     example file:
! files have been generated using the python script
! n(1) n(2) c
! input(1,1)      input(1,2)     ...  input(1,n(2))
! input(2,1)      input(2,2)     ...  input(2,n(2))
! ...
! input(n(1),1)   input(n(1),2)  ...  input(n(1),n(2))
! truth(1,1)      truth(1,2)     ...  truth(1,n(2))
! truth(2,1)      truth(2,2)     ...  truth(2,n(2))
! ...
! truth(n(1),1)   truth(n(1),2)  ...  truth(n(1),n(2))

!       test-cases:
! test001: random real array sorted on first column                                                       (n(1) = 6,  n(2) = 3, c = 1)  expected to PASS
!          The goal is to test the ability of the implementation to sort along specified columns.

! test002: random real array sorted on second column                                                      (n(1) = 6,  n(2) = 3, c = 2)  expected to PASS
!          The goal is to test the ability of the implementation to sort along specified columns.

! test003: random real array sorted on third column                                                       (n(1) = 6,  n(2) = 3, c = 3)  expected to PASS
!          The goal is to test the ability of the implementation to sort along specified columns.

! test004: random integer array with elements from [0,4)                                                  (n(1) = 12, n(2) = 3, c = 1)  expected to FAIL
!          The goal is to test stability of the algorithm. This was done by creating an array of
!          integers with guaranteed repeats. By sorting with a known stable algorithm we can then 
!          see if the implemented algorithm recreates the result. If it does it is stable, if it
!          doesn't, it is unstable

! test005: random real array with elements in first column set to 0.0                                     (n(1) = 12, n(2) = 3, c = 1)  expected to PASS
!          The goal is to test stability when attempting to sort an array that is already sorted.
!          In this case done by setting all elements in the first column to 0.0

! test006: random real array with elements in first column already sorted                                 (n(1) = 12, n(2) = 3, c = 1)  expected to PASS
!          The goal is to test stability when attempting to sort an array that is already sorted.
!          In this case done by pre-sorting the first column

! test007: random real array which is presorted and flipped along first axis                              (n(1) = 12, n(2) = 3, c = 1)  expected to PASS
!          The goal is to test whether it can handle sorted but flipped arrays

! test008: implemented in fortran.                                                                        (n(1) = 12, n(2) = 3)  expected to PASS
!          Tests whether random real array sorted first on column 1, then column 2, then column 1 
!          gives good results

! test009: implemented in fortran.                                                                        (n(1) = 12, n(2) = 3)  expected to FAIL
!          Tests whether random real array sorted first on column 1, then column 2, then column 1 
!          gives good results. In this case the elements are rounded. It is expected that this will 
!          not give stable results for unstable algorithms such as shellsort. So it is expected to fail

module test

    use sortingmodule
    
    public :: test_sort
    private

    contains

    subroutine test_sort()
        call test_sort_from_file()
        call test_in_fort()
    end subroutine

    subroutine test_sort_from_file()
        real*8, allocatable :: input(:,:), truth(:,:)
        integer :: i, j, n(2), c, file_index 
        character(7) :: file
        logical :: exists, faulty

        print *, "starting test suite for subroutine sortingmodule"

        !loop over every testcase in file
        file_index = 1
        do
            !get filename
            if(file_index<10) then
                write(file, "(A,I1)") "test00", file_index
            else if(file_index<100) then
                write(file, "(A,I2)") "test0", file_index
            else 
                write(file, "(A,I3)") "test", file_index
            endif

            !check if the file exists
            inquire(file="src\tests\"//file, exist=exists)
            !close if non-existent
            if(.not.exists) exit

            !open the file
            open(file="src\tests\"//file, unit=21)

            file_index = file_index + 1

            !read in the parameters for the test and allocate arrays
            read(21,*) n(1), n(2), c
            allocate(input(n(1), n(2)), truth(n(1), n(2)))

            !read input one row at a time
            do i = 1, n(1)
                read(21, *) input(i,:)
            enddo

            !read truth
            do i = 1, n(1)
                read(21, *) truth(i,:)
            enddo

            !sort array using written algorithm
            call sort(input, c)

            !compare arrays
            do i = 1, n(1)
                do j = 1, n(2)
                    if (input(i,j).ne.truth(i,j)) faulty = .true.
                enddo
            enddo

            if(.not.faulty) print *, file, "    passed"
            if(faulty) print *, file, "    failed"
            faulty = .false.

            deallocate(input, truth)
        enddo

    end subroutine


    subroutine test_in_fort()
        real*8 :: test1(12,3), test2(12,3)
        logical :: faulty
        integer :: i, j, n(2)=(/12,3/)

        !implement test008
        faulty = .false.

        !randomize test1
        call random_number(test1)
        
        call sort(test1, 1)
        test2 = test1
        call sort(test1, 2)
        call sort(test1, 1)

        !compare arrays
        do i = 1, n(1)
            do j = 1, n(2)
                if (test1(i,j).ne.test2(i,j)) faulty = .true.
            enddo
        enddo

        if(.not.faulty) print *, "test008", "    passed"
        if(faulty) print *, "test008", "    failed"


        !implement test009
        faulty = .false.

        !randomize test1
        call random_number(test1)
        test1 = test1 * 3
        !round them to nint
        do i = 1, n(1)
            do j = 1, n(2)
                test1(i,j) = nint(test1(i,j))
            enddo
        enddo
        
        call sort(test1, 1)
        test2 = test1
        call sort(test1, 2)
        call sort(test1, 1)

        do i = 1, n(1)
            do j = 1, n(2)
                if (test1(i,j).ne.test2(i,j)) faulty = .true.
            enddo
        enddo

        if(.not.faulty) print *, "test009", "    passed"
        if(faulty) print *, "test009", "    failed"

    end subroutine


    







end module