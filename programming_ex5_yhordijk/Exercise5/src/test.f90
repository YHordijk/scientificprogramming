
!module for running tests on the three subroutines implemented in signal.f90

!       data representation:

!   filter and smooth:
!during filtering and smoothing data is created with size size(out_signal) == size(in_signal)
!since this is always the case (in current implementation) we can represent the data as pairs
!with elements corresponding to the test in_signal and ground-truth (as implemented in python)
!we can therefore read both at the same time
!the first row gives the parameters of the data and generation
!n is number of datapoints, m is either the number of smoothing neighbours or the maximum value
!for smoothing and filtering respectively

!   example file:
!n m
!in_1 truth_1
!in_2 truth_2
!in_3 truth_3
!...
!in_n truth_n


!   compress:
!during compression, data is lost, meaning that size(indices) <= size(in_signal)
!where indices is the resulting array from compression and signal is the input to the subroutine
!since the sizes differ, it is wise to use a contiguous stream of data, instead of pairs,
!where first the input signal is written followed by the ground truth data
!the first row now gives three parameters, first n - number of input points, k - number of indices, m - minimum compression value

!   example file:
!n k m
!in_1
!in_3
!in_3
!...
!in_n
!truth_1
!truth_2
!truth_3
!...
!truth_k


module test

use signal_processing

public :: test_compress, test_smooth, test_filter
private 

contains

!compress tests
subroutine test_compress()
    character(7) :: file
    logical :: exists, faulty=.false.
    real*8, allocatable :: input(:)
    integer, allocatable :: truth(:), test(:)
    integer :: file_index, input_n, truth_n
    real*8 :: thresh

    print *, "starting test suite for subroutine signal_processing.compress"

    !loop over every testcase
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
        inquire(file="src\tests\compress\"//file, exist=exists)
        !close if non-existent
        if(.not.exists) exit

        !open the file
        open(file="src\tests\compress\"//file, unit=21)

        file_index = file_index + 1

        !read in the parameters for the test and allocate arrays
        read(21,*) input_n, truth_n, thresh
        allocate(input(input_n))
        allocate(truth(truth_n))

        !read the data
        do i = 1, input_n
            read(21,*) input(i)
        enddo
        
        do i = 1, truth_n
            read(21,*) truth(i)
        enddo

        !calculate the compressed array 
        call compress(input, thresh, test)

        !loop over the test and truth arrays and compare
        do i = 1, truth_n
            if(test(i).ne.truth(i)) then
                ! print *, "Error found in value pair:", test(i), truth(i)
                ! print *, "in file", "src\tests\compress\"//file
                faulty = .true. 
            endif
        enddo

        if(.not.faulty) print *, file, "    passed"
        if(faulty) print *, file, "    failed"
        faulty = .false.
    
        !deallocate the arrays when done
        deallocate(input)
        deallocate(truth)
    enddo
end subroutine


!smoothing tests
subroutine test_smooth()
    character(7) :: file
    logical :: exists, faulty=.false.
    real*8, allocatable :: input(:), truth(:), test(:)
    integer :: file_index, input_n, m

    print *, "starting test suite for subroutine signal_processing.smooth"

    !loop over every testcase
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
        inquire(file="src\tests\smooth\"//file, exist=exists)
        !close if non-existent
        if(.not.exists) exit

        !open the file
        open(file="src\tests\smooth\"//file, unit=21)

        file_index = file_index + 1

        !read in the parameters for the test and allocate arrays
        read(21,*) input_n, m
        allocate(input(input_n))
        allocate(test(input_n))
        allocate(truth(input_n))

        !read the data
        do i = 1, input_n
            read(21,*) input(i), truth(i)
        enddo
        

        !calculate the compressed array 
        call smooth(input, m, test)

        !loop over the test and truth arrays and compare
        !we work with reals, so we must account for precision errors as well
        !(compare difference with epsilon)
        do i = 1, input_n
            if(abs(test(i)-truth(i)) > 1e-8) then
                print *, "Error found in value pair:",i,  test(i), truth(i)
                faulty = .true. 
            endif
        enddo

        if(.not.faulty) print *, file, "    passed"
        if(faulty) print *, file, "    failed"
        faulty = .false.
    
        !deallocate the arrays when done
        deallocate(input)
        deallocate(truth)
        deallocate(test)
    enddo
end subroutine


!filter tests
subroutine test_filter()
    character(7) :: file
    logical :: exists, faulty=.false.
    real*8, allocatable :: input(:), truth(:), test(:)
    integer :: file_index, input_n
    real*8 :: m

    print *, "starting test suite for subroutine signal_processing.filter"

    !loop over every testcase
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
        inquire(file="src\tests\filter\"//file, exist=exists)
        !close if non-existent
        if(.not.exists) exit

        !open the file
        open(file="src\tests\filter\"//file, unit=21)

        file_index = file_index + 1

        !read in the parameters for the test and allocate arrays
        read(21,*) input_n, m
        allocate(input(input_n))
        allocate(test(input_n))
        allocate(truth(input_n))

        !read the data
        do i = 1, input_n
            read(21,*) input(i), truth(i)
        enddo
        

        !calculate the compressed array 
        call filter(input, m, test)

        !loop over the test and truth arrays and compare
        !we work with reals, so we must account for precision errors as well
        !(compare difference with epsilon)
        do i = 1, input_n
            if(abs(test(i)-truth(i)) > 1e-8) then
                print *, "Error found in value pair:",i,  test(i), truth(i)
                faulty = .true. 
            endif
        enddo

        if(.not.faulty) print *, file, "    passed"
        if(faulty) print *, file, "    failed"
        faulty = .false.
    
        !deallocate the arrays when done
        deallocate(input)
        deallocate(truth)
        deallocate(test)
    enddo
end subroutine


end module test