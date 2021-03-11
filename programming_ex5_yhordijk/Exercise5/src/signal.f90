module signal_processing

    public :: smooth, filter, compress
    private 

    contains

    subroutine smooth(signal, m, out)
        integer, intent(in) :: m
        real*8, intent(in) :: signal(:)
        real*8, intent(out) :: out(:)
        real*8, allocatable :: tempsig(:)
        integer :: i, j
        real*8 :: s

        !allocate temporary signal array
        allocate(tempsig(size(signal)))
        tempsig = signal

        if (m.le.0) then
            print *, "m must be an integer larger than zero"
            return
        endif

        do i = 0, size(signal)
            if(m < i .and. i < size(signal)-m) then
                s = 0.
                do j = -m, m
                    s = s + signal(i+j)
                enddo
                tempsig(i) = s/(2*m+1)
            endif
        enddo

        out = tempsig

    end subroutine smooth


    subroutine filter(signal, max_value, out)
        real*8, intent(in) :: max_value
        real*8, intent(in) :: signal(:)
        real*8, intent(out) :: out(:)
        integer :: i

        if (max_value.lt.0) then
            print *, "max_value must be a real larger than zero"
            return
        endif

        do i = 1, size(signal)
            out(i) = min(max_value, max(-max_value, signal(i)))
        enddo

    end subroutine filter


    subroutine compress(signal, threshold, indices)
        real*8, intent(in) :: threshold
        real*8, intent(inout) :: signal(:)
        integer, allocatable, intent(out) :: indices(:)
        integer, allocatable :: temp_indices(:)
        integer :: i, n

        !start by allocating indices the same amount of space as signal
        !resize later in the subroutine
        allocate(temp_indices(size(signal)))

        n = 0
        do i = 1, size(signal)
            if(signal(i) >= threshold) then
                n = n + 1 !keep track of number of points to be kept
                temp_indices(n) = i
                ! print *, temp_indices(n)
            endif
        enddo

        !resize indices array
        allocate(indices(n))
        indices = temp_indices(1:n)

    end subroutine compress

end module 