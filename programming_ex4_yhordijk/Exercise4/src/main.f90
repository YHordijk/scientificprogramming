! program that calculates the integral of a function 
!(implemented in numericalintegration.f90 > evaluatefunction)
!It also finds values for m and n to converge to a certain eps (default 1e-4)
!Strategy is to increase m and n at a fixed rate (1.2) each cycle until convergence. 

!First attempt was steepest descent, but its compuationally expensive and does not work
!properly for integer type independent variables


program integration
    use NumericalIntegration

    implicit none

    integer :: mn(2), maxiter, i
    real*8 :: R, result1, result2, eps, mn_rate(2), delta
    real*8, ALLOCATABLE :: r_i(:), phi_j(:), w_i(:)

    !PARAMETER INITIALIZATION
    !input integration bound R
    write(6, fmt = "(a32)", advance="no") "Please input integration bound: "
    read(5, *) R

    !initialize starting values
    mn              = (/10, 10/)
    mn_rate         = (/1.2, 1.2/)
    eps             = 1e-4
    maxiter         = 100

    allocate(r_i(mn(1)))
    allocate(phi_j(mn(2)))
    allocate(w_i(mn(1)))

    print *, ' '
    print *, "cycle | prev int[f] |    delta    |     m |     n"
    print *, "-------------------------------------------------"

    do i = 1, maxiter
        !resize the arrays to the new m and n
        !we must first deallocate the arrays
        deallocate(r_i)
        deallocate(phi_j)
        deallocate(w_i)
        allocate(r_i(mn(1)))
        allocate(phi_j(mn(2)))
        allocate(w_i(mn(1)))

        !calculate the new grids
        call definegrid(R, mn, r_i, phi_j, w_i)

        !calculate the integral and store in result2
        result2 = evaluateintegral(r_i, phi_j, w_i)

        !check for convergence here and log the current data
        if (i>1) then
            delta = abs(result2 - result1)
            print 10, i, result1, delta, mn(1), mn(2)
            if (delta.le.eps) exit
        endif

        !update m and n
        mn = (/int(mn(1)*mn_rate(1))+1, int(mn(2)*mn_rate(2))+1/)
        result1 = result2
    enddo

    !print results
    print *, "===================="
    if(i.eq.(maxiter+1)) then 
        print *, "Optimization failed: Did not converge within iteration limit"
    else
        print *, "Optimization complete:"
    endif

    print *, "  cycles = ", i
    print *, "  int[f] =  ", result1
    write(6, fmt="(a12 e11.2)") "   delta  = ", delta
    print *, "  m      = ", mn(1)
    print *, "  n      = ", mn(2)
    
10 FORMAT(i6, ' | ', f11.5, ' | ', e11.5, ' | ', I5, ' | ', I5)
end program integration