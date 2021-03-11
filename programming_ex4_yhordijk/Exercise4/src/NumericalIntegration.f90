!module with functions and subroutines for calculating numerical integrals
!using the given method
!checkweights is called in the definegrid subroutine

module NumericalIntegration

    real*8 :: pi = 4 * atan(1.d0)

    public:: definegrid, evaluateintegral

    private 


    contains
    subroutine definegrid(R, mn, r_i, phi_j, w_i)
        real*8, intent(in) :: R
        integer, intent(in) :: mn(2)
        real*8, intent(out) :: r_i(mn(1)), phi_j(mn(2)), w_i(mn(1))
        real*8 :: w_pre, r_pre, phi_pre
        integer :: i, j

        !calculate some prefactors to speed up computation
        w_pre = (pi * R * R) / (mn(1)* mn(1) * mn(2))
        r_pre = R / mn(1)
        phi_pre = (2 * pi) / mn(2) 

        !get the r and w vectors
        do i = 1, mn(1)
            r_i(i) = r_pre * (i - 0.5)
            w_i(i) = w_pre * (2*i - 1)
        enddo

        !get the phi vector
        do j = 1, mn(2)
            phi_j(j) = phi_pre * (j - 0.5)
        enddo

        call checkweights(R, mn, w_i)
    end subroutine definegrid

    subroutine checkweights(R, mn, w_i)
        !the sum should be (pi * R * R / n), right?
        !$$\sum_{i=1}^m w(i) = \sum_{i=1}^m \frac{\pi R^2 (2i-1)}{m^2n} =\frac{\pi R^2 }{m^2n} \sum_{i=1}^m (2i-1) = \frac{\pi R^2 }{n}$$
        !since $$\sum_{i=1}^m (2i-1) = m^2$$

        reaL*8, intent(in) :: R
        real*8, intent(inout) :: w_i(:)
        integer, intent(in) :: mn(2)
        integer :: j
        real*8 :: sum

        sum = 0.0
        do j = 1, size(w_i)
            sum = sum + w_i(j)
        enddo

        !if sum is not correct, normalize the weights
        if (sum.ne.(pi * R * R / mn(2))) then
            ! print *, "Error, weights do not sum to pi * R * R / n: ", sum
            call normalizeweights(w_i, (pi * R * R / mn(2)) / sum)
        endif
    end subroutine checkweights

    subroutine normalizeweights(w_i, scaling)
        !Used for when the weights are not summed up to pi*R*R/n
        !normalizes using scaling factor
        real*8 :: scaling
        reaL*8, intent(inout) :: w_i(:)
        integer :: i

        do i = 1, size(w_i)
            w_i(i) = w_i(i) * scaling
        enddo
    end subroutine normalizeweights

    type(real*8) function evaluatefunction(r, phi)
        real*8, intent(in) :: r, phi

        evaluatefunction = sin(3*phi) * cos(phi/3) * exp((1-r)*(r-1))
    end function

    type(real*8) function evaluateintegral(r_i, phi_j, w_i)
        real*8, intent(in) :: r_i(:), phi_j(:), w_i(:)
        integer :: i, j

        evaluateintegral = 0.0
        do i = 1, size(r_i)
            do j = 1, size(phi_j)
                evaluateintegral = evaluateintegral + evaluatefunction(r_i(i), phi_j(j)) * w_i(i)
            enddo
        enddo
    end function 
end module