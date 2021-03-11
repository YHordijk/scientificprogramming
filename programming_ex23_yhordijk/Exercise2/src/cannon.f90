program cannon
    !import integration subroutine, and a particle custom type
    use pkg

    real*8 :: g, angle, c, v_init, dt, t, ymax

    type(particle) :: part

    !read user input from special file unit 5 (stdin)
    100 print *, "Input an angle between 0 and 90 degrees:" ! in degrees
    read(5,*) angle

    if (0 > angle .or. angle > 90) then 
        print *, "Error: Input outside allowed range"
        go to 100
    end if 

    angle = angle * pi/180.d0 !convert to radians

    !some constants
    g = 9.81  ! m/s^2
    c = 0.001 ! 1/m
    v_init = 200 ! m/s
    dt = 0.001

    t = 0.
    ymax = 0.

    call trajectory(t, ymax, part, angle, dt, g, c, v_init)


!print some statistics
print *, "Velocity on impact", (part%vx * part%vx + part%vy * part%vy)**0.5
print *, "Angle on impact", atan(part%y/part%x) * 180/pi
print *, "Distance covered", part%x
print *, "Time of flight", t
print *, "Highest point reached", ymax


! !optimise angle for maximum distance
! angle = 45. * pi/180.d0
! print *, angle
! call optimize_distance(angle, 10.**-5)


end program

! !find optimal angle for biggest distance travelled
! !use steepest descent
! !could also use newton-raphson since its only 1D (only angle parameter)
! subroutine optimize_distance(angle, eps)
!     use pkg
!     integer :: i
!     real :: angle, eps, delta, disto, distn, dist2, grad

!     type(particle) :: part

!     print *, "Hello"
!     print *, angle, grad, eps
!     call trajectory(t, ymax, part, angle, dt, g, c, v_init)
!     disto = part%x


!     delta = eps + 1 !set delta to some number > eps

!     i = 0
!     do while(abs(delta) > eps .and. i < 500)
!         i = i + 1
!         !calculate the derivative:
!         call trajectory(t, ymax, part, angle+0.001, dt, g, c, v_init)
!         dist2 = part%x
!         grad = (dist2 - disto)/0.001
!         angle = angle + grad

!         print *, angle, grad

!         call trajectory(t, ymax, part, angle, dt, g, c, v_init)
!         distn = part%x

!         delta = disto - distn

        
!     end do

!     print *, angle

! end subroutine optimize_distance