module pkg
   real*8 :: pi = 4.0*datan(1.d0)

   type :: particle
      real*8 :: x, y, vx, vy
   end type particle

   public :: pi , particle, traj_result, result, integrate, trajectory
   private


contains

!subroutine that does one integration cycle
subroutine integrate(p, dt, g, c)
   implicit none
   real*8, intent(in) :: dt, g, c
   type(particle) ::  p
   
   p%x = p%x + p%vx * dt
   p%y = p%y + p%vy * dt
   p%vx = p%vx - c * p%vx * p%vx * dt
   if (p%vy > 0) then
         p%vy = p%vy - g*dt - c * p%vy * p%vy * dt
   else
         p%vy = p%vy - g*dt + c * p%vy * p%vy * dt
   end if

end subroutine integrate


subroutine trajectory(t, ymax, part, angle, dt, g, c, v_init)
   real*8, intent(in) :: g, angle, c, v_init, dt
   real*8, intent(inout) :: t, ymax

   type(particle), intent(out) :: part
 
   !initialize particle type
   part%x = 0
   part%y = 0
   part%vx = v_init * cos(angle)
   part%vy = v_init * sin(angle)

   !do the actual calculation here, continue integrating untill we hit y=0 again
   do while(part%y >= 0)
      !integrate the movement of the particle
      call integrate(part, dt, g, c)
      ymax = max(ymax, part%y)
      t = t + dt
   end do 

end subroutine trajectory

end module pkg