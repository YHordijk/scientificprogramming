program test

  use unit_testing
  use num_der
  implicit none

  integer :: ier
  integer :: nder   
  character(70) :: derstr 
  integer :: ipoint 
  integer :: nprint 
  integer :: lunit  
  integer :: npoints
  real*8,  allocatable, dimension(:) :: xcoor  
  real*8,  allocatable, dimension(:) :: ycoor  
  real*8 :: yerr   
  integer :: ierr   
  real*8,  allocatable, dimension(:,:) :: derres
  real*8,  allocatable, dimension(:,:) :: derres_ref

  allocate( xcoor  (13) , stat=ier ); call errall(ier,'test','xcoor  ')
  allocate( ycoor  (13) , stat=ier ); call errall(ier,'test','ycoor  ')
  allocate( derres(4,3) , stat=ier ); call errall(ier,'test','derres')
  allocate( derres_ref(4,3) , stat=ier ); call errall(ier,'test','derres_ref')

  call get_data_nder   (nder   )
  call get_data_derstr (derstr )
  call get_data_ipoint (ipoint )
  call get_data_nprint (nprint )
  call get_data_lunit  (lunit  )
  call get_data_npoints(npoints)
  call get_data_xcoor  (xcoor  )
  call get_data_ycoor  (ycoor  )
  call get_data_yerr   (yerr   )
  call get_data_ierr   (ierr   )

  call nddrv( &
          nder    = nder    &
         ,derstr  = derstr  &
         ,ipoint  = ipoint  &
         ,nprint  = nprint  &
         ,lunit   = lunit   &
         ,npoints = npoints &
         ,xcoor   = xcoor   &
         ,ycoor   = ycoor   &
         ,yerr    = yerr    &
         ,ierr    = ierr    &
         ,derres = derres &
         )

  call get_data_derres_ref(derres_ref)

  call compare_2vectors(12,derres_ref,derres,1.0d-14)

  deallocate( xcoor   )
  deallocate( ycoor   )
  deallocate( derres )
  deallocate( derres_ref )

  if(all_tests_passed()) write(6,*) ' test result: =====  OK  ====='

contains

  subroutine errall(ier,chsub,charr)
  implicit none

  integer      :: ier
  character(*) :: chsub, charr

  if(ier /= 0)then
    write(6,*)
    write(6,'(a,a,a,a)')' Error in allocation of ',charr,' in subroutine ',chsub
    write(6,'(a,i4)'   )' iostat =',ier
    stop 'ERROR in allocation (see output form more information)'
  endif

  end subroutine

  subroutine get_data_nder   (nder   )
  implicit none
  integer :: nder

  nder    =            4

  end subroutine

  subroutine get_data_derstr (derstr )
  implicit none
  character(*) :: derstr

  derstr  = 'CENTRAL'

  end subroutine

  subroutine get_data_ipoint (ipoint )
  implicit none
  integer :: ipoint

  ipoint  =            7

  end subroutine

  subroutine get_data_nprint (nprint )
  implicit none
  integer :: nprint

  nprint  =            1

  end subroutine

  subroutine get_data_lunit  (lunit  )
  implicit none
  integer :: lunit

  lunit   =            3

  end subroutine

  subroutine get_data_npoints(npoints)
  implicit none
  integer :: npoints

  npoints =           13

  end subroutine

  subroutine get_data_xcoor  (xcoor  )
  implicit none
  real*8, dimension(*) :: xcoor

  xcoor  (           1) =  -0.34211278586564475D+00
  xcoor  (           2) =  -0.28509398822137061D+00
  xcoor  (           3) =  -0.22807519057709652D+00
  xcoor  (           4) =  -0.17105639293282238D+00
  xcoor  (           5) =  -0.11403759528854826D+00
  xcoor  (           6) =  -0.57018797644274130D-01
  xcoor  (           7) =   0.00000000000000000D+00
  xcoor  (           8) =   0.57018797644274130D-01
  xcoor  (           9) =   0.11403759528854826D+00
  xcoor  (          10) =   0.17105639293282238D+00
  xcoor  (          11) =   0.22807519057709652D+00
  xcoor  (          12) =   0.28509398822137061D+00
  xcoor  (          13) =   0.34211278586564475D+00

  end subroutine

  subroutine get_data_ycoor  (ycoor  )
  implicit none
  real*8, dimension(*) :: ycoor

  ycoor  (           1) =  -0.16579876115919999D-16
  ycoor  (           2) =  -0.13907340588530001D-16
  ycoor  (           3) =  -0.11299810924709999D-16
  ycoor  (           4) =  -0.87559347056730007D-17
  ycoor  (           5) =  -0.62753462779890000D-17
  ycoor  (           6) =  -0.38579395472780003D-17
  ycoor  (           7) =  -0.15039675167090001D-17
  ycoor  (           8) =   0.78585261928509998D-18
  ycoor  (           9) =   0.30102693790870000D-17
  ycoor  (          10) =   0.51679229235660000D-17
  ycoor  (          11) =   0.72565503446160005D-17
  ycoor  (          12) =   0.92733414683429999D-17
  ycoor  (          13) =   0.11215232135030001D-16

  end subroutine

  subroutine get_data_yerr   (yerr   )
  implicit none
  real*8 :: yerr

  yerr    =   0.00000000000000000D+00

  end subroutine

  subroutine get_data_ierr   (ierr   )
  implicit none
  integer :: ierr

  ierr    =            0

  end subroutine

  subroutine get_data_derres_ref(derres_ref)
  implicit none
  real*8, dimension(*) :: derres_ref

  derres_ref(           1) =   0.40724732716668484D-16
  derres_ref(           2) =  -0.19714233902489241D-16
  derres_ref(           3) =  -0.58815580857970344D-17
  derres_ref(           4) =  -0.68926121544162399D-16
  derres_ref(           5) =   0.99999999999999699D-22
  derres_ref(           6) =   0.99999999999999464D-21
  derres_ref(           7) =   0.99999999999999429D-18
  derres_ref(           8) =   0.99999999999999899D-17
  derres_ref(           9) =   0.00000000000000000D+00
  derres_ref(          10) =   0.00000000000000000D+00
  derres_ref(          11) =   0.00000000000000000D+00
  derres_ref(          12) =   0.00000000000000000D+00

  end subroutine
end program
