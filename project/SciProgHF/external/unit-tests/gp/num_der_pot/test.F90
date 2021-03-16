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

  ycoor  (           1) =  -0.32064837846948949D+04
  ycoor  (           2) =  -0.32064855623761514D+04
  ycoor  (           3) =  -0.32064870429125444D+04
  ycoor  (           4) =  -0.32064882120769762D+04
  ycoor  (           5) =  -0.32064890549452498D+04
  ycoor  (           6) =  -0.32064895557852524D+04
  ycoor  (           7) =  -0.32064896979442847D+04
  ycoor  (           8) =  -0.32064894638322421D+04
  ycoor  (           9) =  -0.32064888349043949D+04
  ycoor  (          10) =  -0.32064877916116138D+04
  ycoor  (          11) =  -0.32064863133671979D+04
  ycoor  (          12) =  -0.32064843782902053D+04
  ycoor  (          13) =  -0.32064819630582738D+04

  end subroutine

  subroutine get_data_yerr   (yerr   )
  implicit none
  real*8 :: yerr

  yerr    =   0.10000000000000000D-09

  end subroutine

  subroutine get_data_ierr   (ierr   )
  implicit none
  integer :: ierr

  ierr    =            0

  end subroutine

  subroutine get_data_derres_ref(derres_ref)
  implicit none
  real*8, dimension(*) :: derres_ref

  derres_ref(           1) =   0.75355837327563705D-03
  derres_ref(           2) =   0.11571074772130342D+00
  derres_ref(           3) =   0.97402710546654475D-01
  derres_ref(           4) =   0.89378514267103637D-01
  derres_ref(           5) =   0.99999999999999820D-08
  derres_ref(           6) =   0.99999999999999943D-07
  derres_ref(           7) =   0.99999999999999910D-04
  derres_ref(           8) =   0.99999999999999937D-03
  derres_ref(           9) =   0.12041293660092347D-07
  derres_ref(          10) =   0.21479983085420076D-06
  derres_ref(          11) =   0.16403644216615243D-04
  derres_ref(          12) =   0.37546854620346381D-03

  end subroutine
end program
