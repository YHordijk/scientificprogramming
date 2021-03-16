program test

  use unit_testing
  use mcscf_routines
  implicit none

  integer :: ier
  integer :: key   
  integer :: nxy   
  integer :: ix    
  integer :: iy    
  integer :: irepij
  integer :: iprint
  integer :: nz    
  integer :: lupri 
  integer :: norbt 
  integer :: nasht 
  integer :: nnashx
  integer, allocatable, dimension(:,:) :: ipqtoq
  real*8,  allocatable, dimension(:,:,:,:) :: h2xy  
  real*8,  allocatable, dimension(:,:,:,:) :: h2yx  
  real*8,  allocatable, dimension(:,:,:,:,:) :: pv    
  real*8,  allocatable, dimension(:,:,:) :: fq    

  open(60,file='raddfq90.inp_4_3_1',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating raddfq90.inp_4_3_1 file'
    stop
  endif

  read(60,*)
  read(60,*) nxy   
  read(60,*)
  read(60,*) ix    
  read(60,*)
  read(60,*) iy    
  read(60,*)
  read(60,*) irepij
  read(60,*)
  read(60,*) iprint
  read(60,*)
  read(60,*) nz    
  read(60,*)
  read(60,*) lupri 
  read(60,*)
  read(60,*) norbt 
  read(60,*)
  read(60,*) nasht 
  read(60,*)
  read(60,*) nnashx

  allocate( ipqtoq(4,8) , stat=ier ); call errall(ier,'test','ipqtoq')
  read(60,*)
  read(60,*) ipqtoq

  allocate( h2xy  (37,4,2,3) , stat=ier ); call errall(ier,'test','h2xy  ')
  read(60,*)
  read(60,*) h2xy  

  allocate( h2yx  (37,4,2,3) , stat=ier ); call errall(ier,'test','h2yx  ')
  read(60,*)
  read(60,*) h2yx  

  allocate( pv    (4,4,10,2,3) , stat=ier ); call errall(ier,'test','pv    ')
  read(60,*)
  read(60,*) pv    

  allocate( fq    (37,4,2) , stat=ier ); call errall(ier,'test','fq    ')
  read(60,*)
  read(60,*) fq    


  close(60,status='keep')

  key    = 1
  call raddfq90( &
          key    = key    &
         ,nxy    = nxy    &
         ,ix     = ix     &
         ,iy     = iy     &
         ,irepij = irepij &
         ,iprint = iprint &
         ,nz     = nz     &
         ,lupri  = lupri  &
         ,norbt  = norbt  &
         ,nasht  = nasht  &
         ,nnashx = nnashx &
         ,ipqtoq = ipqtoq &
         ,h2xy   = h2xy   &
         ,h2yx   = h2yx   &
         ,pv     = pv     &
         ,fq     = fq     &
         )

  call compare_files(60,'raddfq90.out','raddfq90.out_4_3_1_ref',1.0d-15)

  deallocate( ipqtoq )
  deallocate( h2xy   )
  deallocate( h2yx   )
  deallocate( pv     )
  deallocate( fq     )

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

end program
