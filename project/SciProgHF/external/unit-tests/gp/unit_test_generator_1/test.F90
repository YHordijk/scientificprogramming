program test

  use unit_testing
  use unit_test_generator
  implicit none

  integer :: ier
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

  allocate( ipqtoq(4,8) , stat=ier ); call errall(ier,'test','ipqtoq')
  allocate( h2xy  (37,4,2,3) , stat=ier ); call errall(ier,'test','h2xy  ')
  allocate( h2yx  (37,4,2,3) , stat=ier ); call errall(ier,'test','h2yx  ')
  allocate( pv    (4,4,10,2,3) , stat=ier ); call errall(ier,'test','pv    ')
  allocate( fq    (37,4,2) , stat=ier ); call errall(ier,'test','fq    ')

  open(60,file='INPUT',status='old',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating INPUT file'
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

  read(60,*)
  read(60,*) ipqtoq

  read(60,*)
  read(60,*) h2xy  

  read(60,*)
  read(60,*) h2yx  

  read(60,*)
  read(60,*) pv    

  read(60,*)
  read(60,*) fq    


  close(60,status='keep')

  !------------------------------------------------------------------------------!
  !             Generate test with input data writen to the file                 !
  !------------------------------------------------------------------------------!
  ! if raddfq90.inp_4_3_1 already exist delete it
  open(60,file='raddfq90.inp_4_3_1',status='unknown',form='formatted',access='sequential',iostat=ier)
  if(ier /= 0)then
    write(6,*) 'test: ERROR while creating raddfq90.inp_4_3_1 file'
    stop
  endif
  close(60,status='delete')

  call unit_test__constructor
  call unit_test__set_suffix(nxy)
  call unit_test__set_suffix(ix)
  call unit_test__set_suffix(iy)
  call unit_test__add_input_file ('raddfq90.inp')
  call unit_test__set_routine_name('raddfq90')
  call unit_test__add_module('mcscf_routines')

  call unit_test__add_opt_variable('key   ','1','integer')

  call unit_test__add_inp_variable(nxy   ,'nxy   ')
  call unit_test__add_inp_variable(ix    ,'ix    ')
  call unit_test__add_inp_variable(iy    ,'iy    ')
  call unit_test__add_inp_variable(irepij,'irepij')
  call unit_test__add_inp_variable(iprint,'iprint')
  call unit_test__add_inp_variable(nz    ,'nz    ')
  call unit_test__add_inp_variable(lupri ,'lupri ')
  call unit_test__add_inp_variable(norbt ,'norbt ')
  call unit_test__add_inp_variable(nasht ,'nasht ')
  call unit_test__add_inp_variable(nnashx,'nnashx')
  call unit_test__add_inp_variable(ipqtoq,'ipqtoq')
  call unit_test__add_inp_variable(h2xy  ,'h2xy  ')
  call unit_test__add_inp_variable(h2yx  ,'h2yx  ')
  call unit_test__add_inp_variable(pv    ,'pv    ')

  call unit_test__add_inp_variable(fq    ,'fq    ')
  call unit_test__add_out_variable(fq    ,'fq    ',compare='file-file',file='raddfq90.out')

  call unit_test__generate_test()

  ! compare input file and test
  call compare_text_files(60, 'raddfq90.inp_4_3_1', 'raddfq90.inp_ref')
  call compare_text_files(60, 'test.F90_4_3_1', 'test.F90_4_3_1_ref')

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
