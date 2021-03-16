  module contraction
  use wrapper_sorting
  use symmetry_offset
  use modified_sorting
    implicit none

    public contraction_222
    public contraction_442   
    public contraction_242
    public contraction_422
    public contraction_444
    public contraction_424
    public contraction_244

     private

     contains

   subroutine contraction_222 (left_order,right_order,im_order,leftTensor,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,irrep_left,irrep_right)

!----------calling variables-----------------------------------
    character(2), dimension(2), intent(in)   ::left_order,im_order,right_order
    real(8), intent(in)                      :: rightTensor(:),leftTensor(:)   
    integer, intent(in)                      :: nrep
    real(8), intent(in)                      :: dia_factor, dia_sum_factor             !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)                   :: imTensor(:)
    integer,intent(in),optional              :: irrep_left,irrep_right
!--------------------------------------------------------------
    integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
    character(1)   :: transA,transB
    character(2)   :: freeindex_left(2)
    character(2)   :: freeindex_right(2)
    logical        :: array_equal = .true. 
    real*8,allocatable :: ImTensor_temp(:) 
    real*8,allocatable :: ImTensor_temp_temp(:) 
    integer             :: irrep_left_loc,irrep_right_loc
    integer             :: irrep_target

#include "param.inc"
#include "complex.inc"

    if (.not.present(irrep_left)) then
       irrep_left_loc = 1
    else
      irrep_left_loc = irrep_left
    endif 

    if (.not.present(irrep_right)) then
      irrep_right_loc = 1
    else
      irrep_right_loc = irrep_right
    endif 

    transA = 'N'  
    transB = 'N'  

!  1) fetch the fixed arrays, properly sort to contract with BLAS.
   
    call onebody_for_blas (left_order,im_order,right_order,transA,x_row,x_column,freeindex_left)

    call onebody_for_blas (right_order,left_order,im_order,transB,y_row,y_column,freeindex_right)
   
    array_equal = all((/freeindex_left(1),freeindex_right(2)/)==im_order)

    if (array_equal) then

       call cntrct_new (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor, &
   &                   dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.false.) 

   else

      allocate (ImTensor_temp(dot_product(x_row,y_column)*rcw))

      ImTensor_temp = 0.0d0

      call cntrct_new (transA,transB,x_row,y_column,x_column,dia_factor*A1,leftTensor,rightTensor, &
   &                   dia_sum_factor*A1,imTensor_temp,irrep_left_loc,irrep_right_loc,irrep_target,.false.) 

      allocate (ImTensor_temp_temp(dot_product(x_row,y_column)*rcw))
      ImTensor_temp_temp = 0.0d0

      call SRT1C1_exp (irrep_target,.false.,x_row,y_column,imTensor_temp,imTensor_temp_temp,2)

!     call srt1c1n (nrep,x_row,y_column,imTensor_temp,imTensor_temp_temp)

      call daxpy (size(imTensor),dia_sum_factor*A1,imTensor_temp_temp,1,imTensor,1) 

      deallocate (ImTensor_temp)
      deallocate (ImTensor_temp_temp)

   endif

   end subroutine

   subroutine contraction_442 (left_order,right_order,im_order,ImTensor, &
           &   dia_factor,dia_sum_factor,nrep,LeftTensor,RightTensor,irrep_left,irrep_right)
     
!----------calling variables-----------------------------------
    real(8), intent(in) ::  dia_sum_factor   !dia_factor = \alpha in DGEMM ; dia_sum_factor= \beta in DGEMM
    real(8), intent(in) ::  dia_factor        
    character(2),dimension(4), intent(in) :: left_order, right_order
    character(2),dimension(2), intent(in) :: im_order
    real(8),optional,intent(in)           :: RightTensor(:),LeftTensor(:)   
    integer,intent(in)    :: nrep
    integer,intent(in),optional  :: irrep_left,irrep_right
    integer   :: x_row(nrep),x_column(nrep)
    integer   :: y_column(nrep),y_row(nrep)              
    real(8),intent(inout) :: ImTensor(:)
    integer               :: match_index_R,match_index_L  
    character(1)          :: transA,transB
    real(8),allocatable   :: v_sort(:)
    real(8),allocatable   :: right_sorted(:)    
    character(2)          :: sorted_oparray_left(4)
    character(2)          :: sorted_oparray_right(4)
    logical               :: array_equal = .true. 
    integer               :: asym_factor_R,asym_factor_L,i 
    integer               :: maxdim1, maxdim2 
    integer               :: irrep_left_loc,irrep_right_loc,irrep_target
    real(8),allocatable   :: ImTensor_temp(:) 
    real(8),allocatable   :: ImTensor_temp_temp(:) 
    integer :: im_temp_size,  im_half_size_x,  im_half_size_y
!--------------------------------------------------------------

#include "param.inc"
#include "complex.inc"

      maxdim1 = sorted_dimension(left_order)
      maxdim2 = sorted_dimension(right_order)

 if (.not.present(irrep_left)) then
   irrep_left_loc = 1
 else
   irrep_left_loc = irrep_left
 endif 

 if (.not.present(irrep_right)) then
   irrep_right_loc = 1
 else
   irrep_right_loc = irrep_right
 endif 

  allocate(v_sort(maxdim1*rcw))  
  allocate(right_sorted(maxdim2*rcw))  
 
      transA = 'N'
      transB = 'N'

! 1) fetch 2-e integrals and sort them for BLAS.

     call decide_sorting_routine(13,left_order,right_order,im_order,transA,match_index_L,asym_factor_L, &
    &   .true.,sorted_opArray=sorted_oparray_left)

!     write(*,*)'debug info left array 13 ',transA,match_index_L,asym_factor_L

     if (present(LeftTensor)) then 
      call sorting_for_blas(13,left_order,match_index_L,.false.,transA,v_sort,x_row,x_column,LeftTensor, &
            irrep_input=irrep_left_loc)
     else
      call sorting_for_blas(13,left_order,match_index_L,.false.,transA,v_sort,x_row,x_column, &
            irrep_input=irrep_left_loc)
     endif 

! 2) general sorting routine for fixed arrays. 

     call decide_sorting_routine(31,right_order,im_order,sorted_oparray_left,transB,match_index_R,asym_factor_R, &
    &         .false.,sorted_opArray=sorted_oparray_right)

!     write(*,*)'debug info right array 31 ',transB,match_index_R,asym_factor_R,&
!    & sorted_oparray_right

     array_equal = all((/sorted_oparray_left(1),sorted_oparray_right(4)/)==im_order)

     if (present(RightTensor)) then 
     call sorting_for_blas(31,right_order,match_index_R,.false.,transB,right_sorted,y_row,y_column,rightTensor, &
           irrep_input=irrep_right_loc)
     else
     call sorting_for_blas(31,right_order,match_index_R,.false.,transB,right_sorted,y_row,y_column, &
          irrep_input=irrep_right_loc)
     endif  

! 3) contraction routine for symmetry-packed arrays.

   if (array_equal) then
      call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort, &
   &                right_sorted,dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.false.) 

    else

   im_temp_size = dot_product(x_row,y_column)
!
! aspg: introduced below an additional "cross"-dot product to get the dimension of the work variable below.
!
!       just doing a dot product of x_row and y_column to get the total dimesion of the variable
!       seems fine as far as i can see for groups without inversion symmetry. 
!
!       however, in the case of groups with inversion the dot product cannot allow for  calculating the dimensions of kind
!
!       x_row(irrep_g) * y_column(irrep_u) (and vice-versa) 
!
!       that may occur when one targets a "u" irrep. if the value for the "u" symmetry is smaller than the
!       value coming from the application of the dot product there will be no consequences, otherwise the xgemm 
!       call inside cntrct_new can show a segmentation fault. 
!
!       below, the final work variable dimension is then set to the maximum between the pure dot product
!       and the "cross"-dot product
!
   im_half_size_x = size(x_row,1)/2     
   im_half_size_y = size(y_column,1)/2     
   im_temp_size = max(im_temp_size, &
  &                   (dot_product(x_row(im_half_size_x+1:im_half_size_x*2),y_column(1:im_half_size_y)) &
  &                   +dot_product(x_row(1:im_half_size_x),y_column(im_half_size_y+1:im_half_size_y*2))))

   allocate(ImTensor_temp(im_temp_size*rcw))

      ImTensor_temp = 0.0d0

        call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort, &
    &                 right_sorted,A1,imTensor_temp,irrep_left_loc,irrep_right_loc,irrep_target,.false.) 

      allocate (ImTensor_temp_temp(dot_product(x_row,y_column)*rcw))

      call srt1c1_exp (irrep_target,.false.,x_row,y_column,imTensor_temp,imTensor_temp_temp,2)
!     call srt1c1n (nrep,x_row,y_column,imTensor_temp,imTensor_temp_temp)

      deallocate (ImTensor_temp)

      call daxpy (size(imTensor),dia_sum_factor*A1,imTensor_temp_temp,1,imTensor,1) 

      deallocate (ImTensor_temp_temp)

   endif

   deallocate(v_sort)  
   deallocate(right_sorted)  

   end subroutine

   subroutine contraction_422 (left_order,right_order,im_order,rightTensor,imTensor, &
           &                   dia_factor,dia_sum_factor,nrep,LeftTensor,irrep_left,irrep_right)
     
!----------calling variables-----------------------------------
    character(2),dimension(4),intent(in)   :: left_order
    character(2),dimension(2),intent(in)   :: im_order, right_order
    real(8),intent(in),optional :: LeftTensor(:)                
    real(8),intent(in)     :: rightTensor(:)   
    integer,intent(in)     :: nrep
    real(8),intent(in)     :: dia_sum_factor,dia_factor   !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(inout)  :: imTensor(:)
    integer,intent(in),optional  :: irrep_left,irrep_right
!--------------------------------------------------------------
    integer        :: x_row(nrep), x_column(nrep),y_row(nrep),y_column(nrep)
    character(1)   :: transA,transB
    integer        :: match_index,i 
    integer        :: asym_factor
    integer        :: maxdim, contraction_type
    logical        :: array_equal, totsym_local
    integer        :: irrep_left_loc,irrep_right_loc, irrep_target, nrp(nrep)
    real(8),allocatable        :: v_sort(:)
    real(8),allocatable        :: RightTensor_sort(:)
    character(2)   :: sorted_left(4),freeindx(2)
#include "param.inc"
#include "complex.inc"

!-------------------------------------------------------------
!    use of xgemv is different from dgemm because when I mention 'T', then the dimensionality of the incoming array is also changed.
!-------------------------------------------------------------

!--------------------------------------------------------------------------
!  1> fetch the array which is flexible in terms of indices. sort it to make suitable for BLAS.
!     also determine the row and column index for the matrix-vector multiplication routine.

  maxdim = sorted_dimension(left_order)


 if (.not.present(irrep_left)) then
   irrep_left_loc = 1
 else
   irrep_left_loc = irrep_left
 endif 

 if (.not.present(irrep_right)) then
   irrep_right_loc = 1
 else
   irrep_right_loc = irrep_right
 endif 

   if (irrep_right_loc == 1) then
    totsym_local = .false.
   else
    totsym_local = .false.
   endif

  allocate(v_sort(maxdim*rcw))  

     transA = 'N'

     call decide_sorting_routine(22,left_order,right_order,im_order,transA,match_index,asym_factor,.true., &
    &                        sorted_opArray=sorted_left)

     if (present(LeftTensor)) then

     call sorting_for_blas(22,left_order,match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor, &
                     if_totsym=totsym_local,irrep_input=irrep_left_loc)

     else

     call sorting_for_blas(22,left_order,match_index,.false.,transA,v_sort,x_row,x_column,if_totsym=totsym_local, &
          irrep_input=irrep_left_loc)

     endif

     call onebody_for_blas (right_order,sorted_left(3:4),im_order,transB,y_row,y_column) 

     array_equal = all((/sorted_left(1),sorted_left(2)/)==im_order)

   if ((match_index==2).or.(match_index==3)) then
    contraction_type = 2
   else
    contraction_type = 1
   endif

!      write(*,*)'inside 422: array_equal' ,array_equal 

     nrp = 0
     nrp(1) = 1


     if (transB=='T') then

      allocate(rightTensor_sort(x_column(irrep_right_loc)*rcw))

!      call srt1c1n (nrep,y_column,y_row,rightTensor,rightTensor_sort)
       call srt1c1_exp (irrep_right_loc,.true.,y_column,y_row,rightTensor,rightTensor_sort,1)

!!      if (carith) call conjuga(x_column(1),rightTensor,1)

   select case(contraction_type)

    case(1)

      call cntrct_4 (transa,'N',x_row,nrp,x_column,asym_factor*dia_factor*A1,v_sort,rightTensor_sort, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

     case(2)
      call cntrct_5 (transa,'N',x_row,1,x_column,asym_factor*dia_factor*A1,v_sort,rightTensor_sort, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

  end select 


!     call xgemm (transA,'N',x_row(irrep_right_loc),1,x_column(irrep_right_loc),asym_factor*dia_factor*A1, &
!            v_sort(1+x_row(irrep_right_loc-1)*x_column(irrep_right_loc-1)), &
!    &       x_row(irrep_right_loc),rightTensor_sort,x_column(irrep_right_loc),dia_sum_factor*A1,imTensor,x_row(irrep_right_loc))

!     call xgemv (transA,x_row(1),x_column(1),asym_factor*dia_factor*A1,v_sort,x_row(1),rightTensor_sort, &
!     &            1,dia_sum_factor*A1,imTensor,1) 

      deallocate(rightTensor_sort)
     
      else

   select case(contraction_type)

    case(1)

      call cntrct_4 (transa,'N',x_row,nrp,x_column,asym_factor*dia_factor*A1,v_sort,rightTensor, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

     case(2)
      call cntrct_5 (transa,'N',x_row,1,x_column,asym_factor*dia_factor*A1,v_sort,rightTensor, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

  end select 


!     call xgemm (transA,'N',x_row(irrep_target),1,x_column(irrep_right_loc),asym_factor*dia_factor*A1,v_sort, &
!     &       lda,rightTensor,x_column(irrep_right_loc),dia_sum_factor*A1,imTensor,x_row(irrep_target))

!!       call xgemv (transA,x_row(1),x_column(1),asym_factor*dia_factor*A1,v_sort,x_row(1),rightTensor,1,dia_sum_factor*A1,imTensor,1) 
      endif

      deallocate(v_sort)
    end subroutine

   subroutine contraction_242 (left_order,right_order,im_order,LeftTensor,ImTensor, &
           &                   dia_factor,dia_sum_factor,nrep,RightTensor,irrep_left,irrep_right)


!----------calling variables-----------------------------------
    character(2), dimension(2), intent(in)   :: left_order,im_order
    character(2), dimension(4), intent(in)   :: right_order
    real(8), intent(in)                      :: LeftTensor(:)
    real(8), intent(in), optional            :: RightTensor(:)   
    integer, intent(in)                      :: nrep
    real(8), intent(in)                      :: dia_sum_factor        
    real(8), intent(in)                      :: dia_factor            
    real(8), intent(inout)                   :: ImTensor(:)
    integer,intent(in),optional  :: irrep_left,irrep_right
!--------------------------------------------------------------
    character(2)   :: freeindex_left(2)
    character(2)   :: sorted_oparray_right(4) 
    integer        :: y_row(nrep),y_column(nrep)
    integer        :: x_row(nrep),x_column(nrep)
    character(1)   :: transA,transB,trans1A
    integer        :: match_index,i 
    integer        :: asym_factor
    integer        :: maxdim
    integer        :: irrep_left_loc,irrep_right_loc,irrep_target 
    integer        :: ldb,contraction_type 
    logical        :: array_equal,totsym_local
    real(8),allocatable        :: right_sort(:)
    real(8),allocatable        :: left_sort(:)
!---------------------------------------------------------------    
#include "param.inc"
#include "complex.inc"

    maxdim = sorted_dimension(right_order)

 if (.not.present(irrep_left)) then
   irrep_left_loc = 1
 else
   irrep_left_loc = irrep_left
 endif 

 if (.not.present(irrep_right)) then
   irrep_right_loc = 1
 else
   irrep_right_loc = irrep_right
 endif 

   if (irrep_left_loc == 1) then
    totsym_local = .false.
   else
    totsym_local = .false.
   endif

   allocate(right_sort(maxdim*rcw))  

     transA  =  'N'
     trans1A = 'N'

!** 1> fetch the integral array, properly sort to contract with BLAS.

   call decide_sorting_routine(22,right_order,im_order,left_order,transA,match_index,asym_factor,.true.,&
   &                          sorted_opArray=sorted_oparray_right)


   if (present(RightTensor)) then
  call sorting_for_blas(22,right_order,match_index,.false.,transA,right_sort,y_row,y_column,presorted=RightTensor, &
         if_totsym=totsym_local,irrep_input=irrep_right_loc)
   else
  call sorting_for_blas(22,right_order,match_index,.false.,transA,right_sort,y_row,y_column,if_totsym=totsym_local, &
       irrep_input=irrep_right_loc)
   endif

!**  2> sort onebody integral properly to use in BLAS.

   call onebody_for_blas (left_order,im_order,sorted_oparray_right(1:2),trans1A,x_row,x_column)

     array_equal = all((/sorted_oparray_right(3),sorted_oparray_right(4)/)==im_order)

!        if (transA.eq.'N'.or.transA.eq.'n') then
!           ldb = y_row(irrep_right_loc)
!        else
!           ldb = y_column(irrep_target)
!        endif


   if ((match_index==2).or.(match_index==3)) then
    contraction_type = 2
   else
    contraction_type = 1
   endif

!**  3>  to do B = XA  type operation with xgemv, we are doing B^T = A^T*X^T then transpose(B^T). that implies following modifications.   

   if (trans1A=='T') then

   allocate (left_sort(y_row(irrep_left_loc)*rcw))

!   call srt1c1n (nrep,x_column,x_row,LeftTensor,left_sort)
   call srt1c1_exp (irrep_left_loc,.false.,x_column,x_row,LeftTensor,left_sort,1)

   select case(contraction_type)

    case(1)

      call cntrct_6 ('N',transA,1,y_column,y_row,asym_factor*dia_factor*A1,left_sort,right_sort, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

     case(2)

      call cntrct_7 ('N',transA,1,y_column,y_row,asym_factor*dia_factor*A1,left_sort,right_sort, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

  end select 

!  call xgemm ('N',transA,1,y_column(irrep_left_loc),y_row(irrep_left_loc),asym_factor*dia_factor*A1,left_sort, &
!     &       1,right_sort,ldb,dia_sum_factor*A1,imTensor,1)

   deallocate (left_sort)

   else
  
 select case(contraction_type)

    case(1)

      call cntrct_6 ('N',transA,1,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,right_sort, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

     case(2)

      call cntrct_7 ('N',transA,1,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,right_sort, &
     &            dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.)

  end select 
  
!  call xgemm ('N',transA,1,y_column(irrep_left_loc),y_row(irrep_left_loc),asym_factor*dia_factor*A1,LeftTensor, &
!     &       1,right_sort,ldb,dia_sum_factor*A1,imTensor,1)

   endif
   deallocate (right_sort)

    end subroutine

   subroutine contraction_444 (left_order,right_order,im_order,imTensor, &
           &    dia_factor,dia_sum_factor,nrep,RightTensor,LeftTensor,irrep_left,irrep_right)
    
!----------calling variables-----------------------------------
    character(*), dimension(4), intent(in)   :: left_order, right_order, im_order
    real(8), intent(in), optional :: LeftTensor(:)       
    real(8), intent(in), optional :: RightTensor(:)   
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor,dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout)  :: imTensor(:)
    integer,intent(in),optional  :: irrep_left,irrep_right
!--------------------------------------------------------------
   integer        :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
   character(1)   :: transA, transB
   integer        :: save_index_l,save_index_r,save_index,i
   character(2)   :: sorted_opArray_left(4),sorted_opArray_right(4)
   real(8),allocatable :: v_sort(:)
   real(8),allocatable :: right_sorted(:)
   real(8),allocatable :: im_temp(:)
   real(8),allocatable :: im_temp_temp(:)
   integer             :: asym_factor_R,asym_factor_L
   integer             :: irrep_left_loc,irrep_right_loc
   integer             :: irrep_target
   integer             :: maxdim1, maxdim2
   integer             :: contraction_type

#include "param.inc"
#include "complex.inc"
#include "dcbham.h"

  maxdim1 = max(sorted_dimension(left_order),sorted_dimension(im_order))
  maxdim2 = sorted_dimension(right_order)

 if (.not.present(irrep_left)) then
   irrep_left_loc = 1
 else
   irrep_left_loc = irrep_left
 endif 

 if (.not.present(irrep_right)) then
   irrep_right_loc = 1
 else
   irrep_right_loc = irrep_right
 endif 

  allocate(v_sort(maxdim1*rcw))  
  allocate(right_sorted(maxdim2*rcw))  

   transA = 'N'
   transB = 'N'

!  1) fetch/provide the left tensor and sort properly for BLAS. also extract row and column dimensions of the sorted tensor.

    call decide_sorting_routine(22,left_order,right_order,im_order,transA,save_index_l,asym_factor_L,.true., &
    &  sorted_opArray=sorted_opArray_left)

!    write(*,*)'debug info left array 22:  ',transA,save_index_L,asym_factor_L,sorted_opArray_left

    if (present(LeftTensor)) then 
    call sorting_for_blas(22,left_order,save_index_l,.false.,transA,v_sort,x_row,x_column, presorted=LeftTensor, &
        if_TotSym=.false.,irrep_input=irrep_left_loc)
    else
    call sorting_for_blas(22,left_order,save_index_l,.false.,transA,v_sort,x_row,x_column,if_TotSym=.false., &
          irrep_input=irrep_left_loc)
    endif

!  2) fetch/provide the right tensor and sort properly for BLAS. also extract row and column dimensions of the sorted tensor.

    call decide_sorting_routine(22,right_order,im_order,sorted_opArray_left,transB,save_index_r,asym_factor_R,&
    &          .false.,sorted_opArray=sorted_opArray_right)


!     write(*,*)'debug info right array 22:  ',transB,save_index_r,asym_factor_r,sorted_opArray_right

    if (present(RightTensor)) then 
    call sorting_for_blas(22,right_order,save_index_r,.false.,transB,right_sorted,y_row,y_column,&
   &                  presorted=RightTensor,if_TotSym=.false.,irrep_input=irrep_right_loc)
    else
    call sorting_for_blas(22,right_order,save_index_r,.false.,transB,right_sorted,y_row,y_column, &
        if_TotSym=.false.,irrep_input=irrep_right_loc)
    endif

!  3) Distinguish the contraction cases based on symmetry ordering of the sorted tensors.
!     This disntiction is required since some of them use multb(*,*,1) for final
!     packing and some others multb(*,*,2). As of now I am hard-coding them, until
!     I get a better solution.

    if (((save_index_l==2).or.(save_index_l==3)).and.(save_index_r==1)) then
      contraction_type = 2 
    elseif ((save_index_l==1).and.((save_index_r==2).or.(save_index_r==3))) then
      contraction_type = 3
    elseif (((save_index_l==2).or.(save_index_l==3)).and.((save_index_r==2).or.(save_index_r==3))) then
      contraction_type = 4
    else
      contraction_type = 1
    endif  

!  4) matrix-matrix multiplication routine runs for each irreps.

   if (dia_sum_factor==0.0d0) then

    allocate(im_temp(dot_product(x_row,y_column)*rcw))

    im_temp = 0.0d0

   select case(contraction_type)

    case(1)

    call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort,&
    &    right_sorted,A0,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.true.) 

    case(2)

    call cntrct_1 (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort, &
    &      right_sorted,A0,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.true.) 

    case(3)

    call cntrct_2 (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort, &
    &              right_sorted,A0,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.true.) 
    case(4)

    call cntrct_3 (transA,transB,x_row,y_column,x_column,asym_factor_L*asym_factor_R*dia_factor*A1,v_sort, &
    &           right_sorted,A0,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.true.) 

  end select  

!  5) resort and antisymmetrize the product tensor. sum that up with the particular class of intermediate/residue/sigma etc. 

   transA = 'N'

   v_sort = 0.0d0

   call decide_sorting_routine(22,im_order,left_order,(/sorted_opArray_left(1:2),sorted_opArray_right(3:4)/), &
                       & transA,save_index,asym_factor_R,.false.,sorted_opArray=sorted_opArray_right,inverse=.true.)

! write(*,*)'debug info intermediate array 22:  ' ,transA,save_index,asym_factor_r,sorted_opArray_right

   if (transA=='T') then

!    write(*,*)'possibility of mistake',x_row(1),y_column(1),x_row(2),y_column(2)

   allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))

   im_temp_temp = 0.0d0

!   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)
   call srt1c1_exp (irrep_target,.true.,x_row,y_column,im_temp,im_temp_temp,1)

   call sorting_for_blas(22,im_order,save_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false., &
      irrep_input=irrep_target)

     deallocate(im_temp_temp)
   else

   call sorting_for_blas(22,im_order,save_index,.true.,transA,v_sort,x_row,x_column,im_temp,if_TotSym=.false., &
             irrep_input=irrep_target)

   endif

   deallocate(im_temp)
   call daxpy (size(imTensor),A1*asym_factor_R,v_sort,1,imTensor,1)

   else
    call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor_R*asym_factor_L*dia_factor*A1,v_sort,right_sorted, &
   &              dia_sum_factor*A1,imTensor,irrep_left_loc,irrep_right_loc,irrep_target,.true.) 
   endif

   deallocate(v_sort)
   deallocate(right_sorted)

  end subroutine

   subroutine contraction_424 (left_order,right_order,im_order,ImTensor,dia_factor, &
           &     dia_sum_factor,nrep,LeftTensor,RightTensor,irrep_left,irrep_right)

!----------calling variables-----------------------------------
    character(2), dimension(4), intent(in)   :: left_order,im_order
    character(2), dimension(2), intent(in)   :: right_order
    real(8), intent(in), optional :: LeftTensor(:), RightTensor(:) 
    integer, intent(in) :: nrep
    real(8), intent(in) :: dia_factor,dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8), intent(inout) :: ImTensor(:)
    integer,intent(in),optional  :: irrep_left,irrep_right
    integer      :: x_row(nrep),x_column(nrep),y_row(nrep),y_column(nrep)
    character(1) :: transA, transB
    integer      :: match_index,i
    integer             :: irrep_left_loc,irrep_right_loc
    integer             :: irrep_target
    character(2) :: FreeIndx(2),sorted_opArray(4)
    real(8),allocatable :: v_sort (:)
    integer :: asym_factor 
    integer :: maxdim 
    real(8),allocatable :: im_temp(:)
    real(8),allocatable :: im_temp_temp(:)
    integer :: im_temp_size, im_half_size_x, im_half_size_y

#include "param.inc"
#include "complex.inc"

  maxdim = max (sorted_dimension(left_order),sorted_dimension(im_order))

 if (.not.present(irrep_left)) then
   irrep_left_loc = 1
 else
   irrep_left_loc = irrep_left
 endif 

 if (.not.present(irrep_right)) then
   irrep_right_loc = 1
 else
   irrep_right_loc = irrep_right
 endif 

  allocate(v_sort(maxdim*rcw))  

   transA = 'N'
   transB = 'N'

!  1) fetch/provide the left tensor and sort properly for BLAS. also extract row and column dimensions of the sorted tensor.

   call decide_sorting_routine(31,left_order,right_order,im_order,transA,match_index,asym_factor, &
   &             .true.,sorted_opArray=sorted_opArray)

 if (present(LeftTensor)) then

   call sorting_for_blas(31,left_order,match_index,.false.,transA,v_sort,x_row,x_column,presorted=LeftTensor, &
         irrep_input=irrep_left_loc)
 else

   call sorting_for_blas(31,left_order,match_index,.false.,transA,v_sort,x_row,x_column, &
          irrep_input=irrep_left_loc)
 endif

   call onebody_for_blas (right_order,sorted_opArray,im_order,transB,y_row,y_column,FreeIndx)

!  2) matrix-matrix multiplication routine runs for each irreps.

   im_temp_size = dot_product(x_row,y_column)
!
! aspg: introduced below an additional "cross"-dot product to get the dimension of the work variable below.
!
!       just doing a dot product of x_row and y_column to get the total dimesion of the variable
!       seems fine as far as i can see for groups without inversion symmetry. 
!
!       however, in the case of groups with inversion the dot product cannot allow for  calculating the dimensions of kind
!
!       x_row(irrep_g) * y_column(irrep_u) (and vice-versa) 
!
!       that may occur when one targets a "u" irrep. if the value for the "u" symmetry is smaller than the
!       value coming from the application of the dot product there will be no consequences, otherwise the xgemm 
!       call inside cntrct_new can show a segmentation fault. 
!
!       below, the final work variable dimension is then set to the maximum between the pure dot product
!       and the "cross"-dot product
!
   im_half_size_x = size(x_row,1)/2 
   im_half_size_y = size(y_column,1)/2 
   im_temp_size = max(im_temp_size, &
  & (dot_product(x_row(im_half_size_x+1:im_half_size_x*2),y_column(1:im_half_size_y)) &
  & +dot_product(x_row(1:im_half_size_x),y_column(im_half_size_y+1:im_half_size_y*2))))

   allocate(im_temp(im_temp_size*rcw))
   im_temp = 0.0d0

    if (match_index == 0) then
      if (present(RightTensor)) call cntrct_2 (transA,transB,x_row,y_column,x_column,asym_factor*dia_factor*A1,v_sort, &
    &                                      RightTensor,A1,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.false.)

    else

    if (present(RightTensor)) call cntrct_new (transA,transB,x_row,y_column,x_column,asym_factor*dia_factor*A1,v_sort, &
    &                                      RightTensor,A1,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.false.)
    endif


!  3) re-sort the arrays to align with the desired structure of the intermediate.

   transA = 'N'

   call decide_sorting_routine(31,im_order,left_order,(/sorted_opArray(1:3),FreeIndx(2)/),transA,match_index, &
           &                    asym_factor,.false.,inverse=.true.,sorted_opArray=sorted_opArray)

   if (transA=='T') then

!    write(*,*)'possibility of mistake',x_row(1),y_column(1),x_row(2),y_column(2)

   allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))

   im_temp_temp = 0.0d0

  ! write(*,*)'possibility of mistake'
  ! call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)
   call srt1c1_exp (irrep_target,.false.,x_row,y_column,im_temp,im_temp_temp,2)

   deallocate(im_temp)

   call sorting_for_blas(31,im_order,match_index,.true.,transA,v_sort,x_row,x_column,im_temp_temp,if_TotSym=.false., &
           irrep_input=irrep_target)

   deallocate(im_temp_temp)

   else

   call sorting_for_blas(31,im_order,match_index,.true.,transA,v_sort,x_row,x_column,im_temp,if_TotSym=.false., &
      irrep_input=irrep_target)  

   deallocate(im_temp)

   endif

   call daxpy (size(imTensor),A1*asym_factor*dia_sum_factor,v_sort,1,ImTensor,1)

   deallocate(v_sort)

  end subroutine

  subroutine contraction_244 (left_order,right_order,im_order,ImTensor,dia_factor, &
           &     dia_sum_factor,nrep,LeftTensor,RightTensor,irrep_left,irrep_right)

!----------calling variables-----------------------------------
    character(2),dimension(4),intent(in) :: right_order,im_order
    character(2),dimension(2),intent(in) :: left_order
    real(8),intent(in),optional :: LeftTensor(:), RightTensor(:) 
    integer,intent(in) :: nrep
    real(8),intent(in) ::  dia_sum_factor      !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(in) :: dia_factor           !dia_factor = alpha in DGEMM ; dia_sum_factor= beta in DGEMM
    real(8),intent(inout) :: imTensor(:)
    integer,intent(in),optional  :: irrep_left,irrep_right
    integer :: x_row(nrep), x_column(nrep), y_row(nrep), y_column(nrep)
    character(1) :: transA, transB
    integer :: match_index,i
    character(2) :: FreeIndx(2),sorted_opArray(4)
    real(8),allocatable :: right_sort (:)
    integer :: asym_factor
    integer             :: irrep_left_loc,irrep_right_loc, irrep_target
    integer :: maxdim 
    real(8),allocatable :: im_temp(:)
    real(8),allocatable :: im_temp_temp(:)
    integer :: im_temp_size,  im_half_size_x,  im_half_size_y

#include "param.inc"
#include "complex.inc"

  maxdim = max (sorted_dimension(right_order),sorted_dimension(im_order))

 if (.not.present(irrep_left)) then
   irrep_left_loc = 1
 else
   irrep_left_loc = irrep_left
 endif 

 if (.not.present(irrep_right)) then
   irrep_right_loc = 1
 else
   irrep_right_loc = irrep_right
 endif 

  allocate(right_sort(maxdim*rcw))  

   transA = 'N'
   transB = 'N'
!  1) fetch the array which is flexible in terms of indices. this array is also suitable for BLAS. we are also getting row and column indices for this array here.

   call decide_sorting_routine(13,right_order,im_order,left_order,transB,match_index,asym_factor, &
    &                    .true.,sorted_opArray=sorted_opArray)

 if (present(RightTensor)) then
 
   call sorting_for_blas(13,right_order,match_index,.false.,transB,right_sort,y_row,y_column,RightTensor, &
                       irrep_input=irrep_right_loc)
 
 else

   call sorting_for_blas(13,right_order,match_index,.false.,transB,right_sort,y_row,y_column, &
                    irrep_input=irrep_right_loc)   

 endif

!  2) properly sort the fixed array to contract with BLAS. row and colum indices for this array are also extracted here.

   call onebody_for_blas (left_order,im_order,sorted_opArray,transA,x_row,x_column,FreeIndx)

!  3) matrix-matrix multiplication routine runs over each irreps.

   im_temp_size = dot_product(x_row,y_column)
!
! aspg: introduced below an additional "cross"-dot product to get the dimension of the work variable below.
!
!       just doing a dot product of x_row and y_column to get the total dimesion of the variable
!       seems fine as far as i can see for groups without inversion symmetry. 
!
!       however, in the case of groups with inversion the dot product cannot allow for  calculating the dimensions of kind
!
!       x_row(irrep_g) * y_column(irrep_u) (and vice-versa) 
!
!       that may occur when one targets a "u" irrep. if the value for the "u" symmetry is smaller than the
!       value coming from the application of the dot product there will be no consequences, otherwise the xgemm 
!       call inside cntrct_new can show a segmentation fault. 
!
!       below, the final work variable dimension is then set to the maximum between the pure dot product
!       and the "cross"-dot product
!
   im_half_size_x = size(x_row,1)/2     
   im_half_size_y = size(y_column,1)/2     
   im_temp_size = max(im_temp_size, &
  & (dot_product(x_row(im_half_size_x+1:im_half_size_x*2),y_column(1:im_half_size_y)) &
  & +dot_product(x_row(1:im_half_size_x),y_column(im_half_size_y+1:im_half_size_y*2))))

   allocate(im_temp(im_temp_size*rcw))
   im_temp = 0.0d0

! if (present(LeftTensor)) call cntrct (transA,transB,x_row,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,&
! &        right_sort,A0,im_temp,nrep) 

  if (present(LeftTensor)) call cntrct_new (transA,transB,x_row,y_column,y_row,asym_factor*dia_factor*A1,LeftTensor,&
  &        right_sort,A1,im_temp,irrep_left_loc,irrep_right_loc,irrep_target,.false.) 

!!   right_sort = 0.0d0
!  4) re-sort the arrays to align with the desired structure of the imtermediate.

   transA = 'N'

   call decide_sorting_routine(13,im_order,left_order,(/FreeIndx(1),sorted_opArray(2:4)/),transA,match_index, &
   &                        asym_factor,.false.,inverse=.true.,sorted_opArray=sorted_opArray)

   if (transA=='T') then

!    write(*,*)'possibility of mistake',x_row(1),y_column(1),x_row(2),y_column(2)
  allocate(im_temp_temp(dot_product(x_row,y_column)*rcw))
!  call srt1c1_exp (irrep_target,.false.,x_row,y_column,im_temp,im_temp_temp,2)
   call srt1c1n (nrep,x_row,y_column,im_temp,im_temp_temp)

 deallocate(im_temp)
   call sorting_for_blas(13,im_order,match_index,.true.,transA,right_sort,x_row,x_column,im_temp_temp,if_TotSym=.false., &
                  irrep_input=irrep_target)   

 deallocate(im_temp_temp)
   else
   call sorting_for_blas(13,im_order,match_index,.true.,transA,right_sort,x_row,x_column,im_temp,if_TotSym=.false., &
                  irrep_input=irrep_target)   

 deallocate(im_temp)
   endif
   call daxpy (size(ImTensor),A1*asym_factor*dia_sum_factor,right_sort,1,ImTensor,1)

  deallocate(right_sort)

  end subroutine

!===========================================================================================

      subroutine cntrct_new (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine was
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)  :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)    :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)

      irrepc = multb(irrepa+nrep,irrepb+nrep,1)
      
      shift = 0
      if(shift_for_boson) shift=nrep 
       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do krp = 1, nrep
!      mrp = multb(krp+shift,irrepa+nrep,1)  
       nrp = multb(krp+shift,irrepb+nrep,1)  
       mrp = multb(irrepa+nrep,krp+shift,2)  
!        nrp = multb(irrepb+nrep,krp+shift,2)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twoDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twoDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twoDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twoDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twoDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo

      call dealloc_array(e,f,g)
      end subroutine

      subroutine cntrct_1 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine was
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)  :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)    :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      irrepc = multb(irrepa+nrep,irrepb+nrep,1)

      call alloc_array(e,nrep,f,g)
      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do mrp = 1, nrep
         krp = multb(irrepa+nrep,mrp+shift,2)  
         nrp = multb(krp+shift,irrepb+nrep,2)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twoDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twoDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twononDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twononDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twononDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo
     
      call dealloc_array(e,f,g)

      end subroutine

      subroutine cntrct_2 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine was
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)  :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)    :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)

      irrepc = multb(irrepa+nrep,irrepb+nrep,1)

      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

      do nrp = 1, nrep
         krp = multb(irrepb+nrep,nrp,2)  
         mrp = multb(krp,irrepa+nrep,2)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twononDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twononDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twoDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twoDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twononDirac(mrp,nrp)*rcw
!          off3 = 1 + g%twoDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo

      call dealloc_array(e,f,g)
      end subroutine

      subroutine cntrct_3 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine was
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)  :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)   :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)

      irrepc = multb(irrepa+nrep,irrepb+nrep,1)

      shift = 0
      if(shift_for_boson) shift=nrep 

       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif
       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif
        call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)


      do krp = 1, nrep
         mrp = multb(krp+shift,irrepa+nrep,2)  
         nrp = multb(krp+shift,irrepb+nrep,2)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1 + e%twononDirac(mrp,krp)*rcw
         else 
          off1 = 1 + e%twononDirac(krp,mrp)*rcw
         endif

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twononDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twononDirac(nrp,krp)*rcw
         endif

         off3 = 1 + g%twononDirac(mrp,nrp)*rcw

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      enddo

      call dealloc_array(e,f,g)
      end subroutine

        subroutine cntrct_4 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!  C = A*B operation, where A is a matrix and B is a vector
!  We will use dgemm for this operation as well..
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)  :: c(*)
    complex*16,intent(in) :: alpha,beta

    integer,intent(in)    :: mdim(nrep),kdim(nrep)
    integer,intent(in)    :: ndim(nrep)

!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)    :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,krp,nrp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)

      irrepc = multb(irrepa+nrep,irrepb+nrep,1)

      shift = 0
      if(shift_for_boson) shift=nrep 
       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif

       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif

       call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

         krp = irrepb
         mrp = multb(irrepa+nrep,krp+shift,2)  
!         mrp = multb(krp+shift,irrepa+nrep,1)  
         nrp = multb(krp+shift,irrepb+nrep,1)  

         m = mdim(mrp)
         n = ndim(nrp)
         k = kdim(krp)

         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif

         ldb = k
         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
          off1 = 1  + e%twoDirac(mrp,krp)*rcw
         else 
          off1 = 1  + e%twoDirac(krp,mrp)*rcw
         endif

!         off2 = 1

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 !+ f%twoDirac(krp,nrp)*rcw
         else 
          off2 = 1 !+ f%twoDirac(nrp,krp)*rcw
         endif


         off3 = 1 !+ g%twoDirac(mrp,nrp)*rcw


         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
      
!     if (transa.eq.'N'.or.transa.eq.'n') then
!      call xgemv (transa,m,k,alpha,a(off1),lda,b(off2),1,beta,c(off3),1) 
!     else
!      call xgemv (transa,k,m,alpha,a(off1),lda,b(off2),1,beta,c(off3),1) 
!     endif

      call dealloc_array(e,f,g)
      end subroutine

      subroutine cntrct_5 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!  C = A*B operation, where A is a matrix and B is a vector
!  We will use dgemm for this operation as well..
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)  :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim(nrep),kdim(nrep)
    integer,intent(in)    :: ndim
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)    :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)

      shift = 0

      if(shift_for_boson) shift=nrep 

      irrepc = multb(irrepa+shift,irrepb+shift,1)
       if (transa.eq.'N'.or.transa.eq.'n') then
        call auto_symmetry_offset(e,mdim,kdim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(e,kdim,mdim,shift_for_boson,shift_for_boson)
       endif

  !      call auto_symmetry_offset(g,mdim,ndim,shift_for_boson,shift_for_boson)

         krp = irrepb
         mrp = multb(krp+shift,irrepa+nrep,1)  

         m = mdim(mrp)
         n = ndim
         k = kdim(krp)

         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif

         ldb = k

         ldc = m

         if (transa.eq.'N'.or.transa.eq.'n') then
       off1 = 1 + e%twoDirac(mrp,krp)*rcw
         else 
       off1 = 1 + e%twoDirac(krp,mrp)*rcw
         endif

         off2 = 1 

         off3 = 1 

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif

!     if (transa.eq.'N'.or.transa.eq.'n') then
!      call xgemv (transa,m,k,alpha,a(off1),lda,b(off2),1,beta,c(off3),1) 
!     else
!      call xgemv (transa,k,m,alpha,a(off1),lda,b(off2),1,beta,c(off3),1) 
!     endif

      call dealloc_array(e,f,g)
      end subroutine


      subroutine cntrct_6 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine was
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)     :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim
    integer,intent(in)    :: ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)    :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------

      call alloc_array(e,nrep,f,g)
      shift = 0
      if(shift_for_boson) shift=nrep 

      irrepc = multb(irrepa+shift,irrepb+shift,2)

       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif

         krp = irrepa
         nrp = multb(krp+shift,irrepb+nrep,1)  

         m = mdim
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

          off1 = 1 

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twoDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twoDirac(nrp,krp)*rcw
         endif

         off3 = 1

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
     
      call dealloc_array(e,f,g)

      end subroutine

      subroutine cntrct_7 (transa,transb,mdim,ndim,kdim,alpha,a,b, &
     &                   beta,c,irrepa,irrepb,irrepc,shift_for_boson)

#include "symm.inc"
!---------------description--------------------------------------------
!     generalization of cntrct routine by Luuk. the previous routine was
!     applicable only for the tensors of totally symmetric irrep.
!     now my target is to contract the tensors of any given bosonic irreps.
!     the contraction will generate a tensor belonging to the direct product 
!     irrep of the incoming ones. 
!---------------routines called----------------------------------------

!---------------last modified------------------------------------------
!     author : avijit shee
!---------------calling variables--------------------------------------

    real*8,intent(in)     :: a(*),b(*)
    real*8,intent(inout)     :: c(*)
    complex*16,intent(in) :: alpha,beta
    integer,intent(in)    :: mdim
    integer,intent(in)    :: ndim(nrep),kdim(nrep)
!   integer,intent(in)    :: nrep           ! total number of irreps
    integer,intent(in)    :: irrepa,irrepb  ! irreps to which A & B belong.
    integer,intent(out)   :: irrepc         ! irrep to which C belong.
    logical,intent(in)    :: shift_for_boson! a logical variable indicating whether the shift for bosonic irrep is necessary 
    character*1,intent(in):: transa,transb

!---------------common blocks--------------------------------------

#include "param.inc"
#include "complex.inc"

!---------------local variables--------------------------------------
     integer :: mrp,nrp,krp 
     integer :: lda,ldb,ldc,m,n,k 
     integer :: off1,off2,off3,shift 
     type(Offset) :: e,f,g
!---------------executable code--------------------------------------


      call alloc_array(e,nrep,f,g)
      shift = 0
      if(shift_for_boson) shift=nrep 

      irrepc = multb(irrepa+shift,irrepb+shift,2)

       if (transb.eq.'N'.or.transb.eq.'n') then
        call auto_symmetry_offset(f,kdim,ndim,shift_for_boson,shift_for_boson)
       else
        call auto_symmetry_offset(f,ndim,kdim,shift_for_boson,shift_for_boson)
       endif

         krp = irrepa
         nrp = multb(krp+shift,irrepb+nrep,2)  

         m = mdim
         n = ndim(nrp)
         k = kdim(krp)
         if (transa.eq.'N'.or.transa.eq.'n') then
            lda = m
         else
            lda = k
         endif
         if (transb.eq.'N'.or.transb.eq.'n') then
            ldb = k
         else
            ldb = n
         endif
         ldc = m

          off1 = 1 

         if (transb.eq.'N'.or.transb.eq.'n') then
          off2 = 1 + f%twononDirac(krp,nrp)*rcw
         else 
          off2 = 1 + f%twononDirac(nrp,krp)*rcw
         endif

         off3 = 1

         if (k.eq.0) then
            call xscal (m*n,beta,c(off3),1)
         else
            call xgemm (transa,transb,m,n,k,alpha,a(off1),lda,b(off2),ldb,beta,c(off3),ldc)
         endif
     
      call dealloc_array(e,f,g)

      end subroutine

  end module
