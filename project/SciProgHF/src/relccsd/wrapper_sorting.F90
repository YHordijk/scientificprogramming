  module wrapper_sorting

! This module function as a wrapper around the original RELCCSD sorts that have been generalized to 
! handle also non-symmetric operators (i.e. cases in which the irrep product of the orbital indices
! need not be totally symmetric). 
! Purpose of these sorts is to change the order in which tensor elements are stored such that 
! contractions can be carried out by calls to a BLAS DGEMM or ZGEMM (for complex groups). After a
! contraction, tensors need typically be sorted back to the standard "Dirac" order, which is why
! most of them come with an inverse.

! Module was written by Avijit Shee.

  use modified_sorting
  use symmetry_offset

    implicit none

    public shuffle_indices
    public match_array
    public onebody_for_blas
    public sorting_for_blas
    public decide_sorting_routine 
    public sorted_dimension 

    private
    public indices 
    integer      :: sortType
    logical      :: do_direct, do_exchange
    character(20):: arraytype
    character(2),pointer :: array_original(:)
    character(2)         :: exchange(4)
    integer              :: reorder_left(4) 
    type sorting_list
      integer :: num 
      integer :: trans
      procedure(),pointer, nopass :: unrestricted  
      procedure(),pointer, nopass :: restricted_braket, restricted_bra, restricted_ket 
    endtype sorting_list 

    type indices
      integer, dimension(:,:),pointer :: tuple
      integer, dimension(:),pointer   :: tuple_combInd, orbital_type, full_array 
    endtype indices

    type, extends(indices) :: middle_indices
       type(indices) :: triangular 
    endtype

   type, extends(middle_indices) :: more_indices
      type(middle_indices) :: pqr,qrs,totsym
   end type

    type bit_processing
       integer :: bra,ket,braKet,I,J,K,L,ijk,jkl,ijkl
       integer :: I_reorder,J_reorder,K_reorder,L_reorder,ijkl_reorder
       integer :: bra_reorder,ket_reorder,ijk_reorder,jkl_reorder
    endtype bit_processing 

    type presorting
       real*8,allocatable :: presorted(:) 
       real*8, pointer :: integral(:)
    endtype

!  The list with the original sorters, which will be assigned on basis of the input arrays 
!  and desired resorted order.
   external SRT1SS4
   external SRT1L1
   external SRT1C1
   external SRT16
   external SRT1R1
   external SRT1TT4
   external SRT1TT5
   external SRT1ST4
   external SRT36 
   external SRT22
   external SRT32
   external SRT26
   external SRT20D
   external SRT1T3
   external SRT1T2
   external SRT1S2
   external SRT19
   external SRT1S3
   external SRT7
   external SRT6
   external SRT9
   external SRT1TS4
!  subroutines to retrieve anti-symmetrized integrals from disk storage, used when sorting
!  integrals that do not reside in memory
   external getvoov
   external getoovv
   external getooov
   external getovoo
   external getoovo
   external getvvvo_incore

    contains
    
   subroutine shuffle_indices(reorder_array,no_match,save_index,trans)

! ------- variables --------
    integer,intent(inout)        :: reorder_array(4)
    logical,intent(inout)        :: no_match
    character(1),intent(inout)   :: trans
    integer                      :: shifted_array(4)
    integer,allocatable          :: P(:,:)
    integer                      :: temp,i,j  
    integer, intent(out)         :: save_index
    type(sorting_list)           :: a
! --------------------------

    save_index = 0

    select case(sortType)
    
     case(13)

        shifted_array(1:3) =  reorder_array(2:4)
        temp = reorder_array(1)
        allocate (P(6,3))
        call permutate(shifted_array(1:3),P)

        do i = 1, size(P,1)
          a%num = toDecimal((/temp,P(i,:)/))
          a%trans = toDecimal((/temp,P(i,:)/),.true.)    
          call match_array (a%num, no_match, save_index, trans)  
          if (no_match) call match_array_transpose (a%trans, no_match, save_index, trans)  
          if (.not.no_match) then
             reorder_array(2:4)=P(i,:)
             exit
          endif
        enddo

        deallocate(P)

     case(31) 

        shifted_array(1:3) = reorder_array(1:3)
        temp = reorder_array(4)
        allocate (P(6,3))
        call permutate(shifted_array(1:3),P)

        do i = 1, size(P,1)
          a%num = toDecimal((/P(i,:),temp/))
          a%trans = toDecimal((/P(i,:),temp/),.true.)    
          call match_array (a%num, no_match, save_index, trans)  
          if (no_match) call match_array_transpose (a%trans, no_match, save_index, trans)  
          if (.not.no_match) then
             reorder_array(1:3)=P(i,:)
             exit
           endif
        enddo

        deallocate(P)

     case (22) 

        shifted_array(1:4) = reorder_array(1:4) 
        i = 1
        do
          temp = reorder_array(i)
          reorder_array(i) = reorder_array(i+1)
          reorder_array(i+1) = temp
          a%num = toDecimal(reorder_array)
          a%trans = toDecimal(reorder_array,.true.)    
          call match_array (a%num,no_match,save_index,trans)  
          if (no_match) call match_array_transpose (a%trans,no_match,save_index,trans) 
          if ((.not.no_match).or.(i>3)) exit
          i = i+2
          reorder_array(1:4) = shifted_array(1:4)
        enddo

    end select

   end subroutine

  subroutine match_array (permute,no_match,save_index,trans) 

!-----------descripton---
! choose the sorting type depending on whether the array is square or tringular and pattern of sorting (e.g, (2,2);(1,3))
!------------------------

!-----------calling variable-------------
   integer, intent(in)        :: permute 
   logical, intent(inout)     :: no_match
   integer, intent(out)       :: save_index
   character(1), intent(inout):: trans
!--------------------------------
!---------local variable---------
   integer, allocatable       :: A(:)
   integer                    :: i,dim_A,upbound_A
!--------------------------------

  save_index = 0

  select case(arraytype)

  case('unrestricted')

  select case(sortType)

  case(22) 

  allocate(A(6))
  
  A = (/1234,1324,1432,1342,2134,1243/)

  case(31)

  allocate(A(6))

  A = (/1234,1243,1342,3421,4123,3412/)

  case(13)

  allocate(A(1))

  A = (/1234/)

  end select

  case('restricted_braket')

  select case(sortType) 

  case(22)
  allocate(A(3))
  A = (/1234,1324,1342/)

  case(31)
 
  allocate(A(3))

  A= (/1234,4123,1342/)
  
  case(13) 
  allocate(A(1))

  A= (/1234/)

  end select
 
  case('restricted_bra')

  select case(sortType) 

  case(22)
  allocate(A(4))
  A = (/1234,1324,1432,1342/)

  case(31)
  allocate(A(4))
  A = (/1234,1342,1243,4123/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select
 
  case('restricted_ket')

  select case(sortType)

  case(22)

  allocate(A(3))
  A = (/1234,1324,1432/)

  case(31)

  allocate(A(5))
  A = (/1234,4123,1342,3421,3412/)


  case(13)
  allocate(A(1))
  A = (/1234/)

  end select

  end select

  dim_A = size(A)
  upbound_A = dim_A

   i = 1

   do 
   if (permute == A(i)) then
   save_index = i
   no_match = .FALSE.
   exit
   endif
     i = i+1
   if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine match_array_transpose (permute_transpose, no_match, save_index, trans) 

!-----------variable-------------
   integer,intent(in) :: permute_transpose
   logical,intent(inout):: no_match
   integer,intent(out):: save_index
   character(1),intent(inout) :: trans
!--------------------------------
!---------local variable---------
  integer, allocatable :: A(:)
  integer :: i, dim_A, upbound_A
  integer :: sort_type
!--------------------------------

  save_index = 0 

  if (sortType == 31) sort_type = 13
  if (sortType == 13) sort_type = 31
  if (sortType == 22) sort_type = 22

 select case(arraytype)

  case('unrestricted')

  select case(sort_type)

  case(22) 

  allocate(A(6))

  A = (/1234,1324,1432,1342,2134,1243/)

  case(31)

  allocate(A(6))

  A = (/1234,1243,1342,3421,4123,3412/)

  case(13)

  allocate(A(1))

  A = (/1234/)

  end select

  case('restricted_braket')

  select case(sort_type) 

  case(22)
  allocate(A(3))
  A = (/1234,1324,1342/)

  case(31)
 
  allocate(A(3))

  A= (/1234,4123,1342/)
  
  case(13) 
  allocate(A(1))

  A= (/1234/)

  end select
 
  case('restricted_bra')

  select case(sort_type) 

  case(22)
  allocate(A(4))
  A = (/1234,1324,1432,1342/)

  case(31)
  allocate(A(4))
  A = (/1234,1342,1243,4123/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select
 
  case('restricted_ket')

  select case(sort_type)

  case(22)

  allocate(A(3))
  A = (/1234,1324,1432/)

  case(31)

  allocate(A(5))
  A = (/1234,4123,1342,3421,3412/)

  case(13)
  allocate(A(1))
  A = (/1234/)

  end select

  end select

  dim_A = size(A)
  upbound_A = dim_A

   i = 1
   do 
    if (permute_transpose == A(i)) then 
      save_index = i
      no_match = .FALSE.
!      trans ='C'  
       trans ='T'  
      exit
    endif
     i = i+1
    if (i > upbound_A) exit
   enddo 

  deallocate(A)

  end subroutine

  subroutine decide_sorting_routine(sort_type,opArray_left,opArray_right,opArray_target,if_trans, &
                                &    save_index,asym_factor,do_shuffle,sorted_opArray,inverse)

!-------------------------------------
!    choose the sorting routine. 
!-------------------------------------

!--------------calling variables---------------
    character(2), intent(in),target        :: opArray_left(:)
    character(2), intent(in)               :: opArray_target(:),opArray_right(:)
    integer,      intent(in)               :: sort_type
    logical,      intent(in), optional     :: inverse  
    character(1), intent(out)              :: if_trans
    integer,      intent(out)              :: save_index 
    integer,      intent(out)              :: asym_factor
    logical,      intent(in)               :: do_shuffle
    character(2), intent(out),optional     :: sorted_opArray(size(opArray_left))
!----------------------------------------------

!-------local variables--------
    character(2) :: array_exchange(size(opArray_left))
    integer      :: reorder_v(size(opArray_left))
    integer      :: i,j,k,l 
    logical      :: no_match  
    character(2) :: p,q,r,s

    type(sorting_list)     :: a
!------------------------------

     asym_factor = 1

     sortType = sort_type       !very dirty trick. will soon amend it.

     array_original => opArray_left 

     arraytype = 'unrestricted'

     p = opArray_left(1)
     q = opArray_left(2)
     r = opArray_left(3)
     s = opArray_left(4)

    if (p(1:1) == q(1:1)) then
      if (r(1:1) == s(1:1)) then
       arraytype = 'restricted_braket'
      else
       arraytype = 'restricted_bra'
      endif 
      elseif (r(1:1)==s(1:1)) then
       arraytype = 'restricted_ket'
    endif

    do_direct = .true.
    do_exchange = .false.
    no_match = .true.

   if (present(inverse)) then 

    call align_array(array1=opArray_left,array2=opArray_target,reorder_array=reorder_v)

!findloc feature is not available for intel compilers. so I can't use it now.

    a%num = toDecimal(reorder_v)    
    a%trans = toDecimal(reorder_v,.true.)    

    call match_array (a%num,no_match,save_index,if_trans)

    if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

    endif

    if ((no_match).and.(do_direct)) then

    call align_array(opArray_left,opArray_target,opArray_right,reorder_v)

!findloc feature is not available for intel compilers. so I can't use it now.

    a%trans = toDecimal(reorder_v,.true.)    
    a%num   = toDecimal(reorder_v)    

    call match_array (a%num,no_match,save_index,if_trans)

    if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

!   if (arraytype=='unrestricted') then 

    if ((no_match).and.(do_shuffle)) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

!   endif

     if (no_match) then
!ashee      write (*,*) "====no matching sorting subroutine available for direct mode of contraction.====="
      do_exchange = .true.
     endif

    endif 

   if (do_exchange) then

! In this mode I generate all the antisymmetrized components of an array on-the-fly, wherever the interchange is possible between the same orbital-type (i.e, h or p).
! I have included the change in sign due to anti-symmetrization.

! ashee       write(*,*)'=====contraction will be done in exchange mode====='

       array_exchange = opArray_left

     if (arraytype == 'restricted_braket') then
        k = 1
        i = 1 

     do

       call swap_character(array_exchange,i,i+1)   

       call align_array(array_exchange,opArray_target,opArray_right,reorder_v)

       a%num   = toDecimal(reorder_v)
       a%trans = toDecimal(reorder_v,.true.)    

       call match_array (a%num,no_match,save_index,if_trans)
      
       if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

       if (.not.no_match) exit
       if (k>=3)           exit
       i = i + size(opArray_left)/2 
       if (k==2) then
        array_exchange = opArray_left
        i = 3
       endif
       k = k+1
     enddo

     if(.not.no_match) asym_factor = asym_factor*(-1)**(k)

!      if ((no_match).and.(do_shuffle)) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

   else

       if (arraytype == 'restricted_ket') call swap_character(array_exchange,3,4)   
       if (arraytype == 'restricted_bra') call swap_character(array_exchange,1,2)   

       call align_array(array_exchange,opArray_target,opArray_right,reorder_v)

       a%num   = toDecimal(reorder_v)
       a%trans = toDecimal(reorder_v,.true.)    

       call match_array (a%num,no_match,save_index,if_trans)
      
       if (no_match) call match_array_transpose (a%trans,no_match,save_index,if_trans)  

       if ((no_match).and.(do_shuffle)) call shuffle_indices (reorder_v,no_match,save_index,if_trans)

       if(.not.no_match) asym_factor = -1*asym_factor

    endif
   endif

    exchange = array_exchange

!  if ((no_match).and.(do_shuffle)) then

!  do_exchange = .false.

!  call align_array(opArray_left,opArray_target,opArray_right,reorder_v)

!  call shuffle_indices (reorder_v,no_match,save_index,if_trans)

!  endif 

    if (no_match) then
      call quit('No Sorting subroutine available. code stops here.')
    endif

    if (present(sorted_opArray)) then 
      if (.not.do_exchange) then
        sorted_opArray = get_operator_array(opArray_left,toDecimal(reorder_v))
      else
        sorted_opArray = get_operator_array(array_exchange,toDecimal(reorder_v))
      endif
    endif

   reorder_left = reorder_v 

  if ((if_trans=='T').and.(sortType==31)) then

   reorder_left(1)=reorder_v(4)
   reorder_left(2:4)=reorder_v(1:3)

  endif 

  if ((if_trans=='T').and.(sortType==13)) then

   reorder_left(4)=reorder_v(1)
   reorder_left(1:3)=reorder_v(2:4)

  endif 

  if ((if_trans=='T').and.(sortType==22)) then

   reorder_left(1:2)=reorder_v(3:4)
   reorder_left(3:4)=reorder_v(1:2)

  endif 

  end subroutine

  subroutine sorting_for_blas(Sort_Type,orb_array,save_index,do_inverse,if_trans,sorted,row,column,presorted,if_TotSym,irrep_input) 

!---------------description-------
!        carry-out the sorting.
!---------------------------------

!--------------calling variables---------------
    integer,intent(in)                    :: Sort_Type 
    character(2),intent(in)               :: orb_array(:)
    real(8),intent(in),target,optional    :: presorted(:)
    logical,intent(in),optional           :: if_TotSym
    integer,intent(inout)                 :: save_index
    logical,intent(in)                    :: do_inverse  
    character(1),intent(inout)            :: if_trans
    integer,intent(out)                   :: row(:),column(:)
    real(8),intent(inout)                 :: sorted(:) !think a bit more carefully to replace it by dynamic allocation. not easy.
    integer,intent(in)           :: irrep_input 
!----------------------------------------------

!-------local variables------------------------

    character(20) :: Array_Type
    logical      :: do_TotSym,if_bra,if_ket
    integer      :: i,j,k,l 
    integer      :: sort_type_dummy,AllocateStatus,total_dim

#include "symm.inc"
#include "param.inc"
#include "complex.inc"

    type(sorting_list)     :: a
    type(middle_indices)   :: c(-14:5)!,e  
    type(Offset)           :: e,f      
    type(more_indices)     :: m(0:15)  
    type(bit_processing)   :: d 
    type(presorting),target  :: V 
    type(presorting)       :: B

    type(sorting_list) :: srt_22(0:5)
    type(sorting_list) :: srt_31(0:5)
    type(sorting_list) :: srt_13(0:5)
! (2,2)->(2,2) type of sorting routines         
       srt_22(1)%unrestricted => SRT1SS4_exp   
       srt_22(2)%unrestricted => SRT16_exp     
       srt_22(3)%unrestricted => SRT46_exp 
       srt_22(4)%unrestricted => SRT1L1    
       srt_22(5)%unrestricted => SRT1R1 
       srt_22(1)%restricted_braket => SRT1TT4_exp
       srt_22(2)%restricted_braket => SRT1TT5_exp
!      srt_22(1)%restricted_bra =  SRT1LS1
       srt_22(1)%restricted_bra => SRT1TS4_exp
       srt_22(2)%restricted_bra => SRT36_exp
       srt_22(3)%restricted_bra => SRT56_exp
       srt_22(1)%restricted_ket => SRT1ST4_exp
       srt_22(2)%restricted_ket => SRT26_exp
!      srt_22(3)%restricted_ket => SRT20D

! (2,2)->(3,1) type of sorting routines.
       srt_31(0)%unrestricted => SRT1S3_exp
       srt_31(1)%unrestricted => SRT19_exp
       srt_31(2)%unrestricted => SRT6_exp
       srt_31(3)%unrestricted => SRT9_exp
       srt_31(4)%unrestricted => SRT32_exp
       srt_31(5)%unrestricted => SRT10_exp


       srt_31(0)%restricted_ket => SRT1T3_exp
       srt_31(1)%restricted_ket => SRT22_exp
       srt_31(2)%restricted_ket => SRT6_exp
       srt_31(3)%restricted_ket => SRT9_exp
       srt_31(4)%restricted_ket => SRT10_exp
       srt_31(0)%restricted_braket => SRT1T3_exp
       srt_31(1)%restricted_braket => SRT22_exp
       srt_31(2)%restricted_braket => SRT7_exp


       srt_31(0)%restricted_bra => SRT1S3_exp
       srt_31(1)%restricted_bra => SRT7_exp
       srt_31(2)%restricted_bra => SRT19_exp
       srt_31(3)%restricted_bra => SRT32_exp
    
!(2,2)->(1,3) type of sorting routines.
       srt_13(0)%unrestricted => SRT1S2_exp
       srt_13(0)%restricted_bra => SRT1T2_exp
       srt_13(0)%restricted_ket => SRT1S2_exp
       srt_13(0)%restricted_braket => SRT1T2_exp

!----------------------------------------------

    Array_Type = 'unrestricted'

    if_bra = all(orb_array(1:size(orb_array)/2)(1:1) == orb_array(1)(1:1)) 
    if_ket = all(orb_array((size(orb_array)/2)+1:size(orb_array))(1:1) == orb_array((size(orb_array)/2)+1)(1:1)) 

    if (if_ket) Array_Type = 'restricted_ket' 

    if (if_bra) then
     Array_Type = 'restricted_bra'  
     if (if_ket) Array_Type = 'restricted_braket' 
    endif 

     save_index = save_index - 1

     do_TotSym =.false.

     if (present(if_TotSym)) do_TotSym = if_TotSym

     if (.not.present(presorted)) then

     call input_index(orb_array,c,d)


     call alloc_array(e,nrep,f)
      select case(Array_Type)
 
      case('restricted_braket')

      call auto_symmetry_offset_triangular(e,c(d%I)%orbital_type,c(d%J)%orbital_type)
      call auto_symmetry_offset_triangular(f,c(d%K)%orbital_type,c(d%L)%orbital_type)

       case('restricted_bra')

      call auto_symmetry_offset_triangular(e,c(d%I)%orbital_type,c(d%J)%orbital_type)
      call auto_symmetry_offset(f,c(d%K)%orbital_type,c(d%L)%orbital_type,.false.,.false.)

       case('restricted_ket')

      call auto_symmetry_offset(e,c(d%I)%orbital_type,c(d%J)%orbital_type,.false.,.false.)
      call auto_symmetry_offset_triangular(f,c(d%K)%orbital_type,c(d%L)%orbital_type)

       case('unrestricted')

      call auto_symmetry_offset(e,c(d%I)%orbital_type,c(d%J)%orbital_type,.false.,.false.)
      call auto_symmetry_offset(f,c(d%K)%orbital_type,c(d%L)%orbital_type,.false.,.false.)

       end select

       total_dim = dot_product(e%oneNonDirac(1:nrep),f%oneNonDirac(1:nrep)) 

      call dealloc_array(e,f)

      allocate(V%presorted(total_dim*rcw))

     call integral_fetching_incore(orb_array,V%presorted) 

     B%integral => v%presorted(1:size(v%presorted)) 

     elseif (present(presorted)) then

      B%integral => presorted(1:size(presorted)) 

     endif

     sort_type_dummy = Sort_Type

     if (if_trans == 'T') then
       if (Sort_Type == 31) sort_type_dummy = 13
       if (Sort_Type == 13) sort_type_dummy = 31
     endif

    if (Sort_Type == 22) then

     if (save_index > 0) then
      select case(Array_Type)
        
         case('unrestricted')

       call input_index(orb_array,c,d)

    if (if_trans=='N') then
      call srt_22(save_index)%unrestricted (irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &     c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)
    else
      call srt_22(save_index)%unrestricted (irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &     c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)
    endif

         case('restricted_braket') 

      call input_index(orb_array,c,d)


    if (if_trans=='N') then
     call srt_22(save_index)%restricted_braket(irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type,&
    &    c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)
    else
     call srt_22(save_index)%restricted_braket(irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type,&
    &    c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)
    endif

         case('restricted_bra')

      call input_index(orb_array,c,d)


    if (if_trans=='N') then
      call srt_22(save_index)%restricted_bra(irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     & c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)  
    else
      call srt_22(save_index)%restricted_bra(irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     & c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)  
    endif

      case('restricted_ket')

      call input_index(orb_array,c,d)


    if (if_trans=='N') then
      call srt_22(save_index)%restricted_ket(irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     & c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)  
    else
      call srt_22(save_index)%restricted_ket(irrep_input,nrep,multb,do_totsym,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     & c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)  
    endif

       end select
     else

     if (present(presorted)) then

     call dcopy(size(presorted),B%integral,1,sorted,1) 

     else

     call dcopy(size(V%presorted),B%integral,1,sorted,1) 

     endif

     call input_index(orb_array,c,d)


     call alloc_array(e,nrep,f)
      select case(Array_Type)
 
      case('restricted_braket')

      call auto_symmetry_offset_triangular(e,c(d%I)%orbital_type,c(d%J)%orbital_type)
      call auto_symmetry_offset_triangular(f,c(d%K)%orbital_type,c(d%L)%orbital_type)

       case('restricted_bra')

      call auto_symmetry_offset_triangular(e,c(d%I)%orbital_type,c(d%J)%orbital_type)
      call auto_symmetry_offset(f,c(d%K)%orbital_type,c(d%L)%orbital_type,.false.,.false.)

       case('restricted_ket')

      call auto_symmetry_offset(e,c(d%I)%orbital_type,c(d%J)%orbital_type,.false.,.false.)
      call auto_symmetry_offset_triangular(f,c(d%K)%orbital_type,c(d%L)%orbital_type)

       case('unrestricted')

      call auto_symmetry_offset(e,c(d%I)%orbital_type,c(d%J)%orbital_type,.false.,.false.)
      call auto_symmetry_offset(f,c(d%K)%orbital_type,c(d%L)%orbital_type,.false.,.false.)

       end select

       if (if_trans=='N') then
       row = e%oneNonDirac 
       column = f%oneNonDirac   
       else 
       column = e%oneNonDirac 
       row = f%oneNonDirac   
       endif

     call dealloc_array(e,f)
     endif

     elseif (sort_type_dummy == 31) then
       select case(Array_Type)
      
       case('unrestricted')

       call input_index(orb_array,c,d)


      if (if_trans=='N') then

      call srt_31(save_index)%unrestricted(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

      else

      call srt_31(save_index)%unrestricted(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

      endif

       case('restricted_braket')

       call input_index(orb_array,c,d)


      if (if_trans=='N') then

      call srt_31(save_index)%restricted_braket(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

      else

      call srt_31(save_index)%restricted_braket(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

      endif

        case('restricted_bra')

       call input_index(orb_array,c,d)


      if (if_trans=='N') then

      call srt_31(save_index)%restricted_bra(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

      else

      call srt_31(save_index)%restricted_bra(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

      endif

         case('restricted_ket')

       call input_index(orb_array,c,d)

      if (if_trans=='N') then

      call srt_31(save_index)%restricted_ket(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

      else

      call srt_31(save_index)%restricted_ket(irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
     &               c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

      endif

      end select

     elseif (sort_type_dummy == 13) then
      select case(Array_type)
        
         case('unrestricted')

     call input_index(orb_array,c,d)


      if (if_trans=='N') then

        call srt_13(save_index)%unrestricted (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

     else

        call srt_13(save_index)%unrestricted (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

     endif 

       case('restricted_braket')

      call input_index(orb_array,c,d)


      if (if_trans=='N') then

        call srt_13(save_index)%restricted_braket (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

     else

        call srt_13(save_index)%restricted_braket (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

     endif 

         case('restricted_bra')

         call input_index(orb_array,c,d)


      if (if_trans=='N') then

        call srt_13(save_index)%restricted_bra (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

     else

        call srt_13(save_index)%restricted_bra (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

     endif 

      case('restricted_ket')

     call input_index(orb_array,c,d)

      if (if_trans=='N') then

        call srt_13(save_index)%restricted_ket (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,row,column)

     else

        call srt_13(save_index)%restricted_ket (irrep_input,nrep,multb,do_inverse,c(d%I)%orbital_type,c(d%J)%orbital_type, &
       &                          c(d%K)%orbital_type,c(d%L)%orbital_type,B%integral,sorted,column,row)

     endif 

      end select
   endif
     
     if (allocated(V%presorted)) deallocate(V%presorted)
     if (associated(B%integral)) nullify(B%integral)
 
  end subroutine

  subroutine onebody_for_blas (array1,array2,array3,trans_1b,row,column,FreeIndx) 

!--------purpose-------------
!    arrange onebody routine for blas. get the indices for matrix multiplication as well as the end tensor.
!----------------------------


!------------------------------
   character(*),dimension(:),intent(in)  :: array1,array2,array3 ! follows the same instruction as it is in the align_array.
   character(*),intent(out)              :: trans_1b
   character(2),optional,intent(out)     :: FreeIndx(size(array1))
   integer                               :: reorder_array(size(array1))      
   integer,intent(out)                   :: row(:), column(:) 
   type(middle_indices)                  :: a(-14:2)  
   integer,dimension(2)                  :: bit_reordered
   integer                               :: i
!------------------------------
#include "symm.inc"

  call align_array(array1,array2,array3,reorder_array)   

  if ((toDecimal(reorder_array)) == 21) then
  trans_1b = 'T'
  else
  trans_1b = 'N'
  endif
  
   a(0)%orbital_type(1:nrep) => NO(1:nrep)
   a(1)%orbital_type(1:nrep) => NV(1:nrep)
   a(-12)%orbital_type(1:nrep) => NCONT(1:nrep)

   bit_reordered = bitString(array1,reorder_array) 

   row = a(bit_reordered(1))%orbital_type
   column = a(bit_reordered(2))%orbital_type

   if (present(FreeIndx)) then

   FreeIndx(1:size(array1)) = array1(reorder_array(1:size(array1)))

  endif 

  end subroutine 

  subroutine swap_character(array,k,l)

    character(*),intent(inout) :: array(:)
    integer, intent(in)        :: k,l  
    character(2)               :: temp

    temp = array(k)

    array(k) = array(l)
    array(l) = temp
  end subroutine

  subroutine align_array(array1,array2,array3,reorder_array)

!-------------------------------------------------------------
    character(*),intent(in) :: array1(:), array2(:)                    ! array1 is the array to be reordered
    character(*),intent(in),optional :: array3(:) 
    integer, dimension(size(array1)),intent(out) :: reorder_array      ! array2 is the array from where we will extract the first index of the reordered 
                                                                       ! array i.e, for left array it is always intermediate and for right always    
    integer                                  :: k,i                    ! left array. array3 is the array to extract rest of the indices e.g, for left array
                                                                       ! always right and for right array always intermediate.
!-------------------------------------------------------------
     k = 1
     do i = 1, size(array2)       ! loop over vorder 
        if (findloc(array1,array2(i)) > 0) then 
          reorder_array(k) = findloc(array1, array2(i))
          k=k+1
        endif
     enddo

    if (present(array3)) then 
     do i = 1, size(array3)       ! loop over vorder 
        if (findloc(array1,array3(i)) > 0) then 
          if (k<=4) reorder_array(k) = findloc(array1,array3(i))
          if (k<4) k=k+1
        endif
     enddo
    endif
    
  end subroutine

  subroutine integral_fetching_incore(orb_array,array)
!-------------------------------------------------------------
  integer :: value
  character(2),intent(in) :: orb_array(:)
  integer:: bit_original(size(orb_array))
  real*8,intent(out) :: array(:)
!-------------------------------------------------------------

    bit_original = bitString (orb_array,(/1,2,3,4/)) 

    value = bin2dec(bit_original)

!here you can think of keeping all the VVVO integrals in fast memory. 

   select case(value) 

    case (0)

      call getoooo(array)

    case (1)

      call getooov(array)

    case (2)

      call getoovo(array) 

    case (3)

      call getoovv(array)

   case (4)

      call getovoo(array)

    case (5)

!     allocate(G%presorted(nv4*rcw),stat = AllocateStatus)
!     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     call getovov(G%presorted)

      write(*,*)'not yet implemented'

    case (6)

!     allocate(G%presorted(nv4*rcw),stat = AllocateStatus)
!     IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!     call getovvo(G%presorted)


      write(*,*)'not yet implemented'

    case (7)

      call getovvv_incore(array)

!      write(*,*)'ovvv type array has to be fetched in out-of-core way.'

    case (8)

      call getvooo(array)

    case (9)

      call getvoov(array)

    case (10)

      call getvovo(array)

    case (11)

      call getvovv_incore(array)

    case (12)

      call getvvoo(array)

    case (13)

      call getvvov_incore(array)

    case (14)

      call getvvvo_incore(array)

    case default
      write(*,*)'no suitable array to fetch'
   end select

  end subroutine

  subroutine input_index(orb_array,a,b)
#include "symm.inc"
! This subroutine identifies the dimension of each individual
! orbital-type (i.e, hole, particle, continuous, active etc.) 

!-------------------------------
 character(2),intent(in)                   :: orb_array(:)
 type(middle_indices),intent(out)     :: a(-14:5)  
 type(bit_processing),intent(out)     :: b 
 integer,allocatable :: bit_original(:)
!-------------------------------

   allocate(bit_original(size(orb_array)))

    bit_original  = bitString (orb_array,(/1,2,3,4/)) 

    b%I = bit_original(1) 
    b%J = bit_original(2)
    b%K = bit_original(3)
    b%L = bit_original(4)

    a(0)%orbital_type(1:nrep) => NO(1:nrep)
    a(1)%orbital_type(1:nrep) => NV(1:nrep)
    a(-12)%orbital_type(1:nrep) => NCONT(1:nrep)

   deallocate(bit_original)

  end subroutine

  integer function findloc(array,value)

  character(*), dimension(:),intent(in) :: array
  character(2),intent(in) :: value
  integer      :: i  

     do i = 1, size(array)        
        if (array(i)==value) then
          findloc = i
          exit
          else
          findloc = 0
        endif
     enddo
   end function

  function bitString (generic_array,reorder_index_array)

  character(2), dimension(:),intent(in) :: generic_array
  integer,      dimension(:),intent(in) :: reorder_index_array
  integer      :: i,k
  integer, dimension(size(generic_array)) ::bitString 
  character(2) :: temp

     do i = 1, size(generic_array)
        temp=generic_array(reorder_index_array(i))
        k = ichar(temp(1:1))-ichar('o')
        bitString(i)=k
     enddo

  end function   

  integer function toDecimal(index_array,trans)
   
  integer,dimension(:),intent(in),target :: index_array
  logical,optional, intent(in) :: trans
  integer :: scratch(size(index_array))
  integer,pointer :: scratch_bra(:)
  integer,pointer :: scratch_ket(:)
  logical                       :: do_trans 
  integer                       :: temp, i
  
  scratch = index_array
  if (present(trans)) then
   do_trans = trans
  else
   do_trans = .false.
  endif

  if (do_trans) then

! it should not be case dependent. if possible try to modify it.

  select case (sortType)
     case (22)
     scratch_bra => index_array(3:4) 
     scratch_ket => index_array(1:2) 
     case (31)
     scratch_bra => index_array(4:4)
     scratch_ket => index_array(1:3)
     case (13)
     scratch_bra => index_array(2:4)
     scratch_ket => index_array(1:1)
  end select  

     scratch = (/scratch_bra,scratch_ket/)  
  
  endif   

    toDecimal=0

    do i = size(index_array), 1 , -1  
      toDecimal = toDecimal + scratch(i)*10**(size(index_array)-i)
    enddo 

  end function

  integer function bin2dec(bitstring)

    integer, intent(in) :: bitstring(:)
    integer             :: i 

    bin2dec=0
 
    do i = 1, size(bitstring)
     bin2dec = bin2dec  + bitstring(i)*2**(size(bitstring)-i)
    enddo

  end function

  recursive subroutine permutate(E, P) 
    integer, intent(in)  :: E(:)       ! array of objects 
    integer, intent(out) :: P(:,:)     ! permutations of E 
    integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
       N = size(E); Nfac = size(P,1); 
     do i=1,N                           ! cases with E(i) in front 
          if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
!
!    Fix introduced by Tjerk Straatsma:
!
!    The following is a work around for the PGI compiler error generated
!    by the original forall loop
!
      do k=1,Nfac/N
           P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
      enddo
!
!         forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
!
     end do 
  end subroutine permutate 

  function get_operator_array(character_array,integer_for_array)
  integer,intent(in) :: integer_for_array
  character(*),intent(in) :: character_array(:)  
  character(2)            :: get_operator_array(size(character_array))  
  integer :: y, temp
  integer :: x, i

  temp = integer_for_array
  i = 1

  do  
     y = temp/10**(size(character_array)-i)  
     x = mod(temp, 10**(size(character_array)-i))
     temp = x
     get_operator_array(i) = character_array(y)
     i = i+1
     if (i > size(character_array)) exit
  enddo 
  end function

  character(len=10) function seqstring(seqnum) 
    integer, intent(in) :: seqnum
    character(10)       :: temp
    write (temp,'(I0)') seqnum
    seqstring = trim(temp)
  end function

  function sorted_dimension(generic_array)

#include "symm.inc"

    character(*),intent(in) :: generic_array(:)
    integer :: sorted_dimension
    integer :: k,i,temp(4)

    k = 0 

    do i = 1, size(generic_array)
    
       temp = bitstring(generic_array,(/1,2,3,4/)) 
       if (temp(i) == 1) then
         k = k+1
       endif

    enddo   

    select case(k)

       case(0)
         
       sorted_dimension = joooo(nrep+1)  

       case(1)
         
       sorted_dimension = jooov(nrep+1)  

       case(2)

       sorted_dimension = jvovo(nrep+1)  

       case(3)

       sorted_dimension = jvvov(nrep+1) 

       case default

       write(*,*)'needs very large array, NOT an in-core case'

     end select  

  end function


  end module   
