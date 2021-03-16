!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!
!  new module for symmetry and GASpace handling of mo-spinors in KR/(KU?)-CI
!  written by sknecht, March 2011
!
!
!  notes on notation: a. arrays/values without any particular ending (or 1) 
!                        refer to the set of unbarred spinors
!                     b. arrays/values ending with "2"
!                        refer to the set of barred spinors

module mospinor_info

   use memory_allocator
   implicit none

   public mospinor_info_init
   public mospinor_info_delete

   private

   integer, parameter,   public :: max_num_dbg_irreps = 128 ! max. number of double group irreps (in sync with MXNDGIRR)
   integer, parameter,   public :: max_num_spinors    = 800 ! max. number of spinors  (in sync with MXPORB)
   integer, parameter,   public :: max_num_gaspaces   =   8 ! max. number of GASpaces (in sync with MXPNGAS)

   integer,              public :: norb0, norb1, norb2, norb3, norb4
   integer,              public :: nacob , nocob , ntoob
   integer,              public :: nacob2, nocob2, ntoob2
   
!  allocatable arrays/matrices
!  a. unbarred spinors
   integer, allocatable, public :: nocobs(:)
   integer, allocatable, public :: ntoobs(:)
   integer, allocatable, public :: itoobs(:)

   integer, allocatable, public :: ireots(:)
   integer, allocatable, public :: ismfto(:)
   integer, allocatable, public :: itpfto(:)

   integer, allocatable, public :: nobpt(:)
   integer, allocatable, public :: ngsobt(:)

   integer, allocatable, public :: nobpts(:,:)
   integer, allocatable, public :: iobpts(:,:)
   integer, allocatable, public :: ngsob(:,:)

   integer, allocatable, public :: imosp_dirac_counter1(:)
   integer, allocatable, public :: imosp_dirac_mjub(:)
   integer, allocatable, public :: imosp_luci2dirac1(:)

!  b. barred spinors
   integer, allocatable, public :: nocobs2(:)
   integer, allocatable, public :: ntoobs2(:)
   integer, allocatable, public :: itoobs2(:)

   integer, allocatable, public :: ireots2(:)
   integer, allocatable, public :: ismfto2(:)
   integer, allocatable, public :: itpfto2(:)

   integer, allocatable, public :: nobpt2(:)
   integer, allocatable, public :: ngsobt2(:)

   integer, allocatable, public :: nobpts2(:,:)
   integer, allocatable, public :: iobpts2(:,:)
   integer, allocatable, public :: ngsob2(:,:)

   integer, allocatable, public :: imosp_dirac_counter2(:)
   integer, allocatable, public :: imosp_dirac_mjb(:)
   integer, allocatable, public :: imosp_luci2dirac2(:)
contains

!*******************************************************************************
   subroutine mospinor_info_init()

!>    reset arrays
      call mospinor_info_delete()

!>    allocate arrays

!     1-d arrays with max_num_dbg_irreps   
      call alloc(nocobs ,max_num_dbg_irreps,'occ-spinors1')
      call alloc(ntoobs ,max_num_dbg_irreps,'num-spinors1')
      call alloc(itoobs ,max_num_dbg_irreps,'off-spinors1')
      call alloc(nocobs2,max_num_dbg_irreps,'occ-spinors2')
      call alloc(ntoobs2,max_num_dbg_irreps,'num-spinors2')
      call alloc(itoobs2,max_num_dbg_irreps,'off-spinors2')

!     1-d arrays with max_num_spinors
      call alloc(ireots ,max_num_spinors,'reorder-type-sym1')
      call alloc(ismfto ,max_num_spinors,'reorder-sym1')
      call alloc(itpfto ,max_num_spinors,'reorder-type1')
      call alloc(ireots2,max_num_spinors,'reorder-type-sym2')
      call alloc(ismfto2,max_num_spinors,'reorder-sym2')
      call alloc(itpfto2,max_num_spinors,'reorder-type2')

      call alloc(imosp_dirac_counter1,max_num_spinors,'dirac_mosp-counter1')
      call alloc(imosp_dirac_counter2,max_num_spinors,'dirac_mosp-counter2')
      call alloc(imosp_dirac_mjub,    max_num_spinors,'dirac_mosp-mjub')
      call alloc(imosp_dirac_mjb,     max_num_spinors,'dirac_mosp-mjb')
      call alloc(imosp_luci2dirac1   ,max_num_spinors,'luci2dirac-mosp-index1')
      call alloc(imosp_luci2dirac2   ,max_num_spinors,'luci2dirac-mosp-index2')
      imosp_dirac_counter1 = 0
      imosp_dirac_counter2 = 0
      imosp_dirac_mjub     = 0
      imosp_dirac_mjb      = 0
      imosp_luci2dirac1    = 0
      imosp_luci2dirac2    = 0

!     1-d arrays with max_num_gaspaces
      call alloc(nobpt , max_num_gaspaces,'type-GASpace1')
      call alloc(nobpt2, max_num_gaspaces,'type-GASpace2')
      call alloc(ngsobt, max_num_gaspaces,'type-GASpace3')
      call alloc(ngsobt2,max_num_gaspaces,'type-GASpace4')

!     2-d arrays with max_num_gaspaces,max_num_dbg_irreps
      call alloc(nobpts ,max_num_gaspaces,max_num_dbg_irreps,'GASpace-type-sym1')
      call alloc(iobpts ,max_num_gaspaces,max_num_dbg_irreps,'GASpace-type-sym-off1')
      call alloc(nobpts2,max_num_gaspaces,max_num_dbg_irreps,'GASpace-type-sym2')
      call alloc(iobpts2,max_num_gaspaces,max_num_dbg_irreps,'GASpace-type-sym-off2')

!     2-d arrays with max_num_dbg_irreps,max_num_gaspaces
      call alloc(ngsob  ,max_num_dbg_irreps,max_num_gaspaces,'sym1-GASpace-type')
      call alloc(ngsob2 ,max_num_dbg_irreps,max_num_gaspaces,'sym2-GASpace-type')

   end subroutine mospinor_info_init

!*******************************************************************************

   subroutine mospinor_info_delete()
!-------------------------------------------------------------------------------
!
!  purpose: deallocate all arrays related to symmetry and GASpace handling 
!           of mo-spinors in KR/(KU?)-CI.
!-------------------------------------------------------------------------------

      if (allocated(nocobs )) call dealloc( nocobs )
      if (allocated(ntoobs )) call dealloc( ntoobs )
      if (allocated(itoobs )) call dealloc( itoobs )
      if (allocated(nocobs2)) call dealloc( nocobs2)
      if (allocated(ntoobs2)) call dealloc( ntoobs2)
      if (allocated(itoobs2)) call dealloc( itoobs2)

      if (allocated(ireots )) call dealloc( ireots )
      if (allocated(ismfto )) call dealloc( ismfto )
      if (allocated(itpfto )) call dealloc( itpfto )
      if (allocated(ireots2)) call dealloc( ireots2)
      if (allocated(ismfto2)) call dealloc( ismfto2)
      if (allocated(itpfto2)) call dealloc( itpfto2)

      if (allocated(imosp_dirac_counter1)) call dealloc(imosp_dirac_counter1)
      if (allocated(imosp_dirac_counter2)) call dealloc(imosp_dirac_counter2)
      if (allocated(imosp_dirac_mjub))     call dealloc(imosp_dirac_mjub)
      if (allocated(imosp_dirac_mjb))      call dealloc(imosp_dirac_mjb)
      if (allocated(imosp_luci2dirac1   )) call dealloc(imosp_luci2dirac1)
      if (allocated(imosp_luci2dirac2   )) call dealloc(imosp_luci2dirac2)


      if (allocated(nobpt  )) call dealloc( nobpt )
      if (allocated(ngsobt )) call dealloc( ngsobt)
      if (allocated(nobpt2 )) call dealloc( nobpt2)
      if (allocated(ngsobt2)) call dealloc( ngsobt2)

      if (allocated(nobpts )) call dealloc( nobpts )
      if (allocated(iobpts )) call dealloc( iobpts )
      if (allocated(ngsob  )) call dealloc( ngsob  )
      if (allocated(nobpts2)) call dealloc( nobpts2)
      if (allocated(iobpts2)) call dealloc( iobpts2)
      if (allocated(ngsob2))  call dealloc( ngsob2)

   end subroutine mospinor_info_delete

end module
