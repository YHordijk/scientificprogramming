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

!*******************************************************************************
!
!  new module for symmetry handling in KR-CI
!
!  --> simplification of existing symmetry handling
!  --> implementation of linear symmetry
!
!      written by Stefan Knecht with contributions from 
!      - Hans Joergen Aagaard Jensen
!      - Timo Fleig
!         October 2010 - ?
!
!*******************************************************************************

module symmetry_setup_krci

   use memory_allocator
   implicit none

   public match_spinor_2_boson_irrep
   public symmetry_setup_init
   public symmetry_setup_delete
   public convert_mj_ferm_2_dbg_irrep_linsym_krci

   private

   integer,              public :: lbound_mat = -1
   integer,              public :: ubound_mat = -1
   integer,              public :: nirr_pn    = -1 ! number of spatial irreps
   integer,              public :: nirr_dg    = -1 ! number of double group irreps
   integer,              public :: pntgrp     = -1 ! point  group ID
   integer,              public :: dougrp     = -1 ! double group ID
   integer,              public :: nsmsx      = -1 ! number of symmetries of single excitations
   integer,              public :: nsmdx      = -1 ! number of symmetries of double excitationss
   integer,              public :: nsmst      = -1 ! number of symmetries of strings
   integer,              public :: nsmci      = -1 ! number of symmetries of CI spaces
   integer,              public :: nsmxt      = -1 ! ?
   integer,              public :: nsmob      = -1 ! number of spinor symmetries
   integer,              public :: itssx      = -1 ! total symmetry of single excitation
   integer,              public :: itsdx      = -1 ! total symmetry of double excitation
   integer,              public :: itsxt      = -1 ! ?
!  allocatable arrays/matrices
   integer, allocatable, public :: idbgmult(:,:)
   integer, allocatable, public :: mj2rep(:,:)
   integer, allocatable, public :: iadjsym(:)
   integer, allocatable, public :: invelm(:)

contains

!*******************************************************************************
   subroutine symmetry_setup_init(double_group,number_of_dbg_irreps)

!     --------------------------------------------------------------------------
      integer, intent(in) :: double_group
      integer, intent(in) :: number_of_dbg_irreps
!     --------------------------------------------------------------------------

!>    reset arrays
      call symmetry_setup_delete()

!>    allocate arrays

!     double group multiplication table
      call alloc(idbgmult,number_of_dbg_irreps,number_of_dbg_irreps,'doubleGmult')

!     adjoint and inverse elements
      call alloc(iadjsym,number_of_dbg_irreps,'adjoint_elm')
      call alloc(invelm, number_of_dbg_irreps,'inverse_elm')

!     mj to spinor representation (irrep)
      call alloc(mj2rep,(number_of_dbg_irreps+1),2)

!     set array boundaries
      lbound_mat = 1
      ubound_mat = number_of_dbg_irreps

!>    fill arrays
      call define_dbg_mult_tables_krci(double_group,number_of_dbg_irreps)

   end subroutine symmetry_setup_init

!*******************************************************************************
   subroutine define_dbg_mult_tables_krci(double_group,number_of_dbg_irreps)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: double_group
      integer, intent(in)  :: number_of_dbg_irreps
!     --------------------------------------------------------------------------
      integer              :: i, jj, ii
      integer              :: nr_boson_irp
      integer              :: nr_gerade
!     --------------------------------------------------------------------------

      nr_boson_irp = number_of_dbg_irreps/2
     
!     define multiplication tables for each double group
!     --------------------------------------------------

      select case(double_group)

!       C1 double group
        case(8)

        call quit('C1 symmetry is currently buggy. contact stefan (knecht@ifk.sdu.dk) to ask for help if needed.')

!       adjoints and inverse elements
!       -----------------------------
        do i = 1, number_of_dbg_irreps
          iadjsym(i) = i
          invelm(i)  = i
        end do

!       double group multiplication table
!       ---------------------------------
        idbgmult(1,1) = 1
        idbgmult(1,2) = 2
        idbgmult(2,1) = 2
        idbgmult(2,2) = 1

!       Ci double group
        case(7)


!       adjoints and inverse elements
!       -----------------------------
        do i = 1, number_of_dbg_irreps
          iadjsym(i) = i
          invelm(i)  = i
        end do

!       double group multiplication table
!       ---------------------------------
        idbgmult(1,1) = 1
        idbgmult(1,2) = 2
        idbgmult(1,3) = 3
        idbgmult(1,4) = 4
        idbgmult(2,1) = 2
        idbgmult(2,2) = 1
        idbgmult(2,3) = 4
        idbgmult(2,4) = 3
        idbgmult(3,1) = 3
        idbgmult(3,2) = 4
        idbgmult(3,3) = 1
        idbgmult(3,4) = 2
        idbgmult(4,1) = 4
        idbgmult(4,2) = 3
        idbgmult(4,3) = 2
        idbgmult(4,4) = 1

!       Cs and C2 double groups
        case(6,5)

!       adjoints and inverse elements
!       -----------------------------
!       boson irreps
        do i = 1, number_of_dbg_irreps/2
          iadjsym(i) = i
        end do
!       fermion irreps
        do i = (number_of_dbg_irreps/2)+1,number_of_dbg_irreps
          iadjsym(i) = i + (-1)**(i-1)
        end do

        call icopy(number_of_dbg_irreps,iadjsym,1,invelm,1)

!       double group multiplication table
!       ---------------------------------
        idbgmult(1,1) = 1
        idbgmult(1,2) = 2
        idbgmult(1,3) = 3
        idbgmult(1,4) = 4
        idbgmult(2,1) = 2
        idbgmult(2,2) = 1
        idbgmult(2,3) = 4
        idbgmult(2,4) = 3
        idbgmult(3,1) = 3
        idbgmult(3,2) = 4
        idbgmult(3,3) = 2
        idbgmult(3,4) = 1
        idbgmult(4,1) = 4
        idbgmult(4,2) = 3
        idbgmult(4,3) = 1
        idbgmult(4,4) = 2

!       C2h double group
        case(4)

!       adjoints and inverse elements
!       -----------------------------
!       boson irreps
        do i = 1, number_of_dbg_irreps/2
          iadjsym(i) = i
        end do
!       fermion irreps
        do i = (number_of_dbg_irreps/2)+1,number_of_dbg_irreps
          iadjsym(i) = i + (-1)**(i-1)
        end do

        call icopy(number_of_dbg_irreps,iadjsym,1,invelm,1)

!       double group multiplication table
!       ---------------------------------
        do jj = 1, number_of_dbg_irreps
           do ii= 1, number_of_dbg_irreps
!             both are bosonic
              if(ii <= 4 .and. jj <= 4)then
                if(ii == jj)then
                  idbgmult(ii,jj) = 1
                else if((ii+jj) == 3 .or. (ii+jj) == 7)then
                  idbgmult(ii,jj) = 2
                else if((ii+jj) == 5)then
                  idbgmult(ii,jj) = 4
                else
                  idbgmult(ii,jj) = 3
                end if
!             both are fermionic
              else if(ii > 4 .and. jj > 4)then
                if(ii == jj)then
                  idbgmult(ii,jj) = 2
                else if((ii+jj) == 13)then
                  idbgmult(ii,jj) = 3
                else if((ii+jj) == 11 .or. (ii+jj) == 15)then
                  idbgmult(ii,jj) = 1
                else
                  idbgmult(ii,jj) = 4
                end if
!             mixed fermionic/bosonic
              else
                if((ii+jj) == 9)then
                  idbgmult(ii,jj) = 8
                else if((ii+jj) == 7 .or. (ii+jj) == 11)then
                  idbgmult(ii,jj) = 6
                else if((ii+jj) == 6 .or. (ii+jj) == 12)then
                  idbgmult(ii,jj) = 5
                else if(ii == 6 .and. jj == 2)then
                  idbgmult(ii,jj) = 5
                else if(ii == 2 .and. jj == 6)then
                  idbgmult(ii,jj) = 5
                else if(ii == 3 .and. jj == 7)then
                  idbgmult(ii,jj) = 5
                else if(ii == 7 .and. jj == 3)then
                  idbgmult(ii,jj) = 5
                else
                  idbgmult(ii,jj) = 7
                end if
              end if
           end do
        end do

!       Cinf double group
        case(10)

!       cinf mapped to c32; TODO: increase mapping for linear symmetry
!       to higher double group
!#warning increase mapping for linear symmetry to higher double group

!       adjoints and inverse elements
!       -----------------------------

!       boson irreps

        iadjsym(1) = 1

        do i = 2, nr_boson_irp - 1
          iadjsym(i) = i + (-1)**i
        end do

        iadjsym(nr_boson_irp) = nr_boson_irp

!       fermion irreps
        do i = 1, nr_boson_irp
          iadjsym(nr_boson_irp + i) = nr_boson_irp + i + (-1)**(i-1)
        end do

        call icopy(nr_boson_irp*2,iadjsym,1,invelm,1)

!       double group multiplication table
!       ---------------------------------
        call set_dbg_mult_table_linsym_krci(double_group,            &
                                            idbgmult,                &
                                            mj2rep,                  & 
                                            number_of_dbg_irreps)

!       Cinfh double group
        case(11)

!       cinfh mapped to c16h; TODO: increase mapping for linear symmetry
!       to higher double group
!#warning increase mapping for linear symmetry to higher double group

        nr_gerade = nr_boson_irp / 2

!       adjoints and inverse elements
!       -----------------------------

!       boson irreps

        iadjsym(1)           = 1
        iadjsym(nr_gerade+1) = nr_gerade+1

        do i = 2, nr_gerade - 1
          iadjsym(i)           = i + (-1)**i
          iadjsym(nr_gerade+i) = nr_gerade+i + (-1)**i
        end do

        iadjsym(nr_gerade)    = nr_gerade
        iadjsym(nr_boson_irp) = nr_boson_irp

!       fermion irreps

        do i = 1, nr_boson_irp
          iadjsym(nr_boson_irp+i) = nr_boson_irp + i + (-1)**(i-1)
        end do

        call icopy(nr_boson_irp*2,iadjsym,1,invelm,1)

!       double group multiplication table
!       ---------------------------------
        call set_dbg_mult_table_linsym_krci(double_group,            &
                                            idbgmult,                &
                                            mj2rep,                  & 
                                            number_of_dbg_irreps)

      end select
 
#if defined LUCI_DEBUG
      write(6,*) '   symmetry information from define_dbg_mult_tables_krci'
      write(6,*) '   ====================================================='
      write(6,'(/a,i5)')  '    double group    : ', double_group
      write(6,'(/a,i5/)') '    number of irreps: ', number_of_dbg_irreps
      write(6,*) '   multiplication table '
      write(6,*) '   ==================== '

      if(double_group > 9)then
        do ii = 1, number_of_dbg_irreps
          print *, ' irrep ',ii
          write(6,'(2x,64i3)') (idbgmult(ii,jj),jj=1,number_of_dbg_irreps)
        end do
      else
        do ii = 1, number_of_dbg_irreps
          print *, ' irrep ',ii
          write(6,'(2x,8i3)') (idbgmult(ii,jj),jj=1,number_of_dbg_irreps)
        end do
      end if
!
      write(6,*) '   adjoint symmetry array '
      write(6,*) '   ======================'
      write(6,'(10x,8i3)') (iadjsym(i),i=1,number_of_dbg_irreps)
!
      write(6,*) '   inverse elements '
      write(6,*) '   ================='
      write(6,'(10x,8i3)') (invelm(i),i=1,number_of_dbg_irreps)
#endif

   end subroutine define_dbg_mult_tables_krci

!*******************************************************************************
   subroutine set_dbg_mult_table_linsym_krci(double_group,            &
                                             dbg_mult_tab,            &
                                             mj2rep,                  & 
                                             num_dbg_irreps)
!-------------------------------------------------------------------------------
!
! purpose: define multiplication table for approximate linear symmetry in KR-CI.
!          the present double groups are C32 and C16h, respectively.
!          this allows a maximum M_J value of:
!          --> +- 16   (C32;  double group 10) 
!          --> +- 8g/u (C16h; double group 11)
!
!          the multiplication table follow the notation in 
!          koster, dimmock, wheeler, statz: "properties of the 32 point groups".
!
!          routine inspired by gmultln(moltra/traone.F) originally written by 
!          L. Visscher
!-------------------------------------------------------------------------------
      integer, intent(out) :: dbg_mult_tab(num_dbg_irreps,num_dbg_irreps)
      integer, intent(out) :: mj2rep(-num_dbg_irreps/2:num_dbg_irreps/2,2)
      integer, intent(in)  :: double_group
      integer, intent(in)  :: num_dbg_irreps
!-------------------------------------------------------------------------------
      character (len=4), allocatable :: rep_string(:)
      character (len=1)              :: parity
      integer                        ::  irrep,  jrrep,  ijrrep
      integer                        :: mi_val, mj_val, mij_val
      integer                        :: ifsym
!-------------------------------------------------------------------------------
 
!     allocate scratch memory       
      allocate(rep_string(num_dbg_irreps*4))

      if(double_group == 10)then

!       initialize MJ array
!       @@@@@@@@@@@@@@@@@@@
 
!       according to the group table in our reference we start with the 
!       boson irreps (integer) followed by the fermion (half-integer) irreps.

!       boson irreps
        irrep            = 0
        mj_val           = 0
        irrep            = irrep + 1
        mj2rep(mj_val,1) = irrep
        write(rep_string(irrep),'(i4)') mj_val

        do mj_val = 2,(num_dbg_irreps/2)-2,2
          irrep             = irrep + 1
          mj2rep( mj_val,1) = irrep
          write(rep_string(irrep),'(i4)')  mj_val
          irrep             = irrep + 1
          mj2rep(-mj_val,1) = irrep
          write(rep_string(irrep),'(i4)') -mj_val
        end do

        mj_val            = num_dbg_irreps/2
        irrep             = irrep + 1
        mj2rep( mj_val,1) = irrep
        mj2rep(-mj_val,1) = irrep
        write(rep_string(irrep),'(i4)')  mj_val

!       fermion irreps
        do mj_val = 1,(num_dbg_irreps/2)-1,2
          irrep             = irrep + 1
          mj2rep( mj_val,1) = irrep
          write(rep_string(irrep),'(i4)')  mj_val
          irrep             = irrep + 1
          mj2rep(-mj_val,1) = irrep
          write(rep_string(irrep),'(i4)') -mj_val
        end do

#ifdef LUCI_DEBUG
        do irrep = 1,num_dbg_irreps
          print *, ' *** mij-val and irrep',rep_string(irrep),irrep
        end do
#endif

!       double group multiplication table
!       ---------------------------------

!       =============================================================
!       product irrep = irrep_1 (x) irrep_2 = mj_irrep_1 + mj_irrep_2
!       =============================================================

!       due to the finite range it might be necessary 
!       to add or subtract, e.g. :
!       mj = 33/2 -> mj = -31/2 in the C32 subgroup 
!       <--> MJ = +/- 16 (32/2) transforms according to  the totally symmetric irrep in C32
!        
        do mi_val = -num_dbg_irreps/2, num_dbg_irreps/2, 1
          do mj_val = -num_dbg_irreps/2, num_dbg_irreps/2, 1
            mij_val = mi_val + mj_val 
!           map Cinf to lower Cx group (x_present = 32)
            if(mij_val  <  -num_dbg_irreps/2) mij_val = mij_val + num_dbg_irreps
            if(mij_val  >   num_dbg_irreps/2) mij_val = mij_val - num_dbg_irreps
            if(mij_val  == -num_dbg_irreps/2) mij_val = num_dbg_irreps/2
            irrep                     = mj2rep( mi_val,1)
            jrrep                     = mj2rep( mj_val,1)
            ijrrep                    = mj2rep(mij_val,1)
            dbg_mult_tab(irrep,jrrep) = ijrrep
          end do
        end do

      else ! double_group 11
!
!
!       initialize MJ array (we have inversion symmetry now)
!       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
!       according to the group table in our reference we start with the 
!       boson irreps (integer) followed by the fermion (half-integer) irreps 
!       while sub-dividing into "gerade" and "ungerade".

        irrep = 0

        do ifsym = 1, 2

          if(ifsym == 1)then
            parity = 'g'
          else
            parity = 'u'
          end if

!         boson irreps
          mj_val               = 0
          irrep                = irrep + 1
          mj2rep(mj_val,ifsym) = irrep
          write(rep_string(irrep),'(i3,a1)') mj_val, parity

          do mj_val = 2,(num_dbg_irreps/4)-2,2
            irrep                 = irrep + 1
            mj2rep( mj_val,ifsym) = irrep
            write(rep_string(irrep),'(i3,a1)')  mj_val, parity
            irrep                 = irrep + 1
            mj2rep(-mj_val,ifsym) = irrep
            write(rep_string(irrep),'(i3,a1)') -mj_val, parity
          end do

          mj_val                = num_dbg_irreps/4
          irrep                 = irrep + 1
          mj2rep( mj_val,ifsym) = irrep
          mj2rep(-mj_val,ifsym) = irrep
          write(rep_string(irrep),'(i3,a1)')  mj_val, parity
        end do

        do ifsym = 1, 2

          if(ifsym == 1)then
            parity = 'g'
          else
            parity = 'u'
          end if

!         fermion irreps
          do mj_val = 1,(num_dbg_irreps/4)-1,2
            irrep                 = irrep + 1
            mj2rep( mj_val,ifsym) = irrep
            write(rep_string(irrep),'(i3,a1)')  mj_val, parity
            irrep                 = irrep + 1
            mj2rep(-mj_val,ifsym) = irrep
            write(rep_string(irrep),'(i3,a1)') -mj_val, parity
          end do

        end do

#ifdef LUCI_DEBUG
        do irrep = 1,num_dbg_irreps
          print *, ' *** mij-val and irrep',rep_string(irrep),irrep
        end do
#endif

!       double group multiplication table
!       ---------------------------------

!       =============================================================
!       product irrep = irrep_1 (x) irrep_2 = mj_irrep_1 + mj_irrep_2
!       =============================================================

!       due to the finite range it might be necessary 
!       to add or subtract, e.g. :
!       mj = 17/2g(u) -> mj = -15/2g(u) in the C16h subgroup 
!       <--> MJ = +/- 8g(u) (16/2g(u)) transform according to the symmetric irreps (Ag/Au) in C16h
!        
        do mi_val = -num_dbg_irreps/4, num_dbg_irreps/4, 1
          do mj_val = -num_dbg_irreps/4, num_dbg_irreps/4, 1

!           mij value
            mij_val = mi_val + mj_val 

!           map Cinfh to lower Cxh group (x_present = 16)
            if(mij_val  <  -num_dbg_irreps/4) mij_val = mij_val + num_dbg_irreps/2
            if(mij_val  >   num_dbg_irreps/4) mij_val = mij_val - num_dbg_irreps/2
            if(mij_val  == -num_dbg_irreps/4) mij_val = num_dbg_irreps/4

!           g (x) g = g 
            irrep                     = mj2rep( mi_val,1)
            jrrep                     = mj2rep( mj_val,1)
            ijrrep                    = mj2rep(mij_val,1)
            dbg_mult_tab(irrep,jrrep) = ijrrep

!           u (x) g = u 
            irrep                     = mj2rep( mi_val,2)
            jrrep                     = mj2rep( mj_val,1)
            ijrrep                    = mj2rep(mij_val,2)
            dbg_mult_tab(irrep,jrrep) = ijrrep

!           g (x) u = u 
            irrep                     = mj2rep( mi_val,1)
            jrrep                     = mj2rep( mj_val,2)
            ijrrep                    = mj2rep(mij_val,2)
            dbg_mult_tab(irrep,jrrep) = ijrrep

!           u (x) u = u 
            irrep                     = mj2rep( mi_val,2)
            jrrep                     = mj2rep( mj_val,2)
            ijrrep                    = mj2rep(mij_val,1)
            dbg_mult_tab(irrep,jrrep) = ijrrep
          end do
        end do
 
      end if ! double_group (10 or 11)

!     release scratch matrix
      deallocate(rep_string)

   end subroutine

!*******************************************************************************
   subroutine convert_mj_ferm_2_dbg_irrep_linsym_krci(mj2rep,                  &
                                                      num_dbg_irreps,          &
                                                      active_mj_val,           &
                                                      active_ferm_sym,         &
                                                      active_dbg_irrep)
!-------------------------------------------------------------------------------
!
! purpose: return the active double group irrep for a given mj-value and fermion
!          irrep in approximate linear symmetry KR-CI/MCSCF calculations.
!          the present double groups are C32 and C16h, respectively.
!          this allows a maximum M_J value of:
!          --> +- 16   (C32;  double group 10)
!          --> +- 8g/u (C16h; double group 11)
!
!          the multiplication table follow the notation in
!          koster, dimmock, wheeler, statz: "properties of the 32 point groups".
!
!-------------------------------------------------------------------------------
      integer, intent(out) :: active_dbg_irrep
      integer, intent(in ) :: mj2rep(-num_dbg_irreps/2:num_dbg_irreps/2,2)
      integer, intent(in)  :: num_dbg_irreps
      integer, intent(in)  :: active_mj_val
      integer, intent(in)  :: active_ferm_sym
!-------------------------------------------------------------------------------

      print *, ' evaluating active dbg irrep...'

      print *, ' # dbg irreps               : ',num_dbg_irreps
      print *, ' # active MJ value (doubled): ',active_mj_val
      print *, ' # active fermion sym       : ',active_ferm_sym

      if(abs(active_mj_val) .gt. num_dbg_irreps/2)then
        call quit(' *** error in assigning the active dbg irrep: active MJ-value out of bounds. ***')
      end if

      active_dbg_irrep = mj2rep(active_mj_val,active_ferm_sym)

      print *, ' # active double group irrep: ',active_dbg_irrep

   end subroutine convert_mj_ferm_2_dbg_irrep_linsym_krci
!*******************************************************************************
   subroutine match_spinor_2_boson_irrep(ngsob1,                     &
                                         ngsob2,                     &
                                         boson_symm_array,           &
                                         ngsh,                       &
                                         orb_sym_vec,                &
                                         orb_off,                    &
                                         num_pos_sh,                 &
                                         num_inact_sh,               &
                                         mx_num_dbg_irreps,          &
                                         mx_num_gas_spaces,          &
                                         nr_fsym,                    &
                                         num_gas,                    &
                                         double_group,               &
                                         number_of_dbg_irreps,       &
                                         mj2rep,                     &
                                         imosp_dirac_counter1,       &
                                         imosp_dirac_counter2,       &
                                         imosp_dirac_mjub,           &
                                         imosp_dirac_mjb)

!*******************************************************************************
!
!     double group 9: spinfree d2h
!
!     --------------------------------------------------------------------------
      integer, intent(inout) :: ngsob1(mx_num_dbg_irreps,mx_num_gas_spaces)
      integer, intent(inout) :: ngsob2(mx_num_dbg_irreps,mx_num_gas_spaces)
      integer, intent(inout) :: boson_symm_array(*)
      integer, intent(in)    :: nr_fsym
      integer, intent(in)    :: num_gas
      integer, intent(out)   :: imosp_dirac_counter1(*)
      integer, intent(out)   :: imosp_dirac_counter2(*)
      integer, intent(out)   :: imosp_dirac_mjub(*)
      integer, intent(out)   :: imosp_dirac_mjb(*)
      integer, intent(in)    :: ngsh(2,mx_num_gas_spaces)
      integer, intent(in)    :: orb_sym_vec(*)
      integer, intent(in)    :: orb_off(nr_fsym)
      integer, intent(in)    :: num_pos_sh(nr_fsym)
      integer, intent(in)    :: num_inact_sh(nr_fsym)
      integer, intent(in)    :: mx_num_dbg_irreps
      integer, intent(in)    :: mx_num_gas_spaces
      integer, intent(in)    :: double_group
      integer, intent(in)    :: number_of_dbg_irreps
      integer, intent(in)    :: mj2rep(-number_of_dbg_irreps/2:number_of_dbg_irreps/2,2)
!     --------------------------------------------------------------------------
      integer                :: i, j, is_gas, is_shell
      integer                :: offset
      integer                :: mjval
      integer                :: local_boson_off
      integer                :: mat_address
      integer                :: mat_address2
      integer                :: varphi_1
      integer                :: varphi_2
      integer                :: spin_1
      integer                :: spin_2
      integer                :: orb_counter
      integer                :: deg_ferm    ! determine whether fermion irreps are degenerate (== 0) or not (== 1)
!     --------------------------------------------------------------------------
!#define LUCI_DEBUG
!     initialize
      call izero(ngsob1(1,1),mx_num_dbg_irreps*mx_num_gas_spaces)
      call izero(ngsob2(1,1),mx_num_dbg_irreps*mx_num_gas_spaces)

      deg_ferm    = 1
      if( mod(number_of_dbg_irreps/(2*nr_fsym),2) /= 0) deg_ferm = 0

!     initialize orbital counter
      orb_counter = 0

      do i = 1, nr_fsym

!       set offset in boson_symm_array / orb_sym_vec 
        j = orb_off(i) + num_pos_sh(i) + num_inact_sh(i)

!      miro: fix Index '0' of dimension 1 of array 'boson_symm_array' below lower bound of 1
        if(double_group /= 9 .and. j > 0)                      &
          boson_symm_array(j)      = (i-1)*number_of_dbg_irreps/4

        offset = 0

        do is_gas = 1, num_gas
!         print *, ' fsym, igas ==> # of spinors ',i,is_gas,ngsh(i,is_gas)
          do is_shell = 1, ngsh(i,is_gas)
    
            orb_counter = orb_counter + 1
            offset      = offset      + 1

            if(double_group == 9) local_boson_off = boson_symm_array(j+offset)

!
            if(double_group >= 10)then

!             a. read mj-value
              mjval = orb_sym_vec(j+offset)

!             b. place the spinors into the right position in the ngsobX arrays
!             kr_pair(varphi_1,varphi_2)
              varphi_1 =  mjval
              varphi_2 = -mjval
#ifdef LUCI_DEBUG
              print *,' mjval        = ',mjval
              print *,' i            = ',i
              print *,' mj2rep(..)   = ',mj2rep(varphi_1,i)
              print *,' mj2rep(..)2  = ',mj2rep(varphi_2,i)
#endif
              ngsob1(mj2rep(varphi_1,i),is_gas) = ngsob1(mj2rep(varphi_1,i),is_gas) + 1
              ngsob2(mj2rep(varphi_2,i),is_gas) = ngsob2(mj2rep(varphi_2,i),is_gas) + 1
!             store symmetry order for dirac-ordered mo spinors (e.g., ordered wrt energy)
              imosp_dirac_counter1(orb_counter) = mj2rep(varphi_1,i)
              imosp_dirac_counter2(orb_counter) = mj2rep(varphi_2,i)
              imosp_dirac_mjub(orb_counter)     = varphi_1
              imosp_dirac_mjb(orb_counter)      = varphi_2
            else ! non-linear double groups

!             mat_address = (i-1)*number_of_dbg_irreps/4 + (is_gas-1) * mx_num_dbg_irreps + mx_num_dbg_irreps/2 + 1
              select case(double_group)
!               C2h, Cs, C2
                case(4,5,6)
                  ngsob1((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 1,is_gas) =           & 
                  ngsob1((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 1,is_gas) + 1
                  ngsob2((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 2,is_gas) =           & 
                  ngsob2((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 2,is_gas) + 1
!               Ci, C1
                case(7,8)
                  ngsob1((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 1,is_gas) =           & 
                  ngsob1((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 1,is_gas) + 1
                  ngsob2((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 1,is_gas) =           & 
                  ngsob2((i-1)*number_of_dbg_irreps/4+number_of_dbg_irreps/2 + 1,is_gas) + 1
              end select

!             the offset calculation below has to be verified for Ci/C1
              imosp_dirac_counter1(orb_counter) = (i-1)*number_of_dbg_irreps/4 + number_of_dbg_irreps/2 + 1
              imosp_dirac_counter2(orb_counter) = (i-1)*number_of_dbg_irreps/4 + number_of_dbg_irreps/2 + 1 + 1*deg_ferm
            end if

          end do
        end do
      end do

#ifdef LUCI_DEBUG
      print *, 'printing NGSob1 in match_spinor_2_boson_irrep'
      do is_gas = 1, num_gas
        print *, 'is_gas = ',is_gas
        write(6,'(2x,32i3)') (ngsob1(j,is_gas),j=1,number_of_dbg_irreps)
      end do
      print *, 'printing NGSob2 in match_spinor_2_boson_irrep'
      do is_gas = 1, num_gas
        print *, 'is_gas = ',is_gas
        write(6,'(2x,32i3)') (ngsob2(j,is_gas),j=1,number_of_dbg_irreps)
      end do
      print *, 'printing original NGSH in match_spinor_2_boson_irrep'
      do is_gas = 1, num_gas
        print *, 'is_gas = ',is_gas
        write(6,'(2x,32i3)') (ngsh(j,is_gas),j=1,nr_fsym)
      end do
      print *, 'printing imosp_dirac_counterX in match_spinor_2_boson_irrep'
      do is_gas = 1, orb_counter
        write(6,'(2x,32i3)') imosp_dirac_counter1(is_gas)
        write(6,'(2x,32i3)') imosp_dirac_counter2(is_gas)
      end do
#endif

#undef LUCI_DEBUG
   end subroutine match_spinor_2_boson_irrep
!*******************************************************************************

   subroutine symmetry_setup_delete()
!-------------------------------------------------------------------------------
!
!  purpose: deallocate all symmetry-related arrays for KR-CI.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      if (allocated(  mj2rep)) call dealloc(mj2rep  )
      if (allocated(  invelm)) call dealloc(invelm  )
      if (allocated( iadjsym)) call dealloc(iadjsym )
      if (allocated(idbgmult)) call dealloc(idbgmult)

   end subroutine symmetry_setup_delete

end module
