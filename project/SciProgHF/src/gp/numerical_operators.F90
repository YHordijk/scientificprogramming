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

module numoper_integrals

! 
!     the following sets up the matrix elements for a so-called 
!     "numerical operator" (that is, when one has values of the given 
!     operator on a set of gridpoints ({P}), and contract
!     that with the value of the orbitals (or their derivatives, whatever)
!     and the integration weight to get the final integral: 
!
!     M_ab = < phi^n_a | operator | phi^m_b > = \sum^P_i w_i phi^n_a(i) operator(i) phi^m_b(i) 
!
!     where i is the i-th grid point, w_i the integration weight, and
!     phi^{n,m}_{a,b}(i) is the value of the {n,m}-th order derivative of
!     the (contracted) basis function {a,b} at point i.
!
!     information about the ao-basis blocks active for a given operator is
!     carried in the variable COMB, as indicated in a comment in dirac/dirone.F:
!
!     [...]
!     Chj     CLSCMB and PDOINT is now '+', '-', or '0',
!     Chj     corresponding to a factor of +1, -1, or 0 on that block.
!     C       The four blocks (AO basis): LL, SL, LS, SS
!     [...]
!
!     still to do is to have a general interface to numerical operators
!     and the properties section etc... but one example of such routines
!     can be found at 
!
!        embedding/emb_importPotential.F90 
!

   implicit none
   private
   public NumOper_OneElOpMatrix
 
   interface NumOper_OneElOpMatrix
      module procedure NumOper_OneElOp_SetMatrix_PS
   end interface

   integer, parameter :: max_boson_irreps = 8
   integer, save      :: n_boson_irreps = 1
   integer, save      :: n_ao
   logical, save      :: numoper_cb_initialized = .false.
   logical, save      :: debug = .false.

   type, private :: symmetry_basis_offsets
      integer  :: large_start(0:max_boson_irreps-1) 
      integer  :: small_start(0:max_boson_irreps-1) 
      integer  ::   large_end(0:max_boson_irreps-1) 
      integer  ::   small_end(0:max_boson_irreps-1)
   end type symmetry_basis_offsets 

   type(symmetry_basis_offsets), save :: basis_bs_offt, basis_bu_offt

   contains

      subroutine NumOper_BasisOffsetTables_print(offt)
         implicit none
         integer :: i
         type(symmetry_basis_offsets) :: offt

         do i = 0, n_boson_irreps-1
            write (*,*) 'For boson irrep',(i+1)
            write (*,*) '    Large cmp. basis starts at',offt%large_start(i),  &
                        ', ends at',offt%large_end(i),  &
                        '(',offt%large_end(i)-offt%large_start(i)+1,'elements)'
            write (*,*) '    Smale cmp. basis starts at',offt%small_start(i),  &
                        ', ends at',offt%small_end(i),  &
                        '(',offt%small_end(i)-offt%small_start(i)+1,'elements)'
         enddo

      end subroutine NumOper_BasisOffsetTables_print

      subroutine NumOper_BasisOffsetTables_init(offt)
         implicit none
         integer :: i
         type(symmetry_basis_offsets) :: offt

         do i = 0,(max_boson_irreps-1)
           offt%large_start(i) = 0
           offt%small_start(i) = 0

           offt%large_end(i) = 0
           offt%small_end(i) = 0
         enddo 

      end subroutine NumOper_BasisOffsetTables_init

      subroutine NumOper_CB_init
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "dgroup.h"
#include "symmet.h"

          integer :: i

          if (numoper_cb_initialized) return

          call NumOper_BasisOffsetTables_init(basis_bs_offt)
          call NumOper_BasisOffsetTables_init(basis_bu_offt)

          n_boson_irreps = nbsym
          n_ao           = ntbas(0)

! first collecting the boson unsorted (ie LS(irrep1) LS(irrep2) etc
! radovan: better naming for boson unsorted: hermit sorted

          i = 0 
          basis_bu_offt%large_start(i) = 1
          basis_bu_offt%large_end(i)   = nbbas(0,1)

          basis_bu_offt%small_start(i) = nbbas(0,1) + 1
          basis_bu_offt%small_end(i)   = nbbas(0,1) + nbbas(0,2)
           
          do i = 1, n_boson_irreps-1
               basis_bu_offt%large_start(i) = basis_bu_offt%small_end(i-1) + 1
               basis_bu_offt%large_end(i)   = basis_bu_offt%small_end(i-1) + nbbas(i,1)

               basis_bu_offt%small_start(i) = basis_bu_offt%large_end(i) + 1 
               basis_bu_offt%small_end(i)   = basis_bu_offt%large_end(i) + nbbas(i,2) 
          enddo

! radovan: boson unsorted are actually hermit offsets
          if (debug) then
             write (*,*) ' '
             write (*,*) '=============================='
             write (*,*) 'Boson-unsorted (DIRAC) offsets' 
             write (*,*) '=============================='
             write (*,*) ' '
             call NumOper_BasisOffsetTables_print(basis_bu_offt)
          endif

! and now the boson_sorted (ie: L(irrep1)...L(irrepN)S(irrep1)...S(irrepN)
! radovan: better naming for boson sorted: dirac sorted (and i think the examples are not right)

          do i = 0, n_boson_irreps-1
               basis_bs_offt%large_start(i) =  ibbas(i,1) + 1
               basis_bs_offt%large_end(i)   =  ibbas(i,1) + nbbas(i,1)

               basis_bs_offt%small_start(i) =  ibbas(i,2) + 1
               basis_bs_offt%small_end(i)   =  ibbas(i,2) + nbbas(i,2)
          enddo

! radovan: boson sorted are actually dirac offsets
          if (debug) then
             write (*,*) ' '
             write (*,*) '============================='
             write (*,*) 'Boson-sorted (HERMIT) offsets' 
             write (*,*) '============================='
             write (*,*) ' '
             call NumOper_BasisOffsetTables_print(basis_bs_offt)
          endif
!
! aspg, debugging X2C runs... perhapes it needs to be initialized again and again...
!         numoper_cb_initialized = .true.

      end subroutine NumOper_CB_init

!
! NumOper_OneElOp_SetMatrix_PS
!
! calculates the matrix elements of a one-electron operator given over a
! grid (found in file unit "fromfile", which is assumed to be opened prior 
! to entering this routine.
!
! operates on a packed storage (_PS_) matrix
!

      subroutine NumOper_OneElOp_SetMatrix_PS(oneint,from_file,doint,irrep_oper)
         use memory_allocator
         use interface_ao_specific
         use dirac_ao_eval

         implicit none
         logical :: doint(4)
         integer :: from_file, ipoint, npoint, ndim, intclass
         integer :: a_start, a_end, b_start, b_shift, b_end
         integer :: irrep_oper       ! the boson irrep of the original operator
         real(kind=8) :: oneint(*)
         real(kind=8), allocatable :: ao(:), buffer(:)
         real(kind=8), pointer     :: grid(:,:), operator(:)
         real(kind=8) :: wt_oper
         integer :: irrep, irrep1, irrep2, k, i
         integer, allocatable :: offt_bu(:)

         if (.not.numoper_cb_initialized) call NumOper_CB_init()


         rewind(from_file)
         read (from_file,*) npoint

         if (associated(grid)) nullify(grid)
         call allocp(grid,4,npoint)
         if (associated(operator)) nullify(operator)
         call allocp(operator,npoint)

         do i = 1, npoint
            read(from_file,*) grid(:,i),operator(i)
         enddo

         call interface_ao_write()
         call dirac_ao_eval_init(0, 0, .false.)

         ndim = n_ao
         call alloc(ao,ndim)
         call alloc(offt_bu,ndim)
         call alloc(buffer,ndim)
 
         call resort_so(offt_bu)

! start numerical integration...

         do ipoint = 1, npoint
            call get_ao(1,                  &
                        (/grid(1,ipoint)/), &
                        (/grid(2,ipoint)/), &
                        (/grid(3,ipoint)/), &
                        ao, buffer)

            wt_oper = operator(ipoint)*grid(4,ipoint)            

            do irrep1 = 0,n_boson_irreps-1  ! irrep of the bra 

               irrep2 = ieor(irrep_oper,irrep1) ! symmetry of the product bra * operator gives ket

               do intclass = 1, 4
                  if (doint(intclass)) then 

                     if (intclass.eq.1) then ! LL block
                        a_start = basis_bu_offt%large_start(irrep1)
                        a_end   = basis_bu_offt%large_end(irrep1)

                        b_shift = 0
                        b_start = basis_bu_offt%large_start(irrep2)
                        b_end   = basis_bu_offt%large_end(irrep2)

                     elseif (intclass.eq.2) then ! SL block

                        write (*,*) 'Construction of numerical 1e-operators for '//  &
                                    'the SL block (packed storage) not implemented.'
                        call quit('Numerical 1e-operator cannot be created.')  

                     elseif (intclass.eq.3) then ! LS block
                        a_start = basis_bu_offt%small_start(irrep1)
                        a_end   = basis_bu_offt%small_end(irrep1) 

                        b_shift = basis_bu_offt%large_end(irrep2)
                        b_end   = basis_bu_offt%small_end(irrep2)

                     elseif (intclass.eq.4) then ! SS block
                        a_start = basis_bu_offt%small_start(irrep1)
                        a_end   = basis_bu_offt%small_end(irrep1)

                        b_shift = 0
                        b_start = basis_bu_offt%small_start(irrep2) 
                        b_end   = basis_bu_offt%small_end(irrep2) 

                     endif

                     call NumOper_OneElOp_PS_MatrixUpdateAtGridPoint(oneint,ao,wt_oper,a_start,a_end,b_shift,b_end,offt_bu)

                  endif ! do this class of integrals?
               enddo ! intclass
            enddo ! irrep1
         enddo ! gridpoints

         call deallocp(grid)
         call deallocp(operator)
         call dealloc(ao)
         call dealloc(offt_bu)
         call dealloc(buffer)

      end subroutine NumOper_OneElOp_SetMatrix_PS

!
! NumOper_OneElOp_PS_MatrixUpdateAtGridPoint
!
! updates all matrix elements for a given integral class (LL, LS, SL or SS), by
! adding the values for a new point in the numerical integration.
!
! works for packed storage (_PS_) matrices, operating on the lower diagonal.
!
! the a_start, a_end, b_shift and b_end variables indicate the ranges of indexes
! related to a given integral class, and are currently set outsite... they
! may be perhaps moved inside? thoughts, anyone?

      subroutine NumOper_OneElOp_PS_MatrixUpdateAtGridPoint(oneint,   &
                                                            ao,       & 
                                                            operator, & 
                                                            a_start,  &
                                                            a_end,    &
                                                            b_shift,  &
                                                            b_end,    & 
                                                            offt_bu)

         integer, intent(in) :: a_start, a_end, b_shift,b_end, offt_bu(:)
         integer :: a, b, ab, irrep
         real(kind=8), intent(in) :: ao(*), operator
         real(kind=8) :: oneint(*)

         do a = a_start, a_end
            do b = a+b_shift, b_end
               ab =  a + b*(b-1)/2
               oneint(ab) = oneint(ab) + ao(offt_bu(a))*ao(offt_bu(b))*operator 
            enddo 
         enddo 

      end subroutine NumOper_OneElOp_PS_MatrixUpdateAtGridPoint

      subroutine resort_so(offset_bs_to_bu)
         integer :: i, j, l, k, offset_bs_to_bu(:)

         do l = 0, n_boson_irreps-1 
            k = 0
            do i = basis_bs_offt%large_start(l), basis_bs_offt%large_end(l)
               j = basis_bu_offt%large_start(l) + k
               offset_bs_to_bu(j) = i
               k = k + 1
            enddo
         enddo

         do l = 0, n_boson_irreps-1 
            k = 0
            do i = basis_bs_offt%small_start(l), basis_bs_offt%small_end(l)
               j = basis_bu_offt%small_start(l) + k
               offset_bs_to_bu(j) = i
               k = k + 1
            enddo
         enddo

      end subroutine resort_so

end module numoper_integrals
