#ifdef MOD_LAO_REARRANGED
module london_utils

  use memory_allocator 
  use dirac_cfg
  use num_grid_gen
  use dirac_interface

  implicit none

  public check_symmetry 
  public combine_kappa_T
  public get_m_ai 
  public get_m_ij 
  public get_m_ee 
  public get_m_ep 
  public test_new_routines
  save

  private

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "aovec.h"
#include "blocks.h"
#include "dcbfir.h"
#include "dcbham.h"
#include "iratdef.h"
#include "cbihr2.h"
#include "dcbdhf.h"
#include "shells.h"
#include "dcbxlr.h"

#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
  integer istat(mpi_status_size)
#endif /* ifdef VAR_MPI */



  integer, parameter :: max_nr_xyz_comp = 3
  integer            ::     nr_xyz_comp

  integer, parameter :: tbmo_unit = 71
  logical            :: debug_me


contains

  subroutine get_m_ai(m, ee, ep)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     logical, intent(in) :: ee, ep
     real(8), allocatable :: m_cp(:,:,:)
     integer, allocatable :: indices_virt(:)
     integer, allocatable :: indices_occ(:)
     integer :: isym, iz, i, j, ip
!    ---------------------------------------------------------------------------

     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)

     do isym = 1, nfsym

!       find absolute indices of occuppied orbitals
!       -------------------------------------------
        call alloc(indices_occ, nocc(isym))
        indices_occ = 0
        call get_orbital_indices_from_strings(xlr_indstr(1, isym),         &
                                              indices_occ,                 &
                                              .true., .false., .false.,    &
                                              isym)

!       find absolute indices of virtual orbitals
!       -----------------------------------------
        call alloc(indices_virt, norb(isym) - nocc(isym))
        indices_virt = 0
        call get_orbital_indices_from_strings(xlr_indstr(3, isym),     &
                                              indices_virt,            &
                                              .false., ee, ep,         &
                                              isym)


!       copy only needed elements of matrix m to helper matrix m_cp
!       -----------------------------------------------------------
        do iz = 1, nz
           do j = 1, norb(isym) - nocc(isym)
              do i = 1, nocc(isym)
                 if (indices_occ(i) .ne. 0 .and. indices_virt(j) .ne. 0) then
!                   occ-virt
                    m_cp(indices_occ(i), indices_virt(j), iz) = &
                       m(indices_occ(i), indices_virt(j), iz)
!                   virt-occ
                    m_cp(indices_virt(j), indices_occ(i), iz) = &
                       m(indices_virt(j), indices_occ(i), iz)
                 end if
              end do
           end do
        end do

        call dealloc(indices_occ)
        call dealloc(indices_virt)
     end do


     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ai, m full:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
        write(*, *) 'get_m_ai, m_cp, ai only:'
        call prqmat(m_cp, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if

!    clean
!    -----
     call dzero(m, n2orbxq)
     call dcopy(n2orbxq, m_cp, 1, m, 1)
     call dealloc(m_cp)


  end subroutine


  subroutine get_m_ij(m)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     real(8), allocatable :: m_cp(:,:,:)

     integer :: isym, iz, i, j
     integer, allocatable :: indices_occ(:)
!    ---------------------------------------------------------------------------

     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)

     do isym = 1, nfsym

!       find absolute indices of occuppied orbitals
!       -------------------------------------------
        call alloc(indices_occ, nocc(isym))
        indices_occ = 0
        call get_orbital_indices_from_strings(xlr_indstr(1, isym),         &
                                              indices_occ,                 &
                                              .true., .false., .false.,    &
                                              isym)


!       copy only needed elements of matrix m to helper matrix m_cp
!       -----------------------------------------------------------
        do iz = 1, nz
           do j = 1, nocc(isym)
              do i = 1, nocc(isym)
                 if (indices_occ(i) .ne. 0 .and. indices_occ(j) .ne. 0) then
                    m_cp(indices_occ(i), indices_occ(j), iz) = &
                       m(indices_occ(i), indices_occ(j), iz)
                 end if
              end do
           end do
        end do

        call dealloc(indices_occ)
     end do

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ij, m full:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
        write(*, *) 'get_m_ij, m_cp, ij only:'
        call prqmat(m_cp, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if


!    clean
!    -----
     call dzero(m, n2orbxq)
     call dcopy(n2orbxq, m_cp, 1, m, 1)
     call dealloc(m_cp)

  end subroutine


  subroutine combine_kappa_T(kappa, length, nbtyp, lab1, ih)
!    ---------------------------------------------------------------------------
     real(8), intent(inout)   :: kappa(length, nz)
     integer, intent(in)      :: length, nbtyp, ih
     character*16, intent(in) :: lab1

     real(8), allocatable   :: T(:, :, :)
     real(8), allocatable   :: T_ai_vec(:, :)
     real(8), allocatable   :: kappa_matrix(:)
     integer, allocatable   :: orbital_rotation_indices(:, :)
     integer :: which_icomp, chosen_icomp, irep, iz, i, j
     logical :: sum_ep_rn, ep, ee
     real(8) :: h
     integer :: iht, iq, iqsym(nz), iih

!    this subroutine does:
!    kappa_{ai}^B := kappa_{ai}^B - T_{ai}^B

!    where: kappa_{ai}^B - elements of solution vector
!           T_{ai}^B     - elements of connection matrix (occ-virt blocks)
!                          the 'full' connection matrix is read from TBMO file

!    it is used to calculate the NMR shielding tensor
!    from Eq. 40-42 from JCP 131, 124119
!    ---------------------------------------------------------------------------

!    which component of magnetic field is that:
!    ------------------------------------------
     if (lab1 == 'LAO-XRM1H1') chosen_icomp = 1
     if (lab1 == 'LAO-YRM1H1') chosen_icomp = 2
     if (lab1 == 'LAO-ZRM1H1') chosen_icomp = 3

!    prepare info about solution vector:
!    -----------------------------------
     ep = .false.
     ee = .false.
     if (nbtyp == 1) then
       ee = .true.
       ep = .false.
     else if (nbtyp == 2) then
       ep = .true.
       ee = .false.
     else
        call quit('error in combine_kappa_T')
     end if

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'combine_kappa_T: solution vector for operator = ', lab1
        write(*, *) 'which class of rotations: ee, ep = ', ee, ep
        write(*, *) 'length, nz, ih = ', length, nz, ih
        write(*, *) 'entering kappa:'
        call prbvec(lupri, kappa , 1, length)
     end if

!    get information about orbital rotation indices
!    ----------------------------------------------
     call alloc(orbital_rotation_indices, 2, length)
     call fill_orbital_rotation_indices_vector(ee,     &
                                               ep,     &
                                               length, & 
                                               orbital_rotation_indices)
     

!    read connection matrix from TBMO file:
!    --------------------------------------
     call alloc(T, norbt, norbt, nz)

     open(tbmo_unit,                 &
          file   = 'TBMO',           &
          status = 'old',            &
          form   = 'unformatted',    &
          access = 'sequential',     &
          position = 'rewind',       &
          action = 'read')
    
1    continue
    
     read(tbmo_unit) nr_xyz_comp
     read(tbmo_unit) which_icomp
     read(tbmo_unit) irep
     read(tbmo_unit) sum_ep_rn
     read(tbmo_unit) T

     if (nr_xyz_comp > 1) then
!       when calculating full shielding tensor
!       there will be three matrices on TBMO file: T^Bx, T^By and T^Bz
!       here we need to read the requested one:
        if (chosen_icomp == 2) then
           if (which_icomp == 2) then
              go to 2
           else
              go to 1
           end if
        else if (chosen_icomp == 3) then
           if (which_icomp == 3) then
              go to 2
           else
              go to 1
           end if
        end if
     end if
     
2    continue

     close(tbmo_unit, status = 'keep')


!    leave only occ-virt part of T:
!    ------------------------------
!
     debug_me = .false.
     if (debug_me) then
        write(*, *) 'combine_kappa_t: full T, ee, ep = ', ee, ep
        call prqmat(T, norbt, norbt, norbt, norbt, nz,  &
                    ipqtoq(1,irep), lupri)
     end if

     call get_m_ai(T, ee, ep)

     debug_me = .true.
     if (debug_me) then
        write(*, *) 'combine_kappa_t: T_ai,ee, ep = ', ee, ep
        call prqmat(T, norbt, norbt, norbt, norbt, nz,  &
                    ipqtoq(1,irep), lupri)
     end if

!    create a vector, T_ai_vec, from matrix T
!    the order of indices on T_ai_vec is the same as on kappa vector:
!    gosia fixme: it works now for symcon, but we should be careful, because:
!    *) if symcon then T is antihermitian 
!    *) if natural con. then T is a general matrix
!    and we take only "virt-occ" elements to a vector
!    we want T to "follow" kappas, but in case of natural connection now we 
!    neglect all "occ-virt" (we take only 'virt-occ")
!    i think we should treat 'kappa_ia - T_ia' and 'kappa_ai - T_ai' separately

     call alloc(T_ai_vec, length,nz)
     call dzero(T_ai_vec, length*nz)
     call matrix_to_vector(T_ai_vec, length, T, irep, ih, orbital_rotation_indices)
     debug_me = .false.
     if (debug_me) then
       write(*, *) 'combine_kappa_t: T_ai_vec, length,ee,ep = ', length,ee, ep
       call prqmat(T_ai_vec(1,1), length, 1, length, 1, nz,  &
                   ipqtoq(1,irep), lupri)
     end if

!    substract T from kappa:
     !kappa = kappa - T_ai_vec
     call daxpy(length*nz, -1.0d0, T_ai_vec, 1, kappa, 1)

     call dealloc(T_ai_vec)
     call dealloc(T)
     call dealloc(orbital_rotation_indices)

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'combine_kappa_t: kappa on output, ee, ep = ', ee, ep
        call prqmat(kappa(1,1), length, 1, length, 1, nz,  &
                    ipqtoq(1,irep), lupri)
     end if

  end subroutine


  subroutine matrix_to_vector(vector, length, matrix, irep, ih, orbital_rotation_indices)
     real(8), intent(inout) :: vector(length, nz)
     real(8), intent(inout) :: matrix(norbt, norbt, nz)
     integer, intent(in)    :: length, irep, orbital_rotation_indices(2, length), ih
     integer :: iz, i, s, is

     do iz = 1, nz
        do is = 1, length
           i = orbital_rotation_indices(1, is)
           s = orbital_rotation_indices(2, is)
           vector(is, iz) = matrix(s, i, iz)
        end do
     end do
   
  end subroutine


  subroutine fill_orbital_rotation_indices_vector(ee,     &
                                                  ep,     &
                                                  length, &
                                                  orbital_rotation_indices)
!    ---------------------------------------------------------------------------
     logical, intent(in)    :: ee, ep
     integer, intent(in)    :: length
     integer, intent(inout) :: orbital_rotation_indices(2, *)
     integer :: i1, i, j, isym, ioff
     integer :: ip(2), ap(2), ie(2), ae(2)
     integer, allocatable   :: temp_occ(:), temp_virt(:) 
     integer :: ind, nr_occ_used(2), nr_virt_used(2)
!    ---------------------------------------------------------------------------

     orbital_rotation_indices(1:2, 1:length) = 0

     do isym = 1, nfsym
!       occuppied indices:
!       ------------------
        call alloc(temp_occ, nocc(isym))
        temp_occ = 0
        call get_orbital_indices_from_strings(xlr_indstr(1, isym),         &
                                              temp_occ,                    &
                                              .true., .false., .false.,    &
                                              isym)

        nr_occ_used(isym) = 0
        do i = 1, nocc(isym)
           if (temp_occ(i) .ne. 0) nr_occ_used(isym) = nr_occ_used(isym) + 1
        end do


!       virtual indices:
!       ----------------
        call alloc(temp_virt, norb(isym) - nocc(isym))
        temp_virt = 0
        call get_orbital_indices_from_strings(xlr_indstr(3, isym),     &
                                              temp_virt,               &
                                              .false., ee, ep,         &
                                              isym)

        nr_virt_used(isym) = 0
        do i = 1, norb(isym) - nocc(isym)
           if (temp_virt(i) .ne. 0) nr_virt_used(isym) = nr_virt_used(isym) + 1
        end do

!       combine occuppied-virtual:
!       --------------------------
!       orbital_rotation_indices(1, i) = occ1,      ...      occ1, occ2,       ...      occ2, etc.
!                                        |- 'length_virt' times -| |- 'length_virt' times -|                    
!       
!       orbital_rotation_indices(1, i) = virt1, virt2, ..., virtn,     virt1, virt2, ..., virtn, etc.
!                                        |-      occ1           -|     |-          occ2       -| 


        do i = 1, nr_occ_used(isym)
           do j = 1, nr_virt_used(isym)
              if (temp_occ(i) .ne. 0 .and. temp_virt(j) .ne. 0) then
                ! miro: fix due to runtimecheck
                 if (isym.gt.1) then
                     ioff = (isym - 1)*nr_occ_used(isym - 1)*nr_virt_used(isym - 1) 
                 else
                     ioff = 0
                 endif
                 ind = ioff + (i-1)*nr_virt_used(isym) + j
                 orbital_rotation_indices(1, ind) = temp_occ(i)
                 orbital_rotation_indices(2, ind) = temp_virt(j)
              end if
           end do
        end do

        call dealloc(temp_occ)
        call dealloc(temp_virt)
     end do

     debug_me = .true.
     if (debug_me) then
        do i = 1, length
           write(*, *) 'i, orbital_rotation_indices(1, i), orbital_rotation_indices(2, i) = ', &
                        i, orbital_rotation_indices(1, i), orbital_rotation_indices(2, i) 
        end do
     end if

  end subroutine


  subroutine check_symmetry(mat, ndim, nz)
    integer, intent(in) :: ndim, nz
    real(8), intent(in) :: mat(ndim, ndim, *)
    real(8) :: a
    integer :: i, j, iz


    a = 0.0d0
    do iz = 1, nz
      do i = 1, ndim - 1
        do j = i + 1, ndim
          a = a + (dabs(mat(i, j, iz)) - dabs(mat(j, i, iz)))
        end do
      end do
    end do
    a = a/dfloat(((ndim*ndim - ndim)/2)*nz)
    write(*, *) 'symmetry mat: ', a

  end subroutine


  subroutine get_orbital_indices_from_strings(indices_string, indices_vect, &
                                              occ, e, p, isym)
!    ---------------------------------------------------------------------------
     character*72, intent(inout) :: indices_string(3, 2)
     integer, intent(inout) :: indices_vect(*)
     logical, intent(in)    :: occ, e, p
     integer, intent(in)    :: isym
     integer :: ifrp, ioff, i, ip, nvec, nvirt, ind1, ind2, ishift
!    ---------------------------------------------------------------------------


     if (occ) then

!       occuppied orbitals indices
!       --------------------------

        nvec = 1
        call numlst(indices_string, indices_vect,     &
                    nocc(isym), 1, nocc(isym), isym, nvec)

!       we need absolute indices:
        i = 0
        ioff = npsh(isym) + iorb(isym)
        do while (i < nocc(isym))
           i = i + 1
           if (indices_vect(i) .ne. 0) then
             indices_vect(i) = ioff + indices_vect(i)
           end if
        end do

        debug_me = .false.
        if (debug_me) then
           do i = 1, nocc(isym)
              write(*, *) 'OCC: isym, iorb(isym), i, indices_vect(i) = ', &
                                isym, iorb(isym), i, indices_vect(i)
           end do
        end if

     else

!       virtual orbitals indices
!       ------------------------

        nvirt = norb(isym) - nocc(isym)

        call numlst(indices_string, indices_vect,     &
                    nvirt, -npsh(isym), nesh(isym), isym, nvec)

        ioff = 0
        ip = 0

!       indices_vect contains ALL virtual indices, now, we need to decide what we keep:
!       *)either only positronic indices (p = true, e = false)
!       *)or only electronic indices (p = false, e = true)
!       *)or both (p = true, e = true)

        if (p .and. e) then
           do i = 1, nvirt
              if (indices_vect(i) .lt. 0) then
!                check for positronic virtual indices
                 ioff = npsh(isym) + iorb(isym) + 1
                 indices_vect(i) = ioff + indices_vect(i)
              else if (indices_vect(i) .gt. nocc(isym)) then
!                check for electronic virtual indices
                 ioff = npsh(isym) + iorb(isym)
                 indices_vect(i) = ioff + indices_vect(i)
              end if
           end do

        else if (p .and. .not. e) then
           do i = 1, nvirt
              if (indices_vect(i) .lt. 0) then
!                check for positronic virtual indices
                 ioff = npsh(isym) + iorb(isym) + 1
                 indices_vect(i) = ioff + indices_vect(i)
              else if (indices_vect(i) .gt. nocc(isym)) then
!                if there are any electronic indices, we don't want them here, so put to 0
                 indices_vect(i) = 0
              end if
           end do

        else if (e .and. .not. p) then
           ishift = 0
           do i = 1, nvirt
              if (indices_vect(i) .lt. 0) then
!                check for positronic virtual indices
!                we don't want them here, so just check how many there are (kept in "ishift")
                 ishift = ishift + 1
                 indices_vect(i) = 0
              end if
           end do
           do i = ishift+1, nvirt
              if (indices_vect(i) .ne. 0) then
              !if (indices_vect(i) .gt. nocc(isym)) then
!                the electronic indices we want to write from indices_vect(1:), so:
                 ioff = npsh(isym) + iorb(isym)
                 indices_vect(i-ishift) = ioff + indices_vect(i)
!                zero-out remainings:
                 if (ishift .gt. 0) indices_vect(i) = 0
              end if
           end do
           if (ishift .gt. 0) then
              indices_vect((nvirt - ishift + 1):nvirt) = 0
           end if 

        end if

        debug_me = .false.
        if (debug_me) then
           do i = 1, norb(isym) - nocc(isym)
              write(*, *) 'VIRT: isym, i, indices_vect(i) = ', &
                                 isym, i, indices_vect(i)
           end do
        end if

     end if

  end subroutine

! subroutine copy-pasted from visual.F90:
  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, nz)
    real(8), intent(out) :: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              ¦
!              ¦  column
!              ¦  ¦
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                          - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine



  subroutine get_m_ee(m)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     real(8), allocatable   :: m_cp(:,:,:)
     integer, allocatable :: indices_virt(:)
     integer, allocatable :: indices_occ(:)
     integer :: isym, iz, i, j, ip
! get occ-"virtual-electronic" blocks of matrix m
!    ---------------------------------------------------------------------------

     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)

     do isym = 1, nfsym

!       find absolute indices of occuppied orbitals
!       -------------------------------------------
        call alloc(indices_occ, nocc(isym))
        indices_occ = 0
        call get_orbital_indices_from_strings(xlr_indstr(1, isym),         &
                                              indices_occ,                 &
                                              .true., .false., .false.,    &
                                              isym)

!       find absolute indices of virtual orbitals
!       -----------------------------------------
        call alloc(indices_virt, norb(isym) - nocc(isym))
        indices_virt = 0
        call get_orbital_indices_from_strings(xlr_indstr(3, isym),     &
                                              indices_virt,            &
                                              .false., .true., .false.,         &
                                              isym)


!       copy only needed elements of matrix m to helper matrix m_cp
!       -----------------------------------------------------------
        do iz = 1, nz
           do j = 1, norb(isym) - nocc(isym)
              do i = 1, nocc(isym)
                 if (indices_occ(i) .ne. 0 .and. indices_virt(j) .ne. 0) then
!                   occ-'virt-e'
                    m_cp(indices_occ(i), indices_virt(j), iz) = &
                       m(indices_occ(i), indices_virt(j), iz)
!                   'virt-e'-occ
                    m_cp(indices_virt(j), indices_occ(i), iz) = &
                       m(indices_virt(j), indices_occ(i), iz)
                 end if
              end do
           end do
        end do

        call dealloc(indices_occ)
        call dealloc(indices_virt)
     end do


     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ee, m full:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
        write(*, *) 'get_m_ee, m_cp, ee only:'
        call prqmat(m_cp, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if


!    clean
!    -----
     call dzero(m, n2orbxq)
     call dcopy(n2orbxq, m_cp, 1, m, 1)
     call dealloc(m_cp)

  end subroutine

  subroutine get_m_ep(m)
!    ---------------------------------------------------------------------------
     real(8), intent(inout) :: m(norbt, norbt, nz)
     real(8), allocatable   :: m_cp(:,:,:)
     integer, allocatable :: indices_virt(:)
     integer, allocatable :: indices_occ(:)
     integer :: isym, iz, i, j, ip
!    ---------------------------------------------------------------------------

     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)

     do isym = 1, nfsym

!       find absolute indices of occuppied orbitals
!       -------------------------------------------
        call alloc(indices_occ, nocc(isym))
        indices_occ = 0
        call get_orbital_indices_from_strings(xlr_indstr(1, isym),         &
                                              indices_occ,                 &
                                              .true., .false., .false.,    &
                                              isym)

!       find absolute indices of virtual orbitals
!       -----------------------------------------
        call alloc(indices_virt, norb(isym) - nocc(isym))
        indices_virt = 0
        call get_orbital_indices_from_strings(xlr_indstr(3, isym),     &
                                              indices_virt,            &
                                              .false., .false., .true.,         &
                                              isym)


!       copy only needed elements of matrix m to helper matrix m_cp
!       -----------------------------------------------------------
        do iz = 1, nz
           do j = 1, norb(isym) - nocc(isym)
              do i = 1, nocc(isym)
                 if (indices_occ(i) .ne. 0 .and. indices_virt(j) .ne. 0) then
!                   occ-'virt-p'
                    m_cp(indices_occ(i), indices_virt(j), iz) = &
                       m(indices_occ(i), indices_virt(j), iz)
!                   'virt-p'-occ
                    m_cp(indices_virt(j), indices_occ(i), iz) = &
                       m(indices_virt(j), indices_occ(i), iz)
                 end if
              end do
           end do
        end do

        call dealloc(indices_occ)
        call dealloc(indices_virt)
     end do

     debug_me = .false.
     if (debug_me) then
        write(*, *) 'get_m_ep, m full:'
        call prqmat(m, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
        write(*, *) 'get_m_ep, m_cp, ep only:'
        call prqmat(m_cp, &
                    norbt, norbt, norbt, norbt, &
                    nz, ipqtoq(1, 0), lupri)
     end if


!    clean
!    -----
     call dzero(m, n2orbxq)
     call dcopy(n2orbxq, m_cp, 1, m, 1)
     call dealloc(m_cp)


  end subroutine


  subroutine test_new_routines(m)
     real(8), intent(inout) :: m(norbt, norbt, nz)
     real(8), allocatable   :: m_cp(:,:,:)
     logical :: ee, ep

! test "ij" and "ai"
!     ee = .true.
!     ep = .true.
!     if (xlr_skipep .or. .not. doeprn) ep = .false.
!     if (xlr_skipee) ee = .false.
!
!     call alloc(m_cp, norbt, norbt, nz)
!     call dzero(m_cp, n2orbxq)
!     call dcopy(n2orbxq, m, 1, m_cp, 1)
!
!     call get_m_ij(m_cp)
!     call daxpy(n2orbxq, 1.0d0, m_cp, 1, m, 1)
!
!     call get_m_ai(m, ee, ep)
!     call daxpy(n2orbxq, 1.0d0, m_cp, 1, m, 1)
!
!     call dealloc(m_cp)     

! test "ee" and "ep":

     call alloc(m_cp, norbt, norbt, nz)
     call dzero(m_cp, n2orbxq)
     call dcopy(n2orbxq, m, 1, m_cp, 1)

     call get_m_ee(m_cp)
     call daxpy(n2orbxq, 1.0d0, m_cp, 1, m, 1)

     call get_m_ep(m)
     call daxpy(n2orbxq, 1.0d0, m_cp, 1, m, 1)

     call get_m_ij(m)
     call daxpy(n2orbxq, 1.0d0, m_cp, 1, m, 1)

     call dealloc(m_cp)     



  end subroutine

end module
#endif
