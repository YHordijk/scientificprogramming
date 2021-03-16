  module projectors
  use relcc_cfg
  implicit none

    private

    public create_projector
    public create_t_projector

! the excited determinant derived type holds the information on occupied/virtual spinor energies, irreps and position within
! the irreps, for the spinors making up a given excited state determinant
!
! what's still to be done, if found necessary/useful, is to have functions which inquire the variable of this type and provide
! information on this excited determinant: whether it is a 1h/1p/2h1p/2p1h/1p1h/2p2h etc (which one can get by the dimensions
! of the "occupied_energies" and "virtuals_energies" variables), the irrep to which it belongs to etc
   public excited_determinant
   type excited_determinant
      real*8,  allocatable :: occupied_energies(:)
      integer, allocatable :: occupied_irreps(:)
      integer, allocatable :: occupied_indexes_in_irreps(:)

      real*8,  allocatable :: virtuals_energies(:)
      integer, allocatable :: virtuals_irreps(:)
      integer, allocatable :: virtuals_indexes_in_irreps(:)
   end type excited_determinant

! the mapping derived type is meant to collect the different determinants into an array, and with that the position of the
! determinant on the hamiltonian would be provided by its index in the "det" variable. having a derived type has the advantage of
! being able to keep track of multiple mappings if needed.
   public mapping
   type mapping
       type(excited_determinant), allocatable :: det(:)
   endtype mapping
 

   contains

! function to create the operator to project out certain lines of the trial vectors in the
! davidson procedure
! it returs an array of logical variables, where the entries are set to false if a given line is
! to be suppressed. we determine these lines by creating a mapping between the composition of the 
! excited determinants in terms of the orbital energies and and their position on the similarity
! transformed hamiltonian with a call to mapping_spinors_configurations(), which works quite similarly
! to the routines to print out the solution vectors.
function create_projector(eom_type, ndet, state_sym, eps)
#include "files.inc"
#include "inpt.inc"
   logical, pointer :: create_projector(:) 

   integer, intent(in) :: eom_type
   integer, intent(in) :: ndet
   integer, intent(in) :: state_sym
   real(kind=8), intent(in) :: eps(*)

   type(mapping) :: det_map
   integer :: i, j, ndim_occ, ndim_virt
   logical :: spinors_in_occupied_window(2)
   logical :: spinors_in_virtuals_window(2)

!   if (.not.(relcc_projectors_do_core_valence_separation.or.relcc_projectors_do_restricted_excitation_window)) then
!      create_projector => NULL()
!      return
!
!   else
      call mapping_spinors_configurations (det_map, eom_type, ndet, state_sym, eps)

      allocate(create_projector(ndet))
      create_projector = .false.

      spinors_in_occupied_window = .true.
      spinors_in_virtuals_window = .true. 

! we loop over the composition of each determinant in the "excitation" manifold 
      do i = 1, ndet
         if (allocated(det_map%det(i)%occupied_energies)) then 
            ndim_occ = size(det_map%det(i)%occupied_energies,1)
            do j = 1, ndim_occ
               if ( (det_map%det(i)%occupied_energies(j).le.relcc_projectors_rew_occ_max_energy).and. &
                  & (det_map%det(i)%occupied_energies(j).ge.relcc_projectors_rew_occ_min_energy) ) then
                  spinors_in_occupied_window(j) = .true.
               else
                  spinors_in_occupied_window(j) = .false.
               end if
            end do
            if (ndim_occ.eq.1) spinors_in_occupied_window(2) = spinors_in_occupied_window(1) 

            if ((ndim_occ.eq.2).and.relcc_projectors_rew_remove_double_occupied) then
! if we require to remove doubly excited determinants with both indices in the
! occupied window, force these out
                if (spinors_in_occupied_window(1).and.spinors_in_occupied_window(2)) then
                    spinors_in_occupied_window(1) = .false. 
                    spinors_in_occupied_window(2) = .false. 
                end if
            end if
         end if

         if (allocated(det_map%det(i)%virtuals_energies)) then
            ndim_virt = size(det_map%det(i)%virtuals_energies,1)
            do j = 1, ndim_virt
               if ( (det_map%det(i)%virtuals_energies(j).ge.relcc_projectors_rew_virt_min_energy).and. &
                  & (det_map%det(i)%virtuals_energies(j).le.relcc_projectors_rew_virt_max_energy) ) then
                     spinors_in_virtuals_window(j) = .true.
               else
                     spinors_in_virtuals_window(j) = .false.
               end if
            end do
            if (ndim_virt.eq.1) spinors_in_virtuals_window(2) = spinors_in_virtuals_window(1)
         end if 
         if (relcc_projectors_rew_strict) then 
! strictly enforces that only configurations with ionization/excitations/attachments within the window are kept
            create_projector(i) = spinors_in_occupied_window(1).and.spinors_in_virtuals_window(1)
            create_projector(i) = create_projector(i).and.spinors_in_occupied_window(2).and.spinors_in_virtuals_window(2)
         else
! more loosely enforces ionization/excitations/attachments within the window; requires only that at least one of the idexes in virtual/occupied space belongs to the window.
! for ionization/excitations, if there's no restriction on the virtual space, this is equivalent to cvs 
            create_projector(i) = spinors_in_virtuals_window(1).or.spinors_in_virtuals_window(2)
            create_projector(i) = create_projector(i).and.(spinors_in_occupied_window(1).or.spinors_in_occupied_window(2))
         end if 

         if (iprnt.ge.2) then
            if (create_projector(i)) then
               write (iw,*) 'determinant ',i,' is within  REW'
            else
               write (iw,*) 'determinant ',i,' is outside REW'
            end if
            write (iw,*) '    occupied:',det_map%det(i)%occupied_energies
            write (iw,*) '    virtuals:',det_map%det(i)%virtuals_energies
         end if
      end do
!   end if
end function 

! function to create the operator to project out core contribution in the t amplitudes
! it returs an array of logical variables, where the entries are set to be true if a given line is
! to be suppressed. 
function create_t_projector(ndet, eps)
#include "files.inc"
#include "inpt.inc"
   logical, pointer :: create_t_projector(:) 

   integer, intent(in) :: ndet
   real(kind=8), intent(in) :: eps(*)

   type(mapping) :: det_map
   integer :: i, j, ndim_occ, ndim_virt
   logical :: spinors_in_occupied_window(2)

      call mapping_spinors_configurations_for_t (det_map, ndet, eps)

      allocate(create_t_projector(ndet))
      create_t_projector = .false.
      spinors_in_occupied_window = .false.

! we loop over the composition of each determinant in the "excitation" manifold 
      do i = 1, ndet
         if (allocated(det_map%det(i)%occupied_energies)) then 
            ndim_occ = size(det_map%det(i)%occupied_energies,1)
            do j = 1, ndim_occ
               if ( (det_map%det(i)%occupied_energies(j).le.relcc_projectors_frozen_core_max_energy) ) then
                  spinors_in_occupied_window(j) = .true.
               else
                  spinors_in_occupied_window(j) = .false.
               end if
            end do
            if (ndim_occ.eq.1) spinors_in_occupied_window(2) = spinors_in_occupied_window(1)
         end if
        create_t_projector(i) = (spinors_in_occupied_window(1).or.spinors_in_occupied_window(2))

         if (iprnt.ge.2) then
            if (create_t_projector(i)) then
               write (iw,*) 'determinant ',i,' is core'
            else
               write (iw,*) 'determinant ',i,' is valence'
            end if
            write (iw,*) '    occupied orbital(s) energy:',det_map%det(i)%occupied_energies
         end if
       end do
end function 



   subroutine mapping_spinors_configurations (map, select_eom, ndt, state_sym, eps)

!---------------Description--------------------------------------------
!
!     create a mapping between the index of a determinant (i.e. its position on R/L (eigen)vectors Analyze the excited states and calculate transition moments
!     
!---------------Arguments----------------------------------------------
!     mapping    : derived type storing information on the determinants
!     select_eom : selection of eomtype
!     eps        : orbital energies
!     print_level: print level
!
!---------------Last modified------------------------------------------
!
!     Author : Andre Gomes, adapting from Avijit Shee's evector_analysis 
!
!---------------Calling variables--------------------------------------

#include "complex.inc"
#include "symm.inc"

   real*8,  intent(in) :: eps(*)
   integer, intent(in) :: select_eom 
   integer, intent(in) :: ndt
   integer, intent(in) :: state_sym
   type(mapping), intent(inout) :: map 

   integer            :: nstate_sym(nrep,4)

!---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "ccpar.inc"

!---------------Local variables--------------------------------------
   character(8) :: eomtype
   integer       :: i,ii,irep,j,aoff,a,aa,ai,arep,ssym
   integer       :: ijrp, abrp, jrp, jj, kk, irp, ioff, imin, bb, amin, arp, brp, b

!---------------Executable code--------------------------------------
 
    ii = 1
    allocate(map%det(ndt))
    if (select_eom == 2) eomtype = 'EA'
    if (select_eom == 3) eomtype = 'IP'
    if (select_eom == 4) eomtype = 'EE'
    select case(eomtype)
      case('EA')

         do irep=1,nrep
            arep = multb (irep,1+nrep,1)
            do i = 1, ncont(irep)
               aoff = io(nrep+1) + iv(arep)
               do a = 1, nv(arep)
                  aa = aoff + a

                  allocate(map%det(ii)%virtuals_energies(1))
                  allocate(map%det(ii)%virtuals_irreps(1))
                  allocate(map%det(ii)%virtuals_indexes_in_irreps(1))

                  map%det(ii)%virtuals_energies(1) = eps(aa)
                  map%det(ii)%virtuals_irreps(1) = arep
                  map%det(ii)%virtuals_indexes_in_irreps(1) = a+no(arep)

                  ii = ii + 1
               enddo
            enddo
         enddo

         do irep=1,nrep
            ijrp = irep
            abrp = multb (ijrp+nrep,1+nrep,1)
            do jrp = 1, nrep
               jj = io(jrp)
               irp = multb(jrp,ijrp+nrep,2)
               do j = 1, no(jrp)
                  jj = jj + 1
                  do i = 1, ncont(irp)
                     do brp = 1, nrep
                        bb = iv(brp) + io(nrep+1)
                        arp = multb(brp,abrp+nrep,2)
                        if (arp.lt.brp) cycle
                           aoff = iv(arp) + io(nrep+1)
                           do b = 1, nv(brp)
                              bb = bb + 1
                              amin = 1
                              if (arp.eq.brp) amin = b + 1
                              do a = amin, nv(arp)
                                 aa = aoff + a

                                 allocate(map%det(ii)%occupied_energies(1))
                                 allocate(map%det(ii)%occupied_irreps(1))
                                 allocate(map%det(ii)%occupied_indexes_in_irreps(1))


                                 map%det(ii)%occupied_energies(1) = eps(jj)
                                 map%det(ii)%occupied_irreps(1) = jrp 
                                 map%det(ii)%occupied_indexes_in_irreps(1) = j 

                                 allocate(map%det(ii)%virtuals_energies(2))
                                 allocate(map%det(ii)%virtuals_irreps(2))
                                 allocate(map%det(ii)%virtuals_indexes_in_irreps(2))

                                 map%det(ii)%virtuals_energies(1) = eps(aa)
                                 map%det(ii)%virtuals_energies(2) = eps(bb)

                                 map%det(ii)%virtuals_irreps(1) = arp
                                 map%det(ii)%virtuals_irreps(2) = brp

                                 map%det(ii)%virtuals_indexes_in_irreps(1) = a+no(arp)
                                 map%det(ii)%virtuals_indexes_in_irreps(2) = b+no(brp)

                                 ii = ii + 1
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo

      case('IP')

         jj = 0
         do irep = 1, nrep
            arep = multb (irep,1+nrep,1)
            do i = 1, no(irep)
               jj = jj + 1
               do a = 1, ncont(arep)

                  allocate(map%det(ii)%occupied_energies(1))
                  allocate(map%det(ii)%occupied_irreps(1))
                  allocate(map%det(ii)%occupied_indexes_in_irreps(1))

                  map%det(ii)%occupied_energies(1) = eps(jj)
                  map%det(ii)%occupied_irreps(1) = irep
                  map%det(ii)%occupied_indexes_in_irreps(1) = i

                  ii = ii + 1
               enddo
            enddo
         enddo

         do ijrp = 1, nrep
!           abrp = multb (ijrp+nrep,1+nrep,1)
            abrp = multb (1+nrep,ijrp+nrep,2)
            do jrp = 1, nrep
               jj = io(jrp)
               irp = multb(jrp,ijrp+nrep,2)
               if (irp.lt.jrp) cycle
               ioff = io(irp)
               do j = 1, no(jrp)
                  jj = jj + 1
                  imin = 1
                  if (irp.eq.jrp) imin = j + 1
                  do i = imin, no(irp)
                     kk = ioff + i
                     do brp = 1, nrep
                        bb = iv(brp) + io(nrep+1)
                        arp = multb(brp,abrp+nrep,2)
                        aoff = iv(arp) + io(nrep+1)
                        do b = 1, nv(brp)
                           bb = bb + 1
                           do a = 1, ncont(arp)
                              allocate(map%det(ii)%occupied_energies(2))
                              allocate(map%det(ii)%occupied_irreps(2))
                              allocate(map%det(ii)%occupied_indexes_in_irreps(2))

                              map%det(ii)%occupied_energies(1) = eps(kk)
                              map%det(ii)%occupied_energies(2) = eps(jj)

                              map%det(ii)%occupied_irreps(1) = irp
                              map%det(ii)%occupied_irreps(2) = jrp

                              map%det(ii)%occupied_indexes_in_irreps(1) = i
                              map%det(ii)%occupied_indexes_in_irreps(2) = j

                              allocate(map%det(ii)%virtuals_energies(1))
                              allocate(map%det(ii)%virtuals_irreps(1))
                              allocate(map%det(ii)%virtuals_indexes_in_irreps(1))

                              map%det(ii)%virtuals_energies(1) = eps(bb)
                              map%det(ii)%virtuals_irreps(1) = brp
                              map%det(ii)%virtuals_indexes_in_irreps(1) = b+no(brp)

                              ii = ii + 1
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo


      case('EE')

         jj = 0
         do irep = 1, nrep
            arep = multb (state_sym+nrep,irep,2)
            do i = 1, no(irep)
               jj = jj + 1
               aoff = io(nrep+1) + iv(arep)
               do a = 1, nv(arep)
                  aa = aoff + a

                  allocate(map%det(ii)%occupied_energies(1))
                  allocate(map%det(ii)%occupied_irreps(1))
                  allocate(map%det(ii)%occupied_indexes_in_irreps(1))
                  
                  map%det(ii)%occupied_energies(1) = eps(jj)
                  map%det(ii)%occupied_irreps(1) = irep
                  map%det(ii)%occupied_indexes_in_irreps(1) = i
               
                  allocate(map%det(ii)%virtuals_energies(1))
                  allocate(map%det(ii)%virtuals_irreps(1))
                  allocate(map%det(ii)%virtuals_indexes_in_irreps(1))
      
                  map%det(ii)%virtuals_energies(1) = eps(aa)
                  map%det(ii)%virtuals_irreps(1) = arep
                  map%det(ii)%virtuals_indexes_in_irreps(1) = a+no(arep)

                  ii = ii + 1
               enddo
            enddo
         enddo

         do ijrp = 1, nrep
            abrp = multb (state_sym+nrep,ijrp+nrep,2)
            do jrp = 1, nrep
               jj = io(jrp)
               irp = multb(jrp,ijrp+nrep,2)
               if (irp.lt.jrp) cycle
               ioff = io(irp)
               do j = 1, no(jrp)
                  jj = jj + 1
                  imin = 1
                  if (irp.eq.jrp) imin = j + 1
                  do i = imin, no(irp)
                     kk = ioff + i
                     do brp = 1, nrep
                        bb = iv(brp) + io(nrep+1)
                        arp = multb(brp,abrp+nrep,2)
                        if (arp.lt.brp) cycle
                        aoff = iv(arp) + io(nrep+1)
                        do b = 1, nv(brp)
                           bb = bb + 1
                           amin = 1
                           if (arp.eq.brp) amin = b + 1
                           do a = amin, nv(arp)
                              aa = aoff + a

                              allocate(map%det(ii)%occupied_energies(2))
                              allocate(map%det(ii)%occupied_irreps(2))
                              allocate(map%det(ii)%occupied_indexes_in_irreps(2))

                              map%det(ii)%occupied_energies(1) = eps(kk)
                              map%det(ii)%occupied_energies(2) = eps(jj)

                              map%det(ii)%occupied_irreps(1) = irp
                              map%det(ii)%occupied_irreps(2) = jrp

                              map%det(ii)%occupied_indexes_in_irreps(1) = i
                              map%det(ii)%occupied_indexes_in_irreps(2) = j

                              allocate(map%det(ii)%virtuals_energies(2))
                              allocate(map%det(ii)%virtuals_irreps(2))
                              allocate(map%det(ii)%virtuals_indexes_in_irreps(2))

                              map%det(ii)%virtuals_energies(1) = eps(aa)
                              map%det(ii)%virtuals_energies(2) = eps(bb)

                              map%det(ii)%virtuals_irreps(1) = arp
                              map%det(ii)%virtuals_irreps(2) = brp

                              map%det(ii)%virtuals_indexes_in_irreps(1) = a+no(arp)
                              map%det(ii)%virtuals_indexes_in_irreps(2) = b+no(brp)

                              ii = ii + 1
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

      end select
 
  end subroutine


   subroutine mapping_spinors_configurations_for_t (map, ndt, eps)

!---------------Description--------------------------------------------
!
!     create a mapping between the index of a determinant 
!     
!---------------Arguments----------------------------------------------
!     mapping    : derived type storing information on the determinants
!     eps        : orbital energies
!
!---------------Last modified------------------------------------------
!
!     Author : Andre Gomes, adapting from Avijit Shee's evector_analysis 
!
!---------------Calling variables--------------------------------------

#include "complex.inc"
#include "symm.inc"

   real*8,  intent(in) :: eps(*)
   integer, intent(in) :: ndt
   type(mapping), intent(inout) :: map 

   integer            :: nstate_sym(nrep,4)

!---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "ccpar.inc"

!---------------Local variables--------------------------------------
   integer       :: i,ii,irep,j,aoff,a,aa,ai,arep,ssym
   integer       :: ijrp, abrp, jrp, jj, kk, irp, ioff, imin, bb, amin, arp, brp, b

!---------------Executable code--------------------------------------
 
    ii = 1
    allocate(map%det(ndt))

         jj = 0
         do irep = 1, nrep
            arep = irep !multb (state_sym+nrep,irep,2)
            do i = 1, no(irep)
               jj = jj + 1
               aoff = io(nrep+1) + iv(arep)
               do a = 1, nv(arep)
                  aa = aoff + a

                  allocate(map%det(ii)%occupied_energies(1))
                  allocate(map%det(ii)%occupied_irreps(1))
                  allocate(map%det(ii)%occupied_indexes_in_irreps(1))
                  
                  map%det(ii)%occupied_energies(1) = eps(jj)
                  map%det(ii)%occupied_irreps(1) = irep
                  map%det(ii)%occupied_indexes_in_irreps(1) = i
               
                  allocate(map%det(ii)%virtuals_energies(1))
                  allocate(map%det(ii)%virtuals_irreps(1))
                  allocate(map%det(ii)%virtuals_indexes_in_irreps(1))
      
                  map%det(ii)%virtuals_energies(1) = eps(aa)
                  map%det(ii)%virtuals_irreps(1) = arep
                  map%det(ii)%virtuals_indexes_in_irreps(1) = a+no(arep)

                  ii = ii + 1
               enddo
            enddo
         enddo

         do ijrp = 1, nrep
            abrp = ijrp !multb (state_sym+nrep,ijrp+nrep,2)
            do jrp = 1, nrep
               jj = io(jrp)
               irp = multb(jrp,ijrp+nrep,2)
               if (irp.lt.jrp) cycle
               ioff = io(irp)
               do j = 1, no(jrp)
                  jj = jj + 1
                  imin = 1
                  if (irp.eq.jrp) imin = j + 1
                  do i = imin, no(irp)
                     kk = ioff + i
                     do brp = 1, nrep
                        bb = iv(brp) + io(nrep+1)
                        arp = multb(brp,abrp+nrep,2)
                        if (arp.lt.brp) cycle
                        aoff = iv(arp) + io(nrep+1)
                        do b = 1, nv(brp)
                           bb = bb + 1
                           amin = 1
                           if (arp.eq.brp) amin = b + 1
                           do a = amin, nv(arp)
                              aa = aoff + a

                              allocate(map%det(ii)%occupied_energies(2))
                              allocate(map%det(ii)%occupied_irreps(2))
                              allocate(map%det(ii)%occupied_indexes_in_irreps(2))

                              map%det(ii)%occupied_energies(1) = eps(kk)
                              map%det(ii)%occupied_energies(2) = eps(jj)

                              map%det(ii)%occupied_irreps(1) = irp
                              map%det(ii)%occupied_irreps(2) = jrp

                              map%det(ii)%occupied_indexes_in_irreps(1) = i
                              map%det(ii)%occupied_indexes_in_irreps(2) = j

                              allocate(map%det(ii)%virtuals_energies(2))
                              allocate(map%det(ii)%virtuals_irreps(2))
                              allocate(map%det(ii)%virtuals_indexes_in_irreps(2))

                              map%det(ii)%virtuals_energies(1) = eps(aa)
                              map%det(ii)%virtuals_energies(2) = eps(bb)

                              map%det(ii)%virtuals_irreps(1) = arp
                              map%det(ii)%virtuals_irreps(2) = brp

                              map%det(ii)%virtuals_indexes_in_irreps(1) = a+no(arp)
                              map%det(ii)%virtuals_indexes_in_irreps(2) = b+no(brp)

                              ii = ii + 1
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
 
  end subroutine
   end module 

