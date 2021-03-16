
!    radovan bast <radovan.bast@uit.no> 2011-04-20
!    module for the export/import of basis set information

!    general design
!    --------------
!  - file is formatted
!  - file is read in free format
!  - sections have labels
!  - order of sections does not matter
!  - section symmetry can be omitted if irrelevant (e.g. C1)
!  - normalization is absorbed in contraction coefficients
!    (we can discuss this point)
!    AOs are normalized based on angular momentum and exponent
!    such that each shell has a common normalization
!    with this < AO | AO > =   1 for s, px, py, pz, dxy, dxz, dyz, fxyz
!                              3 for dxx, dyy, dzz, fxxy, ...
!                             15     fxxx, fyyy, fzzz, gxxxy, ...
!                            105     gxxxx, ...
!                            ...     ...
!                              9 for gxxyy, ...
!                            etc     ...
!  - ireps go from 0 to nr_boson_ireps-1

module interface_ao

   use file_units

   implicit none

   public interface_ao_read
   public interface_ao_clear

   integer, parameter, public :: l_max = 20

   logical, public :: interface_ao_set

   integer, public :: ao_fermion_corep(0:7, 2)
   integer, public :: ao_nr(0:7, 2)
   integer, public :: ao_off(0:7, 2)
   integer, public :: ao_off_hermit(8, 2)
   integer, public :: ao_q(0:7, 2)
   integer, public :: iq(0:3)
   integer, public :: irep_of_axes(3)
   integer, public :: irep_of_rotations(3)
   integer, public :: ao_block_start(0:16) !cheat: 0 start is to catch cases where one of AO blocks is empty
   integer, public :: ao_block_nr(0:16)    !       in dens eval these blocks are skipped
   integer, public :: nr_ao_blocks
   integer, public :: llss_block_partner(16, 0:3, 0:7)
   integer, public :: lssl_block_partner(16, 0:3, 0:7)
   real(8), public :: llss_prefactor(16)
   integer, public :: nr_ao
   integer, public :: nr_ao_cartesian
   integer, public :: nr_ao_gerade
   integer, public :: nr_ao_large
   integer, public :: nr_ao_small
   integer, public :: nr_ao_spherical
   integer, public :: nr_ao_ungerade
   integer, public :: nr_boson_ireps
   integer, public :: nr_fermion_coreps
   integer, public :: nr_quaternion_blocks
   integer, public :: nr_shells
   integer, public :: nr_shells_large
   integer, public :: nr_shells_small
   integer, public :: nr_so_ao
   integer, public :: nr_centers
   integer, public :: pq_to_uq(4, 0:7)
   integer, public :: uq_to_pq(4, 0:7)
   logical, public :: basis_is_spherical
   logical, public :: use_only_large
   !real(8), public :: global_gauge_origin(3)

   integer, allocatable, public :: cartesian_deg(:)
   integer, allocatable, public :: cartesian_off(:)
   integer, allocatable, public :: i_so_ao(:, :)
   integer, allocatable, public :: l_quantum_nr(:)
   integer, allocatable, public :: nr_nonzero_contraction_coefs(:)
   integer, allocatable, public :: shell_center(:)
   integer, allocatable, public :: shell_stabilizing_irep(:)
   integer, allocatable, public :: spherical_deg(:)
   integer, allocatable, public :: spherical_off(:)
   integer, allocatable, public :: center_deg(:)
   integer, allocatable, public :: center_stabilizing_irep(:)
   integer, allocatable, public :: ao_center(:)
   real(8), allocatable, public :: ao_center_xyz(:, :)
   real(8), allocatable, public :: ao_distance(:, :, :)
   real(8), allocatable, public :: contraction_coef(:, :)
   real(8), allocatable, public :: f_so_ao(:)
   real(8), allocatable, public :: primitive_exp(:, :)
   real(8), allocatable, public :: center_charge(:)
   real(8), allocatable, public :: center_xyz(:, :)

   private

   integer              :: c2s_index(l_max)
   real(8), allocatable :: c2s_mat(:)

contains

   subroutine interface_ao_clear()

      interface_ao_set           = .false.

      ao_fermion_corep           = 0
      ao_nr                      = 0
      ao_off                     = 0
      ao_off_hermit              = 0
      ao_q                       = 0
      iq                         = 0
      irep_of_axes               = 0
      irep_of_rotations          = 0
      ao_block_start             = 0
      ao_block_nr                = 0
      nr_ao_blocks               = 0
      llss_block_partner         = 0
      lssl_block_partner         = 0
      llss_prefactor             = 0.0d0
      nr_ao                      = 0
      nr_ao_cartesian            = 0
      nr_ao_gerade               = 0
      nr_ao_large                = 0
      nr_ao_small                = 0
      nr_ao_spherical            = 0
      nr_ao_ungerade             = 0
      nr_boson_ireps             = 0
      nr_fermion_coreps          = 0
      nr_quaternion_blocks       = 0
      nr_shells                  = 0
      nr_shells_large            = 0
      nr_shells_small            = 0
      nr_so_ao                   = 0
      nr_centers = 0
      pq_to_uq                   = 0
      uq_to_pq                   = 0
      basis_is_spherical         = .false.
      use_only_large             = .false.
      !global_gauge_origin        = 0.0d0

      if (allocated(cartesian_deg                )) deallocate(cartesian_deg                )
      if (allocated(cartesian_off                )) deallocate(cartesian_off                )
      if (allocated(i_so_ao                      )) deallocate(i_so_ao                      )
      if (allocated(l_quantum_nr                 )) deallocate(l_quantum_nr                 )
      if (allocated(nr_nonzero_contraction_coefs )) deallocate(nr_nonzero_contraction_coefs )
      if (allocated(shell_center                 )) deallocate(shell_center                 )
      if (allocated(shell_stabilizing_irep       )) deallocate(shell_stabilizing_irep       )
      if (allocated(spherical_deg                )) deallocate(spherical_deg                )
      if (allocated(spherical_off                )) deallocate(spherical_off                )
      if (allocated(center_deg                   )) deallocate(center_deg                   )
      if (allocated(center_stabilizing_irep      )) deallocate(center_stabilizing_irep      )
      if (allocated(ao_center                    )) deallocate(ao_center                    )
      if (allocated(ao_center_xyz                )) deallocate(ao_center_xyz                )
      if (allocated(ao_distance                  )) deallocate(ao_distance                  )
      if (allocated(contraction_coef             )) deallocate(contraction_coef             )
      if (allocated(f_so_ao                      )) deallocate(f_so_ao                      )
      if (allocated(primitive_exp                )) deallocate(primitive_exp                )
      if (allocated(center_charge                )) deallocate(center_charge                )
      if (allocated(center_xyz                   )) deallocate(center_xyz                   )

   end subroutine

   subroutine interface_ao_read(use_mpi_if_possible)

      use sigma_prefactor_setting
      use xc_mpi

!     --------------------------------------------------------------------------
      logical, intent(in)  :: use_mpi_if_possible
!     --------------------------------------------------------------------------
      character(120)       :: line
      integer              :: i, j, k, l, m, n
      integer              :: sum_nr_primitives
      integer              :: ioff_c, ioff_s
      integer              :: idummy, irep, cent
      logical              :: symmetry_info_found
      real(8)              :: prim, coef
      character(1)         :: ch
      integer              :: inter_irep(0:3), jxyz, irep1, irep2, ic, nr1, nr2, mc, icount
      integer              :: nr_ao_irep, isym, deg
      integer              :: lvar, mvar, nvar
      integer              :: iao, iso, iso_dirac_sorting
      integer              :: ishell, nshell
      integer              :: nround
      integer              :: ixyz
      integer              :: iblock, jblock
      integer              :: nl, nt, nc, ns, kc, ks, oc, os
      integer              :: lssl_ic(2) = (/2, 1/)
      real(8)              :: f, px, py, pz
      integer, allocatable :: isymao(:, :)
      integer, allocatable :: irepao(:)
      integer, allocatable :: ipind(:)
      integer, allocatable :: kstart(:)
      real(8), allocatable :: cs_mat(:, :)
      real(8), parameter   :: bitstring_parity(0:7) = (/1.0d0, -1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0/)
      integer :: size1_ao_center
      integer :: size1_ao_center_xyz
      integer :: size2_ao_center_xyz
      integer :: size1_ao_distance
      integer :: size2_ao_distance
      integer :: size3_ao_distance
      integer :: size1_cartesian_deg
      integer :: size1_cartesian_off
      integer :: size1_contraction_coef
      integer :: size2_contraction_coef
      integer :: size1_f_so_ao
      integer :: size1_i_so_ao
      integer :: size2_i_so_ao
      integer :: size1_l_quantum_nr
      integer :: size1_nr_nonzero_contraction_coefs
      integer :: size1_primitive_exp
      integer :: size2_primitive_exp
      integer :: size1_shell_center
      integer :: size1_shell_stabilizing_irep
      integer :: size1_spherical_deg
      integer :: size1_spherical_off
      integer :: size1_center_charge
      integer :: size1_center_deg
      integer :: size1_center_stabilizing_irep
      integer :: size1_center_xyz
      integer :: size2_center_xyz
!     --------------------------------------------------------------------------

      symmetry_info_found = .false.

      call interface_ao_clear()

      if (xc_mpi_is_master()) then
         inquire(file = interface_file_name_ao, exist = interface_file_exists)
         if (.not. interface_file_exists) then
            print *, 'interface_ao file not found'
            stop
         end if
         open(interface_file_unit,             &
              file   = interface_file_name_ao, &
              status = 'unknown',              &
              form   = 'formatted',            &
              access = 'sequential')
         rewind(interface_file_unit)

         do while (.true.)
            read(interface_file_unit, '(a120)', end=1) line

            if (word_contains(line, '*** geometry')) then
               read(interface_file_unit, *)
               read(interface_file_unit, *) nr_centers
               read(interface_file_unit, *)
               allocate(center_charge(nr_centers))
               allocate(center_xyz(nr_centers, 3))
               do i = 1, nr_centers
                  read(interface_file_unit, *) center_charge(i), &
                                     center_xyz(i, 1:3)
               end do
            end if

            if (word_contains(line, '*** symmetry')) then
               symmetry_info_found = .true.
               read(interface_file_unit, *)
               allocate(center_deg(nr_centers))
               allocate(center_stabilizing_irep(nr_centers))
               do i = 1, nr_centers
                  read(interface_file_unit, *) idummy, center_deg(i), center_stabilizing_irep(i)
               end do
               read(interface_file_unit, *)
               read(interface_file_unit, *) nr_boson_ireps
               read(interface_file_unit, *)
               read(interface_file_unit, *) nr_fermion_coreps
               read(interface_file_unit, *)
               read(interface_file_unit, *) irep_of_axes(1:3)
               read(interface_file_unit, *)
               read(interface_file_unit, *) irep_of_rotations(1:3)
               read(interface_file_unit, *)
               do irep = 0, nr_boson_ireps - 1
                  read(interface_file_unit, *) idummy,                    &
                                     ao_nr(irep, 1),            &
                                     ao_nr(irep, 2),            &
                                     ao_off(irep, 1),           &
                                     ao_off(irep, 2),           &
                                     ao_q(irep, 1),             &
                                     ao_q(irep, 2),             &
                                     ao_fermion_corep(irep, 1), &
                                     ao_fermion_corep(irep, 2), &
                                     ao_off_hermit(irep+1, 1),  &
                                     ao_off_hermit(irep+1, 2)
               end do
               read(interface_file_unit, *)
               do i = 1, 4
                  read(interface_file_unit, *) pq_to_uq(i, 0:7)
               end do
               do i = 1, 4
                  read(interface_file_unit, *) uq_to_pq(i, 0:7)
               end do
               read(interface_file_unit, *)
               read(interface_file_unit, *) iq(0:3)
            end if

            if (word_contains(line, '*** basis')) then
               read(interface_file_unit, *)
               read(interface_file_unit, *) basis_is_spherical
               read(interface_file_unit, *)
               read(interface_file_unit, *) nr_quaternion_blocks
               read(interface_file_unit, *)
               read(interface_file_unit, *) use_only_large
               read(interface_file_unit, *)
               read(interface_file_unit, *) sum_nr_primitives
               read(interface_file_unit, *)

!              overkill allocation for convenience
!              it will "never become a bottleneck" (2011-04-29)
               allocate(shell_center(sum_nr_primitives))
               allocate(l_quantum_nr(sum_nr_primitives))
               allocate(nr_nonzero_contraction_coefs(sum_nr_primitives))
               allocate(cartesian_deg(sum_nr_primitives))
               allocate(cartesian_off(sum_nr_primitives))
               allocate(spherical_deg(sum_nr_primitives))
               allocate(spherical_off(sum_nr_primitives))
               allocate(primitive_exp(sum_nr_primitives, sum_nr_primitives))
               allocate(contraction_coef(sum_nr_primitives, sum_nr_primitives))

               j = 0
               k = 0
               ioff_c = 0
               ioff_s = 0
               nr_shells_large = 0
               nr_shells_small = 0
               do i = 1, sum_nr_primitives
                  read(interface_file_unit, *) cent,   &
                                     ch,     &
                                     ishell, &
                                     l,      &
                                     prim,   &
                                     coef
                  if (ishell > k) then
                     j = 1
                     shell_center(ishell)           = cent
                     l_quantum_nr(ishell)           = l
                     cartesian_deg(ishell) = (l+1)*(l+2)/2
                     cartesian_off(ishell) = ioff_c
                     ioff_c = ioff_c + cartesian_deg(ishell)
                     spherical_deg(ishell) = 2*l + 1
                     spherical_off(ishell) = ioff_s
                     ioff_s = ioff_s + spherical_deg(ishell)
                     if (ch == 'L') then
                        nr_shells_large = nr_shells_large + 1
                     end if
                     if (ch == 'S') then
                        nr_shells_small = nr_shells_small + 1
                     end if
                  else
                     j = j + 1
                  end if
                  primitive_exp(ishell, j)    = prim
                  contraction_coef(ishell, j) = coef
                  nr_nonzero_contraction_coefs(ishell) = j
                  k = ishell
               end do
            end if

         end do
1        continue

         close(interface_file_unit, status = 'keep')

         nr_shells = nr_shells_large + nr_shells_small

         if (.not. symmetry_info_found) then
            allocate(center_deg(nr_centers))
            allocate(center_stabilizing_irep(nr_centers))
            center_deg              = 1
            center_stabilizing_irep = 0

            nr_boson_ireps    = 1
            nr_fermion_coreps = 1
            irep_of_axes      = 0
            irep_of_rotations = 0

            ao_nr            = 0
            ao_off           = 0
            ao_q             = 0
            ao_fermion_corep = 0
            ao_off_hermit    = 0

            nl = 0
            do ishell = 1, nr_shells_large
               if (basis_is_spherical) then
                  deg = spherical_deg(ishell)
               else
                  deg = cartesian_deg(ishell)
               end if
               nl = nl + deg
            end do

            nt = 0
            do ishell = 1, nr_shells
               if (basis_is_spherical) then
                  deg = spherical_deg(ishell)
               else
                  deg = cartesian_deg(ishell)
               end if
               nt = nt + deg
            end do

            ao_nr(0, 1) = nl
            ao_nr(0, 2) = nt - nl
            ao_off(0, 1) = 0
            ao_off(0, 2) = nl
            ao_q(0, 1) = 1
            ao_q(0, 2) = 1
            ao_fermion_corep(0, 1) = 1
            ao_fermion_corep(0, 2) = 1
            ao_off_hermit(1, 1) = 0
            ao_off_hermit(1, 2) = nl

            pq_to_uq = 0
            do i = 1, 4
               pq_to_uq(i, 0) = i
            end do
            uq_to_pq = 0
            do i = 1, 4
               uq_to_pq(i, 0) = i
            end do

            iq(0) = 1
            iq(1) = 4
            iq(2) = 3
            iq(3) = 2
         end if

         if (use_only_large) then
            mc = 1
         else
            mc = 2
         end if

         i = 0
         llss_prefactor = 1.0d0
         do irep = 0, nr_boson_ireps - 1
            do ic = 1, mc
               if (ao_nr(irep, ic) > 0) then
                  i = i + 1
                  ao_block_nr(i)    = ao_nr(irep, ic)
                  ao_block_start(i) = ao_off(irep, ic) + 1
                  if (ic == 2) then
                     llss_prefactor(i) = sigma_prefactor
                  end if
               end if
            end do
         end do
         nr_ao_blocks = i

!        transform as rotation
         inter_irep(0) = 0
         do jxyz = 1, 3
            inter_irep(jxyz) = irep_of_rotations(jxyz)
         end do
         do irep1 = 0, nr_boson_ireps - 1
            do ic = 1, mc
               do irep = 0, nr_boson_ireps - 1
                  do jxyz = 0, 3
                     irep2 = ieor(irep1, ieor(inter_irep(jxyz), irep))
                     if ((ao_nr(irep1, ic) > 0) .and. (ao_nr(irep2, ic) > 0)) then
                        do i = 1, nr_ao_blocks
                           if (ao_block_start(i) == (ao_off(irep1, ic) + 1)) then
                              iblock = i
                           end if
                        end do
                        do i = 1, nr_ao_blocks
                           if (ao_block_start(i) == (ao_off(irep2, ic) + 1)) then
                              jblock = i
                           end if
                        end do
                        llss_block_partner(iblock, jxyz, irep) = jblock
                     end if
                  end do
               end do
            end do
         end do

!        transform as translation
         inter_irep(0) = 0
         do jxyz = 1, 3
            inter_irep(jxyz) = irep_of_axes(jxyz)
         end do
         do irep1 = 0, nr_boson_ireps - 1
            do ic = 1, mc
               do irep = 0, nr_boson_ireps - 1
                  do jxyz = 0, 3
                     irep2 = ieor(irep1, ieor(inter_irep(jxyz), irep))
                     if ((ao_nr(irep1, ic) > 0) .and. (ao_nr(irep2, lssl_ic(ic)) > 0)) then
                        do i = 1, nr_ao_blocks
                           if (ao_block_start(i) == (ao_off(irep1, ic) + 1)) then
                              iblock = i
                           end if
                        end do
                        do i = 1, nr_ao_blocks
                           if (ao_block_start(i) == (ao_off(irep2, lssl_ic(ic)) + 1)) then
                              jblock = i
                           end if
                        end do
                        lssl_block_partner(iblock, jxyz, irep) = jblock
                     end if
                  end do
               end do
            end do
         end do

         if (basis_is_spherical) then
            allocate(isymao(l_max, 2*l_max + 1))
            isymao = 0
            do l = 0, l_max-1
               do i = 1, 2*l + 1
                  isymao(l + 1, i) = ireplm(l, mdef(l, i))
               end do
            end do
         else
            allocate(isymao(l_max, (l_max+1)*(l_max+2)/2))
            isymao = 0
            do l = 1, l_max
               k = 0
               do m = 1, l
                  do n = 1, m
                     k = k + 1
                     lvar = mod(l-m, 2)*irep_of_axes(1)
                     mvar = mod(m-n, 2)*irep_of_axes(2)
                     nvar = mod(n-1, 2)*irep_of_axes(3)
                     isymao(l, k) = ieor(lvar, ieor(mvar, nvar))
                  end do
               end do
            end do
         end if

         allocate(shell_stabilizing_irep(nr_shells))
         do ishell = 1, nr_shells
            shell_stabilizing_irep(ishell) = center_stabilizing_irep(shell_center(ishell))
         end do

         nr_ao_cartesian = 0
         nr_ao_spherical = 0
         do ishell = 1, nr_shells
            do irep = 0, nr_boson_ireps - 1
               if (iand(irep, shell_stabilizing_irep(ishell)) == 0) then
                  nr_ao_cartesian = nr_ao_cartesian + cartesian_deg(ishell)
                  nr_ao_spherical = nr_ao_spherical + spherical_deg(ishell)
               end if
            end do
         end do

         if (basis_is_spherical) then
            nr_ao = nr_ao_spherical
         else
            nr_ao = nr_ao_cartesian
         end if

         allocate(irepao(nr_ao))
         allocate(ipind(nr_ao))
         irepao = 0
         ipind = 0
         iso = 0
         do irep = 0, nr_boson_ireps - 1
            iao = 0
            do ishell = 1, nr_shells
               if (basis_is_spherical) then
                  deg = spherical_deg(ishell)
               else
                  deg = cartesian_deg(ishell)
               end if
               do k = 1, deg
                  iao = iao + 1
                  if (iand(shell_stabilizing_irep(ishell), ieor(irep, isymao(l_quantum_nr(ishell) + 1, k))) == 0) then
                     iso = iso + 1
                     ipind(iso) = ishell*2**16 + k*2**8 + irepao(iao)
                     irepao(iao) = irepao(iao) + 1
                  end if
               end do
            end do
         end do
         deallocate(irepao)

         allocate(kstart(nr_shells))
         kstart = 0
         do ishell = 2, nr_shells
            if (basis_is_spherical) then
               deg = spherical_deg(ishell-1)
            else
               deg = cartesian_deg(ishell-1)
            end if
            kstart(ishell) = kstart(ishell - 1) + center_deg(shell_center(ishell-1))*deg
         end do

         call setup_c2s_mat()

         allocate(cs_mat(nr_ao_cartesian, nr_ao))
         cs_mat = 0.0d0
         oc = 0
         os = 0
         do ishell = 1, nr_shells
            do irep = 0, nr_boson_ireps - 1
               if (iand(irep, shell_stabilizing_irep(ishell)) == 0) then
                  if (basis_is_spherical) then
                     if (l_quantum_nr(ishell) > 1) then
                        k = c2s_index(l_quantum_nr(ishell))
                        do kc = 1, cartesian_deg(ishell)
                           do ks = 1, spherical_deg(ishell)
                              cs_mat(oc + kc, os + ks) = c2s_mat(k)
                              k = k + 1
                           end do
                        end do
                     else
                        do kc = 1, cartesian_deg(ishell)
                           cs_mat(oc + kc, os + kc) = 1.0d0
                        end do
                     end if
                  else
                     do kc = 1, cartesian_deg(ishell)
                        cs_mat(oc + kc, oc + kc) = 1.0d0
                     end do
                  end if
                  oc = oc + cartesian_deg(ishell)
                  os = os + spherical_deg(ishell)
               end if
            end do
         end do
         deallocate(c2s_mat)

         nr_so_ao = 0
         do nround = 1, 2
!           two rounds: first round gets dimensions
!                       second round allocates and sets arrays

            if (nround == 2) then
               allocate(i_so_ao(2, nr_so_ao))
               allocate(f_so_ao(nr_so_ao))
               nr_so_ao = 0
            end if

            i = 1
            do irep = 0, nr_boson_ireps - 1
               icount = 0
               nr_ao_irep = ao_nr(irep, 1) + ao_nr(irep, 2)
               do iso = i, i + nr_ao_irep - 1

                  ishell = iand(ishft(ipind(iso), -16), 65535)
                  nshell = iand(ishft(ipind(iso),  -8),   255)

                  iao = kstart(ishell) + nshell

                  icount = icount + 1
                  if (icount > ao_nr(irep, 1)) then
!                    is small
                     ic = 2
                  else
!                    is large
                     ic = 1
                  end if

!                 do not plan to AO-SO transform small component functions
!                 in case of zero small-small metric (.LEVY-LEBLOND and .X2C4)
                  if (.not. (use_only_large .and. (ic == 2))) then

!                    shift by difference dirac sorting minus hermit sorting
                     iso_dirac_sorting = iso + ao_off(irep, ic) - ao_off_hermit(irep + 1, ic)

                     do isym = 0, nr_boson_ireps - 1
                        if (iand(isym, shell_stabilizing_irep(ishell)) == 0) then

                           f = bitstring_parity(iand(isym, ieor(irep, isymao(l_quantum_nr(ishell) + 1, nshell))))

                           do kc = 1, nr_ao_cartesian
                              if (dabs(cs_mat(kc, iao)) > 0.0d0) then
                                 nr_so_ao = nr_so_ao + 1
                                 if (nround == 2) then
                                    i_so_ao(1, nr_so_ao) = iso_dirac_sorting
                                    i_so_ao(2, nr_so_ao) = kc
                                    f_so_ao(nr_so_ao)    = f*cs_mat(kc, iao)
                                 end if
                              end if
                           end do

                           if (basis_is_spherical) then
                              iao = iao + spherical_deg(ishell)
                           else
                              iao = iao + cartesian_deg(ishell)
                           end if

                        endif
                     end do
                  end if
               end do
               i = i + nr_ao_irep
            end do
         end do

         deallocate(cs_mat)
         deallocate(isymao)
         deallocate(ipind)
         deallocate(kstart)

         nr_ao_gerade   = 0
         nr_ao_ungerade = 0
         if (nr_fermion_coreps == 1) then
            nr_ao_gerade = nr_ao
         else
            do irep = 0, nr_boson_ireps - 1
               select case (ao_fermion_corep(irep, 1))
                  case (1)
                     nr_ao_gerade = nr_ao_gerade + ao_nr(irep, 1)
                     if (nr_shells_small > 0) then
                        nr_ao_ungerade = nr_ao_ungerade + ao_nr(irep, 2)
                     end if
                  case (2)
                     nr_ao_ungerade = nr_ao_ungerade + ao_nr(irep, 1)
                     if (nr_shells_small > 0) then
                        nr_ao_gerade = nr_ao_gerade + ao_nr(irep, 2)
                     end if
               end select
            end do
         end if

         allocate(ao_center(nr_ao))
         allocate(ao_center_xyz(nr_ao, 3))
         allocate(ao_distance(nr_ao, nr_ao, 3))
         i = 1
         do ishell = 1, nr_shells
            if (basis_is_spherical) then
               deg = spherical_deg(ishell)
            else
               deg = cartesian_deg(ishell)
            end if
            do irep = 0, nr_boson_ireps - 1
               if (iand(irep, shell_stabilizing_irep(ishell)) == 0) then

                  px = bitstring_parity(iand(irep_of_axes(1), irep))*center_xyz(shell_center(ishell), 1)
                  py = bitstring_parity(iand(irep_of_axes(2), irep))*center_xyz(shell_center(ishell), 2)
                  pz = bitstring_parity(iand(irep_of_axes(3), irep))*center_xyz(shell_center(ishell), 3)

                  do j = i, i + deg - 1
                     ao_center_xyz(j, 1) = px
                     ao_center_xyz(j, 2) = py
                     ao_center_xyz(j, 3) = pz
                     ao_center(j) = shell_center(ishell)
                  end do

                  i = i + deg
               end if
            end do
         end do
         do ixyz = 1, 3
            do l = 1, nr_ao
               do k = 1, nr_ao
                  ao_distance(k, l, ixyz) = ao_center_xyz(k, ixyz) - ao_center_xyz(l, ixyz)
               end do
            end do
         end do

         nr_ao_large = 0
         nr_ao_small = 0
         do irep = 0, nr_boson_ireps - 1
            nr_ao_large = nr_ao_large + ao_nr(irep, 1)
            nr_ao_small = nr_ao_small + ao_nr(irep, 2)
         end do

!        limit active nr of shells
         if (use_only_large) then
            nr_shells = nr_shells_large
         end if

         interface_ao_set = .true.
      end if

#ifdef VAR_MPI
      if (use_mpi_if_possible .and. (xc_mpi_get_nr_proc() > 1)) then

         call xc_mpi_bcast(nr_ao             )
         call xc_mpi_bcast(nr_so_ao          )
         call xc_mpi_bcast(nr_centers        )
         call xc_mpi_bcast(sum_nr_primitives )

         if (xc_mpi_is_master()) then
            size1_ao_center                    = size(ao_center)
            size1_ao_center_xyz                = size(ao_center_xyz, 1)
            size2_ao_center_xyz                = size(ao_center_xyz, 2)
            size1_ao_distance                  = size(ao_distance, 1)
            size2_ao_distance                  = size(ao_distance, 2)
            size3_ao_distance                  = size(ao_distance, 3)
            size1_cartesian_deg                = size(cartesian_deg)
            size1_cartesian_off                = size(cartesian_off)
            size1_contraction_coef             = size(contraction_coef, 1)
            size2_contraction_coef             = size(contraction_coef, 2)
            size1_f_so_ao                      = size(f_so_ao)
            size1_i_so_ao                      = size(i_so_ao, 1)
            size2_i_so_ao                      = size(i_so_ao, 2)
            size1_l_quantum_nr                 = size(l_quantum_nr)
            size1_nr_nonzero_contraction_coefs = size(nr_nonzero_contraction_coefs)
            size1_primitive_exp                = size(primitive_exp, 1)
            size2_primitive_exp                = size(primitive_exp, 2)
            size1_shell_center                 = size(shell_center)
            size1_shell_stabilizing_irep       = size(shell_stabilizing_irep)
            size1_spherical_deg                = size(spherical_deg)
            size1_spherical_off                = size(spherical_off)
            size1_center_charge                = size(center_charge)
            size1_center_deg                   = size(center_deg)
            size1_center_stabilizing_irep      = size(center_stabilizing_irep)
            size1_center_xyz                   = size(center_xyz, 1)
            size2_center_xyz                   = size(center_xyz, 2)
         else
            size1_ao_center                    = 0
            size1_ao_center_xyz                = 0
            size2_ao_center_xyz                = 0
            size1_ao_distance                  = 0
            size2_ao_distance                  = 0
            size3_ao_distance                  = 0
            size1_cartesian_deg                = 0
            size1_cartesian_off                = 0
            size1_contraction_coef             = 0
            size2_contraction_coef             = 0
            size1_f_so_ao                      = 0
            size1_i_so_ao                      = 0
            size2_i_so_ao                      = 0
            size1_l_quantum_nr                 = 0
            size1_nr_nonzero_contraction_coefs = 0
            size1_primitive_exp                = 0
            size2_primitive_exp                = 0
            size1_shell_center                 = 0
            size1_shell_stabilizing_irep       = 0
            size1_spherical_deg                = 0
            size1_spherical_off                = 0
            size1_center_charge                = 0
            size1_center_deg                   = 0
            size1_center_stabilizing_irep      = 0
            size1_center_xyz                   = 0
            size2_center_xyz                   = 0
         end if

         call xc_mpi_bcast(size1_ao_center                    )
         call xc_mpi_bcast(size1_ao_center_xyz                )
         call xc_mpi_bcast(size2_ao_center_xyz                )
         call xc_mpi_bcast(size1_ao_distance                  )
         call xc_mpi_bcast(size2_ao_distance                  )
         call xc_mpi_bcast(size3_ao_distance                  )
         call xc_mpi_bcast(size1_cartesian_deg                )
         call xc_mpi_bcast(size1_cartesian_off                )
         call xc_mpi_bcast(size1_contraction_coef             )
         call xc_mpi_bcast(size2_contraction_coef             )
         call xc_mpi_bcast(size1_f_so_ao                      )
         call xc_mpi_bcast(size1_i_so_ao                      )
         call xc_mpi_bcast(size2_i_so_ao                      )
         call xc_mpi_bcast(size1_l_quantum_nr                 )
         call xc_mpi_bcast(size1_nr_nonzero_contraction_coefs )
         call xc_mpi_bcast(size1_primitive_exp                )
         call xc_mpi_bcast(size2_primitive_exp                )
         call xc_mpi_bcast(size1_shell_center                 )
         call xc_mpi_bcast(size1_shell_stabilizing_irep       )
         call xc_mpi_bcast(size1_spherical_deg                )
         call xc_mpi_bcast(size1_spherical_off                )
         call xc_mpi_bcast(size1_center_charge                )
         call xc_mpi_bcast(size1_center_deg                   )
         call xc_mpi_bcast(size1_center_stabilizing_irep      )
         call xc_mpi_bcast(size1_center_xyz                   )
         call xc_mpi_bcast(size2_center_xyz                   )

         if (.not. xc_mpi_is_master()) then
            allocate(ao_center(size1_ao_center))
            allocate(ao_center_xyz(size1_ao_center_xyz, size2_ao_center_xyz))
            allocate(ao_distance(size1_ao_distance, size2_ao_distance, size3_ao_distance))
            allocate(cartesian_deg(size1_cartesian_deg))
            allocate(cartesian_off(size1_cartesian_off))
            allocate(contraction_coef(size1_contraction_coef, size2_contraction_coef))
            allocate(f_so_ao(size1_f_so_ao))
            allocate(i_so_ao(size1_i_so_ao, size2_i_so_ao))
            allocate(l_quantum_nr(size1_l_quantum_nr))
            allocate(nr_nonzero_contraction_coefs(size1_nr_nonzero_contraction_coefs))
            allocate(primitive_exp(size1_primitive_exp, size2_primitive_exp))
            allocate(shell_center(size1_shell_center))
            allocate(shell_stabilizing_irep(size1_shell_stabilizing_irep))
            allocate(spherical_deg(size1_spherical_deg))
            allocate(spherical_off(size1_spherical_off))
            allocate(center_charge(size1_center_charge))
            allocate(center_deg(size1_center_deg))
            allocate(center_stabilizing_irep(size1_center_stabilizing_irep))
            allocate(center_xyz(size1_center_xyz, size2_center_xyz))
         end if

         call xc_mpi_bcast(interface_ao_set              )

         call xc_mpi_bcast(ao_fermion_corep              )
         call xc_mpi_bcast(ao_nr                         )
         call xc_mpi_bcast(ao_off                        )
         call xc_mpi_bcast(ao_off_hermit                 )
         call xc_mpi_bcast(ao_q                          )
         call xc_mpi_bcast(iq                            )
         call xc_mpi_bcast(irep_of_axes                  )
         call xc_mpi_bcast(irep_of_rotations             )
         call xc_mpi_bcast(ao_block_start                )
         call xc_mpi_bcast(ao_block_nr                   )
         call xc_mpi_bcast(nr_ao_blocks                  )
         call xc_mpi_bcast(llss_block_partner            )
         call xc_mpi_bcast(lssl_block_partner            )
         call xc_mpi_bcast(llss_prefactor                )
         call xc_mpi_bcast(nr_ao_cartesian               )
         call xc_mpi_bcast(nr_ao_gerade                  )
         call xc_mpi_bcast(nr_ao_large                   )
         call xc_mpi_bcast(nr_ao_small                   )
         call xc_mpi_bcast(nr_ao_spherical               )
         call xc_mpi_bcast(nr_ao_ungerade                )
         call xc_mpi_bcast(nr_boson_ireps                )
         call xc_mpi_bcast(nr_fermion_coreps             )
         call xc_mpi_bcast(nr_quaternion_blocks          )
         call xc_mpi_bcast(nr_shells                     )
         call xc_mpi_bcast(nr_shells_large               )
         call xc_mpi_bcast(nr_shells_small               )
         call xc_mpi_bcast(pq_to_uq                      )
         call xc_mpi_bcast(uq_to_pq                      )
         call xc_mpi_bcast(basis_is_spherical            )
         call xc_mpi_bcast(use_only_large                )
         !call xc_mpi_bcast(global_gauge_origin           )

         call xc_mpi_bcast(cartesian_deg                 )
         call xc_mpi_bcast(cartesian_off                 )
         call xc_mpi_bcast(i_so_ao                       )
         call xc_mpi_bcast(l_quantum_nr                  )
         call xc_mpi_bcast(nr_nonzero_contraction_coefs  )
         call xc_mpi_bcast(shell_center                  )
         call xc_mpi_bcast(shell_stabilizing_irep        )
         call xc_mpi_bcast(spherical_deg                 )
         call xc_mpi_bcast(spherical_off                 )
         call xc_mpi_bcast(center_deg    )
         call xc_mpi_bcast(center_stabilizing_irep    )
         call xc_mpi_bcast(ao_center                     )
         call xc_mpi_bcast(ao_center_xyz                 )
         call xc_mpi_bcast(ao_distance                   )
         call xc_mpi_bcast(contraction_coef              )
         call xc_mpi_bcast(f_so_ao                       )
         call xc_mpi_bcast(primitive_exp                 )
         call xc_mpi_bcast(center_charge )
         call xc_mpi_bcast(center_xyz    )
      end if
#endif

   end subroutine

   function word_contains(word, substring)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: word
      character(*), intent(in) :: substring
!     --------------------------------------------------------------------------
      logical                  :: word_contains
!     --------------------------------------------------------------------------

      if (index(word, substring) > 0) then
         word_contains = .true.
      else
         word_contains = .false.
      end if

   end function

   function fac(n)

!     --------------------------------------------------------------------------
      integer, intent(in) :: n
      real(8)             :: fac
!     it is real(8) so that expressions don't overflow!
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      fac = 1.0d0
      do i = 1, n
         fac = fac*dble(i)
      end do

   end function

   function fac2(n)

!     --------------------------------------------------------------------------
      integer, intent(in) :: n
      real(8)             :: fac2
!     it is real(8) so that expressions don't overflow!
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      if (n < 0) then
         fac2 = dble(n + 2)
         do i = n + 4, 1, 2
            fac2 = fac2*dble(i)
         end do
         fac2 = 1.0d0/fac2
      else if (n == 0) then
         fac2 = 1.0d0
      else
         fac2 = dble(n)
         do i = n - 2, 1, -2
            fac2 = fac2*dble(i)
         end do
      end if

   end function

   function binom(n, k)

!     taken from: http://blog.plover.com/math/choose.html
!     thanks to Ulf for sending me this link
!     --------------------------------------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: k
      integer             :: binom
!     --------------------------------------------------------------------------
      integer             :: j, m
!     --------------------------------------------------------------------------

      if (k == 0) then
         binom = 1
         return
      end if
      if (n == 0) then
         binom = 0
         return
      end if

      m = n
      binom = 1
      do j = 1, k
         binom = binom*m
         binom = binom/j
         m = m - 1
      end do

   end function

   subroutine setup_c2s_mat()

!     --------------------------------------------------------------------------
      integer              :: nc, ns
      integer              :: i, j, k, l, m, n
      integer              :: ix, ilm
      real(8)              :: cm, cmk, cmki, cmkij
      real(8), allocatable :: tmat(:, :)
      real(8), allocatable :: cossin_coef(:, :)
      real(8), allocatable :: legendre_coef(:)
!     --------------------------------------------------------------------------


      c2s_index = 0
      n = 0
      do l = 2, l_max
         c2s_index(l) = n + 1
         n = n + ((l+1)*(l+2)/2)**2
      end do
      allocate(c2s_mat(n))
      c2s_mat = 0.0d0

      do l = 2, l_max

         nc = (l+1)*(l+2)/2
         ns = 2*l + 1

         allocate(legendre_coef(0:l))
         legendre_coef = 0.0d0
         do k = 0, l/2
            legendre_coef(l - 2*k) = (dble((-1)**k)/dble(2**l))*binom(l, k)*binom(2*(l - k), l)
         end do

         allocate(cossin_coef(0:l, 0:l))
         cossin_coef = 0.0d0
         do m = 0, l
            cossin_coef(m, 0) = 1.0d0
            do k = 1, m
               cossin_coef(m, k) = cossin_coef(m - 1, k - 1)*dble((-1)**(k - 1))
               if (m > k) cossin_coef(m, k) = cossin_coef(m, k) + cossin_coef(m - 1, k)
            end do
         end do

         allocate(tmat(nc, ns))
         tmat = 0.0d0
         do m = 0, l
            cm = sqrt(2.0d0*fac(l - m)/fac(l + m))
            if (m  ==  0) cm = 1.0d0
            cm = cm/dsqrt(dble(fac2(2*l - 1)))
            do k = mod(l - m, 2), l - m, 2
               if (m > 0) legendre_coef(k) = dble((k + 1))*legendre_coef(k + 1)
               cmk = cm*legendre_coef(k)
               do i = 0, (l - k - m)/2
                  cmki = cmk*binom((l - k - m)/2,i)
                  do j = 0, i
                     cmkij = cmki*binom(i, j)
                     do n = 0, m
                        ix = l - 2*j - m + n
                        ix = ix*(ix + 1)/2 + l + 1 - m - 2*i
                        if (mod(n, 2) == 1) then
                           ilm = 1 + l - m
                        else
                           ilm = 1 + l + m
                        end if
                        tmat(ix, ilm) = tmat(ix, ilm) + cmkij*cossin_coef(m, n)
                     end do
                  end do
               end do
            end do
         end do

         k = c2s_index(l)
         do i = 1, nc
            do j = 1, ns
               c2s_mat(k) = tmat(i, j)
               k = k + 1
            end do
         end do

         deallocate(tmat)
         deallocate(cossin_coef)
         deallocate(legendre_coef)

      end do

   end subroutine

   function mdef(l, i)

!     --------------------------------------------------------------------------
      integer, intent(in) :: l
      integer, intent(in) :: i
      integer             :: mdef
!     --------------------------------------------------------------------------

      if (l == 1) then
          select case (i)
             case (1)
                mdef = 1
             case (2)
                mdef = -1
             case (3)
                mdef = 0
          end select
      else
         mdef = i - l - 1
      endif

   end function

   function ireplm(l, m)

!     --------------------------------------------------------------------------
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer             :: ireplm
!     --------------------------------------------------------------------------

      ireplm = 0
      if (mod(l + m,  2) == 1) then
         ireplm = irep_of_axes(3)
      end if
      if (mod(abs(m), 2) == 1) then
         ireplm = ieor(ireplm, irep_of_axes(1))
      end if
      if (m < 0) then
         ireplm = ieor(ireplm, ieor(irep_of_axes(1), irep_of_axes(2)))
      end if

   end function

end module
