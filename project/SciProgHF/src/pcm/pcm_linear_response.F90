module pcm_linear_response

#ifdef HAS_PCMSOLVER
    use, intrinsic :: iso_c_binding

    implicit none

    public pcm_linear_response_initialize
    public pcm_linear_response_finalize

    public pcm_oper_ao_driver
    public compute_mep_asc
    public pcm_oit_asc

    private

    ! If false the interface will refuse to be accessed
    logical                     :: is_initialized = .false.
    type(c_ptr)                 :: context_
    real(c_double), allocatable :: tess_cent(:)
    integer(c_size_t)           :: nr_points
    integer(c_size_t)           :: nr_points_irr
    ! A (maybe clumsy) way of passing LUPRI without using common blocks
    integer                     :: global_print_unit
    ! A counter for the number of LR iterations
    integer, save               :: LR_iteration_counter = 1

contains

    subroutine pcm_linear_response_initialize(print_unit)

        use pcmsolver, only: PCMSOLVER_READER_HOST, PCMSOLVER_READER_OWN, &
                             PCMInput, &
                             pcmsolver_new, &
                             pcmsolver_is_compatible_library, &
                             pcmsolver_get_cavity_size,      &
                             pcmsolver_get_irreducible_cavity_size, &
                             pcmsolver_get_centers
        use pcm_gp, only: collect_nctot, collect_atoms, collect_symmetry_info
        use pcmmod_cfg

        integer, intent(in) :: print_unit

        integer(c_int) :: nr_nuclei
        real(c_double), allocatable :: charges(:)
        real(c_double), allocatable :: centers(:)
        integer(c_int) :: symmetry_info(4)
        type(PCMInput) :: host_input

        global_print_unit = print_unit

        if (.not. pcmsolver_is_compatible_library()) then
            write(global_print_unit, *) 'Error: incompatible version of PCMSolver library'
            stop 1
        end if

        nr_nuclei = collect_nctot()
        allocate(charges(nr_nuclei))
        charges = 0.0_c_double
        allocate(centers(3*nr_nuclei))
        centers = 0.0_c_double
        call collect_atoms(charges, centers)
        call collect_symmetry_info(symmetry_info)

        if (pcmmod_host_provides_input) then
            host_input = pcmsolver_input()
            context_ = pcmsolver_new(PCMSOLVER_READER_HOST, nr_nuclei, charges, centers, symmetry_info, host_input)
        else
            context_ = pcmsolver_new(PCMSOLVER_READER_OWN, nr_nuclei, charges, centers, symmetry_info, host_input)
        end if

        deallocate(charges)
        deallocate(centers)

        nr_points = pcmsolver_get_cavity_size(context_)
        nr_points_irr = pcmsolver_get_irreducible_cavity_size(context_)
        allocate(tess_cent(3*nr_points))
        tess_cent = 0.0d0
        call pcmsolver_get_centers(context_, tess_cent)

        is_initialized = .true.

    end subroutine pcm_linear_response_initialize

    subroutine pcm_linear_response_finalize()

        use pcmsolver, only: pcmsolver_delete

        if (.not. is_initialized) then
            print *, 'Error: pcm_linear_response was never initialized.'
            stop 1
        end if
        ! Free the memory taken from the free store both in Fortran and in C++
        deallocate(tess_cent)
        call pcmsolver_delete(context_)


        is_initialized = .false.

    end subroutine pcm_linear_response_finalize

    subroutine check_if_interface_is_initialized()

        if (.not. is_initialized) then
            print *, 'Error: pcm_linear_response is not initialized.'
            stop 1
        end if

    end subroutine

    subroutine pcm_oper_ao_driver(oper, charge_name, work, lwork)
        !
        ! Polarizable Continuum Model.
        ! Driver routine for Fock matrix contributions, J and X(D).
        !
        ! Roberto Di Remigio, 2012
        !
        use pcm_integrals, only: get_electronic_mep
        use pcmmod_cfg
        use pcmsolver, only: pcmsolver_get_surface_function

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

        ! Passed variables
        ! oper(*) is just the PCM contribution in AO basis
        real(8), intent(out)     :: oper(ntbas(0), ntbas(0), nz)
        character(*), intent(in) :: charge_name
        real(8)                  :: work(*)
        integer                  :: lwork
        ! Local variables
        real(c_double), allocatable :: asc(:)
        real(8)                     :: average
        real(8), allocatable        :: oper_mo(:, :, :)
        integer                     :: i, j

        call check_if_interface_is_initialized

        allocate(asc(nr_points))
        asc = 0.0d0
        call pcmsolver_get_surface_function(context_, nr_points, asc, charge_name)
        ! Contraction of electrostatic potential integrals with polarization charges
        call get_electronic_mep(nr_points, nr_points_irr, oper, tess_cent,   &
            asc, work, lwork, 0, .true., pcmmod_skipss)
        ! Symmetrize Fock matrix contribution
        do i = 1, ntbas(0)
            do j = 1, i
                average = 0.5d0 * (oper(i, j, 1) + oper(j, i, 1))
                oper(i, j, 1) = average
                oper(j, i, 1) = average
            end do
        end do

        if (pcmmod_dospf) then
            allocate(oper_MO(norb(1), norb(1), nz))
            oper_MO = 0.0d0
            print *, "FOCK_PCM before elimination of spin-dependent terms."
            call prqmat(oper, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                nz, ipqtoq(1,0), 6)
            call spfao(oper_MO, oper, work, lwork)
            print *, "FOCK_PCM after elimination of spin-dependent terms."
            call prqmat(oper, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
                nz, ipqtoq(1,0), 6)
            print *, "FOCK_PCM_MO after elimination of spin-dependent terms."
            call prqmat(oper_mo, norb(1), norb(1), norb(1), norb(1), &
                nz, ipqtoq(1,0), 6)
            deallocate(oper_MO)
        end if

        deallocate(asc)

    end subroutine pcm_oper_ao_driver

    subroutine compute_mep_asc(density_matrix, work, lfree)
        !
        ! Calculate the molecular electrostatic potential and
        ! the apparent surface charge at the cavity points.
        !
        ! The user can control via the DIRAC input the following:
        !    * switch between separate and total evaluation of the
        !      nuclear and electronic parts;
        !
        use pcm_integrals, only: get_nuclear_mep, get_electronic_mep, get_mep
        use pcm_write, only: pcm_write_file, pcm_write_file_separate
        use pcmmod_cfg
        use pcmsolver, only: pcmsolver_set_surface_function, &
                             pcmsolver_compute_asc, &
                             pcmsolver_get_surface_function

        real(8), intent(in)    :: density_matrix(*)
        real(8), intent(inout) :: work(*)
        integer                :: lfree
        ! Local variables
        real(c_double), allocatable :: mep(:)
        real(c_double), allocatable :: asc(:)
        real(c_double), allocatable :: nuc_pot(:), nuc_pol_chg(:)
        real(c_double), allocatable :: ele_pot(:), ele_pol_chg(:)
        character(7)                :: potName, chgName
        character(7)                :: potName1, chgName1, potName2, chgName2
        integer                     :: i
        integer(c_int)              :: irrep

        allocate(mep(nr_points))
        mep = 0.0d0
        allocate(asc(nr_points))
        asc = 0.0d0
        ! The irreducible representation, for the moment the totally
        ! symmetric one
        irrep = 0_c_int

        if (.not.(pcmmod_separate)) then
            potName = 'TotMEP'//c_null_char
            chgName = 'TotASC'//c_null_char
            ! Calculate the (total) Molecular Electrostatic Potential
            call get_mep(nr_points, nr_points_irr, tess_cent, mep,        &
                density_matrix, work, lfree, 0, pcmmod_skipss)
            ! Set a cavity surface function with the MEP
            call pcmsolver_set_surface_function(context_, nr_points, mep, potName)
            ! Compute polarization charges and set the proper surface function
            call pcmsolver_compute_asc(context_, potName, chgName, irrep)
            ! Get polarization charges @tesserae centers
            call pcmsolver_get_surface_function(context_, nr_points, asc, chgName)

            ! Print some information
            if (pcmmod_print > 5) then
                write(global_print_unit, '(A, T27, A, T62, A)') "Finite element #", "Total MEP", "Total ASC"
                do i = 1, nr_points
                    write(global_print_unit, '(I6, 2(20X, F15.12))') i, mep(i), asc(i)
                end do
            end if

            ! Write to file MEP and ASC
            call pcm_write_file(nr_points, mep, asc)
        else
            ! Allocation
            allocate(nuc_pot(nr_points))
            nuc_pot = 0.0d0
            allocate(nuc_pol_chg(nr_points))
            nuc_pol_chg = 0.0d0
            allocate(ele_pot(nr_points))
            ele_pot = 0.0d0
            allocate(ele_pol_chg(nr_points))
            ele_pol_chg = 0.0d0

            potName1 = 'NucMEP'//c_null_char
            chgName1 = 'NucASC'//c_null_char
            call get_nuclear_mep(nr_points, nr_points_irr, tess_cent, nuc_pot, 0)
            call pcmsolver_set_surface_function(context_, nr_points, nuc_pot, potName1)
            call pcmsolver_compute_asc(context_, potName1, chgName1, irrep)
            call pcmsolver_get_surface_function(context_, nr_points, nuc_pol_chg, chgName1)

            potName2 = 'EleMEP'//c_null_char
            chgName2 = 'EleASC'//c_null_char
            call get_electronic_mep(nr_points, nr_points_irr, density_matrix, tess_cent,  &
                ele_pot, work, lfree, 0, .false., pcmmod_skipss)
            call pcmsolver_set_surface_function(context_, nr_points, ele_pot, potName2)
            call pcmsolver_compute_asc(context_, potName2, chgName2, irrep)
            call pcmsolver_get_surface_function(context_, nr_points, ele_pol_chg, chgName2)

            ! Print some information
            if (pcmmod_print > 5) then
                write(global_print_unit, '(A, T27, A, T62, A, T97, A, T132, A)') "Finite element #", &
                    "Nuclear MEP", "Nuclear ASC", "Electronic MEP", "Electronic ASC"
                do i = 1, nr_points
                    write(global_print_unit, '(I6, 4(20X, F15.12))') i, nuc_pot(i), nuc_pol_chg(i), ele_pot(i), ele_pol_chg(i)
                end do
            end if

            ! Obtain vector of total MEP
            potName  = 'TotMEP'//c_null_char
            mep(:) = nuc_pot(:) + ele_pot(:)
            call pcmsolver_set_surface_function(context_, nr_points, mep, potName)

            ! Obtain vector of total polarization charges
            chgName  = 'TotASC'//c_null_char
            asc(:) = nuc_pol_chg(:) + ele_pol_chg(:)
            call pcmsolver_set_surface_function(context_, nr_points, asc, chgName)

            ! Write to file MEP and ASC
            call pcm_write_file_separate(nr_points, nuc_pot, nuc_pol_chg, ele_pot, ele_pol_chg)

            deallocate(nuc_pot)
            deallocate(nuc_pol_chg)
            deallocate(ele_pot)
            deallocate(ele_pol_chg)
        end if
        deallocate(mep)
        deallocate(asc)

    end subroutine compute_mep_asc

    subroutine pcm_oit_asc(mep_oit_asc, perturbed_density_matrix, nr_density_matrices, work, lwork)
        !
        ! Polarizable Continuum Model.
        ! Driver routine for "explicit" contributions to linear response.
        ! perturbed_density_matrix is the b-transformed AO density matrix.
        ! oit stands for One-Index Transformed
        !
        ! Roberto Di Remigio, 2012-2013
        !

        use pcm_integrals, only: get_electronic_mep
        use pcmmod_cfg
        use pcmsolver, only: pcmsolver_set_surface_function, &
                             pcmsolver_compute_asc, &
                             pcmsolver_get_surface_function

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

        ! Passed variables
        real(c_double) :: mep_oit_asc(ntbas(0), ntbas(0), nz)
        real(8)        :: perturbed_density_matrix(*)
        integer        :: nr_density_matrices
        real(8)        :: work(*)
        integer        :: lwork
        ! Local variables
        real(c_double), allocatable :: mep_oit(:), asc_oit(:)
        character(8)                :: potName, chgName
        real(8)                     :: average
        integer                     :: i, j
        integer(c_int)              :: irrep

        call check_if_interface_is_initialized

        allocate(mep_oit(nr_points))
        mep_oit = 0.0d0
        potName = 'MEP_OIT'//c_null_char
        allocate(asc_oit(nr_points))
        asc_oit = 0.0d0
        chgName = 'ASC_OIT'//c_null_char
        ! The irreducible representation, for the moment the totally
        ! symmetric one
        irrep = 0_c_int
        ! Compute electronic MEP @tesserae centers
        do i = 1, nr_density_matrices
            call get_electronic_mep(nr_points, nr_points_irr, perturbed_density_matrix(i), tess_cent, &
                mep_oit, work, lwork, 0, .false., pcmmod_skipss)
        enddo
        call pcmsolver_set_surface_function(context_, nr_points, mep_oit, potName)
        ! Compute polarization charges and set the proper surface function
        call pcmsolver_compute_asc(context_, potName, chgName, irrep)
        call pcmsolver_get_surface_function(context_, nr_points, asc_oit, chgName)
        ! Now contract unpertubed potential matrix elements with perturbed ASC
        call get_electronic_mep(nr_points, nr_points_irr, mep_oit_asc, tess_cent,  &
            asc_oit, work, lwork, 0, .true., pcmmod_skipss)
        do i = 1, ntbas(0)
            do j = 1, i
                average = 0.5d0 * (mep_oit_asc(i, j, 1) + mep_oit_asc(j, i, 1))
                mep_oit_asc(i, j, 1) = average
                mep_oit_asc(j, i, 1) = average
            end do
        end do

        deallocate(mep_oit)
        deallocate(asc_oit)

        ! Update the LR iteration counter
        LR_iteration_counter = LR_iteration_counter + 1

    end subroutine
#endif /* HAS_PCMSOLVER */

end module pcm_linear_response
