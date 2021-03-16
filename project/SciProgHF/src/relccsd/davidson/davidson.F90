module davidson


implicit none

public davidson_get_eigenvalue_shift
public davidson_setup
public davidson_cleanup
public davidson_driver
public davidson_print_vectors
public solve_values_vectors
public sort_values_vectors

!public sort_values_vectors
!public set_safe_minimum
!public davidson_create_results
!public free_davidson_results

interface davidson_driver
    module procedure get_values_vectors_davidson
end interface

public davidson_results

private

! derived types

type davidson_results
    real(kind=8), pointer :: eValues(:,:) 
    real(kind=8), pointer :: eVectorsR(:,:)
    real(kind=8), pointer :: eVectorsL(:,:) 
end type


! module parameters

integer, parameter :: max_label_length = 30
integer, parameter :: nr_max_irreps = 32
real(kind=8), parameter :: davidson_eigenvalues_shift_threshold = 1.0e-3

! module variables, should be set with a call to davidson_setup

real(kind=8), save :: safe_minimum = 1.0d-36

real(kind=8), save :: davidson_convergence_threshold 
real(kind=8), save :: davidson_eigenvalues_shift(2)

integer, save      :: davidson_max_subspace_size 
integer, save      :: davidson_max_iterations
integer, save      :: davidson_output_file_unit 
integer, save      :: davidson_subspaces_file_unit 
integer, save      :: davidson_sigmas_file_unit 
integer, save      :: davidson_results_file_unit 
integer, save      :: davidson_trial_vector_refresh_rate
integer, save      :: davidson_maximum_trial_vector_refresh_rate = 10

logical, save      :: davidson_apply_eigenvalues_shift
logical, save      :: davidson_solve_rhs 
logical, save      :: davidson_solve_lhs 
logical, save      :: davidson_solve_values
logical, save      :: davidson_verbose 
logical, save      :: davidson_symmetric 
logical, save      :: davidson_complex_mode_with_reals 
logical, save      :: davidson_save_subspaces 
logical, save      :: davidson_save_sigmas 
logical, save      :: davidson_save_results
logical, save      :: davidson_overlap_sorting
logical, save      :: davidson_debug = .false.
logical, save      :: davidson_use_preconditioner = .true.
logical, save      :: davidson_refresh_trial_vectors = .false.
logical, save      :: davidson_print_config = .true.

! 
! variables for controlling i/o
!
character(len=16), save :: left_sigma_file_name     = 'LSIGMAV.sym'
character(len=16), save :: left_vectors_file_name   = 'LEVECTO.sym'
character(len=16), save :: left_subspace_file_name  = 'LSUBSPA.sym'

character(len=16), save :: right_sigma_file_name    = 'RSIGMAV.sym'
character(len=16), save :: right_vectors_file_name  = 'REVECTO.sym'
character(len=16), save :: right_subspace_file_name = 'RSUBSPA.sym'

integer, save :: right_sigma_file_unit    = 140
integer, save :: right_vectors_file_unit  = 141
integer, save :: right_subspace_file_unit = 142

integer, save :: left_sigma_file_unit     = 240
integer, save :: left_vectors_file_unit   = 241
integer, save :: left_subspace_file_unit  = 242


contains 

!
! create_sigma_vector is the interface used with the diagonalizer. 
! inside it we will handle left and right sigma vector cases, and the sigma_eom module 
! will have to be initialized so that it knows whether to do left or right vectors
!  
   function create_sigma_vector(A, r, old_sigma_vector, state_sym)
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal
      use sigma_eom, only : intermediates, create_sigma_vector_left, create_sigma_vector_right

      real(kind=8), pointer :: create_sigma_vector(:,:)
      real(kind=8), pointer :: old_sigma_vector(:,:)
      type(intermediates) :: A
      real(kind=8), intent(in) :: r(:,:)
      integer, intent(in) :: state_sym

      if (davidson_solve_rhs.and..not.davidson_solve_lhs) then 
!     do right vectors, this is the default
          create_sigma_vector => create_sigma_vector_right(A, r, old_sigma_vector, state_sym)

      else if (davidson_solve_lhs.and..not.davidson_solve_rhs) then
!     do left vectors vectors
          create_sigma_vector => create_sigma_vector_left(A, r, old_sigma_vector, state_sym)

      else
! do both at the same time; this will require some thinking, one option is to return
! a matrix of (N,(Lr+Lf)), L(r/l) the number of rhs/lhs vecturs, containing Lr rhs vectors 
! in the first Lr columns and the transpose of sigma_L in the remaining columns. it is easy
! if there's always the same number of columns flr left and right, otherwise there's no 
! simple way to keep track of the dimensions.
          create_sigma_vector => NULL()
      endif 
   end function


! this should never be used to update an instance, only to create new ones. it's better to have results
! for left and right separately if we determine these in two steps. but we must check whether an 
! element of the derived type has been allocated already... perhaps this should then be changed to a
! function returning a results type, and have the davidson driver return that to its caller
    subroutine davidson_create_results(results, N, nroots, values, left, right)
        type(davidson_results), intent(inout) :: results
        integer, intent(in)    :: N, nroots 
        integer                :: N_local
        logical, intent(in)    :: values, left, right

! we don't use complex variables but rather use twice the amount of real variables
! for symmetric matrices we only need one of the vectors, and we choose it to be the right one
! left vectors, if asked, are ignored
        N_local = N

        results%eValues => NULL()
        results%eVectorsR => NULL()
        results%eVectorsL => NULL()


        if (davidson_symmetric) then
            if (values) then 
                allocate(results%eValues(1,nroots))
                results%eValues   = 0.0d0
            end if
            if (right) then
                allocate(results%eVectorsR(N_local,nroots))
                results%eVectorsR = 0.0d0
            end if
        else
            if (values) then
                allocate(results%eValues(2,nroots))
                results%eValues   = 0.0d0
            end if
            if (right) then
                allocate(results%eVectorsR(N_local,nroots))
                results%eVectorsR = 0.0d0
            else
                allocate(results%eVectorsR(1,nroots))
                results%eVectorsR = 0.0d0
            end if
            if (left) then
                allocate(results%eVectorsL(N_local,nroots))
                results%eVectorsL = 0.0d0
            else
                allocate(results%eVectorsL(1,nroots))
                results%eVectorsL = 0.0d0
            end if
        end if
    end subroutine


    subroutine free_davidson_results(results)
        type(davidson_results), intent(inout) :: results

        if (associated(results%eValues)) deallocate(results%eValues)
        if (associated(results%eVectorsR)) deallocate(results%eVectorsR)
        if (associated(results%eVectorsL)) deallocate(results%eVectorsL)
    end subroutine


    subroutine davidson_setup(convergence_threshold, &
                              max_subspace_size, &
                              max_iterations,    & 
                              unit_output,       & 
                              unit_sigmas,       & 
                              unit_subspaces,    &
                              unit_results,      &
                              solve_values,      &
                              solve_right,       &
                              solve_left,        &
                              symmetric,         &
                              complex_mode,      &
                              use_preconditioner,&
                              refresh_trial_rate,&
                              verbose,           & 
                              overlap_sorting,   &
                              save_results,      &
                              save_sigmas,       &
                              save_subspaces,    &
                              energy_shift,      &
                              print_config)

        real(kind=8), intent(in), optional :: convergence_threshold
        real(kind=8), intent(in), optional :: energy_shift(2)
        integer, intent(in), optional :: max_subspace_size 
        integer, intent(in), optional :: max_iterations
        integer, intent(in), optional :: unit_output
        integer, intent(in), optional :: unit_subspaces
        integer, intent(in), optional :: unit_sigmas
        integer, intent(in), optional :: unit_results
        integer, intent(in), optional :: refresh_trial_rate
        logical, optional :: solve_values
        logical, optional :: solve_right
        logical, optional :: solve_left
        logical, optional :: verbose
        logical, optional :: symmetric
        logical, optional :: overlap_sorting
        logical, optional :: save_sigmas 
        logical, optional :: save_subspaces
        logical, optional :: save_results
        logical, optional :: complex_mode
        logical, optional :: use_preconditioner
        logical, optional :: print_config

        call set_safe_minimum() 

        if (present(print_config)) then
            davidson_print_config = print_config
        else
            davidson_print_config = .true.
        end if
        if (present(convergence_threshold)) then
            davidson_convergence_threshold = convergence_threshold
        else
            davidson_convergence_threshold = 1.0d-8 
        end if

        if (present(energy_shift)) then
            davidson_eigenvalues_shift = energy_shift
            davidson_apply_eigenvalues_shift = .true.

            if (abs(davidson_eigenvalues_shift(1)).lt.davidson_eigenvalues_shift_threshold) &
               davidson_eigenvalues_shift(1) = 0.0d0 
 
            if (abs(davidson_eigenvalues_shift(2)).lt.davidson_eigenvalues_shift_threshold) &
               davidson_eigenvalues_shift(2) = 0.0d0 

            if ((abs(davidson_eigenvalues_shift(1)).lt.davidson_eigenvalues_shift_threshold).and. &
                (abs(davidson_eigenvalues_shift(2)).lt.davidson_eigenvalues_shift_threshold)) then
               davidson_apply_eigenvalues_shift = .false.
            end if
        else
            davidson_eigenvalues_shift = (/0.0d0, 0.0d0/) 
            davidson_apply_eigenvalues_shift = .false.
        end if

        if (present(max_subspace_size)) then
            davidson_max_subspace_size = max_subspace_size  
        else
            davidson_max_subspace_size = 1024
        end if

        if (present(max_iterations)) then
            davidson_max_iterations = max_iterations
        else
            davidson_max_iterations = 50
        end if

        if (present(solve_values)) then
            davidson_solve_values = solve_values
        else
            davidson_solve_values = .true. 
        end if

        if (present(solve_left)) then
            davidson_solve_lhs = solve_left 
        else
            davidson_solve_lhs = .false. 
        end if

        if (present(solve_right)) then
            davidson_solve_rhs = solve_right
        else
            davidson_solve_rhs = .true. 
        end if

        if (present(overlap_sorting)) then
            davidson_overlap_sorting = overlap_sorting
        else
            davidson_overlap_sorting = .false.
        end if

        if (present(use_preconditioner)) then
            davidson_use_preconditioner = use_preconditioner
        else
            davidson_use_preconditioner = .true.
        end if

        if (present(verbose)) then
            davidson_verbose = verbose
            davidson_debug   = verbose
        else
            davidson_verbose = .false.
        end if

        if (present(symmetric)) then
            davidson_symmetric = symmetric
        else
            davidson_symmetric = .false. 
        end if

        if (present(complex_mode)) then
            davidson_complex_mode_with_reals = complex_mode 
        else
            davidson_complex_mode_with_reals = .false. 
        end if

        if (present(refresh_trial_rate)) then
            if (refresh_trial_rate.le.0) then
                davidson_refresh_trial_vectors = .false.
            else 
                davidson_refresh_trial_vectors = .true.
                if (refresh_trial_rate.le.davidson_maximum_trial_vector_refresh_rate) then
                    davidson_trial_vector_refresh_rate = davidson_maximum_trial_vector_refresh_rate 
                else
                    davidson_trial_vector_refresh_rate = refresh_trial_rate
                end if
            end if
        else
            davidson_refresh_trial_vectors = .false.
        end if

        if (present(unit_output)) then
            davidson_output_file_unit = unit_output
        end if

        if (present(save_subspaces)) then
           davidson_save_subspaces = save_subspaces
        else
           davidson_save_subspaces = .false. 
        end if

        if (present(save_sigmas)) then
           davidson_save_sigmas = save_sigmas
        else
           davidson_save_sigmas = .false. 
        end if

        if (present(save_results)) then
           davidson_save_results = save_results
        else
           davidson_save_results = .false. 
        end if

        if (present(unit_subspaces)) then
            davidson_subspaces_file_unit = unit_subspaces
        else
            davidson_subspaces_file_unit = 137 
        end if

        if (present(unit_sigmas)) then
            davidson_sigmas_file_unit = unit_sigmas
        else
            davidson_sigmas_file_unit = 138
        end if

        if (present(unit_results)) then
            davidson_results_file_unit = unit_results
        else
            davidson_results_file_unit = 138
        end if

        if (davidson_print_config) call davidson_show_configurations()
    end subroutine


    subroutine davidson_show_configurations

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "Configuration variables for matrix-free diagonalizer"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  output written to unit        : ", davidson_output_file_unit
        write(davidson_output_file_unit, *) "  verbose output                : ", davidson_verbose

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) " convergence control"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  convergence threshold         : ", davidson_convergence_threshold
        write(davidson_output_file_unit, *) "  maximum subspace size         : ", davidson_max_subspace_size
        write(davidson_output_file_unit, *) "  maximum number of iterations  : ", davidson_max_iterations
        write(davidson_output_file_unit, *) "  refresh trial vectors         : ", davidson_refresh_trial_vectors 
        if (davidson_refresh_trial_vectors) then
            write(davidson_output_file_unit, '(A,I4,A)') &
                         "     set to refresh vectors every ",davidson_trial_vector_refresh_rate," interations"
        end if

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) " restart/data storage"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  save subspaces                :", davidson_save_subspaces
        if (davidson_save_subspaces) &
        write(davidson_output_file_unit, *) "  subspaces written to unit     : ", davidson_subspaces_file_unit
        write(davidson_output_file_unit, *) "  save sigma vectors            : ", davidson_save_sigmas
        if (davidson_save_sigmas) &
        write(davidson_output_file_unit, *) "  sigma vectors written to unit : ", davidson_sigmas_file_unit
        write(davidson_output_file_unit, *) "  save results                  :", davidson_save_results
        if (davidson_save_results) &
        write(davidson_output_file_unit, *) "  results written to unit       : ", davidson_results_file_unit

        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) " diagonalization characteristics"
        write(davidson_output_file_unit, *) ""
        write(davidson_output_file_unit, *) "  solve for right eigenvectors  : ", davidson_solve_rhs
        write(davidson_output_file_unit, *) "  solve for left eigenvectors   : ", davidson_solve_lhs
        write(davidson_output_file_unit, *) "  symmetric eigenproblem        : ", davidson_symmetric
        write(davidson_output_file_unit, *) "  complex mode (real variables) : ", davidson_complex_mode_with_reals
        write(davidson_output_file_unit, *) "  root following via overlap    : ", davidson_overlap_sorting
        write(davidson_output_file_unit, *) "  energy shift in eigval sorting: ",davidson_apply_eigenvalues_shift
!       if (davidson_apply_eigenvalues_shift) then
           write(davidson_output_file_unit, *) "     shift value, real  part    : ",davidson_eigenvalues_shift(1)
           write(davidson_output_file_unit, *) "     shift value, imag. part    : ",davidson_eigenvalues_shift(2)
!       end if

        write(davidson_output_file_unit, *) ""

    end subroutine


    subroutine set_safe_minimum()
        double precision dlamch
        safe_minimum = dlamch('S')
        if (davidson_debug) print *, "debug, safe minimum (LAPACK) is ",safe_minimum
    end subroutine 


    subroutine davidson_cleanup(results)
        type(davidson_results), intent(inout) :: results
        call free_davidson_results(results)
    end subroutine

    function updated_subspace_vectors_R(old_vectors_R, delta, orthogonalize_all,verbose)
        real(kind=8), pointer:: updated_subspace_vectors_R(:,:)

        real(kind=8), intent(inout), pointer :: old_vectors_R(:,:)
        real(kind=8), intent(in)      :: delta(:,:)
        logical, optional, intent(in) :: orthogonalize_all
        logical, optional, intent(in) :: verbose

        integer :: vectors_L, vectors_N
        integer :: vectors_Ld, vectors_Nd
        integer :: new_L
        integer :: i, ortho_from_vector
        logical :: orthogonalize_all_local
        logical :: verbose_local
        real(kind=8), allocatable :: c(:,:)

        vectors_N = size(old_vectors_R,1)
        vectors_L = size(old_vectors_R,2)

        vectors_Nd= size(delta,1)
        vectors_Ld= size(delta,2)

        if (present(orthogonalize_all)) then
            orthogonalize_all_local = orthogonalize_all 
        else
            orthogonalize_all_local = .false.
        end if

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        new_L = vectors_L + vectors_Ld
        allocate(updated_subspace_vectors_R(vectors_N, new_L))
        updated_subspace_vectors_R = 0.0d0

        do i = 1, vectors_L
            updated_subspace_vectors_R(:,i) = old_vectors_R(:,i)
        end do 
        deallocate(old_vectors_R)

        do i = 1, vectors_Ld  
           updated_subspace_vectors_R(:, vectors_L + i) = delta(:,i)
        end do

        if (verbose_local) &
            call davidson_print_vectors(updated_subspace_vectors_R, label="new subspace vectors before OR")

        if (orthogonalize_all_local) then
            call orthonormalize_subspace(updated_subspace_vectors_R, verbose=verbose_local)
        else ! only orthogonalize the new vectors with respect to the old ones
            ortho_from_vector=vectors_L + 1
            call orthonormalize_subspace(updated_subspace_vectors_R, start_from=ortho_from_vector, verbose=verbose_local)
        end if

    end function

    subroutine orthonormalize_subspace(b, start_from, method, verbose)
        real(kind=8), intent(inout)   :: b(:,:)
        integer, optional, intent(in) :: start_from
        logical, optional, intent(in) :: verbose
        integer, optional, intent(in) :: method

        integer :: method_local
        integer :: start_from_local
        logical :: verbose_local

        if (present(start_from)) then
            start_from_local = start_from  
        else
            start_from_local = 1 ! orthonormalize the full set of trial of vectors
        end if

        if (present(method)) then
            method_local = method 
        else
            method_local = 1 ! 1 = Gram-Schmidt procedure, 2 = QR
        end if

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        if (method_local .eq. 1) then
            call gram_schmidt_orthogonalization(b, start_from_local, verbose=verbose_local)
        else if (method_local .eq. 2) then
            call qr_orthogonalization(b, start_from_local)
        end if

    end subroutine


    subroutine qr_orthogonalization(b, start_from)
        real(kind=8), intent(inout)   :: b(:,:)
        integer, intent(in) :: start_from
        integer :: m, n, info, lwork, i, k, m_complex
        real(kind=8), allocatable :: work(:), tau(:), q(:,:)
        real(kind=8) :: norm 
        double precision ddot

        m = size(b,1)
        n = size(b,2)
        k = min(m,n)
 
        allocate(q(m,n))
        q = b
        stop ('orthogonalization via QR decomposition not fully debugged yet')
        if (davidson_complex_mode_with_reals) then
!           stop ('orthogonalization via QR decomposition not implemented yet for complex variables')
            m_complex = m / 2
            k = min(m_complex,n)
            allocate(tau(k))
 
            lwork = -1
            allocate(work(1))
            call dgeqrf(m, n, b, m, tau, work, lwork, info)
            lwork = 2*work(1)
            deallocate(work)

            allocate(work(lwork))
            if (davidson_debug) &
                call davidson_print_vectors(b, label="b before dgeqrf               ")
            call dgeqrf(m, n, b, m, tau, work, lwork, info)
            if (davidson_debug) &
                call davidson_print_vectors(b, label="b after  dgeqrf               ")
            deallocate(work)
 
            lwork = -1
            allocate(work(1))
            call dorgqr(m, n, k, b, m, tau, work, lwork, info)
            lwork = 2*work(1)
            deallocate(work)
            allocate(work(lwork))
            if (davidson_debug) &
                call davidson_print_vectors(b, label="b before dorgqr               ")
            call dorgqr(m, n, k, b, m, tau, work, lwork, info)
            if (davidson_debug) &
                call davidson_print_vectors(b, label="b after  dorgqr               ")
            deallocate(work)
            deallocate(tau)
 
        else
            allocate(tau(k))
            tau = 0.0d0
 
            lwork = -1
            allocate(work(1))
            work = 0.0d0
            call dgeqrf(m, n, b, m, tau, work, lwork, info)
            print *,"debug in dgeqrf, info #1:",info
            lwork = work(1)
            deallocate(work)
 
            allocate(work(2*lwork))
            work = 0.0d0
            if (davidson_debug) &
                call davidson_print_vectors(b, label="b before dgeqrf               ")
            call dgeqrf(m, n, b, m, tau, work, lwork, info)
            if (davidson_debug) then 
                call davidson_print_vectors(b, label="b after  dgeqrf               ")
                print *,"debug in dgeqrf, info #2:",info
            end if
            deallocate(work)
 
            lwork = -1
            allocate(work(1))
            work = 0.0d0
            call dorgqr(m, n, k, b, m, tau, work, lwork, info)
            print *,"debug in dgeqrf, info #3:",info
            lwork = work(1)
            deallocate(work)
 
            allocate(work(2*lwork))
            work = 0.0d0
            if (davidson_debug) &
            call davidson_print_vectors(b, label="b before dorgqr               ")
            call dorgqr(m, n, k, b, m, tau, work, lwork, info)
            if (davidson_debug) then 
                call davidson_print_vectors(b, label="b after  dorgqr               ")
                print *,"debug in dgeqrf, info #4:",info
            end if
            deallocate(work)
            deallocate(tau)
            if (davidson_debug) then
                do i = 1, n
                    norm = ddot(m, b(:,i), 1, b(:,i), 1)
                    print *, 'norm of new orthogonal vector ',i,':', norm
                end do
            end if
        end if

        do i = start_from, n
            b(:,i) = q(:,i)
        end do
        deallocate(q)
    end subroutine
 

    subroutine gram_schmidt_projection_R(projection, u, v)

!LV     NB: we assume that all vectors are already normalized
!LV     because normalizing here is order N^2 in the number of trial vectors

        real(kind=8), intent(inout) :: projection(:)
        real(kind=8), intent(in) :: u(:)
        real(kind=8), intent(in) :: v(:)
        real(kind=8) :: numerator, denominator, factor
        integer :: N, N_complex
        double precision ddot, dznrm2
        real(kind=8) :: factor_complex(2)

        N = size(u,1)

        projection = 0.0d0 

        if (davidson_complex_mode_with_reals) then
            N_complex = N / 2
            call xdotc(factor_complex,N_complex, u, 1, v, 1)
            projection = u
            call zscal(N_complex, factor_complex, projection, 1)
            if (davidson_debug) print *, "GS projection, factor: ",factor_complex

        else
            factor   = ddot(N, v, 1, u, 1)
            projection = u
            call dscal(N, factor, projection, 1)
            if (davidson_debug) print *, "GS projection, factor: ",factor
        end if

    end subroutine 

    subroutine gram_schmidt_orthogonalization(vectors, start_from, verbose)
        real(kind=8), intent(inout) :: vectors(:,:)
        integer, intent(in) :: start_from
        logical, optional, intent(in) :: verbose

        real(kind=8), allocatable   :: projection(:)
        integer :: vectors_L, vectors_N
        integer :: i, j, k, igs
        logical :: verbose_local
        real(kind=8) :: norm, vr, vi, pr, pi

        vectors_L = size(vectors,2)
        vectors_N = size(vectors,1)

        allocate(projection(vectors_N))

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false.
        end if

        do j = start_from, vectors_L
            norm = vector_norm(vectors(:,j))
            if (verbose_local) then
                if (davidson_debug) print *, "debug, Gram-Schmidt orthogonalization"
                print *, "  renormalizing vector v(",j,") : (norm = ",norm,")"
                print *, ""
            end if 
            if (norm .ge. safe_minimum ) then
                vectors(:,j) = vectors(:,j) / norm
            else
                stop ('division by zero in orthogonalization') 
            end if
            do igs = 1, 2 ! Gram-Schmidt needs to be done twice for numerical stability
               do i = 1, j-1
                   call gram_schmidt_projection_R(projection, vectors(:,i),vectors(:,j))
                   vectors(:,j) = vectors(:,j) - projection(:)
               end do
               norm = vector_norm(vectors(:,j))
               if (norm .ge. safe_minimum ) then
                   vectors(:,j) = vectors(:,j) / norm
               else
                   stop ('division by zero in orthogonalization')
               end if
            end do
        end do

        if (verbose_local) then
            call davidson_print_vectors(vectors, label="all vectors                   ")
        end if

        deallocate(projection)
    end subroutine 

    subroutine davidson_print_vectors(vectors, label, nvectors, complex_dimension)
        real(kind=8), intent(in) :: vectors(:,:)
        integer, intent(in), optional :: nvectors
        character(len=max_label_length), intent(in), optional :: label 
        integer, intent(in), optional :: complex_dimension

        integer :: L, N, i, j, N_complex, L_complex, N_step, L_step, kr, ki
        real(kind=8) :: vr, vi
! aspg, for debugging, set this threshold to eliminate small values in printout
!       for both real (and imaginary) parts
        real(kind=8) :: print_vector_threshold = 1.0d-10 
        character(len=max_label_length) :: label_local

        N = size(vectors,1)
        L = size(vectors,2)
        if (present(nvectors)) L = nvectors

        if (present(label)) then
            label_local = label
            print *, ""
            print "(a)", label_local(1:max_label_length)
            print *, ""
            if (davidson_debug) print *, 'debug: size 1 ',N, '  size 2 ',L
            print *, ""
        end if
        if (present(complex_dimension)) then
            if (complex_dimension.eq.1) then
                N_step = 2
                L_step = 1
            else if (complex_dimension.eq.2) then
                N_step = 1
                L_step = 2
            end if
        else
           N_step = 2
           L_step = 1
        end if

        if (davidson_complex_mode_with_reals) then
            N_complex = N / N_step
            L_complex = L / L_step

            do i = 1, N_complex
                print "(2x,i6,$)", i
                do j = 1, L_complex
                    if (N_step.eq.2) then
                        kr = N_step*(i-1) + 1
                        ki = kr + 1
                        vr = vectors(kr,j)
                        vi = vectors(ki,j)
                    else 
                        kr = L_step*(j-1) + 1
                        ki = kr + 1
                        vr = vectors(i, kr)
                        vi = vectors(i, ki)
                    end if
                    if (abs(vi).lt.print_vector_threshold) vi = 0.0d0
                    if (abs(vr).lt.print_vector_threshold) vr = 0.0d0
                    if (vi < 0.0d0) then
                        print "(4x,1e14.6,a,1e13.6,a,$)", vr," -",abs(vi)," i"
                    else
                        print "(4x,1e14.6,a,1e13.6,a,$)", vr," +",abs(vi)," i"
                    end if
                end do
                print *, ""
            end do
        else
            do i = 1, N
                print "(2x,i6,$)", i
                do j = 1, L
                    vr = vectors(i,j)
                    if (abs(vr).lt.print_vector_threshold) vr = 0.0d0
                    print "(4x,1e14.6,$)", vr
                end do
                print *, ""
            end do
        end if 
    end subroutine


    subroutine print_values_vectors_R(output_file_unit, eValues, eVectorsR, eVectorsL, left, right, nroots)
        real(kind=8), intent(in) :: eVectorsL(:,:)
        real(kind=8), intent(in) :: eVectorsR(:,:)
        real(kind=8), intent(in) :: eValues(:,:)
        integer, optional, intent(in) :: nroots
        integer, intent(in) :: output_file_unit

        logical :: left, right
        integer :: vectors_N
        integer :: nroots_local
        integer :: i
        character(len=max_label_length) :: label
        real(kind=8) :: vr, vi
! aspg, for debugging, set this threshold to eliminate small values in printout
!       for both real (and imaginary) parts
        real(kind=8) :: print_values_threshold = 1.0d-10


        vectors_N = size(eValues,1)
        nroots_local = vectors_N
        if (present(nroots)) nroots_local = nroots

        label = "Eigenvalues"
        if (davidson_debug) print *, "debug, evalues inside print_values_vectors_R"
        if (davidson_debug) print *, eValues(1,:)
        if (davidson_debug) print *, eValues(2,:)

        do i = 1, nroots_local 
            vr = eValues(1,i)
            vi = eValues(2,i)
            if (abs(vi).lt.print_values_threshold) vi = 0.0d0
            if (abs(vr).lt.print_values_threshold) vr = 0.0d0

            if (vi .lt. 0.0d0) then
                write (davidson_output_file_unit,'(4x,a,i6,a,1e14.6,a,1e14.6,a)')  &
                      "Eigenvalue ",i," : ", vr," -",abs(vi)," i"
            else if (vi .gt. 0.0d0) then
                write (davidson_output_file_unit,'(4x,a,i6,a,1e14.6,a,1e14.6,a)')  &
                      "Eigenvalue ",i," : ", vr," +",abs(vi)," i"
            else
                write (davidson_output_file_unit,'(4x,a,i6,a,1e14.6)')  &
                      "Eigenvalue ",i," : ", vr
            end if
        end do
        if (left) then
            label = "Left Eigenvectors            "
            call davidson_print_vectors(eVectorsL, label=label) 
        end if
        if (right) then
            label = "Right Eigenvectors           "
            call davidson_print_vectors(eVectorsR, label=label) 
        end if

    end subroutine 

    function vector_norm(v)
        real(kind=8) :: vector_norm
        real(kind=8), intent(in) :: v(:)
        integer :: vector_N, i, vector_N_complex
        double precision dnrm2, dznrm2

        vector_N = size(v,1)
        if (davidson_complex_mode_with_reals) then
            vector_N_complex = vector_N / 2
            vector_norm = dznrm2(vector_N_complex, v, 1)
        else
            vector_norm = dnrm2(vector_N, v, 1)
        end if
    end function

    

    function correction_to_vectors(output_file_unit, nroots, e_small_M, c_small_M, sigma_vector, b, preconditioner, correct_all)
        real(kind=8), pointer     :: correction_to_vectors(:,:)
        integer, intent(in)       :: nroots, output_file_unit
        real(kind=8), intent(in)  :: e_small_M(:,:)
        real(kind=8), intent(in)  :: c_small_M(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        real(kind=8), pointer :: preconditioner(:,:)
        logical, optional, intent(in) :: correct_all

        integer :: vectors_L, vectors_N
        integer :: k, unconverged_roots
        real(kind=8) :: norm_r = 0.0d0
        real(kind=8), allocatable :: r(:,:)
        logical, allocatable :: vectors_to_update(:)
        logical :: correct_all_local

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (present(correct_all)) then
            correct_all_local = correct_all
        else
            correct_all_local = .false. 
        end if

        allocate(vectors_to_update(nroots))
        vectors_to_update = .true.

        allocate(r(vectors_N,nroots))
        r = 0.0d0
        call form_r(r, nroots, e_small_M, c_small_M, sigma_vector, b)

        unconverged_roots = 0
        do k = 1, nroots
            norm_r = vector_norm(r(:,k))
            if (norm_r .lt. davidson_convergence_threshold) then 
               write(output_file_unit, "(4x,a,i6,a,1e14.6,4x,a)") "Root ",k,"  ||q_k,M|| ", norm_r, " <<< converged ! "
               vectors_to_update(k) = .false.
            else
               write(output_file_unit, "(4x,a,i6,a,1e14.6)")  "Root ",k,"  ||q_k,M|| ", norm_r
               vectors_to_update(k) = .true.
               unconverged_roots = unconverged_roots + 1
            end if
        end do
        write (davidson_output_file_unit,*) " "

        if (any(vectors_to_update)) then
            if (correct_all_local) then
                unconverged_roots = nroots
                vectors_to_update = .true.
            end if
            allocate(correction_to_vectors(vectors_N,unconverged_roots))
            correction_to_vectors = 0.0d0
            call form_delta(correction_to_vectors, vectors_to_update, r, preconditioner)
        else
            correction_to_vectors => NULL()
        end if

        deallocate(vectors_to_update)
        deallocate(r)
    end function

    function form_preconditioner(M_diag, eigenvalues)
        real(kind=8), pointer     :: form_preconditioner(:,:)
        real(kind=8), intent(in)  :: eigenvalues(:,:)
        real(kind=8), intent(in)  :: M_diag(:)

        integer :: i, k, l, j, r
        integer :: vectors_N
        integer :: nroots, offset
        real(kind=8) :: x2py2

! careful with dimension 1 of eigenvalues, 
        nroots    = size(eigenvalues,2)
        vectors_N = size(M_diag,1)

        allocate(form_preconditioner(vectors_N,nroots))
        form_preconditioner = 0.0d0
        if (davidson_debug) print *, 'debug, vectors_N in form_preconditioner :',vectors_N

        if (.not.davidson_use_preconditioner) then
            if (davidson_complex_mode_with_reals) then
                offset = 2
            else
                offset = 1
            end if
            do k = 1, nroots
                do j = 1, vectors_N, offset
                    form_preconditioner(j,k) = 1.0d0
                end do
            end do
        else
            if (davidson_complex_mode_with_reals) then
                do k = 1, nroots
                   r = 0
                    do j = 1, vectors_N, 2
                        r = r + 1
!                       if (davidson_debug) print *, 'debug, H_II(re,im) : ',j,M_diag(j),M_diag(j+1)
!                        form_preconditioner(j,k)   = eigenvalues(1,k) - M_diag(r)
! the diagonal appears not to have the complex part, so for debugging purposes i will set it to zero here
! and only construct the real part
 !                      form_preconditioner(j,k)   = eigenvalues(1,k) - M_diag(j)
 !                      form_preconditioner(j+1,k) = eigenvalues(2,k) - M_diag(j+1)
                        
                        form_preconditioner(j,k)   = eigenvalues(1,k) - M_diag(j)
                        form_preconditioner(j+1,k) = 0.0d0

! using the formula for calculating the reciprocal of the nonzero complex number z = x + yi 
!
!{\displaystyle {\frac {1}{z}}={\frac {\bar {z}}{z{\bar {z}}}}
!                             ={\frac {\bar {z}}{x^{2}+y^{2}}}
!                             ={\frac {x}{x^{2}+y^{2}}}-{\frac {y}{x^{2}+y^{2}}}i.} 
!
                        x2py2 = form_preconditioner(j,k)*form_preconditioner(j,k)  &
                              + form_preconditioner(j+1,k)*form_preconditioner(j+1,k)

                        if (dabs(x2py2) .gt. safe_minimum) then
                            form_preconditioner(j,k)   =  form_preconditioner(j,k) / x2py2 
                            form_preconditioner(j+1,k) = -form_preconditioner(j+1,k) / x2py2
                        else
!                           stop ('detected division by zero')
                            write (davidson_output_file_unit,*) 'Warning : numerator is small', &
                                                             (1.0d0 / form_preconditioner(j,k)),& 
&                                    ' element of preconditioner set to 1'
                            form_preconditioner(j,k) = 1.0d0
                            form_preconditioner(j+1,k) = 0.0d0
                        end if
                    end do
                end do
            else
                do k = 1, nroots
                    form_preconditioner(:,k) = eigenvalues(1,k) - M_diag(:)
                    do l = 1, vectors_N
                        if (dabs(form_preconditioner(l,k)) .gt. safe_minimum) then
                            form_preconditioner(l,k) = 1.0d0 / form_preconditioner(l,k)
                        else
!                           stop ('detected division by zero')
                            x2py2=1.0d0 / form_preconditioner(l,k)
                            write (davidson_output_file_unit,*) &
                          'warning : numerator is smaller than', safe_minimum, &
     &                               '(',x2py2,'). element of preconditioner set to 1.'
                            form_preconditioner(l,k) = 1.0d0 
                        end if
                    end do
                end do
            end if
        end if ! use preconditioner
        if (davidson_debug) call davidson_print_vectors(form_preconditioner, label="preconditioner                ")
    end function

    subroutine form_delta(delta, vectors_to_update, r, preconditioner)
        real(kind=8), intent(inout) :: delta(:,:)
        logical, intent(in)         :: vectors_to_update(:)
        real(kind=8), intent(in)    :: preconditioner(:,:)
        real(kind=8), intent(in)    :: r(:,:)

        integer :: nroots
        integer :: i, k, l
        real(kind=8) :: scaling_factor(2), vr, vi
        integer :: vectors_N, vectors_N_complex

        nroots    = size(vectors_to_update,1)

        if (davidson_complex_mode_with_reals) then
            vectors_N = size(delta,1) 
            vectors_N_complex = vectors_N / 2
            scaling_factor(2) = 0.0d0 
 
!!          if (davidson_debug) write (davidson_output_file_unit,*) 'debug, vectors to update',vectors_to_update
            i = 1
            do k = 1, nroots
               if (vectors_to_update(k)) then
                   do l = 1, vectors_N, 2
                       vr = preconditioner(l,k)*r(l,k)   - preconditioner(l+1,k)*r(l+1,k)
                       vi = preconditioner(l,k)*r(l+1,k) + preconditioner(l+1,k)*r(l,k)
                       delta(l,i) = vr
                       delta(l+1,i) = vi
                   end do
                   scaling_factor(1) = 1.0d0 / vector_norm(delta(:,i))
                   call zscal(vectors_N_complex,scaling_factor,delta(:,i),1)
                   i = i + 1
               end if
            end do
        else
            i = 1
            do k = 1, nroots
                if (vectors_to_update(k)) then
                   delta(:,i) = preconditioner(:,k)*r(:,k)
                   delta(:,i) = delta(:,i) / vector_norm(delta(:,i))
                   i = i + 1
               end if
            end do
        end if

! aspg: the call below just fails.. don't quite get why
!       if (davidson_debug) call davidson_print_vectors(delta,          label="deltas                        ")
!       if (davidson_debug) call davidson_print_vectors(preconditioner, label="preconditioner                ")

    end subroutine 

    subroutine form_r(r, nroots,  e_small_M, c_small_M, sigma_vector,  b)
        real(kind=8), intent(inout) :: r(:,:)
        integer, intent(in)       :: nroots
        real(kind=8), intent(in)  :: e_small_M(:,:)
        real(kind=8), intent(in)  :: c_small_M(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        integer :: vectors_L, vectors_N
        integer :: k
        real(kind=8), allocatable :: residue_a(:,:)
        real(kind=8), allocatable :: residue_b(:,:)

        vectors_L = size(b,2)
        vectors_N = size(b,1)

        allocate(residue_a(vectors_N,vectors_L))
        allocate(residue_b(vectors_N,vectors_L))

        residue_a = 0.0d0
        residue_b = 0.0d0

        if (davidson_complex_mode_with_reals) then
            vectors_N = vectors_N / 2
            call zgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                       (1.0d0, 0.0d0), sigma_vector, vectors_N, c_small_M, vectors_L, &
                       (0.0d0, 0.0d0), residue_a, vectors_N)
            call zgemm('N', 'N', vectors_N, vectors_L, vectors_L,                   &
                       (1.0d0, 0.0d0), b, vectors_N, c_small_M, vectors_L,            &
                       (0.0d0, 0.0d0), residue_b, vectors_N)
        else
            call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,              &
                       (1.0d0), sigma_vector, vectors_N, c_small_M, vectors_L, &
                       (0.0d0), residue_a, vectors_N)
            call dgemm('N', 'N', vectors_N, vectors_L, vectors_L,              &
                       (1.0d0), b, vectors_N, c_small_M, vectors_L,            &
                       (0.0d0), residue_b, vectors_N)
        end if

!       if (davidson_debug) then
!           call davidson_print_vectors(residue_a,  label="debug: in form_r, residue a   ", complex_dimension=1)
!           call davidson_print_vectors(residue_b,  label="debug: in form_r, residue b   ", complex_dimension=1)
!           call davidson_print_vectors(e_small_M,  label="debug: in form_r, e_small M   ", complex_dimension=1)
!       end if
! todo: generalize this to the multiplication of the complex e_small with  complex residue_b

        do k = 1, nroots
            if (davidson_complex_mode_with_reals) then
                call zscal(vectors_N,e_small_M(:,k),residue_b(:,k),1)
            else
                call dscal(vectors_N,e_small_M(1,k),residue_b(:,k),1)
            end if
            r(:,k) = residue_a(:,k) - residue_b(:,k)
        end do
!       if (davidson_debug) then
!           call davidson_print_vectors(r,  label="debug: in form_r, r           ", complex_dimension=1)
!       end if

        deallocate(residue_a)
        deallocate(residue_b)
    end subroutine

    function create_reduced_problem(b, sigma_vector)
        real(kind=8), pointer     :: create_reduced_problem(:,:)
        real(kind=8), intent(in)  :: sigma_vector(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        real(kind=8), pointer     :: sigmat(:,:), bt(:,:)

        integer :: vectors_L, vectors_N
        real(kind=8) :: r_alpha_lapack = 1.0d0
        real(kind=8) :: r_beta_lapack = 1.0d0 

        integer :: vectors_L_complex, vectors_N_complex

        complex*16, pointer :: c_b(:,:), c_sigma(:,:), c_reduced(:,:)
        integer :: i, j, k

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (davidson_complex_mode_with_reals) then
            vectors_L_complex = 2*vectors_L
            vectors_N_complex = vectors_N / 2

            allocate(create_reduced_problem(vectors_L_complex,vectors_L))
            create_reduced_problem = 0.0d0

            call zgemm('C', 'N', vectors_L, vectors_L, vectors_N_complex, (1.d0, 0.0d0), b, vectors_N_complex, &
                       sigma_vector, vectors_N_complex, (0.d0, 0.0d0), create_reduced_problem, vectors_L)

            if (davidson_debug) &
            call davidson_print_vectors(create_reduced_problem,  label="debug: reduced problem cmplx 1", complex_dimension=1)
        else
            allocate(create_reduced_problem(vectors_L,vectors_L))
            create_reduced_problem = 0.0d0
            call dgemm('T', 'N', vectors_L, vectors_L, vectors_N, r_alpha_lapack, b, vectors_N, &
                       sigma_vector, vectors_N, r_beta_lapack, create_reduced_problem, vectors_L)
        end if

    end function

    subroutine davidson_transpose_conjugate_vector(v, vt, complex_dimension_v)
        real(kind=8), intent(in)    :: v(:,:)
        real(kind=8), intent(inout) :: vt(:,:)
        integer, intent(in) :: complex_dimension_v
        integer :: vectors_L, vectors_N
        real(kind=8) :: vr, vi, vtr, vti

        integer :: i, j, kr, ki, lr, li, N_step, L_step, vectors_N_complex, vectors_L_complex

        vectors_N = size(v,1)
        vectors_L = size(v,2)

        if (complex_dimension_v.eq.1) then
            N_step = 2 
            L_step = 1
        else if (complex_dimension_v.eq.2) then
            N_step = 1 
            L_step = 2
        end if
        vectors_N_complex = vectors_N / N_step 
        vectors_L_complex = vectors_L / L_step

        do i = 1, vectors_N_complex 
            do j = 1, vectors_L_complex
                if (complex_dimension_v.eq.1) then
                    kr = N_step*(i-1) + 1
                    ki = kr + 1

                    vr = v(kr,j)
                    vi = v(ki,j)
                    write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                                 'i:',i,' j:',j,' kr',kr,' ki:',ki,'  vr:',vr,'  vi:',vi 

                    vt(j,kr) =  vr
                    vt(j,ki) = -vi

                     write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1e14.6,a,1e14.6,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           ' v(',i,',',j,')=(',vr,',',vi,') ->  vt(',j,',',i,')=(',vt(j,kr),',',vt(j,ki),')'
                else if (complex_dimension_v.eq.2) then

                    kr = L_step*(j-1) + 1
                    ki = kr + 1

                    vr = v(i, kr)
                    vi = v(i, ki)
                     write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           'i:',i,' j:',j,' kr',kr,' ki:',ki,'  vr:',vr,'  vi:',vi

                    vt(kr, i) =  vr
                    vt(ki, i) = -vi

                     write (davidson_output_file_unit,'(a,1i4,a,1i4,a,1e14.6,a,1e14.6,a,1i4,a,1i4,a,1e14.6,a,1e14.6)') &
                           ' v(',i,',',j,')=(',vr,',',vi,') ->  vt(',j,',',i,')=(',vt(kr,i),',',vt(ki,i),')'

                end if
            end do
        end do

        if (davidson_debug) call davidson_print_vectors(v,  label="debug:  original vector v     ", complex_dimension=1)
!       print *, v
        if (davidson_debug) call davidson_print_vectors(vt, label="debug:  conj. transpose of v  ", complex_dimension=2)
!       print *, vt

    end subroutine

    subroutine solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left)
        real(kind=8), intent(inout), pointer :: M(:,:)
        real(kind=8), intent(inout), pointer :: eVectorsL(:,:)
        real(kind=8), intent(inout), pointer :: eVectorsR(:,:)
        real(kind=8), intent(inout), pointer :: eValues(:,:)
        logical, intent(in) :: left 

        integer :: vectors_N, vectors_L, i

        integer :: lwork, info
        character(len=1) :: do_left = 'N', do_right = 'V'
        real(kind=8), allocatable :: work(:)
        real(kind=8), allocatable :: rwork(:)

        vectors_N = size(M,1)
        vectors_L = size(M,2)

        if (associated(eVectorsL)) deallocate(eVectorsL)
        if (associated(eVectorsR)) deallocate(eVectorsR)
        if (associated(eValues))   deallocate(eValues)

        if (davidson_debug) call davidson_print_vectors(M, label="debug:  small matrix M #2     ")

        if (left) then
            do_left = 'V'
!           do_right = 'N'
            allocate(eVectorsL(vectors_N,vectors_L))
            eVectorsL = 0.0d0
        else
            allocate(eVectorsL(1,vectors_L))
            eVectorsL = 0.0d0
        end if
        allocate(eVectorsR(vectors_N,vectors_L))
        allocate(eValues(2,vectors_L))
        eVectorsR = 0.0d0
        eValues   = 0.0d0

! first query optimal lwork value
        lwork = -1
        allocate(work(4*vectors_L))
!       print *, "debug, lwork in diagonalization before query:",lwork,4*vectors_L

        if (davidson_complex_mode_with_reals) then
            allocate(rwork(2*vectors_N))
            call zgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                        eValues, eVectorsL, vectors_L, eVectorsR, &
                        vectors_L, WORK, LWORK, rwork, INFO)
            deallocate(rwork)
        else
            call dgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                        eValues(1,:), eValues(2,:), eVectorsL, vectors_L, eVectorsR, &
                        vectors_L, WORK, LWORK, INFO)
        end if

        lwork = work(1)
!       print *, "debug, lwork in diagonalization after query:",lwork,work(1:4)
        deallocate(work)
        allocate(work(2*lwork))
        if (davidson_complex_mode_with_reals) then
            allocate(rwork(2*vectors_N))
            call zgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                        eValues, eVectorsL, vectors_L, eVectorsR, &
                        vectors_L, WORK, LWORK, rwork, INFO)
            deallocate(rwork)
        else
            call dgeev (do_left, do_right, vectors_L, M, vectors_L,      &
                        eValues(1,:), eValues(2,:), eVectorsL, vectors_L, eVectorsR, &
                        vectors_L, WORK, LWORK, INFO)
            call check_for_complex_values_vectors(left, vectors_N, vectors_L, eValues, eVectorsL, eVectorsR)
        end if
!       print *, 'in small_solve_eigenvalue: info',info
        if (info .gt. 0) then
            print *, 'in small_solve_eigenvalue: not all eigenvalues converged' 
            print *, 'converged'
            do i = info+1, vectors_L 
                print *, " real part : ",eValues(1,i), " imaginary part : ",eValues(2,i)
            end do
        end if
        deallocate(work)

        call renormalize_evectors(left, vectors_L, eVectorsR, eVectorsL)

    end subroutine

! for the moment this routine only checks whether we found a complex conjugate pair
! of eigenvalues and prints a warning. still to thing about how to handle the complex
! eigenvectors if we are in real algebra.
    subroutine check_for_complex_values_vectors(left,nrows,nvectors,eV,vL,vR)
        logical, intent(in) :: left
        integer, intent(in) :: nrows, nvectors
        real(kind=8), intent(inout) :: vR(:,:), vL(:,:), eV(:,:)
        integer :: i, j
        real(kind=8) :: norm, tolerance = 1.0d-10
        real(kind=8), allocatable :: complex_eigenvectors(:,:)
        integer, allocatable :: complex_eigenvalue_pairs(:)

        allocate(complex_eigenvalue_pairs(nvectors))
        complex_eigenvalue_pairs = 0

!       allocate(complex_eigenvectors(nrows,2))
!       complex_eigenvectors = 0.0d0

! expected results from dgeev: 
! If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
! If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, then 
!   v(j)   = VR(:,j) + i*VR(:,j+1) 
!   v(j+1) = VR(:,j) - i*VR(:,j+1).

        do i = 1, nvectors - 1
           j = i + 1
           if ((dabs(eV(2,i)).gt.tolerance).and.(complex_eigenvalue_pairs(j).ne.0)) then
              complex_eigenvalue_pairs(i) = j
              complex_eigenvalue_pairs(j) = i
              write (davidson_output_file_unit,'(a,1i4,a,1i4,a,2e14.6,a,2e14.6,a)') &
              "Warning! Eigenvalues",i," and ",j," form c.c. pair: (", &
              eV(1,i),eV(2,i),"), (",eV(1,j),eV(2,j),")"
           end if
        end do

        deallocate(complex_eigenvalue_pairs)
    end subroutine 
 
    subroutine renormalize_evectors(left, nvectors, vR, vL)
        logical, intent(in) :: left
        integer, intent(in) :: nvectors
        real(kind=8), intent(inout) :: vR(:,:), vL(:,:) 
        integer :: vectors_L, i
        real(kind=8) :: norm, tolerance = 1.0d-10

        do i = 1, nvectors
           if (left) then
              norm = vector_norm(vL(:,i))
              if (dabs(norm - 1.0d0).gt.tolerance) then
                 write (davidson_output_file_unit,'(a,1i4,a,1e14.6,a)') &
                 "Warning! LHS subspace eigenvector",i," not normalized(",norm,"), normalizing it."
                 if (norm.gt.safe_minimum) norm = 1.0d0/norm 
                 vL(:,i) = norm*vL(:,i)
              end if
           else
              norm = vector_norm(vR(:,i))
              if (dabs(norm - 1.0d0).gt.tolerance) then
                 write (davidson_output_file_unit,'(a,1i4,a,1e14.6,a)') &
                 "Warning! RHS subspace eigenvector",i," not normalized(",norm,"), normalizing it."
                 if (norm.gt.safe_minimum) norm = 1.0d0/norm
                 vR(:,i) = norm*vR(:,i)
              end if
           end if
        end do
    end subroutine

    subroutine get_values_vectors_davidson(nroots, A, results, b_guess, state_sym, projector)
! bad coupling between modules, to be replaced by an abstract class that the sigma vector module extends
! and that will have associated get_matrix_size, create_sigma_vector and create_diagonal
        use sigma_eom, only : intermediates, get_matrix_size, create_diagonal
        type(intermediates), intent(in)   :: A
! back to general stuff
        type(davidson_results), intent(inout) :: results
        real(kind=8), intent(inout), pointer :: b_guess(:,:)
        logical, intent(in), pointer :: projector(:)
        integer, intent(in) :: nroots
        integer, intent(in) :: state_sym
!
! nullifying all pointers at the start to give them a state, see discussion on http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html
! "[...] When a pointer is declared its status is undefined, and cannot be safely queried with the associated intrinsic. [...]"
! 
        real(kind=8), pointer :: M_diag(:)           => NULL()
        real(kind=8), pointer :: delta(:,:)          => NULL()
        real(kind=8), pointer :: preconditioner(:,:) => NULL()
        real(kind=8), pointer :: sigma_vector(:,:)   => NULL()
        real(kind=8), pointer :: b(:,:)              => NULL()
        real(kind=8), pointer :: M_small(:,:)        => NULL()
        real(kind=8), pointer :: c_small_r(:,:)      => NULL()
        real(kind=8), pointer :: c_small_l(:,:)      => NULL()
        real(kind=8), pointer :: e_small(:,:)        => NULL()

        real(kind=8), pointer :: c_small_r_old(:,:)  => NULL()
        real(kind=8), pointer :: c_small_l_old(:,:)  => NULL()

        integer :: vectors_N, L, min_subspace_size, max_subspace_size, niter
        integer :: i, j, k

        vectors_N = get_matrix_size(A)
!
! davidson_create_results() will allocate storage for results from davidson diagonalization when called a first time. 
! 
        call davidson_create_results(results,vectors_N,nroots,&
                                     davidson_solve_values,   &
                                     davidson_solve_lhs,      &
                                     davidson_solve_rhs)

! todo: replace lines below with a call to a function which also performs the allocation, like
!        M_diag => create_diagonal(A)
        allocate(M_diag(size(A%H_II))) 
        M_diag = A%H_II 

        min_subspace_size = nroots 
        max_subspace_size = max(nroots, davidson_max_subspace_size)
        max_subspace_size = min(vectors_N, max_subspace_size)

        if (associated(projector)) then
            write(davidson_output_file_unit,*) ""
            write(davidson_output_file_unit,"(4x,a)") "Projector (for CVS/REW/... calculations) is active"
            write(davidson_output_file_unit,*) ""
        end if

        write(davidson_output_file_unit,"(4x,a,i6,a)") "Requested ",nroots," eigenvales"
        if ( nroots .gt. vectors_N ) then 
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested larger than matrix order !"
           return 
        else if ( nroots .eq. vectors_N ) then
           write(davidson_output_file_unit, "(4x,a)") "Number of roots requested equal to matrix order !"
        endif

        b => b_guess

        if (davidson_verbose) then
            write(davidson_output_file_unit,"(4x,a,i4,a,i8)") &
                  "Start by mapping M onto a subspace of size ",min_subspace_size, &
                  ", will extend it up to ",max_subspace_size
            call davidson_print_vectors(b, label="Initial subspace vectors      ")
        end if

       if (davidson_overlap_sorting) &
           call setup_overlap_sorting(c_small_r_old, c_small_l_old, nroots, davidson_solve_rhs, davidson_solve_lhs)


        write(davidson_output_file_unit,*) ""
        do niter = 1, davidson_max_iterations
            write(davidson_output_file_unit,*) "Iteration ", niter
            write(davidson_output_file_unit,*) ""

            if (associated(M_small)) deallocate(M_small)

            if (associated(projector)) call apply_projector(projector,b)
            if (davidson_debug) call davidson_print_vectors(b, label="debug:  b     vector          ")

            sigma_vector => create_sigma_vector(A, b, sigma_vector, state_sym)
            if (associated(projector)) call apply_projector(projector,sigma_vector)

            if (davidson_debug) call davidson_print_vectors(sigma_vector, label="debug:  sigma vector          ")

            if (davidson_solve_rhs) M_small => create_reduced_problem(b, sigma_vector)
            if (davidson_solve_lhs) M_small => create_reduced_problem(sigma_vector, b)

            if (davidson_debug) call davidson_print_vectors(M_small, label="debug:  small matrix M        ")
            if (davidson_debug) call dump_to_file(2222, "M_small__.bin", M_small)

            if (davidson_debug) print *,"niter, follow overlap",niter,(davidson_overlap_sorting.and.(niter.gt.1))

            call calculate_eigenvalues_and_eigenvectors(davidson_output_file_unit,    &
                                                        M_small,                      &
                                                        e_small,                      &
                                                        c_small_r, c_small_l,         &
                                                        c_small_r_old, c_small_l_old, &
                                                        davidson_eigenvalues_shift,   &
                                                        follow_overlap=(davidson_overlap_sorting.and.(niter.gt.1)),  &
                                                        apply_evalues_shift=davidson_apply_eigenvalues_shift, &
                                                        left=davidson_solve_lhs,      &
                                                        right=davidson_solve_rhs,     &
                                                        nroots=nroots,                &
                                                        verbose=davidson_verbose)

            if (niter == davidson_max_iterations) exit ! should not extend space anymore because this was the final iteration

            if (davidson_solve_rhs.and.davidson_overlap_sorting) then
                if (davidson_complex_mode_with_reals) then
                    c_small_r_old = c_small_r(1:2*nroots,1:nroots)
                else
                    c_small_r_old = c_small_r(1:nroots,1:nroots)
                end if
            end if
            if (davidson_solve_lhs.and.davidson_overlap_sorting) then
                if (davidson_complex_mode_with_reals) then
                    c_small_l_old = c_small_l(1:2*nroots,1:nroots)
                else
                    c_small_l_old = c_small_l(1:nroots,1:nroots)
                end if
            end if

            if (associated(preconditioner)) deallocate(preconditioner)
            if (associated(delta)) deallocate(delta)

            preconditioner =>  form_preconditioner(M_diag, e_small)

           if (davidson_solve_rhs) delta => correction_to_vectors(davidson_output_file_unit, nroots,  &
                                           e_small, c_small_r, sigma_vector, & 
                                           b, preconditioner)

           if (davidson_solve_lhs) delta => correction_to_vectors(davidson_output_file_unit, nroots,  &
                                           e_small, c_small_l, sigma_vector, & 
                                           b, preconditioner)

            if (.not.associated(delta)) then
                write(davidson_output_file_unit,*) "All roots converged."
                 if (associated(preconditioner)) deallocate(preconditioner)
                exit
            end if

! not converged yet, if asked to refresh the trial vectors do it now, otherwise create
! new trial vectors
            if (davidson_refresh_trial_vectors.and.(mod(niter,davidson_trial_vector_refresh_rate).eq.0)) then
                write(davidson_output_file_unit,*) "Resetting trial vectors at iteration",niter
                call form_eigenvalues_full_problem(results%eValues, nroots, e_small)
                if (davidson_solve_lhs) then
                   write(davidson_output_file_unit,*) "    forming left  eigenvectors for original problem"
                   call form_eigenvector_full_problem(results%eVectorsL, nroots,b,c_small_l,davidson_solve_rhs,davidson_solve_lhs)
                   if (davidson_debug) &
                       call davidson_print_vectors(results%eVectorsL, label="intm  eigenvalues full problem")
                   deallocate(b)
                   allocate(b(vectors_N,nroots))
                   write(davidson_output_file_unit,*) "    copying left  eigenvectors to trial vectors"
                   write(davidson_output_file_unit,*) " "
                   b = results%eVectorsL
                end if
                if (davidson_solve_rhs) then
                   write(davidson_output_file_unit,*) " "
                   write(davidson_output_file_unit,*) "    forming right eigenvectors for original problem"
                   call form_eigenvector_full_problem(results%eVectorsR, nroots,b,c_small_r,davidson_solve_rhs,davidson_solve_lhs)
                   if (davidson_debug) &
                       call davidson_print_vectors(results%eVectorsR, label="intm  eigenvalues full problem")
                   deallocate(b)
                   allocate(b(vectors_N,nroots))
                   write(davidson_output_file_unit,*) "    copying right eigenvectors to trial vectors"
                   write(davidson_output_file_unit,*) " "
                   b = results%eVectorsR
                end if
                deallocate(sigma_vector)
            else

               b => updated_subspace_vectors_R(b, delta, verbose=davidson_debug)

! in the case of a parallel run, at each davidson step synchronize the new trial vectors 
               call davidson_synchronize_trial_vectors(b,size(delta,2))

               if (davidson_verbose) then
                   call davidson_print_vectors(delta, label="delta                         ")
                   call davidson_print_vectors(b, label="Updated subspace vectors      ")
               end if
            end if
 
! now check whether we reached the maximum subspace size
            L = size(b,2)
            if (L .eq. max_subspace_size) then
               if (associated(preconditioner)) deallocate(preconditioner)
               if (associated(delta)) deallocate(delta)
               write(davidson_output_file_unit, *) "Reached maximum subspace size, stopping."
               exit
            end if
        end do ! davidson iterations

        if (davidson_solve_rhs.and.davidson_overlap_sorting) deallocate(c_small_r_old)
        if (davidson_solve_lhs.and.davidson_overlap_sorting) deallocate(c_small_l_old)

        write(davidson_output_file_unit, *) "Final size of reduced subspace : ",L

        call form_eigenvalues_full_problem(results%eValues, nroots, e_small) 

        if (davidson_solve_lhs) then 
            call form_eigenvector_full_problem(results%eVectorsL, nroots,b,c_small_l,davidson_solve_rhs,davidson_solve_lhs)
            if (davidson_debug) call davidson_print_vectors(results%eVectorsL, label="final eigenvalues full problem")
        end if
        if (davidson_solve_rhs) then 
            call form_eigenvector_full_problem(results%eVectorsR, nroots,b,c_small_r,davidson_solve_rhs,davidson_solve_lhs)
            if (davidson_debug) call davidson_print_vectors(results%eVectorsR, label="final eigenvalues full problem")
        end if

        if (associated(b)) deallocate(b)
        deallocate(M_small)
        if (associated(sigma_vector)) deallocate(sigma_vector)
        deallocate(M_diag)
        deallocate(e_small)
        if (associated(c_small_r)) deallocate(c_small_r)
        if (associated(c_small_l)) deallocate(c_small_l)

    end subroutine


    subroutine setup_overlap_sorting(c_r_old, c_l_old, nroots, solve_rhs, solve_lhs)
        real(kind=8), intent(inout), pointer :: c_r_old(:,:) 
        real(kind=8), intent(inout), pointer :: c_l_old(:,:)  
        integer, intent(in) :: nroots
        logical, intent(in) :: solve_rhs, solve_lhs

        integer :: i, factor 

        if (davidson_complex_mode_with_reals) then
            factor = 2 
        else
            factor = 1 
        end if

        if (solve_rhs) then
            allocate(c_r_old(factor*nroots,nroots))
            c_r_old = 0.0d0
            do i = 1, nroots*factor, factor
                c_r_old(i,(i+1)/2) = 1.0d0
            enddo
        endif
        if (solve_lhs) then
            allocate(c_l_old(factor*nroots,nroots))
            c_l_old = 0.0d0
            do i = 1, nroots*factor, factor
                c_l_old(i,(i+1)/2) = 1.0d0
            enddo
        endif
    end subroutine

    subroutine calculate_eigenvalues_and_eigenvectors(output_file_unit, M, eValues, eVectorsR, eVectorsL, &
                                eVectorsR_old, eVectorsL_old, eValuesShift, left, right,    &
                                sorting_order, nroots, verbose, follow_overlap, apply_evalues_shift)
        real(kind=8), intent(inout), pointer   :: M(:,:)
        real(kind=8), intent(inout), pointer   :: eValues(:,:)
        real(kind=8), intent(inout), pointer   :: eVectorsR(:,:), eVectorsL(:,:)
        real(kind=8), intent(in)   :: eVectorsR_old(:,:), eVectorsL_old(:,:) !old vectors
        real(kind=8), intent(in)   :: eValuesShift(2)
        logical, optional, intent(in) :: left, right, verbose, follow_overlap, apply_evalues_shift
        integer, optional, intent(in) :: sorting_order, nroots
        integer, intent(in) :: output_file_unit 

        logical :: left_local, right_local, verbose_local, follow_overlap_local, apply_evalues_shift_local
        integer :: sorting_order_local, M_nroots, nroot, nroots_to_print
        real(kind=8), allocatable :: overlap(:,:)

        if (present(apply_evalues_shift)) then
            apply_evalues_shift_local = apply_evalues_shift
        else
            apply_evalues_shift_local = .false.
        end if

        if (present(follow_overlap)) then
            follow_overlap_local = follow_overlap
        else
            follow_overlap_local = .false.
        end if

        if (present(sorting_order)) then
            sorting_order_local = sorting_order
        else
            sorting_order_local = 1 
        end if

        if (present(right)) then
            right_local = right
        else
            right_local = .false.
        end if

        if (present(left)) then
            left_local = left
        else
            left_local = .false. 
        end if

        if (present(verbose)) then
            verbose_local = verbose
        else
            verbose_local = .false. 
        end if

        M_nroots = size(M,1) 
        if (davidson_complex_mode_with_reals) then
            M_nroots = M_nroots / 2
        end if

        if (present(nroots)) then
            nroots_to_print = nroots
        else
            nroots_to_print = M_nroots
        end if

        call solve_values_vectors(M, eValues, eVectorsR, eVectorsL, left_local)
        if (davidson_debug) then
            call davidson_print_vectors(eValues, label="debug, eValues after solving  ")
            if (davidson_verbose) then 
                call print_values_vectors_R(davidson_output_file_unit,eValues, eVectorsR, eVectorsL, &
                                            left_local, right_local, nroots=nroots_to_print)
            else
                call print_values_vectors_R(davidson_output_file_unit,eValues, eVectorsR, eVectorsL, &
                                            .false., .false., nroots=nroots_to_print)
            end if
        end if

        if (follow_overlap_local) then
            if (davidson_debug) print *, "debug, following overlap !"
            
            if(right_local) then
                call sort_by_overlap_of_subspace(eValues, eVectorsR, nroots)
            else if (left_local) then
                call sort_by_overlap_of_subspace(eValues, eVectorsL, nroots)
            end if
        else
            if (davidson_debug) print *, "not following overlap !"
            call sort_values_vectors(sorting_order_local, eValues, eVectorsR, eVectorsL, &
                                     right_local, left_local, apply_evalues_shift, eValuesShift, &
                                     debug=davidson_debug)
        end if

        if (davidson_verbose) then 
            call print_values_vectors_R(davidson_output_file_unit,eValues, eVectorsR, eVectorsL, &
                                        left_local, right_local, nroots=nroots_to_print)
        else
            call print_values_vectors_R(davidson_output_file_unit,eValues, eVectorsR, eVectorsL, &
                                        .false., .false., nroots=nroots_to_print)
        end if

    end subroutine

    subroutine form_eigenvalues_full_problem(eValues, nroots, e_small)
        real(kind=8), intent(inout) :: eValues(:,:)
        real(kind=8), intent(in)    :: e_small(:,:)
        integer, intent(in) :: nroots

        integer :: i

        do i = 1, nroots
            eValues(:,i) = e_small(:,i)
        end do
    end subroutine

    subroutine form_eigenvector_full_problem(eVectorsFull,nroots,b,eVectorsSmall,right_vector,left_vector) 
        real(kind=8), intent(inout) :: eVectorsFull(:,:)
        real(kind=8), intent(in)    :: eVectorsSmall(:,:)
        real(kind=8), intent(in)    :: b(:,:)
        logical, intent(in) :: left_vector, right_vector
        integer, intent(in) :: nroots

        integer :: i
        real(kind=8), allocatable :: bc(:,:)
        integer :: vectors_L, vectors_N

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (davidson_debug) then
            print *,'debug sizes 1,2 of b',size(b,1),size(b,2)
            print *,'debug sizes 1,2 of eVectorsFull',size(eVectorsFull,1),size(eVectorsFull,2)
            print *,'debug sizes 1,2 of eVectorsSmall',size(eVectorsSmall,1),size(eVectorsSmall,2)
        end if

        if (left_vector) then
            allocate(bc(vectors_L,vectors_N))
            bc = 0.0d0
!           call davidson_print_vectors(eVectorsSmall, label="L eigenvectors small problem  ")
            if (davidson_complex_mode_with_reals) then
                vectors_N = vectors_N / 2
! aspg: verify here whether transpose or conjugate transpose
                print *, 'warning, calling zgemm with T instead of C in form_eigenvector_full_problem'
                call zgemm('T', 'T', vectors_L, vectors_N, vectors_L, (1.0d0, 0.0d0), eVectorsSmall, &
                           vectors_L, b, vectors_N, (1.0d0, 0.0d0), bc, vectors_L)
            else
                call dgemm('T', 'T', vectors_L, vectors_N, vectors_L, (1.0d0), eVectorsSmall, &
                           vectors_L, b, vectors_N, (1.0d0), bc, vectors_L)
            end if 
!           call davidson_print_vectors(bc, label="bc                            ")
            do i = 1, nroots
! aspg, double check where the index i should go?
!               eVectorsFull(:,i) = bc(:,i) 
                eVectorsFull(:,i) = bc(i,:) 
            end do
        end if

        if (right_vector) then 
            allocate(bc(vectors_N,vectors_L))
            bc = 0.0d0
!           call davidson_print_vectors(eVectorsSmall, label="R eigenvectors small problem  ")
            if (davidson_complex_mode_with_reals) then
                vectors_N = vectors_N / 2
                call zgemm('N', 'N', vectors_N, vectors_L, vectors_L, (1.0d0, 0.0d0), b, &
                           vectors_N, eVectorsSmall, vectors_L, (1.0d0, 0.0d0), bc, vectors_N)
            else
                call dgemm('N', 'N', vectors_N, vectors_L, vectors_L, (1.0d0), b, &
                           vectors_N, eVectorsSmall, vectors_L, (1.0d0), bc, vectors_N)
            end if

!           call davidson_print_vectors(bc, label="bc                            ")
            do i = 1, nroots
                eVectorsFull(:,i) = bc(:,i) 
            end do
        end if

        deallocate (bc)
    end subroutine


    function update_sigma_vector(M, b, old_sigma_vector)
        real(kind=8), pointer :: update_sigma_vector(:,:)
        real(kind=8), intent(in)  :: M(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        real(kind=8), pointer :: old_sigma_vector(:,:)

        integer :: vectors_L, vectors_N
        integer :: vectors_L_old
        integer :: new_L, i

        vectors_N = size(b,1)
        vectors_L = size(b,2)
        vectors_L_old = size(old_sigma_vector,2)

        allocate(update_sigma_vector(vectors_N,vectors_L))
        do i = 1, vectors_L_old 
           update_sigma_vector(:, i) = old_sigma_vector(:, i)  
        end do
        deallocate(old_sigma_vector)

        new_L = vectors_L - vectors_L_old
        if (davidson_complex_mode_with_reals) then
            vectors_N = vectors_N / 2
            call zgemm('N', 'N', vectors_N, new_L, vectors_N, (1.0d0, 0.0d0), M, vectors_N, &
       &               b(:, vectors_L_old:vectors_L), vectors_N, (0.0d0, 0.0d0),            &
       &               update_sigma_vector(:,vectors_L_old:vectors_L), vectors_N)
         else
            call dgemm('N', 'N', vectors_N, new_L, vectors_N, 1.0d0, M, vectors_N, &
       &               b(:, vectors_L_old:vectors_L), vectors_N, 0.0d0,            &
       &               update_sigma_vector(:,vectors_L_old:vectors_L), vectors_N)
         end if
    end function


    function create_sigma_vector_full_matrix(M, b)
        real(kind=8), pointer :: create_sigma_vector_full_matrix(:,:)
        real(kind=8), intent(in)  :: M(:,:)
        real(kind=8), intent(in)  :: b(:,:)
        integer :: vectors_L, vectors_N

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        allocate(create_sigma_vector_full_matrix(vectors_N,vectors_L))

        if (davidson_complex_mode_with_reals) then
            vectors_N = vectors_N / 2
            call zgemm('N', 'N', vectors_N, vectors_L, vectors_N, (1.0d0, 0.0d0), M, vectors_N, &
                       b, vectors_N, (0.0d0, 0.0d0), create_sigma_vector, vectors_N)
        else
            call dgemm('N', 'N', vectors_N, vectors_L, vectors_N, 1.0d0, M, vectors_N, &
                       b, vectors_N, 0.0d0, create_sigma_vector, vectors_N)
        end if
    end function

    subroutine sort_values_vectors(order, e_small, c_small_r, c_small_l, right, left, apply_shift, energy_shift, debug) 
        integer, intent(in) :: order ! 1: ascending, -1: descending
        logical, intent(in) :: right ! if true, sort also the right eigenvectors
        logical, intent(in) :: left  ! if true, sort also the left eigenvectors
        logical, intent(in) :: apply_shift ! if true, sort values of |e_small - energy_shift| instead of e_small
        logical, intent(in), optional :: debug
        real(kind=8), intent(inout) :: e_small(:,:)
        real(kind=8), intent(inout) :: c_small_r(:,:)
        real(kind=8), intent(inout) :: c_small_l(:,:)
        real(kind=8), intent(in) :: energy_shift(2)

        real(kind=8), allocatable :: temp_vector(:,:), e_ascending(:,:)
        real(kind=8) :: temp
        integer, allocatable :: indices_sorted_array(:), i_ascending(:)
        integer :: vectors_N, vectors_L, i_sorted, i, j
        logical :: debug_local

        vectors_L = size(e_small,2)
        allocate(indices_sorted_array(vectors_L))

        if (present(debug)) debug_local = debug

        if (abs(order).ne.1) then
           print *,"In sort_values_vectors: variable order given a value of ",order
           print *,"  allowed values:   1 for sort to ascending order"
           print *,"                   -1 for sort to descending order"
           stop "Unknown ordering requested in sort_values_vectors"
        end if

!       all the sort routines sort in ascending order, which is almost always what we want. 
!       when we want descending order, we will proceed with sorting in ascending order, and then
!       just inverse the final arrays, for the energies and indices

        do i = 1, vectors_L 
            indices_sorted_array(i) = i 
        end do

        if (debug_local) print *, "debug sort_values_vectors, L ",vectors_L

        if (debug_local) call davidson_print_vectors(e_small, label="debug, e_small before sort    ")

        if (davidson_complex_mode_with_reals) then
            vectors_N = 2*vectors_L
            if (apply_shift) then
               call quicksort_vector_complex_as_real_shifted(energy_shift, e_small, indices_sorted_array, 1, vectors_L)
            else
               call quicksort_vector_complex_as_real(e_small, indices_sorted_array, 1, vectors_L)
            end if

        else
            vectors_N = vectors_L
            if (apply_shift) then
               call quicksort_vector_real_shifted(energy_shift(1), e_small(1,:), indices_sorted_array, 1, vectors_L)
            else
               call quicksort_vector_real(e_small(1,:), indices_sorted_array, 1, vectors_L)
            end if

       ! we sort the indices that would sort the eigenvalues array, and apply the change to all other arrays
            do i = 1, vectors_L
                i_sorted = indices_sorted_array(i)

                temp = e_small(2,i)
                e_small(2,i) = e_small(2,i_sorted)
                e_small(2,i_sorted) = temp
            end do
        end if
        if (debug_local) call davidson_print_vectors(e_small, label="debug, e_small  after sort    ")
        if (debug_local) print *, "indices, sorted",indices_sorted_array

        if (order.eq.-1) then

           if (debug_local) then
              print *, "ascending order"
              print *, " indices ",indices_sorted_array
              print *, " evalues ",e_small
           end if

           allocate(i_ascending(vectors_L))
           allocate(e_ascending(2,vectors_L))

           i_ascending = indices_sorted_array
           e_ascending = e_small

           indices_sorted_array = i_ascending(vectors_L:1:-1)

           do i = 1, vectors_L
              j = vectors_L - i + 1 
              e_small(1,j) = e_ascending(1,i)
              e_small(2,j) = e_ascending(2,i)
           end do

           deallocate(i_ascending)
           deallocate(e_ascending)

           if (debug_local) then
              print *, "descending order"
              print *, " indices ",indices_sorted_array
              print *, " evalues ",e_small
          end if

        end if

        allocate(temp_vector(vectors_N,vectors_L))
        temp_vector = 0.0d0

        if (right) then
            if (debug_local) call davidson_print_vectors(c_small_r,   label="debug, c_small_r before sort  ")
            do i = 1, vectors_L
             i_sorted = indices_sorted_array(i)
             temp_vector(:,i) = c_small_r(:,i_sorted)
            end do
            if (debug_local) call davidson_print_vectors(temp_vector, label="debug, temp_vector            ")
            c_small_r = temp_vector
            if (debug_local) call davidson_print_vectors(c_small_r, label="debug, c_small_r  after sort  ")
        end if

        if (left) then
            if (debug_local) call davidson_print_vectors(c_small_l, label="debug, c_small_l before sort  ")
            do i = 1, vectors_L
             i_sorted = indices_sorted_array(i)
             temp_vector(:,i) = c_small_l(:,i_sorted)
            end do
            c_small_l = temp_vector
            if (debug_local) call davidson_print_vectors(c_small_l, label="debug, c_small_l  after sort  ")
        end if
        deallocate(temp_vector)
        deallocate(indices_sorted_array)

    end subroutine

    recursive subroutine quicksort_vector_real(v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:)
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, t
        integer :: i, j, ti

        x = v( (first + last) / 2 )
        i = first
        j = last
        do
           do while (v(i) < x)
              i = i + 1
           end do
           do while (x < v(j))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(i)
           v(i) = v(j)
           v(j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) call quicksort_vector_real(v, indices_sorted, first, i-1)
        if (j+1 < last)  call quicksort_vector_real(v, indices_sorted, j+1, last)
    end subroutine

    real(kind=8) function shifted_value_real(unshifted_value, shift)
       real(kind=8), intent(in)  ::  shift, unshifted_value

       shifted_value_real = abs(unshifted_value - shift)
    end function shifted_value_real

    real(kind=8) function shifted_value_complex_as_real(unshifted_value, shift)
       real(kind=8), intent(in)  ::  shift(2), unshifted_value(2)
       real(kind=8)  :: temp(2)

       temp(1) = unshifted_value(1) - shift(1)
       temp(2) = unshifted_value(2) - shift(2)

       shifted_value_complex_as_real = sqrt(temp(1)*temp(1) + temp(2)*temp(2))
    end function shifted_value_complex_as_real

    recursive subroutine quicksort_vector_real_shifted(shift, v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:)
        real(kind=8), intent(in) ::  shift
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, y, t
        integer :: i, j, ti

        x = shifted_value_real(v( (first + last) / 2 ),shift)
        i = first
        j = last
        do
           do while (shifted_value_real(v(i),shift) < x)
              i = i + 1
           end do
           do while (x < shifted_value_real(v(j),shift))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(i)
           v(i) = v(j)
           v(j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) call quicksort_vector_real_shifted(shift, v, indices_sorted, first, i-1)
        if (j+1 < last)  call quicksort_vector_real_shifted(shift, v, indices_sorted, j+1, last)
    end subroutine 

    recursive subroutine quicksort_vector_complex_as_real(v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:,:)
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, t, value
        integer :: i, j, ti, start_index

        start_index = (first + last) / 2 

        x = sqrt(v(1,start_index)*v(1,start_index) + v(2,start_index)*v(2,start_index))

        i = first
        j = last

        do
           do while ((sqrt(v(1,i)*v(1,i) + v(2,i)*v(2,i)) < x))
              i = i + 1
           end do
           do while ((x < sqrt(v(1,j)*v(1,j) + v(2,j)*v(2,j))))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(1,i)
           v(1,i) = v(1,j)
           v(1,j) = t

           t = v(2,i)
           v(2,i) = v(2,j)
           v(2,j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) then 
            call quicksort_vector_complex_as_real(v, indices_sorted, first, i-1)
        end if
        if (j+1 < last)  then
            call quicksort_vector_complex_as_real(v, indices_sorted, j+1, last)
        end if
    end subroutine

    recursive subroutine quicksort_vector_complex_as_real_shifted(shift, v, indices_sorted, first, last)
        real(kind=8), intent(inout) ::  v(:,:)
        real(kind=8), intent(in) :: shift(2) 
        integer, intent(inout) :: indices_sorted(:)
        integer, intent(in)  :: first, last

        real(kind=8) :: x, t, value
        integer :: i, j, ti, start_index

        start_index = (first + last) / 2

        x = shifted_value_complex_as_real(v(:,start_index), shift)

        i = first
        j = last

        do
           do while ( shifted_value_complex_as_real(v(:,i), shift) < x)
              i = i + 1
           end do
           do while ( x < shifted_value_complex_as_real(v(:,j), shift))
              j = j - 1
           end do
           if (i >= j) exit

           t = v(1,i)
           v(1,i) = v(1,j)
           v(1,j) = t

           t = v(2,i)
           v(2,i) = v(2,j)
           v(2,j) = t

           ti = indices_sorted(i)
           indices_sorted(i) = indices_sorted(j)
           indices_sorted(j) = ti

           i = i + 1
           j = j - 1
        end do
        if (first < i-1) then
            call quicksort_vector_complex_as_real_shifted(shift, v, indices_sorted, first, i-1)
        end if
        if (j+1 < last)  then
            call quicksort_vector_complex_as_real_shifted(shift, v, indices_sorted, j+1, last)
        end if
    end subroutine



! I/O for davidson

    subroutine io_davidson_setup(N, L, nroot, complex_mode, left, right, symmetry, node)
        integer, intent(in) :: N, L, nroot, node, symmetry
        logical, intent(in) :: left, right, complex_mode
        integer :: reclen_real, reclen_vector, factor
        real(kind=8) :: dummy

        factor = 1
        if (complex_mode) factor = 2

        inquire (iolength = reclen_real) dummy

        left_sigma_file_name     = 'LSIGMAV.sym'//trim(str(symmetry))//"."//trim(str(node))
        left_vectors_file_name   = 'LEVECTO.sym'//trim(str(symmetry))//"."//trim(str(node))
        left_subspace_file_name  = 'LSUBSPA.sym'//trim(str(symmetry))//"."//trim(str(node))

        right_sigma_file_name    = 'RSIGMAV.sym'//trim(str(symmetry))//"."//trim(str(node))
        right_vectors_file_name  = 'REVECTO.sym'//trim(str(symmetry))//"."//trim(str(node))
        right_subspace_file_name = 'RSUBSPA.sym'//trim(str(symmetry))//"."//trim(str(node))

        if (right) then
            reclen_vector = factor*N*L*reclen_real
            open(unit=right_sigma_file_unit,file=right_sigma_file_name,&
                 form="unformatted", access="direct", recl=reclen_vector)
            open(unit=right_subspace_file_unit, file=right_subspace_file_name, &
                 form="unformatted", access="direct", recl=reclen_vector)

            reclen_vector = factor*N*nroot*reclen_real + 2 
            open(unit=right_vectors_file_unit, file=right_vectors_file_name, &
                 form="unformatted", access="direct", recl=reclen_vector)
        end if

        if (left) then
            reclen_vector = factor*N*L*reclen_real
            open(unit=left_sigma_file_unit, file=left_sigma_file_name, &
                 form='unformatted', access='direct', recl=reclen_vector)
            open(unit=left_subspace_file_unit, file=left_subspace_file_name, &
                 form='unformatted', access='direct', recl=reclen_vector)

            reclen_vector = factor*N*nroot*reclen_real + 2 
            open(unit=left_vectors_file_unit, file=left_vectors_file_name, &
                 form='unformatted', access='direct', recl=reclen_vector)
        end if

    end subroutine

    subroutine io_davidson_finalize(left,right)
        logical, intent(in) :: left, right

        if (right) then
            close(unit=right_sigma_file_unit)
            close(unit=right_vectors_file_unit)
            close(unit=right_subspace_file_unit)
        end if

        if (left) then
            close(unit=left_sigma_file_unit)
            close(unit=left_vectors_file_unit)
            close(unit=left_subspace_file_unit)
        end if

    end subroutine

    subroutine write_davidson_values_vectors(file_unit,valuesRe,valuesIm,vectors,column_start,column_end)
        real(kind=8), intent(in) :: vectors(:,:), valuesRe(:), valuesIm(:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            write(file_unit, rec=i) valuesRe(i), valuesIm(i), vectors(:,i)
        end do

    end subroutine

    subroutine read_davidson_values_vectors(file_unit,valuesRe,valuesIm,vectors,column_start,column_end)
        real(kind=8), intent(inout) :: vectors(:,:), valuesRe(:), valuesIm(:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            read(file_unit, rec=i) valuesRe(i), valuesIm(i), vectors(:,i)
        end do

    end subroutine

    subroutine write_davidson_data_1d(file_unit,data1d)
        real(kind=8), intent(inout) :: data1d(:)
        integer, intent(in) :: file_unit

        write(file_unit, rec=1) data1d(:)
    end subroutine

    subroutine read_davidson_data_1d(file_unit,data1d)
        real(kind=8), intent(inout) :: data1d(:)
        integer, intent(in) :: file_unit

        read(file_unit, rec=1) data1d(:)
    end subroutine

    subroutine write_davidson_data_2d(file_unit,data2d,column_start,column_end)
        real(kind=8), intent(in) :: data2d(:,:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            write(file_unit, rec=i) data2d(:,i)
        end do
    end subroutine

    subroutine read_davidson_data_2d(file_unit,data2d,column_start,column_end)
        real(kind=8), intent(inout) :: data2d(:,:)
        integer, intent(in) :: column_start, column_end, file_unit
        integer :: i

        do i = column_start, column_end
            read(file_unit, rec=i) data2d(:,i)
        end do
    end subroutine

! taken from http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran
    character(len=20) function str(k)
!       "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

    subroutine dump_to_file(fileno, filename, matrix)
        integer, intent(in) :: fileno
        character(len=*), intent(in) :: filename
        real(kind=8) :: matrix(:,:)
        integer :: irec
 
        inquire (iolength = irec) matrix
        open (unit=fileno, file=filename, form="unformatted", access="direct", recl=irec)
        write(unit=fileno, rec=1) matrix
        close(unit=fileno)
    end subroutine

    subroutine sort_by_overlap_of_subspace(e_small, c_small, nroots)
!
! this routine sorts the eigenvectors with respect to the norm of their projection on
! to the set of trial vectors. Sorting is from high overlap to low overlap,
! so that the first nroots vectors can be taken to extend the space.
!
! Lucas Visscher (2019), based on a routine by ASP Gomes.
!
        integer, intent(in)         :: nroots
        real(kind=8), intent(inout) :: c_small(:,:) ! left or right-side reduced subspace eigenvectors (can be real or complex)
        real(kind=8), intent(inout) :: e_small(:,:) ! reduced subspace eigenvalues (complex, as we deal with a non-hermitian eigenvalue problem)

        real(kind=8), allocatable :: overlap_with_subspace(:)
        real(kind=8), allocatable :: sorted_eigenvalues(:,:)
        real(kind=8), allocatable :: sorted_eigenvectors(:,:)

        real(kind=8) :: overlap
        integer :: subspace_dimension, iroot, itrial, i
        integer :: rcw ! 1 for real, 2 for complex

        subspace_dimension = size(c_small,2)
        rcw = size(c_small,1)/subspace_dimension  ! real/imaginary degree of freedom is "hidden" in the first dimension of this square matrix

        allocate(overlap_with_subspace(subspace_dimension))
        overlap_with_subspace = 0.d0
        do iroot = 1, subspace_dimension
           do itrial = 1, nroots
              if (rcw==1) overlap = abs(c_small(itrial,iroot))
              if (rcw==2) overlap = abs(cmplx(c_small(2*itrial-1,iroot),c_small(2*itrial,iroot)))
              overlap_with_subspace(iroot) = overlap_with_subspace(iroot) + overlap**2
           end do
        end do

        allocate (sorted_eigenvalues(2,subspace_dimension))
        allocate (sorted_eigenvectors(rcw*subspace_dimension,subspace_dimension))
        do iroot = 1, subspace_dimension
           i = maxloc(overlap_with_subspace,1)
           sorted_eigenvalues(:,iroot)  = e_small(:,i)
           sorted_eigenvectors(:,iroot) = c_small(:,i)
           overlap_with_subspace(i) = -max(1.d-12,overlap_with_subspace(i))
           ! with this choice we select the next highest value for the next root without loosing the overlap information, the
           ! max construct is there in case we have multiple vectors with zero overlap.
        end do

        ! copy the sorted vectors back to the original arrays
        e_small = sorted_eigenvalues
        c_small = sorted_eigenvectors

        deallocate (sorted_eigenvalues)
        deallocate (sorted_eigenvectors)
        deallocate (overlap_with_subspace)

    end subroutine sort_by_overlap_of_subspace

    subroutine apply_projector(projector,b)
        logical, intent(in) :: projector(:)
        real(kind=8), intent(inout)  :: b(:,:)
        integer :: vectors_L, vectors_N, i, j, k 

! the values of projector mean 
!  .true. : the entry will remain as is
! .false. : the entre will be set to zero

        vectors_N = size(b,1)
        vectors_L = size(b,2)

        if (davidson_complex_mode_with_reals) vectors_N = vectors_N / 2

        do i = 1, vectors_N
           if (projector(i).eqv..false.) then
              if (davidson_complex_mode_with_reals) then
                 j = 2*(i-1) + 1
                 k = j + 1
                 b(j,:) = 0.0d0
                 b(k,:) = 0.0d0
              else 
                 b(i,:) = 0.0d0
              end if
           end if
        end do
    end subroutine

    subroutine davidson_get_eigenvalue_shift(apply_shift,shift)
       logical, intent(inout)      :: apply_shift
       real(kind=8), intent(inout) :: shift(2)

       shift = davidson_eigenvalues_shift
       apply_shift = davidson_apply_eigenvalues_shift 

    end subroutine


   subroutine davidson_synchronize_trial_vectors(b,nr_to_sync)
! this routine synchronizes the matrix of trial vectors (b) from the 
! davidson procedure across nodes, with the process with rank 0, through 
! a brodcast. 
!
! the routine takes as second argument the number of trial vectors to 
! synchronize (nr_to_sync). 
!
! the routine can oprate in two ways: 
! - if nr_to_sync is smaller or equal to zero, this signals all vectors are broadcast.
! - if nr_to_sync is non-zero, nr_to_sync vectors are broadcast, but starting
!   with the most recent ones which are stored in the last columns of b. 
!
#ifdef VAR_MPI
      use interface_to_mpi
#endif
      real(kind=8), intent(inout) :: b(:,:)
      integer, intent(in) :: nr_to_sync 
      integer :: ndim1, ndim2, start_slice, end_slice, i

#ifdef VAR_MPI
      ndim1 = size(b,1)
      if (nr_to_sync .le. 0) then 
          ndim2 = size(b,2) 
          do i = 1, ndim2
             call interface_mpi_bcast(b(1,i),ndim1,0,global_communicator)
          end do
      else
         end_slice = size(b,2)
         start_slice = end_slice - nr_to_sync
          do i = start_slice + 1, end_slice
             call interface_mpi_bcast(b(1,i),ndim1,0,global_communicator)
          end do
      end if
#else
      return
#endif

   end subroutine

end module

