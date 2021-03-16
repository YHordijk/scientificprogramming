program test_davidson

   use davidson

   implicit none

    type intermediates
      real(8),allocatable,dimension(:,:) :: M
    endtype intermediates

! eigenvectors, eigenvalues 
   double precision, pointer :: eval_real(:)
   double precision, pointer :: evec_real(:,:)
   double precision, pointer :: evec_real_l(:,:)
   double precision, pointer :: evec_real_r(:,:)

   double precision, pointer :: eval_r(:)
   double precision, pointer :: eval_i(:)
   complex (kind=8), pointer :: eval_complex(:)
   complex (kind=8), pointer :: evec_complex(:,:)
   complex (kind=8), pointer :: evec_complex_r(:,:)
   complex (kind=8), pointer :: evec_complex_l(:,:)

! matrices to diagolize
   double precision, pointer :: m_real_symmetric_packed(:)
   double precision, pointer :: m_real_symmetric_full(:,:)
   double precision, pointer :: m_real_nonsymmetric_full(:,:)

   complex (kind=8), pointer :: m_complex_hermitian_packed(:)
   complex (kind=8), pointer :: m_complex_hermitian_full(:,:)
   complex (kind=8), pointer :: m_complex_nonsymmetric_full(:,:)
!
   logical :: diagonalize_real_symm_full       = .true.
   logical :: diagonalize_real_nonsymm_full    = .false.
!
   integer :: mdim
   integer :: nroots 
   integer :: sorting_order = 1

   mdim =  128
   nroots = 4

! allocating eigenvectors and eigenvalues storage
   allocate(eval_real(mdim))
   allocate(eval_r(mdim))
   allocate(eval_i(mdim))
   allocate(eval_complex(mdim))
   eval_real = 0.0d0 
   eval_r = 0.0d0 
   eval_i = 0.0d0 
   eval_complex = cmplx(0.0d0, 0.0d0)

   allocate(evec_real(mdim,mdim))
   allocate(evec_real_l(mdim,mdim))
   allocate(evec_real_r(mdim,mdim))

   evec_real = 0.0d0 
   evec_real_l = 0.0d0 
   evec_real_r = 0.0d0 

   allocate(evec_complex(mdim,mdim))
   allocate(evec_complex_l(mdim,mdim))
   allocate(evec_complex_r(mdim,mdim))
   evec_complex = cmplx(0.0d0, 0.0d0)
   evec_complex_l = cmplx(0.0d0, 0.0d0)
   evec_complex_r = cmplx(0.0d0, 0.0d0)

! reference diagonalization for real symmetric full storage matrix

   allocate(m_real_symmetric_full(mdim,mdim))
   call setup_real_symmetric_matrix_full(m_real_symmetric_full, mdim)
   print *, "test, full real symmetric matrix"
!  call print_real_matrix_full(m_real_symmetric_full, mdim, mdim)
   if (diagonalize_real_symm_full) then
      evec_real = m_real_symmetric_full
      print *, "full diagonalization with dgeev"
      call diagonalize_real_general_matrix_full(m_real_symmetric_full,eval_r,eval_i,evec_real_r,evec_real_l)
!     call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, mdim)
      print *, "sorting eigenvectors and eigenvalues in ascending order"
      call sort_values_vectors(sorting_order,eval_r, eval_i, evec_real_r, evec_real_l, .true.)
!     call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, mdim)
      call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, nroots)

      print *, "davidson diagonalization for lowest root(s)" 
      call setup_real_symmetric_matrix_full(m_real_symmetric_full, mdim)
      eval_r      = 0.0d0
      eval_i      = 0.0d0
      evec_real_l = 0.0d0
      evec_real_r = 0.0d0
      call diagonalize_via_davidson(m_real_symmetric_full, nroots, .true., eval_r, eval_i, evec_real_l, evec_real_r)
      call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, nroots)

   endif

!i reference diagonalization for real nonsymmetric full storage matrix

   allocate(m_real_nonsymmetric_full(mdim,mdim))
   call setup_real_nonsymmetric_matrix_full(m_real_nonsymmetric_full, mdim)
   print *, "test, full real non-symmetric matrix"
!  call print_real_matrix_full(m_real_nonsymmetric_full, mdim, mdim)
   if (diagonalize_real_nonsymm_full) then
      print *, "full diagonalization with dgeev"
      call diagonalize_real_general_matrix_full(m_real_nonsymmetric_full,eval_r,eval_i,evec_real_r,evec_real_l)
!     call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, mdim)
      print *, "sorting eigenvectors and eigenvalues in ascending order"
      call sort_values_vectors(sorting_order,eval_r, eval_i, evec_real_r, evec_real_l, .true.)
      call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, nroots)

      print *, "davidson diagonalization for lowest root(s)"
      call setup_real_nonsymmetric_matrix_full(m_real_nonsymmetric_full, mdim)
      eval_r      = 0.0d0
      eval_i      = 0.0d0
      evec_real_l = 0.0d0
      evec_real_r = 0.0d0
      call diagonalize_via_davidson(m_real_nonsymmetric_full, nroots, .false., eval_r, eval_i, evec_real_l, evec_real_r)
      call print_real_eigenvalues_eigenvectors(eval_r, eval_i, evec_real_l, evec_real_r, mdim, nroots)
   end if

contains


   subroutine print_real_eigenvalues_complex_eigenvectors(eval,evec_l,evec_r,mdim, nroots)
      real(kind=8) :: eval(:)
      complex(kind=8) :: evec_l(:,:), evec_r(:,:)
      integer :: mdim
      integer :: nroots

      print *, ""
      print *, "                ***       eigenvalues  ***"
      call print_real_eigenvalues(eval, nroots)
      print *, "                *** left  eigenvectors ***"
      call print_complex_matrix_full(evec_l, mdim, nroots)
      print *, "                *** right eigenvectors ***"
      call print_complex_matrix_full(evec_r, mdim, nroots)
      print *, ""
      print *, ""
      print *, ""

   end subroutine

   subroutine print_real_eigenvalues_eigenvectors(eval_r,eval_i,evec_l,evec_r,mdim, nroots)
      real(kind=8), pointer :: eval_r(:), eval_i(:)
      real(kind=8), pointer :: evec_l(:,:), evec_r(:,:)
      integer :: mdim
      integer :: nroots

      print *, ""
      print *, "                ***       eigenvalues (real part) ***"
      call print_real_eigenvalues(eval_r, nroots)
      if (associated(eval_i)) then
          print *, "                ***       eigenvalues (imaginary part) ***"
          call print_real_eigenvalues(eval_i, nroots)
      end if
      print *, "                *** left  eigenvectors ***"
      call print_real_matrix_full(evec_l, mdim, nroots)
      print *, "                *** right eigenvectors ***"
      call print_real_matrix_full(evec_r, mdim, nroots)
      print *, ""
      print *, ""
      print *, ""

   end subroutine

   subroutine print_complex_eigenvalues_eigenvectors(eval,evec_l,evec_r,mdim, nroots)
      complex(kind=8) :: eval(:)
      complex(kind=8) :: evec_l(:,:), evec_r(:,:)
      integer :: mdim
      integer :: nroots 

      print *, ""
      print *, "                ***       eigenvalues  ***"
      call print_complex_eigenvalues(eval, nroots)
      print *, "                *** left  eigenvectors ***"
      call print_complex_matrix_full(evec_l, mdim, nroots)
      print *, "                *** right eigenvectors ***"
      call print_complex_matrix_full(evec_r, mdim, nroots)
      print *, ""
      print *, ""
      print *, ""

   end subroutine

   complex function delta_complex(i,j)
      integer :: i, j
      if (i .eq. j) then
         delta_complex = cmplx(1.0d0, 0.0d0)
      else
         delta_complex = cmplx(0.0d0, 0.0d0)
      endif
   end function

   double precision function delta_real(i,j)
      integer :: i, j
      if (i .eq. j) then
         delta_real = 1.0d0
      else
         delta_real = 0.0d0
      endif
   end function

   subroutine setup_real_symmetric_matrix_packed(m, mdim)
      integer :: i, j, k, ij
      integer :: ii, jj
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:)

      k = mdim/2

      ij = 0
      m = 0.0d0
      do i = 1, mdim
         ii = i - 1
         do j = i, mdim
            jj = j - 1

            ij=i+(j*(j-1)/2)
            m(ij) = ii*delta_real(ii,jj) + (ii - jj - k*k)
         enddo
      enddo
   end subroutine

   subroutine setup_complex_hermitian_matrix_packed(m, mdim)
      integer :: i, j, k, ij
      integer :: ii, jj
      integer, intent(in) :: mdim
      complex (kind=8), intent(inout) :: m(:)
      real (kind=8) :: real_part, imag_part

      k = mdim/2

      ij = 0
      m = cmplx(0.0d0, 0.0d0)
      do i = 1, mdim
         ii = i - 1
         do j = i, mdim
            jj = j - 1

            ij=i+(j*(j-1)/2)
            real_part = ii*delta_real(i,j) + (ii - jj - k*k)
!           imag_part = (ii*delta_real(i,j) + (ii - jj))/mdim
            imag_part = 0.0d0
            m(ij) = cmplx(real_part, imag_part)
         enddo
      enddo
   end subroutine

       
   subroutine setup_real_symmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:,:)

      k = mdim/2

      m = 0.0d0
      do i = 1, mdim
         ii = i - 1
         do j = i, mdim
            jj = j - 1

            m(i,j) = ii*delta_real(ii,jj) + (ii - jj - k*k)
            m(j,i) = m(i,j)
         enddo
      enddo
   end subroutine

   subroutine setup_complex_hermitian_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      real (kind=8) :: real_part, imag_part
      complex (kind=8), intent(inout) :: m(:,:)

      k = mdim/2

      m = cmplx(0.0d0, 0.0d0)
      do i = 1, mdim
         ii = i - 1
         do j = i, mdim
            jj = j - 1

            real_part = ii*delta_real(i,j) + (ii - jj - k*k)
!           imag_part = (ii*delta_real(i,j) + (ii - jj))/mdim
            imag_part = 0.0d0
            m(i,j) = cmplx(real_part, imag_part)
            if (j .gt. i) m(j,i) = conjg(m(i,j))
         enddo
      enddo
   end subroutine

  subroutine setup_real_nonsymmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:,:)
    
      k = mdim/2 
            
      m = 0.0d0
      do i = 1, mdim 
         ii = i - 1
         do j = 1, mdim 
            jj = j - 1
            if (jj .le. k) then
               m(i,j) = ii*delta_real(i,j) - (ii - jj - k*k)
            end if
            if (jj .le. mdim) then
               m(i,j) = ii*delta_real(i,j) + (ii - jj - k*k)
            end if
         enddo
      enddo

   end subroutine

  subroutine setup_complex_nonsymmetric_matrix_full(m, mdim)
      integer :: i, j, k
      integer :: ii, jj
      integer, intent(in) :: mdim
      real (kind=8) :: real_part, imag_part
      complex (kind=8), intent(inout) :: m(:,:)
    
      k = mdim/2 
            
      m = cmplx(0.0d0,  0.0d0)
      do i = 1, mdim 
         ii = i - 1
         do j = 1, mdim 
            jj = j - 1

            if (jj .le. k) then
               real_part = ii*delta_real(i,j) - (ii - jj - k*k)
!              imag_part = (ii*delta_real(i,j) - (ii - jj))/mdim
               imag_part = 0.0d0
               m(i,j) = cmplx(real_part, imag_part)
            else if ((jj .le. mdim).and.(jj .ge. (k+1))) then
               real_part = ii*delta_real(i,j) + (ii - jj - k*k)
!              imag_part = (ii*delta_real(i,j) - (ii - jj))/mdim
               imag_part = 0.0d0
               m(i,j) = cmplx(real_part, imag_part)
            endif
         enddo
      enddo
   end subroutine


   subroutine print_real_matrix_packed(m, mdim)
      integer :: i, j, ij
      integer, intent(in) :: mdim
      double precision, intent(inout) :: m(:)
      double precision :: my_zero = 0.0d0

      do i = 1, mdim
         do j = 1, mdim
            if (i .le. j) then
               ij = i + j*(j-1)/2 
               print "(4x,1e14.6)", m(ij)
            else
               print "(4x,1e14.6)", my_zero 
            endif
         enddo
         print *, ""
      enddo
   end subroutine

   subroutine print_complex_matrix_packed(m, mdim)
      integer :: i, j, ij
      integer, intent(in) :: mdim
      complex (kind=8), intent(inout) :: m(:)
      complex (kind=8) :: value
      character(4) :: imag_sign = "   +"

      do i = 1, mdim
         do j = 1, mdim
            if (i .le. j) then
               ij = i + j*(j-1)/2

               if (aimag(m(ij)) .ge. 0.0d0) then
                imag_sign = "   +"
               else
                imag_sign = "   -"
               endif
               print "(4x,1e14.6,a,1e14.6,a)", real(m(ij)),imag_sign,abs(aimag(m(ij)))," i"
            else
               ij = j + i*(i-1)/2
               value = conjg(m(ij))
               if (aimag(value) .ge. 0.0d0) then
                imag_sign = "   +"
               else
                imag_sign = "   -"
               endif
               print "(4x,1e14.6,a,1e14.6,a)", real(value),imag_sign,abs(aimag(value))," i"
            endif
         enddo
         print *, ""
      enddo
   end subroutine


   subroutine print_real_matrix_full(m, na, nb)
      integer :: i, j
      integer, intent(in) :: na, nb 
      double precision, intent(inout) :: m(:,:)

      do i = 1, na 
         do j = 1, nb
            print "(4x,1e14.6,$)", m(i,j)
         enddo
         print *, ""
      enddo
   end subroutine

   subroutine print_real_eigenvalues(e, n)
      integer :: i, n
      double precision, intent(in) :: e(:)

      do i = 1, n
         print "(4x,1e14.6,$)", e(i)
      enddo
      print *, ""
   end subroutine

   subroutine print_complex_eigenvalues(e, n)
      integer :: i, n
      complex (kind=8), intent(in) :: e(:)
      character(4) :: imag_sign = "   +"

      do i = 1, n
         if (aimag(e(i)) .ge. 0.0d0) then
             imag_sign = "   +"
         else
             imag_sign = "   -"
         endif

         print "(4x,1e14.6,a,1e14.6,a,$)", real(e(i)),imag_sign,abs(aimag(e(i)))," i"
      enddo
      print *, ""
   end subroutine

   subroutine print_complex_matrix_full(m, na, nb)
      integer :: i, j
      integer, intent(in) :: na, nb 
      complex (kind=8), intent(inout) :: m(:,:)
      character(4) :: imag_sign = "   +"

      do i = 1, na 
         do j = 1, nb 
            if (aimag(m(i,j)) .ge. 0.0d0) then
                imag_sign = "   +" 
            else
                imag_sign = "   -" 
            endif

            print "(4x,1e14.6,a,1e14.6,a,$)", real(m(i,j)),imag_sign,abs(aimag(m(i,j)))," i"
         enddo
         print *, ""
      enddo

   end subroutine

   subroutine diagonalize_real_general_matrix_full(m,eval_r,eval_i,evec_r,evec_l)
      double precision, intent(inout) :: m(:,:)
      double precision, intent(out), pointer :: eval_r(:), eval_i(:)
      double precision, intent(out), pointer :: evec_r(:,:), evec_l(:,:)

      integer :: vectors_L, i

      integer :: lwork, info
      character(len=1) :: do_left = 'V', do_right = 'V'
      real(kind=8), pointer :: work(:)

      vectors_L = size(m,1)
      allocate(evec_l(vectors_L,vectors_L))
      allocate(evec_r(vectors_L,vectors_L))
      allocate(eval_r(vectors_L))
      allocate(eval_i(vectors_L))

! first query optimal lwork value
        lwork = -1
        allocate(work(4*vectors_L))
        call dgeev (do_left, do_right, vectors_L, m, vectors_L,   &
                    eval_r, eval_i, evec_l, vectors_L, evec_r,    &
                    vectors_L, WORK, LWORK, INFO)
        lwork = work(1)
        deallocate(work)
        allocate(work(lwork))
        call dgeev (do_left, do_right, vectors_L, m, vectors_L,   &
                    eval_r, eval_i, evec_l, vectors_L, evec_r,    &
                    vectors_L, WORK, LWORK, INFO)
        if (info .gt. 0) then
            print *, 'in small_solve_eigenvalue: not all eigenvalues converged'
            print *, 'converged'
            do i = info+1, vectors_L
                print *, " real part : ",eval_r(i), " imaginary part : ",eval_i(i)
            end do
        end if


   end subroutine 

   subroutine diagonalize_via_davidson(m, nroots, symmetric, eval_r, eval_i, evec_l, evec_r)
       use davidson

      double precision, intent(inout) :: m(:,:)
      double precision, intent(inout), pointer :: eval_r(:), eval_i(:)
      double precision, intent(inout), pointer :: evec_r(:,:), evec_l(:,:)
      logical, intent(in) :: symmetric 
      integer, intent(in) :: nroots
      logical :: left = .true., right=.true.

      type(davidson_results) :: davResults
      real(kind=8), pointer :: b_guess(:,:)
      type(intermediates) :: B
     
      logical :: verbose = .false.

       call set_safe_minimum()
       call davidson_driver(nroots, m, symmetric, eval_r, eval_i, evec_l, evec_r, left, right, verbose)

       subroutine davidson_driver(nroots, A, symmetric, results, verbose, b_guess)

   end subroutine



end program 
