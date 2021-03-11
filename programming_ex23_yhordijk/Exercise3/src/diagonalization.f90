module Diagonalization

implicit none

private

public :: diagonalize

interface diagonalize
   module procedure diagonalize
end interface

contains

   subroutine diagonalize (matrix,eigenvectors,eigenvalues)

   real*8, intent(in)               :: matrix(:,:)
   real*8, intent(out)              :: eigenvalues(:)
   real*8, intent(out), optional    :: eigenvectors(:,:)
   real*8, allocatable              :: a(:,:),work(:)
   character, parameter             :: uplo='L'
   character                        :: jobz

   integer                          :: n, lda, info, lwork
   integer                          :: ILAENV

!  Determine the calculation mode
   if (present(eigenvectors)) then
      jobz = 'V'
   else
      jobz = 'N'
   end if

!  Error checking
   if (size(matrix,1) /= size(matrix,2)) then
      print*," diagonalize assumes square matrices"
      stop "error in diagonalize routine"
   end if
   if (size(matrix,1) > size(eigenvalues)) then
      print*," dimension of eigenvalue array too small "
      stop "error in diagonalize routine"
   end if
   if (present(eigenvectors)) then
      if (size(eigenvectors,1) /= size(eigenvectors,2)) then
         print*," diagonalize assumes square matrices"
         stop "error in diagonalize routine"
      end if
      if (size(matrix,1) /= size(eigenvectors,1)) then
         print*," diagonalize assumes same size array for matrix and eigenvectors "
         stop "error in diagonalize routine"
      end if
   end if

!  Initialize and diagonalize using lapack's dsyev routine
   n = size(matrix,1)
   lda = n
   lwork = n * ( 2 + ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 ))
   allocate (work(lwork))
   allocate (a(n,n))
   a = matrix

   call DSYEV( JOBZ, UPLO, N, A, LDA, eigenvalues, WORK, LWORK, INFO )

   if (present(eigenvectors)) eigenvectors = a

   end subroutine

end module
