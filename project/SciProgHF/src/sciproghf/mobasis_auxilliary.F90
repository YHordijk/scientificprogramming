module datatypes
       type, public  :: one_el_t
         integer     :: n_spinor=-1           ! number of spinors, used for consistency checks
         real(8)     :: e_core                ! energy of the core electrons plus nuclear repulsion energy
         complex(8), pointer :: h_core(:,:)   ! one-body part of Hamiltonian (kinetic energy + nuclear attraction + core shielding)
       end type one_el_t
end module datatypes

module diagonalization

      contains 

      SUBROUTINE QDIAG(NZ,N,A,LRA,LCA,EIG,MATZ,VEC,LRV,LCV,IERR)
!***********************************************************************
!
!       This is a driver routine for the diagonalization of
!         * real symmetric matrices                     NZ = 1
!         * complex Hermitian matrices                  NZ = 2
!         * quaternion Hermitian matrices               NZ = 4
!
!       INPUT:
!
!       A       matrix to be diagonalized
!
!       Control parameter MATZ:
!               MATZ = 0        Only eigenvalues desired.
!               MATZ = 1        Eigenvalues and eigenvectors
!
!       OUTPUT:
!
!       EIG     eigenvalues in ascending order.
!       VEC     eigenvectors when MATZ=1.
!
!       Error parameter IERR:
!               IERR = 0        Normal completion
!               IERR.NE.0       Erroneous completion
!
!       TEMPORARY STORAGE:      FV1,FV2,TAU
!       Written by T.Saue November 1994 - Odense
!       OMP paralllelization by H. J. Aa. Jensen, F90 conversion L. Visscher
!***********************************************************************
!
      implicit none
      real(8) A(LRA,LCA,NZ),EIG(N),VEC(LRV,LCV,NZ)
      real(8), allocatable :: FV1(:),FV2(:),TAU(:)
      integer :: i, nz, n, lra, lca, matz, lrv, lcv, ierr

      allocate(FV1(N))
      allocate(FV2(N))
      allocate(TAU(N*NZ))
!
!     Householder transform to real symmetric tridiagonal matrix
!     ==========================================================
!
      IF (NZ.EQ.4) THEN
        CALL  QTRIDI_omp(A,N,LRA,LCA,EIG,FV1,FV2,TAU)
      ELSEIF(NZ.EQ.2) THEN
        CALL  HTRIDI_omp(A,N,LRA,LCA,EIG,FV1,FV2,TAU)
      ELSEIF(NZ.EQ.1) THEN
          CALL  TRED2(LRA,N,A(1,1,1),EIG,FV1,VEC(1,1,1))
      ELSE
        STOP 'QDIAG: Illegal value of NZ !'
      ENDIF
 
!     Initialize vectors as unit matrix
      VEC = 0.D0
      DO I = 1, N
         VEC(I,I,1) = 1.D0
      END DO
 
!     Find eigenvalues/vectors of real symmetric
!     tridiagonal matrix
 
      CALL TQL2_omp(LRV,N,EIG,FV1,VEC(1,1,1),IERR)
 
!     Backtransform to find eigenvectors of
!     quaternion Hermitian matrix
 
      IF    (NZ.EQ.4) THEN
        CALL  QTRIBK_omp(A,N,LRA,LCA,TAU,VEC,N,LRV,LCV)
      ELSEIF(NZ.EQ.2) THEN
        CALL  HTRIBK_omp(A,N,LRA,LCA,TAU,VEC,N,LRV,LCV)
      ENDIF

!     Memory deallocation
      deallocate (FV1,FV2,TAU)
 
      RETURN
 
      END

end module diagonalization
