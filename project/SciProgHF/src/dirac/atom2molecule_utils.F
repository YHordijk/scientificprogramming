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
!
!
! this file contains all the functionality required for atom-to-molecule works.
!
! routines written by tsaue march 2011, collected by sknecht sep 2012
!
!***********************************************************************
      SUBROUTINE PUT_AOBLOCK(AMOL,NBAS_R,NBAS_C,INUC,NONT,NDEG,NZC1,
     &                       ATOM,NR,LR,IR_CLS,NC,LC,IC_CLS,
     &                       rbuf, cbuf, ibuf, jbuf)
C***********************************************************************
C
C     Put an atomic block into a molecular matrix
C       AMOL        - molecular matrix
C       NBAS_R      - row dimension of molecular matrix
C       NBAS_C      - column dimension of molecular matrix
C       INUC        - index of first atom of this type in list of
C                     symmetry-independent centers
C       NONT        - number of symmetry-independent centers of this type
C       NDEG(ICENT) - degeneracy of symmetry-independent center ICENT
C
C       ATOM        - atomic matrix with leading row LR and column LC
C       NR          - real rows of atomic matrix
C       NC          - real columns of atomic matrix            
C
C     Written by T. Saue March 31 2011
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION AMOL(NBAS_R,NBAS_C,NZC1),ATOM(LR,LC,NZC1),NDEG(*)
      integer, intent(inout) :: rbuf(*), cbuf(*), ibuf(*), jbuf(*)

      ICENT=INUC
!     print *, 'nbas_r,nbas_c,lr,lc',nbas_r, nbas_c,lr,lc,nont
      DO ITYP = 1,NONT
!     print *, 'ndeg(ityp)',ndeg(ityp)
        DO IDEG = 1,NDEG(ITYP)
C.........Create pointer mol->atom for rows
          CALL LABCOUNT(N1,RBUF,NBAS_R,1,IR_CLS,ICENT,IDEG,-1)
          IF(N1.NE.NR) THEN
            WRITE(6,'(A,2(/A,I5))') 
     &                'PUT_AOBLOCK: Error in row number',
     &                'Expected  : ',NR,
     &                'Found     : ',N1
            CALL QUIT('PUT_AOBLOCK: Error in row number !')
          ENDIF
C.........Create reverse pointer atom->mol
          CALL REVINDEX(RBUF,NBAS_R,IBUF,NR)
C.........Create pointer mol->atom for colums
          CALL LABCOUNT(N2,CBUF,NBAS_C,1,IC_CLS,ICENT,IDEG,-1)
          IF(N2.NE.NC) THEN
            WRITE(6,'(A,2(/A,I5))') 
     &                'PUT_AOBLOCK: Error in column number',
     &                'Expected  : ',NC,
     &                'Found     : ',N2
            CALL QUIT('PUT_AOBLOCK: Error in column number !')
          ENDIF
C.........Create reverse pointer atom->mol
          CALL REVINDEX(CBUF,NBAS_C,JBUF,NC)
          DO IZ = 1,NZC1
            CALL IPMAT(NR,NC,AMOL(1,1,IZ),NBAS_R,ATOM(1,1,IZ),LR,
     &                 IBUF,JBUF)
!     call wrtmatmn(amol(1,1,iz),nbas_r,nbas_c,nbas_r,nbas_c,6)
          ENDDO
        ENDDO
        ICENT = ICENT + 1
      ENDDO
C
      END
!***********************************************************************

      SUBROUTINE LABCOUNT(N,LIST,NBAS,ILAB,I1,I2,I3,I4)
!***********************************************************************
!
!     Do statistics or make pointers using basis function labels
!
!     ILAB = 1: A)AO-labels
!     ---------------------------------------------------------------------
!     I1   ICLS  - Basis set class: large(1), small(2) ..
!     I2   ICENT - index of center in list of symmetry-independent centers
!     I3   IDEG  - for given symmetry-independent center, index of degeneracy
!     I4   ITYP  - basis function type, e.g. s, px, dxy etc.
!
!     ILAB = 2: SO-labels
!     ---------------------------------------------------------------------
!     I1   IC    - Basis set class: large(1), small(2) ....
!     I2   ICENT - index of center in list of symmetry-independent centers
!     I3   IRP   - index of boson irrep
!     I4   ITYP  - basis function type, e.g. s, px, dxy
!
!     T. Saue Mar 31 2011
!
!***********************************************************************
#include "implicit.h"
#include "maxorb.h"
#include "dcblab.h"
      DIMENSION LIST(NBAS)
      integer :: n, list,nbas,ilab,i1,i2,i3,i4
      integer :: ival, itemp
#include "dcbibt.h"

      n = 0
!.....Initialize index array
      CALL ICOPY(NBAS,1,0,LIST,1)
!.....If active, search first value
      IF(I1.GT.0) THEN
        N=0
        DO I = 1,NBAS
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = IGET(ITEMP)
          IF(IVAL.NE.I1) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDDO
      ENDIF
!.....If active, search second value
      IF(I2.GT.0) THEN
        N=0
        DO I = 1,NBAS
        IF(LIST(I).GT.0) THEN
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = JGET(ITEMP)
          IF(IVAL.NE.I2) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDIF
        ENDDO
      ENDIF
!.....If active, search third value
      IF(I3.GT.0) THEN
        N=0
        DO I = 1,NBAS
        IF(LIST(I).GT.0) THEN
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = KGET(ITEMP)
          IF(IVAL.NE.I3) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDIF
        ENDDO
      ENDIF
!.....If active, search third value
      IF(I4.GT.0) THEN
        N=0
        DO I = 1,NBAS
        IF(LIST(I).GT.0) THEN
          ITEMP = IATTR(IPLAB(I,ILAB),ILAB)
          IVAL  = LGET(ITEMP)
          IF(IVAL.NE.I4) THEN
            LIST(I)=-1
          ELSE
            N=N+1
          ENDIF
        ENDIF
        ENDDO
      ENDIF
 
      END
!***********************************************************************

      SUBROUTINE IPMAT(NR,NC,A,LDA,B,LDB,INDR,INDC)
!-----------------------------------------------------------------------
!     This routine will transfer elements from matrix B to matrix A
!     according to row and column pointers INDR and INDC
!
!	Written by T.Saue, March 2011
!
!-----------------------------------------------------------------------
      real(8), intent(out) :: a(lda,*)
      real(8), intent(in)  :: b(ldb,nc)
      integer, intent(in)  :: indr(nr)
      integer, intent(in)  :: indc(nc)
      integer, intent(in)  :: nc
      integer, intent(in)  :: nr
      integer, intent(in)  :: lda
      integer, intent(in)  :: ldb
!-----------------------------------------------------------------------
      integer              :: i, j, ii, jj
!-----------------------------------------------------------------------

      DO J = 1,NC
        JJ = INDC(J)
        DO I = 1,NR
          II = INDR(I)
          A(II,JJ)=B(I,J)
        END DO
      END DO

      END
!***********************************************************************

      SUBROUTINE REVINDEX(ILIST,NI,JLIST,NJ)
!-----------------------------------------------------------------------
!     Make reverse indexing for elements of ILIST.GE.0
!     Written by T. Saue March 31 2011
!-----------------------------------------------------------------------
      integer, intent(out) :: jlist(nj)
      integer, intent(in)  :: ilist(ni)
      integer, intent(in)  :: ni
      integer, intent(in)  :: nj
!-----------------------------------------------------------------------
      integer              :: i, j
!-----------------------------------------------------------------------

      J = 0
      DO I = 1, NI
        IF(ILIST(I) > 0)THEN
          J = J + 1
          IF(J > NJ) CALL QUIT('REVINDEX: Out of bounds !')
          JLIST(J) = I
        ENDIF
      END DO

      END
!***********************************************************************
