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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mpglvc */
      SUBROUTINE MPGLVC(DMO,XMO,WORK,KFREE,LFREE,KFLG)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: Driver routine for the construction of the MP2
C              Lagrangian or construction of active virtual-positronic
C              density matrix (D_ax), depending on flag KFLG
C
C     INPUT:   DMO - second order density matrix, not sym. packed
C              KFLG - 1: Calculate D dependent Lagrangian
C                   - 2: Calculate D_ax
C
C     OUTPUT:  XMO - 1: The Lagrangian, not sym. packed
C                    2: D_ax, not sym packed
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
C
      CALL QENTER('MPGLVC')
C
      CALL TITLER('Density matrix dependent part Lagrangian','*',116)
C
C
C     Calculate density matrix dependent part of Lagrangian
C     =====================================================
C
      CALL LVCDEN(DMO,XMO,WORK,KFREE,LFREE,KFLG)
C
      CALL QEXIT('MPGLVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcden */
      SUBROUTINE LVCDEN(DMO,XMO,WORK,KFREE,LFREE,KFLG)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: generate density matrix dependent part of Lagrangian
C
C     INPUT:   DMO - second order density matrix
C              KFLG - 1: Calculate D dependent Lagrangian
C                   - 2: Calculate D_ax
C
C     OUTPUT:  XMO - 1: The Lagrangian, not sym. packed
C                        2: D_ax
C
C     REMARKS: 
C
C              L_(pq)=-2<kp||jq>DM_(jk) - 2<cp||bq>DM_(bc)
C
C              with: j and k active occupied orbitals
C                    b and c active virtual orbitals
C
C              pq corresponds the following type of pairs:
C
C                    active virtual - active occupied
C                    active virtual - inactive occupied
C
C                    inactive virtual - active occupied
C                    inactive virtual - inactive occupied
C
C                    positronic - active occupied
C                    positronic - inactive occupied
C              
C              using the definitions for the density matrices:
C
C              DA_(yx)(o-o)=C_(yj)(D_(jk))(C+)_(xk) for the 
C              occupied-occupied block
C                 and
C              DA_(yx)(v-v)=C_(yb)(D_(bc))(C+)_(xc) for the
C              virtual-virtual block
C
C              We can write
C
C              L_(pq)=-2[g_(pq,xy)-g_(py,xq)]DA_(yx)
C
C              with DA_(yx)=DA_(yx)(o-o) + DA_(yx)(v-v)
C
C              L_(pq) now has the same form as the G-matrix (i.e. 
C              a Fock type matrix) used in the sigma vector, which is
C              used in RPA (See eq 67 of T.Saue and H.J.Aa Jensen,
C              J. Chem. Phys. 118, 522 (2003)
C
C              -------------------------------------------------
C
C              If we have KFLG = 2:
C              D_(ax)=-<aj"||xv>DM_(vj") - <av||xj">DM_(vj")^(*)
C
C              with: j" all occupied
C                    a  active virtuals
C                    v  all virtuals + positronic
C                    x  positronic
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
#include "blocks.h"
#include "dcbfir.h"
C
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
C
      CALL QENTER('LVCDEN')
C
      KFRSAV = KFREE
C            
C     New NPOS(*) array for HER_PARDRV, reset dcbfir.h to true
C            

      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /)) 
      call SetIntTaskArrayDimension(NPOS,PARCAL)
      if (NPOS.GT.0) THEN
         CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
      else
         KPOS = KFREE
      endif
C
C     Memory allocation
C     ====================================================
      CALL MEMGET('REAL',KCMO   ,N2BBASXQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT  ,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXLAGAO,NZ*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXLAGRAN,NORBT*NORBT*NZ,WORK,KFREE,LFREE)
C
      CALL LVCDEN1(DMO,XMO,WORK(KXLAGRAN),WORK(KCMO),WORK(KDMAT),
     &             WORK(KXLAGAO),WORK(KPOS),WORK,KFREE,LFREE,KFLG)
C
C     Release memory
C     ==============
      CALL MEMREL('LVCDEN',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('LVCDEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lvcden1 */
      SUBROUTINE LVCDEN1(DMO,XMO,XLAGRAN,CMO,DMAT,XLAGAO,NPOS,
     &                   WORK,KFREE,LFREE,KFLG)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: generate density matrix dependent part of Lagrangian
C
C     INPUT:   DMO - second order density matrix
C              KFLG - 1: Calculate D dependent Lagrangian
C                   - 2: Calculate D_ax
C
C     OUTPUT:  XLAGRAN - 1: The Lagrangian, not sym. packed
C                        2: D_ax
C
C     Other arguments:
C              CMO  - spinor (mo) coefficients
C              DMAT - AO density matrix (or better to say density 
C                     matrix transformed to 'orthonormal SO' basis 
C              XLAGAO - Lagrangian transformed to AO basis
C              NPOS - something needed in parallel twofck
C
C     REMARKS:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include "cbihr2.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbdhf.h"
C
      LOGICAL TOBE
      DIMENSION NPOS(*)
      DIMENSION DMO(*)
      DIMENSION XMO(*), XLAGRAN(*)
      DIMENSION CMO(*), DMAT(*), XLAGAO(N2BBASXQ)
      DIMENSION WORK(*)
C
      CALL QENTER('LVCDEN1')
C
C     a) Get coefficients
C     -------------------
C
      LUBUF = 1
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','LVCDEN1')
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUM,IDUM,TOTERG,2)
      CLOSE(LUBUF,STATUS='KEEP')
C
C     b) Generate AO density matrix
C     -----------------------------
C
      CALL DMP2AO(DMO,DMAT,CMO,WORK,KFREE,LFREE)
C
C     c) Generate AO Lagrangian
C     -------------------------
C
      ISYMLAG = 1 ! Lagrangian is totally symmetric
      IHRMOPT = 1  ! Hermitian operator
      IFCKOPT = 1  ! 'Fock' matrix type
      NRFMAT = 1  ! nr of Fock-type matrices
C
      CALL DZERO(XLAGAO,N2BBASXQ)
C
      CALL TWOFCK(ISYMLAG,IHRMOPT,IFCKOPT,
     &            XLAGAO,DMAT,NRFMAT,NPOS,INTFLG,
     &            IPRTWO,WORK(KFREE),LFREE)
C
      IF(KFLG.EQ.1) THEN
C
C     d') Transform Lagrangian to MO basis
C     ------------------------------------
C
      CALL DZERO(XLAGRAN,NORBT*NORBT*NZ)
      CALL LTAOMO(XLAGAO,XLAGRAN,CMO,WORK,KFREE,LFREE)
C
C     add D dependent part of Lagrangian to the rest of the Lagrangian
C
      CALL DAXPY(NORBT*NORBT*NZ,D1,XLAGRAN,1,XMO,1)
C
      ENDIF
C
      IF(KFLG.EQ.2) THEN
C
C     d") Transform Density matrix to MO basis
C     ----------------------------------------
C
      CALL DTAOMO(XLAGAO,XMO,CMO,WORK,KFREE,LFREE)
C
      ENDIF
C
      CALL QEXIT('LVCDEN1')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'LVCDEN1: Coefficient file not found !'
      CALL QUIT('LVCDEN1: Coefficients not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ltaomo */
      SUBROUTINE LTAOMO(XLAGAO,XLAGRAN,CMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C
C     PURPOSE: Transform Lagrangian from AO-basis to MO-basis
C
C     INPUT:   XLAGAO - Lagrangian in AO-basis
C              CMO  - spinor (MO) coefficients
C
C     OUTPUT:  XLAGRAN - The Lagrangian in spinor (MO) basis, not
C                        symmetry packed
C
C     REMARKS:
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1=1.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION XLAGRAN(NORBT,NORBT,NZ)
      DIMENSION XLAGAO(N2BBASXQ), CMO(*)
      DIMENSION WORK(*)
C
      CALL QENTER('LTAOMO')
C
C     +--------------+
C     |   e-e part   |
C     +--------------+
C
      DO IFSYM = 1,NFSYM
C
C        Parameters for 1st index: virtuals
C
         NVIRTL= NORB(IFSYM) - NPSH(IFSYM) - NOCC(IFSYM)
         NS = NVIRTL      ! nr of virtuals for ifsym
         IS = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM) + 1
         JS = ICMOQ(IFSYM) + (NPSH(IFSYM)+NOCC(IFSYM))*NFBAS(IFSYM,0) +1
C
C        Parameters for 2nd index: occupied
C
         NI = NOCC(IFSYM) 
         II = IORB(IFSYM) + NPSH(IFSYM) + 1
         JI = ICMOQ(IFSYM) + NPSH(IFSYM)*NFBAS(IFSYM,0) + 1
C
       IF(NS.NE.0.AND.NI.NE.0) THEN 
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                XLAGAO(I2BASX(IFSYM,IFSYM)+1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                XLAGRAN(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1) ! last one is printlevel
       ENDIF
      ENDDO
C
C     +--------------+
C     |   e-p part   |
C     +--------------+
C
      DO IFSYM = 1,NFSYM
C
C        Parameters for 1st index: positrons
C
         NS = NPSH(IFSYM)      ! nr of positrons for ifsym
         IS = IORB(IFSYM) + 1
         JS = ICMOQ(IFSYM) +1
C
C        Parameters for 2nd index: occupied
C
         NI = NOCC(IFSYM)
         II = IORB(IFSYM) + NPSH(IFSYM) + 1
         JI = ICMOQ(IFSYM) + NPSH(IFSYM)*NFBAS(IFSYM,0) + 1
C
       IF(NS.NE.0.AND.NI.NE.0) THEN
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                XLAGAO(I2BASX(IFSYM,IFSYM)+1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                XLAGRAN(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1) ! last one is printlevel
       ENDIF
      ENDDO
C
      CALL QEXIT('LTAOMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dtaomo */
      SUBROUTINE DTAOMO(DAOAX,DMOAX,CMO,WORK,KFREE,LFREE)
C***********************************************************************
C        
C     Written by J. van Stralen - september 2003
C      
C     PURPOSE: Transform density matrix from AO-basis to MO-basis
C    
C     INPUT:   DAOAX - Density matrix in AO-basis
C              CMO   - spinor (MO) coefficients
C    
C     OUTPUT:  DMOAX - The Density matrix in spinor (MO) basis, not
C                      symmetry packed
C    
C     REMARKS:
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h" 
      PARAMETER(D0 = 0.0D0, D1=1.0D0)
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION DMOAX(NORBT,NORBT,NZ)
      DIMENSION DAOAX(N2BBASXQ), CMO(*)
      DIMENSION WORK(*)
C
      CALL QENTER('DTAOMO')
C
C     +--------------+
C     |   e-e part   |
C     +--------------+
C
      DO IFSYM = 1,NFSYM
C
C        Parameters for 1st index: active virtuals
C
         NS = NAVIR(IFSYM)   ! nr of active virtuals for ifsym
         IS = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM) + 1
         JS = ICMOQ(IFSYM) + (NPSH(IFSYM)+NOCC(IFSYM))*NFBAS(IFSYM,0) +1
C    
C        Parameters for 2nd index: positronic
C     
         NI = NPSH(IFSYM)
         II = IORB(IFSYM)  + 1
         JI = ICMOQ(IFSYM) + 1
C              
       IF(NS.NE.0.AND.NI.NE.0) THEN
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                DAOAX(I2BASX(IFSYM,IFSYM)+1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                DMOAX(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1) ! last one is printlevel
       ENDIF
      ENDDO
C
      CALL QEXIT('DTAOMO')
      RETURN
      END
