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
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck TRDR3T */
      SUBROUTINE TRDR3T(WORK,KFREE,LFREE,IPRINT,INTFLG,                 &
     &                  NSTR,ANTIS,LMP2,INDX,KQ,KE,KIBE,DINTSKP)
!
!     Written by Luuk Visscher, december 1996
!     Removed by Luuk Visscher, july 2002
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      LOGICAL ANTIS, LMP2, NOPV, NODV, TRIAN(2)
      DIMENSION NSTR(2,0:2,4),INDX(3,*)
      DIMENSION KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION WORK(*)
!
      CALL QUIT ("ASKED FOR OBSOLETE STRATEGY 3")
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck MS4INW*/
      SUBROUTINE MS4INW (WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,IOFF4,       &
     &                   IC,IHM,MHTM,IHTM,INDX,INDXAB,IPASS,            &
     &                   HMAT,IRECC,JBOFF,IBOFF,NBUF,NIJBUF,NBFSZ,      &
     &                   RBUF,LBUF1)

!
!     Sort half-transformed integrals
!     This version writes them to disk
!
!     Luuk Visscher
!
!
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
!
      PARAMETER (HALF=0.5D0)
!
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5),IOFF4(0:7)
      DIMENSION IHTM(0:7),IHM(0:7),MHTM(2)
      DIMENSION HMAT(*)
      DIMENSION WORK(*)
      DIMENSION IND(2)
      LOGICAL SUB3,SUB4
!
!     For the buffered output
!
      DIMENSION IRECC(2),IBOFF(2),NBFSZ(2),NBUF(2),NIJBUF(2)
      DIMENSION JBOFF(2,0:*),RBUF(*),LBUF1(*)
!
#include "ibtfun.h"
!
!     Find out the component type of the functions
!
      CALL ICTYPES(IC,IC1,IC2)
!
!     Loop over the functions inside the shell-block.
!     We do it this way rather than just over the boson irreps of the
!     density and then over the pairs within that irrep, because we
!     need to access the information in INDXAB. The alternative gives
!     less jumping through memory, but requires introduction of a new
!     n^2 type array in which indab is stored.
!
      DO 20 INDAR = 1, NINSHA
         DO 10 INDBR = 1, NINSHB
            INDAB  = INDXAB(INDAR,INDBR,1)
            IF (INDAB.EQ.0) GOTO 10
            IF (IPASS.NE.INDXAB(INDAR,INDBR,5)) GOTO 10
            IREPAB = INDXAB(INDAR,INDBR,3)
            IPAB   = INDXAB(INDAR,INDBR,4)
            IREPIJ = JBTOF(IREPAB,1)
            IF (IC.EQ.3) IREPIJ = JBTOF(IREPAB,2)
            IOFF = (IHM(IREPAB)+(IPAB-1)*NFPCK12(IREPIJ))*NZ+1
            CALL IUNPCK(INDAB,2,IND)
            INDA = IND(1)
            INDB = IND(2)
            IREPA = INDX(2,INDA)
            IREPB = INDX(2,INDB)
            INDA1 = INDA - ICOS(IREPA+1,IC1) ! Absolute position
            INDB1 = INDB - ICOS(IREPB+1,IC2) ! Absolute position
            IF (SUB3) INDA1 = INDX(3,INDA)  ! Relative position
            IF (SUB4) INDB1 = IOFF4(IREPB) + INDX(3,INDB)  ! Relative position
!
!           The offset within the untransformed indexes. These were
!           the right-hand side (slowest varying) and are sorted to
!           become the left hand side. Symmetrization is done
!           in the second half-transformation (the function indices are
!           left triangular), but we scale the diagonal elements to 
!           avoid double counting.
!
            JOFFAB = ISPCK34(IREPA,IREPB,IC)                            &
     &             + (INDB1-1)*NBBAS3(IREPA,IC1)+INDA1
            JOFF1 = IHTM(IREPAB) + JOFFAB
            IF (INDA.EQ.INDB)                                           &
     &         CALL DSCAL (NZ*NFPCK12(IREPIJ),HALF,HMAT(IOFF),1)
!
!           Number of pairs in one buffer
!
            NIJ  = NIJBUF(IREPIJ)
            NIJL = NFPCK12(IREPIJ) - NIJ*(NBUF(IREPIJ)-1)
            IF(NIJL.EQ.0) GOTO 10
!
            DO IZ = 1, NZ
               IJOFF = 0
!
!              Distribute the integrals over the output buffers
!              The buffers are kept in 1-dimensional arrays with
!              the offsets stored in JBOFF. This allows a variable
!              buffer length for the different irreps. The amount of
!              integrals in each buffer will be the same because we
!              do not check for zeroes.
!              The number of integrals in the buffers is JBOFF(IREPIJ,0)
!
               DO IBUF = 1, NBUF(IREPIJ)-1
                  JBUF = JBOFF(IREPIJ,IBUF)
                  CALL DCOPY(NIJ,HMAT(IOFF+IJOFF),1,RBUF(JBUF),1)
!
!                 Write the absolute addresses of the integrals
!
                  DO IJ = 1, NIJ
                     LBUF1(JBUF) = JOFF1+(IJOFF+IJ-1)*MHTM(IREPIJ)
                     JBUF = JBUF + 1
                  ENDDO
                  JBOFF(IREPIJ,IBUF) = JBUF
                  IJOFF = IJOFF + NIJ
               ENDDO
               JBOFF(IREPIJ,0) = JBOFF(IREPIJ,0) + NIJ
!
!              Fill the last buffer : the number of integrals may be
!              different
!
               IBUF = NBUF(IREPIJ)
               JBUF = JBOFF(IREPIJ,IBUF)
               CALL DCOPY(NIJL,HMAT(IOFF+IJOFF),1,RBUF(JBUF),1)
               DO IJ = 1, NIJL
                  LBUF1(JBUF) = JOFF1+(IJOFF+IJ-1)*MHTM(IREPIJ)
                  JBUF = JBUF + 1
               ENDDO
               JBOFF(IREPIJ,IBUF) = JBUF
!
               IOFF  = IOFF  + NFPCK12(IREPIJ)
               JOFF1 = JOFF1 + NSPCK34(IREPAB,IC)
            ENDDO
!
!           When the buffers are full we write them to disk
!           In this case JBUF has become equal to the length
!           of the buffer array for this irrep.
!
            IF (JBOFF(IREPIJ,0).GT.NBFSZ(IREPIJ)) THEN
               CALL QUIT ('Buffer overflow for half-transformed ints')
            ELSEIF (JBOFF(IREPIJ,0).EQ.NBFSZ(IREPIJ)) THEN
!
!              Now write the buffers : write always the full length,
!              but write also the actual length in the first spot
!
               IREC = IRECC(IREPIJ)
               JBUFF = IBOFF(IREPIJ)
               DO IBUF = 1, NBUF(IREPIJ)
                  NBUFT = JBOFF(IREPIJ,IBUF) - JBUFF
                  JBUFL = JBUFF + NBFSZ(IREPIJ) - 1
                  WRITE (LUTRA1+IREPIJ,REC=IREC) NBUFT,                 &
     &                  (RBUF(JBUF),JBUF=JBUFF,JBUFL),                  &
     &                  (LBUF1(JBUF),JBUF=JBUFF,JBUFL) 
!
!                 Reset offsets 
!
                  JBOFF(IREPIJ,IBUF) = JBUFF
                  JBUFF = JBUFL + 1
                  IREC = IREC + 1
               ENDDO
               JBOFF(IREPIJ,0) = 0
               IRECC(IREPIJ) = IREC
            ENDIF
  10     CONTINUE
  20  CONTINUE
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck MS4IN2D*/
      SUBROUTINE MS4IN2D(WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR3,NSTR4,   &
     & Q3S,Q4S,Q3T,Q4T,IRECC,NBFSZ,NBUF,NIJBUF,MHTM,INDXKR,             &
     & INDXB12,INDXB34)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!     Written by Luuk Visscher January 1997.
!
!     Driver for the second index transformation
!     Get part of the half-transformed integrals from memory, call
!     MS4IN2
!
!     Input :
!
!     - IPRINT       Print flag
!     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
!     - ICL          Last class of integrals
!     - NSTR3        Number of active spinors for index 3
!     - NSTR4        Number of active spinors for index 4
!     - Q3S          Coefficients for index 3
!     - Q4S          Coefficients for index 4
!     - Q3T          Coefficients for index 3 but for transposed AO ints
!     - Q4T          Coefficients for index 4 but for transposed AO ints
!
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
!
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION Q3S(*),Q4S(*),Q3T(*),Q4T(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
      DIMENSION IRECC(2),MHTM(2),NBUF(2),NIJBUF(2),NBFSZ(2)
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbibt.h"
!
      CALL QENTER('MS4IN2D')
!
!     The integrals are stored with their absolute address, IOFF is 
!     the start of the current batch. The relative address is obtained 
!     by subtracting IOFF.
!
!     NIJOFF is the offset for INDXKR (one list for all irreps)
!     IJS is the offset for the start of the batch within an irrep
!
      IOFF = 0
      NIJOFF = 1
      NKLOFF = 1
      DO IREPIJ = 1, NFSYM
         IJS = 0
         NIJ  = NIJBUF(IREPIJ)
         NIJL = NFPCK12(IREPIJ) - NIJ*(NBUF(IREPIJ)-1)
         NKL  = NFPCK34(IREPIJ)
         CALL MEMGET('REAL',KHMAT,NIJ*MHTM(IREPIJ),WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KRBUF,NBFSZ(IREPIJ)*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KLBUF1,NBFSZ(IREPIJ)*NZ,WORK,KFREE,LFREE)
         DO IBUF = 1, NBUF(IREPIJ)
            IF (IBUF.EQ.NBUF(IREPIJ)) NIJ = NIJL
            CALL DZERO(WORK(KHMAT),NIJ*MHTM(IREPIJ))
            CALL READBUF(IREPIJ,IBUF,IOFF,IRECC(IREPIJ),NBFSZ(IREPIJ),  &
     &                   NBUF(IREPIJ),WORK(KRBUF),WORK(KLBUF1),         &
     &                   WORK(KHMAT))
            CALL MS4IN2P(WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR3,NSTR4,   &
     &                   IREPIJ,IJS,NIJ,Q3S,Q4S,Q3T,Q4T,                &
     &                   WORK(KHMAT),INDXKR(1,NIJOFF),                  &
     &                   INDXB12(1,NIJOFF),INDXB34(1,NKLOFF))
            IOFF = IOFF + NIJ * MHTM(IREPIJ)
            IJS = IJS + NIJ
            NIJOFF = NIJOFF + NIJ
         ENDDO
         NKLOFF = NKLOFF + NIJ
         CALL MEMREL('MS4IN2D',WORK,1,KHMAT,KFREE,LFREE)
      ENDDO
!
      CALL QEXIT('MS4IN2D')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck MS4IN2P*/
      SUBROUTINE MS4IN2P(WORK,KFREE,LFREE,IPRINT,ICS,ICF,NSTR3,NSTR4,   &
     & IREPIJ,IJS,NIJ,Q3S,Q4S,Q3T,Q4T,HMAT,INDXKR,INDXB12,              &
     & INDXB34)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!     Written by Luuk Visscher January 1997.
!
!     PURPOSE : Do 4-index transformation to molecular spinor basis
!               Second index-pair transformation.
!
!     Input :
!
!     - IPRINT       Print flag
!     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
!     - ICL          Last class of integrals
!     - NSTR3        Number of active spinors for index 3
!     - NSTR4        Number of active spinors for index 4
!     - IREPIJ       Fermion symmetry (parity)
!     - IJS          First blocks in this pass
!     - NIJ          Number of blocks in this pass
!     - Q3S          Coefficients for index 3
!     - Q4S          Coefficients for index 4
!     - Q3T          Coefficients for index 3 but for transposed AO ints
!     - Q4T          Coefficients for index 4 but for transposed AO ints
!     - HMAT         Symmetry packed half-transformed integrals
!
!     Output is written directly in the MOLFDIR-type file MDCINT
!
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
!
      INTEGER   NSTR3(2),NSTR4(2)
      DIMENSION WORK(*)
      DIMENSION Q3S(*),Q4S(*),Q3T(*),Q4T(*)
      DIMENSION HMAT(*)
      DIMENSION INDXKR(2,*),INDXB12(2,*),INDXB34(2,*)
      LOGICAL CORRECT_PARITY
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbibt.h"
!
      CALL QENTER('MS4IN2P')
      KFRSAV = KFREE
!
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
!
!     Number of boson symmetry for each fermion irrep (parity)
!
      NBSYMP = NBSYM/NFSYM
!
!     Transform indices 3 and 4 to molecular spinor basis
!
      NKL = NFPCK34(IREPIJ)
      LGMAT = NKL*NZ*NZ*NBSYMP
      CALL MEMGET('REAL',KGMAT,LGMAT,WORK,KFREE,LFREE)
      IF (ISTRAT.EQ.4) THEN
!        for MDSCRU and for NFSYM.gt.1
         CALL MEMGET('REAL',KGTMP,LGMAT,WORK,KFREE,LFREE)
      ELSE IF (NFSYM.GT.1) THEN
         CALL MEMGET('REAL',KGTMP,NKL*NZ,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KGTMP,     0,WORK,KFREE,LFREE)
      END IF
      NR = 0
      NS = 0
      DO IC = ICS,ICF
         CALL ICTYPES(IC,IC3,IC4)
         DO IREPRS = 0,NBSYM-1
            NR = MAX(NR,NBBAS3(IREPRS,IC3))
            NS = MAX(NS,NBBAS4(IREPRS,IC4))
         END DO
      END DO
      CALL MEMGET('REAL',KTRMAT,NR*NS,WORK,KFREE,LFREE)
!
      IOFF = 1
!
!     Loop over the number of right hand side (12) batches.
!
      DO IJ = 1, NIJ
        JGMAT = KGMAT
        IF (NFSYM.EQ.1) CALL DZERO(WORK(KGMAT),NKL*NZ*NZ*NBSYMP)
!
!       Loop over all boson irreps. On the next line we pick out the
!       ones with the right parity. Note that for Gaunt we need to take
!       into account the difference in parity between large and small
!       components by using JBTOF(..,2) instead of JBTOF(..,1).
!
        DO IREPRS = 0, NBSYM-1
          CORRECT_PARITY = IREPIJ.EQ.JBTOF(IREPRS,1)
          IF (ICS.EQ.3) CORRECT_PARITY = IREPIJ.EQ.JBTOF(IREPRS,2)
          IF (CORRECT_PARITY) THEN
!
!           Loop over the quaternion units of the first transformed
!           density (this is IZ1 in my notes).
!
            DO IZ = 1, NZ
               IF (NFSYM.GT.1) CALL DZERO (WORK(KGTMP),NKL*NZ)
!
!              Loop over the compound component (LL/SS) of the second
!              untransformed density.
!
               IF (ICS.EQ.2) IOFF = IOFF + NSPCK34(IREPRS,1)
               DO IC = ICS, ICF
                  CALL ICTYPES(IC,IC3,IC4)
!
!                 Loop over the boson symmetries of the fourth
!                 (untransformed) index. Get the fermion symmetry
!                 of that index and get the boson and fermion symmetry
!                 of the third index.
!
                  DO IREPS = 0, NBSYM-1
                     IREPR = IBTXOR(IREPRS,IREPS)
                     NR = NBBAS3(IREPR,IC3)
                     NS = NBBAS4(IREPS,IC4)
                     IF (IC3.EQ.1) THEN
                        IREPRI = IREPR
                     ELSE
                        IREPRI = IBTXOR(IXYZ,IREPR)
                     ENDIF
                     IF (IC4.EQ.1) THEN
                        IREPSI = IREPS
                     ELSE
                        IREPSI = IBTXOR(IXYZ,IREPS)
                     ENDIF
!
                     IREPK =  JBTOF(IREPR,IC3)
                     IREPL =  JBTOF(IREPS,IC4)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     IF ((NR*NS.NE.0).AND.(NK*NL.NE.0)) THEN
                        KOFF = ICMOQS(IREPK,3) + IBBAS3(IREPR,IC3)      &
     &                       - IBAS3(IREPK)
                        LOFF = ICMOQS(IREPL,4) + IBBAS4(IREPS,IC4)      &
     &                       - IBAS4(IREPL)
                        NRQ3 = NDMOQS(1,IREPK,3)
                        NCQ3 = NDMOQS(2,IREPK,3)
                        NRQ4 = NDMOQS(1,IREPL,4)
                        NCQ4 = NDMOQS(2,IREPL,4)
!
!                       Do the actual transformation. 
!
!                       Because QTRANS will generate them in a different order
!                       for cases with inversion symmetry, we write to the
!                       temporary array KGTMP and sort later on
!                       ELSE: place them directly in array KGMAT
!
                        IF (NFSYM.GT.1) THEN
                           KLOFF = KGTMP + IFPCK34(IREPK,IREPL)*NZ
                        ELSE
                           KLOFF = JGMAT
                        END IF
!
                        CALL QTRANS('AOMO','S',D1,NR,NS,NK,NL,          &
     &                       HMAT(IOFF),NR,NS,1,IPQTOQ(1,0),            &
     &                       WORK(KLOFF),NK,NL,                         &
     &                       NZ,IPQTOQ(1,IREPRS),                       &
     &                       Q3S(KOFF),NRQ3,NCQ3,NZ,IPQTOQ(1,IREPRI),   &
     &                       Q4S(LOFF),NRQ4,NCQ4,NZ,IPQTOQ(1,IREPSI),   &
     &                       WORK(KFREE),LFREE,IPRINT)
#if 0
                        if ( ij .eq. 1 ) then
                           call header('ij = 1',-1 )
                           write(6,*) 'hmat'
                           call prqmat(hmat(ioff),nr,ns,nr,ns,          &
     &                          nz,ipqtoq(1,0),6)
                           write(6,*) 'khtmp'
                           call prqmat(work(khtmp),nk,nl,nk,nl,         &
     &                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q4s'
                           call prqmat(q3s(koff),nr,nl,nrq3,ncq3,       &
     &                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q3s'
                           call prqmat(q4s(loff),ns,nk,nrq4,ncq4,       &
     &                          nz,ipqtoq(1,0),6)
                        end if
#endif
                     ENDIF
!
!                    Repeat this transformation with the transposed
!                    AO-matrix (the half-transformed integrals were
!                    triangular in the AO-indices). Diagonal elements
!                    are scaled before storage to avoid double-counting.
!
                     IREPK =  JBTOF(IREPS,IC4)
                     IREPL =  JBTOF(IREPR,IC3)
                     NK = NSTR3(IREPK)
                     NL = NSTR4(IREPL)
                     IF ((NR*NS.NE.0).AND.(NK*NL.NE.0)) THEN
                        KOFF = ICMOQT(IREPK,3) + IBBAS4(IREPS,IC4)      &
     &                       - IBAS4(IREPK)
                        LOFF = ICMOQT(IREPL,4) + IBBAS3(IREPR,IC3)      &
     &                       - IBAS3(IREPL)
                        NRQ3 = NDMOQT(1,IREPK,3)
                        NCQ3 = NDMOQT(2,IREPK,3)
                        NRQ4 = NDMOQT(1,IREPL,4)
                        NCQ4 = NDMOQT(2,IREPL,4)
!
                        IF (NFSYM.GT.1) THEN
                           KLOFF = KGTMP + IFPCK34(IREPK,IREPL)*NZ
                        ELSE
                           KLOFF = JGMAT
                        END IF
!
                        CALL MTRSP (NR,NS,HMAT(IOFF),NR,WORK(KTRMAT),NS)
                        CALL QTRANS('AOMO','S',D1,NS,NR,NK,NL,          &
     &                       WORK(KTRMAT),NS,NR,1,IPQTOQ(1,0),          &
     &                       WORK(KLOFF),NK,NL,NZ,IPQTOQ(1,IREPRS),     &
     &                       Q3T(KOFF),NRQ3,NCQ3,NZ,IPQTOQ(1,IREPSI),   &
     &                       Q4T(LOFF),NRQ4,NCQ4,NZ,IPQTOQ(1,IREPRI),   &
     &                       WORK(KFREE),LFREE,IPRINT)
#ifdef UNDEF
                        if ( ij .eq. 1 ) then
                           call header('ij = 1',-1 )
                           write(6,*) 'ktrmat'
                           call prqmat(work(ktrmat),ns,nr,ns,nr,        &
     &                          1,ipqtoq(1,0),6)
                           write(6,*) 'khtmp'
                           call prqmat(work(khtmp),nk,nl,nk,nl,         &
     &                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q3t'
                           call prqmat(q3t(koff),ns,nk,nrq3,ncq3,       &
     &                          nz,ipqtoq(1,0),6)
                           write(6,*) 'q4t'
                           call prqmat(q4t(loff),nr,nl,nrq4,ncq4,       &
     &                          nz,ipqtoq(1,0),6)
                        end if
#endif
                     ENDIF
!
                     IOFF = IOFF + NR * NS
                  ENDDO
!                 ... ENDDO IREPS loop
               ENDDO
!              ... ENDDO IC loop
!
               IF (NFSYM.GT.1) THEN
!
!                 Order the integrals with NKL as first and IZ2
!                 as second index
!
                  DO IZ2 = 1, NZ
                     DO IREPK = 1, NFSYM
                        IREPL = MOD(IREPK+IREPIJ,2) + 1
                        NK = NSTR3(IREPK)
                        NL = NSTR4(IREPL)
                        JOFF1 = JGMAT+(IZ2-1)*NKL+IFPCK34(IREPK,IREPL)
                        JOFF2 = KGTMP + IFPCK34(IREPK,IREPL)*NZ         &
     &                                + (IZ2-1)*NK*NL
                        CALL DCOPY(NK*NL,WORK(JOFF2),1,WORK(JOFF1),1)
                     ENDDO
                  ENDDO
               ENDIF
!
               JGMAT = JGMAT + NKL*NZ
            ENDDO
!           ... ENDDO IZ loop
!
          ENDIF
        ENDDO
!       ... ENDDO IREPRS loop
!
!       We can symmetrize the integral here or write it to a
!       scratch file.
!
        IF (ISTRAT.EQ.3) THEN
           IKR = INDXKR(1,IJ)
           JKR = INDXKR(2,IJ)
           CALL SYMFINT(IPRINT,IREPIJ,INDXB12(1,IJ),INDXB34,IJ,IKR,JKR, &
     &                  NSTR3,NSTR4,NKL,WORK(KGMAT))
        ELSEIF (ISTRAT.EQ.4) THEN
           CALL MDSCRU(IREPIJ,IJS+IJ,LGMAT,WORK(KGMAT),WORK(KGTMP))
        ENDIF
      ENDDO
!     ... ENDDO IJ loop
!
      CALL MEMREL('MS4IN2P',WORK,1,KFRSAV,KFREE,LFREE)
!
      CALL QEXIT('MS4IN2P')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck INITBUF */
      SUBROUTINE INITBUF(IRECC,NBFSZ,NBUF,JBOFF)
!
!     Initialize the buffer file and index arrays for the half-
!     transformed integrals
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      DIMENSION NBFSZ(2),NBUF(2),JBOFF(2,0:*),IRECC(2)
      CHARACTER*10 LFNAM,FNODE,FNAME*5
!
      DO IRREP = 1, NFSYM
         NWORDS = NBFSZ(IRREP)
!
!     The record size needs to be nwords real + nwords + 1 integer
!
         IRECL = IRECLEN(NWORDS,NWORDS+1,0)
         LUTRA11 = LUTRA1 + IRREP
         WRITE (FNAME,'(A4,I1)') 'HTIN',IRREP
         FNODE = LFNAM(FNAME)
         OPEN (UNIT=LUTRA11,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
         IRECC(IRREP) = 1
      ENDDO
!
      JBOFF1 = 1
      DO IRREP = 1, NFSYM
         DO IBUF = 1, NBUF(IRREP)
            JBOFF(IRREP,IBUF) = JBOFF1
            JBOFF1 = JBOFF1 + NBFSZ(IRREP)
         ENDDO
         JBOFF(IRREP,0) = 0
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck DELEBUF */
      SUBROUTINE DELEBUF
!
!     Delete the buffer files for the half-transformed integrals
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
!
      DO IRREP = 1, NFSYM
         LUTRA11 = LUTRA1 + IRREP
         CLOSE (UNIT=LUTRA11,STATUS='DELETE')
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck FINABUF */
      SUBROUTINE FINABUF(IRECC,NBFSZ,NBUF,IBOFF,JBOFF,RBUF,LBUF1)
!
!     Write the last set of (half-filled) buffers
!
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      DIMENSION IRECC(2),NBFSZ(2),NBUF(2),JBOFF(2,0:*),IBOFF(2)
      DIMENSION RBUF(*),LBUF1(*)
!
      DO IREPIJ = 1, NFSYM
         IF (JBOFF(IREPIJ,0).NE.0) THEN
!
            IREC = IRECC(IREPIJ)
            JBUFF = IBOFF(IREPIJ)
            DO IBUF = 1, NBUF(IREPIJ)
               JBUFT = JBOFF(IREPIJ,IBUF) - JBUFF
               JBUFL = JBUFF + NBFSZ(IREPIJ) - 1
!
               WRITE (LUTRA1+IREPIJ,REC=IREC) JBUFT,                    &
     &               (RBUF(JBUF),JBUF=JBUFF,JBUFL),                     &
     &               (LBUF1(JBUF),JBUF=JBUFF,JBUFL)
               JBOFF(IREPIJ,IBUF) = JBUFF
               JBUFF = JBUFL + 1
               IREC = IREC + 1
            ENDDO
            IRECC(IREPIJ) = IREC
            JBOFF(IREPIJ,0) = 0
         ENDIF
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck READBUF */
      SUBROUTINE READBUF(IREPIJ,IBUF,IOFF,IRECC,NBFSZ,NBUF,             &
     &                   RBUF,LBUF1,HMAT)
!
!     Read a subset of the halftransformed integrals into memory
!
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      DIMENSION RBUF(NBFSZ),LBUF1(NBFSZ),HMAT(*)
!
      DO IREC = IBUF, IRECC-1, NBUF
         READ (LUTRA1+IREPIJ,REC=IREC) NABFSZ,RBUF,LBUF1
!
!        Distribute the integrals. 
!
         DO JBUF = 1, NABFSZ
            HMAT(LBUF1(JBUF)-IOFF) = RBUF(JBUF)
         ENDDO
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck MS4IN1D*/
      SUBROUTINE MS4IN1D(WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,IOFF4,LGREC, &
     &                   ICS,ICF,IC,                                    &
     &                   TRIAN,NSTR1,NSTR2,                             &
     &                   IJPASS,Q1,Q2,                                  &
     &                   GMAT,NSIZG,MHTM,IHTM,INDX,                     &
     &                   INDXAB,IRECC,JBOFF,IBOFF,NBUF,                 &
     &                   NIJBUF,NBFSZ,RBUF,LBUF1)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!     Written by Luuk Visscher May 1997.
!
!     PURPOSE : Driver for the first index transformation.
!               Loop over batches for a given shell combination and
!               buffer the half-transformed integrals to file
!
!     Input :
!
!     - IPRINT       Print flag
!     - NPASS        Number of passes through the scalar integrals
!     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX)
!     - ICL          Last class of integrals
!     - IC           Component of the right hand : (XX|LL) or (XX|SS)
!     - TRIAN        Store the transformed integrals as a lower triangle
!     - NSTR1        Number of active spinors for index 1
!     - NSTR2        Number of active spinors for index 2
!     - IJPASS       Number of blocks for each boson symmetry
!     - Q1           Coefficients for index 1
!     - Q2           Coefficients for index 2
!     - GMAT         Symmetry packed scalar integrals
!     - NSIZG        Size of GMAT array
!     - MTHM         Number of scalar integral pairs that contribute
!                    to a given fermion ircop (multiplied with NZ)
!     - ITHM         Start address of a boson representation on the
!                    file of half-transformed integrals.
!     - INDX         For each boson function the boson irrep and the
!                    position in the shell
!     - INDXAB       Information for a pair of boson function belonging
!                    to the current shell pair, 
!     - IRECC        Current record on the half-transformed integr files
!     - JBOFF        Number of integrals in the buffers
!     - IBOFF        Starting adress for buffers of a fermion ircop
!     - NBUF         Number of buffers for half-transformed integrals
!     - NIJBUF       Size of buffers for half-transformed integrals
!     - NBFSZ        Size buffer arrays for half-transformed integrals
!     - RBUF         Buffer arrays with value of half-transformed integ
!     - LBUF1        Buffer arrays with label of half-transformed integ
!
!     Output :
!
!     -              Updated values of buffer blocks
!
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
!
      LOGICAL   TRIAN
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(NSIZG)
      DIMENSION IJPASS(0:7,NPASS),INDX(3,*),IOFF4(0:7)
      LOGICAL   SUB3,SUB4
!
!     For the buffered input
!
      DIMENSION LGREC(NPASS)
!
!     For the buffered output
!
      DIMENSION IRECC(2),IBOFF(2),NBFSZ(2),NBUF(2),NIJBUF(2)
      DIMENSION JBOFF(2,0:*),RBUF(*),LBUF1(*)
!
      DIMENSION IHM(0:7)
!
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbibt.h"
!
      CALL QENTER('MS4IN1D')
      KFRSAV = KFREE
!
      DO IPASS = 1, NPASS
!
!        Get the scalar integrals that are processed in this pass
!        The integral are already in GMAT when NPASS = 1
!
          IF (NPASS.GT.1) THEN
             CALL DZERO(GMAT,NSIZG)
             IREC = LGREC(IPASS)
             CALL READGBF(LGFIL,IREC,NGBFSZ,GMAT)
          ENDIF
!
!        Calculate the size of the H-matrix for this pass
!
         NSIZH = 0
         DO IREPAB = 0, NBSYM-1
            IHM(IREPAB) = NSIZH
            IREPIJ = JBTOF(IREPAB,1)
            NSIZH = NSIZH + IJPASS(IREPAB,IPASS)*NFPCK12(IREPIJ)
         ENDDO
!
         CALL MEMGET('REAL',KHMAT,NSIZH*NZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KHMAT),NSIZH*NZ)
!
!        Do first step of 4-index transformation :
!        Transform first pair of indices
! 
         CALL MS4IN1 (WORK,KFREE,LFREE,IPRINT,ICS,ICF,                  &
     &                TRIAN,NSTR1,NSTR2,NDMOQR,                         &
     &                ICMOQR,IJPASS(0,IPASS),Q1,Q2,                     &
     &                GMAT,WORK(KHMAT))
!
!        Scatter half-transformed integrals to the right position
!        Write them out to disk 
!
         CALL MS4INW (WORK,KFREE,LFREE,IPRINT,SUB3,SUB4,IOFF4,          &
     &                IC,IHM,MHTM,IHTM,INDX,INDXAB,IPASS,               &
     &                WORK(KHMAT),IRECC,JBOFF,IBOFF,NBUF,NIJBUF,        &
     &                NBFSZ,RBUF,LBUF1)
         CALL MEMREL('MS4IN1D',WORK,KFRSAV,KHMAT,KFREE,LFREE)
!
      ENDDO
!
      CALL QEXIT('MS4IN1D')
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck READGBF */
      SUBROUTINE READGBF(LGFIL,IREC,NGBFSZ,GMAT)
!
!     Read a batch of scalar integrals into memory
!
!     Written by Luuk Visscher, May 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
!
      Real(8)              :: GMAT(*)
      Real(8), allocatable :: RBUF(:)
      Integer, allocatable :: LBUF(:)
      Integer, intent(in)  :: NGBFSZ
      Integer :: istat

      IF (IREC.EQ.0) RETURN
! miro: provide with checking parameter
      Allocate (RBUF(NGBFSZ),stat=istat)
      If (istat.ne.0) then
        Write(LUPRI,*)                                                  &
     &  "Real*8 RBUF(NGBFSZ) allocation failed, NGBFSZ=",NGBFSZ
           Call Quit('Real*8 RBUF(NGBFSZ) allocation failed !')
      EndIf
      Allocate (LBUF(NGBFSZ),stat=istat)
      If (istat.ne.0) then
        Write(LUPRI,*)                                                  &
     &  "Integer LBUF(NGBFSZ) allocation failed, NGBFSZ=",NGBFSZ
        Call Quit('Integer LBUF(NGBFSZ) allocation failed !')
      EndIf
!
!     Read a batch and the link to the previous record.
!
    1 READ (LGFIL,REC=IREC) IRECN,NABFSZ,RBUF,LBUF
      IREC = IRECN
!
!     Distribute the integrals. 
!
      DO JBUF = 1, NABFSZ
         GMAT(LBUF(JBUF)) = RBUF(JBUF)
      ENDDO
!
!     We read in reverse order, if we have the first block we're done
!
      IF (IREC.NE.0) GOTO 1
      DeAllocate (RBUF,stat=istat)
      DeAllocate (LBUF,stat=istat)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck FINAGBF */
      SUBROUTINE FINAGBF (LGFIL,NPASS,IREC,LGREC,LGBUF,NGBFSZ,          &
     &                    RGBUF,IGBUF)
!
!     Write the last set of (half-filled) buffers of scalar integrals
!
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      DIMENSION LGREC(NPASS),LGBUF(NPASS)
      DIMENSION RGBUF(NGBFSZ,NPASS),IGBUF(NGBFSZ,NPASS)
!
      DO IPASS = 1, NPASS
         IF (LGBUF(IPASS).NE.0) THEN
!
            IREC = IREC + 1
            WRITE (LGFIL,REC=IREC) LGREC(IPASS),LGBUF(IPASS),           &
     &            (RGBUF(JBUF,IPASS),JBUF=1,NGBFSZ),                    &
     &            (IGBUF(JBUF,IPASS),JBUF=1,NGBFSZ)
            LGREC(IPASS) = IREC
            LGBUF(IPASS) = 0
         ENDIF
      ENDDO
!     print*,'Number of buffers written', IREC
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck INITGBF */
      SUBROUTINE INITGBF(LGFIL,NGBFSZ,NPASS,IREC,LGREC,LGBUF)
!
!     Initialize the buffer file and index arrays for the scalar
!     integrals
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      DIMENSION LGREC(NPASS),LGBUF(NPASS)
      CHARACTER*10 FNODE,LFNAM
!
      NWORDS = NGBFSZ
!     print*,'NGBFSZ',NGBFSZ
!
!     The record size needs to be nwords real + (nwords + 2) integer
!
      IRECL = IRECLEN(NWORDS,NWORDS+2,0)
      IF (IRECL .LE. 0) THEN
         write(lupri,*) 'ERROR: bad record length in INITGBF'
         write(lupri,*) 'NGBFSZ,NPASS,IRECL',NGBFSZ,NPASS,IRECL
         call quit('ERROR: bad record length in INITGBF')
      END IF
!
      FNODE = LFNAM('SCLIN')
      OPEN (UNIT=LGFIL,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
!
      IREC = 0
      CALL IZERO(LGREC,NPASS)
      CALL IZERO(LGBUF,NPASS)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck DELGBUF */
      SUBROUTINE DELGBUF(LGFIL)
!
!     Delete the buffer files for the scalar integrals
!     Written by Luuk Visscher, April 1997
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
!
      CLOSE (UNIT=LGFIL,STATUS='DELETE')
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!/* Deck ICTYPES */
      SUBROUTINE ICTYPES(IC,IC1,IC2)
!
!     Get the individual components from the compound IC label
!     Written by Luuk Visscher, July 2002
!
      IF (IC.LE.2) THEN
         IC1 = IC
         IC2 = IC
      ELSEIF (IC.EQ.3) THEN
         IC1 = 2
         IC2 = 1
      ENDIF
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck gauntr */
      SUBROUTINE GAUNTR(FADD,NRAO,NCAO,NRMO,NCMO,                       &
     &                  FAO,LRAO,LCAO,NZAO,IQAO,                        &
     &                  FMO,LRMO,LCMO,NZMO,IQMO,                        &
     &                  TM1,LR1,LC1,NZTM1,IQTM1,                        &
     &                  TM2,LR2,LC2,NZTM2,IQTM2,                        &
     &                  BUF,LBUF,IPRINT)
!*****************************************************************************
!
!     Dedicated routine based on QTRANS to do transformation of a vector
!     of spin matrices as required in the transformation of the Gaunt
!     integrals. Only AO-to-MO-transformation.
!
!     ADDF gives the possibility of adding the results
!
!     The unitary C-matrix has a time-symmetric structure.
!
!     The AO-to-MO-transformation can be set up quaternionically as
!
!       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
!
!       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
!
!     Written by L.Visscher : July 2002
!****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
!
      DIMENSION FAO(LRAO,LCAO,*),IQAO(*),                               &
     &          FMO(LRMO,LCMO,NZMO,3),IQMO(*),                          &
     &          TM1(LR1,LC1,*),IQTM1(*),                                &
     &          TM2(LR2,LC2,*),IQTM2(*),                                &
     &          BUF(LBUF)
      DIMENSION IQW(4),IQWA(4)
!
      CALL QENTER('GAUNTR')
!
!     *******************************************************
!     *****  AO - to - MO -transformation : F' = (C+)FC *****
!     *******************************************************
!
!
!     Look for the most efficient transformation
!       ------------------------------------------
!
      IF(NCAO.GE.NRAO) THEN
!
!       Perform the two steps:
!         1. W  = FC
!         2. F' = (C+)W
!
        CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
        NBUF = NRAO*NCMO*NZW
        IF(NBUF.GT.LBUF) THEN
           WRITE (LUPRI,'(/A,I10/A,I10)')                               &
     & ' >>> GAUNTR error, need work space    ',NBUF,                   &
     & '                   current work space ',LBUF
        ENDIF
!
!       First part of transformation: W = FC
!       ------------------------------------
!
        CALL QGEMM(NRAO,NCMO,NCAO,D1,                                   &
     &           'N','N',IQAO,FAO,LRAO,LCAO,NZAO,                       &
     &           'N','N',IQTM2,TM2,LR2,LC2,NZTM2,                       &
     &      D0,IQW,BUF,NRAO,NCMO,NZW)
!
!       Second part of transformation: H = (C+)W
!       ----------------------------------------
!
        DO IALPHA = 2, 4
!
!         Multiply the one-index transformed matrix
!         by the quaternion unit corresponding to the
!         alpha matrix.
!
          CALL IQPACK(IALPHA,1,IQW,NZW,IQWA,NZWA)
          CALL QGEMM(NRMO,NCMO,NRAO,D1,                                 &
     &       'H','N',IQTM1,TM1 ,LR1,LC1,NZTM1,                          &
     &       'N','N',IQWA,BUF,NRAO,NCMO,NZWA,                           &
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
        ENDDO
      ELSE
!
!       Perform the two steps:
!         1. W  = (C+)F
!         2. F' = WC
!
        CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
        NBUF = NRMO*NCAO*NZW
        IF(NBUF.GT.LBUF) THEN
           WRITE (LUPRI,'(/A,I10/A,I10)')                               &
     & ' >>> GAUNTR error, need work space    ',NBUF,                   &
     & '                   current work space ',LBUF
        ENDIF
!
!       First part of transformation: W = (C+)F
!       ----------------------------------------
!
        CALL QGEMM(NRMO,NCAO,NRAO,D1,                                   &
     &     'H','N',IQTM1,TM1,LR1,LC1,NZTM1,                             &
     &     'N','N',IQAO,FAO,LRAO,LCAO,NZAO,                             &
     &      D0,IQW,BUF,NRMO,NCAO,NZW)
!
!       Second part of transformation: H = WC
!       -------------------------------------
!
        DO IALPHA = 2, 4
          CALL IQPACK(IQW,NZW,IALPHA,1,IQWA,NZWA)
          CALL QGEMM(NRMO,NCMO,NCAO,D1,                                 &
     &       'N','N',IQWA,BUF,NRMO,NCAO,NZWA,                           &
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,                           &
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
!
        ENDDO
      ENDIF
!
      CALL QEXIT('GAUNTR')
!
      RETURN
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

