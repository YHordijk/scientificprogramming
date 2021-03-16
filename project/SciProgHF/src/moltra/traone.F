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

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRAONE */
      SUBROUTINE TRAONE(TOTERG,WORK,KFREE,LFREE,NSTR,NSTRT,NSPC,
     &                  KQ,KQC,KE,KIBE,NSPC2)
C***********************************************************************
C
C     Driver routine for 2-index transformations
C
C     Written by L. Visscher Mar 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "maxorb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
      LOGICAL TOBE,TRASAM
      real*8, intent(in) :: toterg
      DIMENSION WORK(*),NSTR(2,0:2,2),NSTRT(2),KQ(2,2),
     &          NSPC(2,0:2),KQC(2),KE(2,2),KIBE(2,4),
     &          NSPC2(2,0:3)
      Real(8), Allocatable :: FMO(:),ORBE(:)
      Integer, Allocatable :: IBORB(:)
C
      CALL QENTER('TRAONE')
      KFRSAV = KFREE
C
C     Get hold of the slaves if we run in parallel.
C        ( ITASK = 1 for Fock matrices )
C
      IF (PARCAL) CALL DIRAC_PARCTL( HERFCK_PAR )
C
      CALL TRAHI(IPRTRA,2)
C
C     Transform core-fock matrix to active MO basis
C     ---------------------------------------------
C
      IF(NSPC2(1,3).GT.0) THEN
        Allocate(FMO(N2BBASXQ*2)) ! core density matrix
      ELSE
        Allocate(FMO(N2BBASXQ)) ! core density matrix
      ENDIF
!Miro: used parameter DUMMY -no need for full set of MOs
      CALL GTCFCK(FMO,IPRTRA,WORK(KQ(1,1)),WORK(KQ(1,2)),
     &        WORK(KQC(1)),DUMMY,NDMOQR,ICMOQR,NSTR,NDMOQC,ICMOQC,NSPC,
     &        NSPC2,ITRA_INTFL2,WORK,KFREE,LFREE)
C
C     Write MOLFDIR file MRCONEE
C     --------------------------
C
      Allocate (ORBE(NSTRT(1)))
      Allocate (IBORB(NSTRT(1)))
C
      I_KIBE21=KIBE(2,1)
      IF (KIBE(2,1).LE.0) I_KIBE21=1 ! fix value for work
      CALL MRCONE(FMO,ORBE,WORK(KE(1,1)),WORK(KE(2,1)),
     &  IBORB,WORK(KIBE(1,1)),WORK(I_KIBE21),
     &  NSPC,NSTR,NSTRT,IPRTRA,TOTERG,WORK,KFREE,LFREE)
C
C     Release the slaves if we run in parallel.
C        ( ITASK = -1 )
C
      IF (PARCAL) CALL DIRAC_PARCTL( EXIT_NODEMENU )
C
      DeAllocate (FMO)
      DeAllocate (ORBE)
      DeAllocate (IBORB)
      CALL QEXIT('TRAONE')
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gtcfck */
      SUBROUTINE GTCFCK(FMO,IPRINT,Q1,Q2,QCO,CMO,NDMOQR,ICMOQR,NSTR,
     &          NDMOQC,ICMOQC,NSPC,NSPC2,INTFLG2,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Get core Fock-matrix in MO-basis; 
C     use active coefficients for MO-transformation
C
C     Called from traone.F/TRAONE - uses DUMMY instead of CMO
C                 trapam.F/RESOLV - uses CMO
C
C     Based on Trond Saue's GETFCK, Luuk Visscher, Mar 17 1997
C     Modified by Miro Ilias for DFT-COSCI, Auguste 2016
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
C
      LOGICAL TOBE
      DIMENSION FMO(*),WORK(*)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),NSTR(2,0:2,2)
      DIMENSION NDMOQC(2,2),ICMOQC(2),NSPC(2,0:2),NSPC2(2,0:3)
      DIMENSION Q1(*),Q2(*),QCO(*),CMO(*)
      Real(8), Allocatable :: FAO(:),DMAT(:),DMATfull(:)
C
      CALL QENTER('GTCFCK')
      KFRSAV = KFREE
C
C     Set INTFLG in dcbdhf
C
      INTFLG = INTFLG2
C
      Allocate(FAO(N2BBASXQ))
      IF(NSPC2(1,3).GT.0) THEN
        Allocate(DMAT(N2BBASXQ*2)) ! core density matrix
      ELSE
        Allocate(DMAT(N2BBASXQ)) ! core density matrix
      ENDIF
      Allocate(DMATfull(N2BBASXQ))! Miro: full density matrix
C
      CALL GTCFC1(FMO,FAO,DMAT,DMATfull,Q1,Q2,QCO,CMO,NDMOQR,ICMOQR,
     &      NSTR,NDMOQC,ICMOQC,NSPC,QCO2,NSPC2,IPRINT,WORK,KFREE,LFREE)
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
          CALL HEADER('GTCFCK:Total Fock matrix in MO-basis',-1)
          CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
!Miro: printout the core-energy
!TODO: add xc contribution !
      WRITE (LUPRI,1000) DHFERG,ELERGY,E1PART,E2PART
C
      DeAllocate(FAO)
      DeAllocate(DMAT)
      DeAllocate(DMATfull)
      CALL QEXIT('GTCFCK')
C
 1000 FORMAT (//' Core energy (includes nuclear repulsion) :',T50,F20.10
     & /' - Electronic part :',T50,F20.10
     & /'   - One-electron terms :',T50,F20.10
     & /'   - Two-electron terms :',T50,F20.10)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GTCFC1 */
      SUBROUTINE GTCFC1(
     &   FMO,FAO,DMAT,DMATfull,Q1,Q2,QCO,CMO,NDMOQR,ICMOQR,NSTR,
     &          NDMOQC,ICMOQC,NSPC,QCO2,NSPC2,IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     
C     Get core Fock-matrix in MO-basis; 
C     use active coefficients for MO-transformation
C
C     On input: QCO - core MOs 
C               CMO - full set of MOs
C               DMAT - variable for core density mtx
C               DMATfull - variable for full density matrix
C
C     Called from  GTCFCK
C
C     Based on Trond Saue's GETFC1, Luuk Visscher, Mar 17 1997
C     Miro Ilias,2016: for DFT-COSCI, add xcint contribution based of DMAT/DMATfull
C
C*****************************************************************************
         use dirac_cfg
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "aovec.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbfir.h"
#include "dcbdhf.h"
#include "blocks.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "cbihr2.h"
#include "dcbres.h"
      LOGICAL TOBE, SAVEFLAGS(4)
      LOGICAL ::  remove_spinfree_terms=.false.
      LOGICAL :: dirac_cfg_dft_calculation_save
      DIMENSION FMO(*),FAO(*),DMAT(*),DMATfull(*),WORK(*)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2),NSTR(2,0:2,2)
      DIMENSION NDMOQC(2,2),ICMOQC(2),NSPC(2,0:2),NSPC2(2,0:3)
      DIMENSION Q1(*),Q2(*),QCO(*),QCO2(*),CMO(*)

      Integer, Allocatable :: POS(:)
      real(8)              :: f_occ
C
      CALL QENTER('GTCFC1')
      KFRSAV = KFREE
      LUBUF = 22
C
C     Modified Fall 2019 by J. Pototschnig to allow for open shell core
      IF(NSPC2(1,3).GT.0) THEN
C     Initialize with open shell !
        NOPENBUF = NOPEN
        NFMAT = 2
        NOPEN = 1
C       get number of orbitals
        f_occ = 0
        DO IFRP = 1,NFSYM
          f_occ = f_occ + real(NSPC2(IFRP,0))
        ENDDO
        f_occ=real(NSPC2(1,3))/(2*f_occ)
      ELSE
C     Initialize - core is the closed-shell system !
        NOPENBUF = NOPEN
        NFMAT = 1
        NOPEN = 0
      ENDIF

C     Check if we need two-electron part
      NTOT = 0
      DO IFRP = 1,NFSYM
        NTOT = NTOT + NSPC(IFRP,0)
      ENDDO
        
      IF(NTOT.GT.0) THEN
C
C       1. Get two-electron Fock matrix in AO-basis
C       ===========================================
C
      call SaveTaskDistribFlags(saveflags)
      call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
      call SetIntTaskArrayDimension(NPOS,PARCAL)
        Allocate (POS(NPOS))

C       Modified Fall 2019 by J. Pototschnig to allow for open shell core
        IF(NSPC2(1,3).GT.0) THEN
C         with  frozen open-shell
C         ===========================
C         Get density matrix
          CALL CDENSO(DMAT,QCO,NDMOQC,ICMOQC,NSPC,NSPC2,IPRINT) 
C         Totally symmetric operator
          ISYMOP(1) = 1
          ISYMOP(2) = 1
C         Fock matrix type
          IFCKOP(1) = 1
          IFCKOP(2) = 1
C         Hermitian operator
          IHRMOP(1) = 1
          IHRMOP(2) = 1
        ELSE
C         Closed shell Hartree - Fock
C         ===========================
C         Get core density matrix
          CALL CDENS(DMAT,QCO,NDMOQC,ICMOQC,NSPC,IPRINT) 
C         Totally symmetric operator
          ISYMOP(1) = 1
C         Fock matrix type
          IFCKOP(1) = 1
C         Hermitian operator
          IHRMOP(1) = 1
        ENDIF

C       Get the two-electron Fock matrix using the density matrix
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMO,DMAT,NFMAT,
     &       POS,INTFLG,IPRINT,WORK(KFREE),LFREE)

! Miro:  control printout for DFT-COSCI
        IF(IPRINT.GE.5) THEN
          CALL HEADER('GTCFC1: Two-el.Fock matrix in AO-basis',-1)
          CALL PRQMAT(FMO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

#ifdef MOD_UNRELEASED
! Miro: inject xcint density into FMO as is in dirscf.F
        if (DORES.AND.dirac_cfg_dft_calculation.AND..NOT.NOVXC) then
          IF (COREDENS) THEN
            ! take only core-density for Vxc contribution
            call xcint_potential_rks(NTBAS(0),DMAT,FMO)
          ELSE
            ! default - take whole electronic density for Vxc contribution
            NOPEN=NOPENBUF 
            call DENMAT(DMATfull,CMO,IPRINT) 
            NOPEN=0
            call xcint_potential_rks(NTBAS(0),DMATfull,FMO)
          ENDIF

          IF(IPRINT.GE.6) THEN
           CALL HEADER(
     &'GTCFC1: Two-el.Fock matrix with xc contribution in AO-basis',-1)
            IF (COREDENS) THEN
              WRITE(LUPRI,*) 'Only core-density for xc contribution'
            ELSE
              WRITE(LUPRI,*)
     &        'Full electronic density for xc contribution'
            ENDIF
            CALL PRQMAT(FMO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
        endif
#endif

        IF(PARCAL) call SetTaskDistribFlags(saveflags)

        DeAllocate (POS)
      ELSE
        CALL DZERO(FMO, N2BBASXQ)        
        CALL DZERO(DMAT,N2BBASXQ)        
      ENDIF
C
C     3. Get one-electron Fock matrix
C     ===============================
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK1','OLD','GTCFC1')
        CALL REAFCK(LUBUF,FAO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GTCFC1: No 1-Fock found. Regenerating.'
        CALL ONEFCK(FAO,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
C     4. Calculate core-core interaction.
C     ===================================
C
!Miro for DFT-COSCI - using Kohn-Sham orbitals
! if NOVXC=.true. - do not add Exc to the core SCF energy, FMO (2-el Fock mtx) does not contain Vxc
! if NOVXC=.false. - add Exc from core density only to the core SCF energy, FMO (2-el Fock mtx) contains Vxc
!
#ifdef MOD_UNRELEASED
       if (DORES.AND.dirac_cfg_dft_calculation.AND.NOVXC) then
          dirac_cfg_dft_calculation_save=.true.
          dirac_cfg_dft_calculation=.false.   !Miro: do we have interest in absolute energies of multiplet states ?
       else 
          dirac_cfg_dft_calculation_save=.false.
       endif
#endif
      CALL ERGCAL(FAO,FMO,DMAT,WORK(KFREE),LFREE)
!Miro: restore the DFT status
#ifdef MOD_UNRELEASED
       if (dirac_cfg_dft_calculation_save.AND.NOVXC) then
          dirac_cfg_dft_calculation=.true.   !Miro: do we have interest in absolute energies of multiplet states ?
       endif
#endif

C
C     IF LVCON is used we should take the repulsion between all charge while
C     for the core-core interaction we leave out the valence part. Add on the
C     valence-valence SS-repulsion by recalculating CORRLV
C
C     LV+OF 20-8-2002 : To be tested !
C 
C     IF (LVCON) THEN
C        CORRLV_CORE = CORRLV
C        CALL GENDEN(DMAT,CMO,1,IPRSCF)
C        CALL LVCORR(DMAT,WORK,LWORK,IPRSCF)
C        DHFERG = DHFERG + CORRLV - CORRLV_CORE
C     ENDIF
C
C     5. Add 1- and 2-electron matrix and transform to MO-basis
C     =========================================================
C
      IF(NTOT.GT.0) THEN
        CALL DAXPY(N2BBASXQ,D1,FMO,1,FAO,1) 
        IF(NFMAT.GT.1) THEN
C         for now only one froze open-shell included
          CALL DAXPY(N2BBASXQ,f_occ,FMO(1+N2BBASXQ),1,FAO,1) 
        ENDIF
      ENDIF
C
C     Make the matrix spinfree if desired.
C
C     Note that we did not make the matrices spinfree prior to the energy
C     calculation since this would require one extra large (fock-matrix size)
C     array. For closed shell cores we thus use that the first order spin-orbit
C     effect is zero. This assumption should be checked when implementing
C     non-closed shell cores or unrestricted DHF algorithms !!
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTCFC1:Total Fock matrix in AO-basis',-1)
        CALL PRQMAT(FAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
CMI   ... in the BSS mode the spin-free AO matrix is ready on ONEFCK
!     nothing to be done for X2Cmmf (X2CMMF == T)
      if(spinfr) remove_spinfree_terms = .true.
      if(bss)    remove_spinfree_terms = .false.
      if(X2CMMF) remove_spinfree_terms = .false.
 
      if(remove_spinfree_terms)then
        CALL SPFAO (FMO,FAO,WORK(KFREE),LFREE)
      end if
C
      IF(IPRINT.GE.5.AND.SPINFR) THEN
        CALL HEADER('GTCFC1:Total spinfree Fock matrix in AO-basis',-1)
        CALL PRQMAT(FAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
      CALL DZERO(FMO,N2ORBXQ)
#ifdef debug_mmf
!     stefan DEBUG!!!!
      CALL DZERO(FAO,N2BBASXQ)
!     stefan DEBUG!!!!
#endif
      DO 10 I = 1,NFSYM
        IF(NSTR(I,0,1).EQ.0.OR.NSTR(I,0,2).EQ.0) GOTO 10
        NRQ1 = NDMOQR(1,I,1)
        NCQ1 = NDMOQR(2,I,1)
        NRQ2 = NDMOQR(1,I,2)
        NCQ2 = NDMOQR(2,I,2)
        IF(IPRINT.GE.5) THEN
          CALL HEADER('GTCFC1:First transformation matrix',-1)
          CALL PRQMAT(Q1(ICMOQR(I,1)),NRQ1,NCQ1,NRQ1,NCQ1,NZ,
     &                IPQTOQ(1,0),LUPRI)
          CALL HEADER('GTCFC1:Second transformation matrix',-1)
          CALL PRQMAT(Q2(ICMOQR(I,2)),NRQ2,NCQ2,NRQ2,NCQ2,NZ,
     &                IPQTOQ(1,0),LUPRI)
        ENDIF
        CALL QTRANS('AOMO','S',D0,
     &            NFBAS(I,0),NFBAS(I,0),NSTR(I,0,1),NSTR(I,0,2),
     &            FAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            FMO(I2ORBX(I,I)+1),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &            Q1(ICMOQR(I,1)),NRQ1,NCQ1,NZ,IPQTOQ(1,0),
     &            Q2(ICMOQR(I,2)),NRQ2,NCQ2,NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
 10   CONTINUE
C
C     Restore NOPEN
C
      NOPEN = NOPENBUF
C
      CALL QEXIT('GTCFC1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck CDENS */
      SUBROUTINE CDENS(DMAT,QCO,NDMOQC,ICMOQC,NSPC,IPRINT)
C*****************************************************************************
C
C Get the core density matrix into DMAT, from selected MO's in QCO
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DMAT(*),QCO(*)
      DIMENSION NSPC(2,0:2),NDMOQC(2,2),ICMOQC(2)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      CALL DZERO(DMAT,N2BBASXQ)
      DO 10 I = 1,NFSYM
        IF(NSPC(I,0).EQ.0) GOTO 10
        CALL DENST1(DMAT(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,D1,D0,
     &              QCO(ICMOQC(I)),NDMOQC(1,I),NDMOQC(2,I),
     &              1,NSPC(I,0),NFBAS(I,0))
   10 CONTINUE
C
C     Print section
C
      IF (IPRINT.GE.5) THEN
        CALL TITLER('Output from CDENS','*',103)
        DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &    '*** Fermion corep ',I,'/',NFSYM
          IF (IPRINT.GE.7) THEN
            write(LUPRI,*)
     &      'Core molecular orbitals, QCO:',NDMOQC(1,I),NDMOQC(2,I)
            CALL PRQMAT(QCO(ICMOQC(I)),NDMOQC(1,I),NDMOQC(2,I),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          write(LUPRI,*) 'Core density matrix, DMAT:'
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcone */
      SUBROUTINE MRCONE(FMO,ORBE,ORBE1,ORBE2,IBORB,IBORB1,IBORB2,
     &                  NSPC,NSTR,NSTRT,IPRINT,TOTERG,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Write core Fock-matrix and other information in MRCONEE format 
C     A better interface needs to be defined, this will serve for the
C     moment and has the advantage that RELCCSD and DIRRCI need not be
C     modified.
C
C     NSTRT(i) - number of orbitals for bra(i=1) and ket(i=2)
C     NSTR(ifrp,p,i) - number of orbitals for fermion ircop ifrp, 
C                      total (p=0), negative(p=1) or positive(p=2) energy,
C                      for bra(i=1) and ket(i=2)
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
C
      DIMENSION FMO(NORBT,NORBT,NZ)
      DIMENSION ORBE(*),IBORB(*),WORK(*)
      DIMENSION ORBE1(*),ORBE2(*),IBORB1(*),IBORB2(*)
      DIMENSION NSPC(2,0:2),NSTR(2,0:2,2),NSTRT(2)
      REAL*8, intent(in) :: TOTERG
    
      Real(8), Allocatable :: FMOM(:,:,:),ORBMO(:,:)
      Integer, Allocatable :: IRPMO(:,:),IRPAMO(:,:),IOCCUP(:,:)
      Integer, Allocatable :: IBSPI(:,:)
C
      CALL QENTER('MRCONE')
C
      IF (NSTRT(1).NE.NSTRT(2)) THEN
         WRITE (LUPRI,*) ' Ranges for bra and ket are different'
         WRITE (LUPRI,*) ' NSTRT bra and ket: ',NSTRT(1),NSTRT(2)
         WRITE (LUPRI,*) ' MRCONEE file is not written !!!'
         RETURN
      ENDIF
C
      Allocate (FMOM(NSTRT(1),NSTRT(1),8))
      Allocate (ORBMO(NSTRT(1),2))
      Allocate (IRPMO(NSTRT(1),2))
      Allocate (IRPAMO(NSTRT(1),2))
      Allocate (IOCCUP(NSTRT(1),2))
      Allocate (IBSPI(NSTRT(1),2))
C.....Put orbital energies and supersymmetry information of the one or two
C     fermion ircops into single arrays
      II = 0
      DO I = 1, NFSYM
       DO J = 1, NSTR(I,0,1)
         II = II + 1
         IF (I.EQ.1) THEN
            ORBE(II)  = ORBE1(J)
            IBORB(II) = IBORB1(J)
         ELSEIF (I.EQ.2) THEN
            ORBE(II) = ORBE2(J)
            IBORB(II) = IBORB2(J)
         ELSE
            CALL QUIT ('More than 2 fermion irreps not implemented') 
         ENDIF
       ENDDO
      ENDDO
C
      CALL MRCON1(FMO,ORBE,IBORB,FMOM,NSTR,NSTRT,IRPMO,
     &            IRPAMO,IOCCUP,ORBMO,IBSPI,
     &            NSPC,IPRINT,TOTERG,WORK,KFREE,LFREE)
C
      DeAllocate (FMOM)
      DeAllocate (ORBMO)
      DeAllocate (IRPMO)
      DeAllocate (IRPAMO)
      DeAllocate (IOCCUP)
      DeAllocate (IBSPI)
C
      CALL QEXIT('MRCONE')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mrcon1 */
      SUBROUTINE MRCON1(FMO,ORBE,IBORB,FMOM,NSTR,NSTRT,IRPMO,IRPAMO,
     &                  IOCCUP,ORBMO,IBSPI,NSPC,IPRINT,TOTERG,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Luuk Visscher, Mar 18 1997
C
C     Added number of boson symmetry reps (of parent group).
C      Timo Fleig, Mar 28, 2002
C
C     Miro Ilias, Febr.2007   
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbtra.h"
#include "dgroup.h"
C
C     NORDER gives the order of the Abelian subgroup used linear molecules. 
C     NORDER=16 (default) uses C_16_h for D_inf_h and C_32 for C_inf_v.
C
      PARAMETER (NORDER=16)
      LOGICAL BREIT,FORWARD
      CHARACTER*14 REPN(8)
      CHARACTER*4 REPA(4*NORDER)
      DIMENSION FMO(NORBT,NORBT,NZ),ORBE(NSTRT),IBORB(NSTRT)
      DIMENSION FMOM(NSTRT*2,NSTRT*2,2),WORK(*)
      DIMENSION IRPMO(NSTRT*2),IRPAMO(NSTRT*2),IOCCUP(NSTRT*2)
      DIMENSION ORBMO(NSTRT*2),IBSPI(NSTRT*2)
      DIMENSION NSTR(2),IPHS(4),NAEL(8),NSPC(2,0:2)
      DIMENSION MULTB(64,64), MJ_TO_REP(-2*NORDER:2*NORDER,2)
      DATA IPHS/ 1,-1,-1, 1/
      REAL*8, intent(in) :: TOTERG

      CALL QENTER('MRCON1')

      KFRSAV = KFREE
C
C
C     Fixed information. At the moment we cannot treat the Breit 
C     interaction and the fermion representations are limited to 
C     gerade and ungerade.
C     We treat the quaternion case as pseudocomplex to avoid ordering
C     problems in the CCSD case.
C
      BREIT = .FALSE.
      FORWARD = .TRUE.
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN 
C........Quaternion Abelian, no inversion (C1)
         NSYMRP = 2
         NBSYMRP = 1
         REPN(1) = '  1A          '
         REPN(2) = '  2A          '
      ELSEIF (NFSYM.EQ.1.AND.NZ.EQ.2) THEN 
C........Complex Abelian, no inversion (C2, Cs)
         NSYMRP = 2
         NBSYMRP = 2
         REPN(1) = '  1E          '
         REPN(2) = '  2E          '
      ELSEIF (NFSYM.EQ.1.AND.NZ.EQ.1) THEN
C........Real non-Abelian, no inversion (D2, C2v) -> C2
         NSYMRP = 2
         NBSYMRP = 4
         REPN(1) = '   E  1E      '
         REPN(2) = '   E  2E      '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
C........Quaternion Abelian, inversion (Ci)
         NSYMRP = 4
         NBSYMRP = 2
         REPN(1) = ' 1AG          '
         REPN(2) = ' 2AG          '
         REPN(3) = ' 1AU          '
         REPN(4) = ' 2AU          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.2) THEN
C........Complex Abelian, inversion (C2h)
         NSYMRP = 4
         NBSYMRP = 4
         REPN(1) = ' 1Eg          '
         REPN(2) = ' 2Eg          '
         REPN(3) = ' 1Eu          '
         REPN(4) = ' 2Eu          '
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.1) THEN
C........Real non-Abelian, inversion (D2h) -> C2h
         NSYMRP = 4
         NBSYMRP = 8
         REPN(1) = '  Eg 1Eg      '
         REPN(2) = '  Eg 2Eg      '
         REPN(3) = '  Eu 1Eu      '
         REPN(4) = '  Eu 2Eu      '
      ELSE
         PRINT*, 'Warning representation name is not defined'
      ENDIF
C
C     Get multiplication table for the abelian subgroup
C
      IF (SPINFR.AND. .NOT.NOSFMU) THEN
         CALL GMULTSF (NSYMRPA,REPA,MULTB)
      ELSEIF (LINEAR) THEN
         CALL GMULTLN (NSYMRPA,REPA,MJ_TO_REP,MULTB,NORDER)
      ELSE
         CALL GMULTA (NZ,NFSYM,NSYMRPA,REPA,MULTB)
      ENDIF
C
C     Set the number of closed shell electrons that will be correlated.
C     Open shell electrons are not considered, these require user input.
C
      DO IREP = 1, NSYMRP
         IFRP = (IREP+1) / 2
         NAEL(IREP) = NISH(IFRP) - NSPC(IFRP,1)
      ENDDO
C
C     The correlation codes (RELCCSD, DIRRCI, LUCITA) do not use 
C     time reversal symmetry, write Kramers pairs as spinors:
C     For each fermion ircop, all unbarred orbitals come first,  followed by the partners
C
C     For each orbital the following information is collected:
C        IRPMO  - fermion irrep of parent group
C        IRPAMO - fermion irrep of Abelien (sub)group
C        ORBMO  - orbital energy
C        IBSPI  - supersymmetry (if relevant)
C        IOCCUP - occupied(1)/virtual(0)


      II    = 0
      IK    = 0
      IREP  = 0
      IREPA = 0
      DO I = 1, NFSYM ! loop over fermion ircops
         IREP  = IREP  + 1
         IREPA = IREPA + 1
C........Loop over first Kramers partner
         DO ISTR = 1, NSTR(I) ! total number of orbitals of this ircop
            II = II + 1
            IK = IK + 1
            IRPMO(II)  = IREP
            IF (SPINFR.AND. .NOT.NOSFMU) THEN
               IRPAMO(II) = IBORB(IK) + 1
            ELSEIF (LINEAR) THEN
               MJ = IBORB(IK)
               IF (MJ.GT. 2*NORDER/NFSYM) MJ = MJ - 4*NORDER/NFSYM
               IF (MJ.LT.-2*NORDER/NFSYM) MJ = MJ + 4*NORDER/NFSYM
               IRPAMO(II) = MJ_TO_REP(MJ,I)
            ELSE
               IRPAMO(II) = IREPA
            ENDIF
            ORBMO(II)  = ORBE(IK)
            IBSPI(II)  = IBORB(IK)
            IF (ISTR.LE.NAEL(IREP)) THEN
               IOCCUP(II) = 1
            ELSE
               IOCCUP(II) = 0
            ENDIF
         ENDDO
         IREP = IREP + 1
         IF (NZ.NE.4) IREPA = IREPA + 1
C........Loop over second Kramers partner; reset index IK over partners
         IK = IK - NSTR(I)
         DO ISTR = 1, NSTR(I)
            II = II + 1
            IK = IK + 1
            IRPMO(II)  = IREP
            IF (SPINFR.AND. .NOT.NOSFMU) THEN
               IRPAMO(II) = IBORB(IK) + NBSYM + 1
            ELSEIF (LINEAR) THEN
               MJ = -IBORB(IK)
               IF (MJ.GT. 2*NORDER/NFSYM) MJ = MJ - 4*NORDER/NFSYM
               IF (MJ.LT.-2*NORDER/NFSYM) MJ = MJ + 4*NORDER/NFSYM
               IRPAMO(II) = MJ_TO_REP(MJ,I)
            ELSE
               IRPAMO(II) = IREPA
            ENDIF
            ORBMO(II)  = ORBE(IK)
            IBSPI(II)  = IBORB(IK)
            IF (ISTR.LE.NAEL(IREP)) THEN
               IOCCUP(II) = 1
            ELSE
               IOCCUP(II) = 0
            ENDIF
         ENDDO
      ENDDO
C
C     Fill in complete core-Fock matrix in complex form
C
      CALL MOLTODIR(FMO,FMOM,NSTR,NSTRT,FORWARD,IPRTRA)

#ifdef MOD_CAP
      CALL CAPMFDADAPT(CAP,NSTRT,NSTR,FMO,FMOM,WORK(KFREE),LFREE)
#endif
C
C     Write MRCONEE file
C
C     Pass on all the information that the post-dhf
C     electron correlation codes might need.
C
      OPEN (LUMLF1,FILE='MRCONEE',FORM='UNFORMATTED')
C     Total number of spinors, breit active in DHF ?,
C     core energy (inactive energy + nuclear repulsion), inversion symmetry (yes : 2; no : 1),
C     group type (1 real, 2 complex, 4 quat.),spinfree formalism,
C     total number of orbitals (so including frozen or deleted orbitals)
!     total SCF energy
      WRITE (LUMLF1) 2*NSTRT,BREIT,DHFERG,NFSYM,NZ,SPINFR,
     &               NORBT,TOTERG,NSPC,NCORE2
C     NSYMRP        - Number of fermion irreps in parent group
C     REPN          - names of these irreps (gerade, ungerade)
C     NAEL          - number of spinors active in the transformation
C     NSPC          - number of occupied frozen (core) spinors,
C                     0: total
C                     1: positive energy
C                     2: negative energy
C     NESH-NST-NSPC -number of deleted spinors.
C     NCORE2        - spinors missing due to frozen open shell
      WRITE (LUMLF1) NSYMRP,(REPN(IRP),IRP=1,NSYMRP),  
     &                      (NAEL(IRP),IRP=1,NSYMRP),
     &                      (NSTR(IRP),IRP=1,NFSYM),
     &                      (NSPC(IRP,0),IRP=1,NFSYM),
     &                      (NSPC(IRP,1),IRP=1,NFSYM),
     &                      (NSPC(IRP,2),IRP=1,NFSYM),
     &  (NESH(IRP)-NSTR(IRP)-NSPC(IRP,0),IRP=1,NFSYM)
C     Number of fermion irreps in the Abelian subgroup,
C     names of these irreps.
      WRITE (LUMLF1) NSYMRPA,(REPA(IRPA),IRPA=1,NSYMRPA*2)
C     Multiplication table for direct products in the
C     Abelian subgroup.
      WRITE (LUMLF1) ((MULTB(I,J),I=1,2*NSYMRPA),J=1,2*NSYMRPA)
C     Information for each spinor :
C     - Irrep in parent group (1:gerade, 2:ungerade)
C     - Irrep in Abelian subgroup
C     - Orbital energy taken from DHF
C     - Approximate boson irrep identification (needed in LUCIAREL
C       and LUCITA)
C     - Number of boson symmetry reps (for LUCITA)
      WRITE (LUMLF1) (IRPMO(I),IRPAMO(I),ORBMO(I),I=1,NSTRT*2),
     &               (IBSPI(I),I=1,2*NSTRT),(NORB(I),I=1,NFSYM),NBSYMRP
C     Matrix elements over the core Fock operator.
      WRITE (LUMLF1) ((FMOM(I,J,1),FMOM(I,J,2),
     &                 I=1,NSTRT*2),J=1,NSTRT*2)
      CLOSE (LUMLF1)
C
      WRITE (LUPRI,1000) 
C
      IF (TRA_ASCII) THEN
C     Also write some output in formatted form to facilitate interfacing with other codes
       OPEN (LUASCII,FILE='MO_integrals.txt',FORM='FORMATTED')
       WRITE (LUASCII,*) '*** Nuclear repulsion + core energy ***'
       WRITE (LUASCII,*) DHFERG
       WRITE (LUASCII,*) '*** NZ (1:real, 2:complex, 4:quaternion) ***'
       WRITE (LUASCII,*) NZ
       WRITE (LUASCII,*) '*** Number of double group irreps ***'
       WRITE (LUASCII,*) 2*NSYMRPA
       WRITE (LUASCII,*) '*** Representation names ***'
       DO IRP = 1, 2*NSYMRPA
          WRITE (LUASCII,*) IRP, REPA(IRP)
       END DO
       WRITE (LUASCII,*) '*** Multiplication table ***'
       DO JRP = 1, 2*NSYMRPA
          WRITE (LUASCII,'(64I3)') (MULTB(IRP,JRP),IRP=1,2*NSYMRPA)
       END DO
       WRITE (LUASCII,*) '*** Number of spinors ***'
       WRITE (LUASCII,*) 2*NSTRT
       WRITE (LUASCII,*) '*** Spinor, irrep, occupation, energy ***'
       DO I= 1, 2*NSTRT
          WRITE (LUASCII,*) I,
     &                      IRPAMO(I),IOCCUP(I),ORBMO(I)
       END DO
       WRITE (LUASCII,*) '*** Effective 1-e integrals ***'
       DO J = 1, NSTRT*2
          DO I = 1, NSTRT*2
           IF (NZ.EQ.1.AND.ABS(FMOM(I,J,1)).GT.1.E-15) THEN
             WRITE (LUASCII,'(E30.16,4I6)') FMOM(I,J,1),I,J,0,0
           ELSEIF (NZ.EQ.2.AND.(ABS(FMOM(I,J,1)).GT.1.E-15.OR.
     &       ABS(FMOM(I,J,2)).GT.1.E-15)) THEN
             WRITE (LUASCII,'(2E30.16,4I6)')
     &             FMOM(I,J,1),FMOM(I,J,2),I,J,0,0
           END IF
          END DO
       END DO
C     We leave the file open to write also the 2-e integrals, close after that stage
      END IF
C
      CALL QEXIT('MRCON1')
C
 1000 FORMAT (//' MOLFDIR file MRCONEE is written')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmulta */
      SUBROUTINE GMULTA (NZ,NFSYM,NSYMRPA,REPA,MULTB)
C***********************************************************************
C     
C     Hard-wired multiplication table for the abelian subgroups of D2h
C     Defines representation names as well.
C     2*NSYMRPA symmetries in total in MULTB,
C     first NSYMRPA fermion symmetries, then NSYMRPA boson symmetries.
C     
C     All quaternion and complex subgroups are Abelian.
C     The real groups are non-Abelian and so an Abelian subgroup is chosen
C       D2, C2v -> C2
C       D2h     -> C2h
C     NSYMRPA is the number of boson irreps of the subgroup;
C     since it is Abelian the number of fermion irreps will be the same
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER*4 REPA(64)
      DIMENSION MULTB(64,64)
C
      IF (NFSYM.EQ.1.AND.NZ.EQ.4) THEN
C........Quaternion Abelian, no inversion: C1
         NSYMRPA = 1
         REPA(1) = '   A' ! Fermion irrep
         REPA(2) = '   a' ! Boson irrep
         MULTB(1,1) = 2
         MULTB(2,1) = 1
         MULTB(2,2) = 2
      ELSEIF (NFSYM.EQ.1.AND.NZ.LE.2) THEN
C........Complex Abelian, no inversion (C2, Cs)
         NSYMRPA = 2
         REPA(1) = '  1E'
         REPA(2) = '  2E'
         REPA(3) = '   a'
         REPA(4) = '   b'
         MULTB(1,1) = 4
         MULTB(2,1) = 3
         MULTB(3,1) = 1
         MULTB(4,1) = 2
         MULTB(2,2) = 4
         MULTB(3,2) = 2
         MULTB(4,2) = 1
         MULTB(3,3) = 3
         MULTB(4,3) = 4
         MULTB(4,4) = 3
      ELSEIF (NFSYM.EQ.2.AND.NZ.EQ.4) THEN
C........Quaternion Abelian, inversion: Ci
         NSYMRPA = 2
         REPA(1) = '  AG' ! Fermion gerade
         REPA(2) = '  AU' ! Fermion ungerade
         REPA(3) = '  ag' ! Boson gerade
         REPA(4) = '  au' ! Boson ungerade
         MULTB(1,1) = 3
         MULTB(2,1) = 4
         MULTB(3,1) = 1
         MULTB(4,1) = 2
         MULTB(2,2) = 3
         MULTB(3,2) = 2
         MULTB(4,2) = 1
         MULTB(3,3) = 3
         MULTB(4,3) = 4
         MULTB(4,4) = 3
      ELSEIF (NFSYM.EQ.2.AND.NZ.LE.2) THEN
C........Complex Abelian, inversion (C2h)
         NSYMRPA = 4
         REPA(1) = ' 1Eg'
         REPA(2) = ' 2Eg'
         REPA(3) = ' 1Eu'
         REPA(4) = ' 2Eu'
         REPA(5) = '  ag'
         REPA(6) = '  bg'
         REPA(7) = '  au'
         REPA(8) = '  bu'
C        The gg-block
         MULTB(1,1) = 6
         MULTB(2,1) = 5
         MULTB(5,1) = 1
         MULTB(6,1) = 2
         MULTB(2,2) = 6
         MULTB(5,2) = 2
         MULTB(6,2) = 1
         MULTB(5,5) = 5
         MULTB(6,5) = 6
         MULTB(6,6) = 5
C        The uu-block
         MULTB(3,3) = 6
         MULTB(4,3) = 5
         MULTB(7,3) = 1
         MULTB(8,3) = 2
         MULTB(4,4) = 6
         MULTB(7,4) = 2
         MULTB(8,4) = 1
         MULTB(7,7) = 5
         MULTB(8,7) = 6
         MULTB(8,8) = 5
C        The ug-block
         MULTB(3,1) = 8
         MULTB(4,1) = 7
         MULTB(7,1) = 3
         MULTB(8,1) = 4
         MULTB(3,2) = 7
         MULTB(4,2) = 8
         MULTB(7,2) = 4
         MULTB(8,2) = 3
         MULTB(7,5) = 7
         MULTB(8,5) = 8
         MULTB(7,6) = 8
         MULTB(8,6) = 7
C        The gu-block
         MULTB(5,3) = 3
         MULTB(6,3) = 4
         MULTB(5,4) = 4
         MULTB(6,4) = 3
      ELSE
         CALL QUIT('GMULTA: Multiplication table is not defined')
      ENDIF
C
      DO I = 1, 2*NSYMRPA
         DO J = 1, I-1
            MULTB(J,I) = MULTB(I,J)
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmultln */
      SUBROUTINE GMULTLN (NSYMRPA,REPA,MJ_TO_REP,MULTB,NORDER)
C***********************************************************************
C
C     We make the multiplication table for linear molecules where we 
C     have used the mj quantum number to identify blocks. Here the 
C     full point group is the non-Abelian Dinfh or Cinfv and we use 
C     the Abelian Cinfh or Cinf, respectively. Since the Cinfh and 
C     Cinv have an infinite number of irreps and we use a finite 
C     multiplication table to store the coupling information we need 
C     to truncate at some point. We do this at mj=15/2 in Cinfh and 
C     mj=31/2 in Cinf so that we can keep MULTB dimensioned at 64 like 
C     before. In fact this means that we'll be using the subgroups C8h 
C     and C16 respectively. This can be easily extended or decreased 
C     using the input parameter NORDER that indicates how large MULTB 
C     may become.
C
C     Luuk Visscher, Nov 03 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
C
      CHARACTER*4 REPA(4*NORDER)
      CHARACTER PARITY
      DIMENSION MULTB(64,64), MJ_TO_REP(-2*NORDER:2*NORDER,2)
C
#include "dcbibt.h"
C
      IF (4*NORDER.GT.64) CALL QUIT ("increase dimension of multb")
      NSYMRPA = 2*NORDER
C
      IF (NFSYM.EQ.1) THEN
C
C        Start by specifying the names of the irreps
C        and making the array MJ_TO_REP.
C
         IREP = 0
C
C        Fermion (half-integer) irreps
C
         DO MJ = 1, 2*NORDER-1, 2
            IREP = IREP+1
            WRITE (REPA(IREP),1000) MJ
            MJ_TO_REP(MJ,1) = IREP
            IREP = IREP+1
            WRITE (REPA(IREP),1000) -MJ
            MJ_TO_REP(-MJ,1) = IREP
         ENDDO
C
C        Boson (integer) irreps
C
         MJ = 0
         IREP = IREP + 1
         WRITE (REPA(IREP),1000) MJ
         MJ_TO_REP(MJ,1) = IREP
         DO MJ = 2, 2*NORDER-2, 2
            IREP = IREP+1
            WRITE (REPA(IREP),1000) MJ
            MJ_TO_REP(MJ,1) = IREP
            IREP = IREP+1
            WRITE (REPA(IREP),1000) -MJ
            MJ_TO_REP(-MJ,1) = IREP
         ENDDO
         MJ = 2*NORDER
         IREP = IREP + 1
         WRITE (REPA(IREP),1000) MJ
         MJ_TO_REP( MJ,1) = IREP
         MJ_TO_REP(-MJ,1) = IREP
C
C        The multiplication table is straigthforward since 
C        the product irrep is given by the sum of Mj values.
C        The only complication is the finite range of irreps so
C        we may need to add or subtract something, e.g. :
C        MJ = 33/2 -> MJ = -31/2 in the C16 subgroup.
C        (formulated in another way MJ = +/- 16 transforms
C         according to  the totally symmetric irrep in C16)
C        
          DO MJ = -2*NORDER, 2*NORDER, 1
             DO MI = -2*NORDER, 2*NORDER, 1
                MIJ = MI + MJ 
C               Map Cinf to C16
                IF (MIJ.LT.-2*NORDER) MIJ = MIJ + 4*NORDER
                IF (MIJ.GT. 2*NORDER) MIJ = MIJ - 4*NORDER
                IF (MIJ.EQ.-2*NORDER) MIJ = 2*NORDER
                IREP = MJ_TO_REP(MI,1)
                JREP = MJ_TO_REP(MJ,1)
                IJREP = MJ_TO_REP(MIJ,1)
                MULTB(IREP,JREP) = IJREP
             ENDDO
          ENDDO
C
      ELSE
C
C        Same procedure but we now also have inversion symmetry
C
C        Fermion (half-integer) irreps
C
         IREP = 0
         DO IFSYM = 1, NFSYM
            IF (IFSYM.EQ.1) THEN
               PARITY = 'g'
            ELSE
               PARITY = 'u'
            ENDIF
            DO MJ = 1, NORDER-1, 2
               IREP = IREP+1
               WRITE (REPA(IREP),1001) MJ,PARITY
               MJ_TO_REP(MJ,IFSYM) = IREP
               IREP = IREP+1
               WRITE (REPA(IREP),1001) -MJ,PARITY
               MJ_TO_REP(-MJ,IFSYM) = IREP
            ENDDO
         ENDDO
C
C        Boson (integer) irreps
C
         DO IFSYM = 1, NFSYM
            IF (IFSYM.EQ.1) THEN
               PARITY = 'g'
            ELSE
               PARITY = 'u'
            ENDIF
            MJ = 0
            IREP = IREP + 1
            WRITE (REPA(IREP),1001) MJ,PARITY
            MJ_TO_REP(MJ,IFSYM) = IREP
            DO MJ = 2, NORDER-2, 2
               IREP = IREP+1
               WRITE (REPA(IREP),1001) MJ,PARITY
               MJ_TO_REP(MJ,IFSYM) = IREP
               IREP = IREP+1
               WRITE (REPA(IREP),1001) -MJ,PARITY
               IF (MJ.NE.0) MJ_TO_REP(-MJ,IFSYM) = IREP
            ENDDO
C           For MJ = 16 we will not distinguish between positive and 
C           negative MJ
            MJ = NORDER
            IREP = IREP + 1
            WRITE (REPA(IREP),1001) MJ,PARITY
            MJ_TO_REP( MJ,IFSYM) = IREP
            MJ_TO_REP(-MJ,IFSYM) = IREP
         ENDDO
C
C        The multiplication table is rather straigthforward since 
C        the product irrep is the one with the sum of Mj values.
C        The only complication is the finite range of irreps so
C        we may need to add or subtract something, e.g. :
C        MJ = 17/2 -> MJ = -15/2 in the C8h subgroup.
C        (formulated in another way MJ = +/- 8 transforms
C         according to Ag or Au in C8h)
C        
          DO MJ = -NORDER, NORDER, 1
             DO MI = -NORDER, NORDER, 1
                MIJ = MI + MJ 
C               Map Cinfh to C8h
                IF (MIJ.LT.-NORDER) MIJ = MIJ + 2*NORDER
                IF (MIJ.GT. NORDER) MIJ = MIJ - 2*NORDER
                IF (MIJ.EQ.-NORDER) MIJ = NORDER
C               The gg-part
                IREP = MJ_TO_REP(MI,1)
                JREP = MJ_TO_REP(MJ,1)
                IJREP = MJ_TO_REP(MIJ,1)
                MULTB(IREP,JREP) = IJREP
C               The ug-part
                IREP = MJ_TO_REP(MI,2)
                JREP = MJ_TO_REP(MJ,1)
                IJREP = MJ_TO_REP(MIJ,2)
                MULTB(IREP,JREP) = IJREP
C               The gu-part
                IREP = MJ_TO_REP(MI,1)
                JREP = MJ_TO_REP(MJ,2)
                IJREP = MJ_TO_REP(MIJ,2)
                MULTB(IREP,JREP) = IJREP
C               The uu-part
                IREP = MJ_TO_REP(MI,2)
                JREP = MJ_TO_REP(MJ,2)
                IJREP = MJ_TO_REP(MIJ,1)
                MULTB(IREP,JREP) = IJREP
             ENDDO
          ENDDO
      ENDIF
C
 1000 FORMAT (I4)
 1001 FORMAT (I3,A)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck gmultsf */
      SUBROUTINE GMULTSF (NSYMRPA,REPA,MULTB)
C***********************************************************************
C
C     We make the multiplication table as the direct product of the
C     spatial and spin multiplication tables. For the latter we have
C     two irreps a(lpha) and b(eta) that can couple to give mS=0 or
C     mS = +/- 1. The last class is treated as giving one-dimensional
C     irreps, i.e. we do not employ the full spin symmetry group SU(2)
C     but take the subgroup C4 for this to get an Abelian group that
C     can distinguish the different spin-projection up to |mS|=3/2.
C     An easy extension is to take the subgroup C8 so that one can go
C     up to |mS|=7/2, but this increases the amount of irreps by 2 and
C     has to be taken into account in the multiplication tables.
C
C     Luuk Visscher, Mar 18 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "pgroup.h"
C
      CHARACTER*4 REPA(64)
      CHARACTER*1 REPAS(8)
      DIMENSION MULTB(64,64),MULTBS(8,8)
      DATA MULTBS/7,5,8,6,1,3,4,2,
     &            5,8,6,7,2,4,1,3,
     &            8,6,7,5,3,1,2,4,
     &            6,7,5,8,4,2,3,1,
     &            1,2,3,4,5,6,7,8,
     &            3,4,1,2,6,5,8,7,
     &            4,1,2,3,7,8,6,5,
     &            2,3,4,1,8,7,5,6/
C
#include "dcbibt.h"
C
C     Start by making the spin-multiplication table in point group C_4
C
      REPAS(1) = 'a' ! mS = + 1/2
      REPAS(2) = 'b' ! mS = - 1/2
      REPAS(3) = '3' ! mS = + 3/2
      REPAS(4) = '3' ! mS = - 3/2
      REPAS(5) = '0' ! mS = 0
      REPAS(6) = '4' ! mS = +/- 4/2 = +/- 2
      REPAS(7) = '2' ! mS = + 2/2 = + 1
      REPAS(8) = '2' ! mS = - 2/2 = - 1
C
C     The multiplication table of the boson irreps can be
C     obtained by bit operations.
C     Make the direct product table
C 
      NSYMRPA = NBSYM * 4
C
      JSR = 0
      DO JS = 1, 8
         DO JREP = 0, NBSYM-1
            JSR = JSR + 1
            ISR = 0
            DO IS = 1, 8
               DO IREP = 0, NBSYM-1
                  ISR = ISR + 1
                  IJS = MULTBS(IS,JS)
                  IJREP = IBTXOR(IREP,JREP)
                  IJSR = (IJS-1)*NBSYM + IJREP + 1
                  MULTB(ISR,JSR) = IJSR
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     Give the representations a name.
C
      ISR = 0
      DO IS = 1, 8
         DO IREP = 0, NBSYM-1
            ISR = ISR + 1
            REPA(ISR) = REP(IREP)//REPAS(IS)
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck moltodir */
      SUBROUTINE MOLTODIR(FMO,FMOM,NSTR,NSTRT,FORWARD,IPRINT)
C***********************************************************************
C     
C Purpose: Transform DIRAC Fock matrix, FMO, into MOLFDIR format, FMOM
C          (if FORWARD=.true, otherwise transform back)     
C     
C     Luuk Visscher, Mar 26 2002
C         Miro Ilias, June 2007 - added IPRINT
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
C
      LOGICAL FORWARD
      DIMENSION FMO(NORBT,NORBT,NZ)
      DIMENSION FMOM(NSTRT*2,NSTRT*2,2)
      DIMENSION NSTR(2)

      CALL QENTER('MOLTODIR')
C
C
C     Transform from matrix in quaternion DIRAC form to complex MOLFDIR form 
C     ... or backwards if we want to go from RELCCSD back to DIRAC
C
      IF (FORWARD) THEN
         CALL DZERO(FMOM,NSTRT*NSTRT*8)
      ELSE
         CALL DZERO(FMO ,NORBT*NORBT*NZ)
      ENDIF
C
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
        IOFF3 = IOFF3 + NSTR(IFSYM)
        DO J = 1, NSTR(IFSYM)
          IF (FORWARD) THEN
             CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,1),1,
     &                              FMOM(1+IOFF2,J+IOFF2,1),1)
             CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,1),1,
     &                              FMOM(1+IOFF3,J+IOFF3,1),1)
          ELSE
             CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF2,1),1,
     &                              FMO (1+IOFF1,J+IOFF1,1),1)
          ENDIF
        ENDDO
        IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
          DO J = 1, NSTR(IFSYM)
            IF (FORWARD) THEN
               CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                                FMOM(1+IOFF2,J+IOFF2,2),1)
               CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,2),1,
     &                                FMOM(1+IOFF3,J+IOFF3,2),1)
               CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF3,2),1)
             ELSE
               CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF2,2),1,
     &                                FMO (1+IOFF1,J+IOFF1,2),1)
             ENDIF
          ENDDO
          IF (NZ.GE.4.AND.(.NOT.SPINFR)) THEN
            DO J = 1, NSTR(IFSYM)
              IF (FORWARD) THEN
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,3),1,
     &                                 FMOM(1+IOFF2,J+IOFF3,1),1)
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,3),1,
     &                                 FMOM(1+IOFF3,J+IOFF2,1),1)
                CALL DSCAL(NSTR(IFSYM),DM1,FMOM(1+IOFF3,J+IOFF2,1),1)
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,4),1,
     &                                 FMOM(1+IOFF2,J+IOFF3,2),1)
                CALL DCOPY(NSTR(IFSYM),FMO (1+IOFF1,J+IOFF1,4),1,
     &                                 FMOM(1+IOFF3,J+IOFF2,2),1)
               ELSE
                CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF3,1),1,
     &                                 FMO (1+IOFF1,J+IOFF1,3),1)
                CALL DCOPY(NSTR(IFSYM),FMOM(1+IOFF2,J+IOFF3,2),1,
     &                                 FMO (1+IOFF1,J+IOFF1,4),1)
               ENDIF
            ENDDO
          ENDIF
        ENDIF
        IOFF1 = IOFF1 + NORB(IFSYM)
        IOFF2 = IOFF2 + 2*NSTR(IFSYM)
        IOFF3 = IOFF3 + NSTR(IFSYM)
      ENDDO
C
      IF (IPRINT.GE.7) THEN
        IF (FORWARD) THEN
         CALL HEADER(
     &   'MOLTODIR: transformation from DIRAC into MOLFDIR format',-1)
        ELSE
         CALL HEADER(
     &   'MOLTODIR: transformation from MOLFDIR into DIRAC format',-1)
        ENDIF
C       ....  
        WRITE(LUPRI,'(/3X,A)')
     &  'MOLTODIR:  *** matrix in DIRAC format, FMO ***'
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &               IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/3X,A)')
     &  'MOLTODIR:  *** matrix in MOLFDIR format, FMOM ***'
C       CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,MINO(NZ,2),
        CALL PRQMAT(FMOM,NSTRT*2,NSTRT*2,NSTRT*2,NSTRT*2,2,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
    
      CALL QEXIT('MOLTODIR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck moltodir */
      SUBROUTINE RELTODIR(XMO,XMOM,NFVO,IRCW,NREP,NV,NO)
C***********************************************************************
C     
C     Joost van Stralen, april 2002
C
C     Input:  XMOM - symmetry packed matrix in RELCCSD format (using supersymmetry for linear molecules)
C     Output: XMO  - The Lagrangian in Dirac format, this matrix is not sym. packed
C
C     Auxilliaries (all input and unchanged)
C     NFVO - total size of packed VO blocks (is identical to size of T1 amplitude array)
C     IRCW - real or complex arithmetic used in RELCCSD (1:real, 2:complex)
C     NREP - number of irreps used in RELCCSD
C     NV   - for each irrep the number of virtual spinors
C     NO   - for each irrep the number of occupied spinors
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
#include "dcbmp2.h"
C
      INTEGER NREP,NFVO,IRCW 
      DIMENSION XMO(NORBT,NORBT,NZ)
      DIMENSION XMOM(NFVO,IRCW)
      DIMENSION NV(NREP),NO(NREP)

      integer,allocatable :: idxu2g(:,:),idyu2g(:,:) 

      if (spinfr) NZ = 1

      allocate(idxu2g(norbt,nrep),idyu2g(norbt,nrep))

      idxu2g(1:norbt,1:nrep) = 0
      idyu2g(1:norbt,1:nrep) = 0

      call Super_Symmetry(nrep,nv,idxu2g,'particle') 

      call Super_Symmetry(nrep,no,idyu2g,'hole') 

      IF (nz.NE.4) THEN ! The most common case for small molecules "real" or complex" double groups

         if(nfsym == 1)then ! no inversion symmetry
           norbc = 1
           if (spinfr) then
              do j = 1, nbsym
                CALL MATGAT_symm_relcc(xmom(norbc,1),nv(j),no(j),
     &                                 xmo,norbt,norbt,idxu2g(1,j),
     &                                 idyu2g(1,j),nz)
                 norbc = norbc + no(j)*nv(j)
              end do
           else
              do j = 1, nrep,2
                if (nv(j)>0 .and. no(j)>0) 
     &          CALL MATGAT_symm_relcc(xmom(norbc,1),nv(j),no(j),
     &                                 xmo,norbt,norbt,idxu2g(1,j),
     &                                 idyu2g(1,j),nz)
              norbc = norbc + 2*nv(j)*no(j)
              end do
           endif
         else if(nfsym == 2)then ! inversion symmetry
!           LV: this code is incomplete, should be checked for spinfree with inversion
!           LV: non-spinfree code is identical to non-inversion code, can be simplified further
            norbc = 1
            do j = 1, nrep,2
              if (nv(j)>0 .and. no(j)>0) 
     &        CALL MATGAT_symm_relcc(xmom(norbc,1),nv(j),no(j),
     &                               xmo,norbt,norbt,idxu2g(1,j),
     &                               idyu2g(1,j),nz)
              norbc = norbc + 2*nv(j)*no(j)
            end do
         end if

         deallocate(idxu2g,idyu2g)

         IOFF1 = 0
         IOFF2 = 0
         IOFF3 = 0
         IOFF4 = 0

! Complex groups unless when they are run in spinfree (so that the algebra could be chosen as real)
         IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN 

            DO IFSYM = 1, NFSYM
               IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
               IOFF2 = IOFF2 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
                IAI = 1
                DO I = 1, NO(IFSYM+IOFF4)
                    CALL DCOPY(NV(IFSYM+IOFF4),     XMOM(IAI+IOFF3,2),1,
     &                                      XMO(1+IOFF1,I+IOFF2,2),1)
                    IAI = IAI + NV(IFSYM+IOFF4)
                ENDDO
                IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) + NORB(IFSYM)
                IOFF2 = IOFF2 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) 
     &                                      + NORB(IFSYM)
                IOFF3 = IOFF3 + 2*NV(IFSYM)*NO(IFSYM)
                IOFF4 = IOFF4 + 1
            ENDDO
         ENDIF

      ELSE ! branching for NZ, next comes the code for NZ=4 (quaternion groups, at most inversion symmetry)

         IOFF1 = 0
         IOFF2 = 0
         IOFF3 = 0
         DO IFSYM = 1, NFSYM
           IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
           IOFF2 = IOFF2 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
           IAI = 1
           DO I = 1, NAOCC(IFSYM)
C    Re
              CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                                XMO(1+IOFF1,I+IOFF2,1),1)

C    i
              CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                                XMO(1+IOFF1,I+IOFF2,2),1)

              IAI = IAI + 2*NAVIR(IFSYM)
           ENDDO
           IAI = 2*NAVIR(IFSYM)*NAOCC(IFSYM) + 1
           DO I = 1, NAOCC(IFSYM)
C    j
                CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                                  XMO(1+IOFF1,I+IOFF2,3),1)


C    k
                CALL DCOPY(NAVIR(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                                  XMO(1+IOFF1,I+IOFF2,4),1)

             IAI = IAI + 2*NAVIR(IFSYM)
           ENDDO
           IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) + NORB(IFSYM)
           IOFF2 = IOFF2 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) + NORB(IFSYM)
           IOFF3 = IOFF3 + 4*NAVIR(IFSYM)*NAOCC(IFSYM)
         ENDDO
C
      ENDIF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE DOVTODIR(XMO,XMOM,NFVO,IRCW,NREP,NO,NV)
C***********************************************************************
C     
C     Joost van Stralen, april 2002
C
C     XMO - The Lagrangian in Dirac format, this matrix is not sym. packed
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbtr3.h"
#include "dcbmp2.h"
C
      INTEGER NREP,NFVO,IRCW 
      DIMENSION XMO(NORBT,NORBT,NZ)
      DIMENSION XMOM(NFVO,IRCW)
      DIMENSION NV(NREP),NO(NREP)

      integer,allocatable :: idxu2g(:,:),idyu2g(:,:) 

      allocate(idxu2g(norbt,nrep),idyu2g(norbt,nrep))

      idxu2g(1:norbt,1:nrep) = 0
      idyu2g(1:norbt,1:nrep) = 0

      call Super_Symmetry(nrep,no,idxu2g,'hole') 

      call Super_Symmetry(nrep,nv,idyu2g,'particle') 


      IF (nz.NE.4) THEN

         if(nfsym == 1)then
           norbc = 1
       if (spinfr) then

           do j = 1, nbsym
             if (nv(j)>0 .and. no(j)>0) 
     &       CALL MATGAT_symm_relcc(xmom(norbc,1),no(j),nv(j),
     &                              xmo,norbt,norbt,idxu2g(1,j),
     &                              idyu2g(1,j),nz)
           norbc = norbc + no(j)*nv(j)
           end do

       else

           do j = 1, nrep,2
             if (nv(j)>0 .and. no(j)>0) 
     &       CALL MATGAT_symm_relcc(xmom(norbc,1),no(j),nv(j),
     &                              xmo,norbt,norbt,idxu2g(1,j),
     &                              idyu2g(1,j),nz)
           norbc = norbc + 2*nv(j)*no(j)
           end do
       endif

         else if(nfsym == 2)then
           do k = 1, nfsym 
             if(k == 1)then
               norbc = 1
               do j = 1, nrep/2,2
                 if (nv(j)>0 .and. no(j)>0) 
     &           CALL MATGAT_symm_relcc(xmom(norbc,1),no(j),nv(j),
     &                                  xmo,norbt,norbt,idxu2g(1,j),
     &                                  idyu2g(1,j),nz)
                 norbc = norbc + 2*nv(j)*no(j)
               end do
             else if(k == 2)then
               do j = nrep/2+1,nrep,2
                 if (nv(j)>0 .and. no(j)>0) 
     &           CALL MATGAT_symm_relcc(xmom(norbc,1),no(j),nv(j),
     &                                  xmo,norbt,norbt,idxu2g(1,j),
     &                                  idyu2g(1,j),nz)
                 norbc = norbc + 2*nv(j)*no(j)
               end do
             end if
           end do
         end if

         deallocate(idxu2g,idyu2g)


      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      IOFF4 = 0


         IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN

      DO IFSYM = 1, NFSYM
          IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
          IOFF2 = IOFF2 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
          IAI = 1
          DO I = 1, NV(IFSYM+IOFF4)
              CALL DCOPY(NO(IFSYM+IOFF4),     XMOM(IAI+IOFF3,2),1,
     &                                XMO(1+IOFF2,I+IOFF1,2),1)
          IAI = IAI + NO(IFSYM+IOFF4)
          ENDDO
        IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) + NORB(IFSYM)
        IOFF2 = IOFF2 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) + NORB(IFSYM)
!        IOFF3 = IOFF3 + 2*NAVIR(IFSYM)*NAOCC(IFSYM)
        IOFF3 = IOFF3 + 2*NV(IFSYM)*NO(IFSYM)
        IOFF4 = IOFF4 + 1
      ENDDO
         ENDIF
C
      ELSE
C
      IOFF1 = 0
      IOFF2 = 0
      IOFF3 = 0
      DO IFSYM = 1, NFSYM
        IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
        IOFF2 = IOFF2 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
        IAI = 1
        DO I = 1, NAVIR(IFSYM)
C    Re
           CALL DCOPY(NAOCC(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                             XMO(1+IOFF2,I+IOFF1,1),1)
C    i
           CALL DCOPY(NAOCC(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                             XMO(1+IOFF2,I+IOFF1,2),1)
           IAI = IAI + 2*NAOCC(IFSYM)
        ENDDO
        IAI = 2*NAVIR(IFSYM)*NAOCC(IFSYM) + 1
        DO I = 1, NAVIR(IFSYM)
C    j
             CALL DCOPY(NAOCC(IFSYM),     XMOM(IAI+IOFF3,1),1,
     &                               XMO(1+IOFF2,I+IOFF1,3),1)
C    k
             CALL DCOPY(NAOCC(IFSYM),     XMOM(IAI+IOFF3,2),1,
     &                               XMO(1+IOFF2,I+IOFF1,4),1)
             IAI = IAI + 2*NAOCC(IFSYM)
        ENDDO
        IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) + NORB(IFSYM)
        IOFF2 = IOFF2 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) + NORB(IFSYM)
        IOFF3 = IOFF3 + 4*NAVIR(IFSYM)*NAOCC(IFSYM)
      ENDDO
C
      ENDIF
C
      ioff1 = 0
      ioff2 = 0
      RETURN
      END

!========================================================================================

C/* Deck dootodir */
      SUBROUTINE DOOTODIR(DMO,DMOM,NO,NFOO,NREP,IRCW)
C***********************************************************************
C     
C     Joost van Stralen, may 2002
C
C     Transform occ-occ part of 2nd order Density matrix from RELCCSD format 
C     to Dirac format
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbtr3.h"
C
      INTEGER NFOO,NREP,IRCW
      DIMENSION DMO(NORBT,NORBT,NZ)
      DIMENSION DMOM(NFOO,IRCW)
      DIMENSION NO(NREP)
      integer, allocatable :: idxu2g(:,:)
C
C
!!!!       CALL DZERO(DMO,NORBT*NORBT*NZ)

      allocate(idxu2g(norbt,nrep))
      idxu2g(1:norbt,1:nrep) = 0

      call Super_Symmetry(nrep,no,idxu2g,'hole') 

!      write(*,*)'index array',(idxu2g(1,j),j=1,nrep,2)

      IF(NZ.NE.4) THEN
         IOFF1 = 0
         IOFF2 = 0
         IOFF4 = 0
         if(nfsym == 1)then
           norbc = 1

        if (spinfr) then

           do j = 1, nbsym
             CALL MATGAT_symm_relcc(dmom(norbc,1),no(j),no(j),
     &                              dmo,norbt,norbt,idxu2g(1,j),
     &                              idxu2g(1,j),nz)
           norbc = norbc + no(j)**2
           end do

        else

           do j = 1, nrep,2
             if (no(j)>0) 
     &       CALL MATGAT_symm_relcc(dmom(norbc,1),no(j),no(j),
     &                              dmo,norbt,norbt,idxu2g(1,j),
     &                              idxu2g(1,j),nz)
           norbc = norbc + 2*no(j)**2
           end do
        endif

         else if(nfsym == 2)then
           do k = 1, nfsym 
             if(k == 1)then
               norbc = 1
               do j = 1, nrep/2,2
                 if (no(j)>0) 
     &           CALL MATGAT_symm_relcc(dmom(norbc,1),no(j),no(j),
     &                                  dmo,norbt,norbt,idxu2g(1,j),
     &                                  idxu2g(1,j),nz)
                 norbc = norbc + 2*no(j)**2
               end do
             else if(k == 2)then
               do j = nrep/2+1,nrep,2
                 if (no(j)>0) 
     &           CALL MATGAT_symm_relcc(dmom(norbc,1),no(j),no(j),
     &                                  dmo,norbt,norbt,idxu2g(1,j),
     &                                  idxu2g(1,j),nz)
                 norbc = norbc + 2*no(j)**2
               end do
             end if
           end do
         end if

         deallocate(idxu2g)

         IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
               IJ = 1
               DO I = 1, NO(IFSYM+IOFF4)
                  CALL DCOPY(NO(IFSYM+IOFF4),       DMOM(IJ+IOFF2,2),1,
     &                                        DMO(1+IOFF1,I+IOFF1,2),1)
                  IJ = IJ + NO(IFSYM+IOFF4)
               ENDDO
            IOFF1 = IOFF1 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) +NORB(IFSYM)
            IOFF4 = IOFF4 + 1
            IOFF2 = IOFF2 + 2*no(ifsym)*no(ifsym) 
         ENDDO
         ENDIF
      ELSE
         IOFF1 = 0
         IOFF2 = 0
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NDMOQC(2,IFSYM,1)
C     Re       
            IJ = 1
            DO I = 1, (NO(IFSYM)/2)
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,1),1)
C     i
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,2),1)
            IJ = IJ + NO(IFSYM)
            ENDDO
C     j
            IJ = (NO(IFSYM)/2)*NO(IFSYM) + 1
            DO I = 1, (NO(IFSYM)/2)
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,3),1)
C     k
               CALL DCOPY((NO(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,4),1)
            IJ = IJ + NO(IFSYM)
            ENDDO
            IOFF1 = IOFF1 - NPSH(IFSYM) - NDMOQC(2,IFSYM,1) +NORB(IFSYM)
            IOFF2 = IOFF2 + NO(IFSYM)*NO(IFSYM)
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck dvvtodir */
      SUBROUTINE DVVTODIR(DMO,DMOM,NV,NFVV,NREP,IRCW)
C*****************************************************************************
C     
C     Joost van Stralen, may 2002
C
C     Transform vir-vir part of 2nd order Density matrix from RELCCSD format 
C     to Dirac format
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0=0.D0,DM1=-1.D0,D1 = 1.0D0)
C
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbtr3.h"
C
      INTEGER IRCW,NREP,NFVV
      DIMENSION DMO(NORBT,NORBT,NZ)
      DIMENSION DMOM(NFVV,IRCW)
      DIMENSION NV(NREP)
      integer, allocatable :: idxu2g(:,:)

      allocate(idxu2g(norbt,nrep))

      idxu2g(1:norbt,1:nrep) = 0

      call Super_Symmetry(nrep,nv,idxu2g,'particle') 

      IF(NZ.NE.4) THEN
         IOFF1 = 0
         IOFF2 = 0
         IOFF4 = 0
         if(nfsym == 1)then
           norbc = 1

        if (spinfr) then

           do j = 1, nbsym
             if (nv(j)>0) 
     &       CALL MATGAT_symm_relcc(dmom(norbc,1),nv(j),nv(j),
     &                              dmo,norbt,norbt,idxu2g(1,j),
     &                              idxu2g(1,j),nz)
           norbc = norbc + nv(j)**2
           end do

        else

           do j = 1, nrep,2
             if (nv(j)>0) 
     &       CALL MATGAT_symm_relcc(dmom(norbc,1),nv(j),nv(j),
     &                              dmo,norbt,norbt,idxu2g(1,j),
     &                              idxu2g(1,j),nz)
           norbc = norbc + 2*nv(j)**2
           end do
      
        endif


         else if(nfsym == 2)then
           do k = 1, nfsym 
             if(k == 1)then
               norbc = 1
               do j = 1, nrep/2,2
                 if (nv(j)>0) 
     &           CALL MATGAT_symm_relcc(dmom(norbc,1),nv(j),nv(j),
     &                                  dmo,norbt,norbt,idxu2g(1,j),
     &                                  idxu2g(1,j),nz)
                 norbc = norbc + 2*nv(j)**2
               end do
             else if(k == 2)then
               do j = nrep/2+1,nrep,2
                 if (nv(j)>0) 
     &           CALL MATGAT_symm_relcc(dmom(norbc,1),nv(j),nv(j),
     &                                  dmo,norbt,norbt,idxu2g(1,j),
     &                                  idxu2g(1,j),nz)
                 norbc = norbc + 2*nv(j)**2
               end do
             end if
           end do
         end if

         deallocate(idxu2g)

         IF (NZ.GE.2.AND.(.NOT.SPINFR)) THEN
            DO IFSYM = 1, NFSYM
                 IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
                 IJ = 1
                 DO I = 1, NV(IFSYM+IOFF4)
                  CALL DCOPY(NV(IFSYM+IOFF4),       DMOM(IJ+IOFF2,2),1,
     &                                        DMO(1+IOFF1,I+IOFF1,2),1)
                    IJ = IJ + NV(IFSYM+IOFF4)
                 ENDDO
              IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) +NORB(IFSYM)
              IOFF2 = IOFF2 + 2*NV(IFSYM)*NV(IFSYM)
              IOFF4 = IOFF4 + 1
           ENDDO
         ENDIF
      ELSE
         IOFF1 = 0
         IOFF2 = 0
         DO IFSYM = 1, NFSYM
            IOFF1 = IOFF1 + NPSH(IFSYM) + NOCC(IFSYM)
C     Re       
            IJ = 1
            DO I = 1, (NV(IFSYM)/2)
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,1),1)
C     i
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,2),1)
            IJ = IJ + NV(IFSYM)
            ENDDO
C     j
            IJ = (NV(IFSYM)/2)*NV(IFSYM) + 1
            DO I = 1, (NV(IFSYM)/2)
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,1),1,
     &                                  DMO(1+IOFF1,I+IOFF1,3),1)
C     k
               CALL DCOPY((NV(IFSYM)/2),      DMOM(IJ+IOFF2,2),1,
     &                                  DMO(1+IOFF1,I+IOFF1,4),1)
            IJ = IJ + NV(IFSYM)
            ENDDO
            IOFF1 = IOFF1 - NPSH(IFSYM) - NOCC(IFSYM) +NORB(IFSYM)
            IOFF2 = IOFF2 + NV(IFSYM)*NV(IFSYM)
         ENDDO
      ENDIF
C
      RETURN
      END
      SUBROUTINE MATGAT_symm_relcc(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,IDY,
     &                             NZ)
!***********************************************************************
!
!     !> Construct BMAT from AMAT elements using IDX and ibeig arrays.
!
!     Input: amat
!
!     Output: bmat
!
!     Written by S. Knecht - July 2014
!
!***********************************************************************
      implicit none
      real*8, intent(in)  :: AMAT(NRA,NCA,NZ) 
      real*8, intent(out) :: BMAT(NRB,NCB,NZ)
      integer, intent(in) :: IDX(*),IDY(*)
      integer, intent(in) :: nra, nrb, nca, ncb, nz

      real(8), parameter  :: threshold = 1.0d-16
      integer             :: iz, j, k, jj, kk
C
      DO IZ = 1, NZ
         do 5 J = 1, NCA
            JJ = IDY(J) 
            do 10 K = 1, NRA
               KK = IDX(K)
               if(abs(AMAT(K,J,IZ)).lt.threshold)then
                 BMAT(KK,JJ,IZ) = 0.0d0
               else
                 BMAT(KK,JJ,IZ) = AMAT(K,J,IZ)
               end if
 10         continue
  5      continue
      END DO
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine Super_Symmetry(nrep,orb_dim,idxu2g,orbital_type)

!------------------------------------------------------------------
! Make index array to map RELCCSD MO ordering on DIRAC MO ordering
!------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbtr3.h"
   
      integer, intent(in) :: nrep
      integer, intent(in)    :: orb_dim(nrep) 
      character(*), intent(in)    :: orbital_type 
      integer, intent(inout)   :: idxu2g(norbt,nrep)

      integer, allocatable :: ibeig(:)
      integer              :: nash_x(1:nfsym)

!     !> read ibeig (info array on boson irreps) from DFCOEF
      allocate(ibeig(norbt))
      ibeig(1:norbt)         = 0
      nash_x(1:nfsym)        = 0

!     Start by counting the number of spinors in each fermion irrep
      i = 1
      do j = 1, nrep
         if ( nfsym==2 .and.  j > nrep/2) i = 2
         nash_x(i) = nash_x(i) + orb_dim(j)
      end do

      ! but we want nash (meaning number of active shells = active Kramers pairs) so we need to
      ! divide by 2. Note that this will NOT work for unrestricted calculations.
      nash_x = nash_x/2

      if (spinfr) then
        do ifsym = 1,nfsym
          select case(orbital_type)
           case('hole') 
            nash_x(ifsym) = NAOCC(ifsym)
           case('particle')
           nash_x(ifsym) = NAVIR(ifsym)
          end select 
        enddo
      endif


      if( linear .or. spinfr )then
        call opnfil(lucoef,'DFCOEF','OLD','RIBEIG ')
        call reacmo(lucoef,'DFCOEF',dummy,dummy,ibeig,dummy,8)
        close(lucoef,STATUS='KEEP')
      end if


      if(nfsym == 1)then
        indx_loop = 0
        do j = 1, nrep,2
          indx_loop = indx_loop + 1
          if(orb_dim(j) > 0)then
            norbc = 0
            is_symm = 0
            if(linear)then
              !> get mj-value
              is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
            else if(spinfr)then
!              is_symm = mod(indx_loop,2)+1 - 1
             if (indx_loop.le.nbsym) is_symm = indx_loop-1
            end if
            do i = 1, nash_x(1)

            select case (orbital_type)

            case('hole')

                  if(is_symm == ibeig(iorb(1)       +
     &                                npsh(1)       + 
     &                                NDMOQC(2,1,1) +
     &                                i))then
                    norbc = norbc + 1
                    if (spinfr) then
                    if (indx_loop.le.nbsym) 
     &      IDXU2G(norbc,indx_loop) = iorb(1)+npsh(1)+NDMOQC(2,1,1)+i
                    else
                    IDXU2G(norbc,j) = iorb(1)+npsh(1)+NDMOQC(2,1,1)+i
                    endif

                  end if
            case('particle')

              if(is_symm == ibeig(iorb(1)       +
     &                            npsh(1)       + 
     &                            nocc(1)       +
     &                            i))then
                norbc = norbc + 1

                if (spinfr) then
                if (indx_loop.le.nbsym) 
     &          IDXU2G(norbc,indx_loop) = iorb(1)+npsh(1)+nocc(1)+i
                else
                IDXU2G(norbc,j) = iorb(1)+npsh(1)+nocc(1)+i
                endif

              end if
             end select

            end do
          end if
        end do
      else if(nfsym == 2)then
        do k = 1, nfsym 
          indx_loop = 0
          if(k == 1)then
            do j = 1, nrep/2,2
              indx_loop = indx_loop + 1
              if(orb_dim(j) > 0)then
                norbc = 0
                is_symm = 0
                if(linear)then
                  !> get mj-value
                  is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
                else if(spinfr)then
!                 is_symm = mod(indx_loop,2)+1 - 1
                  if (indx_loop.le.nbsym) is_symm = indx_loop-1
                end if
                do i = 1, nash_x(k)

            select case (orbital_type)

            case('hole')

                  if(is_symm == ibeig(iorb(k)       +
     &                                npsh(k)       + 
     &                                NDMOQC(2,k,1) +
     &                                i))then
                    norbc = norbc + 1


                if (spinfr) then
                if (indx_loop.le.nbsym) 
     &  IDXU2G(norbc,indx_loop) = iorb(k)+npsh(k)+NDMOQC(2,k,1)+i
                else
                IDXU2G(norbc,j) = iorb(k)+npsh(k)+NDMOQC(2,k,1)+i
                endif

                  end if
            case('particle')

              if(is_symm == ibeig(iorb(k)       +
     &                            npsh(k)       + 
     &                            nocc(k)       +
     &                            i))then
                norbc = norbc + 1

                if (spinfr) then
                if (indx_loop.le.nbsym) 
     &  IDXU2G(norbc,indx_loop) = iorb(k)+npsh(k)+nocc(k)+i
                else
                IDXU2G(norbc,j) = iorb(k)+npsh(k)+nocc(k)+i
                endif

              end if
             end select

                end do
              end if
            end do
          else if(k == 2)then
            do j = nrep/2+1,nrep,2
              indx_loop = indx_loop + 1
              if(orb_dim(j) > 0)then
                norbc = 0
                is_symm = 0
                if(linear)then
                  !> get mj-value
                  is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
                else if(spinfr)then
!                  is_symm = mod(indx_loop,2)+1 - 1
                  if (indx_loop.le.nbsym) is_symm = indx_loop-1
                end if
                do i = 1, nash_x(k)

            select case (orbital_type)

            case('hole')

                  if(is_symm == ibeig(iorb(k)       +
     &                                npsh(k)       + 
     &                                NDMOQC(2,k,1) +
     &                                i))then
                    norbc = norbc + 1

                if (spinfr) then
                if (indx_loop.le.nbsym) 
     &  IDXU2G(norbc,indx_loop) = iorb(k)+npsh(k)+NDMOQC(2,k,1)+i
                else
                IDXU2G(norbc,j) = iorb(k)+npsh(k)+NDMOQC(2,k,1)+i
                endif

                  end if
            case('particle')

              if(is_symm == ibeig(iorb(k)       +
     &                            npsh(k)       + 
     &                            nocc(k)       +
     &                            i))then
                norbc = norbc + 1

                if (spinfr) then
                if (indx_loop.le.nbsym) 
     &          IDXU2G(norbc,indx_loop) = iorb(k)+npsh(k)+nocc(k)+i
                else
                IDXU2G(norbc,j) = iorb(k)+npsh(k)+nocc(k)+i
                endif

              end if
             end select

                end do
              end if
            end do
          end if
        end do
      end if
      deallocate(ibeig)

      end subroutine

C  /* Deck dmp2ao */
      SUBROUTINE DMP2AO(DMO,DMAT,CMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by J. van Stralen - Feb 2003
C     
C     PURPOSE: Transform MO MP2 densitry matrix to quaternionic
C              symmetry-packed AO matrix 
C              
C     INPUT:   DMO  - MO density matrix
C              CMO  - spinor (MO) coefficients
C              
C     OUTPUT:  DMAT - AO MP2 denstity matrix
C     
C     REMARKS:
C     
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C         
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),DMO(NORBT,NORBT,NZ),WORK(*)
C           
C     
      CALL DZERO(DMAT,N2BBASXQ)
C     
      DO IFSYM = 1,NFSYM
C     
C        Parameters for first index
C
          NS = NORB(IFSYM)
          IS = IORB(IFSYM) + 1
          JS = ICMOQ(IFSYM) + 1
C
C        Parameters for first index
C
          NI  = NORB(IFSYM) 
          II  = IORB(IFSYM) + 1
          JI  = ICMOQ(IFSYM) + 1
C         
          CALL QTRANS('MOAO','S',D0,
     &                NFBAS(IFSYM,0),NFBAS(IFSYM,0),NS,NI,
     &                DMAT(I2BASX(IFSYM,IFSYM)+1,1),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,0),
     &                DMO(IS,II,1),NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),
     &                CMO(JS),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(IFSYM,0),NORB(IFSYM),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,1)
C     
      ENDDO
c       CALL HEADER('DMP2AO: AO density matrix',-1)
c       CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
c    &              NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
c       CALL FLSHFO(LUPRI)
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck CDENSO */
      SUBROUTINE CDENSO(DMAT,QCO,NDMOQC,ICMOQC,NSPC,NSPC2,IPRINT)
C*****************************************************************************
C
C Get the core density matrix into DMAT, from selected MO's in QCO
C Modified version of CDENS to also include one cor eopen shell
C Johann Pototschnig, Fall 2019
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      real(8)            :: DMAT(*),QCO(*)
      integer            :: NSPC(2,0:2),NDMOQC(2,2),ICMOQC(2)
      integer            :: NSPC2(2,0:3)
      real(8), parameter :: d1 = 1.0D0,d0 = 0.0D0
      integer            :: ncspin(2,0:2), i, j, iopen
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
C
      CALL DZERO(DMAT,2*N2BBASXQ)
C     get closed shell denstiy matrix
      do i = 1,NFSYM
        do j=0,2
          ncspin(i,j)=NSPC(i,j)-NSPC2(i,j)
        enddo

        IF(ncspin(i,0).EQ.0) CYCLE
        CALL DENST1(DMAT(1+I2BASX(i,i)),NTBAS(0),NTBAS(0),NZ,d1,d0,
     &              QCO(ICMOQC(i)),NDMOQC(1,i),NDMOQC(2,i),
     &              1,ncspin(i,0),NFBAS(i,0))
      enddo
C     get open-shell density matrix
      iopen=1
C     only use the first (one) open-shell
      do i = 1,NFSYM

        IF(NSPC2(i,0).EQ.0) CYCLE
        CALL DENST1(DMAT(1+N2BBASXQ*iopen+I2BASX(i,i)),
     &              NTBAS(0),NTBAS(0),NZ,d1,d0,
     &              QCO(ICMOQC(i)),NDMOQC(1,i),NDMOQC(2,i),
     &              ncspin(i,0)+1,NSPC2(i,0),NFBAS(i,0))
      enddo 
C     reset values
      do i = 1,NFSYM
            NDMOQC(2,i) = NSPC(i,0)
      enddo 
C     
      IF (IPRINT.GE.7) THEN
        CALL TITLER('Output from CDENSO','*',103)
        DO 20 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 20
          WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &    '*** Fermion corep ',I,'/',NFSYM
          IF (IPRINT.GE.10) THEN
            write(LUPRI,*)
     &      'Core molecular orbitals, QCO:',NDMOQC(1,I),NDMOQC(2,I)
            CALL PRQMAT(QCO(ICMOQC(I)),NDMOQC(1,I),NDMOQC(2,I),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          write(LUPRI,*) 'Core density matrix, DMAT:'
          CALL PRQMAT(DMAT(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
          write(LUPRI,*) 'Open-shell core density matrix, DMAT:'
          CALL PRQMAT(DMAT(I2BASX(I,I)+1+N2BBASXQ*iopen),NFBAS(I,0),
     &            NFBAS(I,0),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   20   CONTINUE
      ENDIF
      RETURN
      END
