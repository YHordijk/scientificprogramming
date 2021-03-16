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
C  /* Deck prpxlr_mpg */
      SUBROUTINE PRPXLR_MPG(DMO,XMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE: Driver routine for the calculation of the V-O and P-O 
C              part of the D^2 matrix. It starts by setting up parameters
C              that the normal response program gets from the input file 
C              or default settings.This includes setting up an operator.
C
C              Also D_ax parts of the density 
C              matrix are calculated
C
C     INPUT:   DMO - second order density matrix 
C              XMO - the Lagrangian
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxlr.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbxrs.h"
#include "dummy.h"
      dimension work(*)
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      logical orbdef(3)
      dimension pfac(3)
      character plabel(3)*8
      character tmp(4)*4
      CHARACTER PNAME*16

      call qenter('PRPXLR_MPG')
CMI   .. add this
      KFRSAV = KFREE

      CALL TITLER('Solving Z-vector equation','*',123)
C
C     Initialize /XCBLR/
C     ==================
C
      CALL SETRSP
      IPRXLR    = 0            !general printlevel
      THCLR     = 1.d-7         !maximum norm of solutionerror vector
      THRCOM    = 1.0d-6         !threshold for compression
      RESXLR    = 1.d3
      CNVXLR(1) = dummy
      CNVXLR(2) = 1.d-2
      MAXRM     = 200           
      ITRXLR    = 100
      INTXLR    = intgen
      ITRILR(1) = 1
      ITRILR(2) = 1
      NLRAPT    = 0             !number of A operators
      NLRBPT    = 0             !number of B operators
      CALL IZERO(NLRAP,8)
      CALL IZERO(NLRBP,8)
      NBFREQ    = 1             !number of frequencies
      BFREQ(1)  = 0.D0
      ALLCMB    = .FALSE.       !generate all operator combinations
      XLR_SKIPEP= .FALSE.       !do include e-p orbital rotations
!     if number of positronic shells is zero skip e-p rotations. hjj+sk - aug 2010
      if(nonrel.or.x2c.or.bss.or.levyle.or.freepj.or.vextpj)
     &  XLR_SKIPEP = .TRUE.
      XLR_SKIPEE= .FALSE.       !do include e-e orbital rotations
c     qedrun    = .TRUE.        !not implemented
      XLR_LSFG(1) = .TRUE.        !Generate G-matrix in MO basis
      XLR_LSFG(2) = .TRUE.        !Generate G-matrix in MO basis
      XLRNRM    = .TRUE.        !normalize solution vector
      XLRDIH    = .TRUE.        !pre-conditioning of trial vectors
C                                
      ORBDEF(1) = .TRUE.
      ORBDEF(2) = .TRUE.
      ORBDEF(3) = .TRUE.
      TRIAB     = .FALSE.       !force triangularity
      XLR_XVCFIL    = 'NORXVC'
      XSTERN    = .FALSE.
C
C     Make sure common block dcbidx is set
C
      CALL SETDCBIDX()
C
C     set work array to zero
C     ======================
C
      CALL DZERO(WORK,LFREE)
C
C     Initialize VIRTUAL orbital strings
C     ==================================
C
C     Virtual input strings
C
      do I=1,nfsym
         call num2str(nocc(I),tmp(1))
         call num2str(npsh(I),tmp(2))
         call num2str(nocc(I)+1,tmp(3))
         call num2str(nesh(I),tmp(4))
         xlr_indstr(1,I)='1..'//tmp(1)
         xlr_indstr(2,I)='0'
CMI      if(levyle)then
         if(levyle.or.bss.or.x2c.or.nonrel) then
            xlr_indstr(3,I)=tmp(3)//'..'//tmp(4)
         else
            xlr_indstr(3,I)='-'//tmp(2)//'..-1,'//tmp(3)//'..'//tmp(4)
         endif

      enddo
C
C     Initialize the operators needed
C     ===============================
C
      indxpr=-1
      isyxpr=-1
      itrxpr=-1
C
      iptyp=1
      npcomp=1
C
C     ---Right side operator
C
      pfac(1)=-1.d0
      plabel(1)='OVERLAP'
      call xprind('Lagrangian rhs. ',iptyp,npcomp,pfac,plabel,
     &     indxpr,isyxpr,itrxpr,iprxlr)
      call op1ind('NLRBPT',ind1op,llrbpu,nlrbpt,indxpr,maxllr)
C
C     Allocate memory
C     ===============
C     
      nrspf=nlrbpt*nbfreq
      if(xlr_skipee)then
         nee=0
      else
         nee=+1
      endif
      if(xlr_skipep)then
         nep=0
      else
         nep=1
      endif
      nci=0
      n2lrf=nlrapt*nrspf*(1+nee+nep+nci)
      call memget('REAL',klrf,n2lrf,work,kfree,lfree)
      call dzero(work(klrf),n2lrf)
      call memget('REAL',kcnv,nrspf,work,kfree,lfree)
C
C
      CALL MEMCHK('prpxlr_mpg 1.00',WORK,1)
C
      call prpxl1_mpg(work(klrf),work(kcnv),DMO,XMO,
     &                work,kfree,lfree)

CMI   CALL MEMCHK('prpxlr_mpg 1.01',WORK,1)
C
C     Calculate D_ax
C     ==============
C
      IF(.NOT.XLR_SKIPEP) THEN
C
        CALL CALCDAX(DMO,XMO,WORK,KFREE,LFREE)
C
      ENDIF
C
CMI   ... it's good to release the memmory after ...
      CALL MEMREL('prpxlr_mpg',WORK,KFRSAV,KFRSAV,KFREE,LFREE)

      call qexit('PRPXLR_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpxl1_mpg */
      SUBROUTINE PRPXL1_MPG(ABLRF,RCNV,DMO,XMO,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     INPUT:   ABLRF  - something for the linear response function,
C                       not needed for the MP2 gradient
C              RCNV   - total norm of residual relative to norm of
C                       solution vector
C              XMO   -  the Lagrangian
C     
C     REMARKS: -
C     
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
C
#include "aovec.h"
#include "maxorb.h"
#include "infpar.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxlr.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
C
#include "blocks.h"
      parameter(gpthrs=1.d-9)
      dimension ablrf(nfreq,nlrapt,nlrbpt,*)
      dimension rcnv(nbfreq,nlrbpt)
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
C
      call qenter('PRPXL1_MPG')
C
      kfrsav=kfree
C
C     Initialize
C     ==========
C
      call setxlr_mpg(work,kfree,lfree)
C
C     Linear response function
C     ========================
C
C
C     Dimensions of reduced system
C     ----------------------------
C
      kxlr=kfree
      nevecr=nredm*nexsim
      call memget('INTE',KIBTYP,2*nredm,work,kfree,lfree)
      call memget('INTE',kibcvc,  nredm,work,kfree,lfree)
      call memget('INTE',kibevc,  nredm,work,kfree,lfree)
      call memget('INTE',kibpvc,  nredm,work,kfree,lfree)
      call memget('REAL',kevalr, nexsim,work,kfree,lfree)
      call memget('REAL',kevecr, nevecr,work,kfree,lfree)
C
C     Transfer frequencies (actually only one frequency needed 0)
C     -----------------------------------------------------------
C
      call dcopy(nfreq,bfreq,1,work(kevalr),1)
C
C     Loop over operators
C     -------------------
C
      do jopsy=1,nfsym
         nfc=nbsym/nfsym
         do is=1,nfc
            jsymop=jfsym(is,jopsy)
            nopb=nlrbp(jsymop)
            jopb=jlrbp(jsymop)
            do 30 iopb=1,nopb
               indbp=llbpsu(jopb+iopb)
               indprb=llrbpu(indbp)
               jtimop=iprptim(indprb)
C
C     Initialize orbital parameters
C     =============================
C
               call xrspar(xlr_indstr,xlr_skipee,xlr_skipep,
     &                     IPRXLR)
C
C     Property gradient
C     =================
C
               kgrad=kfree
               call memget('REAL',kgpci,nzconfq,work,kfree,lfree)
               call memget('REAL',kgpoe,nzxopeq,work,kfree,lfree)
               call memget('REAL',kgpop,nzxoppq,work,kfree,lfree)
               call pamgrd_mpg(indprb,work(kgpoe),
     &                         work(kgpop),jopsy,
     &                         get_orbital_rotation_indices_pp(),
     &                         get_orbital_rotation_indices_pn(),
     &                         nzconf,nzxope,
     &                         nzxopp,XMO,work,kfree,lfree,iprxlr)
C
C     Check norm of gradient
C     ======================
C
               gpcnrm=dnrm2(nzconfq,work(kgpci),1)
               gpenrm=dnrm2(nzxopeq,work(kgpoe),1)
               gppnrm=dnrm2(nzxoppq,work(kgpop),1)
               gptnrm=sqrt(gpcnrm**2+gpenrm**2+gppnrm**2)
c               write(lupri,'(a,f18.10)') 'Gnorm',gptnrm
               if(gptnrm.lt.gpthrs)then
                  goto 30
               endif
C
C     Solve response equation
C     =======================
C
               ncred=0
               nered=0
               npred=0
               nzred=ncred+nered+npred
               RSPLAB = PRPNAM(INDPRB)
               call xrsctl(work(kgpci),work(kgpoe),work(kgpop),
     &                     work(KIBTYP),work(kibcvc),work(kibevc),
     &                     work(kibpvc),rcnv(1,indbp),work(kevalr),
     &                     work(kevecr),work,kfree,lfree)
               call memrel('PRPGR1.XRSCTL',work,1,kgrad,kfree,lfree)
C
C     Evaluate result
C     ===============
C
               call getlrf_mpg(ablrf,work(KIBTYP),work(kibcvc),
     &                         work(kibevc),work(kibpvc),RCNV(1,INDBP),
     &                         work(kevalr),work(kevecr),
     &                         DMO,work,kfree,lfree)
 30         continue
         enddo
      enddo
      call memrel('PRPGR1.ERR',work,1,kxlr,kfree,lfree)
C
C     Final remarks
C     =============
C
      close(luxvc,status='KEEP')

      call memrel('PRPGR1.FINAL',work,1,kfrsav,kfree,lfree)
CMI   call memrel('PRPGR1.FINAL',work,kfrsav,kfrsav,kfree,lfree)
C
      call qexit('PRPXL1_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxlr_mpg */
      SUBROUTINE SETXLR_MPG(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     Initialize variables/Open necessary files etc.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"

#include "dgroup.h"
#include "dcbxlr.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbxrs.h"
#include "cbihr2.h"

      logical tobe
      dimension work(*)
      call qenter('SETXLR_MPG')
      kfrsav=kfree
C
C     Sort operators on symmetries
C     ============================
C
      call op1srt('LLBPSU',llrbpu,llbpsu,nlrbpt,nlrbp,jlrbp,
     &             work(kfree),lfree)
C
C     Transfer information
C     ====================
C
      lineq     = .true.
      lsvcfg(1) = xlr_lsfg(1)
      lsvcfg(2) = xlr_lsfg(2)
      tknorm    = xlrnrm
      diaghe    = xlrdih
      iprxrs    = iprxlr
      thcxrs    = thclr
      resfac    = resxlr
      maxitr    = itrxlr
      nfreq     = nbfreq
      nexsim    = nfreq
      nexstv    = nfreq
      nexcnv    = nfreq
      nredm     = maxrm
      n2redm    = maxrm*maxrm
      loffty    = 0
      static    = .true.
      cnvint(1) = cnvxlr(1)
      cnvint(2) = cnvxlr(2)
      itrint(1) = itrilr(1)
      itrint(2) = itrilr(2)
      intdef    = intxlr
      sternh    = xstern
C
C     Set printlevel in two-electron routines
C     =======================================
C
      iprtwo=0
C
C     initialize file for solution vector
C     ===================================
C
      luxvc=29
      call opnfil(luxvc,'PAMXVC','UNKNOWN','PRPGRA')
      WRITE(LUXVC) 'END_OF_THIS_FILE'
C
      call qexit('SETXLR_MPG')
C
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pamgrd_mpg */
      SUBROUTINE PAMGRD_MPG(INDXPR,GPOE,GPOP,JOPSY,JXOPE,
     &                      JXOPP,NCONF,NXOPE,NXOPP,XMO,WORK,KFREE,
     &                      LFREE,IPRINT)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE:
C
C     INPUT:   NXOPE  - # of e-e orbital rotations
C              NXOPP  - # of e-p orbital rotations
C
C     REMARKS: -
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbxpr.h"
      dimension gpoe(*)
      dimension gpop(*)
      dimension jxope(*)
      dimension jxopp(*)
      DIMENSION XMO(*)
      dimension work(*)
      
      call qenter('PAMGRD_MPG')
C
C     Consistency check
C     =================
C
      isym=iprpsym(indxpr)
      iopsy=jbtof(isym-1,1)
      if(iopsy.ne.jopsy)then
         write(lupri,'(a)') 'Incorrect fermion ircop of property ',
     &        prpnam(indxpr)
         call quit('*** ERROR(PAMGRD_MPG) Incorrect fermion ircop')
      endif
C
C     Orbital part
C     ============
C
      nxopt=nxope+nxopp
      if(nxopt.gt.0)then
         call grdorb_mpg(indxpr,gpoe,gpop,jxope,jxopp,nxope,nxopp,
     &                   XMO,work,kfree,lfree,iprint)
      endif

      call qexit('PAMGRD_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck grdorb_mpg */
      SUBROUTINE GRDORB_MPG(INDXPR,GPOE,GPOP,JXOPE,JXOPP,NXOPE,
     &                      NXOPP,XMO,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE: Construct property gradient from the Lagrangian
C
C              e[1]_Q*(ma)=L_ma
C     
C              m - secondary (positronic and virtual)
C              a - occupied
C     
C     INPUT:   GPOE   - Property gradient,V-O part (this is a vector)
C              GPOP   - Property gradient,V-P part (this is a vector)
C              JXOPE  - Orbital rotation index, V-O
C              JXOPP  - Orbital rotation index, P-O
C              NXOPE  - # of occupied - virtual rotations
C              NXOPP  - # of occupied - positronic rotations
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
      dimension gpoe(*)
      dimension gpop(*)
      dimension jxope(*)
      dimension jxopp(*)
      DIMENSION XMO(NORBT,NORBT,NZ)
      dimension work(*)
C
      call qenter('GRDORB_MPG')
C
      open(lu1int,status='OLD',form='UNFORMATTED',file='AOPROPER')
C
      close(lu1int,status='KEEP')
C
C     Construct gradient
C     ==================
C
C     + ---------------+
C     | electronic part|
C     + ---------------+
C
      if(nxope.gt.0)then
         call xopgat(gpoe,jxope,nxope,XMO,norbt,nz)
         call dscal(nz*nxope,D2,gpoe,1)
      endif
C
C     +-----------------+
C     | positronic part |
C     +-----------------+
C
      if(nxopp.gt.0)then
         call xopgat(gpop,jxopp,nxopp,XMO,norbt,nz)
         call dscal(nz*nxopp,D2,gpop,1)
      endif
C     ====================
C
C     Printsection
C     ============
C     
       IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from GRDORB_MPG',-1)
        IREP = IPRPSYM(INDXPR)-1
        IF(NXOPE.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-e) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOE,NXOPE,1,NXOPE,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
        IF(NXOPP.GT.0) THEN
          WRITE(LUPRI,'(A,A16)') 
     &    '* Gradient (e-p) of property ',PRPNAM(INDXPR)
          CALL PRQMAT(GPOP,NXOPP,1,NXOPP,1,NZ,IPQTOQ(1,IREP),
     &                LUPRI)
        ENDIF
       ENDIF
C
      call qexit('GRDORB_MPG')
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck num2str */
      SUBROUTINE NUM2STR(NUM,RESSTR)
C***********************************************************************
C     
C  PURPOSE: Convert entering integer number NUM (of range 0-9999)
C           into the 4-byte string RESSTR.
C     
C  Written by C.V. Larsen
C  Last modifications: Miro ILIAS/Strasbourg, March 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      character str(4)*1
      character resstr*4
C
      CALL QENTER('NUM2STR')
C
      IF (NUM.LT.0.OR.NUM.GT.9999) THEN
        write(lupri,*)
     &  'NUM2STR: Entering number ',NUM,' is not in range 0-9999 !'
        CALL QUIT('NUM2STR: Wrong entering value of integer number !!!')
      ENDIF

      ifac=1000
      istr=1
      inum=num
 10   ich=inum/ifac
      str(istr)=char(ich+48)
      istr=istr+1
      inum=mod(inum,ifac)
      ifac=ifac/10
      if(ifac.ge.1)goto 10
C
      resstr=str(1)//str(2)//str(3)//str(4)
C     write(lupri,'(2a)') 'Num2Str: ',resstr
C
      CALL QEXIT('NUM2STR')
C
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getlrf_mpg */
      SUBROUTINE GETLRF_MPG(ABLRF,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,EVALR,
     &                      EVECR,DMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by C.V. Larsen
C     Last revision J. van Stralen - June 2002
C
C     PURPOSE: Retrieves the result of the response equation and puts it
C              in the VO and PO block of the 2nd order density matrix
C
C     INPUT:   ABLRF  - linear response function,
C                       not needed for the MP2 gradient
C              IBTYP  - type of the trial vector
C
C                  trial vectors will be of three kinds:
C
C                  C-type: contains only configurational parameters
C                              (WE DON'T NEED THIS ONE NOW!)
C                  E-type: contains only e-e rotations
C                  P-type: contains only e-p rotations
C              IBCVC  - pointer from C-type trial vectors to IBTYP
C              IBEVC  - pointer from E-type trial vectors to IBTYP
C              IBPVC  - pointer from P-type trial vectors to IBTYP
C              RCNV   - total norm of residual relative to norm of
C                       solution vector
C              EVALR  - eigenvalues of the reduced response equation
C              EVECR  - eigenvectors of the reduced response equation
C              DMO    - the D^2 matrix
C
C     REMARKS: -
C
C***********************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbxrs.h"
#include "dcbxlr.h"
#include "dcbibn.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "../moltra/dcbtr3.h"
      PARAMETER(D0=0.0D0,D2=2.0D0)
      logical first
      integer a
      dimension ablrf(nfreq,nlrapt,nlrbpt,*)
      dimension IBTYP(2,*)
      dimension ibcvc(*)
      dimension ibevc(*)
      dimension ibpvc(*)
      dimension evalr(*)
      dimension evecr(nredm,*)
      DIMENSION DMO(NORBT,NORBT,NZ)
      dimension work(*)
C
      CALL QENTER('GETLRF_MPG')
      kfrsav=kfree
C
C     Evaluate solution vectors for all frequencies
C     =============================================
C
      nesim=0
      npsim=0
      ncsim=0
      ifac=1
      if(.not.static) ifac=2
C
C     Get coefficients
C     ================
C
      call memget('REAL',kcmo,ncmotq,work,kfree,lfree)
      lubuf=22
      call opnfil(lubuf,'DFCOEF','OLD','GETLRF')
      IF(SPINFR) THEN
        CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
        CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),DUM,WORK(KIBRP),
     &              TOTERG,11)
      ELSE
        KIBRP = KFREE
        CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),DUM,IDUM,TOTERG,3)
      ENDIF
      close(lubuf,status='KEEP')
C
C     Solution vector (e-e) part
C     ==========================
C
      if(nzxope.gt.0)then
         nesim=nfreq*ifac
         kee=kfree
         call memget('REAL',kgvec,      nzxopeq,work,kfree,lfree)
         call memget('REAL',kxvec,nzxopeq*nfreq,work,kfree,lfree)
         call memget('REAL',kbbuf,      nzxopeq,work,kfree,lfree)
         CALL MEMGET('INTE',KIVEC,        NESIM,WORK,KFREE,LFREE)
CMI      call getlr1(ablrf(1,1,1,1),jbendx,evecr,evalr,RCNV,
         call getlr1(ablrf,jbendx,evecr,evalr,RCNV,
     &               IBTYP,ibevc,
     &               work(kgvec),work(kxvec),work(kbbuf),nzxope,
     &               work(kcmo),work(kibrp),iorbcl,
     &               get_orbital_rotation_indices_pp(),
     &               WORK(KIVEC),IDUM,IDUM,IDUM,IDUM,IDUM,DUM,DUM,
     &               work,kfree,lfree)
C
C     +-----------------------------+
C     | Result can be isolated here |
C     +-----------------------------+
C
C
         K = 0
         DO IZ=1,NZ
            IO1 = 0
            IO2 = 0
            DO IFRP=1,NFSYM
               IO1 = IO1 + NPSH(IFRP)
               IO2 = IO2 + NPSH(IFRP) + NOCC(IFRP)
               DO I = 1, NOCC(IFRP)
                  II = I + IO1
                  CALL DCOPY(NSSH(IFRP),WORK(KXVEC+K),1,
     &                                  DMO(1+IO2,II,IZ),1)
                  K = K + NSSH(IFRP)
               ENDDO
               IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
               IO2 = IO2 - NPSH(IFRP) - NOCC(IFRP) + NORB(IFRP)
            ENDDO
         ENDDO
         CALL MEMREL('GETLRF.e-e',WORK,1,KEE,KFREE,LFREE)
      endif
C
         IO1 = 0
         IO2 = 0
         DO IFRP = 1, NFSYM
            IO1 = IO1 + NPSH(IFRP)
            IO2 = IO2 + NPSH(IFRP) + NOCC(IFRP)
            DO IZ=1,NZ
               DO I=1,NOCC(IFRP)
                  II = I + IO1
                  DO A=1,NSSH(IFRP)
                     IA = A + IO2
c                    write(lupri,'(2i3,f18.10)') IA,II,DMO(IA,II,IZ)
                  ENDDO
               ENDDO
            ENDDO
            IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
            IO2 = IO2 - NPSH(IFRP) - NOCC(IFRP) + NORB(IFRP)
         ENDDO
C
C     Solution vector (e-p) part
C     ==========================
C
      IF(NZXOPP.GT.0)THEN
         npsim=nfreq*ifac
         call memget('REAL',kxvec,nzxoppq*nfreq,work,kfree,lfree)
         call memget('REAL',kgvec,nzxoppq      ,work,kfree,lfree)
         call memget('REAL',kbbuf,nzxoppq      ,work,kfree,lfree)
         CALL MEMGET('INTE',KIVEC,NPSIM        ,WORK,KFREE,LFREE)
CMI      call getlr1(ablrf(1,1,1,1),jbpndx,evecr,evalr,RCNV,IBTYP,
         call getlr1(ablrf,jbpndx,evecr,evalr,RCNV,IBTYP,
     &               ibpvc,work(kgvec),work(kxvec),work(kbbuf),nzxopp,
     &               work(kcmo),work(kibrp),iorbcl,
     &               get_orbital_rotation_indices_pn(),
     &               WORK(KIVEC),IDUM,IDUM,IDUM,IDUM,IDUM,DUM,DUM,
     &               work,kfree,lfree)
C
C    +---------------------------+
C    |Result can be isolated here|
C    +---------------------------+
C
         K = 0
         DO IZ=1,NZ
            IO1 = 0
            IO2 = 0
            DO IFRP=1,NFSYM
               IO1 = IO1 + NPSH(IFRP)
               DO I = 1, NOCC(IFRP)
                  II = I + IO1
                  CALL DCOPY(NPSH(IFRP),WORK(KXVEC+K),1,
     &                                  DMO(1+IO2,II,IZ),1)
                  K = K + NPSH(IFRP)
               ENDDO
               IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
               IO2 = IO2 + NORB(IFRP)
            ENDDO
         ENDDO
      ENDIF
C
         IO1 = 0
         IO2 = 0
         DO IFRP = 1, NFSYM
            IO1 = IO1 + NPSH(IFRP)
            DO IZ=1,NZ
               DO I=1,NOCC(IFRP)
                  II = I + IO1
                  DO A=1,NPSH(IFRP)
                     IA = A + IO2
c                    write(lupri,'(2i3,f18.10)') IA,II,DMO(IA,II,IZ)
                  ENDDO
               ENDDO
            ENDDO
            IO1 = IO1 - NPSH(IFRP) + NORB(IFRP)
            IO2 = IO2 + NORB(IFRP)
         ENDDO
      CALL MEMREL('GETLRF.end',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL MEMCHK('einde rspmpg',WORK,1)
      CALL QEXIT('GETLRF_MPG')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck calcdax */
      SUBROUTINE CALCDAX(DMO,XMO,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by Joost van Stralen september 2003.
C
C     PURPOSE : Calculate D_ax
C
C     Input :   
C
C     Output :
C
C     Remarks: 
C              D_(ax)=[L_(ax)-<aj"||xv>D_(vj") + <av||xj">D_(vj")^*]/DIV
C
C              with
C
C              DIV = eig_x - eig_a
C            
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
      LOGICAL FILEX
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
C
      CALL QENTER('CALCDAX')
C
      KFRSAV=KFREE
C
C     Allocate memory for D_ax and D_vj" density matrix. Note that we
C     allocate for the complete norbt*norbt*nz size because
C     this size is needed by MPGLVC. The other parts will be put to 0
C     ===============================================================
C
      CALL MEMGET('REAL',KDAX,NORBT*NORBT*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVJ,NORBT*NORBT*NZ,WORK,KFREE,LFREE)
C
C     We also need the eigenvalues
C     ============================
C
      INQUIRE(FILE='DFCOEF',EXIST=FILEX)
      IF(.NOT.FILEX) GOTO 1000
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CALCDAX')
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),IDUM,TOTERG,6)
      CLOSE(LUCOEF,STATUS='KEEP')
      CALL MEMCHK('na reacmo',WORK,1)
C
      CALL CALC1_DAX(DMO,XMO,WORK(KDAX),WORK(KDVJ),WORK(KEIG),WORK,
     &               KFREE,LFREE)
      CALL MEMCHK('na calc1_dax',WORK,1)
C
      CALL MEMREL('CALCDAX',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('CALCDAX')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'CALCDAX: Coefficient file not found !'
      CALL QUIT('CALCDAX: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck calc1_dax */
      SUBROUTINE CALC1_DAX(DMO,XMO,DAX,DVJ,EIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by Joost van Stralen september 2003.
C
C     PURPOSE : Calculate D_ax
C     
C     Input :
C
C     Output :
C            
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2 = 2.0D0)
C
#include "dcborb.h"
#include "dgroup.h"
C     
      DIMENSION DMO(*), DAX(*), DVJ(*)
      DIMENSION XMO(*), EIG(*)
      DIMENSION WORK(*)
C
      CALL QENTER('CALC1_DAX')
C
C     put DAX and DVJ to zero
C     =======================
C
      CALL DZERO(DAX,NORBT*NORBT*NZ)
      CALL DZERO(DVJ,NORBT*NORBT*NZ)
C
C     Copy the right parts of DMO to DVJ
C     ==================================
C
C     first the x - j" part
C     ---------------------
C
      DO IZ = 1,NZ
        IOFF1 = (IZ-1)*NORBT*NORBT
        IOFF2 = IOFF1
        DO IFSYM = 1,NFSYM
          IOFF2 = IOFF2 + NORBT*NPSH(IFSYM)
          IOFF3 = IOFF2 + 1
          DO J = 1,NOCC(IFSYM)
            CALL DCOPY(NPSH(IFSYM),DMO(IOFF3),1,DVJ(IOFF3),1)
            do k=1,npsh(ifsym)
            enddo
            IOFF3 = IOFF3 + NORBT
          ENDDO ! j
          IOFF2 = IOFF1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
        ENDDO   ! ifsym
      ENDDO     ! iz
C
C     next the a" - j" part
C     ---------------------
C
      DO IZ = 1,NZ
        IOFF1 = (IZ-1)*NORBT*NORBT
        IOFF2 = IOFF1
        DO IFSYM = 1,NFSYM
          IOFF2 = IOFF2 + NORBT*NPSH(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM)
          IOFF3 = IOFF2 + 1
          NV = NAVIR(IFSYM) + NIVIR(IFSYM)
          DO J = 1,NOCC(IFSYM)
            CALL DCOPY(NV,DMO(IOFF3),1,DVJ(IOFF3),1)
            do k=1,nv
            enddo
            IOFF3 = IOFF3 + NORBT
          ENDDO ! j
          IOFF2 = IOFF1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
        ENDDO   ! ifsym
      ENDDO     ! iz
C
C
C     Calculate the 2-electron dependent part of D_ax
C     ===============================================
C
c     CALL MPGLVC(DVJ,DAX,WORK,KFREE,LFREE,2)
C
C     Add L_ax and divide by spinor energies
C     ======================================
C
C
      DO IZ = 1,NZ
        IOFD1 = (IZ-1)*NORBT*NORBT
        IOFD2 = IOFD1
        DO IFSYM = 1,NFSYM
          IOFD2 = IOFD2 + NPSH(IFSYM) + NOCC(IFSYM)
          IOFD3 = IOFD2
          IOFEGA1 = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM) 
          IOFEGX1 = IORB(IFSYM)
          DO IX =1,NPSH(IFSYM)
            IOFD4 = IOFD3
            IOFEGX2 = IOFEGX1 + IX
            DO IA = 1,NAVIR(IFSYM)
              IOFD4 = IOFD4 + 1
              IOFEGA2 = IOFEGA1 + IA
              DIV = (EIG(IOFEGA2) - EIG(IOFEGX2))
c             DMO(IOFD4) = (D2*XMO(IOFD4) - DAX(IOFD4))/DIV
              DMO(IOFD4) = (D2*XMO(IOFD4))/DIV
            ENDDO !ia
            IOFD3 = IOFD3 + NORBT
          ENDDO   !ix
          IOFD2 = IOFD1 + NORBT*NORB(IFSYM) + NORB(IFSYM) 
        ENDDO     !ifsym
      ENDDO       !iz
C
      CALL QEXIT('CALC1_DAX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck calcdji */
      SUBROUTINE CALCDJI(DMO,XMO,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Written by Joost van Stralen september 2003.
C
C     PURPOSE : Calculate D_ji' and D_b'a
C
C               D_ji' = L_i'j/(eig_j - eig_i')
C
C               D_b'a = L_ab'/(eig_b' - eig_a)
C
C     Input :  XMO - The Lagrangian (L)
C     
C     Output : DMO - D_ji' and D_b'a
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
      LOGICAL FILEX
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
C     
      CALL QENTER('CALCDJI')
C
      KFRSAV=KFREE
C
      IF(NIOCCT.GT.0 .OR. NIVIRT.GT.0) THEN      
C
C       Open the coefficient file and read the eigenvalues
C       ==================================================
C
        INQUIRE(FILE='DFCOEF',EXIST=FILEX)
        IF(.NOT.FILEX) GOTO 1000
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CALCDJI')
        CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)
        CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),IDUM,TOTERG,6)
        CLOSE(LUCOEF,STATUS='KEEP')
C
C       Calculate D_ji' and D_b'a
C       =========================
C
        CALL CALC1_DJI(DMO,XMO,WORK(KEIG),WORK,KFREE,LFREE)
C
C       Create also the Hermitian counterpart of D_ji' and D_b'a
C       ========================================================
C
        CALL QHMRST(DMO,NORBT,NZ,NORBT,NORBT)
C
        CALL MEMREL('CALCDJI',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      ENDIF
C
      CALL QEXIT('CALCDJI')
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'CALCDJI: Coefficient file not found !'
      CALL QUIT('CALCDJI: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck calc1_dji */
      SUBROUTINE CALC1_DJI(DMO,XMO,EIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Written by Joost van Stralen september 2003.
C
C     PURPOSE : Calculate D_ji' and D_b'a
C
C               D_ji' = (L_i'j)^*/(eig_j - eig_i')
C
C               D_b'a = (L_ab')^*/(eig_b' - eig_a)
C
C     Input :  XMO - The Lagrangian (L)
C
C     Output : DMO - D_ji' and D_b'a
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0, D2 = 2.0D0)
C     
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION DMO(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*)
      DIMENSION EIG(*)
C     
C     
      CALL QENTER('CALC1_DJI')
C
C     Generate D_ji'
C     ==============
C
      IF(NIOCCT.GT.0) THEN
C
        IOFD1 = (1-1)*NORBT*NORBT
        IOFX1 = (1-1)*NORBT*NORBT
        IOFD2 = IOFD1
        IOFX2 = IOFX1
        DO IFSYM=1, NFSYM
          IOFD2 = IOFD2 + NORBT*NPSH(IFSYM) + NPSH(IFSYM) + NIOCC(IFSYM)
          IOFD3 = IOFD2
          IOFX2 = IOFX2 + NORBT*(NPSH(IFSYM)+NIOCC(IFSYM)) + NPSH(IFSYM)
          IOFX3 = IOFX2
          IOFEGI1 = IORB(IFSYM) + NPSH(IFSYM)
          IOFEGJ1 = IORB(IFSYM) + NPSH(IFSYM) + NIOCC(IFSYM)
          DO I = 1, NIOCC(IFSYM)
            IOFD4 = IOFD3
            IOFX3 = IOFX3 + 1
            IOFX4 = IOFX3
            IOFEGI2= IOFEGI1 + I
            DO J = 1, NAOCC(IFSYM)
              IOFD4 = IOFD4 + 1
              IOFEGJ2 = IOFEGJ1 + J
              DIV = (EIG(IOFEGJ2) - EIG(IOFEGI2))
              DMO(IOFD4) = XMO(IOFX4)/DIV
              IOFX4 = IOFX4 + NORBT
            ENDDO ! j
            IOFD3 = IOFD3 + NORBT
          ENDDO   ! i'
          IOFD2 = IOFD1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
          IOFX2 = IOFX1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
        ENDDO     ! ifsym
C     
C     In case we have complex or quaternion groups
C     --------------------------------------------
      IF(NZ.GT.1) THEN
      DO IZ =2,NZ 
        IOFD1 = (IZ-1)*NORBT*NORBT
        IOFX1 = (IZ-1)*NORBT*NORBT
        IOFD2 = IOFD1
        IOFX2 = IOFX1
        DO IFSYM=1, NFSYM
          IOFD2 = IOFD2 + NORBT*NPSH(IFSYM) + NPSH(IFSYM) + NIOCC(IFSYM)
          IOFD3 = IOFD2
          IOFX2 = IOFX2 + NORBT*(NPSH(IFSYM)+NIOCC(IFSYM)) + NPSH(IFSYM)
          IOFX3 = IOFX2
          IOFEGI1 = IORB(IFSYM) + NPSH(IFSYM)
          IOFEGJ1 = IORB(IFSYM) + NPSH(IFSYM) + NIOCC(IFSYM)
          DO I = 1, NIOCC(IFSYM)
            IOFD4 = IOFD3
            IOFX3 = IOFX3 + 1
            IOFX4 = IOFX3
            IOFEGI2= IOFEGI1 + I
            DO J = 1, NAOCC(IFSYM)
              IOFD4 = IOFD4 + 1
              IOFEGJ2 = IOFEGJ1 + J
              DIV = (EIG(IOFEGJ2) - EIG(IOFEGI2))
              DMO(IOFD4) = -XMO(IOFX4)/DIV
              IOFX4 = IOFX4 + NORBT
            ENDDO ! j
            IOFD3 = IOFD3 + NORBT
          ENDDO   ! i'
          IOFD2 = IOFD1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
          IOFX2 = IOFX1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
        ENDDO     ! ifsym
      ENDDO       ! iz
      ENDIF
C
      ENDIF
C
C     Generate D_b'a
C     ==============
C
      IF(NIVIRT.GT.0) THEN
C
C       Re
C      ----
        IOFD1 = (1-1)*NORBT*NORBT
        IOFX1 = (1-1)*NORBT*NORBT
        IOFD2 = IOFD1
        IOFX2 = IOFX1
        DO IFSYM=1, NFSYM
          IOFD2 = IOFD2 + NORBT*(NPSH(IFSYM) + NOCC(IFSYM)) 
     &          + NPSH(IFSYM) + NOCC(IFSYM) + NAVIR(IFSYM)
          IOFD3 = IOFD2
          IOFX2 = IOFX2 + NORBT*(NPSH(IFSYM)+NOCC(IFSYM)+NAVIR(IFSYM)) 
     &          + NPSH(IFSYM) + NOCC(IFSYM)
          IOFX3 = IOFX2
          IOFEGA1 = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM)
          IOFEGB1 = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM) 
     &            + NAVIR(IFSYM)
          DO IA = 1, NAVIR(IFSYM)
            IOFD4 = IOFD3
            IOFX3 = IOFX3 + 1
            IOFX4 = IOFX3
            IOFEGA2= IOFEGA1 + IA
            DO IB = 1, NIVIR(IFSYM)
              IOFD4 = IOFD4 + 1
              IOFEGB2 = IOFEGB1 + IB
              DIV = (EIG(IOFEGB2) - EIG(IOFEGA2))
              DMO(IOFD4) = XMO(IOFX4)/DIV
              IOFX4 = IOFX4 + NORBT
            ENDDO ! ib'
            IOFD3 = IOFD3 + NORBT
          ENDDO   ! ia
          IOFD2 = IOFD1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
          IOFX2 = IOFX1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
        ENDDO     ! ifsym
C
C     In case we have complex or quaternion groups
C     --------------------------------------------
      IF(NZ.GT.1) THEN
      DO IZ = 2,NZ
        IOFD1 = (IZ-1)*NORBT*NORBT
        IOFX1 = (IZ-1)*NORBT*NORBT
        IOFD2 = IOFD1
        IOFX2 = IOFX1
        DO IFSYM=1, NFSYM
          IOFD2 = IOFD2 + NORBT*(NPSH(IFSYM) + NOCC(IFSYM))
     &          + NPSH(IFSYM) + NOCC(IFSYM) + NAVIR(IFSYM)
          IOFD3 = IOFD2
          IOFX2 = IOFX2 + NORBT*(NPSH(IFSYM)+NOCC(IFSYM)+NAVIR(IFSYM))
     &          + NPSH(IFSYM) + NOCC(IFSYM)
          IOFX3 = IOFX2
          IOFEGA1 = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM)
          IOFEGB1 = IORB(IFSYM) + NPSH(IFSYM) + NOCC(IFSYM)
     &            + NAVIR(IFSYM)
          DO IA = 1, NAVIR(IFSYM)
            IOFD4 = IOFD3
            IOFX3 = IOFX3 + 1
            IOFX4 = IOFX3
            IOFEGA2= IOFEGA1 + IA
            DO IB = 1, NIVIR(IFSYM)
              IOFD4 = IOFD4 + 1
              IOFEGB2 = IOFEGB1 + IB
              DIV = (EIG(IOFEGB2) - EIG(IOFEGA2))
              DMO(IOFD4) = -XMO(IOFX4)/DIV
              IOFX4 = IOFX4 + NORBT
            ENDDO ! ib'
            IOFD3 = IOFD3 + NORBT
          ENDDO   ! ia
          IOFD2 = IOFD1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
          IOFX2 = IOFX1 + NORBT*NORB(IFSYM) + NORB(IFSYM)
        ENDDO     ! ifsym
      ENDDO       ! iz
      ENDIF
C
      ENDIF
C
      CALL QEXIT('CALC1_DJI')
      RETURN
      END
