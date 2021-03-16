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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck NMRINP */
      SUBROUTINE NMRINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for magnetic properties:
C --------------------------------------------
C  NMR shielding, magnetizabilities with/without London atomic orbitals 
C 
C     Written by T.Enevoldsen - Nov 1997
C
C     Added: MI & HJAaJ, 2002/2003
C
C*****************************************************************************
      use dirac_cfg
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
#include "implicit.h"
#include "priunit.h"
#include "orgcom.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "nuclei.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "dcbxlr.h"
#include "dgroup.h"
C

      PARAMETER (NTABLE = 16)
      LOGICAL :: SET, NEWDEF, RESET, LBIT, USECM
      LOGICAL :: INTNMR_SET
      LOGICAL :: oneexi, twoexi, fckexi
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &          PNAME*16, PLABEL(3)*8
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'xXXXXXX','.NOLOND','.SYMCON','.NOORTH',
     &            '.INTFLG','.NOTWO ','.NOONEI','.FULMAT',
     &            '.REGRAD','.LONDON','.USECM ','.MLLSUS',
     &            'xXXXXXX','.SEPTEP','.EXPPED','xXXXXXX'  /
      DATA SET/.FALSE./
C
      CALL QENTER('NMRINP')
#include "memint.h"
C
C     IF (RESET) GOTO 999

      NEWDEF = (WORD .EQ. '*NMR')
      IF (SET) THEN
         IF (NEWDEF)
     &      CALL QUIT('NMRINP: Only one "*NMR" input section allowed')
C        hjaaj: repeated input sections give infinite loop ...
         IF (RESET) SET = .FALSE.
         GOTO 999
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C     Initialize 
C     ===================
C
CMI   LONDON = SHIELD .OR. MSUSCP .OR. OPTROT
CMI ... by default we use conventional orbitals to ease running DIRAC tests
      LONDON = .FALSE.
      SYMCON = .FALSE.
      NOORTH = .FALSE.
      NOTWOL = .FALSE.
      NOONEI = .FALSE.
      USECM  = .FALSE.
      RDSUSLLMOD = .FALSE.
      SEPTEP = .FALSE.
      FULMAT = .FALSE.
      EXPPED = .FALSE.
      REGRAD = .FALSE.

! The INTNMR parameter here is read before INTXLR (in XLRINP) 
! assign top level default value from dcbgen.h instead of previous INTNMR = 7
      INTNMR = INTGEN
      INTNMR_SET = .FALSE.
C
C     Process input 
C     =========================
C
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in NMRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in NMRINP.')

!
    1          CONTINUE
               GO TO 100

    2          CONTINUE
C ... use conventional atomic orbitals instead of London atomic orbitals
               LONDON = .FALSE.
               GO TO 100
    3          CONTINUE
C ...  use symmetric connection for reorthonormalization terms
C ... default is to use the natural (non-symmetric) connection
               SYMCON = .TRUE.
               GO TO 100
    4          CONTINUE
C ... do not calculate reorthonormalization terms, {T(1),h(m)},
C     for the exp.value term h(Bm)
               NOORTH = .TRUE.
               GO TO 100
    5          CONTINUE
C INTFLG: specify what two-electron integrals to include
C         for the two-electron London contibution g(B)_mnpq
                 READ(LUCMD,*,IOSTAT=IOS) ILLINT,ISLINT,ISSINT
                 IF (IOS.NE.0) THEN
                  WRITE(LUPRI,*)
     &  'ILLINT,ISLINT,ISSINT:',ILLINT,ISLINT,ISSINT
                  CALL QUIT(
     & 'NMRINP: Error in reading ILLINT,ISLINT,ISSINT!')
                 ENDIF
                 INTNMR = ILLINT + 2*ISLINT + 4*ISSINT
                 INTNMR_SET = .TRUE.
               GO TO 100
    6          CONTINUE
C  ... do not calculate two-electron London contribution, g(B)_mnpq
               NOTWOL = .TRUE.
               GO TO 100
    7          CONTINUE
C ... do not calculate  reorthonormalization terms {T(1),H(0)}
C      for the property gradient H(B)
               NOONEI = .TRUE.
               GO TO 100
    8          CONTINUE
!              .FULMAT
               FULMAT = .true.
               GO TO 100
    9          CONTINUE
C ...REGRAD - restart London contributions to property gradient from file
               REGRAD = .TRUE.
!               READ(LUCMD,'(A6)') REGFIL
               GO TO 100
   10          CONTINUE
C&&&& .LONDON: use London atomic orbitals rather than conventional atomic orbitals
               LONDON = .TRUE.
               GO TO 100
   11          CONTINUE
C&&&& .USECM: use center of mass as a default gauge origin
               USECM = .TRUE.
               GO TO 100
   12          CONTINUE
C&&&& .MLLSUS: Activate the nonrelativistic modification of the relativistic RDSUSLL term
!miro: this is required to reproduce Dalton nonrelativistic results with the Levy-Leblond
! Hamiltonian and with the (RKB,c->infinite) Hamiltonian with the Sternheim's diamagnetic term
               RDSUSLLMOD = .TRUE.
               GO TO 100
   13          continue
               go to 100
   14          CONTINUE ! .SEPTEP
!              gosia: SEParate T matrix to electronic and positronic blocks
!              it does not affect shieldings, it is only for analysis.
               SEPTEP = .true.
               GO TO 100
   15          CONTINUE ! .EXPPED
!              gosia: export contributions to perturbed density to files
!              direct-LAO phase factor part write to pertden_direct_lao.FINAL
!              reorthonormalization part write to pertden_reorth_lao.FINAL
!              it is used further in embedding
               EXPPED = .true.
               GO TO 100
   16          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *         '" not recognized in NMRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in NMRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
C ... header: print out only if relevant (controlled by RESET)
      IF (.NOT.RESET) THEN
        CALL PRSYMB(LUPRI,'=',75,0)
        WRITE(LUPRI,'(3X,A)') 'Magnetic properties:'
        CALL PRSYMB(LUPRI,'=',75,0)
      ENDIF

      IF (USECM.AND.GAGORG_SET) THEN
       WRITE(LUPRI,'(/A/A)') 'WARNING:'//
     &    ' You have asked for center of mass to be the gauge origin!',
     &    'WARNING: ...your defined gauge origin is therefore ignored.'
       GAGORG_SET = .FALSE.
      ENDIF

      IF (GAGORG_SET .AND. .NOT.USECM) THEN
         if (maxrep > 0) then
            IF (.NOT.SYMCON.AND.LONDON) WRITE(LUPRI,'(/A/)')
     &     '   WARNING: User defined gauge '//
     &     'may give wrong results when using symmetry !'   
         end if
CMI      .... Don't forget to copy GAGORG into the DIPORG !!!
         DIPORG(:) = GAGORG(:)
      ELSE
       IF (LONDON.OR.USECM) THEN
C ... make the center of mass the gauge origin (default in dalton)
         CALL MEMGET('REAL',KGEOM  ,3*NATOMS ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KMASS  ,  NATOMS ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KNAT   ,  NATOMS ,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KNUMIS ,  NATOMS ,WORK,KFREE,LFREE)
C ... call the routine in abacus/hergp.F
         CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),
     &               WORK(KNUMIS),IPRDEF)
         CALL MEMREL('NMRINP.CMMASS',WORK,KWORK,KGEOM,KFREE,LFREE)
         CALL DCOPY(3,CMXYZ,1,GAGORG,1)
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
         WRITE(LUPRI,'(A,3F12.6)')
     &  ' * Center of mass used for the gauge & dipole '//
     &  'origins (a.u.) :',GAGORG(1:3)
       END IF
      END IF

      IF (LONDON) THEN
       IF (NOONEI) THEN
         WRITE(LUPRI,'(/A)')
     &' : No calculation of {T(1),H(0)} reorthonormalization'//
     &' terms for the H(B) property gradient !'
       ENDIF

#ifdef MOD_LAO_REARRANGED
       IF (NOORTH.AND. (SHIELD.OR.shielding_rearrange
     &     .OR.OPTROT)) THEN
#else
       IF (NOORTH.AND. (SHIELD.OR.OPTROT)) THEN
#endif
         WRITE(LUPRI,'(/,1X,A)')
     &': No calculation of {T(1),H(m_K)} reorthonormalization'//
     &' terms for the H(B m_K) expectation value terms !'
       ENDIF

       IF ((LEVYLE.OR.BSS.or.x2c).AND.(INTNMR.GT.1)) THEN
         INTNMR=1 ! Include only LL integrals...
         IF (INTNMR_SET) THEN
           WRITE(LUPRI,'(2X,A)') 'WARNING: For a '//
     &     'LevyLeblond/BSS/x2c Hamiltonian only LL integrals '//
     &     'for the dg(1,2)/dB terms are needed !'
         ENDIF
       ENDIF

       IF (SYMCON) THEN
         WRITE(LUPRI,'(1X,A)')
     & ': Symmetric (gauge independent) connection used for LAO-PDBS.'
       ELSE
         WRITE(LUPRI,'(1X,A)')
     & ': Natural (gauge dependent) connection used for LAO-PDBS.'
       ENDIF

       IF (NOTWOL) THEN
         WRITE(LUPRI,'(1X,A)')
     & ': No two-electron London contribution '//
     & '  to magnetic property gradient!'
       ELSE
         WRITE(LUPRI,'(1X,A)')
     & ': Include London magnet.prop.grad.contributions from '//
     & 'the following two-electron integrals:'
         IF(LBIT(INTNMR,1)) WRITE(LUPRI,'(6X,A)') '- LL-integrals'
         IF(LBIT(INTNMR,2)) WRITE(LUPRI,'(6X,A)') '- SL-integrals'
         IF(LBIT(INTNMR,3)) WRITE(LUPRI,'(6X,A)') '- SS-integrals'
         IF (.NOT.INTNMR_SET) THEN
            WRITE(LUPRI,"(2X,A)")
     &   "--> These are default values from the Hamiltonian input."
         ELSE
           WRITE(LUPRI,"(2X,A)")
     &   "--> These are accepted values from user's .INTFLG input."
         ENDIF
       ENDIF
      ELSE   !  (LONDON)
C   .... defined conventional atomic orbitals (AO), not London AO
       IF (SHIELD.OR.MSUSCP.OR.OPTROT) THEN
          WRITE(LUPRI,'(A/A/A,3F12.6)')
     &'  Using conventional atomic orbitals, not London orbitals !',
     &'  WARNING: Calculated magnetic properties are dependent'//
     &  ' on the choice of the gauge origin!',
     &'  NMRINP: The gauge & dipole origin :', GAGORG(1:3)
       END IF
      END IF  !  (LONDON)

C =========================================================================
CMI  ... for the magnetic susceptibilities take care the RDSUSLLMOD flag
C =========================================================================
      IF (MSUSCP) THEN
          IF (RDSUSLLMOD) THEN
            WRITE(LUPRI,'(1X,A)') 
     &'WARNING: Used the LS block of the RDSUSLL operator'//
     &' with factor -2 to achieve '//
     &'the Dalton nonrelativistic limit of the magnetizability !'
          ELSE
            WRITE(LUPRI,'(1X,A)') 
     &  'INFO: Used  full relativistic form of the RDSUSLL operator.'
          ENDIF
      END IF

! easy
!       if (.not. symcon) then
!         fulmat = .true.
!       end if 

       IF (FULMAT) THEN
          WRITE(LUPRI,'(A/A/A)')
     &'  FULMAT: when using the natural connection construct'//
     &' the full symmetric connection matrix'
       END IF
!     export perturbed density
      IF (EXPPED) THEN
         IF (.NOT.LONDON) THEN
            CALL QUIT( 
     &    'export of perturbed density only for runs with .LONDON')
         END IF
         WRITE(LUPRI,'(1X,A)')
     & 'perturbed LAO-density exported on file pertden_lao'
      END IF
      if (expped .and. parcal) then
        call quit('expped is tested in serial mode only; quitting...')
      end if
!     if fde and london, then always export perturbed density
! gosia: activate it if freeze&thaw with response...
!      if (dirac_cfg_fde_response .and. london) then
!         EXPPED = .true.
!      end if
      if (REGRAD) then
        IF (.NOT.LONDON) THEN
            CALL QUIT(
     &    'read london property gradient only for runs with .LONDON')
        END IF
         WRITE(LUPRI,'(1X,A,A6)')
     & 'London contribs to property gradient read from files'
        INQUIRE(FILE='ONELON',EXIST=oneexi)
        INQUIRE(FILE='TWOLON',EXIST=twoexi)
        INQUIRE(FILE='FCKLON',EXIST=fckexi)
        if (.not. oneexi) call quit(
     & 'no ONELON file found, remember to copy it by pam')
        if (.not. twoexi) call quit(
     & 'no TWOLON file found, remember to copy it by pam')
        if (.not. fckexi) call quit(
     & 'no FCKLON file found, remember to copy it by pam')
      end if

  999 CONTINUE
!gosia: here we are after reading xlr input, so epreorth has to be assigned here
!otherwise xlr_skipep is not known yet
      EPREORTH = .true.
      if (xlr_skipep) then
        EPREORTH = .false.
      end if

      CALL QEXIT('NMRINP')
      RETURN
      END

#ifdef MOD_LAO_REARRANGED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETCON_septep(ee,ep,occ,ORTO,WORK,KFREE,LFREE,IPRINT)
!gosia: later i will merge it with GETCON, now it is easier....
C***********************************************************************
      use london_helper
      use london_utils
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcborb.h"   
#include "dcbbas.h"   
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbxlr.h"
#include "dummy.h"
      PARAMETER(D2 = 2.0D0,D1 = 1.0D0,D0 = 0.0D0)
      PARAMETER(DM1 = -1.0D0,DM2 = - 2.0D0)
      DIMENSION ORTO(3,*),WORK(*), iqsym(4)
      LOGICAL AH,DOPCTRA_SAVE, ee, ep, occ
      real(8), allocatable:: conmat(:)
C
#include "ibtfun.h"
      CALL QENTER('GETCON')

      KFRSAV = KFREE

C  ...array for MO coefficients
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C boson irreps identification (read with MO coeff.)
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GETCON')
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident.
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,WORK(KBEIG),
     &            TOTERG,IOPT)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,IDUM,
     &            TOTERG,IOPT)
      ENDIF

C Open property file
      CALL OPNFIL(LU1INT,'AOPROPER','OLD','GETCON')

C  ... arrays for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KTMAT,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTHMAT,N2ORBXQ,WORK,KFREE,LFREE)

      IF (IPRINT.GE.5) CALL HEADER(
     &'GETCON: Creating transformed matrix; {T_c, h_mK}',-1)
C --------------------------------------------
C -- Calculate expect. value of {T(1),h(m)} --
C --------------------------------------------
      NSCOOR = 3*NUCDEP
      NOP = 0
      call alloc(conmat, n2orbxq)
      DO N = 1, NSCOOR  ! Run over coordinates m_K
         INDOP = IPSHIELD(N)
         IF (INDOP.GT.0) THEN ! Get h(m_K) in MO basis into CONMAT
           NOP = NOP + 1
           INDXPR = LLRAPU(INDOP) ! Index of h(m_K)
           ISYM  = IPRPSYM(INDXPR)
           IREP = ISYM - 1
           IOPSY = JBTOF(ISYM-1,1)
           CALL DZERO(CONMAT,N2ORBXQ)
           CALL PRPMAT(INDXPR,IOPSY,CONMAT,.TRUE.,WORK,WORK(KCMO),
     &                 WORK(KBEIG),ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
C --------------------------------------------------------------------
C ...  calculate the {T(1),h(m_K)}_ic reorthonormalization terms 
C    as expectation value
C --------------------------------------------------------------------
C .. do not perform the picture change transformation 
C    of connection matrixes !
           DOPCTRA_SAVE=DOPCT
           DOPCT=.FALSE.

           DO ICOMP  = 1, 3 ! Run over x,y,z components
C .... get the connection matrix (symmetric or natural) into WORK(KTMAT)
             INDXPR2 = IPCON(ICOMP) ! Get index of S(1)/W(1)
             ISYM2   = IPRPSYM(INDXPR2)
             ITIM2   = IPRPTIM(INDXPR2)
             IOPSY2  = JBTOF(ISYM2-1,1)
             IREP2   = ISYM2 - 1 ! 0-7
             CALL DZERO(WORK(KTMAT),N2ORBXQ)
#ifdef MOD_LAO_REARRANGED
!            gosia:
             if (lao_lr_rearrange) then
                CALL PRPMAT(INDXPR2,IOPSY2,WORK(KTMAT),.TRUE.,WORK,
     &               WORK(KCMO),WORK(KBEIG),ICMOQ,NORB,
     &               WORK,KFREE,LFREE,IPRINT)
!              gosia: if lao_lr_rearrange, then we asked for 'dS/dB'
!              integrals in def_d1hblond subroutine
             else
#endif
                CALL PRPMAT(INDXPR2,IOPSY2,WORK(KTMAT),.TRUE.,WORK,
     &               WORK(KCMO),WORK(KBEIG),ICMOQ,NORB,
     &               WORK,KFREE,LFREE,IPRINT)

#ifdef MOD_LAO_REARRANGED
             end if ! if (lao_lr_rearrange)
#endif
     
#ifdef MOD_LAO_REARRANGED
             if (shielding_rearrange) then
!              build (T_ij + T_ij^\dagger), as in getmdm subroutine:
               ih = -1
               do iz = 1, nz
                 iq = ipqtoq(iz, irep2)
                 iqsym(iz) = ihqmat(iq, ih)
               end do
               call qgetam(work(ktmat), norbt, nz, norbt, norbt, iqsym)
               call get_m_ij(work(ktmat))
             end if
             if (lao_lr_rearrange) then
               call get_m_ij(work(ktmat))
             end if

             if (ee) then
                call get_m_ee(work(ktmat))
                CALL DZERO(WORK(KTHMAT),N2ORBXQ)
                AH  = .TRUE. ! T_icomp is antihermitian ! 
                CALL TTRA(ISYM2,ISYM,WORK(KTMAT),
     &             CONMAT,WORK(KTHMAT),.true.,.false.,IPRINT)
             else if (ep) then
                call get_m_ep(work(ktmat))
                CALL DZERO(WORK(KTHMAT),N2ORBXQ)
                AH  = .TRUE. ! T_icomp is antihermitian ! 
                CALL TTRA(ISYM2,ISYM,WORK(KTMAT),
     &            CONMAT,WORK(KTHMAT),.true.,.false.,IPRINT)
             else if (occ) then
                call get_m_ij(work(ktmat))
                CALL DZERO(WORK(KTHMAT),N2ORBXQ)
                AH  = .TRUE. ! T_icomp is antihermitian ! 
                if (lao_lr_rearrange .or.shielding_rearrange) then
!                 then all contribution to {T,h} comes from occ-occ blocks
                  CALL TTRA(ISYM2,ISYM,WORK(KTMAT),
     &                CONMAT,WORK(KTHMAT),.true.,.true.,IPRINT)
                else 
                  CALL TTRA(ISYM2,ISYM,WORK(KTMAT),
     &                CONMAT,WORK(KTHMAT),.true.,.false.,IPRINT)
                end if
             end if
#endif

             IREP12 = IBTXOR(IREP2,IREP) ! resulting irrep of T_ic x m_K

      
C ... calculate the expectation value of the {T_ic,h(mK)} matrix
              EXPVAL = D0
              DO IFS = 1, NFSYM
              IF (NOCC(IFS).GT.0) THEN
C ... run only over occupied electronic states
              DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS)
               IADR1  = IORB(IFS)+I
               IADR2  = (IADR1-1)*NORBT+IADR1 ! diagonal element
               EXPVAL = EXPVAL+WORK(KTHMAT+IADR2-1)
              ENDDO
              ENDIF
              ENDDO

C ... store the reorthonormalization contribution into ORTO
              ITMP = IPTAX(ICOMP,2)
              ORTO(ITMP,NOP) = DM2*EXPVAL

           ENDDO ! Of x,y,z connection matrixes...
           DOPCT=DOPCTRA_SAVE ! restore the pict.change transf.
         ENDIF
      ENDDO
       call dealloc(conmat)
C
C close property file
C
      CLOSE(LU1INT,STATUS = 'KEEP')      

C ... release the memmory...
      CALL MEMREL('GETCON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETCON')
      RETURN
      END

#endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETCON */
      SUBROUTINE GETCON(W,CONMAT,DMAT,ORTO,WORK,
     &                  KFREE,LFREE,IPRINT)
C***********************************************************************
C
C   Get the reortonormalization contribution term,
C   {T(1),h(m_K)}, when using LAO. Calculated as expectation value.
C
C   MI/sept.2003: The .DOEPRN parameter decides whether to do
C   reorthonormalization over both p+e shells, or only over
C   the electronic shells. 
C
C***********************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
      use london_utils
#endif
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dgroup.h"
#include "dcborb.h"   
#include "dcbbas.h"   
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbxlr.h"
#include "dummy.h"
      PARAMETER(D2 = 2.0D0,D1 = 1.0D0,D0 = 0.0D0)
      PARAMETER(DM1 = -1.0D0,DM2 = - 2.0D0)
!               Array for {T,h(mK)}_ic product
      DIMENSION W(NORBT,NORBT,NZ,3),
!              Array for h(mK) property matrix
     &         CONMAT(NORBT,NORBT,NZ),
!              Density matrix
     &         DMAT(NORBT,NORBT,NZ),
     &         ORTO(3,*),WORK(*), iqsym(4)
      LOGICAL AH,DOPCTRA_SAVE
C
#include "ibtfun.h"
      CALL QENTER('GETCON')

      KFRSAV = KFREE

C  ...array for MO coefficients
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C boson irreps identification (read with MO coeff.)
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','GETCON')
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident.
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,WORK(KBEIG),
     &            TOTERG,IOPT)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,IDUM,
     &            TOTERG,IOPT)
      ENDIF

C Open property file
      CALL OPNFIL(LU1INT,'AOPROPER','OLD','GETCON')

C  ... arrays for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KTMAT,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTHMAT,N2ORBXQ,WORK,KFREE,LFREE)

      IF (IPRINT.GE.5) CALL HEADER(
     &'GETCON: Creating transformed matrix; {T_c, h_mK}',-1)
C --------------------------------------------
C -- Calculate expect. value of {T(1),h(m)} --
C --------------------------------------------
      NSCOOR = 3*NUCDEP
      NOP = 0
      DO N = 1, NSCOOR  ! Run over coordinates m_K
         INDOP = IPSHIELD(N)
         IF (INDOP.GT.0) THEN ! Get h(m_K) in MO basis into CONMAT
           NOP = NOP + 1
           INDXPR = LLRAPU(INDOP) ! Index of h(m_K)
           ISYM  = IPRPSYM(INDXPR)
           IREP = ISYM - 1
           IOPSY = JBTOF(ISYM-1,1)
           CALL DZERO(CONMAT,N2ORBXQ)
           CALL PRPMAT(INDXPR,IOPSY,CONMAT,.TRUE.,WORK,WORK(KCMO),
     &                 WORK(KBEIG),ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
           IF (IPRINT.GE.5) THEN
            WRITE(LUPRI,'(2A)')
     &      'GETCON: Get MO matrix of property:', PRPNAM(INDXPR)
            CALL PRQMAT(CONMAT,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
           ENDIF
C --------------------------------------------------------------------
C ...  calculate the {T(1),h(m_K)}_ic reorthonormalization terms 
C    as expectation value
C --------------------------------------------------------------------
C .. do not perform the picture change transformation 
C    of connection matrixes !
           DOPCTRA_SAVE=DOPCT
           DOPCT=.FALSE.

           DO ICOMP  = 1, 3 ! Run over x,y,z components
C .... get the connection matrix (symmetric or natural) into WORK(KTMAT)
             INDXPR2 = IPCON(ICOMP) ! Get index of S(1)/W(1)
             ISYM2   = IPRPSYM(INDXPR2)
             ITIM2   = IPRPTIM(INDXPR2)
             IOPSY2  = JBTOF(ISYM2-1,1)
             IREP2   = ISYM2 - 1 ! 0-7
             CALL DZERO(WORK(KTMAT),N2ORBXQ)
!            gosia:
#ifdef MOD_LAO_REARRANGED
             if (lao_lr_rearrange) then
                CALL PRPMAT(INDXPR2,IOPSY2,WORK(KTMAT),.TRUE.,WORK,
     &               WORK(KCMO),WORK(KBEIG),ICMOQ,NORB,
     &               WORK,KFREE,LFREE,IPRINT)
!              gosia: if lao_lr_rearrange, then we asked for 'dS/dB'
!              integrals in def_d1hblond subroutine
             else
#endif
                CALL PRPMAT(INDXPR2,IOPSY2,WORK(KTMAT),.TRUE.,WORK,
     &               WORK(KCMO),WORK(KBEIG),ICMOQ,NORB,
     &               WORK,KFREE,LFREE,IPRINT)
                IF(FULMAT) THEN
                  CALL FULMAT2('S',NORBT,NORBT,WORK(KTMAT))
                ENDIF
#ifdef MOD_LAO_REARRANGED
             end if ! if (lao_lr_rearrange)
#endif
!g             if (shielding_rearrange .or. lao_lr_rearrange) then
!g               call get_m_ij(work(ktmat))
!g             end if

#ifdef MOD_LAO_REARRANGED
             if (shielding_rearrange) then
!              build (T_ij + T_ij^\dagger), as in getmdm subroutine:
               ih = -1
               do iz = 1, nz
                 iq = ipqtoq(iz, irep2)
                 iqsym(iz) = ihqmat(iq, ih)
               end do
               call qgetam(work(ktmat), norbt, nz, norbt, norbt, iqsym)
               call get_m_ij(work(ktmat))
             end if
             if (lao_lr_rearrange) then
               call get_m_ij(work(ktmat))
             end if
#endif
             IREP12 = IBTXOR(IREP2,IREP) ! resulting irrep of T_ic x m_K
             CALL DZERO(WORK(KTHMAT),N2ORBXQ)
             AH  = .TRUE. ! T_icomp is antihermitian ! 
             !if (shielding_rearrange .or. lao_lr_rearrange) then
#ifdef MOD_LAO_REARRANGED
             if (shielding_rearrange.or.lao_lr_rearrange) then
             CALL TTRA(ISYM2,ISYM,WORK(KTMAT),
     &              CONMAT,WORK(KTHMAT),.true.,.true.,IPRINT)
             else
#endif
C ... do the reorthonormalization ...
             CALL TTRA(ISYM2,ISYM,WORK(KTMAT),
     &          CONMAT,WORK(KTHMAT),.true.,.false.,IPRINT)
#ifdef MOD_LAO_REARRANGED
             end if !if (shielding_rearrange .or. lao_lr_rearrange) then
#endif
             IF (IPRINT.GE.5) THEN
              WRITE(LUPRI,'(4A/)')
     &         'GETCON: Transformed term {}', PRPNAM(INDXPR),';',
     &          PRPNAM(IPCON(ICOMP))
              CALL PRQMAT(WORK(KTHMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
             ENDIF
C ... calculate the expectation value of the {T_ic,h(mK)} matrix
              EXPVAL = D0
              DO IFS = 1, NFSYM
              IF (NOCC(IFS).GT.0) THEN
C ... run only over occupied electronic states
              DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS)
               IADR1  = IORB(IFS)+I
               IADR2  = (IADR1-1)*NORBT+IADR1 ! diagonal element
               EXPVAL = EXPVAL+WORK(KTHMAT+IADR2-1)
              ENDDO
              ENDIF
              ENDDO

C ... store the reorthonormalization contribution into ORTO
              ITMP = IPTAX(ICOMP,2)
              ORTO(ITMP,NOP) = DM2*EXPVAL

              IF (IPRINT.GE.4) THEN
                 WRITE(LUPRI,'(4A,1X,E15.8,/)')
     &           'GETCON: Expect. value {Tc,h(m_K)} for ',
     &            PRPNAM(IPCON(ICOMP)),' and ',
     &            PRPNAM(INDXPR),
     &            ORTO(ITMP,NOP)
              ENDIF
           ENDDO ! Of x,y,z connection matrixes...
           DOPCT=DOPCTRA_SAVE ! restore the pict.change transf.
         ENDIF
      ENDDO
C
C close property file
C
      CLOSE(LU1INT,STATUS = 'KEEP')      

C ... release the memmory...
      CALL MEMREL('GETCON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETCON')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck TTRA */
      SUBROUTINE TTRA(ISYMT,ISYMH,TMAT,HMAT,THMAT,
     &                AH,half,IPRINT)
C*****************************************************************************
C
C   PURPOSE:
C ===========
C    Do one-index transformation of entering MO-matrixes, {T,h}.
C    This expression is used for reorthonormalization terms when
C    employing London atomic orbitals.
C
C Definition (Ruud et al., Chem.Phys. 195(1995) 157):
C -----------------------------------------------------
C  {T,h}= THMAT_{mn} = Sum_{o}(TMAT_{om}^*.HMAT_{on} + HMAT_{mo}.TMAT_{no}) 
C
C           in matrix notation:    T^{+}.H + H.T
C
C COMMENT: Helgaker,Jorgensen, JCP 95(1991) 2595 prefer (1 1^*) definition,
C          however, definition (1^* 1) is correct here because 
C          we have the same type
C          of one- and two-electron (1^* 1|2^* 2) matrix elements.
C          We have to obey the convention carefully since we are
C          working with general (quaternion,complex) quantities.
C
C
C On input:
C-----------
C  EPMIX = .true - run o-summation over both positronic and electronic shells,
C   else (EPMIX=.false.) run only over electronic shells.
C gosia: now epmix is replaced by epreorth (which is false only when xlr_skipep)
C
C  AH = .true. - take -1 for Hermitian conjugate, otherwise take +1
C
C  ISYMT,ISYMH ... boson symmetries on entering operators T,h
C
C  Note: The TMAT might also be a general matrix,
C        the entering h is usually anti-/hermitian
C
C  Written by Miro Ilias (based on TEC's routine TECTRA), august 2003, Odense.
C
C*****************************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, DM1 = -1.0D0)
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbnmr.h"
C
      DIMENSION  TMAT(NORBT,NORBT,NZ), 
     &           HMAT(NORBT,NORBT,NZ),
     &          THMAT(NORBT,NORBT,NZ)
      logical :: AH, half
C
#include "ibtfun.h"
      CALL QENTER('TTRA')

      IREPT = ISYMT-1
      IOPSYT = JBTOF(IREPT,1)

      IREPH = ISYMH-1
      IOPSYH = JBTOF(IREPH,1)

      IREPTH = IBTXOR(IREPT,IREPH) ! Symm.product of G(T)xG(H)
      IOPSYTH = JBTOF(IREPTH,1)

      if (AH) then
        PM = DM1
      else
        PM = D1
      end if

      IF(IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: Entering total H MO matrix of boson symmetry ',
     & ISYMH,'/',NBSYM
       CALL PRQMAT(HMAT,NORBT,NORBT,NORBT,
     &             NORBT,NZ,IPQTOQ(1,IREPH),LUPRI)
       WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: Entering total T MO matrix of boson symmetry ',
     & ISYMT,'/',NBSYM
       CALL PRQMAT(TMAT,NORBT,NORBT,NORBT,
     &             NORBT,NZ,IPQTOQ(1,IREPT),LUPRI)
      ENDIF

C ... do matrix multiplications - run over fermion irreps
C ... (I1,I2) - f.irreps of T, (I2,I3) - f.irreps of H

      DO I1 = 1, NFSYM
        IF (epreorth) THEN ! Matrix multiplication in full e+p MO space
#ifdef MOD_LAO_REARRANGED
          if (lao_lr_rearrange) then
!         called only from get_ftb
!         then we do: TH(I1,I3) = T(I1,I2).H(I2,I3) 
!         where T contains S^B matrix (only 'ai' elements)
!         and H has only 'ij' elements
            I2 = MOD(I1+IOPSYT,2) + 1  ! I1,I2 ... T
            I3 = MOD(I2+IOPSYH,2) + 1  ! I2,I3 ... H
            CALL QGEMM(NORB(I1), NORB(I3), NORB(I2), D1,
!                      -----------------------------------
     &                 'N', 'N', IPQTOQ(1,IREPT),
     &                 TMAT(IORB(I1)+1, IORB(I2)+1, 1),
     &                 NORBT,NORBT,NZ,
!                      -----------------------------------
     &                 'N', 'N', IPQTOQ(1,IREPH),
     &                 HMAT(IORB(I2)+1, IORB(I3)+1, 1),
     &                 NORBT,NORBT,NZ, 
!                      -----------------------------------
     &                 D1, IPQTOQ(1,IREPTH),
     &                 THMAT(IORB(I1)+1, IORB(I3)+1, 1),
     &                 NORBT,NORBT,NZ)
          else !if (lao_lr_rearrange) then
#endif
C  ...  do TH(I1,I3) := H(I1,I2).T(I2,I3)
          I2 = MOD(I1+IOPSYH,2) + 1  ! I1,I2 ... H
          I3 = MOD(I2+IOPSYT,2) + 1  ! I2,I3 ... T
          CALL QGEMM(NORB(I1), NORB(I3), NORB(I2), D1,
!                    ---------------------------------
     &               'N', 'N', IPQTOQ(1,IREPH),
     &               HMAT(IORB(I1)+1, IORB(I2)+1, 1),
     &               NORBT, NORBT, NZ,
!                    ---------------------------------
     &               'N', 'N', IPQTOQ(1,IREPT),
     &               TMAT(IORB(I2)+1, IORB(I3)+1, 1), 
     &               NORBT, NORBT, NZ, 
!                    ---------------------------------
     &               D0, IPQTOQ(1,IREPTH),
     &               THMAT(IORB(I1)+1, IORB(I3)+1, 1),
     &               NORBT,NORBT,NZ)
          IF(IPRINT.GE.5) THEN
            WRITE(LUPRI,'(/2X,A,I1,A,I1,F6.2)')
     &    'TTRA: 1.part transformed H.T over e-p :',
     &     IREPTH+1,'/',NBSYM,DSGN
            CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
          ENDIF
          if (.not.half) then
C        ... do TH(I1,I3) := TH(I1,I3) + dsgn*T^{H}(I1,I2).H(I2,I3)
          I2 = MOD(I1+IOPSYT,2) + 1  ! I1,I2 ... T
          I3 = MOD(I2+IOPSYH,2) + 1  ! I2,I3 ... H
          CALL QGEMM(NORB(I1), NORB(I3), NORB(I2), PM,  
!                    -----------------------------------
     &               'H', 'N', IPQTOQ(1,IREPT),
     &               TMAT(IORB(I2)+1, IORB(I1)+1, 1),
     &               NORBT,NORBT,NZ,
!                    -----------------------------------
     &               'N', 'N', IPQTOQ(1,IREPH),
     &               HMAT(IORB(I2)+1, IORB(I3)+1, 1),
     &               NORBT,NORBT,NZ, 
!                    -----------------------------------
     &               D1, IPQTOQ(1,IREPTH),
     &               THMAT(IORB(I1)+1, IORB(I3)+1, 1),
     &               NORBT,NORBT,NZ)
       end if !if (.not.half) then
#ifdef MOD_LAO_REARRANGED
         end if !if (lao_lr_rearrange) then
#endif
        ELSE  ! Transformation in the e-e space only
#ifdef MOD_LAO_REARRANGED
          if (lao_lr_rearrange) then
!         called only from get_ftb
!         then we do: TH(I1,I3) = T(I1,I2).H(I2,I3) 
!         where T contains S^B matrix (only 'ai' elements and only electronic 'a' indices)
!         and H has only 'ij' elements
            I2 = MOD(I1+IOPSYT,2) + 1  ! I1,I2 ... T
            I3 = MOD(I2+IOPSYH,2) + 1  ! I2,I3 ... H
            CALL QGEMM(NESH(I1), NESH(I3), NESH(I2), D1,
!                      -----------------------------------
     &                 'N', 'N', IPQTOQ(1,IREPT),
     &                 TMAT(IORB(I1)+1+NPSH(I1),IORB(I2)+1+NPSH(I2),1),
     &                 NORBT,NORBT,NZ,
!                      -----------------------------------
     &                 'N', 'N', IPQTOQ(1,IREPH),
     &                 HMAT(IORB(I2)+1+NPSH(I2),IORB(I3)+1+NPSH(I3),1),
     &                 NORBT,NORBT,NZ, 
!                      -----------------------------------
     &                 D1, IPQTOQ(1,IREPTH),
     &                 THMAT(IORB(I1)+1+NPSH(I1),IORB(I3)+1+NPSH(I3),1),
     &                 NORBT,NORBT,NZ)
          else
#endif
C  ...  do TH(I1,I3) := H(I1,I2).T(I2,I3)
          I2 = MOD(I1+IOPSYH,2) + 1  ! I1,I2 ... H
          I3 = MOD(I2+IOPSYT,2) + 1  ! I2,I3 ... T
          CALL QGEMM(NESH(I1),NESH(I3),NESH(I2),D1,
!                    -------------------------------------------------
     &               'N', 'N', IPQTOQ(1,IREPH),
     &               HMAT(IORB(I1)+NPSH(I1)+1,IORB(I2)+NPSH(I2)+1,1),
     &               NORBT, NORBT, NZ,
!                    -------------------------------------------------
     &               'N', 'N', IPQTOQ(1,IREPT),
     &               TMAT(IORB(I2)+NPSH(I2)+1,IORB(I3)+NPSH(I3)+1,1),
     &               NORBT, NORBT, NZ,
!                    -------------------------------------------------
     &               D0, IPQTOQ(1,IREPTH),
     &               THMAT(IORB(I1)+NPSH(I1)+1,IORB(I3)+NPSH(I3)+1,1),
     &               NORBT,NORBT,NZ)
          IF(IPRINT.GE.5) THEN
            WRITE(LUPRI,'(/2X,A,I1,A,I1,F6.2)')
     &    'TTRA: 1.part transformed H.T over e:',
     &     IREPTH+1,'/',NBSYM,DSGN
            CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
          ENDIF
          if (.not.half) then
C        ... do TH(I1,I3) := TH(I1,I3) + dsgn*T^{H}(I1,I2).H(I2,I3)
          I2 = MOD(I1+IOPSYT,2) + 1  ! I1,I2 ... T
          I3 = MOD(I2+IOPSYH,2) + 1  ! I2,I3 ... H
          CALL QGEMM(NESH(I1), NESH(I3), NESH(I2), PM, 
!                    ---------------------------------------------- --
     &               'H', 'N', IPQTOQ(1,IREPT),
     &               TMAT(IORB(I2)+NPSH(I2)+1,IORB(I1)+NPSH(I1)+1,1),
     &               NORBT, NORBT, NZ,
!                    -------------------------------------------------
     &               'N', 'N', IPQTOQ(1,IREPH),
     &               HMAT(IORB(I2)+NPSH(I2)+1,IORB(I3)+NPSH(I3)+1,1),
     &               NORBT, NORBT, NZ,
!                    -------------------------------------------------
     &               D1, IPQTOQ(1,IREPTH),
     &               THMAT(IORB(I1)+NPSH(I1)+1,IORB(I3)+NPSH(I3)+1,1),
     &               NORBT,NORBT,NZ)
          end if ! if (.not.half) then
#ifdef MOD_LAO_REARRANGED
          end if !if (lao_lr_rearrange) then
#endif
        ENDIF
      ENDDO
      IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/2X,A,I1,A,I1)')
     &'TTRA: TOTAL transformed {T,H} MO matrix over e of symm ',
     & IREPTH+1,'/',NBSYM
          CALL PRQMAT(THMAT,NORBT,NORBT,NORBT,
     &           NORBT,NZ,IPQTOQ(1,IREPTH),LUPRI)
      ENDIF

      CALL QEXIT('TTRA')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onelon */
      SUBROUTINE ONELON(CMO,IBEIG,WORK,KFREE,LFREE)
C***********************************************************************
C
C   PURPOSE:
C ==========
C  Get the total one-electron magnetic field property gradient dh_{D}(B)/dB
C  in MO basis when using London orbitals and write it to the DA file ONELON.
C
C  Expression:  RM1H1 + RM1H2 + RM1H3(sq) + RM1RN(sq) (x,y,z)
C
C  On output: CMO exports MO coefficients
C             IBEIG and boson irreps info...
C
C  Written by Miro Ilias, august 2003, Odense  
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"  
      PARAMETER(D0 = 0.0D0, D1= 1.0D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dgroup.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dummy.h"

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      DATA LUONELO / 72 /

      CALL QENTER('ONELON')

      KFRSAV = KFREE

      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','ONELON')
C ... distinguish between spin-free (boson symm.) mode
C     and the full relativistic (fermion symm.) mode
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident. 
       CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IBEIG,
     &            TOTERG,IOPT)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IDUM,
     &            TOTERG,IOPT)

      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')

C Open file for the London one-electron contributions, dh(B)/dB_c
      OPEN(LUONELO,FILE='ONELON',FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='UNKNOWN') ! Why 8* and not 3*?

      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')

C Allocate memory for storing components of dh(B)/dB_c

      CALL MEMGET('REAL',KDHB,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDHBIX,N2ORBXQ,WORK,KFREE,LFREE)

      DO ICOMP = 1, 3
       CALL DZERO(WORK(KDHB),N2ORBXQ)
       DO IX = 1, 4 ! run over all 4 integral types RM1...
         CALL DZERO(WORK(KDHBIX),N2ORBXQ)
         INDXPR=IPLONDON(IX,ICOMP)
         ITIM   = IPRPTIM(INDXPR)
         IREP   = IPRPSYM(INDXPR) - 1
         CALL PRPMAT(INDXPR,IOPSY,WORK(KDHBIX),.TRUE.,WORK,
     &               CMO,IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRPRP)
         IF (IPRPRP.GE.5) THEN
           WRITE(LUPRI,'(/A,A16)') 
     &     'ONELON: *MO-matrix of property ',
     &      PRPNAM(INDXPR)
            CALL PRQMAT(WORK(KDHBIX),NORBT,NORBT,
     &             NORBT,NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
         ENDIF

C   ...  accumulate all 1-el. prop. grad. contributions
         CALL DAXPY(N2ORBXQ,D1,WORK(KDHBIX),1,WORK(KDHB),1)
       ENDDO

C ... and write it the file:
       CALL WRTDAC(LUONELO,N2ORBXQ,WORK(KDHB),ICOMP)
       IF (IPRPRP.GE.5) THEN
         CALL HEADER(
     &'ONELON: The total LAO 1 el. mag.field gradient MO matrix:',-1)
         WRITE(LUPRI,'(2X,A,I2)') '** component:',ICOMP
         CALL PRQMAT(WORK(KDHB),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
       ENDIF
      ENDDO

      CLOSE(LUONELO,STATUS='KEEP')
      CLOSE(LU1INT,STATUS='KEEP')

C     .... release the whole block used in this subroutine
      CALL MEMREL('ONELON',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('ONELON')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getmdm */
      SUBROUTINE GETMDM(ISYMT,DMAT,CMO,IBEIG,TMAT,
     &                   WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C  PURPOSE:
C ----------
C  Generate SYMMETRIZED modified density matrixes which are used for
C the additional two-electron reorthonormalization term 
C when employing London orbitals.
C    
C     We define (based on Helgaker,Jorgensen,JCP 95(1991) 2595):
C   -----------------------------------------------------------------
C             DMAT_{ab} = Sum_{i,o}[ c_{ai} T_{oi} c_{bi}^* ]
C             a,b - AO-indices (L+S)
C             o   - secondary MO-index (electronic or both e+p shells)
C             i   - inactive MO-index (=occupied electronic MO's)
C
C     DMAT_{ab} is general matrix, so do symmetrize it 
C
C   On input: EPMIX=.true.=> o-summation runs over e+p shells, 
C             otherwise only over e-shells 
C             ISYMT - boson symmetry (1-NBSYM) of entering T operator matrix
C             CMO,TMAT - MO coeff. + transformation matrix T
C
C  The T is expected to be a general matrix (natural connection), 
C  or an antihermitian (symmetric connection) matrix.
C
C  Written by Miro Ilias (based on TEC's routine TECTR2),
C    august 2003, Odense.
C
C*****************************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,  D2 = 2.0D0)
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbnmr.h"
#include "dcbham.h"
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),IBEIG(*),
     &          TMAT(NORBT,NORBT,NZ),IQSYM(4),WORK(*)

      CALL QENTER('GETMDM')
      KFRSAV = KFREE
C
C     Initialization
C     ==============
C
      IREP = ISYMT-1 ! irre[ pf T(1)_icomp (0 till nirrep-1)
      IOPSY = JBTOF(IREP,1)
      CALL DZERO(DMAT,N2BBASXQ)
C
C     Print section
C
      IF(IPRINT.GE.8) THEN
        WRITE(LUPRI,'(/,A)')
     &  'GETMDM: Entering T(1) MO-connection matrix:'
        CALL PRQMAT(TMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &               IPQTOQ(1,IREP),LUPRI)
      ENDIF
      IF(IPRINT.GE.10) THEN
        WRITE(LUPRI,'(/A)')
     &  'GETMDM: Entering MO coefficients:'
        CALL PRQMAT(CMO,NTBAS(0),NORBT,NTBAS(0),NORBT,NZ,
     &               IPQTOQ(1,0),LUPRI)
      ENDIF
      write(*, *) 'epreorth in getmdm =', epreorth
C
C     Generate MDM = modified density matrix of f.symm I1,I2
C   ===========================================================
C
      DO  I1 = 1, NFSYM
          I2 = MOD(I1+IOPSY,2) + 1
C
C  DMAT_{ab} = Sum_{i,o}[ c_{ao} T_{oi} c_{bi}^* ] ... C T(oi) C^+
C
!gosia: fixme
!g        if (lao_lr_rearrange) then
!g            CALL QTRANS('MOAO','S',D0,NFBAS(I1,0),NFBAS(I2,0),
!g     &                  NOCC(I1),NOCC(I2),
!g     &       DMAT(I2BASX(I1,I2)+1,1),
!g     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
!g!          T(o,i)
!g     &     TMAT(IORB(I1)+NPSH(I1)+1,IORB(I2)+NPSH(I2)+1,1),
!g     &       NORBT,NORBT,NZ,IPQTOQ(1,IREP),
!g!          C(a,o)^+
!g     &     CMO(ICMOQ(I1)+(NPSH(I1)*NFBAS(I1,0))+1),
!g     &       NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
!g!          C(i,b)
!g     &     CMO(ICMOQ(I2)+(NPSH(I2)*NFBAS(I2,0))+1),
!g     &       NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
!g     &     WORK(KFREE),LFREE,IPRINT)
!g
!g        else  ! if not lao_lr_rearrange

        IF (EPREORTH) THEN ! Run 'o' over both (p+e) shell
          IF (NOCC(I2).GT.0.AND.NORB(I1).GT.0) THEN
            CALL QTRANS('MOAO','S',D0,NFBAS(I1,0),NFBAS(I2,0),
     &                  NORB(I1),NOCC(I2),
     &       DMAT(I2BASX(I1,I2)+1,1),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
!          T(o,i)
     &     TMAT(IORB(I1)+1,IORB(I2)+NPSH(I2)+1,1),
     &       NORBT,NORBT,NZ,IPQTOQ(1,IREP),
!          C(a,o)^+
     &     CMO(ICMOQ(I1)+1),
     &       NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
!          C(i,b)
     &     CMO(ICMOQ(I2)+(NPSH(I2)*NFBAS(I2,0))+1),
     &       NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)
          ENDIF

        ELSE ! Run o only over electronic states

        IF (NOCC(I2).GT.0.AND.NESH(I1).GT.0) THEN
          CALL QTRANS('MOAO','S',D0,NFBAS(I1,0),NFBAS(I2,0),
     &             NESH(I1),NOCC(I2),
     &             DMAT(I2BASX(I1,I2)+1,1),
     &             NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
     &             TMAT(IORB(I1)+NPSH(I1)+1,IORB(I2)+NPSH(I2)+1,1),
     &             NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &             CMO(ICMOQ(I1)+(NPSH(I1)*NFBAS(I1,0))+1),
     &             NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &             CMO(ICMOQ(I2)+(NPSH(I2)*NFBAS(I2,0))+1),
     &             NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &             WORK(KFREE),LFREE,IPRINT)
        ENDIF
      ENDIF
!g      end if ! if lao_lr_rearrange !gosia fixme
      ENDDO
      IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'GETMDM: NON-HERMITIAN modified density matrix',-1)
        IF (EPREORTH) THEN
         WRITE(LUPRI,'(2X,A)') "summation over both e+p shells..."
        ELSE
         WRITE(LUPRI,'(2X,A)') "summation over e shells..."
        ENDIF
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF

!     gosia: if lao_lr_rearrange we only do: DMAT = C*(SB)C(+)
#ifdef MOD_LAO_REARRANGED
      if (.not. lao_lr_rearrange) then
#endif
C
C     Take symmetric/antisymmetric combinations to form MDMAT
C    =========================================================
      IH = -1 ! T is TR antisymm.operator !
      DO IZ = 1,NZ
        IQ = IPQTOQ(IZ,IREP) ! Real quatern.number (1-4)
        IM = IHQMAT(IQ,IH) ! 1-S/2-A
        IQSYM(IZ) = IM
      ENDDO
C  ... in gp/gpsaue.F
      CALL QGETAM(DMAT,NTBAS(0),NZ,NTBAS(0),NTBAS(0),IQSYM)
#ifdef MOD_LAO_REARRANGED
      end if ! if (.not. lao_lr_rearrange) 
#endif
      IF(IPRINT.GE.5) THEN
        CALL HEADER(
     &'GETMDM: Total SYMMETRIZED modified density matrix',-1)
        IF (EPREORTH) THEN
         WRITE(LUPRI,'(2X,A)') "summation over both e+p shells..."
        ELSE
         WRITE(LUPRI,'(2X,A)') "summation over e shells..."
        ENDIF
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
      CALL QEXIT('GETMDM')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsklond */
      SUBROUTINE RSKLOND(FMAT,SKLTON,FACTOR,INDEX,IPOINT,
     &                   NDMAT,ISYMDM,IPRINT)
C*****************************************************************************
C
C   Purpose: Get two-electron Fock London matrixes in SA-AO basis
C
C   Written by Thomas Enevoldsen(TEC), 1998
C       and taken as a black-box by MI &HJAaJ, 2002 
C   Last modifications:  MI, March/2006 for pure twocomponent mode
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      INTEGER A, B, AP, BP, X
      DIMENSION FMAT(NBASIS,NBASIS,NDMAT,3), SKLTON(NBASIS,NBASIS),
     &          INDEX(KMAX,0:MAXREP,*), FACTOR(NBASIS,0:MAXREP),
     &          IPOINT(NBASIS,0:MAXREP),ISYMDM(NDMAT)
#include "nuclei.h"
#include "pincom.h"
#include "shells.h"
#include "symmet.h"
#include "ibtfun.h"

      CALL QENTER('RSKLOND')
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER(
     & 'RSKLOND: Input skeleton Fock matrices ',1)
         DO 100 X = 1, 3
            DO N = 1, NDMAT
             WRITE(LUPRI,'(2X,A,I1,A,I1)')
     & '=== Component(X=1,Y=2,Z=3):',X,' quaternion(1-4):',N 
             CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,
     &            NBASIS,NBASIS,1,LUPRI)
            END DO
  100    CONTINUE
      END IF
C
C     Form anti/symmetric skeleton Fock matrix
C     ========================================
C ... in this case Fock matrixes are imaginary=antisymmetric

      SIGN = -D1
      DO N = 1, NDMAT
         DO 200 X = 1, 3
            DO 210 I = 1, NBASIS
            DO 210 J = 1, I
              FMATIJ = DP5*(FMAT(I,J,N,X)+SIGN*FMAT(J,I,N,X))
              FMAT(I,J,N,X) =      FMATIJ
              FMAT(J,I,N,X) = SIGN*FMATIJ
 210        CONTINUE
 200     CONTINUE
         IF (IPRINT .GT. 5) THEN
            CALL HEADER(
     &  'RSKLOND: Skeleton Fock matrixes',1)
            DO 300 X = 1, 3
               WRITE(LUPRI,'(/,1X,A,I1,A,I1,A,I1)')
     &'N(1-NDMAT)/NDMAT=',N,'/',NDMAT,'  of component(x-1/y-2/z-3)=',X 
               CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,NBASIS,NBASIS
     &              ,1,LUPRI)
 300        CONTINUE
         END IF
         SIGN = D1
      END DO ! End of N

C===============================================================
C
C     Form AO Fock matrix from symmetric skeleton matrix
C     ==================================================
C
C===============================================================
      IF (MAXREP .GT. 0) THEN
C
C        First construct index and symmetry factor arrays
C
         IORB = 0
         DO 400 ISHELL = 1, KMAX
            DO 410 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)).EQ.0) THEN
               DO 420 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  INDEX(ISHELL,ISYMOP,ICOMP) = IORB
  420          CONTINUE
            ELSE
               ISYMF = IBTXOR(ISYMOP,IBTAND(ISYMOP,ISTBAO(ISHELL)))
               DO 430 ICOMP  = 1, KHKT(ISHELL)
                  INDEX(ISHELL,ISYMOP,ICOMP) = INDEX(ISHELL,ISYMF,ICOMP)
  430          CONTINUE
            END IF
  410       CONTINUE
  400    CONTINUE
C
         IORB = 0
         DO 500 ISHELL = 1, NLRGSH
            NHKTA = NHKT(ISHELL)
            DO 510 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 520 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 530 ISYM = 0, MAXREP
                     ISYMPR = IBTXOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) =
     &                        PT(IBTAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  530             CONTINUE
  520          CONTINUE
            END IF
  510       CONTINUE
  500    CONTINUE

C
C        The next paragraph concerns relativistic calculations only !
C        IPARity is invoked since the operation of inversion 
C        receives a minus sign when acting on small components....
C
         DO ISHELL = NLRGSH+1,KMAX
            NHKTA = NHKT(ISHELL)
            DO ISYMOP = 0, MAXREP
               IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
                  DO ICOMP  = 1, KHKT(ISHELL)
                     IORB = IORB + 1
                     DO ISYM = 0, MAXREP
                        ISYMPR = IBTXOR(ISYMOP,ISYM)
                        IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                        FACTOR(IORB,ISYM) = IPAR(ISYM)*
     &                       PT(IBTAND(ISYM,ISYMAO(NHKTA,ICOMP)))
                     END DO
                  END DO
               END IF
            END DO
         END DO
C
C        END of relativistic code .....
C
C
C        Construct full matrix from skeleton
C
         DO N = 1, NDMAT
            DO  X = 1, 3
               CALL DCOPY(NBASIS*NBASIS,FMAT(1,1,N,X),1,SKLTON,1)
               DO 600 ISYMOP = 1, MAXREP
                  FP = PT(IBTAND(ISYMOP,IBTXOR(ISYMDM(N),ISYMAX(X,2))))
                  DO 620 A = 1, NBASIS
                     AP = IPOINT(A,ISYMOP)
                     FA = FACTOR(A,ISYMOP)
                     DO 630 B = 1, NBASIS
                        BP = IPOINT(B,ISYMOP)
                        FB = FACTOR(B,ISYMOP)
                        FMAT(AP,BP,N,X)=FMAT(AP,BP,N,X)
     &                       +FA*FB*FP*SKLTON(A,B)
 630                 CONTINUE
 620              CONTINUE
 600           CONTINUE
            END DO
         END DO
C
         FAC = D1/dble(MAXREP + 1)
         NN = NDMAT*3*NBASIS*NBASIS
         CALL DSCAL(NN,FAC,FMAT,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("AO Fock matrices in RSKLOND",1)
            DO N = 1, NDMAT
               DO 700 X = 1, 3
                WRITE(LUPRI,'(/,1X,A,I1,A,I1,A,I1)')
     &'N(1-NDMAT)/NDMAT=',N,'/',NDMAT,'  of component(x-1/y-2/z-3)=',X 
                  CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,
     &                 NBASIS,NBASIS,1,LUPRI)
 700           CONTINUE
            END DO
         END IF
C
C      Check for relativistic calculation....
C
C     IF(KMAX.GT.NLRGSH) RETURN
CMI   IF(KMAX.GT.NLRGSH) GOTO 134
CMI   IF(KMAX.GE.NLRGSH) GOTO 134
      IF(.NOT.(KMAX.GE.NLRGSH)) THEN
CMI  ... seems that this body is never called....
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(3*(NBASIS**2),FMAT,1,SKLTON,1)
         CALL DZERO(FMAT,3*NBASIS*NBASIS)
         ISTRA = 1
         DO 800 IREPA = 0, MAXREP
            NORBA = NAOS(IREPA+1)
            DO 810 I = ISTRA,ISTRA + NORBA - 1
               IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
               NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
               NHKTA = NHKT(IA)
               KHKTA = KHKT(IA)
               MULA  = ISTBAO(IA)
               INDA  = KSTRT(IA) + NA - KHKTA
               IVARA = IBTXOR(IREPA,ISYMAO(NHKTA,NA))
               DO 820 ISYMA = 0, MAXOPR
               IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
                  FA    = PT(IBTAND(ISYMA,IVARA))
                  INDA  = INDA + KHKTA
                  ISTRB = 1
                  DO 830 IREPB = 0, MAXREP
                     NORBB = NAOS(IREPB+1)
                     DO 840 J = ISTRB,ISTRB + NORBB - 1
                        IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                        NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                        NHKTB = NHKT(IB)
                        KHKTB = KHKT(IB)
                        MULB  = ISTBAO(IB)
                        INDB  = KSTRT(IB) + NB - KHKTB
                        IVARB = IBTXOR(IREPB,ISYMAO(NHKTB,NB))
                        DO 850 ISYMB = 0, MAXOPR
                        IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                           INDB = INDB + KHKTB
                           FAB  = FA*PT(IBTAND(ISYMB,IVARB))
                           DO 860 X = 1, 3
                           IF (IBTXOR(IREPA,IREPB).EQ.ISYMAX(X,2)) THEN
                              FMAT(I,J,N,X) = FMAT(I,J,N,X)
     &                                    + FAB*SKLTON(INDA,INDB)
                           END IF
  860                      CONTINUE
                        END IF
  850                   CONTINUE
  840                CONTINUE
                     ISTRB = ISTRB + NORBB
  830             CONTINUE
               END IF
  820          CONTINUE
  810       CONTINUE
            ISTRA = ISTRA + NORBA
  800    CONTINUE

         IF (IPRINT .GT. 5) THEN
            CALL HEADER("SO Fock matrix in RSKLOND",1)
            DO 900 X = 1, 3
               CALL OUTPUT(FMAT(1,1,N,X),1,NBASIS,1,NBASIS,
     &                       NBASIS,NBASIS,1,LUPRI)
  900       CONTINUE
         END IF
      END IF

      ENDIF

C134  CONTINUE ! Needed this statement...

      CALL QEXIT('RSKLOND')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qmphase */
      SUBROUTINE QMPHASE(TYP,IM,JSYMOP,ONEINT,DOINT)
C*****************************************************************************
C
C     Quaternion phase insertion in lower triangular matrix
C       Phase: (e_K1*)(e_IM)(e_K2)
C
C     Theory: T.Saue &  H.J.Aa.Jensen
C     Written by T.Saue October 1994
C     Last revision: May 27 1998 - tec
C                    Added TYP 'S' for Square matrices, 'L' for lower triangle
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DM1 = -1.0D0)
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "symmet.h"
C
#include "dgroup.h"
#include "dcbbas.h"
      LOGICAL DOINT(2,2),DIAG,SAME
      DIMENSION ONEINT(*)
      CHARACTER TYP*1
C
#include "ibtfun.h"
C
      IL = 2
      IR = 1
      I1 = 0
      IF (TYP.EQ.'L') THEN
         DO 10 IRP1 = 0,NBSYM-1
            DO 20 IC1 = 1,2
               NB1 = NBBAS(IRP1,IC1)
               K1  = JQBAS(IRP1,IC1)
               I2 = 0
               DO 30 IRP2 = 0,IRP1
                  SAME = IRP1.EQ.IRP2
                  ICMX = 2
                  IF(SAME) ICMX = IC1
                  DO 40 IC2 = 1,ICMX
                     NB2 = NBBAS(IRP2,IC2)
                     K2  = JQBAS(IRP2,IC2)
                     DIAG = SAME.AND.(IC1.EQ.IC2)
                     IF(DOINT(IC1,IC2).AND.(IBTXOR(IRP1,IRP2).EQ.JSYMOP)
     &                    ) THEN
                        KFAC = IQPHASE(K1,K2,IM)*IQPH(K1,IL)*IQPH(K2,IR)
                        IF(KFAC.EQ.-1) THEN
                           DO 50 J1 = 1,NB1
                              JMX = NB2
                              IF(DIAG) JMX = J1
                              JOFF = (I1+J1)*(I1+J1-1)/2 + I2
                              DO 60 J2 = 1,JMX
                                 ONEINT(JOFF+J2) = -ONEINT(JOFF+J2)
 60                           CONTINUE
 50                        CONTINUE
                        ENDIF
                     ENDIF
                     I2 = I2 + NB2
 40               CONTINUE
 30            CONTINUE
               I1   = I1 + NB1
 20         CONTINUE
 10      CONTINUE
      ELSE IF (TYP.EQ.'S') THEN
         DO IRP1 = 0,NBSYM-1
            ISY1 = IRP1 + 1
            DO IC1 = 1,2
               NB1 = NBBAS(IRP1,IC1)
               K1  = JQBAS(IRP1,IC1)
               I2 = 0
               DO IRP2 = 0,NBSYM-1
                  ISY2 = IRP2 + 1
                  DO IC2 = 1,2
                     NB2 = NBBAS(IRP2,IC2)
                     K2  = JQBAS(IRP2,IC2)
                     IF(DOINT(IC1,IC2).AND.
     &                    (IBTXOR(IRP1,IRP2).EQ.JSYMOP)) THEN
                        KFAC = IQPHASE(K1,K2,IM)*IQPH(K1,IL)*IQPH(K2,IR)
                        IF(KFAC.EQ.-1) THEN
                           IOFF = I2COSX(ISY1,ISY2,IC1,IC2)
                           DO ICOL = 1,NB2
                              CALL DSCAL(NB1,DM1,ONEINT(IOFF+1),1)
                              IOFF = IOFF + NTBAS(0)
                           END DO
                        ENDIF
                     ENDIF
                     I2 = I2 + NB2
                  END DO
               END DO
               I1   = I1 + NB1
            END DO
         END DO
      ELSE
         CALL QUIT('QMPHASE: Unknown TYP !')
      ENDIF
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtcnvd */
      SUBROUTINE GTCNVD(DMAT,CMO,IBEIG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     (c) 1997/06/09 by J. Thyssen
C
C     Get contravariant denstity matrix without folding
C     Copied from GTCNV1
C     
C     On input: CMO - MO coeficcients
C               IBEIG - info on boson irreps      
C     
C     Last revisions: 1997/06/10 jth
C                     1998  TEC
C                     jan.2003 MI 
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WORK(*),DMAT(N2BBASX,4,NFMAT),CMO(*),IBEIG(*)
      LOGICAL FOLD
      PARAMETER ( D1 = 1.0 D00, D2 = 2.0 D00 )
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dummy.h"
C
#include "ibtfun.h"
      IQOFF(IZ,IMAT) = (IZ-1) * N2BBASX + (IMAT-1) * N2BBASXQ
C
#include "memint.h"
C
C
      CALL QENTER('GTCNVD')
C 
      CALL MEMGET('REAL',KTMP,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Contruct density matrix
C
      DO IMAT = 1, NFMAT ! Run over closed/open-shells
         IBIT = IBTSHL(1,IMAT-1)
C      ... routine in ../dirac/dirden.F
         CALL GENDEN(WORK(KTMP),CMO,IBIT,IPRINT)
C
C        Scale inactive density matrix with D2 (= occupation)
C
         IF (IMAT.EQ.1)
     &     CALL DSCAL(N2BBASXQ,D2,WORK(KTMP),1)
C
C        Transform to unsorted basis
C
         IREP = ISYMOP(IMAT)-1
         CALL BSTOBU(WORK(KTMP),NZ,WORK(KFREE),LFREE)
         IF(IPRINT.GE.5) THEN
           CALL HEADER('GTCNVD: '//
     &     'DMAT After BSTOBU to unsorted basis',-1)
           WRITE(LUPRI,'(A,I3)') '*** Matrix no.: ',IMAT
           CALL PRQMAT(WORK(KTMP),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         ENDIF
C
C        Insert half-phases when NZ < 4 (real and complex groups)
C     =============================================================
C
         IREP = ISYMOP(IMAT)-1
         IF(NZ.LT.4) THEN
            DO IZ = 1,NZ
              IQ = IPQTOQ(IZ,IREP)
              CALL Q2PHASE('D',IQ,1,WORK(KTMP+N2BBASX*(IZ-1)))
            ENDDO
            IF(IPRINT.GE.5) THEN
              CALL HEADER('GTCNVD: '//
     &  'Unsorted DMAT with half-phases inserted',-1)
              CALL PRQMAT(WORK(KTMP),NTBAS(0),NTBAS(0),NTBAS(0),
     &             NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
            ENDIF
         ENDIF
C
C        Transform density matrix from QO to AO
C        ======================================
C   .... analogous to sr. DIRSCF...
C
         CALL DZERO(DMAT(1,1,IMAT),4*N2BBASX)
         IREP = ISYMOP(IMAT)-1
         IPAR = JBTOF(IREP,1)
         DO IZ = 1,4
           IREPD = IRQMAT(IZ,IREP)
           IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
           IPQ   = IQTOPQ(IQ,IREP)
           IF(IPQ.GT.NZ) CALL QUIT('GTCNVD.DENSITY: IPQ.GT.NZ !')
           IF(IPQ.LE.0 ) CALL QUIT('GTCNVD.DENSITY: IPQ.LE.0 !')
           IADR = (IPQ-1) * N2BBASX

           CALL DTSOAO(WORK(KTMP+IADR),DMAT(1,IZ,IMAT),
     &                 NTBAS(0),IREPD,IPRINT)
         ENDDO
      ENDDO

      IF (NASHT.GT.0) THEN ! When open-shell are included
         CALL DAXPY(N2BBASX*4,D1,DMAT(1,1,2),1,DMAT(1,1,1),1)
      END IF
C
C        Final print section
C
      IF (IPRINT .GE. 4) THEN
       CALL HEADER('GTCNVD: Final quaternion density matrices.',-1)
       DO IMAT = 1,NFMAT
         WRITE(LUPRI,'(A,I1,A,I1,A)')
     &    '  **** ----- Matrix number(1-cs,2-os1,3-os2...): ',
     &    IMAT,'/',NFMAT,' ----- ****'
         DO IZ = 1, 4
          WRITE(LUPRI,'(/A,I1)')
     &    '  =====  quaternion number (1-4): ',IZ
          CALL OUTPUT(DMAT(1,IZ,IMAT),1,NTBAS(0),1,NTBAS(0),
     &                  NTBAS(0),NTBAS(0),1,LUPRI) 
         END DO
       END DO
      END IF


C     Dispose used memory
      CALL MEMREL('GTCNVD',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GTCNVD')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EXPVSUSLO(WORK,LWORK)
C*****************************************************************************
C
C  Routine that calculates expectation value terms for magnetic
C  susceptibilities with LAO - two-electron term and both higher and lower order
C  reorthonormalization terms.
C
C  Written by Miroslav Ilias, Prievidza, june 2003.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dcbprp.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbgen.h"
#include "dummy.h"
C
      DIMENSION WORK(*)

      CALL QENTER('EXPVSUSLO')
#include "memint.h" 
C 0. Allocate memory and read coefficients etc.
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE) ! MO coefficients
      CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','EXVSUSLO')
C ... distinguish between spin-free (boson symm.) mode
C     and the full relativistic (fermion symm.) mode
      IF (SPINFR) THEN
       IOPT = 10 ! read MO coeff.+boson irreps ident.
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,WORK(KBEIG),
     &            TOTERG,IOPT)
      ELSE
       IOPT = 2 ! read only MO coefficients
       CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),DUM,IDUM,
     &            TOTERG,IOPT)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')
C  1. Get tensor of two-electron expectation value terms 
      CALL SUSC2EL(WORK(KCMO),WORK(KBEIG),WORK,KFREE,LFREE)
C  2. Compute highest order reorthonormalization term {T(2),h(0)+g(0)}
      CALL HORDRNT(WORK(KCMO),WORK(KBEIG),WORK,KFREE,LFREE)
C  3. Compute lowest order reorthonormalization term  ...many more terms...
      CALL LORDRNT(WORK(KCMO),WORK(KBEIG),WORK,KFREE,LFREE)
C
      CALL MEMREL('EXPVSUSLO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('EXPVSUSLO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck SUSC2EL  */
      SUBROUTINE SUSC2EL(CMO,IBEIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C   Routine that calculates two-electron expectation value contribution
C to magnetic susceptibility when using London atomic orbitals.
C
C     This is the pure London-contribution to the second B-derivative.
C
C     Resulting tensor is stored in SUS2EL of dcbsusc.h
!     (gosia: SUS2EL is evaluated in LNDOUT subroutine in abacus/her2ave.F)
C
C     Written by Miroslav Ilias, Prievidza, june 2003.
C
C*****************************************************************************
      use quaternion_algebra
      use dirac_cfg
      use dft_cfg
      use xcint_main
      use memory_allocator
      use num_grid_gen
      use fde_input
      use fde_mag_cfg
!
      use fde_mod
      use fde_data
      use fde_dirac_matrices_integration
      use fde_evaluators_dirac
#include "implicit.h"
#include "priunit.h"  
#include "dummy.h"   
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "dcbbas.h"
#include "shells.h"
#include "ccom.h"
#include "cbihr2.h"
#include "symmet.h"
#include "dgroup.h"
#include "iratdef.h"
C
#include "aovec.h"
#include "blocks.h"
#include "dcbgen.h"
#include "dcbnmr.h"
#include "dcbprp.h"
#include "dcbfir.h"
C
#include "dcbsusc.h"
C
      DIMENSION CMO(*),IBEIG(*),WORK(*)
      DIMENSION ISYMOP(6),IHRMOP(6),IFCKOP(6),MTOTTK(3)
      LOGICAL   LBIT,NODV,NOPV,NOCONT,RETUR,
     &          FCKDDR,DDFOCK,FIRST,SUSCEP
      real(8), allocatable :: dmat_0(:)
      real(8) :: sus2el_xc(3, 3)
      real(8) :: sus2el_fde(3, 3)
      real(8) :: sus2el_fde_1(3, 3), sus2el_fde_2(3,3)
      real(8) :: sus2el_fde_3(3, 3)
!
      real(8), allocatable   :: dmat(:, :), buffer(:)
      type(fde_import) :: itmp

      CALL QENTER('SUSC2EL')
      KFRSAV=KFREE
C
      NDMAT = 4 ! Number of dens.matrixes (closed shell)
      NCOMP = 6 ! tensor components 
      CALL MEMGET('REAL',KDMAT,N2BBASX*NDMAT,WORK,KFREE,LFREE) ! Contravariant density matrixes
      CALL MEMGET('INTE',KIREP,NDMAT*NCOMP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFCT,NDMAT*NCOMP,WORK,KFREE,LFREE)
C     get contavariant density matrixes
      CALL GTCNVD(WORK(KDMAT),CMO,IBEIG,IPRPRP,WORK(KFREE),LFREE)
C set up for second derivative g12 over B
      DO I=1,NCOMP
         ISYMOP(I) = 1
         IHRMOP(I) = 1
         IFCKOP(I) = 1
      END DO
C prepare IFCTYP and IREPDM
      CALL SETFCK(WORK(KIFCT),WORK(KIREP),NCOMP,NDMAT,
     &            ISYMOP,IHRMOP,IFCKOP,IPRPRP)
C call TWOINT and modify routines inside her2ave
C ...Integral type - London orbitals and without EXPECT in her2ave !!!
      ITYPE = -5  
      MAXDIF = 2  ! Second derivative...=> SUSCEPTIBILITIES !!
      IATOM = 0
      RETUR  = .FALSE.
      RELCAL = .TRUE.
      DDFOCK = .TRUE.
      FCKDDR = .TRUE.
      LONDON = .TRUE.  
      SUSCEP = .TRUE.
      FIRST = .TRUE.
      SCRDUM=-1.0D0
C
C     Check what type of two-electron integrals to include
C
#if defined (VAR_MPI)
      MTOTTK(1) = NLRGBL*(NLRGBL+1)/2 ! MI...hm..some length...
      MTOTTK(2) = NLRGBL*(NLRGBL+1)/2
      MTOTTK(3) = NSMLBL*(NSMLBL+1)/2
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
#endif
      CALL GETTIM(CPU1,WALL1)
!gosia:
! i2typ = 1 (LL|LL)
!         2 (SS|LL)
!         3 (SS|SS) integrals
      DO I2TYP = 1,3
         IF (LBIT(INTNMR,I2TYP)) THEN
#if defined (VAR_MPI)
         IF (PARCAL) THEN
           CALL MEMGET('INTE',KNPOS,MTOTTK(I2TYP),WORK,KFREE,LFREE)
           CALL DIRAC_PARCTL( HERFCK_PAR )
           CALL HER_PARDRV(WORK(KFREE),LFREE,DUMMY,WORK(KDMAT),
     &          NDMAT,WORK(KIREP),WORK(KIFCT),
     &          ITYPE,MAXDIF,IATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,WORK(KNPOS),MTOTTK(I2TYP),
     &          I2TYP,IDUMMY,SCRDUM,DUMMY,DUMMY,DUMMY,DUMMY)
           CALL MEMREL('SUSC2EL.par',WORK,1,KNPOS,KFREE,LFREE)
         ELSE
#endif
           CALL MEMGET('REAL',KFMAT,N2BBASX*NDMAT*NCOMP,
     &          WORK,KFREE,LFREE)
           CALL TWOINT(WORK(KFREE),LFREE,WORK(KFMAT),WORK(KDMAT),
     &          NDMAT,WORK(KIREP),WORK(KIFCT),DUMMY,IDUMMY,IDUMMY,
     &          ITYPE,MAXDIF,IATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,IPRPRP,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY,
     &          I2TYP,IDUMMY,SCRDUM,DUMMY,DUMMY,DUMMY,DUMMY,RELCAL,
     &          .FALSE.,IDUMMY,DUMMY)
           CALL MEMREL('SUSC2EL.seq',WORK,1,KFMAT,KFREE,LFREE)
#if defined (VAR_MPI)
         ENDIF
#endif
         ENDIF
      END DO

      if (iprprp .ge. 1) then
!       sus2el before adding dft contributions:
        CALL HEADER('Two-electron expectation values, ' //
     &              'before adding XC terms',-1)
        CALL POLPRI(SUS2EL,'   ',-2)
        AVERAG = (SUS2EL(1,1) + SUS2EL(2,2) + SUS2EL(3,3))/3.0d0
        WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
      end if

!     gosia: add xc contributions:
!     ============================
!     here inside integrate_xc we calculate each component of sus2el_xc(i, j)
!     and add it to sus2el(i, j) from dcbsusc.h
!
      if (dirac_cfg_dft_calculation) then
        call alloc(dmat_0, n2bbasx*nz)
        call dzero(dmat_0, n2bbasx*nz)
        call genden(dmat_0, cmo, 1, iprprp)

        call generate_num_grid(dmat_0)
#ifdef VAR_MPI
        if (parcal) call dirac_parctl(XCINT_PAR)
#endif

!       charge density contribution dependent on XC potential
        sus2el_xc = 0.0d0
        call integrate_xc(xc_mat_dim   = ntbas(0),
     &                    xc_nz        = nz,
     &                    xc_dmat_0    = dmat_0,
     &                    xc_nr_dmat   = 0,
     &                    xc_nr_fmat   = 0,
     &                    xc_do_london_susc2el_der1 = .true.)
        sus2el_xc = get_susc2el_integrated()

        if (iprprp .ge. 1) then
          CALL HEADER('XC potential contribs to two-electron '//
     &               ' expectation values (charge density)',-1)
          CALL POLPRI(sus2el_xc,'   ',-2)
          AVERAG = (sus2el_xc(1,1) + sus2el_xc(2,2)
     &            + sus2el_xc(3,3))/3.0d0
          WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
        end if
        do i = 1, 3
           do j = 1, 3
              sus2el(j, i) = sus2el(j, i) + sus2el_xc(j, i)
           end do
        end do

        if (.not. dft_cfg_no_sdft) then
!         spin density contribution dependent on XC kernel
#ifdef VAR_MPI
          if (parcal) call dirac_parctl(XCINT_PAR)
#endif
          sus2el_xc = 0.0d0
          call integrate_xc(xc_mat_dim   = ntbas(0),
     &                      xc_nz        = nz,
     &                      xc_dmat_0    = dmat_0,
     &                      xc_nr_dmat   = 0,
     &                      xc_nr_fmat   = 0,
     &                      xc_do_london_susc2el_der2 = .true.)
          sus2el_xc = get_susc2el_integrated()

          if (iprprp .ge. 1) then
            CALL HEADER('XC kernel contribs to two-electron '//
     &                  ' expectation values (spin density)',-1)
            CALL POLPRI(sus2el_xc,'   ',-2)
            AVERAG = (sus2el_xc(1,1) + sus2el_xc(2,2)
     &              + sus2el_xc(3,3))/3.0d0
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
          end if
          do i = 1, 3
             do j = 1, 3
                sus2el(j, i) = sus2el(j, i) + sus2el_xc(j, i)
             end do
          end do
        end if

        call dealloc(dmat_0)

      end if


!     gosia: add fde contributions:
!     =============================
!
      if (fde_magn_expval_lao) then
        call fde_get_import_info(itmp)

        call alloc(dmat_0, n2bbasx*nz)
        call dzero(dmat_0, n2bbasx*nz)
        call genden(dmat_0, cmo, 1, iprprp)

        sus2el_fde   = 0.0d0
        sus2el_fde_1 = 0.0d0
        sus2el_fde_2 = 0.0d0
        sus2el_fde_3 = 0.0d0

!       embedding potential contribution
!       --------------------------------
        if (.not. fde_lao_magn_expval_no_embpot) then
          write(*, *) 'FDE 2-el contribution to expectation value ' //
     &                'dependent on embedding potential included'
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( FDE_PAR )
#endif
          if (itmp%im_vemb) then
            call fde_dirac_emb_matrices_via_integration(
     &                    fde_mat_dim            = ntbas(0),
     &                    fde_nz                 = nz,
     &                    fde_dmat_0             = dmat_0,
     &                    fde_nr_dmat            = 0,
     &                    fde_nr_fmat            = 0,
     &                    fde_use_potential      = .true.,
     &                    fde_do_london_susc2el_der1 = .true.)
          else
            call fde_dirac_emb_matrices_via_integration(
     &                      fde_mat_dim            = ntbas(0),
     &                      fde_nz                 = nz,
     &                      fde_dmat_0             = dmat_0,
     &                      fde_nr_dmat            = 0,
     &                      fde_nr_fmat            = 0,
     &                      fde_do_london_susc2el_der1  = .true.)
          end if
          sus2el_fde_1 = fde_susc2el_integrated()

          if (iprprp .ge. 1) then
            CALL HEADER('FDE 2-el contribution dependent on ' //
     &                  'embedding potential', -1)
            CALL POLPRI(sus2el_fde_1,'   ',-2)
            AVERAG = (sus2el_fde_1(1,1) + sus2el_fde_1(2,2)
     &              + sus2el_fde_1(3,3))/3.0d0
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
          end if

        end if

!       embedding kernel contribution (no coupling)
!       -------------------------------------------
        if (.not.fde_lao_magn_expval_no_uncoup_embker
     &    .and..not. fde_cfg_no_sdft) then
          write(*, *) 'FDE 2-el contribution to expectation value ' //
     &       'dependent on embedding kernel included (uncoupled)'
#ifdef VAR_MPI
          if (parcal) call dirac_parctl( FDE_PAR )
#endif     
          call fde_dirac_emb_matrices_via_integration(
     &                    fde_mat_dim            = ntbas(0),
     &                    fde_nz                 = nz,
     &                    fde_dmat_0             = dmat_0,
     &                    fde_nr_dmat            = 0,
     &                    fde_nr_fmat            = 0,
     &                    fde_do_london_susc2el_der2  = .true.)
          sus2el_fde_2 = fde_susc2el_integrated()

          if (iprprp .ge. 1) then
            CALL HEADER('FDE 2-el contribution dependent on ' //
     &                'embedding kernel (uncoupled)', -1)
            CALL POLPRI(sus2el_fde_2,'   ',-2)
            AVERAG = (sus2el_fde_2(1,1) + sus2el_fde_2(2,2)
     &              + sus2el_fde_2(3,3))/3.0d0
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
          end if

        end if

!       embedding coupling-kernel contribution
!       --------------------------------------
!       this is the coupling-kernel contribution
!       which involves the nonadditive xc and kinetic energy functionals
!       There is no contribution from the Coulomb term from LAO phase factor
!       in closed-shell systems (only the reorthonormalization term,
!       but this is taken care of in LORDRNT subroutine which uses
!       Coulomb term calculated for property gradient and stored on FCKLON file)
        if (.not.fde_lao_magn_expval_no_coupl_nonadd_embker
     &        .and. .not. fde_cfg_no_sdft) then
          if (fde_rsp_mag_lao_import) then
            write(*, *) 'FDE 2-el contribution to expectation value ' //
     &       'dependent on embedding kernel and coupling '          //
     &       'two subsystems included'
#ifdef VAR_MPI
            if (parcal) call dirac_parctl( FDE_PAR )
#endif    
            call fde_dirac_emb_matrices_via_integration(
     &                      fde_mat_dim            = ntbas(0),
     &                      fde_nz                 = nz,
     &                      fde_dmat_0             = dmat_0,
     &                      fde_nr_dmat            = 0,
     &                      fde_nr_fmat            = 0,
     &                      fde_do_london_susc2el_coupling  = .true.)
            sus2el_fde_3 = fde_susc2el_integrated()

            if (iprprp .ge. 1) then
              CALL HEADER('FDE 2-el contribution dependent on ' //
     &                'coupling embedding kernel (no Coulomb)', -1)
              CALL POLPRI(sus2el_fde_3,'   ',-2)
              AVERAG = (sus2el_fde_3(1,1) + sus2el_fde_3(2,2)
     &                + sus2el_fde_3(3,3))/3.0d0
              WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
            end if

          else ! not fde_rsp_mag_lao_import
!gosia todo:  put this quit before, now it is already after costly response equations
            call quit('You asked for coupling kernel term without '//
     &                'importing perturbed density') 
          end if
        end if

        call dealloc(dmat_0)

!       Calculate total value of sus2el_fde:
        do i = 1, 3
           do j = 1, 3
              sus2el_fde(j, i) = sus2el_fde(j, i) 
     &                         + sus2el_fde_1(j, i)
     &                         + sus2el_fde_2(j, i)
     &                         + sus2el_fde_3(j, i)
           end do
        end do

        if (iprprp .ge. 1) then
!         sus2el - fde contributions:
          CALL HEADER('Sum of FDE contribs to two-electron '//
     &                ' expectation value',-1)
          CALL POLPRI(sus2el_fde,'   ',-2)
          AVERAG = (sus2el_fde(1,1) + sus2el_fde(2,2)
     &            + sus2el_fde(3,3))/3.0d0
          WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG

        end if

!       Calculate total value of sus2el:
        do i = 1, 3
           do j = 1, 3
              sus2el(j, i) = sus2el(j, i) + sus2el_fde(j, i)
           end do
        end do

      end if !if(fde_magn_expval_lao)

!     print total sus2el:
      if (iprprp .ge. 1) then
        CALL HEADER('Two-electron expectation values (total)',-1)
        CALL POLPRI(SUS2EL,'   ',-2)
        AVERAG = (SUS2EL(1,1) + SUS2EL(2,2) + SUS2EL(3,3))/3.0d0
        WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
      end if

      CALL GETTIM(CPU2,WALL2)
C
      CALL MEMREL('SUSC2EL',WORK,1,KFRSAV,KFREE,LFREE)
      CPUTOT = CPU2 - CPU1
      WLLTOT = WALL2 - WALL1
      WRITE(LUPRI,'(/A,F13.4,A,F13.4,A/)') 
     &'Cpu time',CPUTOT,'Wall time',WLLTOT,
     &'in call to TWOINT(SUSC2EL).'
      CALL FLSHFO(LUPRI)

      CALL QEXIT('SUSC2EL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck hordrnt  */
      SUBROUTINE HORDRNT(CMO,IBEIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C   PURPOSE:
C  ----------
C   Compute the highest order expectation value
C   reorthonormalization term, {T(2), H(0)},
C   (where {T(2),H(0)} = {T(2),F(0)}) for the magnetic susceptibility 
C   with London atomic orbitals.
C   In the Levy-Leblond case this is identical
C   with the nonrelativistic DALTON counterpart.
C
C   Resulting tensor is stored in SUSFS of dcbsusc.h
C
C   Definition:
C -------------
C   The exact form (Ruud et al, JCP...) is -1/2{S(2),H(0)} 
C   (for both connections because the S(2) operator is unique)
C
C   On input:  CMO - MO coefficients, IBEIG - boson irreps info
C
C   Called from: EXPVSUSLO (prp/pammag.F)
C
C   Written by Miro Ilias,  jun 2003,  Prievidza.
C   Last revisions:
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
      PARAMETER(D0=0.0D0, D2=2.0D0, D1 = 1.0D0)
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbsusc.h"

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      LOGICAL HERM,AH,DOPCTRA_SAVE

      CALL QENTER('HORDRNT')

      KFRSAV = KFREE

C   Get the Fock matrix in MO basis
!   gosia: XC contributions are added in GETFCK subroutine, so no additional programming for DFT-LAO here...

      CALL MEMGET('REAL',KFMAT,N2ORBXQ,WORK,KFREE,LFREE)
      IPRINT = IPRPRP
      CALL GETFCK(WORK(KFMAT),IPRINT,WORK,KFREE,LFREE)
      IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &'HORDRNT: Total Fock matrix in MO basis',-1)
        CALL PRQMAT(WORK(KFMAT),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')

C ... arrays for transformed F - {T(2),F(0)}_comp 
      CALL MEMGET('REAL',KTFMAT,N2ORBXQ,WORK,KFREE,LFREE)

C  ... array for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KT2MAT,N2ORBXQ,WORK,KFREE,LFREE)

      IPRINT = IPRPRP

C ....  DO NOT PICT.CHANGE TRANSFORMATIN OF CONNECTION MATRIXES!
      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.
 
      DO ICOMP= 1, 6 ! XX,XY,XZ,YY,YZ,ZZ
       INDXPR = IPCON23(ICOMP) ! Get indx of T(2)
       ISYM   = IPRPSYM(INDXPR)
       ITIM   = IPRPTIM(INDXPR)
       IOPSY  = JBTOF(ISYM-1,1)
       IREP   = ISYM - 1
       CALL DZERO(WORK(KT2MAT),N2ORBXQ)
C    --- Get the T(2) symm.connection matrix in MO basis ---
       CALL PRPMAT(INDXPR,IOPSY,WORK(KT2MAT),.TRUE.,WORK,
     &             CMO,IBEIG,ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

C ... Call TTRA 
        AH = .FALSE. ! S(2) is hermitian ...
        ISYMT = ISYM
        ISYMH = 1 ! It's the H(0) term
        CALL TTRA(ISYMT,ISYMH,WORK(KT2MAT),
     &    WORK(KFMAT),WORK(KTFMAT),.false.,.false.,IPRINT)

        IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &'HORDRNT: Transformed Fock matrix, {T(2),F(0)}',-1)
         WRITE(LUPRI,'(2X,A,I1)')
     & '*** tensor component(1-6): ',ICOMP
         CALL PRQMAT(WORK(KTFMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP),LUPRI )
        ENDIF

C----------------------------------------------------------------
C  Calculate the expectation value for the higher order
C reorthonormalization term, SUSFS(3,3).
C----------------------------------------------------------------
       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KTFMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO
       EXPVAL = EXPVAL*D2

C       ... control print-out ...
       IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/2X,A,I2,1X,A,D15.9)')
     &  'HORDRNT: Contribution to the highest order reorthonorm. term'
     &  //' -exp.value {S(2),H(0)} of component ',ICOMP,' is:',EXPVAL
       ENDIF

       IF (ICOMP.EQ.1) THEN
         SUSFS(IPTAX(1,2),IPTAX(1,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.2) THEN
         SUSFS(IPTAX(1,2),IPTAX(2,2)) = EXPVAL
         SUSFS(IPTAX(2,2),IPTAX(1,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.3) THEN
         SUSFS(IPTAX(1,2),IPTAX(3,2)) = EXPVAL
         SUSFS(IPTAX(3,2),IPTAX(1,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.4) THEN
         SUSFS(IPTAX(2,2),IPTAX(2,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.5) THEN
         SUSFS(IPTAX(2,2),IPTAX(3,2)) = EXPVAL
         SUSFS(IPTAX(3,2),IPTAX(2,2)) = EXPVAL
       ELSE IF (ICOMP.EQ.6) THEN
         SUSFS(IPTAX(3,2),IPTAX(3,2)) = EXPVAL
       ELSE
         CALL QUIT('HORDRNT: ICOMP exceeds 6!')
       ENDIF

      ENDDO ! of XX,XY,XZ,YY,YZ,ZZ .....
      DOPCT=DOPCTRA_SAVE
C     .... release the whole memory block used in this subroutine
      CALL MEMREL('HORDRNT',WORK,1,KFRSAV,KFREE,LFREE)
C     ... close the AOPROPER file
      CLOSE(LU1INT,STATUS = 'KEEP')
      CALL QEXIT('HORDRNT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lordrrnt  */
      SUBROUTINE LORDRNT(CMO,IBEIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C   PURPOSE:
C  ----------
C   Compute the lowest order expectation value
C   reorthonormalization terms needed for calculating
C   the magnetic susceptibility with London atomic orbitals.
C   In the Levy-Leblond case these are identical
C   with the DALTON counterparts.
C
C   DEFINITIONS OF TRANSFORMED OPERATORS FOR THE EXPECTATION VALUE:
C     (Ruud et. al, Chem.Phys. Vol. 195 (1995) 157,
C                 J.Chem.Phys. Vol. 99  (1993) 3847 )
C
C        Using the symmetric connection :
C -------------------------------------------------------
C   1/2{ S(1)S(1),F(0)} - {S(1),F(1)} + 1/4{S(1),S(1),F(0)}
C
C        Using the natural connection:
C -------------------------------------------------------
C   { S(1)W(1), F(0)} - {W(1),F(1)} + {W(1),W(1),F(0)}
C
C
C Resulting tensor is stored in SUSFSY of dcbsusc.h 
C
C On input:  CMO - MO coefficients, IBEIG - boson irreps info
C
C   Written by Miro Ilias,  august 2003,  Odense.
C
C   Last revisions:
C
!gosia: both GETFCK and GETDHB contain XC contributions, so no additional
! programming for DFT-LAO needed here...
C*****************************************************************************
      use dirac_cfg
      use memory_allocator
      use num_grid_gen
      use fde_input
      use fde_mag_cfg
      use fde_mod
      use fde_data
      use fde_dirac_matrices_integration
      use fde_evaluators_dirac
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
      PARAMETER(D0=0.0D0, DI2=0.5D0, DM1=-1.0D0, D2=2.0D0, D1 = 1.0D0)
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbprp.h"
#include "dcbnmr.h"
#include "dcbxpr.h"
#include "dcbsusc.h"

      DIMENSION CMO(*),IBEIG(*),WORK(*)
      LOGICAL HERM,AH,DOPCTRA_SAVE
      real(8), allocatable :: tb_dmat(:,:), tb(:), dmat_0(:)
      real(8) :: susfsy_xc(3,3), susfsy_fde_1(3,3)
      type(fde_import) :: itmp
#include "ibtfun.h"

      CALL QENTER('LORDRNT')

      KFRSAV = KFREE
C ... allocate space for the F(0) matrix and for the MO coefficients

C   Get the Fock matrix in MO basis
      CALL MEMGET('REAL',KFMAT,N2ORBXQ,WORK,KFREE,LFREE)
      IPRINT = IPRPRP
      CALL GETFCK(WORK(KFMAT),IPRINT,WORK,KFREE,LFREE)
      IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &'LORDRNT: Total Fock matrix in MO basis',-1)
        CALL PRQMAT(WORK(KFMAT),NORBT,NORBT,NORBT,
     &            NORBT,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
C  ... arrays for T(2)_comp connection matrixes in MO basis
      CALL MEMGET('REAL',KSWMAT,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT ,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KWMAT ,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXMAT ,N2ORBXQ,WORK,KFREE,LFREE)

C****************************************************************
C. Get  { S(1)W(1), H(0) } = { S(1)W(1), F(0) } expect. value term
C or 1/2{S(1)S(1),H(0)} = 1/2{S(1)S(1),F(0)} for symm.connection
C****************************************************************
      IPRINT = IPRPRP
C ... define e-p/e-shells for reorthonormalization
C DO NOT PERFORM THE PICT.CHANGE TRANSFORMATIN OF CONNECTION MATRIXES!
      DOPCTRA_SAVE=DOPCT
      DOPCT=.FALSE.
      DO IS1 = 1, 3
       INDXPR1 = IPCON21(IS1) ! Get index of S(1)
       ISYM1   = IPRPSYM(INDXPR1)
       ITIM1   = IPRPTIM(INDXPR1)
       IOPSY1  = JBTOF(ISYM1-1,1)
       IREP1   = ISYM1 - 1
       CALL DZERO(WORK(KSMAT),N2ORBXQ)
C    --- Get the S(1) connection matrix in MO basis ---
       IF (ITIM1.NE.-1)
     &  CALL QUIT('LORDRNT:Wrong TRsymm. of S(1)!')
       CALL PRPMAT(INDXPR1,IOPSY1,WORK(KSMAT),.TRUE.,WORK,CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)

       DO IW1 = 1, 3
C       DO IW1 = IS1, 3
        IF (.NOT.SYMCON) THEN
C ... it's the natural connection
         INDXPR2 = IPCON22(IW1) ! Get index of W(1)
         ISYM2   = IPRPSYM(INDXPR2)
         ITIM2   = IPRPTIM(INDXPR2)
         IOPSY2  = JBTOF(ISYM2-1,1)
         IREP2   = ISYM2 - 1 ! 0-7
         CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis ---
         CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),.TRUE.,WORK,CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ELSE
C ... it's the symmetric connection
         INDXPR2 = IPCON21(IW1) ! Get index of S(1)
         ISYM2   = IPRPSYM(INDXPR2)
         ITIM2   = IPRPTIM(INDXPR2)
         IOPSY2  = JBTOF(ISYM2-1,1)
         IREP2   = ISYM2 - 1 ! 0-7
         CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the S(1) connection matrix in MO basis ---
         IF (ITIM2.NE.-1) CALL QUIT('LORDRNT: Wrong TR symm. of S(1)!')
         CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),.TRUE.,WORK,CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ENDIF
C  ... make a product of operators S(1)_ir1 W(1)/S(1)_ir2 - 
C     to deterimine resulting irrep, use irreps in range 0-7,
C     otherwise use MULTD2H
C ...it's based on order of irreps in DALTON (invented by J.Almlof,1970)
C       so logical function XOR works nice...
!        Direct product of two boson irreps
         IREP12 = IBTXOR(IREP1,IREP2)
         CALL QGEMM(NORBT,NORBT,NORBT,D1,
     &    'N','N',IPQTOQ(1,IREP1),
     &     WORK(KSMAT),NORBT,NORBT,NZ,
     &    'N','N',IPQTOQ(1,IREP2),WORK(KWMAT),NORBT,NORBT,NZ,
     &     D0,IPQTOQ(1,IREP12),WORK(KSWMAT),NORBT,NORBT,NZ)

         IF (IPRINT.GE.5) THEN
          IF (.NOT.SYMCON) THEN 
           CALL HEADER(
     &'LORDRNT: Transformation matrix product, S(1)*W(1)',-1)
          ELSE
           CALL HEADER(
     &'LORDRNT: Transformation matrix product, S(1)*S(1)',-1)
          ENDIF
        WRITE(LUPRI,'(2X,A,2I2)')
     & '*** tensor components: ',IS1,IW1
           CALL PRQMAT(WORK(KSWMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
         ENDIF

C ... call TTRA for {S(1)W(1), H(0)} = {S(1)W(1), F(0)}
C ... call TTRA for {S(1)S(1), H(0)} = {S(1)S(1), F(0)}

        AH  = .FALSE. ! product of S(1)W(1) is hermitian 
        ISYM12 = IREP12 + 1 ! Total symmetry of S(1)*W(1)/S(1) operators
        ISYMF=1 ! F(0) is fully symmetric
        CALL TTRA(ISYM12,ISYMF,WORK(KSWMAT),
     &        WORK(KFMAT),WORK(KXMAT),.false.,.false.,IPRINT)

        IF (IPRINT.GE.5) THEN
         IF (.NOT.SYMCON) THEN 
          CALL HEADER(
     &'LORDRNT: Transformed Fock matrix, {S(1)W(1),F(0)}',-1)
         ELSE
          CALL HEADER(
     &'LORDRNT: Transformed Fock matrix, {S(1)S(1),F(0)}',-1)
         ENDIF
        WRITE(LUPRI,'(2X,A,2I2)')
     &'*** tensor component: ',IS1,IW1
         CALL PRQMAT(WORK(KXMAT),NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
        ENDIF

C  ... save the expectation values for the lower order rn term
C      into SUSFSY(3,3)

       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KXMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO

       IF (.NOT.SYMCON) THEN
        EXPVAL = -EXPVAL*D2 ! Because of -S(1)W(1)
       ELSE
        EXPVAL = -EXPVAL ! Take care of -1/2 S(1)S(1)
       ENDIF

       IF (IPRINT.GE.5) THEN
        IF (.NOT.SYMCON) THEN  
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {S(1)W(1),H(0)} of component ',IS1,IW1,EXPVAL
        ELSE
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {S(1)S(1),H(0)} of component ',IS1,IW1,EXPVAL
        ENDIF
       ENDIF

       SUSFSY(IPTAX(IS1,2),IPTAX(IW1,2)) = EXPVAL

       ENDDO
      ENDDO

C===================================================
C   Proceed to the second term,
C -------------------------------------------------
C    - {W(1)/S(1),H(1)} = {S(1)/W(1),F(1)},
C           where F(1)=h(1)+u(1)...
C===================================================

C     ====   Get F(1)=h(1)+u(1) MO terms ====

      DO IW1 = 1, 3  ! Components of W(1)
       IF (.NOT.SYMCON) THEN
        INDXPR2 = IPCON22(IW1) ! Get index of +W(1)
        ISYM2   = IPRPSYM(INDXPR2)
        ITIM2   = IPRPTIM(INDXPR2)
        IOPSY2  = JBTOF(ISYM2-1,1)
        IREP2   = ISYM2 - 1
        CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis into WORK(KWMAT) ---
        CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),.TRUE.,WORK,CMO,IBEIG,
     &                 ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
      ELSE
        INDXPR2 = IPCON21(IW1) ! Get index of +S(1)
        ISYM2   = IPRPSYM(INDXPR2)
        ITIM2   = IPRPTIM(INDXPR2)
        IOPSY2  = JBTOF(ISYM2-1,1)
        IREP2   = ISYM2 - 1
        CALL DZERO(WORK(KWMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis into WORK(KWMAT) ---
        IF (ITIM2.NE.-1) CALL QUIT('LORDRNT:Wron TR of s(1)!')
        CALL PRPMAT(INDXPR2,IOPSY2,WORK(KWMAT),.TRUE.,WORK,CMO,IBEIG,
     &                 ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
      ENDIF
C ... get components of F(1)=h(1)+u(1) in MO basis into WORK(KSMAT)
!gosia: if we do fde, then fde terms included here are determined by flags 
!made for property gradient (and stored on TWOLON and FCKLON files)
!all the coupling terms are on FCKLON file (next GETDHB call)
        DO IH1 = 1, 3 ! Components of H(1)
          IB = 3 ! 0011 ... get ONELON + TWOLON
          INDXPR1 = IPLONDON(1,IH1) ! RM1H1_comp
          ISYM1 = IPRPSYM(INDXPR1) 
          IREP1 = ISYM1 - 1
          IREP12 = IBTXOR(IREP2,IREP1) 
          CALL DZERO(WORK(KSMAT),N2ORBXQ)
          CALL GETDHB(IB,D0,WORK(KSMAT),
     &           INDXPR1,IOPSY,WORK,KFREE,LFREE,IPRINT)

C ... call TTRA for {W(1)/S(1),F(1)} .. get all into KSWMAT
          AH  = .TRUE. ! T=W =>evidently antisymm...
          CALL TTRA(ISYM2,ISYM1,WORK(KWMAT),
     &      WORK(KSMAT),WORK(KSWMAT),.true.,.false.,IPRINT)

          IF (IPRINT.GE.5) THEN
           IF (.NOT.SYMCON) THEN
            CALL HEADER(
     &  'LORDRNT: Transformed Fock-LAO matrix,nat.conn.{W(1),F(1)}',-1)
           ELSE
            CALL HEADER(
     & 'LORDRNT: Transformed Fock-LAO matrix,symm.conn. {S(1),F(1)}',-1)
           ENDIF
           WRITE(LUPRI,'(2X,A,2I2)')
     & '*** tensor component: ',IW1,IH1
           CALL PRQMAT(WORK(KSWMAT),NORBT,NORBT,NORBT,
     &                 NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
          ENDIF

!
C  ... save the expectation values for the lower order rn term
C      into SUSFSY(3,3)

       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KSWMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO

       IF (.NOT.SYMCON) THEN
         EXPVAL = D2*D2*EXPVAL ! Times 2*
       ELSE
         EXPVAL = D2*EXPVAL ! Standard..
       ENDIF

       IF (IPRINT.GE.5) THEN
        IF (.NOT.SYMCON) THEN
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT:Contribution to the lowest order reorthonorm. term'
     &//' {W(1),H(1)} of component ',IW1,IH1,EXPVAL
        ELSE
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     & 'LORDRNT:Contribution to the lowest order reorthonorm. term'
     &//' {S(1),H(1)} of component ',IW1,IH1,EXPVAL
        ENDIF
       ENDIF

        ! store
        SUSFSY(IPTAX(IW1,2),IPTAX(IH1,2)) = 
     &  SUSFSY(IPTAX(IW1,2),IPTAX(IH1,2)) + EXPVAL
         
      ENDDO

      ENDDO

C=================================================================
C
C Proceed to the third term, {W(1),W(1),H(0)}={W(1),{W(1),H(0)}}
C
C Get F'=({W(1),F(0)}+W) contribution and  then make {W(1),F'}
C
C The same for S(1)
C
C=================================================================

      DO I1  =   1, 3
       IF (.NOT.SYMCON) THEN
C ...  get W(1)_icomp(I1) into WORK(KSMAT)
        INDXPR1 = IPCON22(I1) ! Get index of +W(1)
        ISYM1   = IPRPSYM(INDXPR1)
        ITIM1   = IPRPTIM(INDXPR1)
        IOPSY1  = JBTOF(ISYM1-1,1)
        IREP1   = ISYM1 - 1 ! 0-7
        CALL DZERO(WORK(KSMAT),N2ORBXQ)
C    --- Get the W(1) connection matrix in MO basis into WORK(KSMAT) ---
        CALL PRPMAT(INDXPR1,IOPSY1,WORK(KSMAT),.TRUE.,WORK,CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ELSE
C ...  get S(1)_icomp(I1) into WORK(KSMAT)
        INDXPR1 = IPCON21(I1) ! Get index of +S(1)
        ISYM1   = IPRPSYM(INDXPR1)
        ITIM1   = IPRPTIM(INDXPR1)
        IOPSY1  = JBTOF(ISYM1-1,1)
        IREP1   = ISYM1 - 1 ! 0-7
        CALL DZERO(WORK(KSMAT),N2ORBXQ)
C    --- Get the S(1) connection matrix in MO basis into WORK(KSMAT) ---
        IF (ITIM1.NE.-1) CALL QUIT('LORDRNT: Wrong TR symm. of S(1)!')
        CALL PRPMAT(INDXPR1,IOPSY1,WORK(KSMAT),.TRUE.,WORK,CMO,IBEIG,
     &             ICMOQ,NORB,WORK,KFREE,LFREE,IPRINT)
       ENDIF

        DO I2  =  1, 3
C ...  get FCKLON_icomp(I2) into WORK(KWMAT)
         INDXPR2 = IPLONDON(1,I2) ! Get index of RM1H1_comp
         ISYM2   = IPRPSYM(INDXPR2)
         IREP2   = ISYM2 - 1 ! 0-7
         IREP12  = IBTXOR(IREP1,IREP2) ! Final symm.
         CALL DZERO(WORK(KWMAT),N2ORBXQ)
         IB = 4 ! 0100B -> get F'=FCKLON... {T(1),F(0)}+W into WORK(KWMAT)
         CALL GETDHB(IB,D0,WORK(KWMAT),
     &          INDXPR2,IOPSY,WORK,KFREE,LFREE,IPRINT)

C ...  Do transform  F''={W(1)/S(1)_ic1,F'_ic2} into WORK(KSWMAT)
         AH  = .TRUE. ! T=W=>evidently antiherm...
         CALL DZERO(WORK(KSWMAT),N2ORBXQ)
         CALL TTRA(ISYM1,ISYM2,WORK(KSMAT),
     &    WORK(KWMAT),WORK(KSWMAT),.true.,.false.,IPRINT)

         IF (IPRINT.GE.5) THEN
          IF (.NOT.SYMCON) THEN
           CALL HEADER(
     & 'LORDRNT: TOTAL transformed matrix {W(1),W(1),F(0)}',-1)
          ELSE
           CALL HEADER(
     & 'LORDRNT: TOTAL transformed matrix {S(1),S(1),F(0)}',-1)
          ENDIF
          WRITE(LUPRI,'(2X,A,2I2)')
     & '*** tensor component: ',I1,I2
          CALL PRQMAT(WORK(KSWMAT),NORBT,NORBT,NORBT,
     &               NORBT,NZ,IPQTOQ(1,IREP12),LUPRI )
         ENDIF
C ...  Store expect. values from resulting matrix WORK(KSWMAT)
       EXPVAL = D0
       DO IFS = 1, NFSYM
       IF (NOCC(IFS).GT.0) THEN
       DO I = NPSH(IFS)+1, NPSH(IFS)+NISH(IFS) 
         IADR1 = IORB(IFS)+I
         IADR2 = (IADR1-1)*NORBT+IADR1
         EXPVAL=EXPVAL+WORK(KSWMAT+IADR2-1)
       ENDDO
       ENDIF
       ENDDO

       IF (.NOT.SYMCON) THEN
        EXPVAL = D2*EXPVAL
       ELSE
        CONTINUE  ! Consider the 1/2
       ENDIF

       IF (IPRINT.GE.5) THEN
        IF (.NOT.SYMCON) THEN
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     &'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {W(1),W(1),H(0)} of component ',I1,I2,EXPVAL
        ELSE
         WRITE(LUPRI,'(/2X,A,2I2,1X,D15.9)')
     &'LORDRNT: Contribution to the lowest order reorthonorm. term'
     &//' {S(1),S(1),H(0)} of component ',I1,I2,EXPVAL
        ENDIF
       ENDIF

       SUSFSY(IPTAX(I1,2),IPTAX(I2,2)) = 
     &  SUSFSY(IPTAX(I1,2),IPTAX(I2,2)) + EXPVAL

      ENDDO
      ENDDO

! !gosia: uncoupled FDE contributions are already present (either through GETFCK or GETDHB subroutines)
! !but we need additional code for fde coupling-kernel contributions
!       if (fde_magn_expval_lao) then
!         if (.not.fde_lao_magn_expval_no_coupl_nonadd_embker .and. 
!      &           fde_rsp_mag_lao_import) then
! 
!           write(*, *) 'FDE reorthonormalization contribution '      //
!      &       'to expectation value ' //
!      &       'dependent on embedding kernel and coupling '          //
!      &       'two subsystems included'
! 
!           call fde_get_import_info(itmp)
!         
!           call alloc(dmat_0, n2bbasx*nz)
!           call dzero(dmat_0, n2bbasx*nz)
!           call genden(dmat_0, cmo, 1, iprprp)
!           call generate_num_grid(dmat_0)
! 
! !         get T(1) matrices from AOPROPER
!           call alloc(tb_dmat, n2bbasxq, 3)
!           tb_dmat = 0.0d0
!           call alloc(tb, n2orbxq)
!           tb = 0.0d0
!           
!           do icomp = 1, 3
!             tb = 0.0d0
!             indxpr1 = ipcon(icomp)
!             isym1   = iprpsym(indxpr1)
!             irep1   = isym1 - 1
!             iopsy1  = jbtof(irep1, 1)
!           
!             call prpmat(indxpr1,iopsy1,tb,cmo,ibeig,icmoq,norb,
!      &                  work,kfree,lfree, iprint)
!           
!             call getmdm(isym1,tb_dmat(1, icomp),cmo,ibeig,tb,
!      &                  work,kfree,lfree,iprint)
!           
!           end do
! 
! #ifdef VAR_MPI
!           if (parcal) call dirac_parctl( FDE_PAR )
! #endif     
! !         inside fde_dirac_emb_matrices_via_integration subroutine
! !         we will calculate two coupling contributions
! !         - one is Omega^T(I)*Omega^T(II)
! !         - second is (R \times r)Omega(I)*Omega^T(II) + (R \times r)Omega(II)*Omega^T(I)
! !         with:
! !         Omega^T = MO(a)*MO(b)*(cTc + c.c.) with first-order connection matrix T
! !         (R \times r)Omega - LAO phase factor
!           susfsy_fde_1 = 0.0d0
!           call fde_dirac_emb_matrices_via_integration(
!      &                    fde_mat_dim            = ntbas(0),
!      &                    fde_nz                 = nz,
!      &                    fde_dmat_0             = dmat_0,
!      &                    fde_dmat               = tb_dmat,
!      &                    fde_nr_dmat            = 3,
!      &                    fde_nr_fmat            = 0,
!      &                    fde_do_london_suscreo_coupling  = .true.)
! 
!           susfsy_fde_1 = fde_susc2el_integrated()
! 
!           call dealloc(dmat_0)
!           call dealloc(tb_dmat)
!           call dealloc(tb)
! 
!           if (iprprp .ge. 1) then
!               CALL HEADER('FDE reorthonormalization contribution ' //
!      &                'dependent on ' //
!      &                'embedding kernel, coupling two subsystems', -1)
!             CALL POLPRI(susfsy_fde_1,'   ',-2)
!             AVERAG = (susfsy_fde_1(1,1) + susfsy_fde_1(2,2)
!      &              + susfsy_fde_1(3,3))/3.0d0
!             WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
! 
!           end if
!         end if
! 
! 
! !n        if (.not.fde_lao_magn_expval_no_coupl_coulomb_embker .and. 
! !n     &           fde_rsp_mag_lao_import) then
! !n
! !n          write(*, *) 'FDE contribution '      //
! !n     &       'to expectation value ' //
! !n     &       'dependent on Coulomb term in kernel and coupling '    //
! !n     &       'two subsystems included'
! !n
! !n          call fde_get_import_info(itmp)
! !n        
! !n          call alloc(dmat_0, n2bbasx*nz)
! !n          call dzero(dmat_0, n2bbasx*nz)
! !n          call genden(dmat_0, cmo, 1, iprprp)
! !n          call generate_num_grid(dmat_0)
! !n
! !n!         get T(1) matrices from AOPROPER
! !n          call alloc(tb_dmat, n2bbasxq, 3)
! !n          tb_dmat = 0.0d0
! !n          call alloc(tb, n2orbxq)
! !n          tb = 0.0d0
! !n          
! !n          do icomp = 1, 3
! !n            tb = 0.0d0
! !n            indxpr1 = ipcon(icomp)
! !n            isym1   = iprpsym(indxpr1)
! !n            irep1   = isym1 - 1
! !n            iopsy1  = jbtof(irep1, 1)
! !n          
! !n            call prpmat(indxpr1,iopsy1,tb,cmo,ibeig,icmoq,norb,
! !n     &                  work,kfree,lfree, iprint)
! !n          
! !n            call getmdm(isym1,tb_dmat(1, icomp),cmo,ibeig,tb,
! !n     &                  work,kfree,lfree,iprint)
! !n          
! !n          end do
! !n
! !n          call fde_magn_coulomb(ntbas(0), nz, 3, dmat)
! !n
! !n          call dealloc(dmat_0)
! !n          call dealloc(tb_dmat)
! !n          call dealloc(tb)
! !n        end if
! 
! !       Calculate total value of susfsy:
!         do i = 1, 3
!            do j = 1, 3
!               susfsy(j, i) = susfsy(j, i) + susfsy_fde_1(j, i)
!            end do
!         end do
! 
!       end if !if (fde_magn_expval_lao)


C Symmetrize the off-diagonal terms of SUSFSY !!

      DO I1 = 1, 2
      DO I2 = I1+1, 3
       PS =
     & (SUSFSY(IPTAX(I1,2),IPTAX(I2,2))+
     &  SUSFSY(IPTAX(I2,2),IPTAX(I1,2)))/D2
       SUSFSY(IPTAX(I1,2),IPTAX(I2,2))=PS
       SUSFSY(IPTAX(I2,2),IPTAX(I1,2))=PS
      ENDDO
      ENDDO
C     .... release the whole block used in this subroutine
      CALL MEMREL('LORDRNT',WORK,1,KFRSAV,KFREE,LFREE)
      CLOSE(LU1INT,STATUS = 'KEEP')
      DOPCT=DOPCTRA_SAVE
      CALL QEXIT('LORDRNT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck def_d1hblond */
       SUBROUTINE DEF_D1HBLOND(IMP,IPRINT)
C*****************************************************************************
C
C   PURPOSE:
C ===============
C     Define all one-electron operators
C     for magnetic property
C     gradient dH/dB|B=0 (A,B,or both)
C     when using London orbitals for the specified magnetic dependent
C     property IMP:
C
C     IMP = 1 .... NMR shielding, optical rotation (TEC)
C     IMP = 2 .... magnetic susceptibility (MI)
C     IMP = 3 ...  perturbation from one component of magnetic field
C                  (gosia: temporarily i need it in visual module)
C     IMP = 4 ...  and more (to be added later)
C
C  Written by:  MI, jan.2003, based on TEC's routines
C
C*****************************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER(DI2 = +0.5D00, DMI2 = -0.5D0, DM1 = -1.0D0,
     &          D2 = 2.0D0, D1 = 1.0D0)
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxlr.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbsusc.h"
#include "dcbnmr.h"
      CHARACTER PNAME*16, PLABEL(3)*8, LABELC*7
      DIMENSION PFAC(3),IOP(3)

      CALL QENTER('DEF_D1HBLOND')

!miro: why do we reassign the entering print level ?
      IPRINT = IPRPRP

      IF (IPRINT.GE.2) THEN
      IF (IMP.EQ.1) THEN
#ifdef MOD_LAO_REARRANGED
        IF(SHIELD.OR.shielding_rearrange) then
#else
        IF(SHIELD) then
#endif
         CALL HEADER(
     &'All individual operators for NMR shieldings with LAO',-1)
        end if
        IF(OPTROT) CALL HEADER(
     &'All individual operators for optical rotation with LAO',-1)
      ELSE IF (IMP.EQ.2) THEN
       CALL HEADER(
     &'All individual operators for magnetic susceptibility'//
     &' (magnetizability) with LAO',-1)
      ELSE IF (IMP.EQ.3) THEN
        CALL HEADER(
     &'Individual operators for perturbation from '//
     &' one component of magnetic field with LAO',-1)
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & 'DEF_D1HBLOND: Entering magnetic property IMP '//
     & ' has to be 1 or 2 or 3!'
       CALL QUIT('DEF_D1HBLOND: wrong IMP value!')
      ENDIF
      ENDIF
C------------------------------------------------
C          X component of the B-field
C------------------------------------------------
C     RM1H1 operator ... i Q_MN <r beta' c^2 > 
C
      PNAME    = 'LAO-XRM1H1'
      IPTYP    = 1
      NPCOMP   = 1
      PFAC(1)  = -D2*CVAL*CVAL
      PLABEL(1)= 'dS/dBX'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &           INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(1,1) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(1,IPLONDON(1,1))

      IF (IMP.EQ.1) THEN
        CALL OP1IND('XSHIELB',IPSHIELD(1+MXCOOR),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF (IMP.EQ.2) THEN
        CALL OP1IND('MSUSCX1',IPMSUSC(1,1),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('MSUSCX2',IPMSUSC(2,1),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF ((IMP.EQ.3) .AND. BXLAO) THEN
        CALL OP1IND('BLAOX1',IP_BLAO(1,1),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('BLAOX2',IP_BLAO(2,1),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE
        CONTINUE
      ENDIF

C--------------------------------------------------------------
C       RM1H2 operator ... i Q_MN < r_vec Vnuc >
C--------------------------------------------------------------
      PNAME    = 'LAO-XRM1H2'
      IPTYP    = 1
      NPCOMP   = 1
      PFAC(1)  = D1
      PLABEL(1)= 'XRM1H2'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(2,1) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(2,IPLONDON(2,1))
C--------------------------------------------------------------
C       RM1H3  operator ... i c Q_MN < r_vec (alp.p) >
C--------------------------------------------------------------
      PNAME     = 'LAO-XRM1H3'
      IPTYP     = 8
      NPCOMP    = 3
      PFAC(1)  = CVAL
      PFAC(2)  = CVAL
      PFAC(3)  = CVAL
      PLABEL(1) = 'XRM1H3X '
      PLABEL(2) = 'XRM1H3Y '
      PLABEL(3) = 'XRM1H3Z '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(3,1) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(2,IPLONDON(3,1))
C
C  RM1RN-X operator of +(c/2).< r_n x alpha >
C
      PNAME    = 'LAO-XRM1RN'
      IPTYP    = 5
      NPCOMP   = 2
      PFAC(1)  = 0.5d0*cval
      PFAC(2)  = 0.5d0*cval
      PLABEL(1) = 'ZRM1RN'
      PLABEL(2) = 'YRM1RN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(4,1) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(2,IPLONDON(4,1))

#if !defined (XH4)
C-------------------------------------------------------------
C            Y component of B-field
C-------------------------------------------------------------
C
C         RM1H1-Y operators
C
      PNAME     = 'LAO-YRM1H1'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = -D2*CVAL*CVAL
      PLABEL(1) = 'dS/dBY'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &       INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(1,2) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(1,IPLONDON(1,2))

      IF (IMP.EQ.1) THEN
       CALL OP1IND('XSHIELB',IPSHIELD(2+MXCOOR),LLRBPU,NLRBPT,INDXPR
     &             ,MAXLLR)
      ELSE IF (IMP.EQ.2) THEN
       CALL OP1IND('MSUSCY1',IPMSUSC(1,2),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
       CALL OP1IND('MSUSCY2',IPMSUSC(2,2),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF ((IMP.EQ.3) .AND. BYLAO) THEN
        CALL OP1IND('BLAOY1',IP_BLAO(1,2),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('BLAOY2',IP_BLAO(2,2),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE
       CONTINUE
      ENDIF
C
C        RM1H2-Y operator
C
      PNAME     = 'LAO-YRM1H2'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'YRM1H2'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(2,2) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(2,IPLONDON(2,2))
C
C       RM1H3-Y operator
C
      PNAME    = 'LAO-YRM1H3'
      IPTYP    = 8
      NPCOMP   = 3
      PFAC(1)  = CVAL
      PFAC(2)  = CVAL
      PFAC(3)  = CVAL
      PLABEL(1) = 'YRM1H3X '
      PLABEL(2) = 'YRM1H3Y '
      PLABEL(3) = 'YRM1H3Z '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &         INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(3,2) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(3,IPLONDON(3,2))
C
C          RM1RN-Y operator
C
      PNAME     = 'LAO-YRM1RN'
      IPTYP     = 6
      NPCOMP    = 2
      PFAC(1)   = 0.5d0*cval
      PFAC(2)   = 0.5d0*cval
      PLABEL(1) = 'XRM1RN'
      PLABEL(2) = 'ZRM1RN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(4,2) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(4,IPLONDON(4,2))

C-------------------------------------------------------------
C     Z component of B-field
C-------------------------------------------------------------
C
C         RM1H1-Z operator
C
      PNAME     = 'LAO-ZRM1H1'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = -D2*CVAL*CVAL
      PLABEL(1) = 'dS/dBZ'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &      INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(1,3) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(1,IPLONDON(1,3))

      IF (IMP.EQ.1) THEN
        CALL OP1IND('XSHIELB',IPSHIELD(3+MXCOOR),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF (IMP.EQ.2) THEN
        CALL OP1IND('MSUSCZ1',IPMSUSC(1,3),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('MSUSCZ2',IPMSUSC(2,3),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE IF ((IMP.EQ.3) .AND. BZLAO) THEN
        CALL OP1IND('BLAOZ1',IP_BLAO(1,3),LLRAPU,NLRAPT,INDXPR
     &     ,MAXLLR)
        CALL OP1IND('BLAOZ2',IP_BLAO(2,3),LLRBPU,NLRBPT,INDXPR
     &     ,MAXLLR)
      ELSE
       CONTINUE
      ENDIF
C
C         RM1H2-Z operator
C
      PNAME     = 'LAO-ZRM1H2'
      IPTYP     = 1
      NPCOMP    = 1
      PFAC(1)   = D1
      PLABEL(1) = 'ZRM1H2'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(2,3) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(2,IPLONDON(2,3))
C
C         RM1H3-Z operator
C
      PNAME     = 'LAO-ZRM1H3'
      IPTYP     = 8
      NPCOMP    = 3
      PFAC(1)   = CVAL
      PFAC(2)   = CVAL
      PFAC(3)   = CVAL
      PLABEL(1) = 'ZRM1H3X '
      PLABEL(2) = 'ZRM1H3Y '
      PLABEL(3) = 'ZRM1H3Z '
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &     INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(3,3) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(3,IPLONDON(3,3))
C
C        RM1RN-Z operator of +(c/2).< r_n x alpha >
C
      PNAME     = 'LAO-ZRM1RN'
      IPTYP     = 7
      NPCOMP    = 2
      PFAC(1)   = 0.5d0*cval
      PFAC(2)   = 0.5d0*cval
      PLABEL(1) = 'YRM1RN'
      PLABEL(2) = 'XRM1RN'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPLONDON(4,3) = INDXPR
      IF (IPRINT.GE.2) CALL WRIXPR(4,IPLONDON(4,3))
#endif

C----------------------------------------------------------------
C   Define connection matrix T for the first derivate dH/dB
C----------------------------------------------------------------
! for imp=3 always write the correct indxpr to ipcon(1), it is easier for a moment
      IPTYP  = 1
      NPCOMP = 1
!     gosia: if connection-independent formulation of shieldings, then
!     we want "dS/dB" integrals, we keep the '-' sign.
#ifdef MOD_LAO_REARRANGED
      if (lao_lr_rearrange) then
          PFAC(1) = -1.0d0
          LABELC='dS/dB'
          ILABELC = 5
      else
#endif
        IF (SYMCON) THEN
          PFAC(1) = DMI2
          LABELC='dS/dB'
          ILABELC = 5
        ELSE
          PFAC(1) = -D1
          LABELC='d|S>/dB'
          ILABELC = 7
        END IF
#ifdef MOD_LAO_REARRANGED
      end if
#endif
      PNAME     = LABELC(1:ILABELC)//'X'
      PLABEL(1) = LABELC(1:ILABELC)//'X'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPCON(1) = INDXPR
      if ((imp .eq. 3) .and. bxlao) then
        ipcon(1) = indxpr
      end if
      IF (IPRINT.GE.2) CALL WRIXPR(1,IPCON(1))

      PNAME     = LABELC(1:ILABELC)//'Y'
      PLABEL(1) = LABELC(1:ILABELC)//'Y'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPCON(2) = INDXPR
      if ((imp .eq. 3) .and. bylao) then
        ipcon(1) = indxpr
      end if
      IF (IPRINT.GE.2) CALL WRIXPR(2,IPCON(2))

      PNAME     = LABELC(1:ILABELC)//'Z'
      PLABEL(1) = LABELC(1:ILABELC)//'Z'
      CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
      IPCON(3) = INDXPR
      if ((imp .eq. 3) .and. bzlao) then
        ipcon(1) = indxpr
      end if
      IF (IPRINT.GE.2) CALL WRIXPR(3,IPCON(3))

! if (imp == 3) only ipcon(1) is important, lines below are to catch errors 
      if (imp .eq. 3) then
        ipcon(2) = 0
        ipcon(3) = 0
      end if

C----------------------------------------------------------------
C   Define connection matrixes T(1),T(2) for 
C   the second derivate d2H/dB2 - for magnetic suscpetibilities
C----------------------------------------------------------------
      IF (IMP.EQ.2) THEN

C ... Symmetry connection term T(1)=S(1)
       IPTYP  = 1
       NPCOMP = 1
C      PFAC(1) = DMI2
       PFAC(1) = D1
       LABELC='dS/dB'
       ILABELC = 5
       PNAME     = LABELC(1:ILABELC)//'X'
       PLABEL(1) = LABELC(1:ILABELC)//'X'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON21(1) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(1,IPCON21(1))

       PNAME     = LABELC(1:ILABELC)//'Y'
       PLABEL(1) = LABELC(1:ILABELC)//'Y'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON21(2) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(2,IPCON21(2))

       PNAME     = LABELC(1:ILABELC)//'Z'
       PLABEL(1) = LABELC(1:ILABELC)//'Z'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON21(3) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(3,IPCON21(3))

C ... Natural connection term T(1)=W(1)
C      PFAC(1) = -D1
       PFAC(1) =  D1
       LABELC='d|S>/dB'
       ILABELC = 7
       PNAME     = LABELC(1:ILABELC)//'X'
       PLABEL(1) = LABELC(1:ILABELC)//'X'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON22(1) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(1,IPCON22(1))

       PNAME     = LABELC(1:ILABELC)//'Y'
       PLABEL(1) = LABELC(1:ILABELC)//'Y'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON22(2) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(2,IPCON22(2))

       PNAME     = LABELC(1:ILABELC)//'Z'
       PLABEL(1) = LABELC(1:ILABELC)//'Z'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON22(3) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(3,IPCON22(3))

C Define d2S/dBXX,XY,...
C Specify latter PFAC(1)!!!
       IPTYP  = 1
       NPCOMP = 1
       PFAC(1) = DMI2 ! For the  "-1/2 {S(2),H(0)}" term
       LABELC='dS/dB2'
       ILABELC =6 
       PNAME     = LABELC(1:ILABELC)//'XX'
       PLABEL(1) = LABELC(1:ILABELC)//'XX'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON23(1) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(1,IPCON23(1))

       PNAME     = LABELC(1:ILABELC)//'XY'
       PLABEL(1) = LABELC(1:ILABELC)//'XY'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON23(2) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(2,IPCON23(2))

       PNAME     = LABELC(1:ILABELC)//'XZ'
       PLABEL(1) = LABELC(1:ILABELC)//'XZ'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON23(3) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(3,IPCON23(3))

       PNAME     = LABELC(1:ILABELC)//'YY'
       PLABEL(1) = LABELC(1:ILABELC)//'YY'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON23(4) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(4,IPCON23(4))

       PNAME     = LABELC(1:ILABELC)//'YZ'
       PLABEL(1) = LABELC(1:ILABELC)//'YZ'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON23(5) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(5,IPCON23(5))

       PNAME     = LABELC(1:ILABELC)//'ZZ'
       PLABEL(1) = LABELC(1:ILABELC)//'ZZ'
       CALL XPRIND(PNAME,IPTYP,NPCOMP,PFAC,PLABEL,
     &        INDXPR,ISYXPR,ITRXPR,IPRINT)
       IPCON23(6) = INDXPR
       IF (IPRINT.GE.2) CALL WRIXPR(6,IPCON23(6))
      ENDIF

      CALL QEXIT('DEF_D1HBLOND')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sum_d1hblo */
       SUBROUTINE GETDHB(IB,FAC2LO,PRPMO,INDXPR,
     &              IOPSY,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C
C    PURPOSE:
C  ==============
C   Accumulate selected terms (depends on IB)
C   of the (dH/dB)|B=0 property gradient
C   in MO-basis into the entering PRPMO MO matrix
C   when using London orbitals.
C
C   The gradient is needed for calculating the NMR shieldings and
C   the magnetizabilities.
C
C *** Terms:
C
C       dH/dB|B=0|_c = ONELON + FCKLON + (FAC2LO*)TWOLON
C
C  where ONELON = (RM1RN(sq) + RM1H1+ RM1H2 + RM1H3(sq))
C
C  Bit flags:
C ---------------
C
C    1.bit(1)  -  ONELON
C    2.bit(2)  -  TWOLON
C    3.bit(4)  -  FCKLON
C    4.bit(8)  -  FAC2LO*TWOLON  (+ 2.bit ON!)
C
C  Examples:
C    IB = 0111 (7)  - get all ONELON + FCKLON + TWOLON (full mag.field.grad.)
C    IB = 1011 (11) - get ONELON + FAC2LO*TWOLON (FAC2LO is entering)
C    IB = 1111 (15) - get ONELON + FCKLON + FAC2LO*TWOLON (FAC2LO is entering)
C    IB = 0001 (1)  - get only ONELON 
C    IB = 0101 (5)  - get only ONELON + FCKLON
C
C  The entering INDXPR for this case corresponds to some component of the
C  RM1H1 operator (X/Y/Z). 
C
C  The resulting PRPMO_c OUTPUT matrix is TR antisymmetric, what is 
C  is checked inside this routine.
C
C   Called from: GOPGET
C                GRDOR1
C
C   Written by Miro Ilias, August 2003, Odense
C 
C*****************************************************************************
#ifdef MOD_LAO_REARRANGED
      use london_helper
#endif
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,  D1 = 1.0D0)
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
C
#include "mxcent.h"
#include "dcbnmr.h"
C
!gosia, i need bxlao, bylao, bzlao
#include "dcbprp.h"
      DIMENSION PRPMO(NORBT,NORBT,NZ),WORK(*)
      INTEGER ICONVXYZ
      CHARACTER*1 COMP
      DATA THRS / 1.0D-7 /
      DATA LUONELO / 72 /
#include "ibtfun.h"

      CALL QENTER('GETDHB')

      KFRSAV = KFREE

      CALL DZERO(PRPMO, N2ORBXQ)

C     ... PRPNAM(INDXPR) is the antisymmetric 'LAO-cRM1H1'
      IF (PRPNAM(INDXPR)(6:10).NE.'RM1H1') THEN
       WRITE(LUPRI,'(2X,A,I3,A)')
     & 'GETDHB: Entering INDXPR, PRPNAM(INDXPR):',  
     &   INDXPR, PRPNAM(INDXPR)
       CALL QUIT('GETDHB: Entering wrong INDXPR!')
      ENDIF
      READ (PRPNAM(INDXPR),'(4X,A1)') COMP
      ICOMP = ICONVXYZ(COMP)

      IREP = IPRPSYM(INDXPR)-1

      CALL MEMGET('REAL',KPRPMOLO,N2ORBXQ,WORK,KFREE,LFREE)
C   ... 1. accumulate ONELON into WORK(KPRPMOLO) if 1.bit is ON
      IF (IBTAND(IB,1).EQ.1) THEN
         OPEN(LUONELO,FILE='ONELON',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='OLD')
         CALL READAC(LUONELO,N2ORBXQ,WORK(KPRPMOLO),ICOMP)
         CLOSE(LUONELO,STATUS = 'KEEP')
         IF (IBTAND(IB,8).NE.8) THEN
          CALL DAXPY(N2ORBXQ,D1,WORK(KPRPMOLO),1,PRPMO,1)
         ELSE
           WRITE(LUPRI,'(2X,A,D12.7)')
     &    'GETDHBLO: Added ONELON with FAC2LO=',FAC2LO
          CALL DAXPY(N2ORBXQ,FAC2LO,WORK(KPRPMOLO),1,PRPMO,1)
         ENDIF
         CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A)')'GETDHBLO: Added ONELON'
         ENDIF
      ENDIF
C   ... 2. accumulate TWOLON into WORK(KPRPMOLO) if 2.bit is ON
      IF (IBTAND(IB,2).EQ.2.AND.(.NOT.NOTWOL)) THEN
         OPEN(LUTWOL,FILE='TWOLON',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='OLD')
         !CALL READAC(LUTWOL,N2ORBXQ,WORK(KPRPMOLO),ICOMP)
         CALL READAC(LUTWOL,N2ORBXQ,WORK(KPRPMOLO),ICOMP)
C   ...   if 4.bit is ON, acc TWOLON with some factor FAC2LO
         IF (IBTAND(IB,8).EQ.8) THEN
          CALL DAXPY(N2ORBXQ,FAC2LO,WORK(KPRPMOLO),1,PRPMO,1)
         ELSE
          CALL DAXPY(N2ORBXQ,D1,WORK(KPRPMOLO),1,PRPMO,1)
         ENDIF
         CLOSE(LUTWOL,STATUS = 'KEEP')
         CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A)')'GETDHBLO: Added TWOLON'
          IF (IBTAND(IB,8).EQ.8) THEN
           WRITE(LUPRI,'(2X,A,D12.7)')
     &    'GETDHBLO: Added TWOLON with FAC2LO=',FAC2LO
          ENDIF
         ENDIF
      ENDIF
C   ... 3. accumulate FCKLON into WORK(KPRPMOLO) if 3.bit is ON
      if (bxlao .or. bylao .or. bzlao) then
!       we always write to component "1"
        icomp_reorth = 1
      else
        icomp_reorth = icomp
      end if
      IF (IBTAND(IB,4).EQ.4.AND.(.NOT.NOONEI)) THEN
         OPEN(LUFCKL,FILE='FCKLON',FORM='UNFORMATTED',
     &        ACCESS='DIRECT',RECL=8*N2ORBXQ,STATUS='OLD')
!gosia fixme: if here i do:
!        (1) WORK(KPRPMOLO) = 0.0d0
!        (2) CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
! i get different results.
! Results of (1) and (2) should be the same and the same as of putting 'NOONEI' in input
! while only (2) gives the same results as 'NOONEI'
         CALL READAC(LUFCKL,N2ORBXQ,WORK(KPRPMOLO),icomp_reorth)
         CALL DAXPY(N2ORBXQ,D1,WORK(KPRPMOLO),1,PRPMO,1)
         CLOSE(LUFCKL,STATUS = 'KEEP')
         CALL DZERO(WORK(KPRPMOLO),N2ORBXQ)
         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A)')'GETDHBLO: Added FCKLON'
         ENDIF
      ENDIF
C ... print the total matrix and CHECK if it is anti-symmetric !
        IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(/A,A,A,I3)')
     & 'GETDHBLO: For *** LAO-PDBS ***'
     &//' the total accumulated MO-matrix of property ',
     &  PRPNAM(INDXPR), 'with flag',IB
        CALL PRQMAT(PRPMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &              IPQTOQ(1,IREP),LUPRI)
        ENDIF
C ... check the symmetry of the total property
C      gradient matrix, PRPMO for LAO-PDBS
       DSX=SYMCHECK(PRPMO,NORBT,NORBT,NORBT,NZ)
       IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(4X,A,/,A,A,A,1X,D12.5)')
     &'*** GETDHB: The average symmetry of the final dH/dB',
     &' MO-property matrix ',PRPNAM(INDXPR),
     &' with LAO is: ',DSX
      ENDIF
      IF (DABS(DSX).GE.THRS) THEN
        WRITE(LUPRI,'(4X,A,/,A,A,A,1X,D12.5)')
     &'*** GETDHB: The average symmetry of the final dH/dB',
     &' MO-property matrix ',PRPNAM(INDXPR),
     &' with LAO is: ',DSX
        WRITE(LUPRI,'(A)')
     &'WARNING from GETDHB: The total property gradient'
     &//'matrix dH/dB for LAO-PDBS is not (fully?) symmetric !'
      ENDIF

      CALL MEMREL('GETDHB',WORK,1,KFRSAV,KFREE,LFREE)

      CALL QEXIT('GETDHB')
      RETURN
      END

      INTEGER FUNCTION ICONVXYZ(C)
C******************************************************
C 
C   Integer function returning 1,2 or 3,
C  depending on input(x,y,z)
C 
C  Written by M.Ilias, Odense, august 2003
C
C******************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*1 C

      IF (C.EQ.'X'.OR.C.EQ.'x') THEN
        ICONVXYZ=1
      ELSE IF (C.EQ.'Y'.OR.C.EQ.'y') THEN
        ICONVXYZ=2
      ELSE IF (C.EQ.'Z'.OR.C.EQ.'z') THEN
        ICONVXYZ=3
      ELSE
        WRITE(LUPRI,'(A,A)') 'ICONVXYZ: The input character:',C
        WRITE(LUPRI,'(A)') 'ICONVXYZ: and it has to be x/y/z...!!!'
        CALL QUIT('ICONVXYZ: Not x/y/z input!')
      ENDIF
      RETURN
      END


