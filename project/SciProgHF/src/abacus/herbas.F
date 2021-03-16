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

C
C
C FILE : herbas.F
C
C
C*************************************************************
C* This is the procedures that controls the reading from     *
C* file, and calls all the other subroutines.                *
C* It returns a bunch of variables. It returns:              *
C* - intent(out):
C*     IQM - (Highest angular quantum number) + 1            *
C*     JCO - Redundant variable that gives how many blocks   *
C*           the primitives and contraction coefficients are *
C*           given in (used if SPLIT is true) (array)        *
C*     NUC - Number of primitives for a given                *
C*           (quantum number + 1)                            *
C*     NRC - Number of columns with contraction coefficients *
C*           for a given (quantum number + 1).               *
C*     SEG - Gives whether the contraction coefficients are  *
C*           segmented (for a given (quantum number + 1),    *
C*           (array).                                        *
C*     ALPHA - Array with the primitives (for a given        *
C*             (quantum number + 1)).                        *
C*     CPRIMU - Matrix (3D) with the NOT normalized          *
C*              contraction coefficients (for a given        *
C*              (quantum number + 1)).                       *
C*     CPRIM - Like CPRIMU, but this time the contraction    *
C*             coefficients are normalized.                  *
C*     NBLOCK - total number of blocks; = sum(JCO(1:IQM).    *
C*              The same as IQM if not DOSPLIT.              *
C*     BASREF - documentation lines for basis set            *
C* - intent(in):
C*     KAOVEC - max NBLOCK, typically = MXSHEL.              *
C*     KPRIM  - max number of primitives.                    *
C*     Q    is nuclear charge, used to identify atom         *
C*     QEFF is eff. nuclear charge which is .le. Q,          *
C*             Q-QEFF will be described with ECP             *
C*     DSM    - highest number to be accepted as zero.       *
C*     BASNAM - name of basis set to find, or "HUCKEL"       *
C*     IPRBAS - print level in BASLIB and sub procedures     *
C*************************************************************
C
      SUBROUTINE BASLIB (IQM, JCO, NUC, NRC, SEG, ALPHA, CPRIM,
     &                   CPRIMU, NBLOCK, KAOVEC, KPRIM,
     &                   Q, QEFF, DSM, UNCONT, BASNAM, BASREF, IPRBAS)
#include "implicit.h"
#include "priunit.h"
      LOGICAL   UNCONT, SEG, READ_NU_CALLED
      DIMENSION JCO(*), NUC(KAOVEC), NRC(KAOVEC), SEG(KAOVEC),
     +          ALPHA(KPRIM, KAOVEC), CPRIMU(KPRIM, KPRIM, KAOVEC),
     +          CPRIM(KPRIM,KPRIM, KAOVEC)
! infpar.h : MYTID, MPARID
#include "infpar.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbpsi.h"
C*************************************************************
C* Variable declarations:                                    *
C*     BASNAM - A character variable that contains the name  *
C*              of the file with the basis set.              *
C*     NEWEL - Logical variable, that is returned from the   *
C*             subroutine find_pos, gives whether the next   *
C*             line in basis-file is a new element or not.   *
C*     SEGEJ - This gives whether the contraction            *
C*             coefficients are segmented or not (for        *
C*             (angular quantum-number + 1))                 *
C* External procedures:                                      *
C*     SEGORB - A procedure that checks if the contraction   *
C*              coefficients (for a given (angular quantum   *
C*              number + 1)) are segmented or not.           *
C*     NRMORB - Normalizing a matrix with contraction        *
C*              coefficients                                 *
C*************************************************************
C
      CHARACTER*(*) BASNAM
      CHARACTER*80  BASTMP, BASSAV, BASREF(10)
      LOGICAL NEWEL, SEGEJ, ANO, SADLEJ, POLFUN, FOUND, NQVD, EMSL_TYPE
      logical :: is_pople_basis_set
#ifdef PRG_DIRAC
C     Define a maximum number of primitives that can go in one block
C     for UNCONTRACTED basis sets. The value given is for s-functions
C     that generate 3 small component p-functions. The maximum for
C     other types of functions is scaled down by KHK(IQM+1) that
C     gives the number of functions in the small component for that type.
      PARAMETER (MAX_IN_BLOCK=3*15)
      LOGICAL SPLIT, DO_UNCONT, DO_SPLIT
      integer :: jcount = 0
#endif /* PRG_DIRAC */
      integer :: LU_NOT_USED
C
      CALL QENTER('BASLIB')
C
C     control print in BASLIB, usually IPRBAS = IPREAD = IPRUSR + 1 = 1
C     from calling routine
      IPRINT = IPRBAS
      IF (BASNAM(1:1).eq.' ') IPRINT = 5
      INTQ   = NINT(Q)
      IF (abs(Q-INTQ) .gt. 1.d-8) then
         WRITE (lupri,*) 'BASLIB: Q, QEFF, INTQ',Q,QEFF,INTQ
         WRITE (lupri,*) 'BASLIB: BASNAM ',BASNAM
         write (lupri,*) 'Fatal error, nuclear charge Q must be'//
     &      ' an integer when using a named basis set.'
         CALL QUIT('BASLIB ERROR, nuclear charge not an integer')
      end if
      IF (IPRINT .GE. 5) THEN
         WRITE (lupri,*) 'BASLIB: Q, QEFF, INTQ',Q,QEFF,INTQ
         WRITE (lupri,*) 'BASLIB: BASNAM ',BASNAM
      END IF
      IF (BASNAM(1:1).eq.' ') THEN
         CALL QUIT('BASLIB ERROR, basis set name starts with blank')
      END IF
      len_BASNAM = LNBLNK(BASNAM)

      BASTMP = ' '
      BASSAV = BASNAM
      LUBAS  = LU_NOT_USED(20)
      IQM    = 0
      NBLOCK = 0
      IAUG   = 0
      ANO    = .FALSE.
      SADLEJ = .FALSE.
      NQVD   = .FALSE.
      POLFUN = .FALSE.

C
C     turn on splitting of primitive shells only for methods that
C     perform an integral transformation to the spinor basis
C
      DO_SPLIT = (DOTRA .OR. DOKRMC .OR. DOMP2 .OR. DOMP2NO)
      IF (PARCAL) DO_SPLIT = .true. ! if parallel, then also split for direct SCF

CDIRAC:..by default, allow contraction only for elements Z.LE.36
C        except for 2c NR and RECP calculations
      DO_UNCONT = INTQ.GT.36 .AND. .NOT.(NONREL.OR.LEVYLE.OR.ECPCALC)
!     DO_UNCONT = DO_UNCONT .OR. (TWOCOMP.AND.(.NOT.TWOCOMPBSS))
!     hjaaj: is this still relevant after introducing NONREL above?
      DO_UNCONT = UNCONT .OR. DO_UNCONT
      IF(DO_UNCONT.AND.(.NOT.UNCONT)) THEN
!     print only from one process
      if (mytid == MPARID) then
        WRITE(LUPRI,'(/A/A,I5/2A/A/A)')
     &   '*** INFO (BASLIB) : Decontracting basis set',
     &   ' - nuclear charge: ',INTQ,
     &   ' - basis file    : ',BASNAM,
     &   'Reason: DIRAC cannot create proper small component basis',
     &   '        for contracted large component basis when Z > 36'
      end if
      ENDIF

      ! for windows we keep "-star" basis set files and not "*" files
      ! now we will translate a possible "*" basis set input
      ! first find out whether basis set is a pople basis
      is_pople_basis_set = .false.
      if (index(basnam, '3-21') > 0) then
         is_pople_basis_set = .true.
      end if
      if (index(basnam, '6-31') > 0) then
         is_pople_basis_set = .true.
      end if
      if (is_pople_basis_set) then
         ! replace possible "*" by "-star"
         index_star = index(basnam, '*')
         if (index_star > 0) then
            do i = len(basnam), index_star+5, -1
               basnam(i:i) = basnam(i-4:i-4)
            end do
            basnam(index_star:index_star+4) = '-star'
         end if
         ! replace possible second "*" by "-star"
         index_star = index(basnam, '*')
         if (index_star > 0) then
            do i = len(basnam), index_star+5, -1
               basnam(i:i) = basnam(i-4:i-4)
            end do
            basnam(index_star:index_star+4) = '-star'
         end if
      end if

      IF (BASNAM(1:6) .EQ. 'HUCKEL') THEN

C        by default do not print which basis set directory
C        is used for huckel (IPRBAS .EQ. 1 by default) /hjaaj
         IPRINT = IPRBAS - 1

C
C hjaaj jan 2000:
C        New Huckel based on ano-4 instead of STO-3G.
C        - Advantages: atomic orbitals are orthonormal,
C                      better description of atomic shells
C        - The number of contracted functions is based
C        on the available Huckel parameters in huckel.F
C        (it would be good to extend with 2p parameters for Li, Be etc.)
C        BASNAM is needed here but contains true basis set name,
C        this is restored from BASSAV just before RETURN.
C panor/johhe 2005:
C        Huckel guess now based on ANO-DK3 instead because it
C        includes all elements in the periodic table.
C hjaaj mar 2006: use ano-4 for Z.le.36 to include
C        2p for Li, Be; 3p for Na, Mg; 4p for K, Ca
C        (these basis fu. does not exist in the minimum basis ANO-DK3)

         ANO = .TRUE.
         IF (INTQ .LE. 36) THEN
            BASTMP = 'ano-4   for Huckel'
         ELSE
            BASTMP = 'ANO-DK3 for Huckel'
         END IF

C        NOTE that ECP Huckel start guess is implemented
C        in READ_ANO through the Q-QEFF parameter.
         ISTART = 7
         IF (INTQ .LE. 2) THEN
C           ... H and He
            BASNAM = 'HUCKEL    1 0 0 0'
         ELSE IF (INTQ .LE. 4) THEN
C           ... Li, Be
CHJ         BASNAM = 'HUCKEL    2 0 0 0'
            BASNAM = 'HUCKEL    2 1 0 0'
         ELSE IF (INTQ .LE. 10) THEN
C           ... B, C, N, O, F, Ne
            BASNAM = 'HUCKEL    2 1 0 0'
         ELSE IF (INTQ .LE. 12) THEN
C           ... Na, Mg
CHJ         BASNAM = 'HUCKEL    3 1 0 0'
            BASNAM = 'HUCKEL    3 2 0 0'
         ELSE IF (INTQ .LE. 18) THEN
C           ... Al, Si, P, S, Cl, Ar
CHJ-jun06   BASNAM = 'HUCKEL    3 2 0 0'
CHJ-jun06   added 3d
CKR-march 07 This is dangerous if you run with a minimal basis. Reverted.
            BASNAM = 'HUCKEL    3 2 1 0'
         ELSE IF (INTQ .LE. 20) THEN
C           ...  K, Ca
CHJ         BASNAM = 'HUCKEL    4 2 0 0'
            BASNAM = 'HUCKEL    4 3 0 0'
         ELSE IF (INTQ .LE. 30) THEN
C           ...  Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            BASNAM = 'HUCKEL    4 2 1 0'
         ELSE IF (INTQ .LE. 36) THEN
C           ...  Ga, Ge, As, Se, Br, Kr
            BASNAM = 'HUCKEL    4 3 1 0'
         ELSE IF (INTQ .LE. 38) THEN
C           ...  Rb, Sr
            BASNAM = 'HUCKEL    5 3 1 0'
         ELSE IF (INTQ .LE. 45) THEN
C           ...  Rh
            BASNAM = 'HUCKEL    5 3 2 0'
         ELSE IF (INTQ .LE. 46) THEN
C           ...  Pd
            BASNAM = 'HUCKEL    4 3 2 0'
         ELSE IF (INTQ .LE. 48) THEN
C           ...  Cd
            BASNAM = 'HUCKEL    5 3 2 0'
         ELSE IF (INTQ .LE. 54) THEN
C           ...  Xe
            BASNAM = 'HUCKEL    5 4 2 0'
         ELSE IF (INTQ .LE. 56) THEN
C           ...  Ba
            BASNAM = 'HUCKEL    6 4 2 0'
         ELSE IF (INTQ .LE. 57) THEN
C           ...  La
            BASNAM = 'HUCKEL    6 4 3 0'
         ELSE IF (INTQ .LE. 58) THEN
C           ...  Ce
            BASNAM = 'HUCKEL    6 4 3 1'
         ELSE IF (INTQ .LE. 63) THEN
C           ...  Eu
            BASNAM = 'HUCKEL    6 4 2 1'
         ELSE IF (INTQ .LE. 64) THEN
C           ...  Gd
            BASNAM = 'HUCKEL    6 4 3 1'
         ELSE IF (INTQ .LE. 70) THEN
C           ...  Yb
            BASNAM = 'HUCKEL    6 4 2 1'
         ELSE IF (INTQ .LE. 80) THEN
C           ...  Hg
            BASNAM = 'HUCKEL    6 4 3 1'
         ELSE IF (INTQ .LE. 86) THEN
C           ...  Rn
            BASNAM = 'HUCKEL    6 5 3 1'
C Manu beg 12-02-2007
         ELSE IF (INTQ .LE. 88) THEN
C           ...  Ra
            BASNAM = 'HUCKEL    7 5 3 1'
         ELSE IF (INTQ .LE. 90) THEN
C           ...  Th
            BASNAM = 'HUCKEL    7 5 4 1'
         ELSE IF (INTQ .LE. 93) THEN
C           ...  Np
            BASNAM = 'HUCKEL    7 5 4 2'
         ELSE IF (INTQ .LE. 95) THEN
C           ...  Am
            BASNAM = 'HUCKEL    7 5 3 2'
         ELSE IF (INTQ .LE. 103) THEN
C           ...  Lr
            BASNAM = 'HUCKEL    7 6 4 2'
C Manu end

         ELSE
            CALL QUIT('BASLIB ERROR: HUCKEL not defined for Z > 103')
         END IF

      ELSE

C djw Feb 2005 - notation changed from 'd' to 'd-' for consistency with EMSL
         IF (BASNAM(3:10) .EQ. 'aug-cc-p') THEN
            BASTMP(1:78) = BASNAM(3:80)

            IF (BASNAM(1:2) .EQ. 'd-' .OR.
     &          BASNAM(1:2) .eq. '2-') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:2) .EQ. 't-' .OR.
     &               BASNAM(1:2) .eq. '3-') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:2) .EQ. 'q-' .OR.
     &               BASNAM(1:2) .eq. '4-') THEN
               IAUG = 3

!radovan: added some nonstandard augmentation prefixes
!         i need them to test high-order response
!         use at own risk
            else if (basnam(1:2) .eq. '5-') then
               iaug = 4
            else if (basnam(1:2) .eq. '6-') then
               iaug = 5
            else if (basnam(1:2) .eq. '7-') then
               iaug = 6

            ELSE
               WRITE (LUPRI,'(/A2,2A)') BASNAM(1:2),' is an unknown'//
     &              ' augmentation level for ',BASNAM(1:len_BASNAM)
               CALL QUIT('BASLIB ERROR: too high augmentation level')
            END IF

         else if (basnam(3:15) == 'aug-Turbomole' .or.
     &            basnam(3:13) == 'a-Turbomole') then
           if (basnam(3:4) == 'a-') then
             bastmp(1:) = basnam(5:)
           else
             bastmp(1:) = basnam(7:)
           end if
           if      (basnam(1:2) == 's-') then
             iaug = 1
           else if (basnam(1:2) == 'd-') then
             iaug = 2
           else if (basnam(1:2) == 't-') then
             iaug = 3
           else if (basnam(1:2) == 'q-') then
             iaug = 4
           else
             write(lupri,'(//A/A/A)')
     &          'Illegal aug-Turbomole augmentation in BASLIB',
     &          ' BASNAM:',BASNAM,
     &          ' BASTMP:',BASTMP
             call quit('Illegal aug-Turbomole augmentation in BASLIB')
           end if
         else if (basnam(3:11) == 'aug-dyall') then
           bastmp(1:74) = basnam(7:80)
           if      (basnam(1:2) == 's-') then
             iaug = 1
           else if (basnam(1:2) == 'd-') then
             iaug = 2
           else if (basnam(1:2) == 't-') then
             iaug = 3
           else if (basnam(1:2) == 'q-') then
             iaug = 4
           else
             call quit('illegal aug-dyall augmentation in BASLIB')
           end if

C for backward compatibility the old notation is retained
         ELSE IF (BASNAM(2:9) .EQ. 'aug-cc-p') THEN
            BASTMP(1:79) = BASNAM(2:80)
            IF (BASNAM(1:1) .EQ. 'd') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:1) .EQ. 't') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:1) .EQ. 'q') THEN
               IAUG = 3
            ELSE
               WRITE (LUPRI,'(/A1,2A)')BASNAM(1:1),' is an unknown'//
     &              ' augmentation level for ',BASNAM(1:len_BASNAM)
               CALL QUIT('BASLIB ERROR: '//
     &            'Too high augmentation level for aug-cc-pXXX')
            END IF
         ELSE IF (BASNAM(1:4) .EQ. 'ano-') THEN
            !ano-1, ano-2, ano-3, ano-4 from MOLCAS
            ANO = .TRUE.
            BASTMP(1:len_BASNAM) = BASNAM(1:len_BASNAM)
            ISTART = len_BASNAM + 2
         ELSE IF (BASNAM(1:7) .EQ. 'ANO-RCC') THEN
            ANO = .FALSE. ! ANO-RCC is not written in MOLCAS format
            BASTMP(1:len_BASNAM) = BASNAM(1:len_BASNAM)
            ISTART = len_BASNAM + 2
         ELSE IF (BASNAM(1:7) .EQ. 'ANO-DK3') THEN
            ANO = .TRUE.
            BASTMP(1:len_BASNAM) = BASNAM(1:len_BASNAM)
            ISTART = len_BASNAM + 2
         ELSE IF (BASNAM(1:3) .EQ. 'raf') THEN
            ANO = .TRUE.
            BASTMP(1:len_BASNAM) = BASNAM(1:len_BASNAM)
            ISTART = len_BASNAM + 2
         ELSE IF (BASNAM(1:7) .EQ. 'Sadlej-') THEN
            SADLEJ = .TRUE.
            BASTMP = 'sadlej'
         ELSE IF (BASNAM(1:4) .EQ. 'NQvD') THEN
            NQVD = .TRUE.
            BASTMP(1:len_BASNAM) = BASNAM(1:len_BASNAM)
            ISTART = len_BASNAM + 2
         ELSE
            BASTMP(1:len_BASNAM) = BASNAM(1:len_BASNAM)
            ISTART = len_BASNAM + 2
         END IF
      END IF

      iaug_save = iaug
C
C     Determine if there are any user added basis functions/polarization
C     functions
C
      IPOLST = INDEX(BASNAM,'Pol')
      IF (IPOLST .GT. 0) THEN
         POLFUN = .TRUE.
         IPOLST = IPOLST + 3
      ELSE
         POLFUN = .FALSE.
      END IF

      IF (IPRINT.GE.10) THEN
        CALL HEADER('Output from BASLIB',-1)
        WRITE(LUPRI,*) 'ISTART =',ISTART
        WRITE(LUPRI,'(2A)') ' BASTMP = ',BASTMP
        WRITE(LUPRI,*) 'IPOLST =',IPOLST
        WRITE(LUPRI,*) 'POLFUN =',POLFUN
        WRITE(LUPRI,*) 'ANO    =',ANO
        WRITE(LUPRI,*) 'SADLEJ =',SADLEJ
        WRITE(LUPRI,*) 'NQVD   =',NQVD
      ENDIF
C
C Finds the right element in the file.
C
      IF (NQVD) THEN
         CALL FIND_NQD(BASNAM,INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,NRC,
     &       SEG,KPRIM,KAOVEC,DSM,POLFUN,IPOLST,UNCONT,LUBAS,IPRINT)
      ELSE
         CALL FIND_ELEMENT(BASTMP,INTQ,LUBAS,IPRINT,EMSL_TYPE)
C
C Finds the number of primitives and columns of contraction coefficients
C
         CALL FIND_POS ( NEWEL, INTEXP, INTORB, INTISG, LUBAS,
     &      IPRINT, EMSL_TYPE)
C
 10      CONTINUE
C
         IF (DO_UNCONT) THEN
            INTORB = -INTORB
            IF (ANO .OR. SADLEJ) THEN
               WRITE (LUPRI,'(/A)') ' FATAL ERROR: '//
     &            '.UNCONT not implemented for ano and sadlej yet.'
               CALL QUIT(
     &            '.UNCONT not implemented for ano and sadlej yet')
Chj-aug99: MAERKE find out how to implement UNCONT for ano or sadlej !
            END IF
         END IF
C
         IF ( .NOT. NEWEL) THEN

            READ_NU_CALLED=.FALSE.

            IQM = IQM + 1 ! counting s, p, d, ... total shells
C
C           Compute size of blocks in case we need to split large blocks of primitives.
C           Compute maximum size (should of course always be one, hence the MAX in the
C           determination step) and also check if this fits the size for the fixed dimensioned
C           alpha and cprim arrays. We have two separate hardcoded parameters, this is to make
C           sure no conflicts arise (terrible spaghetti code, but have no time to fix everything, LV, Oct 2020).
C
            KHK_SC = (IQM+1)*(IQM+2)/2 ! number of basis function components for small component l+1
            MAXPRIM_IN_BLOCK = MAX(1,MAX_IN_BLOCK/KHK_SC)
            MAXPRIM_IN_BLOCK = MIN(KPRIM,MAXPRIM_IN_BLOCK)
            N_SPLIT = (INTEXP-1)/MAXPRIM_IN_BLOCK + 1
            IF (N_SPLIT .GT. 1) MAXPRIM_IN_BLOCK = (INTEXP-1)/N_SPLIT+1 ! spread load evenly

            JCO(IQM) = 1
            ! this is an intricate loop construction that should be rewritten by someone who understands programming (LV, 2020)
  11        CONTINUE
CMI           ... to have INTORB fixed for repeating calls of READ_NU
            NBLOCK = NBLOCK + 1
            IF (NBLOCK .GT. KAOVEC) THEN
               CALL QUIT('Too many AO shell blocks; increase MXSHEL')
            END IF
C
            IF (ANO .OR. SADLEJ .OR. INTORB.GT.0) THEN
               ! splitting with ANO type sets is not possible
               ! also the intorb > 0 case is problematic for splitting
               SPLIT = .FALSE.
            ELSEIF (INTEXP.GT.KPRIM) THEN
               ! splitting is always necessary, otherwise the fixed array size is too small
               SPLIT = .TRUE.
            ELSEIF (INTEXP.GT.MAXPRIM_IN_BLOCK .AND. DO_SPLIT) THEN
               ! splitting is desired by the user to reduce memory consumption
               SPLIT = .TRUE.
            ELSE
               ! splitting is not necessary (or we are in the last block after splitting)
               SPLIT = .FALSE.
            ENDIF
           
            IF (SPLIT) THEN
               INTEXP_2 = INTEXP - MAXPRIM_IN_BLOCK
               INTEXP = MAXPRIM_IN_BLOCK
               JCO(IQM) = JCO(IQM) + 1
            ENDIF

C
C Setting ALPHA, CPRIMU AND CPRIM to zero.
C
            CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
            CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
            CALL DZERO(CPRIM (1, 1, NBLOCK), KPRIM*KPRIM)
C
C Reading the primitives and contraction coefficients from file in READ_NU
C
            IF (ANO .OR. SADLEJ) THEN
               IQCORE = NINT(Q-QEFF)
               CALL READ_ANO(INTEXP,INTORB, NBLOCK,KAOVEC, CPRIMU,ALPHA,
     +                     KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                     UNCONT,LUBAS,IQCORE,IPRINT)
               IF (INTORB .EQ. 0) THEN
                  NBLOCK = NBLOCK - 1
                  CALL FIND_POS ( NEWEL, INTEXP, INTORB,
     &                            INTISG, LUBAS, IPRINT, EMSL_TYPE)
                  GOTO 10
               END IF
            ELSE

               IF (READ_NU_CALLED) THEN
C                INTEXP = INTEXP_SAVE
                 INTORB = INTORB_SAVE
                 IF (IPRINT.GE.7) THEN
                 write(lupri,*) 'restored INTORB=',INTORB
                 ENDIF
               END IF

               if (split) then
!                 block has been split and this is not the last part;
!                 only augment last part; therefore
!                 do not augment this part
                  iaug = 0
               else
                  iaug = iaug_save
               end if

               INTEXP_SAVE = INTEXP
               INTORB_SAVE = INTORB
               CALL READ_NU(INTEXP,INTORB,NBLOCK, ALPHA,CPRIMU,
     &                      KAOVEC,KPRIM, IAUG,POLFUN,IPOLST,
     &                      LUBAS,BASNAM,IPRINT)

               READ_NU_CALLED = .TRUE.

            END IF
C
C Checking if the matrix, with the contraction coefficients, is segmented
C
            SEGEJ  = .TRUE.
            CALL SEGORB(SEGEJ,INTEXP,INTORB,CPRIMU(1,1,NBLOCK),KPRIM,
     +                  DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)

!#define PRINT_GAUSSIAN_INP
#ifdef PRINT_GAUSSIAN_INP
            CALL print_gaussian_inp(jcount, IQM, ALPHA(1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)
#endif
C
C Normalizing CPRIMU for (angular quantum number + 1).
C
             CALL NRMORB(IQM, NRC(NBLOCK), NUC(NBLOCK),
     +                   ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK),
     +                   CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
C
C Go back to read additional exponents in case of block splitting
C
            IF (SPLIT) THEN
               INTORB = 0
               INTEXP = INTEXP_2
               GOTO 11
            ENDIF
C
C Define the type of the next relevant line in file
C
            CALL FIND_POS( NEWEL, INTEXP, INTORB, INTISG,
     &         LUBAS, IPRINT, EMSL_TYPE)
            GOTO 10
         END IF ! .not. NEWEL
      END IF
C
      IF (POLFUN) THEN
 32      CONTINUE
         INTEXP = 0
         INTORB = 0
         NBLOCK = NBLOCK + 1
         FOUND = .FALSE.
         ISTART = IPOLST
 33      CONTINUE
         len_BASNAM = LNBLNK(BASNAM(ISTART:))
         IF (len_BASNAM .GT. 0) THEN
            CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT')
            CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA')
            IF (IQUANT .EQ. IQM) THEN
               FOUND  = .TRUE.
               INTEXP = INTEXP + 1
               INTORB = INTORB + 1
               ALPHA(INTEXP,NBLOCK) = EXPON
               CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
            END IF
            GOTO 33
         END IF
         IF (FOUND) THEN
            SEGEJ  = .TRUE.
            CALL SEGORB( SEGEJ, INTEXP, INTORB, CPRIMU(1,1,NBLOCK),
     +                   KPRIM, DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
            CALL NRMORB(IQM, NRC(NBLOCK), NUC(NBLOCK),
     +                  ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK),
     +                  CPRIMU(1,1, NBLOCK), KPRIM, NBLOCK)
            GOTO 32
         ELSE
            NBLOCK = NBLOCK - 1
         END IF
      END IF
C
      IF (NBLOCK .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I5,A/3X,A)')
     &   'ERROR No basis functions found for element',INTQ,' in basis:',
     &   BASNAM(1:len_BASNAM)
         IF (ANO) WRITE(LUPRI,'(/A/A)')
     &   'ERROR NB! For ano basis sets, user must specify number of',
     &   'ERROR     AOs in each shell (example: "basis=ano-1 3 2 1")'
         CALL QUIT('No basis functions found for this element')
      END IF
C
C Find reference information for the basis set.
C
      CALL FIND_REF (BASREF,LUBAS,BASSAV)
C
C Close basis set file
C

      BASNAM = BASSAV
      CLOSE (LUBAS, STATUS = 'KEEP')
      CALL QEXIT('BASLIB')
      RETURN
C
      END ! subroutine BASLIB
C
C*************************************************************
C* This is the subroutine that takes care of the reading of  *
C* the primitives and the contraction coefficients from the  *
C* file (for a given format). The variables that are         *
C* transferred out are:                                      *
C*     CPRIMU - The matrix that the contraction coefficients *
C*              are put in (not normalized).                 *
C*     ALPHA - The matrix where the primitives are put in    *
C*     (Both of the variables are given for a given (quantum *
C*     number + 1))                                          *
C*************************************************************
C
C/* Deck read_nu */
      SUBROUTINE READ_NU(INTEXP,INTORB,NBLOCK, ALPHA,CPRIMU,
     &                   KAOVEC,KPRIM, IAUG,POLFUN,IPOLST,
     &                   LUBAS,BASNAM,IPRINT)
C
C Reading the exponents and contraction coefficients
C for standard Dalton-type basis sets (i.e. not MOLCAS ANO format)
C
C Rewritten May 2016 hjaaj to free format input from LUBAS
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      REAL*8    ALPHA(KPRIM, KAOVEC), CPRIMU( KPRIM, KPRIM, KAOVEC)
      LOGICAL   POLFUN
      CHARACTER BASNAM*(*)

      IF (IPRINT.GE.5) THEN
       CALL HEADER('Output from READ_NU',-1)
       WRITE(lupri,*) 'INTEXP,INTORB,NBLOCK:',INTEXP,INTORB,NBLOCK
       WRITE(lupri,*) 'IPRINT=',IPRINT
      ENDIF

      IF (intexp .GT. KPRIM) THEN
         CALL QUIT('ERROR, too many primitives. Increase .MAXPRI')
      END IF

      INTORB_dim = ABS(INTORB)
      DO J = 1, INTEXP
          IF (INTORB .EQ. 0) THEN
C            uncontracted basis set
             READ (LUBAS, * ,IOSTAT=IOS) ALPHA(J,NBLOCK)
             CPRIMU(J,J,NBLOCK) = 1.0D0
             IF (IPRINT.GE.5) write(lupri,*)
     &         'READ_NU read 1 - NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)',
     &         NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)
          ELSE IF (INTORB .LT. 0) THEN
C            forced uncontracted with .UNCONT
             READ (LUBAS, * ,IOSTAT=IOS) ALPHA(J,NBLOCK),
     &          (DUM, I=1,INTORB_dim)
             IF (IOS.NE.0) CALL QUIT('Error in reading alpha')
             CPRIMU(J,J,NBLOCK) = 1.0D0

             IF (IPRINT.GE.5) write(lupri,*)
     &         'READ_NU read 1 - NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)',
     &         NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)

          ELSE
             READ (LUBAS, * ,IOSTAT=IOS) ALPHA(J,NBLOCK),
     &          (CPRIMU(J,I,NBLOCK), I = 1, INTORB)
          END IF
          IF (IOS.NE.0) THEN
             WRITE(LUPRI,*) 'BASLIB READ_NU: Error in reading basis set'
             WRITE(LUPRI,*) 'IOSTAT code ',IOS
             WRITE(LUPRI,*) 'BASNAM : ',TRIM(BASNAM)
             WRITE(lupri,*) 'INTEXP,INTORB,NBLOCK:',INTEXP,INTORB,NBLOCK
             WRITE(lupri,*) 'Exponent no. ',J
             CALL QUIT('BASLIB READ_NU: Error in reading basis set')
          END IF
          IF (IPRINT.GE.6) THEN
             WRITE(lupri,*)
     &          'READ_NU read - J,NBLOCK,ALPHA(J,NBLOCK)',
     &          J,NBLOCK,ALPHA(J,NBLOCK)
             IF (INTORB .gt. 0) WRITE(lupri,*)
     &          'READ_NU read - CPRIMU:',
     &          (CPRIMU(J,I,NBLOCK),I=1,INTORB)
          END IF
      END DO ! J
C
C     For uncontracted basis sets, update INTORB.
C
      IF ( INTORB .LE. 0 ) INTORB = INTEXP
C
C
C If this is some kind of augmented basis set, we augment it here, kr-96
C
         IF (IAUG .GT. 0) THEN
            do kaug = 1, iaug

!             find lowest nonzero alpha
!             note that they do not have to be monotonically decreasing!
              alpha_lowest = maxval(alpha(1:intexp, nblock))
              do iprim = 1, intexp
                if (
     &              (alpha(iprim, nblock) < alpha_lowest)
     &              .and.
     &              (alpha(iprim, nblock) > 0.0d0)
     &             ) then
                  alpha_lowest = alpha(iprim, nblock)
                end if
              end do

!             find second lowest nonzero alpha
              alpha_second_lowest = maxval(alpha(1:intexp, nblock))
              do iprim = 1, intexp
                if (
     &              (alpha(iprim, nblock) < alpha_second_lowest)
     &              .and.
     &              (alpha(iprim, nblock) > alpha_lowest)
     &             ) then
                  alpha_second_lowest = alpha(iprim, nblock)
                end if
              end do

              if (dabs(alpha_second_lowest
     &                      - alpha_lowest) < 1.0d-20) then

!               if there is only one primitive in the original set
!               then alpha_second_lowest equals alpha_lowest
!               in this case alpha_second_lowest is multiplied by 3.5
!               (factor 3.5 is my choice, no really good reason for it)
                alpha_second_lowest = 3.5d0*alpha_lowest
                write(lupri,*) 'INFO: non-standard augmentation '//
     &          'because only one primitive in block;'//
     &          ' factor 3.5 used.'
                write(lupri,*) 'spdf... block: ',NBLOCK
              end if

              intexp = intexp + 1
              intorb = intorb + 1

              IF (intexp .GT. KPRIM) THEN
               CALL QUIT('ERROR, too many primitives. Increase .MAXPRI')
              END IF
              alpha(intexp, nblock)
     &            = alpha_lowest*alpha_lowest/alpha_second_lowest
              cprimu(intexp, intorb, nblock) = 1.0d0
            end do ! kaug = 1, iaug

            IF (IPRINT.GE.5) THEN
              write(lupri,*) 'READ_NU aug:  INTEXP,NBLOCK,ALPHA',
     &         INTEXP,NBLOCK,ALPHA(INTEXP,NBLOCK)
              write(LUPRI,*) 'READ_NU: CPRIMU(INTEXP,INTORB,NBLOCK) =',
     &         CPRIMU(INTEXP,INTORB,NBLOCK)
            END IF
         END IF

         IF (POLFUN) THEN
            ISTART = IPOLST
 33         CONTINUE
            len_BASNAM = LNBLNK(BASNAM(ISTART:))
            IF (len_BASNAM .GT. 0) THEN
               CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT')
               CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA')
               IF (IQUANT .EQ. NBLOCK) THEN
                  INTEXP = INTEXP + 1
                  INTORB = INTORB + 1
                  IF (intexp .GT. KPRIM) THEN
               CALL QUIT('ERROR, too many primitives. Increase .MAXPRI')
                  END IF
                  ALPHA(INTEXP,NBLOCK) = EXPON
                  CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0

               IF (IPRINT.GE.5) THEN
                  write(lupri,*) 'POLFUN INTEXP,NBLOCK,ALPHA',
     &            INTEXP,NBLOCK,ALPHA(INTEXP,NBLOCK)
               END IF

               END IF
               GOTO 33
            END IF
         END IF

         IF (IPRINT.GE.6) THEN
            WRITE(lupri,*) 'READ_NU...NBLOCK,INTEXP=',NBLOCK,INTEXP
            WRITE(lupri,*) (ALPHA(IX,NBLOCK),IX=1,INTEXP)
         ENDIF

C
      RETURN
      END ! subroutine READ_NU
C/* Deck read_ano */
      SUBROUTINE READ_ANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA,
     +                  KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                  UNCONT,LUBAS,IQCORE,IPRINT)
C
C     IQCORE is charge of core electrons (.eq. zero if not ECP)
C     NBLOCK = l_quantum_number + 1 (1 for s, 2 for p, etc.)
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C
C Reading the exponents. Format from MOLCAS files
C
      LOGICAL   ANO, SADLEJ, POLFUN, UNCONT
      REAL*8    CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
      CHARACTER BASNAM*(*)
C
C Determine number of contracted functions to be read
C
      IF (ANO) THEN
Chjaaj aug00: prepared for adding max # of contracted orbitals
C             in ano-* files (for this consistency check).
         INTORBMAX = INTORB
         IF (INTORBMAX .LE. 0) INTORBMAX = INTEXP
         CALL FREFRM(BASNAM,ISTART,INTORB,DUMMY,'INT',IERR)
         IF (INTORB .GT. INTORBMAX) THEN
            WRITE (LUPRI,'(/2A/A,I5/A,I5)')
     &      'Fatal error in READ_ANO for basis ',BASNAM,
     &      '- max number of contracted functions',INTORBMAX,
     &      '- you asked for                     ',INTORB
            CALL QUIT('BASLIB: '//
     &         'You asked for too many contracted functions for ANO')
         END IF
      END IF

      IF (IPRINT.GE.8) THEN
        CALL HEADER('Output from READ_ANO',-1)
        wRITE(LUPRI,'(A,A)') 'BASNAM=',BASNAM
        wRITE(LUPRI,*) 'ISTART=',ISTART
        wRITE(LUPRI,*) 'INTEXP,INTORB=',INTEXP,INTORB
        WRITE(LUPRI,*) 'NBLOCK=',NBLOCK
        WRITE(LUPRI,*) ' KPRIM=',KPRIM
      ENDIF

      IF (INTEXP .GT. KPRIM) THEN
         CALL QUIT('ERROR, too many primitives. Increase .MAXPRI')
      END IF

!
!      ISTART_ANO=10
!     IF (ANO) CALL FREFRM(BASNAM,ISTART,INTORB,DUMMY,'INT')
!     IF (ANO) CALL FREFRM(BASNAM,ISTART_ANO,INTORB,DUMMY,'INT')

!     IF (IPRINT.GE.8) THEN
!       WRITE(LUPRI,'(/,8x,a)') 'READ_ANO - After FREFRM:'
!       WRITE(LUPRI,'(A,A)') 'BASNAM=',BASNAM
!       WRITE(LUPRI,*) 'ISTART=',ISTART
!       WRITE(LUPRI,*) 'INTEXP,INTORB=',INTEXP,INTORB
!     ENDIF

C
C If ECP calculate number of core orbitals to exclude for each l-shell
C This time we have l-shell number NBLOCK
C
      CALL ECP_LCORE(IQCORE,NBLOCK,ILOFF)
C
C Reading the primitive and contracted coefficients
C
      IF (INTORB .GT. 0) THEN
         READ (LUBAS,*,IOSTAT=IOS) (ALPHA(K,NBLOCK), K = 1, INTEXP)
         IF (IOS.NE.0) CALL QUIT('BASLIB: '//
     &      'Error in reading ANO exponents from basis set file')

         IF (UNCONT .AND. BASNAM(1:6) .NE. 'HUCKEL') THEN
C        ... force uncontracted
            DO I = 1, INTEXP
               READ (LUBAS,*)
               CPRIMU(I,I,NBLOCK) = 1.0D0
            END DO
            INTORB = INTEXP
         ELSE IF (ILOFF .GT. 0) THEN
C           skip core orbitals for Huckel guess for ECP
            INTORB = INTORB - ILOFF
            DO I = 1, INTEXP
               READ (LUBAS,*,IOSTAT=IOS) (XXX, J = 1,ILOFF),
     &          (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
            END DO
         ELSE
            DO I = 1, INTEXP
               READ (LUBAS,*,IOSTAT=IOS)
     &          (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
               IF (IOS.NE.0) CALL QUIT('BASLIB: Error in reading '//
     &          'ANO contraction coefficients from basis set file')
            END DO
         END IF

      ELSE
         READ (LUBAS,*) (XXX, K = 1, INTEXP)
         DO I = 1, INTEXP
            READ (LUBAS,*)
         END DO
      END IF
      IF (POLFUN) THEN
         IREAD = IPOLST
 33      CONTINUE
         len_BASNAM = LNBLNK(BASNAM(IREAD:))
         IF (len_BASNAM .GT. 0) THEN
            CALL FREFRM(BASNAM,IREAD,IQUANT,DUMMY,'INT')
            CALL FREFRM(BASNAM,IREAD,IDUMMY,EXPON,'REA')
            IF (IQUANT .EQ. NBLOCK) THEN
               INTEXP = INTEXP + 1
               INTORB = INTORB + 1
               IF (intexp .GT. KPRIM) THEN
               CALL QUIT('ERROR, too many primitives. Increase .MAXPRI')
               END IF
               ALPHA(INTEXP,NBLOCK) = EXPON
               CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
            END IF
            GOTO 33
         END IF
      END IF

      IF (IPRINT.GE.10) THEN
        WRITE(LUPRI,*) 'READ_ANO-read coefficients...'
        WRITE(LUPRI,*) 'NBLOCK,INTEXP:',NBLOCK,INTEXP
        WRITE(LUPRI,*) 'INTORB:',INTORB
        DO I = 1, INTEXP
          WRITE(LUPRI,*)
     &   (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
        ENDDO
      ENDIF

      RETURN
      END
C*********************************************************
C* This is the subroutine that searches through the file *
C* and finds the element in question.                    *
C*********************************************************
C/* Deck find_element */
      SUBROUTINE FIND_ELEMENT(BASNAM,INTQ,LUBAS,IPRINT,EMSL_TYPE)
#include "implicit.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     in:  BASNAM - The name of the basis set file      *
C*     out: EMSL_TYPE - true if EMSL type basis set      *
C*     local:                                            *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST  - Logical variable that helps inquire if    *
C*             there exists a file with that name.       *
C*     SLASH - Variable for slash in paths (/ or \\)     *
C*     WIN   - Check if we work on Linux or on Windows   *
C*********************************************************
C
      CHARACTER*(*) BASNAM
      LOGICAL       EMSL_TYPE
      CHARACTER*300 STRING  ! long because contains the full directory path
      CHARACTER SIGN
      LOGICAL EXST, EXST2, USE_BASNAM
! infpar.h : MYTID, MPARID
#include "inftap.h"
#include "infpar.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif

#if defined (SYS_WINDOWS)
      CHARACTER, PARAMETER :: SLASH = '\\' ! Specify trailing slash: \\ for Windows
#else
      CHARACTER, PARAMETER :: SLASH = '/'  ! Specify trailing slash: / for Linux
#endif

      CHARACTER*5   EMSL_NAME_CHECK
      CHARACTER*12  EMSL_ELEMENT_NAME(96)
      DATA EMSL_ELEMENT_NAME(1:96) /
     & 'HYDROGEN','HELIUM',
     & 'LITHIUM','BERYLLIUM',
     & 'BORON','CARBON','NITROGEN','OXYGEN','FLUORINE','NEON',
     & 'SODIUM','MAGNESIUM',
     & 'ALUMINUM','SILICON','PHOSPHOROUS','SULFUR','CHLORINE','ARGON',
     & 'POTASSIUM','CALCIUM',
     & 'SCANDIUM','TITANIUM','VANADIUM','CHROMIUM','MANGANESE',
     &   'IRON','COBALT','NICKEL','COPPER','ZINC','GALLIUM',
     & 'GERMANIUM','ARSENIC','SELENIUM','BROMINE','KRYPTON',
     & 'RUBIDIUM','STRONTIUM',
     & 'YTTRIUM','ZIRCONIUM','NIOBIUM','MOLYBDENUM','TECHNETIUM',
     &   'RUTHENIUM','RHODIUM','PALLADIUM','SILVER','CADMIUM',
     & 'INDIUM','TIN','ANTIMONY','TELLURIUM','IODINE','XENON',
     & 'CESIUM','BARIUM',
     & 'LANTHANUM','CERIUM','PRASEODYMIUM','NEODYMIUM',
     &   'PROMETHIUM','SAMARIUM','EUROPIUM',
     &   'GADOLIUM','TERBIUM','DYSPROSIUM','HOLMIUM','ERBIUM',
     &   'THULIUM','YTTERBIUM',
     & 'LUTETIUM','HAFNIUM','TANTALUM','TUNGSTEN','RHENIUM',
     &   'OSMIUM','IRIDIUM','PLATINUM','GOLD','MERCURY',
     & 'THALLIUM','LEAD','BISMUTH','POLONIUM','ASTATINE','RADON',
     & 'FRANCIUM','RADIUM',
     & 'ACTINIUM','THORIUM','PROTACTINIUM','URANIUM','NEPTUNIUM',
     &   'PLUTONIUM','AMERICIUM','CURIUM'/
C
!     CALL FIND_BASFIL(BASNAM,LUBAS,IPRINT)
C
C     Linux:
C     The environment variable is a : separated string, e.g.:
C     .;/programs/Dirac/basis:/programs/Dalton/basis:.
C
C     Windows:
C     The environment variable is a ; separated string, e.g.:
C     .;C:\programs\Dirac\basis;C:\programs\Dalton\basis
C
C     We loop over these directories in search for the basis set and
c     element given.
C
C     Depending on the machine, the string has trailing blanks
C     (ASCII char=32) or null elements (ASCII=0).
      ISTART = 1
      IBLANK=INDEX(BASDIR,CHAR(32))
      INULL =INDEX(BASDIR,CHAR(0))
      ISTOP=MAX(IBLANK,INULL)
      IF(ISTOP.EQ.0) THEN
C       If string is empty add current directory
        BASDIR(1:1) = '.'
        ISTOP = 1
      ELSE
        ISTOP=ISTOP-1
      ENDIF

!     Step through the specified directories for search of
!     basis set. If basis set not found in this directory;
!     the code returns to statement number 10 to find the
!     next directory to search.
 10   CONTINUE
C
#if defined (SYS_WINDOWS)
C     Find next semicolon in BASDIR
          ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:ISTOP), ';' )
#else
C     Find next colon in BASDIR
          ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:ISTOP), ':' )
#endif
C
C     If there are no more colons (: or ;), we are the last directory in BASDIR
      IF ( ICOLON .EQ. ISTART - 1 ) ICOLON = ISTOP + 1
      IEND = ICOLON - 1
      STRING = BASDIR(ISTART:IEND)//' '
C
      IEND = INDEX ( STRING(1:), ' ') - 1
!     stefan: make check with idb indicates a possible access penalty if
!     BASDIR is not set...
      if(iend.lt.1)then
        write(lupri,'(/a)') '  *** error in FIND_ELEMENT: '//
     &        'environment variable BASDIR is not set'
        call quit('*** error in FIND_ELEMENT: BASDIR is not set')
      end if
      IF ( STRING(IEND:IEND) .NE. SLASH ) THEN
         IEND = IEND + 1
         STRING(IEND:IEND) = SLASH
      END IF
      STRING = STRING(1:IEND)//BASNAM
      IEND = INDEX(STRING(1:),' ') - 1
Cluuk Activate when debugging basis set search (should not be on by default)
CLuuk WRITE(LUPRI,'(A/,A1,A,A1)')
CLuuk&     ' Searching basis : trying file: ','"',STRING(1:IEND),'"'
C
C     Inquire if "/path/basis-set" exists.
C
      INQUIRE (FILE = STRING, EXIST = EXST)
CMI  ... on certain systems (like PC-linux,OpenMPI-gfortran) string
C     beginning with '/' is NOT the proper file name ! therefore
C     continue in inquiring:
      EXST2=EXST
      USE_BASNAM=.FALSE.
      IF (.NOT.EXST2) THEN
        INQUIRE (FILE = BASNAM, EXIST = EXST)
        IF (EXST) THEN
          if (mytid == MPARID) then
          ! print only from one process
             WRITE(LUPRI,"(/2A/2A)")
     &     "Info: used proper basis set name: ",BASNAM,
     &     "Info: instead of                : ",STRING
          end if
          USE_BASNAM=.TRUE.
        ENDIF
      ENDIF

      IF (EXST) THEN
        IF (.NOT.USE_BASNAM) THEN
         OPEN(LUBAS, FILE=STRING)
        ELSE
         OPEN(LUBAS, FILE=BASNAM)
        ENDIF
C
C
C Searching the file for the element.
C Traditional Dalton basis set layout with and 'a' or 'A' in first column
C to identify next element.
C
      EMSL_TYPE = .FALSE.
 20   CONTINUE
         READ(LUBAS,'(A)', IOSTAT = IOERR, ERR = 2000, END = 100)
     +       STRING
         READ (STRING, '(A1)',IOSTAT=IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN       ! <- traditional Dalton basis set file
            READ (STRING, '(BN,A1,I4)',IOSTAT=IOS) SIGN, NUCEL
            IF (IOS.NE.0) THEN
               CALL QUIT('FIND_ELEMENT: Error in reading NUCEL')
            ENDIF
            IF (INTQ .EQ. NUCEL) RETURN
         END IF
      GOTO 20
C
C Searching the file for the element again. This time with
C the EMSL Dalton basis set layout without an 'a' or 'A' in first column.
C Element must be identified by name from the table EMSL_ELEMENT_NAME.
C
  100 IF (INTQ .GT. 0 .AND. INTQ .LE. 96) THEN
         EMSL_NAME_CHECK = EMSL_ELEMENT_NAME(INTQ)(1:5)
      ELSE
         GOTO 200
      END IF
      EMSL_TYPE = .TRUE.
      REWIND LUBAS
      NUCEL = -1
 30   CONTINUE
         READ(LUBAS,'(A)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         IF (INDEX(STRING,EMSL_NAME_CHECK) .GT. 0) THEN  ! <- basis set file made by EMSL
            ! we do not know how many lines the element name is repeated in,
            ! therefore we continue until it is not found in order that FIND_POS can work correctly
            ! (EMSL always print "! s functions" for the first block, so this is OK)
 35         READ(LUBAS,'(A)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +          STRING
            IF (INDEX(STRING,EMSL_NAME_CHECK) .GT. 0) GO TO 35
            RETURN
         END IF
      GOTO 30
C
      END IF ! EXST
!
! Basis set not found in this basis set directory,
! find the next directory and search there.
! (note: basis set directories are saved as "dir1:dir2:dir3:etc")
!
      ISTART = ICOLON + 1
      IF(ISTART.GT.ISTOP) GOTO 300
      GOTO 10
 300  CONTINUE
      LSTR = 1
      DO J = LEN(BASNAM), 1, -1
         IF ( BASNAM(J:J) .NE. ' ' ) THEN
            LSTR = J
            GOTO 301
         END IF
      END DO
 301  CONTINUE
      WRITE (LUPRI,'(/3A)')
     &'ERROR: Basis set "',BASNAM(1:LSTR),'" not found.'//
     &' Check spelling or whether DIRAC can locate basis set directory.'
      IF (BASNAM(1:6) .EQ. 'HUCKEL') THEN
         WRITE (LUPRI,'(A)')
     &     'NB! This basis set must be available to use Huckel guess!'
      END IF
      WRITE(LUPRI,'(2X,A,A)') 'Basis name, STRING=',STRING
      WRITE(LUPRI,'(2X,A,A)') 'Basis name, BASNAM=',BASNAM
      CALL QUIT('Basis set not found. Check whether DIRAC '//
     &          'can locate basis set libraries.')
C
C Error messages
C
 2000 CONTINUE
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FIND_ELEMENT: Error when reading from basis file ',
     &   STRING(1:IEND),' IOSTAT =',IOERR
      CALL QUIT('I/O error in BASLIB(FIND_ELEMENT)')
 200  CONTINUE
      len_BASNAM = LNBLNK(BASNAM)
      WRITE (LUPRI,'(//A,I5,2A/A//)')
     &   ' ERROR: Nuclear charge',INTQ,
     &   ' is an unsupported element for basis ',BASNAM(1:len_BASNAM),
     &' ERROR: You must choose another basis set, or enter it manually.'
      CALL QUIT('BASLIB: '//
     &   'Unsupported element in specified basis set file')
      END
C/* Deck find_nqd */
      SUBROUTINE FIND_NQD(BASNAM, INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,
     &                    NRC,SEG,KPRIM,KAOVEC,DSM,POLFUN,IPOLST,
     &                    UNCONT,LUBAS,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST - Logical variable that helps inquire if     *
C*            there exists a file with that name.        *
C*********************************************************
C
C
      CHARACTER*(*) BASNAM
      CHARACTER*300 STRING
      CHARACTER SIGN
      LOGICAL   POLFUN, SEG, SEGEJ, UNCONT, EXST
      DIMENSION ALPHA(KPRIM,KAOVEC), CPRIMU(KPRIM,KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC), NUC(KAOVEC), NRC(KAOVEC),
     &          SEG(KAOVEC)
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     Determine primitive basis
C
      ISTART = 5
      CALL FREFRM(BASNAM,ISTART,ISPRIM,DUMMY,'INT')
      IF (INTQ .GT. 4) THEN
         CALL FREFRM(BASNAM,ISTART,IPPRIM,DUMMY,'INT')
      END IF
      CALL FREFRM(BASNAM,ISTART,ISCONT,DUMMY,'INT')
      IF (INTQ .GT. 4) THEN
         CALL FREFRM(BASNAM,ISTART,IPCONT,DUMMY,'INT')
      END IF
      STRING = BASDIR(1:LENBAS)//'NQvD'
      INQUIRE (FILE = STRING, EXIST = EXST)
      IF (EXST) THEN
         OPEN(LUBAS, FILE = STRING)
C
C
C Searching the file for the element.
C
 20   CONTINUE
         READ(LUBAS,'(A)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)', IOSTAT = IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            ISTART = 2
            CALL FREFRM(STRING,ISTART,NUCEL,DUMMY,'INT')
            CALL FREFRM(STRING,ISTART,IS,DUMMY,'INT')
            IF (NUCEL .EQ. INTQ .AND. IS .EQ. ISPRIM) THEN
               IF (NUCEL .GT. 4) THEN
                  CALL FREFRM(STRING,ISTART,IP,DUMMY,'INT')
                  IF (IP .NE. IPPRIM) GO TO 20
               END IF
               NBLOCK = NBLOCK + 1
               IF (NBLOCK .GT. KAOVEC) THEN
                  CALL QUIT('FIND_NQD: Too many AO blocks'//
     &               ' for this program version; increase WRKMEM')
               END IF
C
C Setting ALPHA, CPRIMU AND CPRIM to zero.
C
               CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
               CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
               CALL DZERO(CPRIM (1,1, NBLOCK), KPRIM*KPRIM)
               IF (NUCEL .GT. 4) THEN
                  CALL DZERO(ALPHA(1, NBLOCK + 1), KPRIM)
                  CALL DZERO(CPRIMU(1,1, NBLOCK + 1), KPRIM*KPRIM)
                  CALL DZERO(CPRIM (1,1, NBLOCK + 1), KPRIM*KPRIM)
               END IF
               IF (NUCEL .LE. 4) THEN
                  ICONT = 1
                  DO ISLOOP = 1, ISPRIM
                     READ (LUBAS,*,IOSTAT=IOS)
     &                ALPHA(ISLOOP,NBLOCK), CONTS
                     IF (IOS.NE.0) THEN
                       CALL QUIT('FIND_NQD: Error reading ALPHA, CONTS')
                     ENDIF
                     IF (NUCEL .EQ. 1) ALPHA(ISLOOP,NBLOCK) =
     &                    ALPHA(ISLOOP,NBLOCK)*1.44D0
                     IF (ICONT .LT. ISCONT) THEN
                        CPRIMU(ISLOOP,ICONT,NBLOCK) = 1.0D0
                        ICONT = ICONT + 1
                     ELSE
                        CPRIMU(ISLOOP,ICONT,NBLOCK) = CONTS
                     END IF
                  END DO
               ELSE
                  ICONTS = 1
                  ICONTP = 1
                  DO IPLOOP = 1, IPPRIM
                     READ (LUBAS,*,IOSTAT=IOS)
     &                     ALPHA(IPLOOP,NBLOCK), CONTS,
     &                     ALPHA(IPLOOP,NBLOCK + 1), CONTP
                     IF (IOS.NE.0) THEN
                       CALL QUIT('FIND_NQD Error reading ALPHA,CONTS..')
                     ENDIF
                     IF (ICONTP .LT. IPCONT) THEN
                        CPRIMU(IPLOOP,ICONTP,NBLOCK + 1) = 1.0D0
                        ICONTP = ICONTP + 1
                     ELSE
                        CPRIMU(IPLOOP,ICONTP,NBLOCK + 1) = CONTP
                     END IF
                     IF (ICONTS .LT. ISCONT) THEN
                        CPRIMU(IPLOOP,ICONTS,NBLOCK) = 1.0D0
                        ICONTS = ICONTS + 1
                     ELSE
                        CPRIMU(IPLOOP,ICONTS,NBLOCK) = CONTS
                     END IF
                  END DO
                  DO ISLOOP = IPLOOP, ISPRIM
                     READ (LUBAS,*,IOSTAT=IOS)
     &                 ALPHA(ISLOOP,NBLOCK), CONTS
                     IF (IOS.NE.0) THEN
                        CALL QUIT('FIND_NQD Error reading ALPHA,CONTS')
                     ENDIF
                     IF (ICONTS .LT. ISCONT) THEN
                        CPRIMU(ISLOOP,ICONTS,NBLOCK) = 1.0D0
                        ICONTS = ICONTS + 1
                     ELSE
                        CPRIMU(ISLOOP,ICONTS,NBLOCK) = CONTS
                     END IF
                  END DO
               END IF
               IF (POLFUN) THEN
                  ISTART = IPOLST
 33               CONTINUE
                  len_BASNAM = LNBLNK(BASNAM(ISTART:))
                  IF (len_BASNAM .GT. 0) THEN
                     CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT')
                     CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA')
                     IF (IQUANT .EQ. NBLOCK) THEN
                        ISPRIM = ISPRIM + 1
                        ISCONT = ISCONT + 1
                        ALPHA(ISPRIM,NBLOCK) = EXPON
                        CPRIMU(ISPRIM,ISCONT,NBLOCK) = 1.0D0
                     ELSE IF (IQUANT .EQ. (NBLOCK + 1)
     &                       .AND. NUCEL .GT. 4) THEN
                        IPPRIM = IPPRIM + 1
                        IPCONT = IPCONT + 1
                        ALPHA(IPPRIM,NBLOCK + 1) = EXPON
                        CPRIMU(IPPRIM,IPCONT,NBLOCK + 1) = 1.0D0
                     END IF
                     GOTO 33
                  END IF
               END IF
C
C Checking if the matrix, with the contraction coefficients, is segmented
C
               SEGEJ  = .TRUE.
               CALL SEGORB(SEGEJ, ISPRIM, ISCONT, CPRIMU(1,1,NBLOCK),
     +                     KPRIM, DSM)
               NUC(NBLOCK) = ISPRIM
               NRC(NBLOCK) = ISCONT
               SEG(NBLOCK) = SEGEJ
C
C              Reorder primitive orbitals
C
               CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),
     &                     NUC(NBLOCK),NRC(NBLOCK),SEG(NBLOCK),
     &                     KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
               CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK),
     +                     ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK),
     +                     CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               IF (NUCEL .GT. 4) THEN
C
C Checking if the matrix, with the contraction coefficients, is segmented
C
                  NBLOCK = NBLOCK + 1
                  SEGEJ = .TRUE.
                  CALL SEGORB(SEGEJ, IPPRIM, IPCONT,
     +                        CPRIMU(1,1,NBLOCK),KPRIM, DSM)
                  NUC(NBLOCK) = IPPRIM
                  NRC(NBLOCK) = IPCONT
                  SEG(NBLOCK) = SEGEJ
C
C                 Reorder primitive orbitals
C
                  CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),
     &                        NUC(NBLOCK),NRC(NBLOCK),SEG(NBLOCK),
     &                        KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
                  CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK),
     +                        ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK),
     +                        CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               END IF
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE
            GOTO 20
         END IF
      ELSE
         LSTR=1
         DO J = 140,1,-1
            IF (STRING(J:J) .NE. ' ') THEN
               LSTR = J
               GO TO 300
            END IF
         END DO
  300    WRITE (LUPRI,'(/,3A)')
     +      'Basis "',STRING(1:LSTR),'" doesn''t exist'
      CALL QUIT('Basis set not found. Check whether DIRAC '//
     &          'can locate basis set libraries.')
      END IF
C
C Error messages
C
 2000 CONTINUE
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FIND_NQD: Error when reading from basis file ',
     &   BASDIR(1:LENBAS)//'NQvD',
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FIND_NQD')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in HERBAS')
      END
C/* Deck find_ref */
      SUBROUTINE FIND_REF(BASREF, LUBAS, BASNAM)
C*********************************************************
C* This is the subroutine that searches through the file *
C* and finds reference info if available.                *
C* Based on FIND_POS, LV, 14-5-2003                      *
C*********************************************************
#include "implicit.h"
#include "priunit.h"
#ifdef VAR_MPI
#include "infpar.h"
#endif

      CHARACTER BASNAM*(*), BASREF(10)*80
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*********************************************************
      CHARACTER SIGN, LINE*80, STRING*100
C
      REWIND (LUBAS)
C
C Initialising IREF (number of documentation lines), set to one
C once we find a line containing the string REFERENCE
C
      IREF = 0
C
C     The loops are written f77 style since f90 may still not be
C     available everywhere, better to use while loops here..
C
 20   CONTINUE
      READ (LUBAS, '(A)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      IF (IOERR .NE. 0) GO TO 200
      READ (STRING, '(A1)',IOSTAT=IOS) SIGN
C
      IF (SIGN .EQ. '$' .OR. SIGN .EQ. '!') THEN
         IF (IREF.EQ.0) THEN
C           Check for start of documentation
            IF (INDEX(STRING,'REFERENCE').GT.0) IREF = 1
         ELSE
C           We have documentation, omit blank lines and store in BASREF
            len_STRING = LNBLNK(STRING(2:))
            IF (len_STRING .GT. 0) THEN
               BASREF(IREF) = STRING(2:)
               IREF = IREF + 1
            ENDIF
         END IF
C
C        Maximum number of documentation lines is 10
C
         IF (IREF.GT.10) RETURN
      ELSE
C        We return if we found a REFERENCE line but no more $/! lines.
         IF (IREF.GT.0) RETURN
      END IF
C
      GOTO 20
C
C No documentation when end of file is found
C
 200  CONTINUE
 1000 CONTINUE
      END
C
C
C*********************************************************
C* This is the subroutine that finds out whether it is a *
C* new block with primitives and contraction             *
C* coefficients, or if it is a new element. It returns   *
C* three variables:                                      *
C*          NEWEL -Logical variable that gives whether it*
C*                 is a new element or not.              *
C*          INTEXP, INTORB - Number of primitives and    *
C*                          columns of contraction       *
C*                          coefficients                 *
C*********************************************************
C
C/* Deck find_pos */
      SUBROUTINE FIND_POS(NEWEL, INTEXP, INTORB, INTISG,
     &   LUBAS, IPRINT)
#include "implicit.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*********************************************************
      LOGICAL   NEWEL
      CHARACTER SIGN
      CHARACTER*300 STRING

C
C Initializing NEWEL, because one is innocent until proven guilty
C
      NEWEL = .FALSE.
C
 20   CONTINUE
      READ (LUBAS, '(A)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      READ (STRING, '(A1)',IOSTAT=IOS) SIGN
C
      IF (SIGN .EQ. ' ') THEN
C
         len_STRING = LNBLNK(STRING)
         IF (len_STRING .LE. 0) THEN
            GOTO 20
         ELSE
            READ (STRING, '(BN,3I5)',IOSTAT=IOS) INTEXP, INTORB, INTISG
            IF (IOS.NE.0) THEN
               CALL QUIT('ERROR in reading basis; '//
     &         'check your mol file formatting')
            END IF
            IF (IPRINT.GE.10) THEN
             write(lupri,'(/,2X,A)') 'Output from FIND_POS...reading'
             write(lupri,*) 'INTEXP, INTORB, INTISG',
     &          INTEXP, INTORB, INTISG
            ENDIF
            RETURN
         END IF
      ELSE IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
         NEWEL = .TRUE.
         IF (IPRINT.GE.10) THEN
            WRITE(lupri,'(/A)') '  Output from FIND_POS: New element'
            WRITE(lupri,*) 'INTEXP, INTORB, INTISG',
     &         INTEXP, INTORB, INTISG
         END IF
         RETURN
      ELSE
         GOTO 20
      END IF
C
C No more orbitals when end of file, so the orbitals for the last element
C are read
C
 200  CONTINUE
      NEWEL = .TRUE.
C
C NEWEL must be .TRUE. to break the if-loop.
C
      IF (IPRINT.GE.10) THEN
        WRITE(lupri,'(/,2X,A)') 'Output from FIND_POS'
        WRITE(lupri,*) 'INTEXP, INTORB, INTISG',
     &          INTEXP, INTORB, INTISG
      END IF

      RETURN
C
 1000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'ERROR in reading basis set file, '//
     &     'your basis is not complete. IOSTAT code:',IOERR
      CALL QUIT('BASLIB FIND_POS: Incomplete basis set file')
      END

      SUBROUTINE print_gaussian_inp(icount, iqm,ALPHA,CPRIMU,NPRI,
     &                              NRCI,SEG,KPRIM,DSM)
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
      LOGICAL SEG
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
      character(1) :: ch
      if (iqm == 1) then
         icount = icount + 1
         if (icount > 1) then
            write(*, '(a)') '****'
         end if
         write(*, '(a2, i2)') namn(icount), 0
      end if
C
C     *****************************
C     *** Segmented contraction ***
C     *****************************
C
      IDONE = 0
      IF (SEG) THEN
         DO 100 ICONTR = 1, NRCI
            ISTART = IDONE + 1
            NLEFT  = NPRI - IDONE
C
C           Find first primitive
C           ====================
C
            IMXA = ISTART + IDAMAX(NLEFT,ALPHA(ISTART),1) - 1
            CALL DSWAP(1,ALPHA(ISTART),1,ALPHA(IMXA),1)
            CALL DSWAP(NRCI,CPRIMU(ISTART,1),KPRIM,CPRIMU(IMXA,1),
     &                 KPRIM)
C
C           Find corresponding contracted function
C           ======================================
C
            IMXC = IDAMAX(NRCI,CPRIMU(ISTART,1),KPRIM)
            CALL DSWAP(NPRI,CPRIMU(1,ICONTR),1,CPRIMU(1,IMXC),1)
C
C           Collect other primitives contributing to same contracted
C           ========================================================
C
            IPRI = 1
            DO 200 I = ISTART + 1, NPRI
               IF (ABS(CPRIMU(I,ICONTR)) .GT. DSM) THEN
                  CALL DSWAP(1,ALPHA(I),1,ALPHA(ISTART+IPRI),1)
                  CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                            CPRIMU(ISTART+IPRI,1),KPRIM)
                  IPRI = IPRI + 1
               END IF
  200       CONTINUE
C
C           Sort primitives
C           ===============
C
            IF (IPRI .GT. 2) THEN
               DO 300 I = ISTART + 1, ISTART + IPRI - 2
                  DO 400 J = I + 1, ISTART + IPRI - 1
                     IF (ALPHA(J) .GT. ALPHA(I)) THEN
                        CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
                        CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                                  CPRIMU(J,1),KPRIM)
                     END IF
  400             CONTINUE
  300          CONTINUE
            END IF
C
            IDONE = IDONE + IPRI
  100    CONTINUE
C
C     ***************************
C     *** General contraction ***
C     ***************************
C
      ELSE
         DO 500 I = 1, NPRI - 1
            DO 600 J = I + 1, NPRI
            IF(ALPHA(J) .GT. ALPHA(I)) THEN
               CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
               CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,CPRIMU(J,1),KPRIM)
           END IF
  600      CONTINUE
  500   CONTINUE
      END IF

      select case (iqm)
         case (1)
            ch = 'S'
         case (2)
            ch = 'P'
         case (3)
            ch = 'D'
         case (4)
            ch = 'F'
         case (5)
            ch = 'G'
         case (6)
            ch = 'H'
         case default
            call quit('print_gaussian_inp: too high l quantum number')
      end select

      do i = 1, npri
         write(*, '(a2, a)') ch, '   1 1.00       0.000000000000'
         write(*, '(e22.10, e18.10)') alpha(i), 1.0d0
      end do

      END
C/* Deck ECP_LCORE */
      SUBROUTINE ECP_LCORE(IQCORE,IQM,ILOFF)
C
C     Nov. 2003, Johan, Patrick and Hans Joergen.
C
C     Purpose: calculate ILOFF = number of core orbitals described by ECP
C     for l_quantum_number = IQM - 1
C     when IQCORE electrons are described with the ECP.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ILCORE(4)
C
      IF (IQCORE .GT. 0 .AND. IQM .LE. 4) THEN
C        only s, p. d, f (never g, h, ... core orbitals)
         IF (MOD(IQCORE,2) .NE. 0) THEN
            WRITE(LUPRI,'(/A,I5)')
     &      'ECP_LCORE ECP: odd number of core electrons',IQCORE
            CALL QUIT('BASLIB ECP_LCORE ECP error: '//
     &         'odd number of core electrons')
         END IF
         DO I = 1, 4
            ILCORE(I) = 0
         END DO
         IF (IQCORE .EQ. 2) THEN
C           1s ECP
            ILCORE(1) = 1
         ELSE IF (IQCORE .EQ. 10) THEN
C           1s,2s,2p ECP
            ILCORE(1) = 2
            ILCORE(2) = 1
         ELSE IF (IQCORE .EQ. 18) THEN
C           1s,2s,2p,3s,3p ECP
            ILCORE(1) = 3
            ILCORE(2) = 2
         ELSE IF (IQCORE .EQ. 28) THEN
C           1s,2s,2p,3s,3p,3d ECP
            ILCORE(1) = 3
            ILCORE(2) = 2
            ILCORE(3) = 1
         ELSE IF (IQCORE .EQ. 36) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p ECP
            ILCORE(1) = 4
            ILCORE(2) = 3
            ILCORE(3) = 1
         ELSE IF (IQCORE .EQ. 46) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p,4d ECP
            ILCORE(1) = 4
            ILCORE(2) = 3
            ILCORE(3) = 2
         ELSE IF (IQCORE .EQ. 60) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p,4d,4f ECP
            ILCORE(1) = 4
            ILCORE(2) = 3
            ILCORE(3) = 2
            ILCORE(4) = 1
         ELSE IF (IQCORE .EQ. 78) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d ECP
            ILCORE(1) = 5
            ILCORE(2) = 4
            ILCORE(3) = 3
            ILCORE(4) = 1
         ELSE
            WRITE(LUPRI,'(/A,I5)') 'READ_ANO FATAL ERROR: '//
     &         'ECP Huckel not implemented for this core charge:',IQCORE
            CALL QUIT(
     &         'READ_ANO: ECP not implemented for this core charge')
         END IF
         ILOFF = ILCORE(IQM)
      ELSE
         ILOFF = 0
C        not ECP or g, h, ...
C        i.e. no core orbitals to exclude
      END IF
      RETURN
      END
C  -- end of herbas.F --
