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
***********************************************************************

      subroutine c_to_ckk(C,CKK,IWAY,NJA,NJB,NKA,NKB,NTA,NTB,
     &                    IKJA,XKJA,IKJB,XKJB)
*
* IWAY = 1 : C(Ja,Jb) => C(Ka,Kb,ITa,ITb)
* IWAY = 2 : C(Ka,KB,Ita,Itb) => C(Ja,Jb) + C(Ka,Kb,Ita,Itb)
*
* Jeppe Olsen, Summer of 99
* Version where matrix notation is eliminated
* Notice : K and T is in general batched. This is not visible here
!
! simplification an f90-fied by Stefan Knecht - April 2011
*
      implicit none

!----------------------------------------------------------------------
      integer, intent(in)    :: IWAY
      integer, intent(in)    :: NJA
      integer, intent(in)    :: NJB
      integer, intent(in)    :: NKA
      integer, intent(in)    :: NKB
      integer, intent(in)    :: NTA
      integer, intent(in)    :: NTB
      integer, intent(in)    :: IKJA(NTA*NKA) 
      integer, intent(in)    :: IKJB(NTB,NKB)
      real(8), intent(in)    :: XKJA(NTA*NKA) 
      real(8), intent(in)    :: XKJB(NTB,NKB)
      real(8), intent(inout) :: C(NJA*NJB) 
      real(8), intent(inout) :: CKK(*)
!----------------------------------------------------------------------
      integer                :: ita, itb
      integer                ::  ka,  kb
      integer                ::  ja,  jb
      integer                :: kakbitab, kakbitabx
      integer                :: itaka
      integer                :: itab
!======================================================================

!     print '(a,4i8)', ' NJA,NJB,NKA,NKB ==> ',nja,njb, nka, nkb
!     print '(a,2i8)', ' NTA,NTB         ==> ',nta,ntb
!     print '(a   )', ' vector IKJA'
!     CALL IWRTMA(IKJA,1,NTA*NKA,1,NTA*NKA)
!     print '(a   )', ' vector IKJB'
!     CALL IWRTMA(IKJB,NTB,NKB,NTB,NKB)
!     print '(a,i8)', ' nta*ntb*nka*nkb  ==> ',nta*ntb*nka*nkb
 
!     initialize
      if(iway.eq.1)then 
        call dzero(ckk,nta*ntb*nka*nkb)
      end if

      DO ITB = 1, NTB
        DO ITA = 1, NTA
          KAKBITAB = ((ITB-1)*NTA + ITA - 1) * NKA * NKB
          DO KB = 1, NKB         
            JB = IKJB(ITB,KB)
            if(JB /= 0)then
              KAKBITABx = KAKBITAB
              ITAKA     = ITA-NTA
              do ka = 1, nka           
                KAKBITABx = KAKBITABx + 1
                ITAKA     = ITAKA + NTA
                JA        = IKJA(ITAKA)
                if(ja /= 0)then
                  if(iway.eq.1)then
                    CKK(KAKBITABx)     = (JB/abs(JB)) * (JA/abs(JA))
!                   CKK(KAKBITABx)     = XKJA(ITAKA) * XKJB(ITB,KB) 
     &                                 * C((abs(JB)-1)*NJA + abs(JA))
                  else
                    C((abs(JB)-1)*NJA  + abs(JA)) 
!    &                                 = XKJA(ITAKA) * XKJB(ITB,KB)
     &                                 = (JB/abs(JB)) * (JA/abs(JA))
     &                                 * CKK(KAKBITABx)
     &                                 + C((abs(JB)-1)*NJA + abs(JA))
                  end if ! iway switch
                end if ! ja > 0
              end do ! ka loop
            end if ! jb > 0
            KAKBITAB  =  KAKBITAB + NKA
          end do ! kb loop
        end do ! nta loop
      end do ! ntb loop
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,*) ' Output from C_TO_CKK, IWAY =  ', IWAY
      WRITE(6,*) ' ======================================='
      WRITE(6,*)
      WRITE(6,*) ' C(Ja,Jb) '
      CALL WRTMAT(C,NJA,NJB,NJA,NJB)
      WRITE(6,*) ' C(Ka, Kb, ITa, Itb) '
      ITAB = 0 
      DO ITB = 1, NTB
      DO ITA = 1, NTA
        WRITE(6,*) ' C(Ka, Kb, Ita, Itb) for Ita,Itb =', Ita, Itb
        CALL WRTMAT(CKK(1+ITAB*NKA*NKB),NKA,NKB,NKA,NKB)
        ITAB = ITAB + 1
      END DO
      END DO
#endif
!#undef LUCI_DEBUG
*
      END 
***********************************************************************

      SUBROUTINE C_TO_CKK_old(C,CKK,IWAY,NJA,NJB,NKA,NKB,NTA,NTB,
     &                        IKJA,XKJA,IKJB,XKJB)
*
* IWAY = 1 : C(Ja,Jb) => C(Ka,Kb,ITa,ITb)
* IWAY = 2 : C(Ka,KB,Ita,Itb) => C(Ja,Jb) + C(Ka,Kb,Ita,Itb)
*
* Jeppe Olsen, Summer of 99
*
* Version where matrix notation is eliminated
*
* Notice : K and T is in general batched. This is not visible here
*
*. Input
#include "implicit.inc"
      INTEGER    IKJA(NTA*NKA), IKJB(NTB,NKB)             
      DIMENSION  XKJA(NTA*NKA), XKJB(NTB,NKB)             
*. Input/Output 
C     DIMENSION C(NJA,NJB),CKK(NKA,NKB,*)
      DIMENSION C(NJA*NJB),CKK(*)

!     WRITE(LUWRT,*)'NJA,NJB, NKA, NKB',NJA,NJB, NKA, NKB
!     WRITE(LUWRT,*)'NTA,NTB ',NTA,NTB
!     WRITE(LUWRT,*)'vector IKJA'
!     CALL IWRTMA(IKJA,1,NTA*NKA,1,NTA*NKA)
!     WRITE(LUWRT,*)'vector IKJB', IKJB
!     CALL IWRTMA(IKJB,NTB,NKB,NTB,NKB)
 
      DO ITB = 1, NTB
        DO ITA = 1, NTA
          ITAB     = (ITB-1)*NTA+ITA-1
          KAKBITAB = ITAB*NKA*NKB
          WRITE(6,*)'ITAB',ITAB
          IF(IWAY.EQ.1) THEN
            DO KB = 1, NKB         
              JB = IKJB(ITB,KB)
              IF(JB.NE.0) THEN
                WRITE(6,*)'JB',JB
                SIGNB = XKJB(ITB,KB)
                JAJB0 = (JB-1)*NJA 
                ITAKA = ITA-NTA
                DO KA = 1, NKA           
                  ITAKA = ITAKA + NTA
                  JA    = IKJA(ITAKA)
                  WRITE(6,*)'JA, ITAKA, NTA',JA, ITAKA, NTA
                  KAKBITAB = KAKBITAB + 1
                  IF(JA.NE.0) THEN
                    JAJB = JAJB0 + JA
                    SIGN = XKJA(ITAKA)*SIGNB
              WRITE(6,*)'KAKBITAB, JAJB',KAKBITAB, JAJB, JAJB0, JA
                    CKK(KAKBITAB) = SIGN*C(JAJB)
                  ELSE 
                    CKK(KAKBITAB) = 0.0D0
                  END IF
                END DO
              ELSE !IF (JB.EQ.0) THEN
                DO KAKBITAB_EFF =  KAKBITAB+1,  KAKBITAB + NKA  
                  CKK(KAKBITAB_EFF) = 0.0D0
                END DO
                KAKBITAB  =  KAKBITAB + NKA
              END IF ! JB .NE. 0
            END DO
          ELSE  ! iway == 2
            DO KB = 1, NKB         
              JB = IKJB(ITB,KB)
              IF(JB.NE.0) THEN
                JAJB0 = (JB-1)*NJA
                SIGNB = XKJB(ITB,KB)
                ITAKA = ITA-NTA
                DO KA = 1, NKA           
                  ITAKA    = ITAKA + NTA
                  JA       = IKJA(ITAKA)
                  KAKBITAB = KAKBITAB + 1
                  IF(JA.NE.0) THEN
                    JAJB    = JAJB0 + JA
                    SIGN    = XKJA(ITAKA)*SIGNB
                    C(JAJB) = C(JAJB) + SIGN*CKK(KAKBITAB)
                  END IF
                END DO
              ELSE
                KAKBITAB =  KAKBITAB + NKA
              END IF ! if JB = 0
            END DO
          END IF ! IWAY switch
        END DO
      END DO
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,*) ' Output from C_TO_CKK, IWAY =  ', IWAY
      WRITE(6,*) ' ======================================='
      WRITE(6,*)
      WRITE(6,*) ' C(Ja,Jb) '
      CALL WRTMAT(C,NJA,NJB,NJA,NJB)
      WRITE(6,*) ' C(Ka, Kb, ITa, Itb) '
      ITAB = 0 
      DO ITB = 1, NTB
      DO ITA = 1, NTA
        WRITE(6,*) ' C(Ka, Kb, Ita, Itb) for Ita,Itb =', Ita, Itb
        CALL WRTMAT(CKK(1+ITAB*NKA*NKB),NKA,NKB,NKA,NKB)
        ITAB = ITAB + 1
      END DO
      END DO
#endif
!#undef LUCI_DEBUG
*
      END 
***********************************************************************

      SUBROUTINE COMPRS2LST_B(I1I,XI1I,N1,I2I,XI2I,N2,NKIN,NKOUT,
     &                        I1O,XI1O,I2O,XI2O)
*
* Two lists of excitations/annihilations/creations are given.
* Compress to common nonvanishing entries
*
* Jeppe Olsen, July 2000 from COMPRS2LST 
*
* Compared to COMPRS2LST : Order of arrys in I1 and I2 interchanged
*                          Output lists differs from input lists
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input lists
!     DIMENSION I1I(N1,NKIN),XI1I(N1,NKIN)
!     DIMENSION I2I(N2,NKIN),XI2I(N2,NKIN)
      DIMENSION I2I(N2,NKIN), I1I(N1,NKIN)
*. Output lists
!     DIMENSION I1O(N1,NKIN),XI1O(N1,NKIN)
!     DIMENSION I2O(N2,NKIN),XI2O(N2,NKIN)
      DIMENSION I2O(N2,NKIN), I1O(N1,NKIN)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1I(I,K) /= 0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2I(I,K) /= 0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
            DO I = 1, N1
               I1O(I,NKOUT) = I1I(I,K)!*XI1I(I,K)
!             XI1O(I,NKOUT) =XI1I(I,K)
            END DO
            DO I = 1, N2
               I2O(I,NKOUT) = I2I(I,K)!* XI2I(I,K)
!             XI2O(I,NKOUT) =XI2I(I,K)
            END DO
        END IF
      END DO
*
      END
***********************************************************************

      SUBROUTINE K_TO_J_TOT(IKJ,XKJ,KSM,IK_B,IK_E,IT,NTOP,
     &                      IT_B,IT_E,IM,XM,IBM,NK,LTOP,IZERO_MAP,IUBL)
*
* Obtain Total map !J> = T-Oper !K>
* T-oper : string of elementary operators 
* !K> String of sym KSM
*
* !K> are of sym KSM and are restricted to IK_B to IK_E ( within sym)
* T_operators are restricted to IT_B to IT_E
*
* Map for each elementary operator is provided by 
*  IM with offset mat IBM
*
* Jeppe Olsen, Summer of 99
* Double group symmetry incorporated for relativistic calculations
* IUBL flag for unbarred/barred operator strings introduced.
*  Timo Fleig, July 2001, January 2002
*
      use mospinor_info
      use symmetry_setup_krci
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      INTEGER IT(LTOP,*),IM(*),
     &        IBM(MXNDGIRR,MXNDGIRR,LTOP),NK(MXPOBS,LTOP)
      DIMENSION XM(*)
*. Output
      INTEGER IKJ(NTOP,(IK_E-IK_B+1))
!     DIMENSION XKJ(NTOP,(IK_E-IK_B+1))
      DIMENSION XKJ(*)
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
#include "ctcc.inc"
      WRITE(6,*) ' Information from K_TO_J.... '
      WRITE(6,*) ' =========================== '
      WRITE(6,*) ' LTOP, NTOP = ', LTOP,NTOP
      WRITE(6,*) ' IK_B,IK_E,IT_B,IT_E = ', IK_B,IK_E,IT_B,IT_E
      WRITE(6,*) ' IUBL, KSM = ',IUBL,KSM
      do i = 1, maxlen_i1
        if(IM(i).ne.0) print *,' nonzero im at i ==> ',i,im(i)
      end do
#endif
!#undef LUCI_DEBUG
*
      IZERO_MAP = 1
      IF(LTOP.GT.0) THEN
        DO KSTR = IK_B, IK_E
!          WRITE(6,*) ' KSTR = ', KSTR
          DO ITOP = IT_B, IT_E
C?          WRITE(6,*) ' ITOP = ', ITOP
            KNSM = KSM
            SIGN = 1.0D0
            KNSTR = KSTR
            DO IOP = 1, LTOP
!              WRITE(6,*) ' Info for IOP/ITOP = ', IOP,ITOP
              JOB_ABS = IT(IOP,ITOP)
              if (IUBL.eq.1) then
                JOB_SM = ISMFTO(JOB_ABS)
                JOB_TP = ITPFTO(JOB_ABS)
              else if (IUBL.eq.2) then
                JOB_SM = ISMFTO2(JOB_ABS)
                JOB_TP = ITPFTO2(JOB_ABS)
              end if
!              WRITE(6,*) ' JOB_ABS, IBTS, JOB_SM ',
!     &        JOB_ABS, JOB_TP, JOB_SM
              if (IUBL.eq.1) then
                JOB_REL = JOB_ABS - IOBPTS(JOB_TP,JOB_SM) + 1
              else if (IUBL.eq.2) then
                JOB_REL = JOB_ABS - IOBPTS2(JOB_TP,JOB_SM) + 1
              end if
              LK = NK(KNSM,IOP)

!             stefan - apr 2011: restricting assignment in accordance
!             with loop restrictions in MAP_EXSTR.
              if(lk.gt.0)then
                IADR  = IBM(JOB_SM,KNSM,IOP) -1 + (JOB_REL-1)*LK+KNSTR
                KNSTR = iabs(IM(IADR))
!               WRITE(6,*) ' LK ==> ',LK
!               WRITE(6,*) ' IADR, JOB_SM, KNSM, KNSTR, JOB_REL-1', 
!     &                      IADR, JOB_SM, KNSM, KNSTR, JOB_REL-1
!               SIGN = SIGN*XM(IADR) 
                if(iabs(IM(IADR)) > 0)then
                  SIGN = SIGN*(IM(IADR)/iabs(IM(IADR)))
                else
                  sign = 0.0d0
                end if
              else
                KNSTR = 0
              end if

              IF(KNSTR.EQ.0) exit
              KNSM = IDBGMULT(KNSM,JOB_SM)
            END DO
*
            IF(KNSTR.NE.0) IZERO_MAP = 0
!            print *,' assign KNSTR to ITOP-IT_B+1,KSTR-IK_B+1',
!     &                       KNSTR,   ITOP-IT_B+1,KSTR-IK_B+1
            IKJ(ITOP-IT_B+1,KSTR-IK_B+1) = KNSTR * SIGN
!           XKJ(ITOP-IT_B+1,KSTR-IK_B+1) = SIGN
          END DO
        END DO
      ELSE
*. No K-operators => Identity map
        IZERO_MAP = 0
        DO KSTR = IK_B,  IK_E
          IKJ(1,KSTR-IK_B+1) = KSTR
!         XKJ(1,KSTR-IK_B+1) = 1.0D0
        END DO
      END IF
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,*) ' Operators '
      print*,'IT_B and IT_E are ==> ',IT_B, IT_E
      print*,'IK_B and IK_E are ==> ',IK_B, IK_E
      CALL IWRTMA(IT(1,IT_B),LTOP,IT_E-IT_B+1,LTOP,IT_E-IT_B+1)
      WRITE(6,*) ' Output from K_TO_J_TOT ' 
      WRITE(6,*) 'IKJ array'
      LT = IT_E-IT_B+1
      LK = IK_E-IK_B+1
      CALL IWRTMA(IKJ,LT,LK,LT,LK)
#endif
!#undef LUCI_DEBUG
      END
***********************************************************************

      SUBROUTINE MXRESC_REL(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                      NSMST,NSTFSMSPGP,NSTFSMSPGP2,
     &                      MXNDGIRR,
     &                      NSMOB,MXPTOB,NTPOB,NTSOB,NTSOB2,
     &                      NTESTG,MXPKA,
     &                      NEL1234,
     &                      MXCJ,MXCIJA,
     &                      MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
*
* Largest block of single excitations MXSXBL

*. Input
* IAB :allowed combination of alpha and beta supergroups
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups

      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXNDGIRR,*),NSTFSMSPGP2(MXNDGIRR,*)
      DIMENSION NTSOB(MXPTOB,NSMOB),NTSOB2(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
*. Note : Only done for alpha-strings, problems when transposing
*         constructing C(J,Ia,Kb)
      MXCJ = 0
      MXADKBLK = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ka,Ib)'
        write(6,*)
      end if
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB = MAX(MXB,NSTFSMSPGP2(ISM,IBTPABS))
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB = ', MXB
            DO 300 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
!               DO 320 ISMOB = 1, NSMOB
                DO 320 ISMOB = NSMOB/2 + 1,NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
*
* matrix C(j,Ia,Kb)
*
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ia,Kb)'
        write(6,*)
      end if
*
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            DO ISM = 1, NSMST
              MXA = MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IBTPABS,1,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO KSM = 1, NSMST
                  MXKB = MAX(MXKB,NSTFSMSPGP2(KSM,KBTP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
!               DO ISMOB = 1,NSMOB
                DO ISMOB = NSMOB/2 + 1,NSMOB
                  MXSOB = MAX(MXSOB,NTSOB2(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
C               LCJBLK = MXSOB*MXKB*MXB
                LCJBLK = MXSOB*MXKB*MXA
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Ib)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP2(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IATPABS,1,IOBTP,1,K1ATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
!               DO ISMOB = 1, NSMOB
                DO ISMOB = NSMOB/2 + 1,NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP 
*  from K1ATP
                  CALL NEWTYP_REL(K1ATP,1,JOBTP,1,KATP)
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
!                   DO JSMOB = 1, NSMOB
                    DO JSMOB = NSMOB/2 + 1,NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,Kb)
* both Ka and Ib blocked
*
      IF(NTEST.GE.100) WRITE(6,*) ' MXCIJB under development '
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ia,Kb)'
        write(6,*)
      end if
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP 
*  from IBTP
              CALL NEWTYP_REL(IBTPABS,1,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
!               DO ISMOB = 1, NSMOB
                DO ISMOB = NSMOB/2 + 1,NSMOB
                  MXISOB = MAX(MXISOB,NTSOB2(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP 
*  from K1ATP
                  CALL NEWTYP_REL(K1BTP,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP2(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
!                   DO JSMOB = 1, NSMOB
                    DO JSMOB = NSMOB/2 + 1,NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB2(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,Kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Kb)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP 
*  from IATP
              CALL NEWTYP_REL(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
!               DO ISMOB = 1, NSMOB
                DO ISMOB = NSMOB/2 + 1,NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP 
*  from IBTP
                  CALL NEWTYP_REL(IBTPABS,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP2(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
!                   DO JSMOB = 1, NSMOB
                    DO JSMOB = NSMOB/2 + 1,NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB2(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is
* obtained as single annihilations
      MXSXBL = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'alpha strings'
        write(6,*)
      end if
*
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP_REL(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
!           DO JSMOB = 1, NSMOB
            DO JSMOB = NSMOB/2 + 1,NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP_REL(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
!               DO ISMOB = 1, NSMOB
                DO ISMOB = NSMOB/2 + 1,NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'beta strings'
        write(6,*)
      end if
*
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP2(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP_REL(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
!           DO JSMOB = 1, NSMOB
            DO JSMOB = NSMOB/2 + 1,NSMOB
               MXJOB = MAX(MXJOB,NTSOB2(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP_REL(KBTP,2,IOBTP,1,JBTP)
              IF(JBTP.GT.0) THEN
                MXIOB = 0
!               DO ISMOB = 1, NSMOB
                DO ISMOB = NSMOB/2 + 1,NSMOB
                  MXIOB = MAX(MXIOB,NTSOB2(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIB)
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL
      END IF
*
      END

***********************************************************************
      SUBROUTINE NEWTYP_REL(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
*
* an input supergroup is given.
* apply a string of elementary operators to this supergroup and
* obtain supergroup number of new string
*
* Jeppe Olsen, October 1993
* GAS-version : July 95
*
* ------
* Input
* ------
*
* INSPGP  : input super group ( given occupation in each AS )
* IACOP(I) = 1 : operator I is an annihilation operator
*          = 2 : operator I is a  creation   operator
* ITPOP(I) : orbitals space of operator I
* NOP : Number of operators
*
* Output
* ------
* OUTSPGP  : supergroup of resulting string
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      INTEGER ITPOP(*),IACOP(*)
*. Number of active spaces  (NGAS )
#include "mxpdim.inc"
#include "cgas.inc"
#include "strbas.inc"
#include "wrkspc.inc"
*. Local scratch
      DIMENSION IEL(MXPNGAS)
*. output
      INTEGER OUTSPGP
*
      INEW_OR_OLD = 1
      IF(INEW_OR_OLD.EQ.1) THEN
        CALL NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &       NGAS,WORK(KSPGPAN),WORK(KSPGPCR),OUTSPGP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
      ELSE IF(INEW_OR_OLD.EQ.2) THEN
*. Number of electrons in AS1,AS2, .... for input supergroup
        CALL GTSPGP_REL(IEL,INSPGP,2)
*
        IDELTA = 0
        DO IOP = 1, NOP
*. Change in number of orbitals
          IF(IACOP(IOP).EQ.1) THEN
            IDELTA = IDELTA - 1
            IDELTAI = -1
          ELSE
            IDELTA = IDELTA + 1
            IDELTAI = +1
          END IF
          IEL(ITPOP(IOP)) = IEL(ITPOP(IOP)) + IDELTAI
        END DO
*. output supergroup
        CALL GTSPGP_REL(IEL,OUTSPGP,1)
      END IF
*     /\ End  of NEW_OR_OLD switch
*
#if defined LUCI_DEBUG
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NEWTYP ,  OUTSPGP ', OUTSPGP
      END IF
#endif
*
      END
***********************************************************************

      SUBROUTINE GTSPGP_REL(IEL,ISPGP,IWAY)
*
*
* Relation between number of electrons in AS1, AS2 ... and
* supergoup number
*
* IWAY = 1 :
* Get ISPGP : Supergroup of strings
*             with IEL(*)  electrons in each AS
* IWAY = 2 :
* GET IEL(*)  : Number of electrons in each AS for supergroup ISPGP
*
*
* Jeppe Olsen, Another lonely night in Lund
*               GAS version July 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Generel input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. input(IWAY = 2 ), output (IWAY = 1 )
      INTEGER IEL(*)
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Number
        ISPGP = -1
        DO JSPGP = 1, NTSPGP
          IF(ISPGP.EQ.-1) THEN
            IEQUAL = 1
            DO IGAS = 1, NGAS
              IF(NELFSPGP(IGAS,JSPGP).NE.IEL(IGAS))  IEQUAL= 0
            END DO
            IF(IEQUAL.EQ.1) ISPGP = JSPGP
          END IF
        END DO
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Number => Occupation
        DO IGAS = 1, NGAS
         IEL(IGAS) = NELFSPGP(IGAS,ISPGP)
        END DO
      END IF
*
#if defined LUCI_DEBUG
      NTEST  = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from GTSPGP_REL '
        WRITE(6,*)
     &   ' IWAY ISPGP IEL ', IWAY,ISPGP,(IEL(IGAS),IGAS = 1, NGAS)
      END IF
#endif
*
      RETURN
      END
