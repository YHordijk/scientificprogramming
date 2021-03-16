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

      SUBROUTINE DBG_INT_LISTS(NL1,IUB1,NL2,IUB2,ISM2,IDENSI,
     &                         ISPINFREE)
*
* Number of double group integral lists to be constructed
*
*     Jeppe Olsen, Feb. 98
*
* Reordered lists for convenience: We only need the list
* (uu|uu) in an initial or plain SOCI calculation. So this
* is going to be list number 1. Allocation modified as well.
* Modifications in GETINT_DBG incorporated (and elsewhere).
* We are furthermore using real wave functions. This implies
* That complex conjugation symmetry can also be used in the
* quaternion matrix groups C1 and Ci, because their irrep
* basis functions are not generally complex. But if we
* generalize the code to e.g. atomic complex 4- (or 2-)spinors
* or the like, these things have to be distinguished.
* So we will have 3 cases:
*
*   1) Initial or plain SOCI: 1 list
*   2) SOCI + spinor optimization in some way (INS or MCSCF): 3 lists
*   3) Complex wave function: 6 lists
*
* Here, time-reversal symmetry has already been accounted for in
* all cases, so these are non-redundant integral lists.
* At the moment, the first two cases are assumed to be implemented,
* so we will distinguish them by the input keyword DENSI, i.e.
* whether density matrices are to be calculated or not.
*
*     Timo Fleig, June 1999
*
* revised for DIRAC environment
*     Timo Fleig, December 2000
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.output
      INTEGER IUB1(2,10),IUB2(4,10),ISM2(3,10)
*
*. ======================================
*. 1e Lists : all four lists constructed
*. ======================================
*
      NL1 = 4
*. 1 : unbarred, unbarred
      IUB1(1,1) = 1
      IUB1(2,1) = 1
*. 2 : unbarred, barred
      IUB1(1,2) = 1
      IUB1(2,2) = 2
*. 3 : barred, unbarred
      IUB1(1,3) = 2
      IUB1(2,3) = 1
*. 4 : barred, barred
      IUB1(1,4) = 2
      IUB1(2,4) = 2
*
*. =========
*. 2e Lists
*. =========
*
*  For the moment, no complex wave function:
      IWFCOMP = 0
*
* The lists of integrals are
*
* 1 : (u u | u u)
*--------------------
* 2 : (u b | u u)
* 3 : (b u | u b)
* 4 : (b u | b u)
*--------------------
* 5 : (u b | u b)
* 6 : (b u | u u)
*--------------------
*
*=====================
* List 1 ( u u ! u u )
*=====================
*
      IUB2(1,1) = 1
      IUB2(2,1) = 1
      IUB2(3,1) = 1
      IUB2(4,1) = 1
*. No Symmetry between 12 and 34
      ISM2(1,1) = 0
      ISM2(2,1) = 0
      ISM2(3,1) = 0
*
      if (IDENSI.ne.0.or.ISPINFREE.eq.0) then
*
*=====================
* List 2 ( u b ! u u )
*=====================
*
        IUB2(1,2) = 1
        IUB2(2,2) = 2
        IUB2(3,2) = 1
        IUB2(4,2) = 1
*. No Symmetry between 12 and 34
        ISM2(1,2) = 0
        ISM2(2,2) = 0
        ISM2(3,2) = 0
*
*=====================
* List 3 ( b u ! u b )
*=====================
*
        IUB2(1,3) = 2
        IUB2(2,3) = 1
        IUB2(3,3) = 1
        IUB2(4,3) = 2
*. No Symmetry between 12 and 34
        ISM2(1,3) = 0
        ISM2(2,3) = 0
        ISM2(3,3) = 0
*
*=====================
* List 4 ( b u ! b u )
*=====================
*
        IUB2(1,4) = 2
        IUB2(2,4) = 1
        IUB2(3,4) = 2
        IUB2(4,4) = 1
*. No Symmetry between 12 and 34
        ISM2(1,4) = 0
        ISM2(2,4) = 0
        ISM2(3,4) = 0
*
        if (IWFCOMP.eq.1) then
*
*=====================
* List 5 ( u b ! u b )
*=====================
*
          IUB2(1,5) = 1
          IUB2(2,5) = 2
          IUB2(3,5) = 1
          IUB2(4,5) = 2
*. Symmetry between 12 and 34
          ISM2(1,5) = 0
          ISM2(2,5) = 0
          ISM2(3,5) = 0
CTEMP     ISM2(3,5) = 1
*
*=====================
* List 6 ( b u ! u u )
*=====================
*
          IUB2(1,6) = 2
          IUB2(2,6) = 1
          IUB2(3,6) = 1
          IUB2(4,6) = 1
*. No Symmetry between 12 and 34
          ISM2(1,6) = 0
          ISM2(2,6) = 0
          ISM2(3,6) = 0
*
        end if
      end if
*
      if (IDENSI.eq.0.and.ISPINFREE.eq.1) then
        NL2 = 1
      else
        NL2 = 4
        if (IWFCOMP.eq.1) NL2 = 6
      end if
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DETTYP(IPRNT)
*
* Allowed set of determinants : Distribution of
* electrons in alpha(up) -and beta(down) spinors
*
* Defined in terms of
* NMS2VAL : Number of MS2 values (MK2 values )
* MS2VAL  : The allowed MS2 values
*
* Jeppe Olsen, July 97
* Timo Fleig, August '97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
*. Output
#include "cstate.inc"
*
#if defined LUCI_DEBUG
      write(6,*) ' -------------------------------- '
      write(6,*) '   DETTYP speaking  '
      write(6,*) ' -------------------------------- '
      write(6,'(/A17,I4,I4)') 'MK2REF_CI,MK2DEL_CI are',
     &                         MK2REF_CI,MK2DEL_CI
#endif
*
      IF(ISPINFREE.EQ.1) THEN
*. Single MS2 space
        NMS2VAL   = 1
        MS2VAL(1) = MK2REF_CI
      ELSE
*. Potentially several MS2/MK2 values
        MK2MAX = MK2REF_CI + MK2DEL_CI
*. MK2MAX gives highest MK2 value for possibly coupled strings (?)
        if (MK2MAX.gt.NACTEL) then
          MK2MAX = NACTEL
        else
          call evenodd(IEVODMX,MK2MAX)
          call evenodd(IEVODN,NACTEL)
          if (IEVODN.eq.1.and.IEVODMX.eq.2) MK2MAX = MK2MAX - 1
        end if
        MK2MIN = MK2REF_CI - MK2DEL_CI
        IF(MK2MIN.LT.-NACTEL) MK2MIN = -NACTEL
* Pairing ? Ensure symmetry between +MK and - MK spaces
        IPAIR = 1
        call evenodd(IEVOD,NACTEL)
        IF(IPAIR.EQ.1) THEN
          IF(MK2MIN.GT.0) THEN
            MK2MINP = -MK2MIN
          ELSE
            if (IEVOD.eq.1) then
               MK2MIN  =  1
               MK2MINP = -1
            else if (IEVOD.eq.2) then
               MK2MIN  =  0
               MK2MINP = -2
            else
              write(6,*) 'Problem in evenodd'
              call quit(
     &       '*** error in dettyp (ipair==1): evenodd problem.***')
            end if
          END IF
          MK2MAXP = -MK2MAX
        ELSE
          if (IEVOD.eq.1) then
             MK2MINP = -1
             MK2MAXP =  1
          else if (IEVOD.eq.2) then
             MK2MINP = -2
             MK2MAXP =  0
          else
             write(6,*) 'Problem in evenodd'
             call quit(
     &       '*** error in dettyp (ipair==0): evenodd problem.***')
          end if
        END IF
*. And construct the MS2 spaces
        NMS2VAL = 0
        DO MS2 = MK2MIN,MK2MAX,2
          NMS2VAL         = NMS2VAL + 1
          MS2VAL(NMS2VAL) = MS2
        END DO
        DO MS2 = MK2MAXP,MK2MINP,2 !sic
          NMS2VAL         = NMS2VAL+1
          MS2VAL(NMS2VAL) = MS2
        END DO
      END IF
*
#if defined LUCI_DEBUG
      WRITE(6,'(/A,I2/)') ' Number of MS2/MK2 values', NMS2VAL
      WRITE(6,'(A)') ' Included MS2/MK2 values '
      WRITE(6,'(A)') ' ======================= '
      DO IMS2 = 1, NMS2VAL
        WRITE(6,'(I4)') MS2VAL(IMS2)
      END DO
#endif
*
      END
***********************************************************************

      SUBROUTINE GASANA_REL(C,NBLOCK,IBLOCK,IBLTP,LUC,LUSCR,ICISTR,IRC)
*
*
* Analyze CI vector
*
* Jeppe Olsen, August 1995
* Driven By IBLOCK, May 1997
*           String occupations added, Feb. 98
*
*
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*
      use memory_allocator
      use symmetry_setup_krci
      use mospinor_info
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "cands.inc"
#include "parluci.h"
#include "dgroup.h"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
      dimension C(*)
      integer               :: mj_statistics(-32:32)
      integer, allocatable  :: mat1(:,:)
      integer, allocatable  :: astring(:,:)
      integer, allocatable  :: bstring(:,:)
      integer, allocatable  :: occupation_class(:,:)
      real(8), allocatable  :: mat2(:,:)
      character (len=2)     :: ferm_lab
      character (len=1)     :: inv_symbol
*
      CALL QENTER('ANACI')
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRWFC)
*
!     MJ statistics
      do i = -32,32
       mj_statistics(i) = 0
!      write(6,'(a,i3,a,i2)') '   initialize - mj_statistics(',i,') = ',
!    & mj_statistics(i)
      end do
      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call itods(-1,1,-1,LUSCR)
      inv_symbol = ' '
*
* Loop over real and imaginary parts
      do IRIR = 1,IRC,1
        if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
        if (IRIR.EQ.2) then
*
*. Skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if
        do IMK2_R = 1,NMS2VAL,1
          MK2_R = MS2VAL(IMK2_R)
          IACTP = IST_FOR_DT(1,IMK2_R)
          IBCTP = IST_FOR_DT(2,IMK2_R)
*
*  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2(IMK2_R-1)
          end if
*
          call REWINE(LUSCR,-1)
          NBL_C = NBLK_MS2(IMK2_R)
          CALL COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          CALL ITODS(-1,1,-1,LUSCR)
#if defined LUCI_DEBUG
          if (NTEST.ge.50) then
            WRITE(6,*)
     &      ' ********************************************* '
            WRITE(6,*)
     &      ' *                                           * '
            WRITE(6,*)
     &      ' * Start of contribution from given MK_2     * '
            WRITE(6,*)
     &      ' *                                           * '
            WRITE(6,*)
     &      ' ********************************************* '
            WRITE(6,*)
            WRITE(6,'(4X,A,2I3)')
     &      ' MK2/IRI vector of C :  MK2_R, IRIR =', MK2_R,IRIR
            WRITE(6,*)
          end if
#endif
          write(6,'(/A/A,I4,I3/A)')
     &         ' ----------------------------------------------',
     &         ' Contrib. Kramers proj., real/im  ',MK2_R,IRIR,
     &         ' ----------------------------------------------'
          NAEL = NELEC(IACTP)
          NBEL = NELEC(IBCTP)
          call alloc(astring,NAEL,MXNSTR, id='astring')
          call alloc(bstring,NBEL,MXNSTR, id='bstring')

!         determine label for print after the analysis
          if(mod(nael+nbel,2).ne.0)then
            ferm_lab = '/2'
          else
            ferm_lab = '  '
          end if
 
*
          NOCTPA   = NOCTYP(IACTP)
          NOCTPB   = NOCTYP(IBCTP)
*
          IOCTPA   = IBSPGPFTP(IACTP)
          IOCTPB   = IBSPGPFTP(IBCTP)
*
          IASPGPTP = IACTP
          IBSPGPTP = IBCTP
*
#if defined LUCI_DEBUG
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' ================'
            WRITE(6,*) ' GASANA speaking '
            WRITE(6,*) ' ================'
            WRITE(6,'(A26)') ' IACTP   IBCTP   NAEL   NBEL'
            WRITE(6,'(3X,4I6)')  IACTP,IBCTP,NAEL,NBEL
            WRITE(6,'(A14,2I6)') ' NOCTPA NOCTPB ', NOCTPA,NOCTPB
          END IF
#endif
*
**. Info on block structure of space
*
*
*. Number of terms to be printed
          IF (NTEST.GE.100) THEN
            THRES  = 5.0D-10
            MAXTRM = 10000
          ELSE
            THRES  = 1.0D-2
            MAXTRM = 100
          END IF
*
*. Number of occupation classes
          IWAY = 1
          NEL = NAEL + NBEL
          CALL OCCLS_REL(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &                   IGSOCC(1,1),IGSOCC(1,2))
*
*. and then the occupation classes
          call alloc(occupation_class,ngas,noccls,id='occ_class')
          call alloc(mat1,10,NOCCLS, id='mat1')
          call alloc(mat2,10,NOCCLS, id='mat2')
          IWAY = 2
          CALL OCCLS_REL(IWAY,NOCCLS,occupation_class,NEL,NGAS,
     &                   IGSOCC(1,1),IGSOCC(1,2))
*
*. Occupation of strings of given sym and supergroup

          CALL GASANAS_REL(C,LUSCR,NBL_C,WORK(KNSTSO(IACTP)),
     &                     WORK(KNSTSO2(IBCTP)),
     &                     NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &                     NBLOCK,IBLOCK(1,IBOFF),
     &                     THRES,MAXTRM,NAEL,NBEL,
     &                     astring, bstring,
     &                     IBLTP,NSMST,
     &                     mat1, mat2,
     &                     NELFSPGP,
     &                     NOCCLS,NGAS,occupation_class,
     &                     IASPGPTP,IBSPGPTP,
     &                     ICISTR,mj_statistics,inv_symbol,
     &                     IPRWFC)
*

          call dealloc(mat2)
          call dealloc(mat1)
!         NAEL/NBEL may be zero and therefore no allocation of string arrays...
          if(allocated(bstring)) call dealloc(bstring)
          if(allocated(astring)) call dealloc(astring)
          call dealloc(occupation_class)
        end do
      end do

!
!     print statistics
      if(linear)then
        write(6,'(/a))') 
     &  '   *************** MJ statistics ***************'
        do mjval = -32,32
        if(mj_statistics(mjval).gt.0)then
         xpercent = mj_statistics(mjval)*100.0D0/L_COMBI
         if(ferm_lab.eq.'  ')then
         write(6,'(a,i4,a1,a,i10))')  
     &  '   number of dets with MJ-value',
     &   mjval,inv_symbol,': ',mj_statistics(mjval)
         write(6,'(a,f10.2,a))')  
     &  '   percentage of total CI space     :',
     &   xpercent,'%'
         else
         write(6,'(a,i4,a2,a1,a,i8))')  
     &  '   number of dets with MJ-value',
     &   mjval,ferm_lab,inv_symbol,': ',mj_statistics(mjval)
         write(6,'(a,f8.2,a))')  
     &  '   percentage of total CI space       :',
     &   xpercent,'%'
         end if
        end if
        end do
        write(6,'(a/))') 
     &  '   *********************************************'
      end if
 
      CALL QEXIT('ANACI')
*
      END
***********************************************************************

      SUBROUTINE GASANAS_REL(C,LUC,NBL_C,NSSOA,NSSOB,
     &                       NOCTPA,NOCTPB,
     &                       MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                       THRES,MAXTRM,NAEL,NBEL,
     &                       IASTR,IBSTR,IBLTP,NSMST,
     &                       NCPMT,WCPMT,NELFSPGP,NOCCLS,
     &                       NGAS,IOCCLS,
     &                       IASPGPTP,IBSPGPTP,
     &                       ICISTR,mj_statistics,inv_symbol,
     &                       IPRNT)
*
* Analyze CI vector :
*
*      1) Print atmost MAXTRM  combinations with coefficients
*         larger than THRES
*         Currently the corresponding dets are not GIVEN !!
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given
*         occupation class
*
* Jeppe Olsen , Jan. 1989 ,
*               Aug 1995 : GAS version
*               May 1997 : BLOCK driven
*
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*
*
      use mospinor_info
*  Orbitals
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION C(*)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBL_C)
*. Output
      DIMENSION NCPMT(10,NOCCLS)
      DIMENSION WCPMT(10,NOCCLS)
#include "dcborb.h"
#include "dgroup.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
!     scratch
      integer, allocatable              :: mjorb_vec(:)
      integer, intent(inout)            :: mj_statistics(-32:32)
      character (len=1) , intent(inout) :: inv_symbol
      character (len=1)                 :: inv_symbol_orb
      character (len=4),  allocatable   :: conf_string_a(:)
      character (len=4),  allocatable   :: conf_string_b(:)
      character (len=10), allocatable   :: stringa(:)
      character (len=10), allocatable   :: stringb(:)
      logical fndlab, ex, ex_open
*
      NTESTL = 0000
      NTEST = max(IPRNT,NTESTL)
*
      MINPRT    = 0
      ITRM      = 0
      IDET      = 0
      IIDET     = 0
      ILOOP     = 0
      NCIVAR    = 0
      CNORM     = 0.0D0

!     allocate array for mj-values and read info from file
      if(linear)then 
        allocate(conf_string_a(4*100))
        allocate(conf_string_b(4*100))
        allocate(stringa(10*100))
        allocate(stringb(10*100))
        allocate(mjorb_vec(norbt))
        mjorb_vec = 0
        luxxx = -1
        INQUIRE(FILE='KRMCSCF',opened=EX_open,number=luxxx)
        IF( .not. EX_open)THEN
          luxxx = 98
          CALL OPNFIL(luxxx,'KRMCSCF','UNKNOWN','ana mj')
        END IF
        rewind(luxxx)
        call ireakrmc(luxxx,'MJVEC   ',mjorb_vec,norbt)
        IF( .not. EX_open) close(luxxx,status="keep")
      end if
!     do i = 1, nfsym
!       write(6,'(/a,i3)') 
!    &     '   mj-values of orbitals per fermion sym',i
!       do j = 1, norb(i)
!         write(6,'( a,i5,a,i5,a)') '     orbital',j+iorb(i),':',
!    &                                    mjorb_vec(j+iorb(i)),'/2'
!       end do
!     end do

      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
 2001 CONTINUE
      IF( ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      IIDET = 0
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
         XMAX = 1.0D0
         XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
         XMAX = XMIN
         XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNT.EQ.2) THEN
*. Print in one shot
        XMAX     = 3006.1956
        XMIN     = -0.001D0
        ipot_max = 12
      END IF
!     restrict to > 10-6
      ipot_max = 5
      IDET     = 0
*
      WRITE(6,'(/A,E10.4,A,E10.4/A/)')
     & '  Printout of coefficients in interval  ',XMIN,' to ',XMAX,
     & '  ========================================================='
*
#if defined LUCI_DEBUG
      if (NTEST.ge.10) then
         WRITE(6,'(A,I6)') ' GASANAS : NBL_C = ',NBL_C
         WRITE(6,'(A,I6)') ' GASANAS : IASPGPTP = ',IASPGPTP
         WRITE(6,'(A,I6)') ' GASANAS : IBSPGPTP = ',IBSPGPTP
      end if
#endif
*
      DO JBLOCK = 1, NBL_C
*
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
#if defined LUCI_DEBUG
         if (NTEST.ge.5) write(6,'(A,1X,I6,I6,I6,I6)')
     &        ' IATP,IBTP,IASM,IBSM = ',IATP,IBTP,IASM,IBSM
#endif
*. Obtain alpha strings of sym IASM and type IATP
         CALL GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,
     &                              NASTR1,IASTR)
*. Obtain beta strings of sym IBSM and type IBTP
         CALL GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,
     &                              NBSTR1,IBSTR)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR_ANA = 1
         ELSE
            IRESTR_ANA = 0
         END IF
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
#if defined LUCI_DEBUG
         if (NTEST.ge.5) write(6,'(A,1X,I9,I9)')
     &        ' NIA,NIB = ',NIA,NIB
#endif
*
         IMZERO = 0
         if (ICISTR.ge.2) then
*. Read in a Type-Type-symmetry block
*? Perhaps use GSTTBL here?
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
#if defined LUCI_DEBUG
            if (NTEST.ge.10) then
              write(6,'(A,I8)') ' Number of elements read in ',IDET
            end if
            if (NTEST.ge.1000) then
              write(6,*)
              write(6,*) 'Partition of C vector is: '
              call wrtmat(C,IDET,1,IDET,1)
            end if
#endif
            IDET = 0
         END IF
C?       write(6,*) ' IMZERO in gasanas:  ',IMZERO
!        if(linear) IMZERO = 0
         IF(IMZERO.NE.1) THEN
*
            IBBAS = 1
            IABAS = 1
*
            DO IB = IBBAS,IBBAS+NIB-1
               IF(IRESTR_ANA.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO  IA = MINIA,IABAS+NIA-1
*
                  IF (ILOOP .EQ. 1 ) NCIVAR = NCIVAR + 1
                  IDET = IDET + 1
                  IF (XMAX .GE. ABS(C(IDET)) .AND.
     &                 ABS(C(IDET)).GT. XMIN ) THEN
                    if(linear)then
!                     calculate MJ value for this determinant 
!                     MJ = mj_orb1 (x) mj_orb2 (x) mj_orb3 (x) ... 
!                     = mj_orb1 + mj_orb2 + mj_orb3 + ...
                      mjval = 0
                      inversion_counter = 0
                      do iel = 1, nael
                        nxtorb = iastr(iel,ia)
!                       write(6,*) 'unbarred: nxtorb is',nxtorb
                        nxtorb = 
     &                  imosp_luci2dirac1(ireots (iastr(iel,ia)))
!                      write(6,*) 'unbarred: reordered nxtorb is',nxtorb
                        iorb_tmp = 0
                        do j = 1, nfsym
                          do i = 1, ngas
                            iorb_idx_gas   = 0
                            do ij = 1, ngsh(j,i)
                              iorb_tmp     = iorb_tmp + 1
                              iorb_idx_gas = iorb_idx_gas + 1
                              if(iorb_tmp.eq.nxtorb) then
                                is_ferm_irrep  = j
                                is_gas_shell   = i
                                is_orb_idx_gas = iorb_idx_gas
                                if(nfsym.eq.1)then
                                  inv_symbol_orb = ' '
                                else
                                  inv_symbol_orb = 'g'
                                  inversion_counter = 
     &                            inversion_counter + 2
                                  if(j.eq.2)then 
                                    inv_symbol_orb = 'u'
                                    inversion_counter = 
     &                              inversion_counter + 1
                                  end if
                                end if
                                exit
                              end if
                            end do
                            if(iorb_tmp.eq.nxtorb) exit
                          end do
                            if(iorb_tmp.eq.nxtorb) exit
                        end do
!                       calculate dirac orbital index
                        ioff_idx = 0
                        if(is_ferm_irrep.gt.1) 
     &                    ioff_idx  = ioff_idx + nash(1)
                          itmp_val = 0
                          do i = 1, is_gas_shell-1
                            itmp_val = 
     &                      itmp_val + ngsh(is_ferm_irrep,i)
                          end do
                          iorbidx_dc_act = 
     &                    ioff_idx + itmp_val + is_orb_idx_gas
                          iorbidx_dc     = IDXU2G(iorbidx_dc_act)
!                         write(6,*) 'alph: iorbidx_dc is',iorbidx_dc
!                         get orbital mj-value
                          mjval_orb = mjorb_vec(iorbidx_dc)
                          write(conf_string_a(iel),'(i3,a1)') 
     &                    mjval_orb,inv_symbol_orb
                          mjval     = mjval + mjval_orb
                      end do
                      do iel = 1, nbel
                        nxtorb = ibstr(iel,ib)
!                       write(6,*) 'barred: nxtorb is',nxtorb
                        nxtorb = 
     &                  imosp_luci2dirac2(ireots2(ibstr(iel,ib)))
!                       write(6,*) 'barred: reordered nxtorb is',nxtorb
                        iorb_tmp = 0
                        do j = 1, nfsym
                          do i = 1, ngas
                            iorb_idx_gas = 0
                            do ij = 1, ngsh(j,i)
                              iorb_tmp     = iorb_tmp + 1
                              iorb_idx_gas = iorb_idx_gas + 1
                              if(iorb_tmp.eq.nxtorb) then
                                is_ferm_irrep  = j
                                is_gas_shell   = i
                                is_orb_idx_gas = iorb_idx_gas
                                if(nfsym.eq.1)then
                                  inv_symbol_orb = ' '
                                else
                                  inv_symbol_orb = 'g'
                                  inversion_counter = 
     &                            inversion_counter + 2
                                  if(j.eq.2)then 
                                    inv_symbol_orb = 'u'
                                    inversion_counter = 
     &                              inversion_counter + 1
                                  end if
                                end if
                                exit
                              end if
                            end do
                            if(iorb_tmp.eq.nxtorb) exit
                          end do
                            if(iorb_tmp.eq.nxtorb) exit
                        end do
!                       calculate dirac orbital index
                        ioff_idx = 0
                        if(is_ferm_irrep.gt.1) 
     &                    ioff_idx  = ioff_idx + nash(1)
                          itmp_val  = 0
                          do i = 1, is_gas_shell-1
                            itmp_val = 
     &                      itmp_val + ngsh(is_ferm_irrep,i)
                          end do
                          iorbidx_dc_act = 
     &                    ioff_idx + itmp_val + is_orb_idx_gas
                          iorbidx_dc     = IDXU2G(iorbidx_dc_act)
!                         write(6,*) 'beta: iorbidx_dc is',iorbidx_dc
!                         get orbital mj-value
                          mjval_orb = -mjorb_vec(iorbidx_dc)
                          write(conf_string_b(iel),'(i3,a1)') 
     &                    mjval_orb,inv_symbol_orb
                          mjval     =  mjval + mjval_orb
                      end do
                      if(mod(nael+nbel,2).eq.0)then
                        mj_statistics(mjval/2) =
     &                  mj_statistics(mjval/2) + 1
                      else
                        mj_statistics(mjval/1) =
     &                  mj_statistics(mjval/1) + 1
                      end if
                    end if
                    ITRM = ITRM + 1
                    IIDET = IIDET + 1
                    IF( ITRM .LE. MAXTRM ) THEN
                       CNORM = CNORM + C(IDET) ** 2
                       WRITE(6,'(A,I8,A,E14.8/A)')
     &                  '  Coefficient of combination ',IDET,' is ',
     &                      C(IDET),
     &                  '  Corresponding alpha- and beta-strings'
                       if(linear)then

                         do IEL = 1, NAEL
                           write(stringa(iel),'(I4,a1,a4,a1)') 
     &                     imosp_luci2dirac1(ireots(iastr(iel,ia))),
     &                     '(',conf_string_a(iel),')'
                         end do
                         WRITE(6,'("  alpha",20(a))')
     &                       (stringa(iel), IEL = 1, NAEL )


                         do IEL = 1, NbEL
                           write(stringb(iel),'(I4,a1,a4,a1)') 
     &                     imosp_luci2dirac2(ireots2(ibstr(iel,ib))),
     &                     '(',conf_string_b(iel),')'
                         end do
                         WRITE(6,'("   beta",20(a))')
     &                       (stringb(iel), IEL = 1, NBEL )
                         inv_symbol = ' '
                         if(nfsym.gt.1)then
                           if(mod(inversion_counter,2).eq.0)then
                             inv_symbol = 'g'
                           else
                             inv_symbol = 'u'
                           end if
                         end if
                         if(mod(nael+nbel,2).eq.0)then
                           write(6,'("  MJ =",i5,a1)')  mjval/2,
     &                                                       inv_symbol
                         else 
                           write(6,'("  MJ =",i5,"/2",a1)')  mjval,
     &                                                       inv_symbol
                         end if
                       else
                       WRITE(6,'("  alpha",20I4)')
     &                      (imosp_luci2dirac1(ireots (iastr(iel,ia))),
     &                       IEL = 1, NAEL )


                       WRITE(6,'("   beta",20I4)')
     &                      (imosp_luci2dirac2(ireots2(ibstr(iel,ib))),
     &                       IEL = 1, NBEL )
                       end if
                    END IF
                  END IF
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of if statement for nonvanishing blocks
      END DO
*     ^ End of loop over blocks
      IF(IIDET .EQ. 0 ) WRITE(6,*) '   ( no coefficients )'
      IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
*
      WRITE(6,'(/A,E15.8/)')
     &     '  Norm of printed CI vector .. ', CNORM


*
*.Size of CI coefficients
*
*
      IDET = 0
      IF(ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      CALL IZERO(NCPMT,10*NOCCLS)
      CALL DZERO(WCPMT,10*NOCCLS)
      DO JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)
*
         IF(IBLTP(IASM).EQ.2) THEN
            IRESTR_ANA = 1
         ELSE
            IRESTR_ANA = 0
         END IF
*. Occupation class corresponding to given occupation
         JOCCLS = 0
         DO JJOCCLS = 1, NOCCLS
            IM_THE_ONE = 1
            DO IGAS = 1, NGAS
               IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &          NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &          IM_THE_ONE = 0
            END DO
            IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
         END DO
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
*
         IMZERO = 0
         IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
            CALL IFRMDS(IDET,1,-1,LUC)
            CALL FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            IDET = 0
         END IF
         IF(IMZERO.EQ.0) THEN

            DO IB = IBBAS,IBBAS+NIB-1
               IF(IRESTR_ANA.EQ.1.AND.IATP.EQ.IBTP) THEN
                  MINIA = IB - IBBAS + IABAS
               ELSE
                  MINIA = IABAS
               END IF
               DO IA = MINIA,IABAS+NIA-1
                  IDET = IDET + 1
                  DO IPOT = 1, ipot_max
                     IF(10.0D0 ** (-IPOT+1).GE.ABS(C(IDET)).AND.
     &                 ABS(C(IDET)).GT. 10.0D0 ** ( - IPOT )) THEN
                       NCPMT(IPOT,JOCCLS)=
     &                       NCPMT(IPOT,JOCCLS)+ 1
                       WCPMT(IPOT,JOCCLS)=
     &                       WCPMT(IPOT,JOCCLS)+ C(IDET) ** 2
                     END IF
                  END DO
*     ^ End of loop over powers of ten
               END DO
*     ^ End of loop over alpha strings
            END DO
*     ^ End of loop over beta strings
         END IF
*     ^ End of test for novanishing blocks
      END DO
*     ^ End of loop over blocks
*
      WRITE(6,'(A)') '   Magnitude of CI coefficients '
      WRITE(6,'(A)') '  =============================='
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, ipot_max
         W = 0.0D0
         N = 0
         DO 290 JOCCLS = 1, NOCCLS
            N = N + NCPMT(IPOT,JOCCLS)
            W = W + WCPMT(IPOT,JOCCLS)
 290     CONTINUE
         WACC = WACC + W
         NACC = NACC + N
         WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &        '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
 300  CONTINUE
*
      WRITE(6,'(a,i2,a,i10)') ' Number of coefficients less than  10-',
     &     ipot_max+1,' IS  ',NCIVAR - NACC
*
      IF(NOCCLS.NE.1) THEN
         WRITE(6,'(A)')
         WRITE(6,'(A)')
     &    '   Magnitude of CI coefficients for each excitation level '
         WRITE(6,'(A)')
     &    '  ========================================================='
         DO 400 JOCCLS = 1, NOCCLS
            N = 0
            DO 380 IPOT = 1, ipot_max
               N = N + NCPMT(IPOT,JOCCLS)
 380        CONTINUE
            IF(N .NE. 0 ) THEN
               WRITE(6,*)
               WRITE(6,'(A,15I3)')'       Occupation of active sets :',
     &              (IOCCLS(IGAS,JOCCLS),IGAS=1, NGAS)
               WRITE(6,'(A,I2,a,i9)')
     &          '         Number of coefficients larger than 10-',
     &                    ipot_max+1,':', N
               WACC = 0.0D0
               DO 370 IPOT = 1, ipot_max
                  N =  NCPMT(IPOT,JOCCLS)
                  W =  WCPMT(IPOT,JOCCLS)
                  WACC = WACC + W
                  WRITE(6,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &                 '  10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
 370           CONTINUE
            END IF
 400     CONTINUE
*
*. Total weight and number of dets per excitation level
*
         WRITE(6,'(/a,i2,a)')
     &        '   Total weight and number of SD''s (> 10 ** -,',
     &        ipot_max+1,')  : '
         WRITE(6,'(A)')
     &        '  ================================================='
       WRITE(6,*) '        N      Weight      Acc. Weight   Occupation '
       WRITE(6,*) ' ==================================================='
         WACC = 0.0D0
         DO 500 JOCCLS = 1, NOCCLS
            N = 0
            W = 0.0D0
            DO 480 IPOT = 1, ipot_max
               N = N + NCPMT(IPOT,JOCCLS)
               W = W + WCPMT(IPOT,JOCCLS)
 480        CONTINUE
            WACC = WACC + W
            IF(N .NE. 0 ) THEN
               WRITE(6,'(1X,I9,3X,E9.4,7X,E9.4,2X,16(1X,I2))')
     &              N,W,WACC,(IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
            END IF
 500     CONTINUE
      END IF

!     release scratch memory
      if(linear)then 
        deallocate(mjorb_vec)
        deallocate(stringb)
        deallocate(stringa)
        deallocate(conf_string_b)
        deallocate(conf_string_a)
      end if

      END
***********************************************************************

      SUBROUTINE GET_3BLKS_REL(KVEC1,KVEC2,LSCR2,CIRUN)
C**********************************************************************
C
C     allocate two blocks: VEC1 and VEC2
C     dimension for resolution block is computed:
C                                    - two partitioned blocks of
C                                      length LSCR2 can be stored
C
C      Jeppe Olsen, Jan 1998 - Relativistic version, Feb. 98
C
C      S. Knecht - Jul 2007: memory saving version enabled
C      S. Knecht - Nov 2008: modified for CI property runs
C
C**********************************************************************
C
      use symmetry_setup_krci
      use mospinor_info
      IMPLICIT REAL*8(A-H,O-Z)
C
      character*6 CIRUN
C
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
#endif
#include "parluci.h"
#include "krciprop.h"
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
C
#include "cgas.inc"
#include "cands.inc"
C
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',KDUM,'GET_3B')
C
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRDIA)
C
      ISPC = ISSPC
      LSCR2 = 0
C
      IF( CIRUN .ne. 'PROP1 ' ) IXSYMACT = 1
C      
      DO I = 1, IXSYMACT
C
        IJSYM = ISSM
        IF( CIRUN .eq. 'PROP1 ' ) IJSYM = IXSYMLIST(I)
C
C       loop over MS2 spaces for "active" symmetry irrep
        DO IMS2 = 1, NMS2VAL
C
          IATP = IST_FOR_DT(1,IMS2)
          IBTP = IST_FOR_DT(2,IMS2)
C
          NAEL = NELEC(IATP)
          NBEL = NELEC(IBTP)
C
          NOCTPA = NOCTYP(IATP)
          NOCTPB = NOCTYP(IBTP)
C
          IOCTPA = IBSPGPFTP(IATP)
          IOCTPB = IBSPGPFTP(IBTP)
C
          call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
C
          CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C         largest block of strings in zero order space
          MXSTBL0 = MXNSTR
C         this must be generalized ! alpha and beta strings with an 
C         electron removed
          IATPM1 = ITPM1(1,IMS2)
          IBTPM1 = ITPM1(2,IMS2)
C         alpha and beta strings with two electrons removed
          IATPM2 = ITPM2(1,IMS2)
          IBTPM2 = ITPM2(2,IMS2)
C         largest number of strings of given symmetry and type
          MAXA = MXNSTR
          IF(NAEL.GE.1) THEN
            MAXA1 = IMNMX_REL(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),
     &                        2)
            MAXA = MAX(MAXA,MAXA1)
          END IF
          IF(NAEL.GE.2) THEN
            MAXA1 = IMNMX_REL(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),
     &                        2)
            MAXA = MAX(MAXA,MAXA1)
          END IF
          MAXB = 0
          IF(NBEL.GE.1) THEN
            MAXB1 =IMNMX_REL(WORK(KNSTSO2(IBTPM1)),NSMST*NOCTYP(IBTPM1),
     &                       2)
            MAXB = MAX(MAXB,MAXB1)
          END IF
          IF(NBEL.GE.2) THEN
            MAXB1 =IMNMX_REL(WORK(KNSTSO2(IBTPM2)),NSMST*NOCTYP(IBTPM2),
     &                       2)
            MAXB = MAX(MAXB,MAXB1)
          END IF
          MXSTBL = MAX(MAXA,MAXB)
C
C         largest number of resolution strings and spectator strings
C         that can be treated simultaneously
          MAXK = MIN(MXINKA,MXSTBL)
C
C         scratch space for CJKAIB resolution matrices
C         size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
          CALL MXRESC_REL(WORK(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                    NSMST,NSTFSMSPGP,NSTFSMSPGP2,MXNDGIRR,NSMOB,
     &                    MXPNGAS,NGAS,NOBPTS,NOBPTS2,IPRCIX,MAXK,
     &                    NELFSPGP,MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,
     &                    MXADKBLK)
          LSCR2 = MAX(LSCR2,MXCJ,MXCIJA,MXCIJB)
        END DO
C       ^ end loop over MK2 values
      END DO
C     ^ end loop over active symmetry irreps
C       
C     L2BLOCK and LBLOCK are stored in crun.inc; dimensioning in Z_BLK_*
      L0BLOCK = MXSOOB
CSK...6-jul-07      LBLOCK = MAX(L2BLOCK,L0BLOCK)
      L3BLOCK = MAX(L2BLOCK,L0BLOCK)
CSK...6-jul-07 new
      L3BLOCK = MIN(LBLOCK,L3BLOCK)
CSK...6-jul-07      LSCR12 = MAX(LBLOCK,2*LSCR2)
!u2T  LSCR12 = MAX(LBLOCK,2*LSCR2)
!     LSCR12 = MAX(LBLOCK,2*LSCR2,lcsblk)
!     this ought to be fine ; no lcsblk
      LSCR12 = MAX(LBLOCK,2*LSCR2)
      LCJRES_SAVE = LSCR12
C     
#if defined LUCI_DEBUG
      write(luwrt,*) 'L3BLOCK is',L3BLOCK
      write(luwrt,*) 'LBLOCK is',LBLOCK
      write(luwrt,*) 'LSCR2 is',LSCR2
      write(luwrt,*) 'LSCR12 is',LSCR12
      write(luwrt,*) 'LCJRES_SAVE is',LCJRES_SAVE
#endif
C
      IDUM = 0
      call memmar(KDUM,IDUM,'FLUSM ',KDUM,'GET_3B')
C
      IF(CIRUN /= 'REFVEC')THEN
      WRITE(LUWRT,'(/3X,A)') 
     & '=============================================================='
      WRITE(LUWRT,'(3X,A)') 
     & '==> allocation of two CI vectors and one resolution vector <=='
      write(LUWRT,*) ' '
      WRITE(LUWRT,'(3X,A,1X,I15)')
     & 'current available free memory in double words:',LMEMFREE_PTR
      WRITE(LUWRT,'(3X,A,1X,I15)')
     & 'allocate two CI vectors each of length:       ',LBLOCK
      WRITE(LUWRT,'(3X,A,1X,I15)')
     & 'allocate resolution vector of length:         ',LCJRES_SAVE
      WRITE(LUWRT,'(3X,A/)') 
     & '=============================================================='
      END IF
C
C     static allocation
      call memmar(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      call memmar(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
C
      END
***********************************************************************

      SUBROUTINE IAIBCM_REL(ICISPC,IATP_TF,IBTP_TF,IAIB)
*
* obtain allowed combination of alpha- and beta- supergroups
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*
*  Corrected for rel. GASCI
*       Timo Fleig, 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "stinf.inc"
#include "strinp.inc"
*. Output
      INTEGER IAIB(*)
*
      NAEL = NELEC(IATP_TF)
      NBEL = NELEC(IBTP_TF)
*
      NOCTPA = NOCTYP(IATP_TF)
      NOCTPB = NOCTYP(IBTP_TF)
*
      IOCTPA = IBSPGPFTP(IATP_TF)
      IOCTPB = IBSPGPFTP(IBTP_TF)
*
C?    write(6,*) ' IAIB ::::::'
C?    write(6,*) ' LCMBSPC, ICISPC, ICMBSPC '
C?    WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
C?    WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),
     &                IGSOCCX,NOCTPA,
     &                NOCTPB,ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                NELFGP,MXPNGAS,NGAS,IAIB,IPRDIA)
*
      RETURN
      END
***********************************************************************

      SUBROUTINE INTDIM_REL(IPRNT)
*
* Number of integrals and storage mode
*
*. Last modifications : Jan 98, relativistic integrals added
*
      use symmetry_setup_krci
      use mospinor_info
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
#include "parluci.h"
*.CSMPRD
*. Local arrays
      integer, allocatable   :: nnorb(:,:)
*
* =======
*. Output
* =======
*. Atmost 10 double group integral lists
C     COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,
C    &              NL1D,IUB1D(2,10),L1D(10),
C    &              NL2D,IUB2D(4,10),ISM2D(3,10),L2D(10)
#include "cintfo.inc"


      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)

      allocate(nnorb(nirr_dg,4))
*
* =======================
* Double group integrals
* =======================
*
*. Number of double group integral lists
      CALL DBG_INT_LISTS(NL1D,IUB1D,NL2D,IUB2D,ISM2D,IDENSI,
     &                   ISPINFREE)
!     print *, ' in INTDIM: NL1D, NL2D',NL1D, NL2D
      DO IL1 = 1, NL1D
        DO INDEX = 1, 2
          IF(IUB1D(INDEX,IL1).EQ.1) THEN
            CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
          ELSE
            CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
          END IF
        END DO
        L1D(IL1) = NSXFSM2(NNORB(1,1),NNORB(1,2),ITSSX,1,IPRNT)
      END DO
*
      DO IL2 = 1, NL2D
        DO INDEX = 1, 4
          IF(IUB2D(INDEX,IL2).EQ.1) THEN
            CALL ICOPVE(NTOOBS,NNORB(1,INDEX),NSMOB)
          ELSE
            CALL ICOPVE(NTOOBS2,NNORB(1,INDEX),NSMOB)
          END IF
        END DO
        L2D(IL2) = NDXFSM2(
     &             NNORB(1,1),NNORB(1,2),NNORB(1,3),NNORB(1,4),
     &             ISM2D(1,IL2),ISM2D(2,IL2),ISM2D(3,IL2),
     &             IPRNT )
      END DO
*. Largest block of two-electron integrals
      LMAXE2LST = NINT2
      DO IL2 = 1, NL2D
        LMAXE2LST = MAX(LMAXE2LST,L2D(IL2))
      END DO
*
#ifdef LUCI_DEBUG
      write(6,'(/A,I10)')
     & '   The largest block of two-electron integrals is',LMAXE2LST
#endif

!     release scratch memory
      deallocate(nnorb)

      END
***********************************************************************

      integer function nsxfsm2(no1ps,no2ps,isxsm,isym,iprnt)

      use symmetry_setup_krci
      implicit none
*
*     number of single excitations of symmetry ISXSM
*
* ISYM = 0 : All symmetry allowed excitations
* ISYM = 1 : Only excitations a+iaj with I.ge.J   ! this is the hardwired default.
* ISYM =-1 : Only excitations a+iaj with I.gt.J
*
* Compared to NSXFSM : Reference to symmetries with call to SYMCOM
*
      integer, intent(in) :: no1ps(*), no2ps(*)
      integer, intent(in) :: isxsm, isym, iprnt
      integer             :: msxfsm, io1sm, io2sm, io2sma
*
      MSXFSM = 0
!     DO IO1SM = 1,NSMOB
      DO IO1SM = NSMOB/2+1,NSMOB

        CALL SYMCOM_REL(2,IO1SM,IO2SMA,ISXSM)
!       IO2SMA is symmetry of annihilation operator, find symmetry of creation op
        IO2SM = iadjsym(IO2SMA)

        IF(IO1SM.EQ.IO2SM)THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)+1)/2
        END IF

      end do
*
      NSXFSM2 = MSXFSM

#if defined LUCI_DEBUG
      WRITE(6,'(A32,I8,A1,I8)')
     &' Number of single excitations of symmetry ',
     &  ISXSM,',',NSXFSM2
#endif
*
      END
***********************************************************************

      integer function ndxfsm2(no1ps,no2ps,no3ps,no4ps,
     &                         is12,is34,is1234,iprnt)
*
*     purpose: compute number of double excitations 
!              with total symmetry IDXSM
*
*     IS12 (0,1,-1)   : Permutational symmetry between index 1 and 2
*     IS34 (0,1,-1)   : Permutational symmetry between index 3 and 3
*     IS1234 (0,1,-1) : permutational symmetry between index 12 and 34
*
      use symmetry_setup_krci
      implicit none
      integer, intent(in) :: no1ps(*),no2ps(*),no3ps(*),no4ps(*)
      integer, intent(in) :: is12, is34, is1234, iprnt
      integer             :: mdx, idxsm
      integer             :: i12sm, i1sm, i2sm, i2sma
      integer             :: i34sm, i3sm, i4sm, i4sma
      integer             :: i12num, n12
      integer             :: i34num, n34
      integer             :: i

#if defined LUCI_DEBUG
      print *,'printing NO1PS'
      write(6,'(2x,8i4)')(NO1PS(i),i=1,NSMOB)
      print *,'printing NO2PS'
      write(6,'(2x,8i4)')(NO2PS(i),i=1,NSMOB)
      print *,'printing NO3PS'
      write(6,'(2x,8i4)')(NO3PS(i),i=1,NSMOB)
      print *,'printing NO4PS'
      write(6,'(2x,8i4)')(NO4PS(i),i=1,NSMOB)
 
      print *,'IS12, IS34, IS1234 and NSMOB',
     &         IS12, IS34, IS1234, NSMOB
#endif

      MDX   = 0
      IDXSM = 1

!     do i1sm = 1, nsmob
      do i1sm = nsmob/2+1, nsmob
        if(NO1PS(I1SM).gt.0)then
!         do i12sm = 1, nsmob
          do i12sm = 1, nsmob

            CALL SYMCOM_REL(2,I1SM,I2SMA,I12SM)
            I2SM   = iadjsym(I2SMA)
            N12    = NO1PS(I1SM)*NO2PS(I2SM)

            if(n12.gt.0)then
              CALL SYMCOM_REL(2,I12SM,I34SM,IDXSM)
!             do i3sm = 1, nsmob
              do i3sm = nsmob/2+1, nsmob

                CALL SYMCOM_REL(2,I3SM,I4SMA,I34SM)
                I4SM   = iadjsym(I4SMA)
                N34    = NO3PS(I3SM)*NO4PS(I4SM)

                MDX    = MDX + N12 * N34
C?              WRITE(6,*) ' I1SM I2SM I3SM I4SM MDX '
C?              WRITE(6,*)   I1SM,I2SM,I3SM,I4SM,MDX
              end do
            end if
          end do
        end if
      end do

      NDXFSM2 = MDX

#if defined LUCI_DEBUG
      WRITE(6,'(A,I10)')
     &' Number of double excitations obtained ', NDXFSM2
#endif
*
      END
***********************************************************************

      SUBROUTINE LCISPC_REL(LCSBLK,IPRNT,CIRUN)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*               July 1997, multiple Ms spaces allowed
*
* GAS VERSION
*
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* ===================
*.Input common blocks
* ===================
*
#include "mxpdim.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "strbas.inc"
*./CSM/ : NSMST is used
#include "stinf.inc"
*./CSTATE/ : IDC is used
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.inc"
#include "parluci.h"
      character (len=6), intent(in) :: cirun
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRNT)
*
      CALL QENTER('LCISP')
*. Number of spaces
      NICISP = NCMBSPC
*
*.    set mark for local memory
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'LCISP ')
*
*.  Largest NOCTPA*NOCTPB
      MXNOCAB = 0
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
        NOCAB = NOCTYP(IATP)*NOCTYP(IBTP)
        MXNOCAB = MAX(MXNOCAB,NOCAB)
      END DO
C     WRITE(6,*) ' MXNOCAB = ',MXNOCAB
      call memmar(KLIOIO,MXNOCAB,   'ADDL  ',2,'KLIOIO')
*
      call memmar(KLBLTP,NSMST,'ADDL  ',2,'KLBLTP')
*. Array defining symmetry combinations of internal strings
      CALL SMOST_REL(NSMST,NSMCI,MXNDGIRR,ISMOST)

!     create batches
      ITTSS_ORD = 2
      MXSB      = 0
      MXSOOB    = 0
      MXNTTS    = 0
      call dzero(XISPSM,MXNDGIRR*MXPICI)
!     calculate number of internal dets for each symmetry
!     ---------------------------------------------------

!     loop over CI spaces
      DO ICI = 1, NICISP
!     loop over symmetry irreps
        DO ISYM = 1, NSMCI

          CALL ZBLTP_REL(ISMOST(1,ISYM),NSMST,WORK(KLBLTP))
          IOFFBLK = 0

!         loop over MS2 spaces
          DO IMS2 = 1, NMS2VAL
            IATP   = IST_FOR_DT(1,IMS2)
            IBTP   = IST_FOR_DT(2,IMS2)
!
            NOCTPA =  NOCTYP(IATP)
            NOCTPB =  NOCTYP(IBTP)
!
            IOCTPA = IBSPGPFTP(IATP)
            IOCTPB = IBSPGPFTP(IBTP)

!           allowed combination of types
            CALL IAIBCM_GAS(LCMBSPC(ICI),ICMBSPC(1,ICI),
     &                      IGSOCCX,NOCTPA,NOCTPB,ISPGPFTP(1,IOCTPA),
     &                      ISPGPFTP(1,IOCTPB),NELFGP,MXPNGAS,NGAS,
     &                      WORK(KLIOIO),IPRNT)

!           a. number of determinants (ci space / symmetry / ms2 value)
            CALL NGASDT_REL(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),
     &                      NGAS,ISYM,NSMST,NOCTPA,NOCTPB,
     &                      WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &                      ISPGPFTP(1,IBSPGPFTP(IATP)),
     &                      ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &                      NCOMB,XNCOMB,MXS,MXSOO,WORK(KLBLTP),NTTSBL,
     &                      LCOL,WORK(KLIOIO))

!           b. number of ttss blocks per (ci space / symmetry / ms2 value)
            call determine_ttss_block_num(WORK(KLBLTP),
     &                                    WORK(KNSTSO(IATP)),
     &                                    WORK(KNSTSO2(IBTP)),NOCTPA,
     &                                    NOCTPB,IOFFBLK,
     &                                    NSMST,WORK(KLIOIO),
     &                                    ISMOST(1,ISYM),
     &                                    ITTSS_ORD)
!
            XISPSM(ISYM,ICI) = XISPSM(ISYM,ICI) + XNCOMB
            MXSOOB           = MAX(MXSOOB,MXSOO)
            MXSB             = MAX(MXSB,MXS)
            MXNTTS           = MAX(MXNTTS,IOFFBLK)
          END DO ! loop over MS2 spaces
        END DO ! loop over symmetries
      END DO ! loop over CI spaces
!
!     Setting size of resolution matrices
      LCSBLK = IMAXLBLKSZ
!      
!     Check for minimum size of resolution matrices:
      if(MXSB.gt.LCSBLK)then
        write(LUWRT,'(/a,i12,a,i12/)') 
     & ' *** warning: resetting LCSBLK from',LCSBLK,' to:',MXSB
        LCSBLK = MXSB
      end if
!
      if(cirun.eq.'NDET  ')then
        write(luwrt,'(/a)')
     &'  Number of determinants per double group irrep '
        write(luwrt,'(a)')
     &'  ============================================= '
        do ici = 1, ncmbspc
          write(luwrt,'(1p,5i20)') (nint(xispsm(ii,ici)),ii=1,nsmci)
        end do
      end if

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(LUWRT,'(/a,i8)') 
     &' Largest symmetry block                          ==>',
     &  MXSB
      WRITE(LUWRT,'(a,i8)')  
     &' Largest type-type-symmetry (TTS) block          ==>',
     &  MXSOOB
      WRITE(LUWRT,'(a,i8/)')  
     &' Largest number of TTS subblocks in CI expansion ==>',
     &  MXNTTS
#endif
!#undef LUCI_DEBUG
!
!     eliminate local memory
      IDUM = 0
      call memmar(KDUM ,IDUM,'FLUSM ',2,'LCISP ')

      CALL QEXIT('LCISP')
      END
***********************************************************************

      SUBROUTINE NGASDT_REL(IOCCMN,IOCCMX,NGAS,ITOTSM,
     &                      NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,
     &                      IAOCC,IBOCC,MXPNGAS,NELFGP,
     &                      NCOMB,XNCOMB,MXSB,MXSOOB,
     &                      IBLTP,NTTSBL,LCOL,IOCOC)
*
* Number of combinations (NCOMB) with symmetry ITOTSM and
* occupation between IOCCMN and IOCCMX
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*),NELFGP(*)
*. block types
      DIMENSION IBLTP(*)
*
#if defined LUCI_DEBUG
      NTEST = 010
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' NSMST = ',NSMST        
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
*
      MXSB   = 0
      MXSOOB = 0
      NCOMB  = 0
      NTTSBL = 0
      LCOL   = 0
      XNCOMB = 0.0d0

      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        CALL SYMCOM_REL(2,IASM,IBSM,ITOTSM)
!       write(6,*) ' Iasm ibsm ', iasm,ibsm
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
             MXBTP = IATP
           ELSE
             MXBTP = NOCTPB
           END IF
           LTSSA = NSSOA(IASM,IATP)
           XLTSSA = LTSSA
!          if(XLTSSA .gt. 0.0d0) print *, ' non-zero alpha part for st',
!    &     IASM, IATP, NSSOA(IASM,IATP)
           if(XLTSSA .gt. 0.0d0) then
             DO 100 IBTP = 1, MXBTP
*
#if defined LUCI_DEBUG
               IF(NTEST.GE.10) THEN
                 WRITE(6,*) ' Alpha super group and beta super group'
                 CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
                 CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
               END IF
#endif
*
               IF(IOCOC(IATP,IBTP).EQ.1) THEN
*. Size of unpacked block
                 LTTSUP =  LTSSA*NSSOB(IBSM,IBTP)
*. Size of packed block
                 IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                   LTTSBL = LTSSA*NSSOB(IBSM,IBTP)
                   XNCOMB = XNCOMB + XLTSSA*FLOAT(NSSOB(IBSM,IBTP))
!                  print *, 'path 1 for xnc...'
                 ELSE
                   LTTSBL = LTSSA *(LTSSA +1)/2
                   XNCOMB = XNCOMB + XLTSSA*(XLTSSA+1.0D0)*0.5
!                  print *, 'path 2 for xnc...'
                 END IF
                 NCOMB = NCOMB + LTTSBL
                 LSB = LSB + LTTSUP
                 MXSOOB = MAX(MXSOOB,LTTSUP)
                 NTTSBL = NTTSBL + 1
                 LCOL = LCOL + NSSOB(IBSM,IBTP)
               END IF
  100        CONTINUE
          end if !XLTSSA .gt. 0.0d0
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
  300 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB ,NTTSBL',
     &               NCOMB,XNCOMB,NTTSBL
      END IF
#endif
*
      END
***********************************************************************

      SUBROUTINE OCCLS_REL(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and
* IOCCLS = allowed distributions of electrons
*
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS)
*. Output
      DIMENSION  IOCCLS(NGAS,*)
*. Local scratch
#include "mxpdim.inc"
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
*
      NTEST = 00
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
#endif
*
      ISKIP = 1
      NOCCLS = 0
*. start with smallest allowed number
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: '
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space
        NEGA=0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: '
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
#if defined LUCI_DEBUG
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : '
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
#endif
            CALL ICOPY(NGAS,IOC,1,IOCCLS(1,NOCCLS),1)
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes '
           CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
#endif
*
      END
***********************************************************************

      SUBROUTINE OCC_INCLUDED(IOCC,INCLUDED)
*
* A supergroup is given as occupation for each gas space
*
* Check if all these groups are included, and return 
* answer in INCLUDED
*
* Jeppe and Timo, Oct 2001 
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Local scratch
*. Specific input
      INTEGER IOCC(NGAS)
#include "ipoist8.inc"
*
      NTEST = 000
*
      INCLUDED = 1
      DO IOBTP = 1, NGAS
        JJGRP = 0
        DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
          IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
        END DO
        IF(JJGRP.EQ.0 ) INCLUDED = 0
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation  '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        WRITE(6,*) ' Included = ',INCLUDED
      END IF
#endif
*
      END 
***********************************************************************

      SUBROUTINE ORBINF_REL(LUOUT,IPRNT)
*
* Obtain information about orbitals from shell information
*
* ======
* Output
* ======
* Orbital/spinor information in /ORBINP/
*
* Jeppe Olsen, Winter of 1991
*              Updated, July 97 : Double group symmetry added

*
      use memory_allocator
      use mospinor_info
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
#include "parluci.h"
*
      integer, allocatable :: imosp_dirac_counterx(:)

      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)

!     step 1: setup NGSOBTx
!     ---------------------
      call izero(ngsobt ,mxpngas)
      call izero(ngsobt2,mxpngas)

!     write(luwrt,*)  'nirr_dg ==> ',nirr_dg
      do isym = nirr_dg/2+1, nirr_dg 
        do igas = 1, ngas
!         write(luwrt,*)  'igas, isym ==> ',igas, isym
          ngsobt(igas)  = ngsobt(igas)  +  ngsob(isym,igas)
          ngsobt2(igas) = ngsobt2(igas) + ngsob2(isym,igas)
        end do
      end do
 
!     step 2: Number of active/occupied unbarred spinors
!     --------------------------------------------------
      call izero(ntoobs,nirr_dg)
      call izero(nocobs,nirr_dg)
      call izero(ntoobs2,nirr_dg)
      call izero(nocobs2,nirr_dg)
 
      NTOOB  = 0
      NTOOB2 = 0
      NACOB  = 0
      NACOB2 = 0
      NOCOB  = 0
      NOCOB2 = 0

      DO IGAS = 1, NGAS
!       total number of unbarred spinors
        CALL IVCSUM(NTOOBS,NTOOBS,NGSOB(1,IGAS),1,1,nirr_dg)
!       occupied/active unbarred spinors
        CALL IVCSUM(NOCOBS,NOCOBS,NGSOB(1,IGAS),1,1,nirr_dg)
        NTOOB = NTOOB + NGSOBT(IGAS)
        NOCOB = NTOOB
        NACOB = NTOOB
        CALL IVCSUM(NTOOBS2,NTOOBS2,NGSOB2(1,IGAS),1,1,nirr_dg)
        CALL IVCSUM(NOCOBS2,NOCOBS2,NGSOB2(1,IGAS),1,1,nirr_dg)
        NTOOB2 = NTOOB2 + NGSOBT2(IGAS)
        NOCOB2 = NTOOB2
        NACOB2 = NTOOB2
      END DO

!     step 3. offsets for spinors of given symmetry
!     ---------------------------------------------
      ITOOBS(1)  = 1
      ITOOBS2(1) = 1
      DO ISMOB = 2, nirr_dg
        ITOOBS(ISMOB)  =  ITOOBS(ISMOB-1) +  NTOOBS(ISMOB-1)
        ITOOBS2(ISMOB) = ITOOBS2(ISMOB-1) + NTOOBS2(ISMOB-1)
      END DO

!     debug print
!     ===========
#ifdef LUCI_DEBUG
      WRITE(LUOUT,'(/A)') ' Number of unbarred spinors per symmetry :'
      WRITE(LUOUT,'( A)') ' ========================================='
      WRITE(LUOUT,*) ' Total number of spinors ', NTOOB
      WRITE(LUOUT,*) ' Total number of occupied spinors ', NOCOB
      print *, ' *** printing NTOOBS array'
      WRITE(LUOUT,'(2x,8i4)')(NTOOBS(i),i=1,nirr_dg)

      WRITE(LUOUT,'(/A)') ' Number of barred spinors per symmetry :'
      WRITE(LUOUT,'(/A)') ' ======================================='
      WRITE(LUOUT,*) ' Total number of spinors ', NTOOB2
      WRITE(LUOUT,*) ' Total number of occupied spinors ', NOCOB2
      print *, ' *** printing NTOOBS2 array'
      WRITE(LUOUT,'(2x,8i4)')(NTOOBS2(i),i=1,nirr_dg)

      WRITE(6,*) ' Offsets for orbital of given symmetry '
      CALL IWRTMA(ITOOBS,1,nirr_dg,1,nirr_dg)
#endif
 
!     step 4: reorder array for unbarred spinors
!     ------------------------------------------
      call alloc(imosp_dirac_counterx,mxporb)
      imosp_dirac_counterx = 0

!     if we have Kramers-pairs then ngsobt and ngsobt2 are identical
!     for unrestricted code - check the global (krci-wide) use of nobpt.
!     i introduce nobpt2 here for completeness but it is not used yet...
!     s. knecht march 2011
      call icopy(ngas,ngsobt, 1,nobpt, 1)
      call icopy(ngas,ngsobt2,1,nobpt2,1)

!     reorder array for unbarred spinors
      call icopy(mxporb,imosp_dirac_counter1,1,imosp_dirac_counterx,1)

      CALL ORBORD_GAS_REL_opt(nirr_dg,mxndgirr,MXPNGAS,NGAS,NGSOB,
     &                        ITOOBS,NTOOB,IREOTS,ISMFTO,
     &                        NOBPTS,IOBPTS,ITPFTO,
     &                        imosp_dirac_counterx,imosp_luci2dirac1)
!     reorder array for barred spinors
      call icopy(mxporb,imosp_dirac_counter2,1,imosp_dirac_counterx,1)

      CALL ORBORD_GAS_REL_opt(nirr_dg,mxndgirr,MXPNGAS,NGAS,NGSOB2,
     &                        ITOOBS2,NTOOB2,IREOTS2,ISMFTO2,
     &                        NOBPTS2,IOBPTS2,ITPFTO2,
     &                        imosp_dirac_counterx,imosp_luci2dirac2)
 
      call dealloc(imosp_dirac_counterx)
      END
!***********************************************************************

      SUBROUTINE ORBORD_GAS_REL_opt(nirrep,mxnirrep,MXPNGAS,NGAS,NGSOB,
     &                              ITOOBS,NTOOB,IREOTS,ISFTO,
     &                              NOBPTS,IOBPTS,ITFTO,
     &                              imosp_dirac_counter,
     &                              imosp_luci2dirac)
*
* Obtain Reordering arrays for orbitals
* ( See note below for assumed ordering )
*
*
* GAS version
*
* =====
* Input
* =====
*  nirrep : Number of double group irreps
* mxnirrep: Max number of orbital symmetries resp. double group irreps
*  mxpngas: Max number of GAS spaces allowed by program
*  ngas   : Number of GAS spaces
*  ngsob  : Number of GAS orbitals per symmetry and space
*  itoobs : offsets of orbitals per symmetry,all types
*
* ======
* Output
* ======
*  IREOTS : Reordering array type     => symmetry
*  ISFTO  : Symmetry array for type ordered orbitals
*  NOBPTS : Number of orbitals per subtype and symmetry
*  IOBPTS : Offsets for orbitals of given subtype and symmetry
*           ordered according to input integrals
*
* ITFTO  : Type of orbital, type ordering
*
* Jeppe Olsen, Winter 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NGSOB(mxnirrep,MXPNGAS)
      DIMENSION ITOOBS(*), imosp_dirac_counter(*)
*. Output
      DIMENSION IREOTS(*),ISFTO(*)
      DIMENSION ITFTO(*), imosp_luci2dirac(*)
      DIMENSION NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
!     scratch
      logical, allocatable :: symmetry_closed_ring(:)
      integer, allocatable :: tmp_space(:)

!
!     stefan: added symmetry check for "double counting" of spinor
!     symmetries, e.g. boson symmetry "16" x "1" (up/alpha) yields "33" 
!     which is for Cinf/Cinfh equal to "1" x "1"  due to a closed 
!     representation of the approximated infinite multiplication table 
!     in C16, C8h.
!
* ==========================
* Note on order of orbitals
* ==========================
*
* The orbitals are supposed to be imported ordered symmetry-type
* ordered as
*
* Loop over symmetries of orbitals
*  Loop over GAS spaces
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of Loop over Gas spaces
* End of loop over symmetries
*
* Internally the orbitals are reordered to type symmetry order
* where the outer loop is over types and the inner loop is
* over symmetries, i.e.
*
* Loop over GAS spaces
*  Loop over symmetries of orbitals
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of loop over symmetries
* End of Loop over Gas spaces
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,*) ' Input for ORBORD '
      WRITE(6,*) ' ================='
      WRITE(6,*) ' itoobs array '
      CALL IWRTMA(itoobs(nirrep/2+1),1,nirrep/2,1,nirrep/2)
      WRITE(6,*) ' ngsob array '
      CALL IWRTMA(ngsob(nirrep/2+1,1),ngas,nirrep/2,ngas,nirrep/2)
#endif
!#undef LUCI_DEBUG

!
!     construct ISFTO, ITFTO, IREOTS, NOBPTS, IOBPTS
!     ----------------------------------------------
      call izero(nobpts,mxpngas*mxnirrep)
      call izero(iobpts,mxpngas*mxnirrep)
      call izero(itfto, mxpngas)
      call izero(isfto, mxpngas)

      allocate(symmetry_closed_ring(nirrep))
      ITSOFF = 1
      do igas = 1, ngas
!       print *, ' loop for igas ==>',igas
!       initialize (.false.) 
        symmetry_closed_ring = .false.
        ibssm = 1
        do isym = nirrep/2+1, nirrep
!         print *, ' isym, igas ==> ngsob',isym,igas,
!    &               ngsob(isym,igas)
          if(.not.symmetry_closed_ring(isym))then
!           if(ngsob(isym,igas).gt.0)then
              ibssm = itoobs(isym)
!             print '(a,2i4)',' itoobs(isym), ibssm',
!    &                          itoobs(isym), ibssm
              NPREV = 0
              DO JGAS = 1, IGAS-1
                NPREV = NPREV + NGSOB(ISYM,JGAS)
              END DO
              IADD = 0
              NOBPTS(IGAS,ISYM) = NGSOB(ISYM,IGAS)
              IOBPTS(IGAS,ISYM) = ITSOFF
!             print '(a,4i4)','ITSOFF, ISYM, NGSOB(ISYM,IGAS)-1 ',
!    &                         ITSOFF, ISYM, NGSOB(ISYM,IGAS)-1
              DO IORB = ITSOFF,ITSOFF+NGSOB(ISYM,IGAS)-1
                IADD = IADD + 1
!            print '(a,4i4)','IBSSM-1+NPREV+IADD, NPREV, IADD ==> IORB',
!    &                        IBSSM-1+NPREV+IADD, NPREV, IADD,    IORB
!            print '(a,2i4)',' IORB ==> imosp_dirac_counter',
!    &                         IORB, imosp_dirac_counter(IORB)
                i = 0
                do 
                  i = i + 1
                  if(imosp_dirac_counter(i).eq.isym) exit
                end do
                imosp_luci2dirac(IBSSM-1+NPREV+IADD) = i
!            print '(a,2i4)',' IORB ==> imosp_luci2dirac',
!    &       IORB, imosp_luci2dirac(IBSSM-1+NPREV+IADD)
                IREOTS(IORB)               = IBSSM-1+NPREV+IADD
                ITFTO(IORB)                = IGAS
                ISFTO(IORB)                = ISYM
                imosp_dirac_counter(i)     = -1
              END DO
              ITSOFF                     = ITSOFF + NGSOB(ISYM,IGAS)
              symmetry_closed_ring(isym) = .true.
!           end if ! spinor symmetry has not been accounted for previously...
          end if ! "orbital" shell is not empty...
        end do
      end do

!     release scratch memory
      deallocate(symmetry_closed_ring)

!#define LUCI_DEBUG
#if defined LUCI_DEBUG
      WRITE(6,*) ' Output from ORBORD '
      WRITE(6,*) ' ==================='
      WRITE(6,*) ' Symmetry of orbitals , type ordered '
      CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
      WRITE(6,*) ' Type => symmetry reordering array '
      CALL IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
      WRITE(6,*) ' luci => dirac reordering array '
      CALL IWRTMA(imosp_luci2dirac,1,NTOOB,1,NTOOB)
#endif
!#undef LUCI_DEBUG
#if defined LUCI_DEBUG
*
      WRITE(6,*) ' NOBPTS '
      CALL IWRTMA(NOBPTS,NGAS,nirrep,MXPNGAS,mxnirrep)
      WRITE(6,*) ' IOBPTS '
      CALL IWRTMA(IOBPTS,NGAS,nirrep,MXPNGAS,mxnirrep)
*
      WRITE(6,*) ' ITFTO array : '
      CALL IWRTMA(ITFTO,1,NTOOB,1,NTOOB)
#endif

      end
***********************************************************************

      SUBROUTINE PART_CIV3(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,MXLNG,IOCOC,ISMOST,
     &                     NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                     ITTSS_ORD)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     PART_CIV3     '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(6,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
*
*. block  zero
*
      IB = 1
      IA = 1
      IASM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, IASM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,IASM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(IASM).EQ.0) GOTO 1000
      IF(IBLTP(IASM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(IASM)
*. can this block be included
      IBSM = ISMOST(IASM)
      NSTA = NSSOA(IASM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IBLTP(IASM).EQ.1.OR.(IBLTP(IASM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(IASM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C?    write(6,*) ' IA IB IASM LBLOCK ', IA,IB,IASM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = IASM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV3'
        WRITE(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        call quit(' Insufficient buffer detected in PART_CIV3.')
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
#if defined LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(6,*) 'Output from PART_CIV3'
        WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lengths of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
#endif
*
      END
***********************************************************************

      SUBROUTINE PART_CIV4(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,MXLNG,IOCOC,ISMOST,
     &                     NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                     ITTSS_ORD)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
*
* IMPORTANT change compared to PART_CIV3:
*
* MXLNG is no longer assumed to be equals to the total number 
* of dets. Therefore we have to partition the CI vector being aware  
* of putting symmetry blocks of same TT into different batches.
*
*  
* REASON: SIGMA-vector / DENSITY matrix calculation is done wrt loops
* over NSMST! 
* 
* S. Knecht - July 02 2007
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(NSMST)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      I_SAVE_IBLOCK = 0
      N_TRY = 0
*
#ifdef LUCI_DEBUG
      WRITE(LUWRT,*) ' =================='
      WRITE(LUWRT,*) '     PART_CIV4     '
      WRITE(LUWRT,*) ' =================='
      WRITE(LUWRT,*) ' IDC = ', IDC
      WRITE(LUWRT,*) ' NSSOA array ( input ) '
      CALL IWRTMAMN(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA,LUWRT)
      WRITE(LUWRT,*) ' NSSOB array ( input ) '
      CALL IWRTMAMN(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB,LUWRT)
      WRITE(LUWRT,*) ' IBLTP array ( input ) '
      CALL IWRTMAMN(IBLTP,1,NSMST,1,NSMST,LUWRT)
#endif
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH          = NBATCH + 1
      LBATCH(NBATCH)  = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH          = 0
      LENGTHP         = 0
      NBLOCK          = 0
      IFRST           = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
!     print *, 'current block types'
!     CALL IWRTMAMN(IBLTP,1,NSMST,1,NSMST,LUWRT)
!     if(IBLTP(ISM).EQ.0)then 
!       print *, 'block type zero for ISM',ISM
!     end if
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
!     print *, 'IA, IB, IOCOC(IA,IB)',IA,IB,IOCOC(IA,IB)
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
!     write(6,*) ' ISM, IBSM, NSTA, NSTB', ISM, IBSM, NSTA, NSTB
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV4'
        WRITE(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        call quit(' Insufficient buffer detected in PART_CIV4.')
      ELSE
*
*     this batch is finished, goto next batch
*
*       is this the TT-1 block?
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(LUWRT,*) ' Not enough scratch space to include'//
     &                 ' two TT-blocks of different symmetry in'//
     &                   'a batch'
            WRITE(LUWRT,*) ' Since I cannot procede I will stop '
            WRITE(LUWRT,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV4'
            WRITE(LUWRT,*) ' block LENGTH, actual batch length ',
     &                       LENGTH, LEBATCH(NBATCH)
            WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from',MXLNG
            call quit(' N_TRY .ge. 1. in PART_CIV4.')
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,6)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
*
#ifdef LUCI_DEBUG
      IF(NTEST.NE.0) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV4'
        WRITE(LUWRT,*) '====================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '  Number of blocks included ', LBATCH(JBATCH)
          WRITE(LUWRT,*) '  TTSS and offsets and lengths of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
#endif
*
      END
***********************************************************************

      SUBROUTINE PART_CIV3_PAR(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         IOFFBLK,IOFFBTC,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,IDEBUGPRNT)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* Stefan Knecht, May 2007
*
* further modified from PART_CIV3 for parallel purposes
*
* Compared to PART_CIV3 the changes are
*
* 1 : partitioning with respect to block <--> node distribution
*
* 2 : to be continued ...
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) '     PART_CIV3_PAR    '
        WRITE(LUWRT,*) ' ====================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
C
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
C       only s-blocks corresponding to distribution
C*causes error?!?        IF( IBLOCKD(IBLOCK) .ne. MYPROC ) IBATCH(1,IBLOCK) = 0
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                  single Block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV3'
        WRITE(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH, LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        CALL QUIT( ' In PART_CIV3_PAR. ')
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV3_PAR'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lengths of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************

      SUBROUTINE PART_CIV3_SPC(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         IOFFBLK,IOFFBTC,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,IDEBUGPRNT)
*
* Partition a CI vector with given MS2 into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete ci vector is constructed
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type , relative to start of supergroup
*   IBATCH(2,*) : Beta type , relative to start of supergroup
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expanded form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, Jan. 1998
*
* Modified from PART_CIV2 for relativistic purposes
*
* Compared to PART_CIV2 the changes are
* 1 : The type numbers are absolute, no need for additional
*     reference to offset for given type
* 2 : Blocks filled from a given offset
*
* Stefan Knecht, May 2007
*
* further modified from PART_CIV3 for parallel purposes
*
* Compared to PART_CIV3 the changes are
*
* 1 : partitioning with respect to block <--> node distribution
*
* 2 : to be continued ...
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "parluci.h"
*.Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ========================='
        WRITE(LUWRT,*) '     PART_CIV3_PAR_SPC    '
        WRITE(LUWRT,*) ' ========================='
        WRITE(LUWRT,*) ' IDC = ', IDC
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' NSSOA array ( input ) '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        WRITE(LUWRT,*) ' NSSOB array ( input ) '
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
#endif
      I_SAVE_IBLOCK = 0
*
*. block  zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = IOFFBTC-1
      IBLOCK = IOFFBLK-1
      IFINI = 0
      N_TRY = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
*. Next block : Order is currently : IB, IA, ISM  (leftmost inner loop )
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' PART_CIV3 IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCK = 0
C
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(LUWRT,*) ' Not enough scratch space to include a 
     &                  single Block'
        WRITE(LUWRT,*) ' Since I cannot procede I will stop '
        WRITE(LUWRT,*) ' Insufficient buffer detected in PART_CIV3_SPC'
        WRITE(LUWRT,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from ', MXLNG
        CALL RMEMCHK
        call quit(' Insufficient buffer detected in PART_CIV3_SPC.')
      ELSE
*. This batch is finished, goto next batch
*
*       is this the TT-1 block?
*
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(LUWRT,*) ' Not enough scratch space to include a '//
     &                     ' single Block'
            WRITE(LUWRT,*) ' Since I cannot procede I will stop '
            WRITE(LUWRT,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV3_SPC'
            WRITE(LUWRT,*) '  LENGTH, LBLOCK ',LENGTH, LBLOCK
            WRITE(LUWRT,*) ' Alter GAS space of raise Buffer from',MXLNG
            CALL RMEMCHK
            CALL QUIT( ' In PART_CIV3_SPC because of N_TRY .ge. 1. ')
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,LUWRT)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
*
      IF( IDEBUGPRNT .ne. 0 ) THEN
        WRITE(LUWRT,*) 'Output from PART_CIV3_SPC'
        WRITE(LUWRT,*) '========================='
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' Number of added batches ', NBATCH-IOFFBTC+1
        DO JBATCH = IOFFBTC, NBATCH
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Info on batch ', JBATCH
          WRITE(LUWRT,*) ' *********************** '
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) '      Number of blocks included ', 
     &                          LBATCH(JBATCH)
          WRITE(LUWRT,*) '      TTSS and offsets and lengths of 
     &                          each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(LUWRT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************

      subroutine determine_ttss_block_num(IBLTP,NSSOA,NSSOB,
     &                                    NOCTPA,NOCTPB,IOFFBLK,NSMST,
     &                                    IOCOC,ISMOST,ITTSS_ORD)
!
! Partition a CI vector with given MS2 into batches of blocks.
! The length of a batch must be atmost MXLNG
!
! output: 
!         number of TTSS blocks ==> nblocks
!
! based on part_civ2 originally written by Jeppe Olsen
!
      implicit none

!-------------------------------------------------------------------------------
      integer, intent(inout) :: ioffblk
      integer, intent(in)    :: nsmst
      integer, intent(in)    :: noctpa, noctpb
      integer, intent(in)    :: nssoa(nsmst,*), nssob(nsmst,*)
      integer, intent(in)    :: iococ(noctpa,noctpb)
      integer, intent(in)    :: ibltp(nsmst)
      integer, intent(in)    :: ismost(*)
      integer, intent(in)    :: ittss_ord
!-------------------------------------------------------------------------------
      integer                :: ib, ia, ism
      integer                :: lblock, length, iblock, ifrst, ifini
      integer(8), parameter  :: max_batch_number = 1000
      integer(8), parameter  :: max_batch_length = 1000000000 * 
     &                          max_batch_number
!-------------------------------------------------------------------------------

!     block  zero
      ib     = 1
      ia     = 1
      ism    = 1
      iblock = ioffblk
      ifini  = 0
!     loop over batches of blocks
 2000 CONTINUE
      length = 0
      ifrst  = 1
!     loop over blocks in batch
 1000 continue
!     next block : order is: IB, IA, ISM (leftmost inner loop)
      if(ifrst .eq. 0)then
        call nxt_tts(ittss_ord,ia,ib,ism,ifini,noctpa,noctpb,nsmst)
      end if
      ifrst = 0
      if(ifini .eq. 1) goto 2002
      if(ibltp(ism) .eq. 0) goto 1000
      if(ibltp(ism) .eq. 2 .and. ia .lt. ib) goto 1000
      if(iococ(ia,ib) .eq. 0) goto 1000

!     can this block be included
      lblock = nssoa(ism,ia) * nssob(ismost(ism),ib)

      if(length+lblock .le. max_batch_length)then
        iblock = iblock + 1
        length = length + lblock
#ifdef LUCI_DEBUG
        print '(a,i6)', 'set ttss block counter ==> ',IBLOCK
#endif
        goto 1000
      end if ! length of batch < MXLNG?
!
 2002 ioffblk = iblock 

      end
***********************************************************************

      SUBROUTINE RSMXMN_REL(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     &                      MIN1_T,MAX1,MIN3,MAX3_T,NTEST)
*
* Construct accumulated MAX and MIN arrays for a RAS set of strings
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION MINEL(*),MAXEL(*)
*
      NORB = NORB1 + NORB2 + NORB3
*. accumulated max and min in each of the three spaces
*. ( required max and min at final orbital in each space )
COLD  MIN1A = MIN1_T
      MIN1A = MAX(MIN1_T,NEL-MAX3_T-NORB2)
      MAX1A = MAX1
*
      MIN2A = NEL - MAX3_T
      MAX2A = NEL - MIN3
*
      MIN3A = NEL
      MAX3A = NEL
*
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF (NORB1.GT.0) MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF (NORB1+NORB2.GT.0)
     &        MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
  100 CONTINUE
*
      IF( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from RSMXMN '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE SMOST_REL(NSMST,NSMCI,MXNDGIRR,ISMOST)
*
* ISMOST(ISYM,ITOTSM) : Symmetry of an internal state is ITOTSM
*                       if symmetry of 1 string is ISYM, the
*                       symmetry of the other string is
*                       ISMOST(ISYM,ITOTSM)
*
* Jeppe Olsen , Spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISMOST(MXNDGIRR,MXNDGIRR)
*
      DO 1000 ITOTSM = 1, NSMCI
       DO 900 ISTSM  = 1, NSMST
        CALL SYMCOM_REL(2,ISTSM,JSTSM,ITOTSM)
        ISMOST(ISTSM,ITOTSM) = JSTSM
  900  CONTINUE
 1000 CONTINUE
*
#ifdef LUCI_DEBUG
      WRITE(6,*) ' ==============='
      WRITE(6,*) ' Info from SMOST '
      WRITE(6,*) ' ==============='
      DO 1010 ITOTSM = 1, NSMCI
        WRITE(6,*) ' ISMOST array for ITOTSM = ', ITOTSM
        CALL IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
 1010 CONTINUE
#endif
*
      END
***********************************************************************

      SUBROUTINE SYMCOM_REL(ITASK,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1 I2 = I12
* IF(ITASK = 2 ) I1 and I12 are known, find I2
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* If obtained symmetry I1 or I2 is outside bounds,
* zero is returned.
*
* Jeppe Olsen , Spring of 1991
*
      use symmetry_setup_krci
#include "implicit.inc"
*
*
#if defined LUCI_DEBUG
      write(6,*) 'symcom test:'
      write(6,*) 'ITASK : ',ITASK
      write(6,*) 'I1    : ',I1
      write(6,*) 'I2    : ',I2
      write(6,*) 'I12   : ',I12
#endif
*
!     stefan: be careful if one of the input parameters i1, i2 or i12 
!     is zero (for whatever reason). this may lead to any result. 
!     in the previous implementation (idbgmult on a common block) 
!     it just worked (why, i do not know). in my new "modular way" 
!     i had to define these magic NULL positions otherwise the program
!     simply crashed as it ought to be. - Oct 2010 -
!     
!     boundary condition check - S. Knecht Oct 2010
*
      if(.not.(i1.lt.lbound_mat.or.i1.gt.ubound_mat))then
        if(itask.eq.2)then
          i2  = 0
          if(.not.(i12.lt.lbound_mat.or.i12.gt.ubound_mat))
     &    i2  = idbgmult(invelm(i1),i12)
!         print *, 'itask 2: I2 = ',i2
        else
          i12  = 0
          if(.not.(i2.lt.lbound_mat.or.i2.gt.ubound_mat))
     &    i12 = idbgmult(i1,i2)
!         print *, 'itask 3: I12 = ',i12
        end if
      else
        print *, ' *** failed boundary check for itask = ',itask
        if(itask.eq.2) i2   = 0
        if(itask.eq.3) i12  = 0
      end if
*
      END
***********************************************************************

      SUBROUTINE Z_BLKFO_REL(IDC,NMS2VAL,ISPC,ISM,KPCLBT,KPCLEBT,
     &                       KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK,
     &                       NBLK_MS2,IBLK_MS2,NBAT_MS2,IBAT_MS2,
     &                       NBLK_MS2_C,CALC_MS2_C,I_SET_L2BLOCK,
     &                       I_USE_PC,NPARBLOCK)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block
* KPCIBT  : Info on each block
* KPCBLTP : Block type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      use symmetry_setup_krci
      use krci_cfg
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
#include "krciprop.h"
*
      LOGICAL CALC_MS2_C
*.output
      INTEGER NBLK_MS2(*), IBLK_MS2(*), NBLK_MS2_C(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
      NTEST = 000
#if defined LUCI_DEBUG
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif
*
*. Pointers to output arrays
      NTTS = MXNTTS
      call memmar(KPCLBT ,MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,MXNTTS,'ADDL  ',1,'CI1BT ')
      call memmar(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ')
      call memmar(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*.    ^ These should be preserved after exit so put mark for 
*       flushing after
*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLKF')
*
      IDUM = 0
*
*     allowed length of each batch
CSK   IF( ICISTR.EQ.3) THEN
      L0BLOCK = MXSOOB
CSK   END IF
*
      IF( I_SET_L2BLOCK .ge. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
csk        WRITE(LUWRT,'(1X,A,1X,I20)')
csk     & '  Current available free memory in double words:',L2BLOCK
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
*
csk        WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
csk     &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1000 ) THEN
        L2BLOCK = 100 000
      END IF
*     ^ I_SET_L2BLOCK 
*
*     test if L2BLOCK has become too large...
      LTEST_BLOCK = L_COMBI_MAX
      IF( I_USE_PC .ne. 0)THEN
CSK     WRITE(LUWRT,*) ' LTEST_BLOCK, L2BLOCK, L_COMBI',
CSK  &                   LTEST_BLOCK, L2BLOCK, L_COMBI
*
        IF( L2BLOCK .gt. LTEST_BLOCK )THEN
*         reset to L_COMBI because that is already enough
          WRITE(LUWRT,*) ' *** restricting memory usage ***'
          WRITE(LUWRT,*) ' reset L2BLOCK =',L2BLOCK,' to LTEST_BLOCK =',
     &                     LTEST_BLOCK
          WRITE(LUWRT,*) ' *** restricting memory usage ***'
          L2BLOCK = LTEST_BLOCK
*
        END IF
        LTEST_BLOCK = MIN(L2BLOCK,L_COMBI_MAX)
      END IF
*
*     ... set LBLOCK value
      LBLOCK  = MAX(L0BLOCK,L2BLOCK,LTEST_BLOCK)
      IF( I_SET_L2BLOCK .gt. 1 .or. I_USE_PC .eq. 0 ) THEN
       LBLOCK = MAX(LBLOCK,LCSBLK)
      END IF
!     add restriction (if desired) for a max LBLOCK (to get some U2 calcs running...)
      if(krci_cfg_max_vec_block > 0 )then
        LBLOCK = krci_cfg_max_vec_block
      end if
      IF( I_SET_L2BLOCK .ge. 1 )
     &  WRITE(LUWRT,*) '  LBLOCK set to:',LBLOCK
#if defined LUCI_DEBUG
      WRITE(LUWRT,*) ' inside Z_BLKFO_REL: reset LBLOCK to',LBLOCK
      WRITE(LUWRT,*) ' LBLOCK, L0BLOCK, L2BLOCK, LCSBLK, LTEST_BLOCK',
     &                 LBLOCK, L0BLOCK, L2BLOCK, LCSBLK, LTEST_BLOCK
#endif
C
C     loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
C
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
C
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
C       info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KPCBLTP))
C
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
C       number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
C
C       batches of C vector
        ITTSS_ORD = 2
#if defined (VAR_MPI2)
        IF( I_USE_PC .eq. 1 ) THEN
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT)
          IDEBUGPRNT = 000
*
        ELSE IF( I_USE_PC .eq. 2 )THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_SPC(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                       WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                       IOFFBLK,IOFFBTC,
     &                       NSMST,LBLOCK,WORK(KLCIOIO),
     &                       ISMOST(1,ISM),
     &                       NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                       WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                       NPARBLOCK,IDEBUGPRNT)
          IDEBUGPRNT = 000
*
        ELSE
#endif
*
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,ISM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
#if defined (VAR_MPI2)
        END IF
#endif
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

#if defined LUCI_DEBUG
      IF(NTEST.GE.1) THEN
       WRITE(LUWRT,*) ' Total number of batches', NBATCH
       WRITE(LUWRT,*) ' Total number of blocks ', NBLOCK
*
       WRITE(LUWRT,*) ' Number of blocks and batches per MS2 value '
       WRITE(LUWRT,*) ' ==========================================='
       WRITE(LUWRT,*)
       WRITE(LUWRT,*) ' IMS2    MS2      NBLOCK      NBATCH '
       WRITE(LUWRT,*) ' ===================================='
       DO IMS2 = 1, NMS2VAL
         IATP = IST_FOR_DT(1,IMS2)
         IBTP = IST_FOR_DT(2,IMS2)
         IIMS2 = NELEC(IATP)-NELEC(IBTP)
         WRITE(LUWRT,'(4(3X,I4))') 
     &              IMS2,IIMS2,NBLK_MS2(IMS2),NBAT_MS2(IMS2)
       END DO
*
      END IF
#endif
*
* Transfer information to common block for KRMC - LUCIAREL
* communication:
      if (INTIMP.eq.6.or.INTIMP.eq.7) then
        IF( MYPROC .eq. MASTER ) THEN
          call civcinf_mc(WORK(KPCIBT),WORK(KPCLBT),WORK(KPCI1BT),
     &                    WORK(KPCLEBT),NBAT_MS2,NMS2VAL,NBLOCK)
        END IF
      end if
C
C     this part is only relevant for sigma calculations as part of 
C     the property module
C
      IF( .NOT. CALC_MS2_C ) GOTO 9999
C
C     if we run a property calculation and the call to this 
C     routine occurs from inside SIGDEN_REL we need to define 
C     atleast NBLK_MS2_C which may differ from NBLK_MS2 above 
C     since either ISSM == ICSM (fine) or ISSM != ICSM (this loop
C     is needed)
C
C     hidden information: IABS(I_SET_L2BLOCK) = ICSM
C
      ICSM_SCR = IABS(I_SET_L2BLOCK)
C     WRITE(LUWRT,*) ' ISM, ICSM_SCR is',ISM, ICSM_SCR
C
      CALL MEMMAR(KPCLBTX ,  MXNTTS,'ADDL  ',1,'CLBTX ')
      CALL MEMMAR(KPCLEBTX,  MXNTTS,'ADDL  ',1,'CLEBTX')
      CALL MEMMAR(KPCI1BTX,  MXNTTS,'ADDL  ',1,'CI1BTX')
      CALL MEMMAR(KPCIBTX ,8*MXNTTS,'ADDL  ',1,'CIBTX ')
      CALL MEMMAR(KPCBLTPX,  NSMST,'ADDL  ',2,'CBLTPX')
C
C     loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
C
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
C
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
C       info needed for generation of block info
        call memmar(KLCIOIOX,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIOX))
C
        CALL ZBLTP_REL(ISMOST(1,ICSM_SCR),NSMST,WORK(KPCBLTPX))
C
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
C       number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
C
C       batches of C vector
        ITTSS_ORD = 2
        CALL PART_CIV4(IDC,WORK(KPCBLTPX),WORK(KNSTSO(IATP)),
     &                 WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                 IOFFBLK,IOFFBTC,
     &                 NSMST,LBLOCK,WORK(KLCIOIOX),
     &                 ISMOST(1,ICSM_SCR),
     &                 NBATCH,WORK(KPCLBTX),WORK(KPCLEBTX),
     &                 WORK(KPCI1BTX),WORK(KPCIBTX),0,ITTSS_ORD)
C
C       number of BLOCKS so far
        NBLOCK = IFRMR_REL(WORK(KPCI1BTX),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBTX),1,NBATCH) - 1
*
C       WRITE(LUWRT,*) ' IMS2, NBLK_MS2_C(IMS2) is',
C    &                   IMS2,NBLOCK - NBLOCKP
        NBLK_MS2_C(IMS2) = NBLOCK - NBLOCKP
*
      END DO
*     ^ End of loop over MS2 values
C
9999  call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLKF')
*
      END
***********************************************************************

      SUBROUTINE Z_BLKFO_REL_1(IDC,NMS2VAL,ISPC,ISM,KPCIBT,KPCBLTP,
     &                         KNODELIST,KBLOCKLIST,KRCCTOS,
     &                         NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                         NBAT_MS2,IBAT_MS2,I_SET_L2BLOCK,
     &                         I_USE_PC,NPARBLOCK,IRC,IT_TTPL,
     &                         IGROUPLIST)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT  : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      use symmetry_setup_krci
      use krci_cfg
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      INTEGER ISTAT(MPI_STATUS_SIZE)
#endif
*
*.output
      DIMENSION IT_TTPL(*)
      INTEGER NBLK_MS2(*), IBLK_MS2(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
#if defined LUCI_DEBUG
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ========================== '
        WRITE(LUWRT,*) ' Output from Z_BLKFO_REL_1  '
        WRITE(LUWRT,*) ' ========================== '
        WRITE(LUWRT,*)
        WRITE(LUWRT,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif
*
*     allowed length of each batch
CSK   IF( ICISTR.EQ.3) THEN
      L0BLOCK = MXSOOB
CSK   END IF
*
      IF( I_SET_L2BLOCK .ge. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
CSK     WRITE(LUWRT,'(1X,A,1X,I20)')
CSK  & '  Current available free memory in double words:',L2BLOCK
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK_SAVE_V = 0
        L2BLOCK_SAVE_V = L2BLOCK
        L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
*       check if L2BLOCK is a positive integer
        IF( L2BLOCK .le. 0 )THEN
*         set back to default...
          WRITE(LUWRT,*) ' WARNING WARNING WARNING!!!'
          WRITE(LUWRT,*) ' MEMFAC too large, reset to default == 9'
          ISMEMFAC = 9
          L2BLOCK  = ( L2BLOCK_SAVE_V - ( ISMEMFAC * 1 000 000 ) ) / 4
        END IF
*
csk     WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
csk  &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1 ) THEN
        L2BLOCK = 100 000
      END IF
*     ^ I_SET_L2BLOCK 
*
*     test if L2BLOCK has become too large...
*
      LTEST_BLOCK = 0
      LTEST_BLOCK = MIN(L2BLOCK,L_COMBI)
*
      IF( L2BLOCK .gt. LTEST_BLOCK )THEN
*       reset to L_COMBI because that is already enough
        L2BLOCK = LTEST_BLOCK
*
      END IF
*
*     ... set LBLOCK value
*
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
!     add restriction (if desired) for a max LBLOCK (to get some U2
!     calcs running...)
      if(krci_cfg_max_vec_block > 0 )then
        LBLOCK = krci_cfg_max_vec_block
      end if
*
      IF( I_SET_L2BLOCK .eq. 2 ) RETURN
*
*     pointers to output arrays
      NTTS = MXNTTS
      CALL MEMMAR(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ')
      CALL MEMMAR(KPCBLTP,   NSMST,'ADDL  ',1,'CBLTP ')
*.    ^ These should be preserved after exit so put mark for 
*       flushing after
*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLKS')
*
*     moved here since they are not in use for GASCI_REL
      call memmar(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BT ')
*
      IDUM = 0
*
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*. Info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
*
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KPCBLTP))
*
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
*. Number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
*. Batches of C vector
        ITTSS_ORD = 2
#if defined (VAR_MPI2)
        IF( I_USE_PC .eq. 1 ) THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT)
*
        ELSE IF( I_USE_PC .eq. 2 )THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_SPC(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                       WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                       IOFFBLK,IOFFBTC,
     &                       NSMST,LBLOCK,WORK(KLCIOIO),
     &                       ISMOST(1,ISM),
     &                       NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                       WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                       NPARBLOCK,IDEBUGPRNT)
*
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
        ELSE
#endif
*
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,ISM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
#if defined (VAR_MPI2)
        END IF
#endif
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLKS')
*
*     start second part...
*
      NUM_BLOCKS  = NBLOCK
      NUM_BLOCKS2 = NBLOCK * IRC
*
      CALL MEMMAR(KNODELIST,  NUM_BLOCKS,'ADDS  ',1,'NODELS')
      CALL MEMMAR(KBLOCKLIST, NUM_BLOCKS,'ADDS  ',1,'BLCKLS')
      CALL MEMMAR(KRCCTOS,   NUM_BLOCKS2,'ADDS  ',1,'RCCTOS')
*.    ^ These should be preserved after exit so put mark for
*       flushing after
*
      MIN2  = - 2
*
      CALL ISETVC(WORK(KNODELIST),MIN2,NUM_BLOCKS)
      CALL IZERO(WORK(KBLOCKLIST),NUM_BLOCKS)
      CALL IZERO(WORK(KRCCTOS),NUM_BLOCKS2)
*     initialize block array
      CALL INI_BLOCKL(WORK(KBLOCKLIST),NBLOCK,WORK(KPCIBT))
*
      
#if defined (VAR_MPI2)
      COMPDISTL     = .TRUE.
      COMPRCCTOS    = .TRUE.
      CALL BLOCK_DISTR_DRV_REL(WORK(KNODELIST),WORK(KBLOCKLIST),
     &                         WORK(KBLOCKLIST),WORK(KRCCTOS),
     &                         IT_TTPL,IGROUPLIST)
#endif
      END
***********************************************************************

      SUBROUTINE Z_BLKFO_REL_2(IDC,NMS2VAL,ISPC,ISM,KPCIBT,KPCBLTP,
     &                         NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                         NBAT_MS2,IBAT_MS2,I_SET_L2BLOCK,
     &                         I_USE_PC,NPARBLOCK,IRC)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM.
*
* Output is given in the form of pointers to vectors in WORK
* and  NBLK_MS2,IBLK_MS2 which should be dimensioned outside
*
* KPCLBT  : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* NBLK_MS2 : Number of blocks with a given MS2 values
* IBLK_MS2 : First block with a given MS2 value
*
* Jeppe Olsen, Feb. 98
*
* For relativistic program
      use symmetry_setup_krci
      use krci_cfg
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
* new
#include "parluci.h"
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      INTEGER ISTAT(MPI_STATUS_SIZE)
#endif
*
      INTEGER NBLK_MS2(*), IBLK_MS2(*)
      INTEGER NBAT_MS2(*), IBAT_MS2(*)
*
#if defined LUCI_DEBUG
      WRITE(LUWRT,*)
      WRITE(LUWRT,*) ' ========================== '
      WRITE(LUWRT,*) ' Output from Z_BLKFO_REL_2  '
      WRITE(LUWRT,*) ' ========================== '
      WRITE(LUWRT,*)
      WRITE(LUWRT,*) ' ISM, ISPC = ', ISM,ISPC
#endif
*
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'Z_BLK2')
*     moved here since they are not in use for GASCI_REL
      call memmar(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBT  ')
      call memmar(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBT ')
      call memmar(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BT ')
*
      IDUM = 0
*
*     allowed length of each batch
CSK   ELSE IF( ICISTR.EQ.3) THEN
      L0BLOCK = MXSOOB
CSK   END IF
*
      IF( I_SET_L2BLOCK .eq. 1 ) THEN
*
*       L2BLOCK = max memory for c-vec and sigma-vec from
*       ( curr. total mem - 3 mio for string info etc.) / 4
*
        L2BLOCK = 0
*
        CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
        L2BLOCK        = LMEMFREE_PTR
*
*       we want to keep three blocks in memory at the same time
*       CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
*       division by a factor of 4 = safety!
*
        L2BLOCK_SAVE_V = 0
        L2BLOCK_SAVE_V = L2BLOCK
        L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
*       check if L2BLOCK is a positive integer
        IF( L2BLOCK .le. 0 )THEN
*         set back to default...
          WRITE(LUWRT,*) ' WARNING WARNING WARNING!!!'
          WRITE(LUWRT,*) ' MEMFAC too large, reset to default == 9'
          ISMEMFAC = 9
          L2BLOCK  = ( L2BLOCK_SAVE_V - ( ISMEMFAC * 1 000 000 ) ) / 4 
        END IF
*
csk        WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK ',
csk     &                    L0BLOCK,L2BLOCK,LCSBLK
*
        L2BLOCK = MIN(LCSBLK,L2BLOCK)
*
      ELSE IF (I_SET_L2BLOCK .eq. -1 ) THEN
        L2BLOCK = 100 000
      END IF
*     ^ I_SET_L2BLOCK 
*
*     test if L2BLOCK has become too large...
*
      LTEST_BLOCK = 0
      LTEST_BLOCK = MIN(L2BLOCK,L_COMBI)
*
      IF( L2BLOCK .gt. LTEST_BLOCK )THEN
*       reset to L_COMBI because that is already enough
        L2BLOCK = LTEST_BLOCK
*
      END IF
*
*     ... set LBLOCK value
*
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
!     add restriction (if desired) for a max LBLOCK (to get some U2
!     calcs running...)
      if(krci_cfg_max_vec_block > 0 )then
        LBLOCK = krci_cfg_max_vec_block
      end if
*
CACTU      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
CTEST      LBLOCK  = MIN(L0BLOCK,L2BLOCK)
*          old versions works without L0BLOCK and L2BLOCK
CVOLD     LBLOCK = MAX(LBLOCK,LCSBLK)
*
      NBATCH = 0
      NBLOCK = 0
*. Loop over MS2 spaces
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NOCTPA =  NOCTYP(IATP)
        NOCTPB =  NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*. Info needed for generation of block info
        call memmar(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
*
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KPCBLTP))
*
        IF(IMS2.EQ.1) THEN
          IOFFBLK = 1
          IOFFBTC = 1
        ELSE
          IOFFBLK = NBLOCK + 1
          IOFFBTC = NBATCH + 1
        END IF
        IBLK_MS2(IMS2) = IOFFBLK
        IBAT_MS2(IMS2) = IOFFBTC
*. Number of batches and blocks obtained until now
        IF(IMS2.EQ.1) THEN
          NBATCHP = 0
          NBLOCKP = 0
        ELSE
          NBATCHP = NBATCH
          NBLOCKP = NBLOCK
        END IF
*. Batches of C vector
        ITTSS_ORD = 2
#if defined (VAR_MPI2)
        IF( I_USE_PC .eq. 1 ) THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_PAR(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                     WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                     IOFFBLK,IOFFBTC,
     &                     NSMST,LBLOCK,WORK(KLCIOIO),
     &                     ISMOST(1,ISM),
     &                     NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                     WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                     NPARBLOCK,IDEBUGPRNT)
*
        ELSE IF( I_USE_PC .eq. 2 )THEN
*
*         special routine that uses only node-blocks to determine 
*         batches of a CI-vector
*
          IDEBUGPRNT = 000
*
          CALL PART_CIV3_SPC(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                       WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                       IOFFBLK,IOFFBTC,
     &                       NSMST,LBLOCK,WORK(KLCIOIO),
     &                       ISMOST(1,ISM),
     &                       NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                       WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD,
     &                       NPARBLOCK,IDEBUGPRNT)
*
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
        ELSE
#endif
*
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,ISM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*     use this line as line 4 for testing more batches ...
*     &                       NSMST,L0BLOCK,WORK(KLCIOIO),
#if defined (VAR_MPI2)
        END IF
#endif
*
*. Number of BLOCKS hitherto
        NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &         + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
*
        NBLK_MS2(IMS2) = NBLOCK - NBLOCKP
        NBAT_MS2(IMS2) = NBATCH - NBATCHP
*
      END DO
*     ^ End of loop over MS2 values

#if defined LUCI_DEBUG
      WRITE(LUWRT,*) ' Total number of batches', NBATCH
      WRITE(LUWRT,*) ' Total number of blocks ', NBLOCK
*
      WRITE(LUWRT,*) ' Number of blocks and batches per MS2 value '
      WRITE(LUWRT,*) ' ==========================================='
      WRITE(LUWRT,*)
      WRITE(LUWRT,*) ' IMS2    MS2      NBLOCK      NBATCH '
      WRITE(LUWRT,*) ' ===================================='
      DO IMS2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
        IIMS2 = NELEC(IATP)-NELEC(IBTP)
        WRITE(LUWRT,'(4(3X,I4))') 
     &             IMS2,IIMS2,NBLK_MS2(IMS2),NBAT_MS2(IMS2)
      END DO
#endif
*
*     transfer information to common block for KRMC - LUCIAREL
*     communication:
      if (INTIMP.eq.6.or.INTIMP.eq.7) then
        IF( MYPROC .eq. MASTER ) THEN
          call civcinf_mc(WORK(KPCIBT),WORK(KPCLBT),WORK(KPCI1BT),
     &                    WORK(KPCLEBT),NBAT_MS2,NMS2VAL,NBLOCK)
        END IF
      end if
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'Z_BLK2')
*
      END
C***********************************************************************

      SUBROUTINE Z_BLKFO_XPROP(IDC,NMS2VAL,ISPC,KNODELIST,KBLOCKLIST,
     &                         KRCCTOS,KXSYMDISTR,IXSYMBLK,IXSYMPAIRS,
     &                         IRC,IGROUPLIST)
C***********************************************************************
C
C     construct information about batch and block structure of all 
C     CI spaces defined by ISPC(=1) and ISM used in property run.
C
C     Output is given in the form of pointers to vectors in WORK
C
C     based on Z_BLKFO_REL. 
C
C     Written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************
C
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "parluci.h"
#include "infpar.h"
#include "krciprop.h"
#if defined (VAR_MPI2)
#include "mpif.h"
      INTEGER ISTAT(MPI_STATUS_SIZE)
#endif
      DIMENSION IXSYMBLK(*)
C
C     min length of each batch
      L0BLOCK = MXSOOB
C
C     L2BLOCK = max memory for c-vec and sigma-vec from
C     ( curr. total mem - 3 mio for string info etc.) / 4
C
      L2BLOCK = 0
C
      CALL MEMMAR(L2BLOCK,0,'SFREEM',2,'SEEFRM')
      L2BLOCK        = LMEMFREE_PTR
C
C     we want to keep three blocks in memory at the same time
C     CB,SB,VEC3(=C2). estimated scratch memory: 9 000 000 real*8
C     division by a factor of 4 = safety!
C
      L2BLOCK_SAVE_V = 0
      L2BLOCK_SAVE_V = L2BLOCK
      L2BLOCK = ( L2BLOCK - ( ISMEMFAC * 1 000 000 ) ) / 4
C     check if L2BLOCK is a positive integer
      IF( L2BLOCK .le. 0 )THEN
C       set back to default...
        WRITE(LUWRT,*) ' WARNING WARNING WARNING!!!'
        WRITE(LUWRT,*) ' MEMFAC too large, reset to default == 9'
        ISMEMFAC = 9
        L2BLOCK  = ( L2BLOCK_SAVE_V - ( ISMEMFAC * 1 000 000 ) ) / 4
      END IF
C
C
      L2BLOCK = MIN(LCSBLK,L2BLOCK)
C
C     test if L2BLOCK has become too large...
      LTEST_BLOCK = 0
      LTEST_BLOCK = MIN(L2BLOCK,L_COMBI_MAX)
C
      IF( L2BLOCK .gt. LTEST_BLOCK )THEN
C       reset to L_COMBI because that is already enough
C       WRITE(LUWRT,*) '  L2BLOCK,LTEST_BLOCK,L_COMBI_MAX ',
C    &                    L2BLOCK,LTEST_BLOCK,L_COMBI_MAX
        L2BLOCK = LTEST_BLOCK
      END IF
C
C     ... set LBLOCK value
      LBLOCK  = MAX(L0BLOCK,L2BLOCK)
C     WRITE(LUWRT,*) '  L0BLOCK,L2BLOCK,LCSBLK,LBLOCK ',
C    &                  L0BLOCK,L2BLOCK,LCSBLK,LBLOCK
C
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK  ',IDUM,'ZBXPRP')
C
      NTTS = MXNTTS
      CALL MEMMAR(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBTX ')
      CALL MEMMAR(KPCBLTP,   NSMST,'ADDL  ',1,'CBLTPX')
      CALL MEMMAR(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBTX ')
      CALL MEMMAR(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBTX')
      CALL MEMMAR(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BTX')
C
      NBLOCK_MAX = 0
      IXSYMACT   = 0
C
C     WRITE(LUWRT,*) ' ... before IJSYM loop: LBLOCK = ', LBLOCK
      DO 100 IJSYM = 1, MXPROPKRCI_SYM
C
        IF( IXSYMLIST(IJSYM) .eq. 0 ) GOTO 100
C       loop over MS2 spaces for "active" symmetry irrep
        DO IMS2 = 1, NMS2VAL
          IATP = IST_FOR_DT(1,IMS2)
          IBTP = IST_FOR_DT(2,IMS2)
C
          NOCTPA =  NOCTYP(IATP)
          NOCTPB =  NOCTYP(IBTP)
C
          IOCTPA = IBSPGPFTP(IATP)
          IOCTPB = IBSPGPFTP(IBTP)
C
C         info needed for generation of block info
          CALL MEMMAR(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CXIOIO')
          CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
          CALL ZBLTP_REL(ISMOST(1,IJSYM),NSMST,WORK(KPCBLTP))
C
          IF(IMS2.EQ.1) THEN
            IOFFBLK = 1
            IOFFBTC = 1
          ELSE
            IOFFBLK = NBLOCK + 1
            IOFFBTC = NBATCH + 1
          END IF
C         number of batches and blocks obtained until now
          IF(IMS2.EQ.1) THEN
            NBATCHP = 0
            NBLOCKP = 0
          ELSE
            NBATCHP = NBATCH
            NBLOCKP = NBLOCK
          END IF
C         batches of C vector
          ITTSS_ORD = 2
C
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,IJSYM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
C
C         number of BLOCKS so far
          NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &           + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
C
        END DO
C       ^ End of loop over MS2 values
        IXSYMBLK(IJSYM)     = NBLOCK
        NBLOCK_MAX          = MAX(NBLOCK_MAX,NBLOCK)
        IXSYMACT            = IXSYMACT + 1
        IXSYMLIST(IXSYMACT) = IJSYM
 100  CONTINUE
C     ^ End of loop over different symmetries
C     WRITE(LUWRT,*) ' ... after IJSYM loop: LBLOCK = ', LBLOCK

      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'ZBXPRP')
C
C     start second part...
C
      IF( NBLOCK_MAX .le. 0 )THEN
          WRITE(LUWRT,*) ' *** ERROR in Z_BLKFO_XPROP: max number of'// 
     &                   ' TTSS blocks is zero.'
          CALL QUIT('*** Z_BLKFO_XPROP: max number of TTSS blocks is
     &               zero ***')
      END IF
C
      NUM_BLOCKS  = NBLOCK_MAX
      NUM_BLOCKS2 = NBLOCK_MAX * IRC
C
      CALL MEMMAR(KNODELIST, NUM_BLOCKS *IXSYMACT  ,'ADDS  ',1,'NODELS')
      CALL MEMMAR(KBLOCKLIST,NUM_BLOCKS *IXSYMACT  ,'ADDS  ',1,'BLCKLS')
      CALL MEMMAR(KXSYMDISTR,            IXSYMACT  ,'ADDS  ',1,'IXSYMD')
      CALL MEMMAR(KRCCTOS,   NUM_BLOCKS2*NXSYMPAIRS,'ADDS  ',1,'RCCTOS')
C     ^ keep them after exit
C
      CALL ISETVC(WORK(KNODELIST),-2,NUM_BLOCKS *IXSYMACT)
      CALL IZERO(WORK(KBLOCKLIST),   NUM_BLOCKS *IXSYMACT)
      CALL IZERO(WORK(KXSYMDISTR),               IXSYMACT)
      CALL IZERO(WORK(KRCCTOS),      NUM_BLOCKS2*NXSYMPAIRS)
C
C     mark for local memory
      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'MARK  ',IDUM,'ZBXPR2')
C
      NTTS = MXNTTS
      CALL MEMMAR(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBTX2')
      CALL MEMMAR(KPCBLTP,   NSMST,'ADDL  ',1,'CBLTX2')
      CALL MEMMAR(KPCLBT ,  MXNTTS,'ADDL  ',1,'CLBTX2')
      CALL MEMMAR(KPCLEBT,  MXNTTS,'ADDL  ',1,'CLEBX2')
      CALL MEMMAR(KPCI1BT,  MXNTTS,'ADDL  ',1,'CI1BX2')
C
C     WRITE(LUWRT,*) ' ... before IXSYMACT loop: LBLOCK = ', LBLOCK
      DO 200 I = 1, IXSYMACT
C
        IJSYM = IXSYMLIST(I)
C       loop over MS2 spaces for "active" symmetry irrep
        DO IMS2 = 1, NMS2VAL
          IATP = IST_FOR_DT(1,IMS2)
          IBTP = IST_FOR_DT(2,IMS2)
C
          NOCTPA =  NOCTYP(IATP)
          NOCTPB =  NOCTYP(IBTP)
C
          IOCTPA = IBSPGPFTP(IATP)
          IOCTPB = IBSPGPFTP(IBTP)
C
C         info needed for generation of block info
          CALL MEMMAR(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
          CALL IAIBCM_REL(ISPC,IATP,IBTP,WORK(KLCIOIO))
C
          CALL ZBLTP_REL(ISMOST(1,IJSYM),NSMST,WORK(KPCBLTP))
C
          IF(IMS2.EQ.1) THEN
            IOFFBLK = 1
            IOFFBTC = 1
          ELSE
            IOFFBLK = NBLOCK + 1
            IOFFBTC = NBATCH + 1
          END IF
C         number of batches and blocks obtained until now
          IF(IMS2.EQ.1) THEN
            NBATCHP = 0
            NBLOCKP = 0
          ELSE
            NBATCHP = NBATCH
            NBLOCKP = NBLOCK
          END IF
C         batches of C vector
          ITTSS_ORD = 2
C
C         WRITE(LUWRT,*) 'IJSYM is',IJSYM
C
          CALL PART_CIV4(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO2(IBTP)),NOCTPA,NOCTPB,
     &                   IOFFBLK,IOFFBTC,
     &                   NSMST,LBLOCK,WORK(KLCIOIO),
     &                   ISMOST(1,IJSYM),
     &                   NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &                   WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
C
C         number of BLOCKS so far
          NBLOCK = IFRMR_REL(WORK(KPCI1BT),1,NBATCH)
     &           + IFRMR_REL(WORK(KPCLBT),1,NBATCH) - 1
C
        END DO
C       ^ End of loop over MS2 values
C
C       initialize symmetry irrep block array
        CALL INI_BLOCKL_PRP(WORK(KBLOCKLIST),NUM_BLOCKS,I,
     &                      NBLOCK,WORK(KPCIBT))
 200  CONTINUE
C     ^ End of loop over different symmetries
C     WRITE(LUWRT,*) ' ... after IXSYMACT loop: LBLOCK = ', LBLOCK

      IDUM = 0
      CALL MEMMAR(KDUM,IDUM,'FLUSM ',IDUM,'ZBXPR2')
C
#if defined (VAR_MPI2)
C     WRITE(LUWRT,*) ' ... before BLOCK_DISTR_PROP: LBLOCK = ', LBLOCK
      CALL BLOCK_DISTR_PROP(WORK(KNODELIST),WORK(KBLOCKLIST),
     &                      WORK(KRCCTOS),IXSYMPAIRS,
     &                      WORK(KXSYMDISTR),IGROUPLIST)
C     WRITE(LUWRT,*) ' ... after BLOCK_DISTR_PROP: LBLOCK = ', LBLOCK
#endif
      END
#if defined (VAR_MPI2)
C***********************************************************************      

      SUBROUTINE BLOCK_DISTR_PROP(NPARBLOCK_X,IBLOCKL_X,RCCTOS_X,
     &                            IXSYMPAIRS,ISYM_DISTR,IGROUPLIST)
C***********************************************************************      
C
C     driver routine for symmetry irrep dependent block distribution 
C
C     written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************
C
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "krciprop.h"
#include "mpif.h"
#include "parluci.h"
#include "infpar.h"
#include "cands.inc"
#include "mxpdim.inc"
      INTEGER   ISTAT(MPI_STATUS_SIZE), ISYM_DISTR(IXSYMACT)
      INTEGER   NPARBLOCK_X(NUM_BLOCKS,*),IBLOCKL_X(NUM_BLOCKS,*)
      INTEGER   IXSYMPAIRS(MXPROPKRCI_SYM,*)
      INTEGER   RCCTOS_X(NUM_BLOCKS2,*)
C
      IXPAIR = 0
C
C     loop over all "active" symmetry irrep pairs
C     for simplicity we assume the same even distribution for a given 
C     sigma irrep no matter from which T_X x C pair it will be generated
C
C     COMPDISTL and COMPRCCTOS on common block in include/parluci.h.
C
C     WRITE(LUWRT,*) ' block distribution for sum_(act. irreps) =',
C    &                 IXSYMACT
C
      DO I = 1, IXSYMACT
        IISYM = IXSYMLIST(I)
        DO J = 1, IXSYMACT
           IJSYM = IXSYMLIST(J)
           IF( IXSYMPAIRS(IISYM,IJSYM) .eq. 1 )THEN
             IXPAIR = IXPAIR + 1
C
C            reset IXSYMPAIRS to current pair number
C            ---------------------------------------
             IXSYMPAIRS(IISYM,IJSYM) = IXPAIR
             IF( ISYM_DISTR(J) .eq. 0 ) THEN
               COMPDISTL     = .TRUE.
               COMPRCCTOS    = .TRUE.
               ISYM_DISTR(J) = 1
             ELSE
               COMPDISTL  = .FALSE.
               COMPRCCTOS = .TRUE.
             END IF
C            set correct symmetry for C and sigma
             ICSM = IISYM
             ISSM = IJSYM
C
C            set T operator arrays
             ISYM_T    = IDBGMULT(INVELM(ICSM),ISSM)
             IPRHAM_CI = 0
             NO_TTWO   = 1
C	 WRITE(LUWRT,*) ' inside loop I,J (1): LBLOCK',I,J,LBLOCK
             CALL SET_HOP_DBG(IHTYPE_X,NO_TTWO,ISYM_T,IPRHAM_CI)
C	 WRITE(LUWRT,*) ' inside loop I,J (2): LBLOCK',I,J,LBLOCK
C
             CALL BLOCK_DISTR_DRV_REL(NPARBLOCK_X(1,J),IBLOCKL_X(1,J),
     &                                IBLOCKL_X(1,I),RCCTOS_X(1,IXPAIR),
     &                                IDUMMY,IGROUPLIST)
C	 WRITE(LUWRT,*) ' inside loop I,J (3): LBLOCK',I,J,LBLOCK
           END IF
        END DO
      END DO
      END
#endif
C***********************************************************************      

      SUBROUTINE GET_SC_PAIRS(IXSYMPAIRS,IXSYM,NX_OP,ISYMOPX,
     &                        IXROOTS,IXROOTSYM)
C***********************************************************************      
C
C     obtain all relevant sigma/C symmetry pairs
C
C     written by S. Knecht - Oct 2008
C
C     Last revision :
C
C***********************************************************************      
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "parluci.h"
#include "mxpdim.inc"
#include "krciprop.h"
      DIMENSION IXSYMPAIRS(IXSYM,IXSYM), ISYMOPX(*),IXROOTSYM(*)
C
      IONE = 1
C
      DO 100 I = 1, NX_OP
         IF( ISYMOPX(I) .eq. 0 ) GOTO 100
         DO J = 1, IXROOTS
           ICSM = IXROOTSYM(J)
           ISSM = IDBGMULT(ICSM,ISYMOPX(I))
           IXSYMPAIRS(ICSM,ISSM) = IONE
         END DO
 100  CONTINUE
C
C     number of symmetry pairs
      NXSYMPAIRS = 0
      DO I = 1, IXSYM
        DO J = 1, IXSYM
          IF( IXSYMPAIRS(I,J).eq. IONE ) NXSYMPAIRS = NXSYMPAIRS + 1
        END DO
      END DO
c     WRITE(LUWRT,*) ' GET_SC_PAIRS: number of symmetry pairs = ',
c    &                 NXSYMPAIRS
C
       END
***********************************************************************

      SUBROUTINE ZBLTP_REL(ISMOST,MAXSYM,ICBLTP)
*
* Generate vector ICBLTP giving type of each block
*
* ICBLTP gives type of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
*. Input
      DIMENSION ISMOST(*)
*. Output
      DIMENSION ICBLTP(*)
*
      call izero(icbltp,maxsym)
      do iasym = 1, maxsym
        icbltp(iasym) = 1
      end do

      END
***********************************************************************

      SUBROUTINE ZNELFSPGP_REL(NTESTG)
*
* Generate for each supergroup the number of electrons in each active
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cgas.inc"
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
#include "gasstr.inc"
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP  = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP)
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
#if defined LUCI_DEBUG
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
#endif
*
      END
