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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck dftgrd */
      SUBROUTINE DFTGRD(DMAT,IPRINT)
C                                                                      C
C     Interface routine for grid points                                C
C     T. Helgaker                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        use num_grid_cfg
        use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"
      PARAMETER (D0=0.0D0)

#include "dftwrk.h"

      DIMENSION DMAT(*)
      real(8), allocatable :: kzan(:)
      real(8), allocatable :: kc(:)
      real(8), allocatable :: ksp(:)
      real(8), allocatable :: kx8(:)
      real(8), allocatable :: ky8(:)
      real(8), allocatable :: kz8(:)
      real(8), allocatable :: kw8(:)
      real(8), allocatable :: kpsmu(:)
      real(8), allocatable :: krj(:)
      real(8), allocatable :: kaccum(:)
      real(8), allocatable :: kxmujn(:)
      real(8), allocatable :: kxmuj2(:)
      real(8), allocatable :: krij(:)  
      real(8), allocatable :: kaij(:)  
      real(8), allocatable :: krbc(:)  
      real(8), allocatable :: krvec(:) 
      real(8), allocatable :: kchimt(:)
      real(8), allocatable :: kaa(:)   
      integer, allocatable :: knuco(:)

#include "ibtfun.h"

      TIMSTR = SECOND()

C
C     ********************************
C     ***** Set up common blocks *****
C     ********************************
C
C     INFOA has been replaced
C     =====

      call alloc(kzan, nucdep)
      call alloc(kc,   nucdep*3)
      CALL SET_GRID_INFO(NAT,NUM,kzan,kc)
C
C     DFTWRK
C     ======
C
      RNDMAX  = 20.0D0
      RADERR  = num_grid_cfg_RADINT
      NSCHEME = num_grid_cfg_ANGINT
      LEBMIN  = num_grid_cfg_ANGMIN

C
C     *******************************************
C     ***** Prepare arguments for CONSTRUCT *****
C     *******************************************
C
      IBIG    = 1000000
      IVECL   = 100
      NNAT    = NAT*(NAT-1)/2
      N2AT    = NAT*NAT 

!     radovan: workaround to get nonzero length array
      if (nnat == 0) nnat = 1

      call alloc(KSP   ,4*NUM    )
      call alloc(KX8   ,IBIG     )
      call alloc(KY8   ,IBIG     )
      call alloc(KZ8   ,IBIG     )
      call alloc(KW8   ,IBIG     )
      CALL alloc(KPSMU ,IVECL*NAT)
      CALL alloc(KRJ   ,IVECL*NAT)
      CALL alloc(KACCUM,IVECL    )
      CALL alloc(KXMUJN,IVECL    )
      CALL alloc(KXMUJ2,IVECL    )
      CALL alloc(KRIJ  ,NNAT     )
      CALL alloc(KAIJ  ,NNAT     )
      CALL alloc(KRBC  ,N2AT     )
      CALL alloc(KRVEC ,3*IVECL  )
C..added for density-based atomic size estimates
      CALL alloc(KCHIMT,NAT*NAT)
C..added for the MOLCAS scheme:
      NDIM=2*NHTYP*NUCIND
      CALL alloc(KAA  ,2*NDIM)
      kaa = 0.0d0
      CALL alloc(KNUCO,NDIM  )


      if (num_grid_cfg_import_grid) go to 999

      lunit_grid = 44
      OPEN(lunit_grid,
     &     FORM   = 'FORMATTED',
     &     STATUS = 'UNKNOWN',
     &     ACCESS = 'SEQUENTIAL',
     &     FILE   = 'numerical_grid')
      REWIND lunit_grid

C
C     Calculate abscissas and weights
C     ===============================
C
      CALL CONSTRUC(KSP,KX8,KY8,KZ8,KW8,
     &              KPSMU,KRJ,KACCUM,
     &              KXMUJN,KXMUJ2,
     &              KRIJ,KAIJ,IVECL,IBIG,
     &              KRBC,KRVEC,
     &              KNUCO,KAA,DMAT,KCHIMT,
     &              NAT,NUM,KZAN,KC,
     &              IPRINT)


!     check accuracy of numerical integration
      if (num_grid_cfg_integration_check_level > 0) then
        call grtest()
      end if

      close(lunit_grid, status='keep')

      TIMGRD = SECOND() - TIMSTR
      WRITE (LUPRI,'(1X)')
      CALL TIMTXT('>>>> Time used in DFTGRD:',TIMGRD,LUPRI)
      WRITE (LUPRI,'(1X)')

 999  CONTINUE

      call dealloc(kzan)
      call dealloc(kc)
      call dealloc(ksp)
      call dealloc(kx8)
      call dealloc(ky8)
      call dealloc(kz8)
      call dealloc(kw8)
      call dealloc(kpsmu)
      call dealloc(krj)
      call dealloc(kaccum)
      call dealloc(kxmujn)
      call dealloc(kxmuj2)
      call dealloc(krij)
      call dealloc(kaij)
      call dealloc(krbc)
      call dealloc(krvec)
      call dealloc(kchimt)
      call dealloc(kaa)
      call dealloc(knuco)

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck construc */
      SUBROUTINE CONSTRUC(WRKSP,X8,Y8,Z8,WT8,PSMU,RJ,ACCUM,XMUIJN,
     &                    XMUIJ2,RIJ,AIJ,IVECL,IBIG,RBC,RVEC,
     &                    NUCORB,AA,DMAT,CHIMT,
     &                    NAT,NUM,ZAN,C,IPRINT)
C                                                                      C
C      Creates a quadrature grid and writes it to disk                 C
C      Original version by C.W.Murray. Rewritten by                    C
C      A.M.Lee and D.J.Tozer.                                          C
C                                                                      C
C      Adapted for Dalton by T. Helgaker                               C
C      Adapted for DIRAC by T.Saue                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       use num_grid_cfg
       use interface_grid

#include "implicit.h" 
#include "codata.h"
#include "priunit.h"
#include "consts.h"

      LOGICAL SWITCH,ADDSML
      DIMENSION X8(IBIG),Y8(IBIG),Z8(IBIG),WT8(IBIG)


#include "dftwrk.h"

#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"
#include "dcbbas.h"
#include "dcbham.h"

      DIMENSION ZAN(*), C(3,NAT)
      DIMENSION WRKSP(4*NUM),
     &          ACCUM(IVECL),XMUIJN(IVECL),XMUIJ2(IVECL),
     &          PSMU(IVECL,NAT),RJ(IVECL,NAT)
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND),
     &     DMAT(*),CHIMT(NAT,NAT)
C
C     CARE: Some equivalencing through argument list!
C
      DIMENSION RIJ(NAT*(NAT-1)/2),AIJ(NAT*(NAT-1)/2)
      DIMENSION RBC(NAT,NAT)
      DIMENSION RVEC(IVECL,3)
      DIMENSION RADWT(2000), RADNDE(2000)
      DIMENSION ANGWT(9000), XANG(9000), YANG(9000), ZANG(9000)
      DIMENSION XPASC(20)
#include "ibtfun.h"
C
      ADDSML = num_grid_cfg_force_4c_grid.OR.
     &         (MC.GT.1.AND.(.NOT.LEVYLE))
C
C     Get estimates for atomic size ratios for use in the Becke
C     partitioning scheme.
C
      CALL TITLER('Slater-Bragg atomic size estimates','*',122)
C
C     Initialize CHIMT with Bragg-Slater data
C
      DO I = 1,NAT
         DO J = 1,NAT
            CHIMT(I,J) = BRAGG(NINT(ZAN(I)))/BRAGG(NINT(ZAN(J)))
         END DO
      END DO
      CALL HEADER('Slater-Bragg relative atomic sizes',-1)
      CALL OUTPUT(CHIMT,1,NAT,1,NAT,NAT,NAT,1,LUPRI)
C
C     Get density-based atomic radii if requested AND coefficients 
C     exist.
C
      IF(num_grid_cfg_estimate_radii) THEN
         CALL TITLER('Density based atomic size estimates','*',122)
         CALL GETCHI(DMAT,CHIMT,NAT,C,IPRINT)
         CALL HEADER('Density based relative atomic sizes',-1)
         CALL OUTPUT(CHIMT,1,NAT,1,NAT,NAT,NAT,1,LUPRI)
         WRITE(LUPRI,*)
      ELSE
         WRITE(LUPRI,'(/A)')
     &        'DFTGRD: ATSIZE not requested or no '//
     &        'trial vectors found. Slater-Bragg '//
     &        'radii will be used.'
         WRITE(LUPRI,*)
      ENDIF

C
C     USE INPUTTED INFORMATION IN COMMON BLOCK DFTWRK TO DEFINE THE TYPE
C     OF QUADRATURE
C
C     Experimental value of NTRANS
C
      NTRANS=10
C
C     Form Pascals triangle in XPASC for fuzzy Voronoi polyhedra code
C
      ISIGN=-1
      DO 5 I=NTRANS,1,-1
         ISIGN=-ISIGN
         XPASC(I+1) = ISIGN*FACULT(NTRANS)/(FACULT(I)*FACULT(NTRANS-I))
         XPASC(1)   = 1.0D0
5     CONTINUE

      APASC=0.0D0
      DO 6 I=1,NTRANS+1
         XPASC(I)=XPASC(I)/DBLE(2*I-1)
         APASC=APASC+XPASC(I)
6     CONTINUE
      APASC=0.5D0/APASC

      ITEMP=0
      RBC = 0.0d0
      DO 7 INA=1,NAT
      DO 7 JNA=1,INA-1
          ITEMP=ITEMP+1
          RIJ(ITEMP) = 1.0D0/(DSQRT((C(1,INA)-C(1,JNA))**2
     &                             +(C(2,INA)-C(2,JNA))**2
     &                             +(C(3,INA)-C(3,JNA))**2))
          RBC(INA,JNA) = 1.0D0/RIJ(ITEMP)
          RBC(JNA,INA) = 1.0D0/RIJ(ITEMP)
          CHI = CHIMT(INA,JNA)
          TEMP=(CHI-1.0D0)/(CHI+1.0D0)
          AIJ(ITEMP)=TEMP/(TEMP*TEMP-1.0D0)
          IF (AIJ(ITEMP).GT.0.5D0) AIJ(ITEMP)=0.5D0
          IF (AIJ(ITEMP).LT.-0.5D0) AIJ(ITEMP)=-0.5D0
7     CONTINUE
C
C     Get information about basis sets
C
      CALL NUCBAS(NUCORB,AA,IPRINT)
C
C     Loop over atoms
C     ===============
C 
      MXOPR = 0 
      NTOTAL = 0
      NALL = 0
      NATOM = 1
      WRITE(LUPRI,'(A)')
     &  'Atom  Deg   Rmin        Rmax        Step size     #rp     #tp',
     &  '============================================================='
      DO 10 IATOM = 1, NUCIND
      DO 10 ISYMOP = 0, MXOPR         
      IF (IBTAND(ISYMOP,ISTBNU(IATOM)).EQ.0) THEN 
         NDEG  = MULT(ISTBNU(IATOM))
         IF(NOORBT(IATOM)) GOTO 10
         IF(NAMN(IATOM)(1:2).EQ.'Gh') GOTO 10
         MULA  = ISTBNU(IATOM)
         NDEG  = MULT(MULA)
         IPT   = 0
         INDEX = 0
         NDEX  = 1
         NTHIS = 0
         NTHAT = 0
C
C     Radial quadrature 
C     =================
C     
C     As proposed by 
C     Roland Lindh, Per-Aake Malmqvist and Laura Gagliardi
C     Theor. Chem. Acc. (2001) 106:178 
C
C     Radial point:         
C     RADNDE(K) = GRDC(EXP[K*H] - 1)
C               = (RADNDE(K-1)+GRDC)*EXP[H]-GRDC
C
C     where appears the constant
C
C          GRDC = RL/(EXP[H] - 1),         
C                   
C     calculated from the step size H and the inner point RL.
C
C     Radial weight:         
C
C     RADWT(K) = (RADNDE(K)+GRDC)*RADNDE(K)*RADNDE(K)*H
C           
         CALL LMGGRID(H,RL,RH,RADERR,IATOM,ADDSML,NUCORB,AA,IPRINT)
         EPH  = EXP(H)
         GRDC = RL/(EPH-D1)
         NR   = NINT(LOG(D1+(RH/GRDC))/H)
         RADNDE(1) = RL
         RADWT(1)  = (RL+GRDC)*RL*RL*H
         DO IR = 2,NR
           RADNDE(IR) = (RADNDE(IR-1)+GRDC)*EPH-GRDC
           RADWT(IR) = (RADNDE(IR)+GRDC)*RADNDE(IR)*RADNDE(IR)*H
         ENDDO
C 
C        loop over radial points
C        =======================
C
         SWITCH = .FALSE.
         RBRAGG = BRAGG(NINT(ZAN(NATOM)))/(5.0D0*XTANG)
         DO 20 IR = 1, NR
            RWT   = RADWT(IR)
            RNODE = RADNDE(IR)
COLAV            IF (RNODE.GT.RNDMAX) GOTO 20
C 
C           angular abscissas and weights for this radial point
C           ===================================================
C
            IF (.NOT.SWITCH) THEN
               CALL SLEB(XANG,YANG,ZANG,ANGWT,NSCHEME,NANG,NATOM,RNODE,
     &                   RBRAGG)
               IF (RNODE.GT.RBRAGG) SWITCH = .TRUE.
            END IF
            IF (NR*NANG.GT.IBIG) THEN
               WRITE(LUPRI,*) 'Have:',IBIG,' need:',NR*NANG,NANG
               CALL QUIT('Storage error in CONSTRUCT')
            END IF
            IF(IPRINT.GE.3) THEN
              WRITE(LUPRI,'(A,I5,F16.8,I5,I5,5X,L1)') 
     &          '* Grid:',IR,RNODE,NANG,NSCHEME,SWITCH
            ENDIF
C
C           loop over angular points
C           ========================
C
            DO 30 iy = 1, NANG
               IPT   = IPT + 1
               INDEX = INDEX + 1
               X8(INDEX)  = XANG(iy)*RNODE+C(1,NATOM)
               Y8(INDEX)  = YANG(iy)*RNODE+C(2,NATOM)
               Z8(INDEX)  = ZANG(iy)*RNODE+C(3,NATOM)
               WT8(INDEX) = RWT*ANGWT(iy)
C
C              Becke partitioning
C              ==================
C
               IF (IPT.EQ.IVECL) THEN
                  CALL BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,
     &                       X8(NDEX),Y8(NDEX),Z8(NDEX),WT8(NDEX),
     &                       NAT,NUC,ZAN,C,
     &                       ACCUM,XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,
     &                       NATOM)
                  NDEX   = NDEX  + IPT
CVT               NTOTAL = NTOTAL+ IPT
                  NTHIS  = NTHIS + IPT
                  INDEX  = INDEX - IVECL + IPT
                  IPT    = 0
               END IF
30          CONTINUE
20       CONTINUE
C
C        Becke partitioning
C        ==================
C
         IF (IPT.NE.0) THEN
            JPT = IPT
            CALL BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,
     &                 X8(NDEX),Y8(NDEX),Z8(NDEX),WT8(NDEX),
     &                 NAT,NUC,ZAN,C,ACCUM,
     &                 XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,NATOM)
            NDEX   = NDEX   + IPT
CVT         NTOTAL = NTOTAL + IPT
            NTHIS  = NTHIS  + IPT
            INDEX  = INDEX  - JPT + IPT
            IPT = 0
         END IF
         IF(NDEG.GT.0) THEN
C          Scale weights with multiplicity of center
           FAC = dble(NDEG) 
           CALL DSCAL(NTHIS,FAC,WT8,1)
         ENDIF
         NTHAT = NTHIS


!        drop points that are redundant and rescale
         if (num_grid_cfg_zipgrid) then
           if (maxrep > 0) then
             call zip_grid(iatom, x8, y8, z8, wt8, nthis)
           end if
         end if


C
C        Write to file
C        ====================
C
         lunit_grid = 44
         call num_grid_write(x8, y8, z8, wt8, 
     &                        lunit_grid, nthis)

         NTOTAL = NTOTAL + NTHIS
         NALL = NALL + NTHAT
         NATOM = NATOM + NDEG
         WRITE(LUPRI,'(A4,I5,3E12.3,2I8)') NAMN(IATOM),NDEG,
     &         RL,RH,H,NR,NTHIS

      END IF
10    CONTINUE

      NLAST = -1
      call num_grid_write(x8, y8, z8, wt8, 
     &                     lunit_grid, nlast)

      WRITE (LUPRI,'(/,2X,A,I10,A,F5.1,A/)') 
     &     ' Number of grid points in quadrature:',NTOTAL,
     &     ' (',FLOAT(100*NTOTAL)/FLOAT(NALL),'%)'


      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck nleb */                                                     C
      FUNCTION NLEB(NSCHEME)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 

      IF(NSCHEME.GT.131)THEN
         WRITE(6,*) 'Requested accuracy:        ',NSCHEME
         WRITE(6,*) 'Maximum accuracy available:',131
         CALL QUIT('NLEB: Requested accuracy not available !')
      ELSEIF(NSCHEME.LE.131.AND.NSCHEME.GT.125) THEN
         NTOT = 5810
      ELSEIF(NSCHEME.LE.125.AND.NSCHEME.GT.119) THEN
         NTOT = 5294
      ELSEIF(NSCHEME.LE.119.AND.NSCHEME.GT.113) THEN
         NTOT = 4802
      ELSEIF(NSCHEME.LE.113.AND.NSCHEME.GT.107) THEN
         NTOT = 4334
      ELSEIF(NSCHEME.LE.107.AND.NSCHEME.GT.101) THEN
         NTOT = 3890
      ELSEIF(NSCHEME.LE.101.AND.NSCHEME.GT.95) THEN
         NTOT = 3470
      ELSEIF(NSCHEME.LE.95.AND.NSCHEME.GT.89) THEN
         NTOT = 3074
      ELSEIF(NSCHEME.LE.89.AND.NSCHEME.GT.83) THEN
         NTOT = 2702
      ELSEIF(NSCHEME.LE.83.AND.NSCHEME.GT.77) THEN
         NTOT = 2354
      ELSEIF(NSCHEME.LE.77.AND.NSCHEME.GT.71) THEN
         NTOT = 2030
      ELSEIF(NSCHEME.LE.71.AND.NSCHEME.GT.65) THEN
         NTOT = 1730
      ELSEIF(NSCHEME.LE.65.AND.NSCHEME.GT.59) THEN
         NTOT = 1454
      ELSEIF(NSCHEME.LE.59.AND.NSCHEME.GT.53) THEN
         NTOT = 1202
      ELSEIF(NSCHEME.LE.53.AND.NSCHEME.GT.47) THEN
         NTOT =  974
      ELSEIF(NSCHEME.LE.47.AND.NSCHEME.GT.41) THEN
         NTOT =  770
      ELSEIF(NSCHEME.LE.41.AND.NSCHEME.GT.35) THEN
         NTOT =  590
      ELSEIF(NSCHEME.LE.35.AND.NSCHEME.GT.31) THEN
         NTOT =  434
      ELSEIF(NSCHEME.LE.31.AND.NSCHEME.GT.29) THEN
         NTOT =  350
      ELSEIF(NSCHEME.LE.29.AND.NSCHEME.GT.23) THEN
         NTOT =  302
      ELSEIF(NSCHEME.LE.23.AND.NSCHEME.GT.21) THEN
         NTOT =  194
      ELSEIF(NSCHEME.LE.21.AND.NSCHEME.GT.19) THEN
         NTOT =  170
      ELSEIF(NSCHEME.LE.19.AND.NSCHEME.GT.17) THEN
         NTOT =  146
      ELSEIF(NSCHEME.LE.17.AND.NSCHEME.GT.15) THEN
         NTOT =  110
      ELSEIF(NSCHEME.LE.15.AND.NSCHEME.GT.11) THEN
         NTOT =   86
      ELSEIF(NSCHEME.LE.11.AND.NSCHEME.GT.9) THEN
         NTOT =   50
      ELSEIF(NSCHEME.LE.9.AND.NSCHEME.GT.7) THEN
         NTOT =   38
      ELSEIF(NSCHEME.LE.7.AND.NSCHEME.GT.5) THEN
         NTOT =   26
      ELSEIF(NSCHEME.LE.5.AND.NSCHEME.GT.3) THEN
         NTOT =   14
      ELSEIF(NSCHEME.LE.3) THEN
         NTOT =    6
      END IF

      NLEB = NTOT

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck nlebset */                                                  C
      FUNCTION NLEBSET(NTOT)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 
 
      IF (NTOT.GE.1453) THEN
         NSCHEME=64
      ELSEIF (NTOT.GE.1201) THEN
         NSCHEME = 59
      ELSEIF (NTOT.GE.973) THEN
         NSCHEME = 53
      ELSEIF (NTOT.GE.769) THEN
         NSCHEME=47
      ELSEIF (NTOT.GE.589) THEN
         NSCHEME=41
      ELSEIF (NTOT.GE.433) THEN
         NSCHEME=35
      ELSEIF (NTOT.GE.349) THEN
         NSCHEME=31
      ELSEIF (NTOT.GE.301) THEN
         NSCHEME=29
      ELSEIF (NTOT.GE.193) THEN
         NSCHEME=23
      ELSEIF (NTOT.GE.169) THEN
         NSCHEME=21
      ELSEIF (NTOT.GE.145) THEN
         NSCHEME=19
      ELSEIF (NTOT.GE.109) THEN
         NSCHEME=17
      ELSEIF (NTOT.GE.85) THEN
         NSCHEME=15
      ELSEIF (NTOT.GE.73) THEN
         NSCHEME=13
      ELSEIF (NTOT.GE.49) THEN
         NSCHEME=11
      ELSE   
         NSCHEME=9
      ENDIF
C
      NLEBSET=NSCHEME
 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE LMGGRID(H,RL,RH,RADERR,IATOM,ADDSML,NUCORB,AA,IPRINT)
C***********************************************************************
C
C     Find parameters definining the radial grid defined here
C     Roland Lindh, Per-Aake Malmqvist and Laura Gagliardi
C     Theor. Chem. Acc. (2001) 106:178
C
C     H  - grid spacing      
C     RL - inner grid point
C     RH - outer grid point
C
C     Input variable:
C       RADERR - acceptable discretization error
C       IATOM  - index of symmetry-independent nucleus
C       ADDSML - include small component basis functions in the generation of
C                grid parameters (4c grid)      
C       NUCORB(LL,2,IATOM) - number of large(1)/small(2) component basis functions
C                            of angular momentum L = LL - 1 on atom IATOM      
C       AA(2,LL,IC,NUCIND) - highest(1)/lowest(2) exponent of angular momentum 
C                            L = LL - 1 of large(IC=1)/small(IC=2) components
C                            on atom IATOM
C     
C     Radial point:         
C     RADNDE(K) = GRDC(EXP[K*H] - 1)
C               = (RADNDE(K-1)+GRDC)*EXP[H]-GRDC
C
C     where appears the constant
C
C          GRDC = RL/(EXP[H] - 1),         
C                   
C     calculated from the step size H and the inner point RL.
C
C     Radial weight:         
C
C     RADWT(K) = (RADNDE(K)+GRDC)*RADNDE(K)*RADNDE(K)*H
C           
C***********************************************************************
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"      
#include "dummy.h"      
      PARAMETER(D0=0.0D0,D2=2.0D0,D3=3.0D0)
C     
      CHARACTER SPDCAR*1
#include "ccom.h"
#include "nuclei.h"
      LOGICAL ADDSML
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)      
C     Grid spacing....
      H = DUMMY
      DO LL = 1,NHTYP
        L = LL-1
        NBAS = NUCORB(LL, 1, IATOM)
!       radovan:
!       if Levy-Leblond or X2C, ignore small component
!       unless user forces 4-c grid
        IF (ADDSML) THEN
          NBAS = NBAS + NUCORB(LL, 2, IATOM)
        END IF
        IF(NBAS.GT.0) THEN
          HTMP = DISERR(RADERR,L) ! discretization error
          H    = MIN(H,HTMP)
          IF(IPRINT.GE.3) THEN
            WRITE(LUPRI,'(3X,A1,A,F6.3)')
     &          SPDCAR(L),'-orbitals --> ',HTMP
          ENDIF
        ENDIF
      ENDDO
C...  Inner grid point
      AH = AA(1, 1, 1, IATOM)
!     radovan:
!     if Levy-Leblond or X2C, ignore small component
!     unless user forces 4-c grid
      if (ADDSML) then
        AH = MAX(AH, AA(1, 1, 2, IATOM))
      end if
      AH = D2*AH
      RL = ((1.9D0+LOG(RADERR))/D3)-(LOG(AH)/D2)
      RL = EXP(RL)
C...  Outer point
      RH = D0
      DO LL = 1,NHTYP
        L = LL-1
        AL= DUMMY
        IF(NUCORB(LL,1,IATOM).GT.0) AL=AA(2,LL,1,IATOM)
!       radovan:
!       if Levy-Leblond or X2C, ignore small component
!       unless user forces 4-c grid
        IF (ADDSML) THEN
          IF(NUCORB(LL,2,IATOM).GT.0) AL=MIN(AL,AA(2,LL,2,IATOM))
        END IF
        IF(AL.LT.DUMMY) THEN
          AL    = AL+AL
          RHTMP = OUTERR(AL,L,RADERR)           
          RH    = MAX(RH,RHTMP)
          IF(IPRINT.GE.3) THEN
            WRITE(LUPRI,'(3X,A1,A,F6.3)')
     &            SPDCAR(L),'-orbitals --> ',RHTMP
           ENDIF
         ENDIF
       ENDDO
       RETURN
       END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck sleb */                                                     C
      SUBROUTINE SLEB(X,Y,Z,WT,NSCHEME,NCOUNT,NATOM,RNODE,RBRAGG)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       use num_grid_cfg

#include "implicit.h" 
#include "priunit.h" 
#include "pi.h"

      PARAMETER (MAXNUM=1455)
      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)

C
C     Lebedev schemes programmed by Cristoph van Wullen
C     and used by kind permission.
C
C     For principal references, see file Lebedev-Laikov.F
C 
      NTOT=NLEB(NSCHEME)

C     Chose the Lebedev grid to use as a function of the radius (rnode)
C     closer to the nuclei scale the number of angular points as

      lebmin = num_grid_cfg_angmin
      NOLD=NSCHEME
      IF(NSCHEME.GT.LEBMIN.AND.RNODE.LT.RBRAGG.AND.
     &    .NOT.num_grid_cfg_no_pruning) THEN
         IANG=NINT(NTOT*RNODE/RBRAGG)
         NTOT=MIN(NTOT,IANG)
         NSCHEME=NLEBSET(NTOT)
         NTOT=NLEB(NSCHEME)
C
C     If the pruned value of NSCHEME is less than LEBMIN, we set 
C     NSCHEME = LEBMIN, map NSCHEME to NTOT and then we map NTOT
C     back to a Lebedev scheme that exists on the list.
C
         IF(NSCHEME.LT.LEBMIN) THEN
            NSCHEME = LEBMIN
            NTOT = NLEB(NSCHEME)
            NSCHEME = NLEBSET(NTOT)
         END IF 
      ENDIF

      DO K=1,NTOT
         X(K)=0.0D0
         Y(K)=0.0D0
         Z(K)=0.0D0
         WT(K)=0.0D0
      ENDDO

      IF(NSCHEME.GT.64)THEN
         CALL LD1454(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.64.AND.NSCHEME.GT.59) THEN
         CALL LD1454(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.59.AND.NSCHEME.GT.53) THEN
         CALL LD1202(X,Y,Z,WT,NCOUNT) 
      ELSEIF(NSCHEME.LE.53.AND.NSCHEME.GT.47) THEN
         CALL LD0974(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.47.AND.NSCHEME.GT.41) THEN
         CALL LD0770(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.41.AND.NSCHEME.GT.35) THEN
         CALL LD0590(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.35.AND.NSCHEME.GT.31) THEN
         CALL LD0434(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.31.AND.NSCHEME.GT.29) THEN
        CALL LD0350(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.29.AND.NSCHEME.GT.23) THEN
        CALL LD0302(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.23.AND.NSCHEME.GT.21) THEN
        CALL LD0194(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.21.AND.NSCHEME.GT.19) THEN
        CALL LD0170(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.19.AND.NSCHEME.GT.17) THEN
        CALL LD0146(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.17.AND.NSCHEME.GT.15) THEN
        CALL LD0110(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.15.AND.NSCHEME.GT.11) THEN
        CALL LD0086(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.11.AND.NSCHEME.GT.9) THEN
        CALL LD0050(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.9) THEN
        CALL LD0038(X,Y,Z,WT,NCOUNT)
      END IF

      DO J=1,NCOUNT
         WT(J)=WT(J)*4.0D0*PI
      ENDDO

      IF (NTOT.NE.NCOUNT) THEN
         WRITE(LUPRI,*)'Want',NTOT,', but get',ncount
         CALL QUIT('Error 1 in subroutine_SLEB, get Aaron.')
      ENDIF
       
C     Check values of nodes, and weights; calculate absolute errors.
C     WTERROR returns a value of 0 if any of the weights is zero.
C     EABS is the sum of all absolute deviations of point from
C     the surface of the unit sphere. It is then divided by NCOUNT.

      EABS=0.0D0
      WTERROR=1.0D0
      DO I=1,NCOUNT
         XDEV=DABS(1-DSQRT(X(I)**2 + Y(I)**2 + Z(I)**2))
         EABS=EABS+XDEV
         IF (XDEV.GT.1D-8) THEN
            WRITE(LUPRI,*) 'NODE ERROR(1) AT I',I,' deviation ',xdev
            CALL QUIT('Error 2 in subroutine_SLEB, get Aaron.')
         ENDIF
         IF (DABS(WT(I)).LT.1D-5) THEN
            WTERROR=0
            WRITE(LUPRI,*) 'BWGHT ERROR(1) AT I',I,' weight ',wt(I)
            CALL QUIT('Error 3 in subroutine_SLEB, get Aaron.')
         ENDIF
      ENDDO
      EABS=EABS/NCOUNT
      IF (EABS.LT.1D-15) EABS=1.0D-15

C     Check values of nodes, and weights; calculate absolute errors.
C     WTERROR returns a value of 0 if any of the weights is zero.
C     EABS is the sum of all absolute deviations of point from
C     the surface of the unit sphere. It is then divided by NCOUNT.

      EABS=0.0D0
      WTERROR=0.0D0
      DO I=1,NCOUNT
         XDEV=DABS(1.0D0 - DSQRT(X(I)**2 + Y(I)**2 + Z(I)**2))
         EABS=EABS+XDEV
         WTERROR=WTERROR+WT(I)
         IF (XDEV.GT.1D-8) THEN
            WRITE(LUPRI,*) 'NODE ERROR(2) AT I',I,' deviation ',xdev
            CALL QUIT('NODE ERROR')
         ENDIF
      ENDDO
      XDEV=DABS(WTERROR/4.0D0/PI-1.0D0)
      IF (XDEV.GT.1D-9) THEN
         WRITE(LUPRI,*) 'SUM OF WEIGHTS NOT EQUAL TO 1',XDEV
         CALL QUIT('SUM OF WEIGHTS IN LEBEDEV NOT EQUAL TO 1')
      ENDIF
      EABS=EABS/NCOUNT
      IF (EABS.LT.1D-15) EABS=1.0D-15
      
      NSCHEME=NOLD

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bweght */                                                   C
      SUBROUTINE BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,X,Y,Z,WT,
     &                 NAT,NUC,ZAN,C,
     &                 ACCUM,XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,NATOM)
C                                                                      C
C      Written by C. W. Murray                                         C
C      BWGHT calculates the weights associated with                    C
C      the Becke partitioning amongst the atoms.                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 

      DIMENSION RJ(IVECL,NAT),
     &          PSMU(IVECL,NAT),
     &          RIJ(NAT*(NAT-1)/2),
     &          AIJ(NAT*(NAT-1)/2),
     &          X(IVECL),Y(IVECL),Z(IVECL),WT(IVECL),
     &          XMUIJN(IVECL),XMUIJ2(IVECL),ACCUM(IVECL),XPASC(20),
     &          ZAN(*),C(3,*)

C
C     WEIGHT FUNCTION BY BECKE PARTITIONING
C
      DO 100 INA=1,NAT
      DO 100 M=1,IPT
         RJ(M,INA)=DSQRT((C(1,INA)-X(M))**2
     &                  +(C(2,INA)-Y(M))**2
     &                  +(C(3,INA)-Z(M))**2)
         PSMU(M,INA)=1.0D0
  100 CONTINUE

      ITEMP=0
      DO 200 INA=1,NAT
      DO 200 JNA=1,INA-1
         ITEMP=ITEMP+1
         DO 210 M=1,IPT
            XMU=(RJ(M,INA)-RJ(M,JNA))*RIJ(ITEMP)
            XMUIJ=XMU+AIJ(ITEMP)*(1-XMU*XMU)
            XMUIJ2(M)=XMUIJ*XMUIJ
            XMUIJN(M)=XMUIJ
            ACCUM(M)=0.0D0
  210    CONTINUE
         DO 220 I=1,NTRANS+1
         DO 220 M=1,IPT
            ACCUM(M)=ACCUM(M)+XPASC(I)*XMUIJN(M)
            XMUIJN(M)=XMUIJN(M)*XMUIJ2(M)
  220    CONTINUE
         DO 230 M=1,IPT
            PSMU(M,INA)=PSMU(M,INA)*(0.5D0-APASC*ACCUM(M))
            PSMU(M,JNA)=PSMU(M,JNA)*(0.5D0+APASC*ACCUM(M))
  230    CONTINUE
  200 CONTINUE

      ACCUM = 0.0d0
      DO 300 INA=1,NAT
      DO 300 M=1,IPT
         ACCUM(M)=ACCUM(M)+PSMU(M,INA)
  300 CONTINUE
C
C     Contract ACCUM into PSMU
C
      DO 400 INA=1,NAT
      DO 400 M=1,IPT
         PSMU(M,INA)=PSMU(M,INA)/ACCUM(M)
  400 CONTINUE

      DO 500 M=1,IPT
         WT(M)=WT(M)*PSMU(M,NATOM)
  500 CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck cmprsq */                                                   C
      SUBROUTINE CMPRSQ(X,Y,Z,WT,NPOINT,KPOINT)                  
C                                                                      C
C     T. Helgaker                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 

      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)


      J = 0
      DO I = 1, NPOINT
         NDEG = MLTPNT(X(I),Y(I),Z(I))
         IF (NDEG.GT.0) THEN
            J = J + 1 
            X(J) = X(I)
            Y(J) = Y(I)
            Z(J) = Z(I)
            WT(J) = dble(NDEG)*WT(I) 
         END IF
      END DO
      KPOINT = J

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C /* Deck mltpnt */                                                    C
      FUNCTION MLTPNT(PX,PY,PZ)
C                                                                      C
C     Symmetry multiplicity of a point in space                        C
C                                                                      C
C     T. Helgaker  Feb 01                                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"

      PARAMETER (A0 = 1.0D-8)
      LOGICAL LBTAX

#include "symmet.h"
#include "ibtfun.h"

      LBTAX(I,J) = IBTAND(2**(I-1),ISYMAX(J,1)) .GT. 0

      IF (MAXREP.EQ.0) NSYMOP = 0 
      IF (MAXREP.EQ.1) NSYMOP = 1 
      IF (MAXREP.EQ.3) NSYMOP = 2 
      IF (MAXREP.EQ.7) NSYMOP = 3 
      ISTAB = 0
      DO I = 1, NSYMOP
         IF(LBTAX(I,1) .AND. ABS(PX).GT.A0) GOTO 100
         IF(LBTAX(I,2) .AND. ABS(PY).GT.A0) GOTO 100
         IF(LBTAX(I,3) .AND. ABS(PZ).GT.A0) GOTO 100
            ISTAB = ISTAB + 2**(I-1)
 100     CONTINUE 
      END DO
      MLTPNT = MULT(ISTAB)
      DO I = 1, NSYMOP
         IF(LBTAX(I,1) .AND. PX.LT.-A0) MLTPNT = 0
         IF(LBTAX(I,2) .AND. PY.LT.-A0) MLTPNT = 0
         IF(LBTAX(I,3) .AND. PZ.LT.-A0) MLTPNT = 0
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      FUNCTION DISERR(RD,L)
C                                                                      C
C     Provide grid spacing h for given angular momentum L              C
C     and discretization error RD                                      C
C                                                                      C
C     Based on eqs. (17) and (18) of                                   C
C       R. Lindh, P.-Aa. Malmqvist and L. Gagliardi                    C 
C       "Molecular integrals by numerical quadrature",                 C
C       Theor. Chem. Acc. 106 (2001) 178-187                           C
C                                                                      C
C     The array CF(4,L) contains coefficients of a 3rd order           C
C     polynomial fit to provide start values for the                   C
C     determination of H by a Newton-Raphson search.                   C
C                                                                      C
C     Written by T. Saue July 2002                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "pi.h"
      PARAMETER(ACC=1.0D-5)
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,TI=1.0D1)
      PARAMETER(MXIT=20)
      DIMENSION CF(4,0:4)
      DATA CF/0.91570D0,0.78806D-1,0.28056D-2,3.4197D-05,
     &        0.74912D0,0.61502D-1,0.21558D-2,2.6100D-05,
     &        0.65449D0,0.52322D-1,0.18217D-2,2.2004D-05,
     &        0.59321D0,0.46769D-1,0.16261D-2,1.9649D-05,
     &        0.55125D0,0.43269D-1,0.15084D-2,1.8270D-05/
C
C     Initialization
C
      FAC  = SQRT(D2)*D2*D2
      IFAC = 1
      DO I = 1,L
        FAC   = FAC*D2
        IFAC  = IFAC*(2*I+1)
      ENDDO
      FAC = FAC/IFAC
      LM = MIN(L,4)
      RDLOG = LOG(RD)
      DISERR = POLVAL(3,CF(1,LM),RDLOG)
      HTLOG = LOG(DISERR)
C     Newton-Raphson search
      DO IT = 1,MXIT
        PIH  = PI/DISERR
        PIHL = PIH
        PIEX = PI*PIH/D2
        DO I = 1,L
          PIHL = PIHL*PIH
        ENDDO
        U0   = FAC*PIHL*EXP(-PIEX)
        U1   = U0*((PIEX/DISERR)-(L+1)/PIH)
        F0   = LOG(U0)-RDLOG
        F1   = DISERR*U1/U0
        DX = F0/F1
        HTLOG = HTLOG - DX
        DISERR = EXP(HTLOG)
        IF(ABS(DX).LT.ACC) GOTO 10
      ENDDO
      CALL QUIT('Error in DISERR')

 10   RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      FUNCTION OUTERR(AL,L,RD)
C                                                                      C
C     Provide outer grid point for given angular momentum L            C
C     outer exponent AL and discretization error RD                    C
C                                                                      C
C     Based on eq. (19) of                                             C
C       R. Lindh, P.-Aa. Malmqvist and L. Gagliardi                    C
C       "Molecular integrals by numerical quadrature",                 C
C       Theor. Chem. Acc. 106 (2001) 178-187                           C
C                                                                      C
C     The variable U = AL*R*R is found by a Newton-Raphson search.     C 
C                                                                      C
C     Written by T. Saue July 2002                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "pi.h"
      PARAMETER(ACC=1.0D-6)
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,TI=1.0D1)
      PARAMETER(MXIT=8)
C
C     Initialization
C
      TOLEN = D2
      FAC = D1
      DO I = 1,L
        TOLEN = TOLEN*D2
        FAC   = FAC*(2*I+1)
      ENDDO
      EXPL = (2*L+1)/D2
      A = SQRT(PI)*FAC/TOLEN
      ALN = LOG(A)
      RLN = LOG(RD)
      U = 35.0D0
C     Newton-Raphson search
      DO IT = 1,MXIT
        F0HLN = ALN+EXPL*LOG(U)-U-RLN
        F1HLN = EXPL/U-D1
        DX = F0HLN/F1HLN
        U = U - DX
        IF(ABS(DX).LT.ACC) THEN
          OUTERR = SQRT(U/AL)
          GOTO 10
        ENDIF
      ENDDO
      CALL QUIT('Error in OUTERR')

 10   RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1wt  */
      SUBROUTINE GETCHI(DMAT,CHIMT,NAT,C,IPRINT)
C***********************************************************************
C
C     Get relative atomic size estimates by analyzing atomic
c     contributions to the density along bonds. 
C
C     Written by Olav Fossgaard - June 2003
C
C***********************************************************************

      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION DMAT(*)
      DIMENSION C(3,NAT)
#include "dcbgen.h"
C
      DIMENSION CHIMT(NAT,NAT)
      real(8), allocatable :: kcoor(:)
      real(8), allocatable :: kdmtao(:)
      real(8), allocatable :: kgao(:)
      real(8), allocatable :: ksmat(:)
      real(8), allocatable :: kasmat(:)
      integer, allocatable :: kiclab(:)
C
#include "dcbbas.h"
#include "dcbdhf.h"
C
      CALL QENTER('GETCHI')
C
C
C     Quit if open shell
C
      IF (NFMAT .GT. 1) THEN
        WRITE(LUPRI,'(/A)')
     &  '--> Not implemented for open shell yet, nothing done.'
        CALL QUIT('GETCHI: Not implemented for open shell')
      END IF
C
C     Memory allocation
C
      call alloc(KCOOR  , 6                 )
      call alloc(KDMTAO , N2BBASX           )
      call alloc(KGAO   , NTBAS(0)          )
      call alloc(KSMAT  , NTBAS(1)*NTBAS(1) )
      call alloc(KASMAT , NAT*NAT           )
      call alloc(KICLAB , NTBAS(0)          )
C
C     Transform density matrix to unsorted basis, then to AO basis
C
      CALL BSTOBU_no_work(DMAT,1)
      CALL DTSOAO(DMAT,KDMTAO,NTBAS(0),0,IPRINT)
C
C     Get atomic size estimates
C
      CALL CHI1(CHIMT,NAT,C,KDMTAO,KGAO,KSMAT,
     &     KASMAT,KICLAB,KCOOR,IPRINT)
C
C     Backtransform density matrix to sorted basis
C
      CALL BUTOBS_no_work(DMAT,1)
C
C     Memory deallocation
C
      call dealloc(kcoor)
      call dealloc(kdmtao)
      call dealloc(kgao)
      call dealloc(ksmat)
      call dealloc(kasmat)
      call dealloc(kiclab)
C
 9999  CALL QEXIT('GETCHI')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1rh1 */
      SUBROUTINE CHI1(CHIMT,NAT,C,DMAT,GAO,SMAT,ASMAT,ICLAB,COORD,
     &     IPRINT)
C
C     Get atomic size estimates by analyzing pairwise partial densities
C     in one dimension, along bonds.
C
C     Written by Trond Saue May 9 2000
C     Modified for CHI1 by Olav Fossgaard - June 2003
C
C***********************************************************************
#include "implicit.h"

      PARAMETER(MAXCEN = 100)

#include "priunit.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "dcbbas.h"
      DIMENSION C(3,NAT)
#if defined (VAR_G77) || defined (SYS_SX)
C     g77 fails in the complicated dimension statement that was used
C     in the original version, make things easier for the compiler with an
C     equivalence statment
      EQUIVALENCE (NBAST,NTBAS(0)),(NBASL,NTBAS(1))
      DIMENSION COORD(3,2),CHIMT(NAT,NAT),DMAT(NBAST,NBAST),
     &     GAO(NBAST),SMAT(NBASL,NBASL),ASMAT(NAT,NAT),
     &     ICLAB(NBAST)
#else
      DIMENSION COORD(3,2),CHIMT(NAT,NAT),DMAT(NTBAS(0),NTBAS(0)),
     &     GAO(NTBAS(0)),SMAT(NTBAS(1),NTBAS(1)),ASMAT(NAT,NAT),
     &     ICLAB(NTBAS(0))
#endif
#include "ibtfun.h"
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(8X,A)') 'X  : Coordinate along the bond' 
      WRITE(LUPRI,'(8X,A)') 'Y  : Sum of atomic contributions to '//
     &     'the large component density'
      WRITE(LUPRI,'(8X,A)') 'RHO: Total large component density '
     
      DO 10 I = 1,NAT
         DO 20 J = 1,I-1
            CHIMT(I,J) = 1.0D0
            COORD(1,1) = C(1,I)
            COORD(2,1) = C(2,I)
            COORD(3,1) = C(3,I)
            COORD(1,2) = C(1,J)
            COORD(2,2) = C(2,J)
            COORD(3,2) = C(3,J)
            CALL CHI2(CHIMT,NAT,DMAT,GAO,SMAT,ASMAT,ICLAB,I,J,COORD,
     &              IPRINT)
 20      END DO
 10   END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CHI2(CHIMT,NAT,DMAT,GAO,SMAT,ASMAT,ICLAB,ICENTA,ICENTB,
     &     COORD,IPRINT)
C                                                                      C
C     Get atomic size estimates by analyzing pairwise partial densities
C     in one dimension, along bonds.

C     Written by Olav Fossgaard - June 2003
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "consts.h"
#include "dcbbas.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
C     g77 fails in the complicated dimension statement that was used
C     in the original version, make things easier for the compiler with an
C     equiavelence statment
      EQUIVALENCE (NBAST,NTBAS(0)),(NBASL,NTBAS(1))
      DIMENSION CHIMT(NAT,NAT),COORD(3,2),DMAT(NBAST,NBAST),
     &     GAO(NBAST),SMAT(NBASL,NBASL),ASMAT(NAT,NAT),
     &     ICLAB(NBAST)
#else
      DIMENSION CHIMT(NAT,NAT),COORD(3,2),DMAT(NTBAS(0),NTBAS(0)),
     &     GAO(NTBAS(0)),SMAT(NTBAS(1),NTBAS(1)),ASMAT(NAT,NAT),
     &     ICLAB(NTBAS(0))
#endif
C
C     Arrays for 2nd order polynomial fit:
C
      DIMENSION A(3,3),B(3),COFF(3),D(3),X(3),Y(3)
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A,I2,A,I2,A)') '** Region of density minimum '//
     &      'for atom pair (', ICENTA,',', ICENTB,')'
      WRITE(LUPRI,*)
C
C     Calculate distance between centers
C
      DSTEP = 0.01D0

      DX = (COORD(1,2)-COORD(1,1))
      DY = (COORD(2,2)-COORD(2,1))
      DZ = (COORD(3,2)-COORD(3,1))
      DISTAN = SQRT(DX*DX+DY*DY+DZ*DZ)
      DX     = DX/DISTAN
      DY     = DY/DISTAN
      DZ     = DZ/DISTAN
      ASTEP  = DSTEP/XTANG
      NPOINT = INT(DISTAN/ASTEP)
      IF(MOD(NPOINT,2).EQ.1) NPOINT = NPOINT+1
      ASTEP  = DISTAN/dble(NPOINT)

      RHOA   = 0.0D0
      RHOB   = 0.0D0
      TOTRHO = 0.0D0
      SUM    = 0.0D0
      SUMOLD = 0.0D0
C
C     Start between nuclei and get sum of partial densities for
C     atom pair (ICENTA,ICENTB).

      IP  = NPOINT/2
      PP  = dble(IP)*ASTEP
      PX  = COORD(1,1) + PP*DX
      PY  = COORD(2,1) + PP*DY
      PZ  = COORD(3,1) + PP*DZ
C
      CALL CHI3(PX,PY,PZ,RHOA,RHOB,TOTRHO,GAO,DMAT,ICENTA,ICENTB,SMAT,
     &     ASMAT,ICLAB,NAT,IPRINT)

      PPOLD = PP
      SUMOLD = RHOA +RHOB
C
C     Take one step towards ICENTB
C
      IP  = IP + 1
      PP  = dble(IP)*ASTEP
      PX  = COORD(1,1) + PP*DX
      PY  = COORD(2,1) + PP*DY
      PZ  = COORD(3,1) + PP*DZ
C
      CALL CHI3(PX,PY,PZ,RHOA,RHOB,TOTRHO,GAO,DMAT,ICENTA,ICENTB,SMAT,
     &     ASMAT,ICLAB,NAT,IPRINT)

      SUM = RHOA + RHOB
C
C     Decide which way to go. Start walking until a minimum is found
C     on the grid.
C
      IF(SUM.LT.SUMOLD) THEN
         ISTART = NPOINT/2 + 1
         IEND   = NPOINT
         INCR   =  1
         IBAK   = -1
      ELSE
         ISTART = NPOINT/2 - 1
         IEND   =  0
         INCR   = -1
         IBAK   =  1
      ENDIF

      DO 1 IP = ISTART,IEND,INCR
         PP  = dble(IP)*ASTEP
         PX  = COORD(1,1) + PP*DX
         PY  = COORD(2,1) + PP*DY
         PZ  = COORD(3,1) + PP*DZ
C
         CALL CHI3(PX,PY,PZ,RHOA,RHOB,TOTRHO,GAO,DMAT,ICENTA,ICENTB,
     &        SMAT,ASMAT,ICLAB,NAT,IPRINT)

         SUM = RHOA + RHOB
         IF(SUM.LE.SUMOLD) THEN
            SUMOLD = SUM
            PPOLD = PP
            GO TO 1
         ELSEIF(SUM.GT.SUMOLD) THEN
            IPMIN = IP + IBAK
            PPMIN = PPOLD
            SUMMIN = SUMOLD
            GO TO 2
         END IF
 1    CONTINUE
C
C     Fit a second order polynomial around minimum point on the grid
C
 2    K = 1
      DO IP  = IPMIN-1,IPMIN+1
         PP  = dble(IP)*ASTEP
         PX  = COORD(1,1) + PP*DX
         PY  = COORD(2,1) + PP*DY
         PZ  = COORD(3,1) + PP*DZ
C     
         CALL CHI3(PX,PY,PZ,RHOA,RHOB,TOTRHO,GAO,DMAT,ICENTA,ICENTB,
     &        SMAT,ASMAT,ICLAB,NAT,IPRINT)
C
         X(K) = PP
         Y(K) = RHOA + RHOB
C
         WRITE(LUPRI,'(8X,A,I1,A,F12.8,4X,A,I1,A,F12.8,4X,A,F12.8)') 
     &        'X(',K,') = ',X(K),'Y(',K,') =',Y(K),'RHO = ',
     &        TOTRHO
C
         K = K +1
      ENDDO

      CALL POLSVD(3,3,A,B,X,Y,COFF,D,CHISQ,ISKIP)
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(8X,A)') 'Coefficients for 2nd order '//
     &     'interpolating polynomial:'
      WRITE(LUPRI,'(12X,A,I1,A,1P,E14.6)') ('c(',(I-1),'):  ',
     &      COFF(I),I=1,3)
      WRITE(LUPRI,*)

C
C     Minimum of 2. order interpolating polynomial is -b/2a
c
      PPMIN = -COFF(2)/2.0D0/COFF(3)
      
      WRITE(LUPRI,'(8X,A,I1,A,F12.8,A)') 'Radius of atom ',ICENTA,': ',
     &     PPMIN*XTANG,' A'
      WRITE(LUPRI,'(8X,A,I1,A,F12.8,A)') 'Radius of atom ',ICENTB,': ',
     &     (DISTAN - PPMIN)*XTANG, ' A'

      CHIMT(ICENTA,ICENTB) = PPMIN/(DISTAN - PPMIN)
      CHIMT(ICENTB,ICENTA) = 1.0D0/CHIMT(ICENTA,ICENTB)

      RETURN
C
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      SUBROUTINE CHI3(X,Y,Z,RHOA,RHOB,TOTRHO,GAO,DMAT,ICENTA,ICENTB,
     &     SMAT,ASMAT,ICLAB,NAT,IPRINT)

C     Get atomic size estimates by analyzing pairwise partial densities
C     in one dimension, along bonds.
C                                                                      C
C      Written by Olav Fossgaard - June 2003                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"

      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
C
#if defined (VAR_G77) || defined (SYS_SX)
C     g77 fails in the complicated dimension statement that was used
C     in the original version, make things easier for the compiler with an
C     equiavelence statment
      EQUIVALENCE (NBAST,NTBAS(0)),(NBASL,NTBAS(1))
      DIMENSION GAO(NBAST),DMAT(NBAST,NBAST),
     &     SMAT(NBASL,NBASL),ASMAT(NAT,NAT),
     &     ICLAB(NBAST)
#else
      DIMENSION GAO(NTBAS(0)),DMAT(NTBAS(0),NTBAS(0)),
     &     SMAT(NTBAS(1),NTBAS(1)),ASMAT(NAT,NAT),
     &     ICLAB(NTBAS(0))
#endif
      GAO   = 0.0d0
      SMAT  = 0.0d0
      ASMAT = 0.0d0
      iclab = 0
C
C     GET AO's
C
      CALL GETAOS(GAO,TEMP,TEMP,TEMP,TEMP,X,Y,Z,NTBAS(0),
     &     0,.FALSE.,IPRINT)
C
C     Get "overlap" matrix
C
      CALL DGEMM('N','T',NTBAS(1),NTBAS(1),1,ONE,GAO,NTBAS(1),
     &     GAO,NTBAS(1),ZERO,SMAT,NTBAS(1))
C
C     Multiply with density matrix elements
C
      DO I = 1,NTBAS(1)
         DO J = 1,NTBAS(1)
            SMAT(I,J) = SMAT(I,J)*0.5D0*DMAT(I,J)
         END DO
      END DO
C
C     Map basis functions onto atomic centers
C
      CALL LSQMA1(ICLAB)
C
C     Group SMAT matrix elements onto atomic centers. Overlap
c     elements are ignored.
C
      CALL GATMAT(1,NTBAS(1),NTBAS(1),SMAT,ASMAT,ICLAB,NAT)
C
C     Calculate partial densities.
C
      RHOA = ASMAT(ICENTA,ICENTA)
      RHOB = ASMAT(ICENTB,ICENTB)
C
C     Calculate total density for comparison with RHO1
C
      TOTRHO = 0.0D0
      DO I = 1,NAT
         DO J = 1,NAT
            TOTRHO = TOTRHO + ASMAT(I,J)
         END DO
      END DO
      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bragg */
      FUNCTION BRAGG(NCHARGE)
#include "implicit.h"
C     Returns Slater-Bragg radius in Angstrom;
C     Returns radius of MAXCHR when data not available.
C 
C     The below data gives atomic radii in Angstroms and stems from table I of 
C     J.C.Slater: "Atomic Radii in Crystals"
C     J.Chem.Phys. 41(1964) 3199-3204
C     Values for elements marked with an asterisk has been
C     guessed/interpolated
C
      PARAMETER (MAXCHR=103)
      DIMENSION RAD(0:MAXCHR)
      DATA (RAD(I), I = 0, MAXCHR)/0.75D0,
C       H      He*   
     &  0.35D0,0.35D0,
C       Li     Be     B      C      N      O      F      Ne*
     &  1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
C       Na     Mg     Al     Si     P      S      Cl     Ar*
     &  1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0,
C       K      Ca     Sc     Ti     V      Cr     Mn     Fe     Co 
     &  2.20D0,1.80D0,1.60D0,1.40D0,1.35D0,1.40D0,1.40D0,1.40D0,1.35D0,
C       Ni     Cu     Zn     Ga     Ge     As     Se     Br     Kr*
     &  1.35D0,1.35D0,1.35D0,1.30D0,1.25D0,1.15D0,1.15D0,1.15D0,1.10D0,
C       Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru     Rh
     &  2.35D0,2.00D0,1.80D0,1.55D0,1.45D0,1.45D0,1.35D0,1.30D0,1.35D0,
C       Pd     Ag     Cd     In     Sn     Sb     Te     I      Xe*
     &  1.40D0,1.60D0,1.55D0,1.55D0,1.45D0,1.45D0,1.40D0,1.40D0,1.40D0,
C       Cs     Ba     La    
     &  2.60D0,2.15D0,1.95D0,
C       Ce     Pr     Nd     Pm     Sm     Eu     Gd
     &  1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.80D0,
C       Tb     Dy     Ho     Er     Tm     Yb     Lu
     &  1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
C       Hf     Ta     W      Re     Os     Ir     Pt     Au     Hg
     &  1.55D0,1.45D0,1.35D0,1.30D0,1.30D0,1.35D0,1.35D0,1.35D0,1.50D0,
C       Tl     Pb*    Bi     Po     At*    Rn*
     &  1.90D0,1.75D0,1.60D0,1.90D0,1.50D0,1.50D0,
C       Fr*    Ra     Ac     
     &  2.15D0,2.15D0,1.95D0,
CTROND rad(U): 1.75 --> 1.37D0
C       Th     Pa     U      Np     Pu     Am     Cm*     
     &  1.80D0,1.80D0,1.37D0,1.75D0,1.75D0,1.75D0,1.75D0,
CTROND       Th     Pa     U      Np     Pu     Am     Cm*     
CTROND     &  1.80D0,1.80D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
C       Bk*    Cf*    Es*    Fm*    Md*    No*    Lw*
     &  1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0/ 
C
      IF(NCHARGE.GT.MAXCHR) THEN
        BRAGG = RAD(MAXCHR)
      ELSEIF(NCHARGE.LT.0) THEN
        BRAGG = -1.0D0
      ELSE
        BRAGG = RAD(NCHARGE)
      ENDIF
      RETURN
      END


      SUBROUTINE SET_GRID_INFO(NAT,NUM,ZAN,C)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"
      DIMENSION ZAN(*), C(3,NUCDEP)
#include "ibtfun.h"
      NAT = 0
      NUM = NBASIS
      DO ICENT = 1, NUCIND
        MULCNT = ISTBNU(ICENT)
        ZCNT   = CHARGE(ICENT)
        IF(NAMN(ICENT)(1:2).EQ.'Gh') CYCLE
        IF(NOORBT(ICENT)) CYCLE
        DO ISYMOP = 0, MAXOPR
        IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
          NAT = NAT + 1
          ZAN(NAT) = ZCNT
          C(1,NAT)= PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
          C(2,NAT)= PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
          C(3,NAT)= PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
        ENDIF
        ENDDO
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* ZIP_GRID */
      SUBROUTINE ZIP_GRID(IATOM,X,Y,Z,W,NR_POINTS)
C=======================================================================
C                                                           radovan bast
C                                           last revision: february 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

#include "nuclei.h"


#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"

      logical :: lsymop(0:7)
      DIMENSION X(*),
     &          Y(*),
     &          Z(*),
     &          W(*)

      lsymop = .false.

      do i = 1, maxrep
        lsymop(jsop(i)) = .true.
      end do

      IF(lsymop(1)) THEN
        IF(CORD(1,IATOM) .EQ. D0) THEN
          CALL ZIP_GRID_SUB(X,Y,Z,W,NR_POINTS)
        ENDIF
      ENDIF
      IF(lsymop(2)) THEN
        IF(CORD(2,IATOM) .EQ. D0) THEN
          CALL ZIP_GRID_SUB(Y,Z,X,W,NR_POINTS)
        ENDIF
      ENDIF
      IF(lsymop(4)) THEN
        IF(CORD(3,IATOM) .EQ. D0) THEN
          CALL ZIP_GRID_SUB(Z,X,Y,W,NR_POINTS)
        ENDIF
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* ZIP_GRID_SUB */
      SUBROUTINE ZIP_GRID_SUB(X,Y,Z,W,NR_POINTS)
C=======================================================================
C                                                           radovan bast
C                                           last revision: february 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION X(*),
     &          Y(*),
     &          Z(*),
     &          W(*)

      IWRITE = 1
      DO IREAD = 1,NR_POINTS
        IF(X(IREAD) .GT. D0) THEN
          X(IWRITE) = X(IREAD)
          Y(IWRITE) = Y(IREAD)
          Z(IWRITE) = Z(IREAD)
          W(IWRITE) = W(IREAD)
          IWRITE    = IWRITE + 1
        ELSEIF(X(IREAD) .LT. D0) THEN
!         drop these points
        ELSE
          X(IWRITE) = X(IREAD)
          Y(IWRITE) = Y(IREAD)
          Z(IWRITE) = Z(IREAD)
          W(IWRITE) = W(IREAD)/D2
          IWRITE    = IWRITE + 1
        ENDIF
      ENDDO
      NR_POINTS = IWRITE - 1

      CALL DSCAL(NR_POINTS,D2,W,1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck Lebedev */
      SUBROUTINE LEBEDEV(NSCHEME,X,Y,Z,WT,NCOUNT)
C***********************************************************************
C     Wrapper routine for Lebedev quadrature
C
C***********************************************************************      
#include "implicit.h"
      DIMENSION X(*),Y(*),Z(*),WT(*)
      IF(NSCHEME.GT.131)THEN
         WRITE(6,*) 'Requested accuracy:        ',NSCHEME
         WRITE(6,*) 'Maximum accuracy available:',131
         CALL QUIT('LEBEDEV: Requested accuracy not available !')
      ELSEIF(NSCHEME.LE.131.AND.NSCHEME.GT.125) THEN
         CALL LD5810(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.125.AND.NSCHEME.GT.119) THEN
         CALL LD5294(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.119.AND.NSCHEME.GT.113) THEN
         CALL LD4802(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.113.AND.NSCHEME.GT.107) THEN
         CALL LD4334(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.107.AND.NSCHEME.GT.101) THEN
         CALL LD3890(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.101.AND.NSCHEME.GT.95) THEN
         CALL LD3470(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.95.AND.NSCHEME.GT.89) THEN
         CALL LD3074(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.89.AND.NSCHEME.GT.83) THEN
         CALL LD2702(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.83.AND.NSCHEME.GT.77) THEN
         CALL LD2354(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.77.AND.NSCHEME.GT.71) THEN
         CALL LD2030(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.71.AND.NSCHEME.GT.65) THEN
         CALL LD1730(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.65.AND.NSCHEME.GT.59) THEN
         CALL LD1454(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.59.AND.NSCHEME.GT.53) THEN
         CALL LD1202(X,Y,Z,WT,NCOUNT) 
      ELSEIF(NSCHEME.LE.53.AND.NSCHEME.GT.47) THEN
         CALL LD0974(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.47.AND.NSCHEME.GT.41) THEN
         CALL LD0770(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.41.AND.NSCHEME.GT.35) THEN
         CALL LD0590(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.35.AND.NSCHEME.GT.31) THEN
         CALL LD0434(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.31.AND.NSCHEME.GT.29) THEN
        CALL LD0350(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.29.AND.NSCHEME.GT.23) THEN
        CALL LD0302(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.23.AND.NSCHEME.GT.21) THEN
        CALL LD0194(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.21.AND.NSCHEME.GT.19) THEN
        CALL LD0170(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.19.AND.NSCHEME.GT.17) THEN
        CALL LD0146(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.17.AND.NSCHEME.GT.15) THEN
        CALL LD0110(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.15.AND.NSCHEME.GT.11) THEN
        CALL LD0086(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.11.AND.NSCHEME.GT.9) THEN
        CALL LD0050(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.9.AND.NSCHEME.GT.7) THEN
        CALL LD0038(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.7.AND.NSCHEME.GT.5) THEN
        CALL LD0026(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.5.AND.NSCHEME.GT.3) THEN
        CALL LD0014(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.3) THEN
        CALL LD0006(X,Y,Z,WT,NCOUNT)
      END IF
      RETURN
      END
