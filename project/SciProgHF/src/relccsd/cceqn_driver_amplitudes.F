      SUBROUTINE cceqn_driver_amplitudes(MAXDIM,MAXIT,EPS,FOO,FVO,FVV,
     &                 T1,S1,T2,S2,AZ,H,HO,HV,HOV,GO,GV,
     &                 BUF1,BUF2,BUF3,CONV,NBUF2,NBUF3)

      use interface_to_mpi
      use xmlout
      implicit none

C---------------Description--------------------------------------------
C
C     Solves CCSD amplitude or lambda equations.
C
C---------------Authors------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 CONV                ! Desired convergence
      INTEGER MAXIT              ! Maximum number of iterations that is allowed
      INTEGER MAXDIM             ! Maximum size of DIIS expansion space
      INTEGER NBUF2,NBUF3        ! Size of the buffer arrays

C     Pre-allocated arrays used for intermediaries etc.
      REAL*8 EPS(*),FOO(*),FVO(*),FVV(*),T1(*),S1(*),T2(*),S2(*)
      REAL*8 AZ(*),H(*),HO(*),HV(*),HOV(*),GO(*),GV(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)


C---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "ihm.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      integer ISTAT(df_MPI_STATUS_SIZE)
#endif

C---------------Local variables--------------------------------------

      REAL*8 RMS                             ! Actual convergence
      real(8), parameter :: rms_limit = 1.D4 ! divergence criterium
      REAL*8 DIAG                            ! T1-diagnostic 
      REAL*8 CPU,CPU4V,CPUH1,CPUH2           ! used in timing
      LOGICAL     FOCKSP,PRENE
      CHARACTER*2 DELS1                      ! specifies S1 blocks to be deleted in denomf
      CHARACTER*4 DELS2                      ! specifies S2 blocks to be deleted in denomf
      REAL*8 DDOT
      Complex(8) ZDOTC
      INTEGER it,itc,nit,IERR

!     Arrays needed by DIIS routine
      REAL*8  CC1(MAXDIM),BB1(MAXDIM+1,MAXDIM+2),BB2(MAXDIM+1,MAXDIM+2)

C---------------Executable code--------------------------------------
C
C
C     Get initial values of timers for the most important contractions
C     This is important in Fock space calculations where we enter this
C     routine more than once. We do it in this manner rather than resetting
C     the timer because we want the overall timing at the end.

      CALL XTIMRD (CPU4V,'-- BINTM                       ')
      CALL XTIMRD (CPUH1,'--- HINTM: VVOO contribution   ')
      CALL XTIMRD (CPUH2,'--- T2EQNS: HINTM*T2           ')

      CALL cceqn_driver_init(MAXIT,FVO,T1,S1,T2,S2,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,
     &                 NBUF2,NBUF3,RMS,FOCKSP,DELS1,DELS2)

C
C  START ITERATIVE PROCEDURE
C
      IT = 0
      ITC = 0
      NIT = 0

      do
      NIT = NIT + 1
      IF (NIT.GT.MAXIT) THEN
         maxit = maxit + 1
         NIT = NIT - 1
         WRITE(IW,6004) RMS
         exit
      ENDIF
C
C  EVALUATE T1 EQUATIONS
C  CHECK FOR CASES WITH NO NON_ZERO T1s TO AVOID CRASHES
C
      CALL XTIME(5,1,'- T1 equations                  ')
      CALL HOINTM (FOO,EPS,T1,T2,BUF1,BUF2,BUF3,HO)
      CALL HVINTM (FVV,EPS,T1,T2,BUF1,BUF2,BUF3,HV)
C      print*,"T1 amplitudes ",zdotc(ndimt1,t1,1,t1,1)
C      print*,"T2 amplitudes ",zdotc(ndimt2,t2,1,t2,1)
C      print*,"TAU intermediate ",zdotc(NO(1)*NVVOT(1),BUF2,1,BUF2,1)
           
      IF (NDIMT1.GT.0) THEN
          CALL HOVINTM (FVO,T1,BUF1,BUF2,HOV)
C          print*,"HOV intermediate ",zdotc(nfvo,hov,1,hov,1)
          CALL amplitude_equation_T1 (FVO,T1,T2,HO,HV,HOV,
     &                 BUF1,BUF2,BUF3,NBUF3,S1)
      ENDIF
      CALL XTIME(5,2,'- T1 equations                  ')

C
C     In case of NOCCD we can skip the calculation of T2 amplitudes
C
      IF (        NOCCD) GOTO 100
           
C
C  EVALUATE T2 INTERMEDIATES
C
      CALL XTIME(5,1,'- T2 equations                  ')
      CALL XTIME(0,1,'-- GOINTM                        ')
      CALL XCOPY (NFOO,HO,1,GO,1)
      CALL GOINTM (FVO,T1,BUF1,BUF2,GO)
      CALL XTIME(0,2,'-- GOINTM                        ')
      CALL XTIME(0,1,'-- GVINTM                        ')
      CALL XCOPY (NFVV,HV,1,GV,1)
      CALL GVINTM (FVO,T1,BUF3,NBUF3,BUF2,NBUF2,GV)
      CALL XTIME(0,2,'-- GVINTM                        ')
      CALL XTIME(0,1,'-- AINTM                         ')
      CALL AINTM (T1,T2,BUF1,BUF2,AZ)
      CALL XTIME(0,2,'-- AINTM                         ')
      CALL XTIME(0,1,'-- HINTM                         ')
      CALL HINTM (T1,T2,BUF1,BUF2,BUF3,NBUF3,H)
      CALL XTIME(0,2,'-- HINTM                         ')
C
C  INTERMEDIATES FORMED --- START T2 EQUATIONS
C
      CALL XTIME(5,1,'-- T2 EQNS                     ')
      CALL amplitude_equation_T2 (T1,T2,GO,GV,HOV,H,AZ,
     &        BUF1,BUF2,BUF3,NBUF3,S2)
      CALL XTIME(5,2,'-- T2 EQNS                     ')
C
 100  CONTINUE
C
C     Delete the Pi->Q excitations in case of Int. Ham. CCFS scheme 2
      IF (FOCKSP.AND.DOIH.AND.IHSCHEME.EQ.2) THEN
         CALL DENOMF (EPS,S1,S2,S1,S2,3)
      ENDIF
C
C     Add the first order wave function contribution
C
      CALL S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)
C
C  ADD PARTIAL S1 AND S2 (PARALLEL RUNS)
C
      CALL XTIME(4,1,'- adding partial T1/T2 amplitudes ')
      CALL PARTS (1,RMS,S1,S2,BUF1,BUF2)
      CALL XTIME(4,2,'- adding partial T1/T2 amplitudes ')
C
      IF (FOCKSP) THEN
C       Adding folded diagrams
        CALL FOLDED(ICURFSS,T1,T2,EPS,FVO,
     &              BUF1,BUF2,BUF3,S1,S2)
      ENDIF
C
C  Make new T1s and T2 by dividing S1s and S2 by eigenvalues
C
      CALL XCOPY (NDIMT1,S1,1,T1,1)
      CALL XCOPY (NDIMT2,S2,1,T2,1)
C
      IF (FOCKSP) THEN
C       Remove wave function amplitudes from higher sectors
        CALL DELFCK ('VO',DELS1,T1)
        CALL DELINT ('VVOO',DELS2,T2,0,0,0)
C       Scale with denominators
        CALL DENOMF (EPS,T1,T2,T1,T2,1)
C       Remove the Hamiltonian elements, T1 and T2 should only contain
C       the wave function (amplitudes) appropriate for this sector.
        CALL DENOMF (EPS,T1,T2,T1,T2,2)
      ELSE
        CALL DENOM (EPS,T1,T2,T1,T2)
      ENDIF
C
C     Zero amplitudes if the NOCCX options are specified
C
      IF (NOCCS) CALL XCOPY (NDIMT1,A0,0,T1,1)
      IF (NOCCD) CALL XCOPY (NDIMT2,A0,0,T2,1)
C
      CALL ZCORE (T1,T2)
      !mlv
      CALL FCORE(T1,T2,EPS)
C
C  EXTRAPOLATE THE AMPLITUDES
C
      CALL XTIME(5,1,'- DIIS extrapolation            ')
      CALL DIISX(NDIMT1,NDIMT2,T1,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
      CALL XTIME(5,2,'- DIIS extrapolation            ')
C
C  SYNCHRONIZE T1 AND T2 (PARALLEL RUNS)
C
      CALL XTIME(4,1,'- synchronizing T1 & T2 amplitudes')
      CALL PARTS (2,RMS,T1,T2,BUF1,BUF2)
      CALL XTIME(4,2,'- synchronizing T1 & T2 amplitudes')
C
C  EVALUATE CORRELATION ENERGY AND WRITE OUT RESULTS
C
      IF (IPRNT.GE.3) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
      CALL CCENG  (T1,T2,FVO,BUF1,BUF2,ECCSD,ECCSDIM)
C
      CALL T1DIAG (DIAG,T1,BUF1)
      IF (IPRNT.GE.1) WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      CALL FLSHFO(IW)
C
      IF(RMS.LT.CONV) exit
      IF(RMS.GT.rms_limit) THEN
        WRITE (IW,*) ' WARNING: aborted CC iterations (divergence)'
        exit
      END IF

      end do
C
      CALL T1DIAG (DIAG,T1,BUF1)
      WRITE(IW,6071)
C
C     Only print energy when doing amplitudes for the reference
C     wave function.
C
      IF (.NOT.FOCKSP) THEN
         PRENE =  .TRUE.
      ELSE
         PRENE = ICURFSS .EQ. 1
      ENDIF
C
      IF (PRENE) THEN
         ETOT = ESCF + ECCSD
         WRITE(IW,6072) ESCF
         IF (.NOT.(NOCCD.OR.NOCCS)) THEN
            WRITE(IW,6074) ECCSD
         ELSE
            IF (NOCCS) WRITE(IW,6084) ECCSD
            IF (NOCCD) WRITE(IW,6085) ECCSD
         ENDIF
         WRITE(IW,6076) ETOT
         IF (IPRNT.GE.2) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)
      ENDIF
C
      WRITE(IW,6077) DIAG
      WRITE(IW,6078) RMS
      WRITE(IW,6079) NIT
      if (doxml) then
         call xml_quantity('CC iterations used',NIT)
         call xml_quantity('CC convergence',RMS,' ')
         call xml_quantity('T1 diagnostic',DIAG,' ')
      end if
C
C     Read out timers for the most important contractions
C
      CALL XTIMRD (CPU,'-- BINTM                       ')
      CPU4V = CPU - CPU4V
      CALL XTIMRD (CPU,'--- HINTM: VVOO contribution   ')
      CPUH1 = CPU - CPUH1
      CALL XTIMRD (CPU,'--- T2EQNS: HINTM*T2           ')
      CPUH2 = CPU - CPUH2
C
C     Compute & print the performance in Mega/Giga/TeraFlops
C
      IF (TIMING) CALL BFLOPS (CPU4V,CPUH1,CPUH2,NIT)
C
C     Close and delete DIIS scratch files
C
      IT = -2
      CALL DIISX(NDIMT1,NDIMT2,T1,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)
C
 111  FORMAT(A20,I4,2f25.15)
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6002 FORMAT(3X,'MP2   =',F20.15/3X,'EMP2  =',F20.15/)
 6003 FORMAT(//3X,'NIT',10X,'ENERGY',20X,'RMS',7X,'T1-DIAGN')
 6004 FORMAT(//' WARNING: Requested convergence not reached in CCSD :',
     &         'RMS = ',F20.15)
 6071 FORMAT(//,'  CCSD results',/)
 6072 FORMAT(' SCF energy :',T40,F25.15)
 6073 FORMAT(' MP2 correlation energy :',T40,F25.15)
 6074 FORMAT(' CCSD correlation energy :',T40,F25.15)
 6084 FORMAT(' CCD correlation energy :',T40,F25.15)
 6085 FORMAT(' CCS correlation energy :',T40,F25.15)
 6075 FORMAT(' Total MP2 energy :',T40,F25.15)
 6076 FORMAT(' Total CCSD energy :',T40,F25.15)
 6077 FORMAT(' T1 diagnostic :',T40,F25.15)
 6078 FORMAT(' Convergence :',T40,F25.15)
 6079 FORMAT(' Number or iterations used :',T60,I5)
      RETURN
      END
