C  /* Deck getbed */                                                  
      SUBROUTINE GET_BED(OSCBED,OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Build representation matrix of full light-matter interaction
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC, to be opened with unit number LUXVC
C     Reference:
C       N.H. List, T.R.L Melin, M. van Horn, T. Saue, J. Chem. Phys. 152 (2020) 184110
C                                                                      
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,DP5=0.5D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),OSCBED(MAXEXC),
     &          WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      real(8), allocatable              :: TMOM(:,:,:)
      character (len=8),  allocatable   :: LABINT(:)
      character (len=5)                 :: TPAR(3)
      integer, allocatable              :: INTREP(:),IPSYM2INT(:,:)
      integer, allocatable              :: IVECS(:)
      real(8)                   :: PBUF(3,2),EVECR(NREDM,*)
      integer                   :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      logical, allocatable              :: first(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
      NOP = 2
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(8))
      allocate(INTREP(8))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      CALL PR1DIR('CXIKR   ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C      
C     Find what symmetries are available
C
      allocate(IPSYM2INT(NBSYM,NOP))
      NCOMP=NPRPINT/2
      CALL ICOPY(NOP*NBSYM,-1,0,IPSYM2INT,1)
      DO I = 1,NCOMP
        ISYM = INTREP(I) + 1
        IPSYM2INT(ISYM,1) = I
        II = NCOMP + I
        ISYM = INTREP(II) + 1
        IPSYM2INT(ISYM,2) = II
      ENDDO
C
C     Operator setup
C
      IREP   =  JSYMOP-1
      ITIM   = -JTIMOP  ! Time antisymmetric operator
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      IOPTYP  = 8                ! iA.P           * dot-product of alpha and vector operator
      NMAT    = MCMP(IOPTYP)     ! number of components, should be three
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &   ' ** ERROR in GET_BED **  Redimension MAXPRPS to ',NPRPS+NOP
        CALL QUIT('GET_BED: Too may properties !')
      ENDIF
      DO IOP = 1,NOP
        WRITE(PRPNAM(NPRPS+IOP)(1:8),'(A3,I1,A1,A3)')
     &        'BED',IOP,'_',REP(IREP)
        IPRPSYM(NPRPS+IOP) = JSYMOP ! Total symmetry of operator
        IPRPTIM(NPRPS+IOP) = ITIM
        DO IMAT = 1,NMAT
          IAIND = JM4(IMAT,IOPTYP)
          IAREP = JM4REP(IAIND)
          IBSYM = IBTXOR(IAREP,IREP)+1
          IF(IPSYM2INT(IBSYM,IOP).EQ.-1) THEN
            PBUF(IMAT,IOP) = D0
          ELSE
            PBUF(IMAT,IOP)          = UPOL(IMAT) ! factor 0.5 not included
            IPRPLBL(IMAT,NPRPS+IOP) = IPSYM2INT(IBSYM,IOP)+NPRPLBL
          ENDIF
        ENDDO
        IF(DSUM(NMAT,PBUF(1,IOP),1).EQ.D0) THEN
           IPRPTYP(NPRPS+IOP) = -1  ! Means that there is no contribution
        ELSE
           IPRPTYP(NPRPS+IOP) = IOPTYP
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
      CALL HEADER('Full light-matter interaction: Anisotropic case',-1)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      NPAR = 0
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-e '
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-p) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-p '
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
      TPAR(NPAR+1) = 'Total'
C      
C     Loop over excitation energies
C     =============================
C
      allocate(TMOM(NOP,NPAR+1,NEXCNV))
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      DO 10 IEXC = 1,NEXCNV
        IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
           WRITE(LUPRI,*)
     &     'GET_BED: Skipping excitation. Frequency too small: ',
     &     OMEGA(IEXC)
           GOTO 10
        ENDIF
C
C       Length of wave vector
C 
        FAC=OMEGA(IEXC)/CVAL
        DO IOP = 1,NOP
           WRITE(PRPNAM(NPRPS+IOP)(9:16),'(A4,I4)') '_exc',IEXC
           DO IMAT = 1,NMAT
              FACPRP(IMAT,NPRPS+IOP) = PBUF(IMAT,IOP)/FAC
           ENDDO
        ENDDO
C
C       Debug output
C
        IF(IPRXPP.GE.2) THEN
          DO IOP = 1,NOP
            WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP),FAC
            IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
              WRITE(6,*) ' - skipped'
            ELSE
              DO IMAT = 1,NMAT
                II = IPRPLBL(IMAT,NPRPS+IOP)
                WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C     
C         Generate integrals
C
        IF(BEDCHK) THEN
           CALL DCOPY(3,D0,0,WAVEVEC,1)
        ELSE
          DO I = 1,3
            WAVEVEC(I) = FAC*UWAVE(I)
          ENDDO
        ENDIF
C       
        CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &     NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &     TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &     0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
        IF(IPRXPP.GE.5) THEN
          DO ICOMP = 1,NPRPINT
            WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &      '. Symmetry: ',REP(INTREP(ICOMP))
            CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
          ENDDO
        ENDIF
C
        OSCBED(IEXC)        = D0
        TMOM(1,NPAR+1,IEXC) = D0
        TMOM(2,NPAR+1,IEXC) = D0
        IPAR = 0
C
C       Gradient - orbital (e-e) part
C       ====================================
        IF(NZXOPE.GT.0) THEN
          IPAR = IPAR + 1
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
              CALL GPGET(INDXPR,GPOE(1,IOP),
     &             get_orbital_rotation_indices_pp(),
     &             NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &             WORK,KFREE,LFREE,IPRXPP)
              TMOM(IOP,IPAR,IEXC)=D2*
     &             DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
              IF(IPRXPP.GE.3) THEN
                WRITE(6,'(A16,2X,A/)') PRPNAM(INDXPR),' (e-e) part'
                CALL PRINT_GP_XP(GPOE,XPOE,
     &             get_orbital_rotation_indices_pp(),NZXOPE,NZ)
              ENDIF
            ELSE
              TMOM(IOP,IPAR,IEXC)=D0
            ENDIF
            TMOM(IOP,NPAR+1,IEXC)=TMOM(IOP,NPAR+1,IEXC)
     &                           +TMOM(IOP,IPAR,IEXC)
          ENDDO
        ENDIF
C
C       Gradient - orbital (e-p) part
C       ====================================
        IF(NZXOPP.GT.0) THEN       
          IPAR = IPAR + 1
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
              CALL GPGET(INDXPR,GPOP(1,IOP),
     &             get_orbital_rotation_indices_pn(),
     &             NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &             WORK,KFREE,LFREE,IPRXPP)
              TMOM(IOP,IPAR,IEXC)=D2*
     &             DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
              IF(IPRXPP.GE.3) THEN
                WRITE(6,'(A16,2X,A/)') PRPNAM(INDXPR),' (e-p) part'
                CALL PRINT_GP_XP(GPOP,XPOP,
     &             get_orbital_rotation_indices_pn(),NZXOPP,NZ)
              ENDIF
            ELSE
              TMOM(IOP,IPAR,IEXC)=D0
            ENDIF
            TMOM(IOP,NPAR+1,IEXC)=TMOM(IOP,NPAR+1,IEXC)
     &                           +TMOM(IOP,IPAR,IEXC)
          ENDDO
        ENDIF
C
C       Oscillator strength
C
        TMP = D0
        DO IOP = 1,NOP
          TMP = TMP + TMOM(IOP,NPAR+1,IEXC)*TMOM(IOP,NPAR+1,IEXC)
        ENDDO
        OSCBED(IEXC) = D2*OMEGA(IEXC)*TMP
 10   CONTINUE
C
C     Print section
C
      WRITE(LUPRI,'(3X,A,3F18.10)')                                    
     &     '* Unit wave vector   :', (UWAVE(I),I=1,3),                  
     &     '* Polarization vector:', (UPOL(I),I=1,3)                   
      WRITE(6,'(/A)') '* Transition moments:'
      WRITE(6,'(/A4,1X,A18,24X,3(3X,A5,13X))')     
     &  'Exc.','Frequency (au)    ',(TPAR(IPAR),IPAR=1,NPAR+1)
       DO IEXC = 1,NEXCNV
         WRITE(6,'(I4,1X,F18.10,A23,3(3X,ES18.10))')
     &   IEXC,OMEGA(IEXC),' real part (sin[kr]) : ',
     &   (TMOM(2,IPAR,IEXC),IPAR = 1,NPAR+1)
         WRITE(6,'(23X,A23,3(3X,ES18.10))')
     &   ' imag part (cos[kr]) : ',
     &   (-TMOM(1,IPAR,IEXC),IPAR = 1,NPAR+1)
      ENDDO
      WRITE(6,'(/A)') '* Oscillator strengths:'
      WRITE(6,'(/A4,1X,A18)')     
     &  'Exc.','Frequency (au)    '
      DO IEXC = 1,NEXCNV
        WRITE(6,'(I4,1X,F18.10,(3X,ES14.6))')
     &   IEXC,OMEGA(IEXC),OSCBED(IEXC)
      ENDDO
C     
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(TMOM)
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(IPSYM2INT)
      deallocate(first)
C     
      RETURN
      END
C
C  /* Deck geted_length */                                                  
      SUBROUTINE GET_ED_LENGTH(OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Test routine for electric dipole approximation in length gauge.
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC, to be opened with unit number LUXVC
C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D2=2.0D0,D1=1.0D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "pgroup.h"
#include "dgroup.h"
#include "nuclei.h"      
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      integer, allocatable              :: IVECS(:)
      real(8)                           :: EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      logical, allocatable              :: first(:),doatom(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(3))
      allocate(INTREP(3))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      allocate(doatom(nucind))
      CALL PR1DIR('DIPLEN  ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,DOATOM,NUCIND,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '+00+'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.TRUE.
      DOINT(2,1)=.FALSE.
      DOINT(1,2)=.FALSE.
      DOINT(2,2)=.TRUE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C
C     Operator setup
C
      NOP = 3
      IREP   = JSYMOP-1
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      IOPTYP  = 1                ! iA.P           * dot-product of alpha and vector operator
      NMAT    = MCMP(IOPTYP)     ! number of components, should be three
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &        ' ** ERROR in GET_ED_LENGTH **  Redimension MAXPRPS to ',
     &        NPRPS+NOP
        CALL QUIT('GET_ED_LENGTH: Too may properties !')
      ENDIF
      DO IOP = 1,NOP
        WRITE(PRPNAM(NPRPS+IOP),'(A4,A3,A1,A8)')
     &        'ED__',REP(IREP),'_',LABINT(IOP)
        IPRPSYM(NPRPS+IOP) = JSYMOP ! Total symmetry of operator
        IPRPTIM(NPRPS+IOP) = JTIMOP ! Time antisymmetric operator
        IF(INTREP(IOP).EQ.IREP) THEN
           IPRPTYP(NPRPS+IOP)   = IOPTYP
           IPRPLBL(1,NPRPS+IOP) = IOP + NPRPLBL
           FACPRP(1,NPRPS+IOP)  = D1
        ELSE
           IPRPTYP(NPRPS+IOP) = -1
           FACPRP(1,NPRPS+IOP) = D0
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      ITIM = JTIMOP
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))         
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)

         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-e) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C      
C     Loop over excitation energies
C     =============================
C     
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      WRITE(6,'(/A18,10X,A9,8X,A5)')
     &   'Transition moments','Freq(au)','Value'
      DO 10 IEXC = 1,KEXCNV(JSYMOP)
        IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
           WRITE(LUPRI,*)
     &     'GET_ED_LENGTH: Skipping excitation. Frequency too small: ',
     &     OMEGA(IEXC)
           GOTO 10
        ENDIF
C
C       Debug output
C
        IF(IPRXPP.GE.2) THEN
          DO IOP = 1,NOP
            WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP)
            IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
              WRITE(6,*) ' - skipped'
            ELSE
              DO IMAT = 1,1
                II = IPRPLBL(IMAT,NPRPS+IOP)
                WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C     
C       Generate integrals
C
        CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DOATOM,
     &       NUCIND,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
        IF(IPRXPP.GE.5) THEN
          DO ICOMP = 1,NPRPINT
            WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &        '. Symmetry: ',REP(INTREP(NCOMP))
            CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-e) part
C     ====================================
        IF(NZXOPE.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-e :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-p) part
C     ====================================
        IF(NZXOPP.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-p :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
 10   CONTINUE
C
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(first)
      deallocate(doatom)
C     
      RETURN
      END
C  /* Deck geted_length */                                                  
      SUBROUTINE GET_ED_VELOCITY(OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Test routine for electric dipole approximation in velocity gauge.
C     Scalar integrals are sorted on symmetries (max. 8)
C
C     Upon entering this routine solution vectors have been written
C     to unformatted file PAMXVC, to be opened with unit number LUXVC
C     References are:
C     1. Nanna Holmgaard List, Joanna Kauczor, Trond Saue,
C        Hans Jørgen Aagaard Jensen and Patrick Norman, J. Chem. Phys. 142 (2015) 244111
C                               
C                                                                      
C***********************************************************************
      use orbital_rotation_indices
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D2=2.0D0,D1=1.0D0)
#include "dcbgen.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "chrxyz.h"
#include "pgroup.h"
#include "dgroup.h"
#include "nuclei.h"      
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
C#include "cbiher.h"
      LOGICAL DOINT(2,2)
      INTEGER LWORK, NCOMP, NOPTP,IBEIG(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),CMO(*),WORK(*)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:)
      integer, allocatable              :: IVECS(:)
      real(8)                           :: EVECR(NREDM,*)
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      logical, allocatable              :: first(:),doatom(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE
C
#include "ibtfun.h"
      NOP = 3
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
C     
C     Allocate placeholders for integral labels and irreps 
C
      allocate(LABINT(1))
      allocate(INTREP(1))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      allocate(doatom(nucind))
      CALL PR1DIR('OVERLAP ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,DOATOM,NUCIND,NBAST,NELMNT,IPRINT)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_ED_VELOCITY:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C
C     Operator setup
C
      NOP    = 3
      IREP   = JSYMOP-1
      ITIM   = -JTIMOP ! Time antisymmetric operator
      IF((NPRPS+NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,A,I5)') 
     &        ' ** ERROR in GET_ED_VELOCITY **  ',
     &        'Redimension MAXPRPS to ',
     &        NPRPS+NOP
        CALL QUIT('GET_ED_LENGTH: Too may properties !')
      ENDIF
      DO IOP = 1,NOP
        WRITE(PRPNAM(NPRPS+IOP),'(A4,A3,2A1,A7)')
     &        'ED__',REP(IREP),'_',CHRXYZ(IOP),'dipvel_'
        IOPTYP             = 1 + IOP
        IAIND = JM4(1,IOPTYP)
        IAREP = JM4REP(IAIND)
        IPRPSYM(NPRPS+IOP) = JSYMOP ! Total symmetry of operator
        IPRPTIM(NPRPS+IOP) = ITIM   
        IF(IAREP.EQ.IREP) THEN
           IPRPTYP(NPRPS+IOP)   = IOPTYP
           IPRPLBL(1,NPRPS+IOP) = 1 + NPRPLBL
           FACPRP(1,NPRPS+IOP)  = CVAL
        ELSE
           IPRPTYP(NPRPS+IOP)  = -1
           FACPRP(1,NPRPS+IOP) = D0
        ENDIF
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))         
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-e) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C      
C     Loop over excitation energies
C     =============================
C     
      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      WRITE(6,'(/A18,10X,A9,8X,A5)')
     &   'Transition moments','Freq(au)','Value'
      DO 10 IEXC = 1,KEXCNV(JSYMOP)
        IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
           WRITE(LUPRI,*)
     &     'GET_ED_LENGTH: Skipping excitation. Frequency too small: ',
     &     OMEGA(IEXC)
           GOTO 10
        ENDIF
C
C       Debug output
C
        IF(IPRXPP.GE.2) THEN
          DO IOP = 1,NOP
            WRITE(6,*) PRPNAM(NPRPS+IOP),IPRPTYP(NPRPS+IOP)
            IF(IPRPTYP(NPRPS+IOP).EQ.-1) THEN
              WRITE(6,*) ' - skipped'
            ELSE
              DO IMAT = 1,1
                II = IPRPLBL(IMAT,NPRPS+IOP)
                WRITE(6,*) PRPLBL(II),FACPRP(IMAT,NPRPS+IOP)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C     
C       Generate integrals
C
        CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DOATOM,
     &       NUCIND,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
        IF(IPRXPP.GE.5) THEN
          DO ICOMP = 1,NPRPINT
            WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &        '. Symmetry: ',REP(INTREP(ICOMP))
            CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-e) part
C     ====================================
        IF(NZXOPE.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-e :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
C
C     Gradient - orbital (e-e) part
C     ====================================
        IF(NZXOPP.GT.0) THEN       
          DO IOP = 1,NOP
            INDXPR = NPRPS+IOP
            IF(IPRPTYP(INDXPR).NE.-1) THEN
               CALL GPGET(INDXPR,GPOP(1,IOP),
     &              get_orbital_rotation_indices_pn(),
     &              NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)
               TMP=D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
               WRITE(6,'(A,1X,A16,F18.10,3X,F14.10)')
     &            'e-p :',PRPNAM(NPRPS+IOP),OMEGA(IEXC),TMP
            ENDIF
          ENDDO
        ENDIF
 10   CONTINUE
C
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(first)
      deallocate(doatom)
C     
      RETURN
      END
C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module shared_bedled_var_module   
      LOGICAL NEWTASK

      INTEGER LWORK, NCOMP, NOPTP,IIP
      character (len=5)                 :: TPAR(3)
      real(8)             :: BCPU, PBUF(3,2),EEWT(6), TMP
      logical, allocatable              :: first(:)
      LOGICAL ANTI, TRASPH, SOLVNT, TRIANG, SQUARE         
      integer  POINT_PER_PROC, MOD_POINT, ITASK, SNODE, PPPP, NPOINTS
      integer  PSTART, PEND, PDIF, LINDX, MMAXEXC, MPIGO
      LOGICAL DOINT(2,2)

      integer, allocatable              :: POINTS(:)
      real(8), allocatable              :: PRPINT(:,:),BBUF(:)
      real(8), allocatable              :: TMOM(:,:,:),OSC(:)
      character (len=8),  allocatable   :: LABINT(:)
      integer, allocatable              :: INTREP(:),IPSYM2INT(:,:)
      integer, allocatable              :: IVECS(:)
      real(8), allocatable              :: GPOE(:,:),GPOP(:,:)
      real(8), allocatable              :: XPOE(:,:),XPOP(:,:)
      real(8), allocatable              :: X(:),Y(:),Z(:),WT(:)      
      integer, allocatable :: norot(:)
      LOGICAL GASWF
      end module shared_bedled_var_module
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bedleb */                                                   C
      SUBROUTINE BEDLEB(OSCBED,OMEGA,ISSYM,CMO,IBEIG,EVECR,IBTYP,
     &                   IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***********************************************************************
      use orbital_rotation_indices
      use interface_to_mpi
      use shared_bedled_var_module
C
C      use integer_kind_mpilib
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,DP5=0.5D0)
#include "maxash.h"
#include "dcborb.h"
C
#include "dcbgen.h"
#include "chrxyz.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
#include "maxorb.h"
#include "dcbidx.h"
#if defined (VAR_MPI)
#include "mpif.h"
#include "infpar.h"
#endif
C#include "cbiher.h"
C     ************************************************************
C     * Most variables have been moved to shared_bedled_var_module
C     * for sharing with MPI slave processes
C     ************************************************************
      CHARACTER*6 MXFORM, FMT
      CHARACTER SECTID*12, TOTTID*12, WALLTID*12, CPUTID*12 
      CHARACTER ALLCPUTID*12
#if defined (VAR_MPI)
      INTEGER :: ISTAT(df_mpi_status_size)
#endif
      INTEGER IBEIG(*)
      DIMENSION CMO(*),WORK(*)
      DIMENSION OMEGA(MAXEXC),ISSYM(MAXEXC),OSCBED(MAXEXC)
      real(8)             :: EVECR(NREDM,*), CPUTIME
      integer             :: IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*)
      integer   KFRSAV
C
#include "ibtfun.h"

      CALL QENTER('BEDLEB')
C
#if defined (VAR_MPI)
      MPARID = 0
      MASTER = MPARID
      call MPI_COMM_RANK(global_communicator, MYPROC, ierr)
      call MPI_COMM_SIZE(global_communicator, NMPROC, ierr)

C     ****************************************************
C     * All initialization will be performed by the 
C     * master node from this point 
C     ****************************************************     
C     
C     Setup timing of the BEDLEB routine 
C
      CPUTIME = 0
      CALL GETTIM(BCPU,WALL1)
      IF (MYPROC.eq.MASTER) THEN

#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
      NEXCNV = KEXCNV(JSYMOP)
      IF (NEXCNV .EQ. 0) RETURN
      NPOINTS = NLEB(NROTAV)
C     
C     Allocate placeholders for integral labels and irreps 
C
      NC = 2
      allocate(LABINT(8))
      allocate(INTREP(8))
      allocate(first(nz))
C
C     Integral directives
C
      IPRINT=0
      CALL PR1DIR('CXIKR   ',INTTYP,NPRPINT,INTREP,
     &            ANTI,SQUARE,DUMMY,LABINT,TRIANG,TRASPH,
     &            SOLVNT,0,.FALSE.,0,NBAST,NELMNT,IPRINT)
C     NPRPINT is by construction an even number,
C     the first and second halfs corresponding to the real (eq. 44) and imaginary parts (eq. 45)
      IF(NPRPLBL+NPRPINT.GT.MAXPRPLBL) THEN
        CALL QUIT('GET_BED:Redimension MAXPRPLBL !')
      ENDIF
      DO ICOMP = 1,NPRPINT
        PRPLBL (NPRPLBL+ICOMP) = LABINT(ICOMP)
        IPRLREP(NPRPLBL+ICOMP) = INTREP(ICOMP)
        PDOINT (NPRPLBL+ICOMP) = '0++0'
        IPRLTYP(NPRPLBL+ICOMP) =  1        
      ENDDO
      DOINT(1,1)=.FALSE.
      DOINT(2,1)=.TRUE.
      DOINT(1,2)=.TRUE.
      DOINT(2,2)=.FALSE.
C.....allocate array for holding integrals (on lower triangular form)
      allocate(PRPINT(NELMNT,NPRPINT))
      IORDER=0
      NPQUAD=0
C      
C     Find what symmetries are available
C
      allocate(IPSYM2INT(NBSYM,NC))
      NCOMP=NPRPINT/2
      CALL ICOPY(NC*NBSYM,-1,0,IPSYM2INT,1)
      DO I = 1,NCOMP
        ISYM = INTREP(I) + 1
        IPSYM2INT(ISYM,1) = I
        II = NCOMP + I
        ISYM = INTREP(II) + 1
        IPSYM2INT(ISYM,2) = II
      ENDDO
C
C     Operator setup
C
      NOP    = 3
      IREP   = JSYMOP-1
      ITIM   = -JTIMOP ! Time antisymmetric operator
      IF((NPRPS+2*NOP).GT.MAXPRPS) THEN
        WRITE(LUPRI,'(A,I5)') 
     &        ' ** ERROR in BEDLEB **  Redimension MAXPRPS to ',
     &        NPRPS+2*NOP
        CALL QUIT('BEDLEB: Too may properties !')
      ENDIF
      II = 0
      DO IC = 1,NC
        DO IOP = 1,NOP
          II = II + 1
          IOPTYP = 1 + IOP
          IAIND  = JM4(1,IOPTYP)
          IAREP  = JM4REP(IAIND)
          IBSYM  = IBTXOR(IAREP,IREP)+1
          WRITE(PRPNAM(NPRPS+II),'(A3,A1,I1,A7)')
     &        'BED',CHRXYZ(IC),IC,REP(IREP)
          IPRPSYM(NPRPS+II) = JSYMOP ! Total symmetry of operator
          IPRPTIM(NPRPS+II) = ITIM   
          IF(IPSYM2INT(IBSYM,IC).EQ.-1) THEN
            PBUF(IOP,IC) = D0
            IPRPTYP(NPRPS+II) = -1 ! Means that there is no contribution
          ELSE
            PBUF(IOP,IC)         = D1 ! factor 0.5 not included
            IPRPTYP(NPRPS+II)    = IOPTYP
            IPRPLBL(1,NPRPS+II)  = IPSYM2INT(IBSYM,IC)+NPRPLBL
          ENDIF
        ENDDO
      ENDDO 
      CALL PRSYMB(LUPRI,'=',80,0)
      CALL HEADER('Full light-matter interaction: Isotropic case',-1)
      WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',JSYMOP,' : ',REP(IREP)
      WRITE(6,*) 'BED rotational average. '
      WRITE(6,*) '  Lmax    = ',NROTAV
      WRITE(6,*) '  Npoints = ',NPOINTS
      WRITE(6,*) '  Nprocs  = ',NMPROC
C
C     Prepare solution vectors
C     ========================
C
C      - orbital (e-e) part      
C
      NPAR = 0
      IF(NZXOPE.GT.0) THEN
         allocate(XPOE(NZXOPEQ,NEXCNV))
         allocate(BBUF(NZXOPEQ))
         allocate(IVECS(NESIM))
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-e '
         OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBENDX,XPOE,EVECR,NEXCNV,IBTYP,
     &        IBEVC,IVECS,BBUF)
         CLOSE(LUBOE,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
C
C     - orbital (e-p) part      
C
      IF(NZXOPP.GT.0) THEN
         allocate(XPOP(NZXOPPQ,NEXCNV))
         allocate(BBUF(NZXOPPQ))
         allocate(IVECS(NPSIM))         
         NPAR = NPAR + 1
         TPAR(NPAR) = ' e-p '
         OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
         CALL XRSXV1(ITIM,JBPNDX,XPOP,EVECR,NEXCNV,IBTYP,
     &        IBPVC,IVECS,BBUF)
         CLOSE(LUBOP,STATUS = 'KEEP')
         deallocate(BBUF)
         deallocate(IVECS)         
      ENDIF
      TPAR(NPAR+1) = 'Total'
C     
C     Loop over quadrature points
C     =============================
C
      allocate(TMOM(NOP,NC,NEXCNV))
      allocate(OSC(NEXCNV))
      CALL DCOPY(NEXCNV,D0,0,OSC,1)

      IF(NZXOPE.GT.0) allocate(GPOE(NZXOPEQ,NOP))
      IF(NZXOPP.GT.0) allocate(GPOP(NZXOPPQ,NOP))      
      allocate( X(NPOINTS))
      allocate( Y(NPOINTS))
      allocate( Z(NPOINTS))
      allocate(WT(NPOINTS))
      CALL LEBEDEV(NROTAV,X,Y,Z,WT,NCOUNT)
      IF(NCOUNT.NE.NPOINTS) CALL QUIT('BEDLEB: Grid inconsistency.')

      allocate(POINTS(NPOINTS))
C    List of indicies
      POINT_PER_PROC = NPOINTS
      DO N = 1,NPOINTS
        POINTS(N) = N
      ENDDO

#ifdef VAR_MPI
C     ====================================================
C     MPI : Call slaves and broadcast data to Slaves 
C     ====================================================
      ENDIF
C     Initialize the Slaves (if there are any other than Master)
      KFRSAV = KFREE
      IF (NMPROC.gt.1) THEN
        CALL bedleb_master_distr(NOP,NC,NPAR,
     &       OMEGA,ISSYM,IBEIG,EVECR,IBTYP, IPRINT,
     &       IBCVC,IBEVC,IBPVC, INTTYP, NBAST,
     &       NELMNT,NPRPINT, CMO, WORK,KFREE,LFREE)
      ENDIF
#endif

C     ====================================================
C     LOOP OVER LEBDEV POINTS AND INTEGRATE 
C     ====================================================
      DO IIP = 1,POINT_PER_PROC
        IP = POINTS(IIP)

C.......Set up chi-average over polarization directions 
        EEWT(1) =  DP5*(Y(IP)*Y(IP)+Z(IP)*Z(IP)) ! (1,1)
        EEWT(2) = -X(IP)*Y(IP)                   ! (2,1) ! factor 0.5 cancelled by triangular sum
        EEWT(3) = -X(IP)*Z(IP)                   ! (3,1) ! factor 0.5 cancelled by triangular sum
        EEWT(4) =  DP5*(X(IP)*X(IP)+Z(IP)*Z(IP)) ! (2,2)
        EEWT(5) = -Y(IP)*Z(IP)                   ! (3,2) ! factor 0.5 cancelled by triangular sum
        EEWT(6) =  DP5*(X(IP)*X(IP)+Y(IP)*Y(IP)) ! (3,3)
        CALL DCOPY(NOP*NC*NEXCNV,D0,0,TMOM,1)
        DO 10 IEXC = 1,NEXCNV
          IF(ABS(OMEGA(IEXC)).LT.1.0D-3) THEN
             WRITE(LUPRI,*)
     &       'BEDLEB: Skipping excitation. Frequency too small: ',
     &       OMEGA(IEXC)
             GOTO 10
          ENDIF
C
C         Length of wave vector
C
          FAC=OMEGA(IEXC)/CVAL
          II = 0
          DO IC = 1,NC
            DO IOP = 1,NOP
              II = II + 1
              WRITE(PRPNAM(NPRPS+II)(9:16),'(A4,I4)') '_exc',IEXC
              FACPRP(1,NPRPS+II) = PBUF(IOP,IC)/FAC
            ENDDO
          ENDDO 
C     
C         Generate integrals
C
          IF(BEDCHK) THEN
            CALL DCOPY(3,D0,0,WAVEVEC,1)
          ELSE
            WAVEVEC(1) = FAC*X(IP)
            WAVEVEC(2) = FAC*Y(IP)
            WAVEVEC(3) = FAC*Z(IP)
          ENDIF
          CALL PR1IN2(WORK(KFREE),LFREE,PRPINT,INTTYP,
     &       NPRPINT,INTREP,ANTI,SQUARE,DUMMY,LABINT,TRIANG,
     &       TRASPH,.FALSE.,SOLVNT,0,0,DUMMY,
     &       0,NBAST,NELMNT,IPRINT,DOINT,.FALSE.)
          IF(IPRXPP.GE.5) THEN
            DO ICOMP = 1,NPRPINT
              WRITE(6,*) '* Integrals of ',LABINT(ICOMP),
     &          '. Symmetry: ',REP(INTREP(ICOMP))
              CALL OUTPAK(PRPINT(1,ICOMP),NBAST,1,LUPRI)
            ENDDO
          ENDIF
C
C         Accumulate transition moments
C         ============================= 
C
C         
C         - orbital (e-e) part
C         ====================================
          IF(NZXOPE.GT.0) THEN
            II = 0
            DO IC = 1,NC
              DO IOP = 1,NOP
                II = II + 1
                INDXPR = NPRPS+II
                IF(IPRPTYP(INDXPR).NE.-1) THEN
            
                  CALL GPGET(INDXPR,GPOE(1,IOP),
     &              get_orbital_rotation_indices_pp(),
     &              NZXOPE,CMO,IBEIG,JBENDX,.FALSE.,PRPINT,
     &              WORK,KFREE,LFREE,IPRXPP)

                  TMOM(IOP,IC,IEXC)=TMOM(IOP,IC,IEXC)
     &                + D2*DDOT(NZXOPEQ,GPOE(1,IOP),1,XPOE(1,IEXC),1)
                ENDIF
              ENDDO
            ENDDO
          ENDIF ! IF(NZXOPE.GT.0) THEN
C
C          - orbital (e-p) part
C         ====================================
          IF(NZXOPP.GT.0) THEN
            II = 0
            DO IC = 1,NC
              DO IOP = 1,NOP
                II = II + 1
                INDXPR = NPRPS+II
                IF(IPRPTYP(INDXPR).NE.-1) THEN
                  CALL GPGET(INDXPR,GPOP(1,IOP),
     &               get_orbital_rotation_indices_pn(),
     &               NZXOPP,CMO,IBEIG,JBPNDX,.FALSE.,PRPINT,
     &               WORK,KFREE,LFREE,IPRXPP)
                  TMOM(IOP,IC,IEXC)=TMOM(IOP,IC,IEXC)
     &              +D2*DDOT(NZXOPPQ,GPOP(1,IOP),1,XPOP(1,IEXC),1)
                ENDIF
              ENDDO
            ENDDO
          ENDIF ! IF(NZXOPP.GT.0) THEN
C
C         Accumulate oscillator strengths
C         ===============================
C
          TMP = D0
          DO IC = 1,NC
             TMP = TMP
     &     + EEWT(1)*TMOM(1,IC,IEXC)*TMOM(1,IC,IEXC)
     &     + EEWT(2)*TMOM(2,IC,IEXC)*TMOM(1,IC,IEXC)
     &     + EEWT(3)*TMOM(3,IC,IEXC)*TMOM(1,IC,IEXC)
     &     + EEWT(4)*TMOM(2,IC,IEXC)*TMOM(2,IC,IEXC)
     &     + EEWT(5)*TMOM(3,IC,IEXC)*TMOM(2,IC,IEXC)
     &     + EEWT(6)*TMOM(3,IC,IEXC)*TMOM(3,IC,IEXC)
          ENDDO
          OSC(IEXC) = OSC(IEXC) + WT(IP)*TMP
 10    CONTINUE                 ! DO 10 IEXC = 1,NEXCNV
      ENDDO ! DO IIP = 1,POINT_PER_PROC


#ifdef VAR_MPI
C     ********************************************************
C     If running parallel MPI - gather all DATA here
C     ********************************************************
      IF (NMPROC.gt.1) THEN
        call bedleb_gather_slave_data(NPAR, WORK,KFREE,LFREE,KFRSAV,
     &                                CPUTIME)
C       SLAVES WILL EXIT BEDLED HERE
        IF (MYPROC.ne.MASTER) THEN
          CALL QEXIT('BEDLEB')
          RETURN
        ENDIF
      ENDIF
C     
C     End timing of the BEDLEB routine 
C
      CALL GETTIM(CPU2,WALL2)
      CPU    = CPU2 - BCPU
      WALL   = WALL2 - WALL1
      TOTWAL = TOTWAL + WALL
      CPUTIME = CPUTIME + CPU
C
      TOTTID = SECTID(TOTWAL)
      WRITE(LUPRI,'(A,A12)')
     &    '>>>> Total wall time used in BEDLEB           : ',
     &       TOTTID
      CPUTID = SECTID(CPU)
      ALLCPUTID = SECTID(CPUTIME)
      WRITE(LUPRI,'(A,A12,A1,A12)')
     &    '>>>> Master CPU/Total CPU time used in BEDLEB : ',
     &       CPUTID,'/',ALLCPUTID
C     ********************************************************
#endif
C
C     Final values for oscillator strengths
C
      DO IEXC = 1,NEXCNV
        OSCBED(IEXC) = D2*OMEGA(IEXC)*OSC(IEXC)
      ENDDO
C
C     Print section
C
      WRITE(6,'(/A)') '* Oscillator strengths:'
      WRITE(6,'(/A4,1X,A18)')     
     &  'Exc.','Frequency (au)    '
      DO IEXC = 1,NEXCNV
        WRITE(6,'(I4,1X,F18.10,3X,ES14.6)')
     &   IEXC,OMEGA(IEXC),OSCBED(IEXC)
      ENDDO
C            
C     Memory deallocation
      IF(NZXOPE.GT.0) THEN
        deallocate(GPOE)
        deallocate(XPOE)
      ENDIF
      IF(NZXOPP.GT.0) THEN
        deallocate(GPOP)
        deallocate(XPOP)         
      ENDIF
      deallocate(TMOM)
      deallocate(OSC)
      deallocate(PRPINT)
      deallocate(LABINT)
      deallocate(INTREP)
      deallocate(IPSYM2INT)
      deallocate(first)
      deallocate(POINTS)
      deallocate(X)
      deallocate(Y)
      deallocate(Z)
      deallocate(WT)
      CALL QEXIT('BEDLEB')
C     
      RETURN
      END
#if defined (VAR_MPI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bedleb */                                                   C
C     * This is the slave interface to the BEDLEB subroutine below     C
C     * Will receive the necessary information to begin working        C
C     * on the BEDLEB routine - allocating the needed arrays           C
C     * for entering BEDLED                                            C
      SUBROUTINE bedleb_launch_slave_procs()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use memory_allocator                                 
      use interface_to_mpi                                 
#include "implicit.h"                                      
C                                                          
C                                                          
c TODO added due to NESH
#include "dcborb.h"
#include "dcbbas.h"
#include "priunit.h"                                       
#include "mxcent.h"                                        
#include "maxorb.h"                                        

#include "aovec.h"                                         
C TODO commented out due to dcborb.h                                                          
c#include "inforb.h"                                        
#include "infpar.h"                                        
#include "dcbxpp.h"
#include "molinp.h"                                        
#include "blocks.h"                                        
#if defined (PRG_DIRAC)                                    
#include "dcbgen.h"                                        
#include "dcbgrd.h"                                        
#else                                                      
#include "gnrinf.h"                                        
#include "energy.h"                                        
#endif                                                     
C
C                                                          
      DIMENSION DINTSKP(14*3*3)                            
      INTEGER LWORK, NCOMP, NOPTP, MYPROC
      LOGICAL   NODV,NOPV,NOCONT,RETUR,TKTIME,NEWGEO,FINISH
C                                                          
      real(8), allocatable :: WORK(:), CMO(:), OMEGA(:), TMPISSYM(:)
      real(8), allocatable :: TMPOSCBED(:)
      integer, allocatable :: IBEIG(:)
      INTEGER MAXEXC, I
      call legacy_lwork_get(LWORK)                         
      call alloc(WORK,LWORK,id='WORK in BEDLEB')           
#include "memint.h"

      call MPI_COMM_RANK(global_communicator, MYPROC, ierr)
C                                                          
      CAll interface_mpi_BCAST(MAXEXC,1,MPARID,
     &                  global_communicator)
      allocate(TMPOSCBED(MAXEXC))
      CAll interface_mpi_BCAST(NREDM,1,MPARID,
     &                  global_communicator)

C     * OMEGA
      allocate(OMEGA(MAXEXC))
      CAll interface_mpi_bcast_r1_work_f77(OMEGA,MAXEXC,MPARID,
     &                  global_communicator)
C     * ISSYM

      allocate(TMPISSYM(MAXEXC))
      CAll interface_mpi_bcast_r1_work_f77(TMPISSYM,MAXEXC,MPARID,
     &                  global_communicator)
C     * N2BBASXQ for Allocation of CMO
      call interface_mpi_bcast(N2BBASXQ, 1, MPARID,
     &              global_communicator)
C     * IBEIG
      call interface_mpi_bcast(NTBAS, 3, MPARID,
     &              global_communicator)

      allocate(CMO(N2BBASXQ))
      allocate(IBEIG(NTBAS(0)))
      CALL BEDLEB(TMPOSCBED,OMEGA,TMPISSYM,CMO,IBEIG,EVECR,IBTYP,
     &            IBCVC,IBEVC,IBPVC,WORK,KFREE,LFREE)
      deallocate(TMPISSYM)
      deallocate(OMEGA)
      deallocate(CMO)
      deallocate(IBEIG)
      deallocate(TMPOSCBED)
      RETURN

      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bedleb */                                                   C
C     * This master routine initially starts/calls all the slaves      C
C     * by sending the iTask ID to main and then follows up            C
C     * by broadcasting all the initial allocatable variables          C
C                                                                      C
      SUBROUTINE bedleb_gather_slave_data(NPAR,WORK,KFREE,LFREE,KFRSAV,
     &                                    CPUTIME)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use orbital_rotation_indices
      use interface_to_mpi
      use shared_bedled_var_module
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,DP5=0.5D0)
      real(8)             :: TMPCPU
#include "maxash.h"
#include "dcborb.h"
C
#include "dcbgen.h"
#include "chrxyz.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
#include "maxorb.h"
#include "dcbidx.h"
#include "mpif.h"
      integer MPARID, MASTER, NMPROC, MYPROC
C
#include "ibtfun.h"

      call MPI_COMM_SIZE(global_communicator, NMPROC, ierr)
      call MPI_COMM_RANK(global_communicator, MYPROC, ierr)

      CALL GETTIM(CPU2,WALL2)
      CPU    = CPU2 - BCPU

      MPARID = 0
      MASTER = MPARID
C     Wait for all processes before doing any Send/Recv
      call interface_mpi_barrier(global_communicator)

C     If MPI - this remaining part will only be  done by master
C     *********************************************************     
C     * Send All Slave results back to the Master
C     *********************************************************     
      IF (MYPROC.ne.MASTER) THEN
C       Send back the OSC Results
        CALL interface_MPI_RECV(MPIGO, 1,
     &                          MASTER, 27,global_communicator)
C       *************** Sends OSC ***************
        DO IEXC = 1,NEXCNV
          tmp = OSC(IEXC)
          CALL interface_MPI_SEND(tmp, 1,
     &                            MASTER, 28,global_communicator)
        ENDDO
        CALL interface_MPI_SEND(CPU, 1,
     &                          MASTER, 28,global_communicator)
C       *************** Deallocate ***************

C       Deallocate all data for the Slaves
        deallocate(TMOM)
        deallocate(OSC)
        deallocate(X)
        deallocate(Y)
        deallocate(Z)
        deallocate(WT)
        deallocate(POINTS)
        deallocate(PRPINT)
        deallocate(INTREP)
        deallocate(LABINT)
        IF(NZXOPE.GT.0) THEN
          deallocate(GPOE)
          deallocate(XPOE)
        ENDIF
        IF(NZXOPP.GT.0) THEN
          deallocate(GPOP)
          deallocate(XPOP)
        ENDIF
C     *********************************************************     
C     * Receive All Slave results
C     *********************************************************     
      ELSE
C       Receive results from each Node
C       And each IEXC
C       ****************************************************
C       * Same as splitting - now receiving OSC val from
C       * each slave node
C       ****************************************************     
C       Offset by POINT_PER_PROC - since master will also work
        PEND = POINT_PER_PROC 
C
        DO SNODE = 1, NMPROC-1
C       ************ Send go to Slave ************
        CALL interface_MPI_SEND(1, 1,
     &                          SNODE, 27,global_communicator)
C         ------------------------------
C         Calc Point distrib
C         ------------------------------
C         Offset by POINT_PER_PROC - since master will also work
          PSTART = 1 + PEND 
          PEND   = PSTART + POINT_PER_PROC - 1

C         Once the the remainder is equally divisible - make it so
          if ((SNODE).eq.MOD_POINT) then
              POINT_PER_PROC = POINT_PER_PROC - 1
          endif

C         Make sure that the 'end' is max the length if the list
          if (PEND.gt.NPOINTS) PEND = NPOINTS

C         ------------------------------
C         Receive OSC results
C         ------------------------------
          DO IEXC = 1,NEXCNV
            CALL interface_MPI_RECV(TMP,1,
     &                            SNODE, 28,global_communicator)
            OSC(IEXC) = OSC(IEXC) + TMP
          ENDDO
          CALL interface_MPI_RECV(TMPCPU, 1,
     &                          SNODE, 28,global_communicator)
          CPUTIME = CPUTIME + TMPCPU
C       end SNODE
        ENDDO
      ENDIF
c     ********************************************************** 
c     * END OF GATHERING SLAVE DATA
c     ********************************************************** 
      call interface_mpi_barrier(global_communicator)

      CALL MEMREL('BEDLEB',WORK,1,KFRSAV,KFREE,LFREE)
c     ********************************************************** 
c     * END OF MPI PART 
c     ********************************************************** 
      Return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck bedleb */                                                   C
C     * This master routine initially starts/calls all the slaves      C
C     * by sending the iTask ID to main and then follows up            C
C     * by broadcasting all the initial allocatable variables          C
C                                                                      C
      SUBROUTINE bedleb_master_distr(NOP,NC,NPAR,
     &       OMEGA,ISSYM,IBEIG,EVECR,IBTYP, IPRINT,
     &       IBCVC,IBEVC,IBPVC,INTTYP, NBAST,
     &       NELMNT,NPRPINT, CMO, WORK,KFREE,LFREE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use orbital_rotation_indices
      use interface_to_mpi
      use shared_bedled_var_module
C     
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,DP5=0.5D0)
#include "maxash.h"
#include "dcborb.h"
C
#include "dcbgen.h"
#include "chrxyz.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "orgcom.h"
#include "dcbxpr.h"
#include "dcbprl.h"            
#include "dcbibn.h"
#include "dcbwav.h"
#include "maxorb.h"
#include "dcbidx.h"
#include "mpif.h"
      integer MPARID, MASTER, NMPROC, MYPROC
C
#include "ibtfun.h"

      call MPI_COMM_SIZE(global_communicator, NMPROC, ierr)
      call MPI_COMM_RANK(global_communicator, MYPROC, ierr)

      MPARID = 0
      MASTER = MPARID
C     * =================================================
C     *       MASTER ONLY
C     * =================================================
      IF (MYPROC.eq.MASTER) THEN
C       ****************************************************
C       *     Scatter the points to calculate to the slaves
C       ****************************************************     
C       Reduce NPOINTS to the number of points per process
C       Calculate the number of points per process 
        POINT_PER_PROC = NPOINTS / NMPROC
        MOD_POINT = MOD(NPOINTS, NMPROC)
        IF (MOD_POINT.gt.0) POINT_PER_PROC = POINT_PER_PROC + 1

C       ****************************************************
C       * Call upon the slaves in main/dirac.f (DIRNOD) 
C       * and begin to distribute the data created    
C       ****************************************************     
        ITASK = 21  
        DO SNODE = 1, NMPROC-1
          CALL interface_MPI_SEND(ITASK, 1, SNODE, 27,
     &                            global_communicator)
        ENDDO
C       ****************************************************
C       * Distribute the data needed for the slaves to enter
C       * BEDLED in the first place : 
C       *     -> bedleb_launch_slave_procs() 
C       ****************************************************     
C       Sizes for OMEGA
        
        CAll interface_mpi_BCAST(MAXEXC,1,MPARID,
     &                    global_communicator)
        CAll interface_mpi_BCAST(NREDM,1,MPARID,
     &                    global_communicator)
C       OMEGA      
        CAll interface_mpi_bcast_r1_work_f77(OMEGA,MAXEXC,MPARID,
     &                    global_communicator)
C       ISSYM
        CAll interface_mpi_bcast_r1_work_f77(ISSYM,MAXEXC,MPARID,
     &                    global_communicator)
C       Sizes for CMO and IBEIG
        call interface_mpi_bcast(N2BBASXQ, 1, MPARID,
     &                global_communicator)
        call interface_mpi_bcast(NTBAS, 3, MPARID,
     &                global_communicator)

C       *****************************************************
C       * Split and distribute the points to each slave 
C       * process to calculate - each slave will receive
C       * this information in below in the SLAVE ONLY SECTION
C       *****************************************************     
        PEND = POINT_PER_PROC 
        PPPP = POINT_PER_PROC
        DO SNODE = 1, NMPROC-1
C         Offset by POINT_PER_PROC - since master will also work
          PSTART = 1 + PEND 
          PEND   = PSTART + POINT_PER_PROC - 1

C         Once the the remainder is equally divisible - make it so
          if ((SNODE).eq.MOD_POINT) POINT_PER_PROC = POINT_PER_PROC - 1

C         Make sure that the 'end' is max the length if the list
          if (PEND.gt.NPOINTS) PEND = NPOINTS

          PDIF = PEND - PSTART + 1
C         ----------------------------------------------------------- 
C         - Send a sub array of points to the target slave node
C         ----------------------------------------------------------- 
          CALL interface_MPI_SEND(PDIF, 1, SNODE, 27,
     &                            global_communicator)
          CALL interface_MPI_SEND(POINTS(PSTART:PEND), 
     &                            PDIF, SNODE, 27,global_communicator)
        ENDDO
C       If 1 was subtracted - put it back
        POINT_PER_PROC = PPPP
        POINTS = POINTS(1:POINT_PER_PROC)
      ELSE
C       * =================================================
C       *       SLAVES ONLY
C       * =================================================
        CALL interface_MPI_RECV(POINT_PER_PROC,1,MASTER,27,
     &                            global_communicator)
        allocate(POINTS(POINT_PER_PROC))
        CALL interface_MPI_RECV(POINTS,POINT_PER_PROC,MASTER,27,
     &                            global_communicator)
C       * =================================================
      ENDIF


CCCC     ****************************************************
CCCC     * Broadcast all the information needed to 
CCCC     * perform the BEDLEB integral for the Slaves 
CCCC     ****************************************************
      call interface_mpi_bcast(LFREE, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPOINTS, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NEXCNV, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPAR, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NELMNT, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NBAST, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NC, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NOP, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast_r2_work_f77(PBUF, 3,2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPRPINT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NZXOPE, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NZXOPEQ, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NZXOPPQ, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NZXOPP, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NZHOPE, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(norbt, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(SKIPEE, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(SKIPEP, 1,  MASTER,
     &              global_communicator)
    
      IF (MYPROC.ne.MASTER) THEN
        KFREE = 1
        CALL MEMGET2('REAL','LINDX',LINDX,NELMNT,
     &               WORK,KFREE,LFREE)
        
        allocate(TMOM(NOP,NC,NEXCNV))
        allocate(OSC(NEXCNV))
        allocate(X(NPOINTS))
        allocate(Y(NPOINTS))
        allocate(Z(NPOINTS))
        allocate(WT(NPOINTS))
        allocate(PRPINT(NELMNT,NPRPINT))
        allocate(INTREP(8))
        allocate(LABINT(8))
        IF(NZXOPE.GT.0) THEN
          allocate(GPOE(NZXOPEQ,NOP))
          allocate(XPOE(NZXOPEQ,NEXCNV))
          call DZERO(GPOE, NZXOPEQ*NOP)
        ENDIF
        IF(NZXOPP.GT.0) THEN
          allocate(GPOP(NZXOPPQ,NOP))      
          allocate(XPOP(NZXOPPQ,NEXCNV))
          call DZERO(GPOP, NZXOPPQ*NOP)
        ENDIF
      ENDIF

      call DZERO(TMOM, NOP*NC*NEXCNV)
      call interface_mpi_bcast(OSC(1), NEXCNV, 
     &             MASTER, global_communicator)
      call interface_mpi_bcast(X(1), NPOINTS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(Y(1), NPOINTS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(Z(1), NPOINTS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(WT(1), NPOINTS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(N2BBASXQ, 1, MPARID,
     &              global_communicator)
      call interface_mpi_bcast_r1_work_f77(CMO,N2BBASXQ, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(ICMOQ, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(I2BASX(1,1), 2*2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(I2BAST, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(ICMO, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NCMOTQ, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NCMOT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NCMOQ, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NCMO, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(N2ORB, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NORBT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(N2ORBT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(N2ORBX, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NISH, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NFRO, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NASH, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NESH, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPSH, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(INTTYP, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPRPINT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(INTREP, 8*1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(ANTI, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(SQUARE, 1,  MASTER,
     &              global_communicator)
      call interface_mpi_bcast(LABINT(1), 8*8, MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(TRIANG, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(TRASPH, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l0(SOLVNT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast_l1(DOINT, 4*1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(DUMMY, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IPRINT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(JOPSY, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(N2BBASX, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NFBAS, 6, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NTBAS, 3, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NBBAS, 8*3, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IBBAS, 8*3, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(MXFBAS, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(CVAL, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPRPS, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IOBTYP, MXCORB, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IFSMO, MXCORB, MASTER,
     &              global_communicator)
C    ** bcast file name labels **
      call interface_mpi_bcast(IPRPLBL(1,1), 3*MAXPRPS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(PRPLBL(1), 8*MAXPRPLBL, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(PRPNAM(1), 16*MAXPRPS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IPRPTYP, MAXPRPS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IPRPSYM, MAXPRPS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IPRPTIM, MAXPRPS, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IPRLTYP, MAXPRPLBL, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IPRLREP, MAXPRPLBL, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(PDOINT, 4*MAXPRPLBL, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPMO, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NORB, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(IORB, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NFSYM, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPMOT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NNBBASX, 2, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(NPRPLBL, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(LU1INT, 1, MASTER,
     &              global_communicator)
      call interface_mpi_bcast(INDSTR(1,1),6*72, MASTER,
     &              global_communicator)
C     All slaves must setup orbital informations
      IF (MYPROC.ne.MASTER) THEN
        allocate(norot(norbt*1))
        call IZERO(norot, norbt*1)
        GASWF = .FALSE.
        CALL SETXOP(INDSTR,SKIPEE,SKIPEP,GASWF,JOPSY,NOROT,
     &            NZXOPE,NZHOPE,NZXOPP,
     &            IPRINT)
        deallocate(norot)
      ENDIF
C
      call interface_mpi_bcast(NZXOPP, 1, MASTER,
     &              global_communicator)
      IF(NZXOPP.GT.0) then
        call interface_mpi_bcast(XPOP(1,1), NZXOPPQ*NEXCNV, 
     &             MASTER, global_communicator)
      endif
      call interface_mpi_bcast(NZXOPE, 1, MASTER,
     &              global_communicator)
      if (NZXOPE.GT.0) then
        call interface_mpi_bcast(XPOE(1,1), NZXOPEQ*NEXCNV, 
     &             MASTER, global_communicator)
      endif

      RETURN
      END
#endif      
CCCCCCCCCCCCCCCCCCCCCCCCC
