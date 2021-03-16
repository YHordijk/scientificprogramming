!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE property_handler(ioch,ioch2,lavail,nsp,irepspi)

      use polprp_cfg   !make some input variables available

      IMPLICIT none
!
! this routine reads in the MO-transformed
! property integrals from MDPROP via RDPROP_ADC maintaining
! complex arrays throughout.
! RDPROP_ADC is a modified version of the old relccsd/rdprop routine.
! The property integrals are analyzed for the CARITH=F cases where we still
! can have complex wave functions.
!
!---------------Calling variables--------------------------------------
!
      integer                          ::  ioch,ioch2  !file handles
      logical                          ::  lavail      !flag for property integrals
      integer                          ::  nsp         !number of active spinors
      integer, dimension(:)            ::  irepspi     !spinor translation array
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!provides MXREP,NXREP etc.
#include  "../relccsd/complex.inc"
!provides CARITH,RCW
!
!---------------Local variables--------------------------------------
!
      logical      :: existing,success
      integer      :: npaoo, npaov, npavo, npavv  !specific prop array dimensions

      integer      :: i,iprop
      integer      :: MDPROP

      character*8  :: nameofprop
      real*8, allocatable   :: CBUF(:)
      real*8, allocatable   :: CPVO(:), CPOV(:), CPOO(:), CPVV(:)
      complex*16            :: cphase(4)
      character*8,dimension(3)   :: naoplist
!
!---------------Executable code --------------------------------------
!
      call PST('Entering property handler+')
!
!**  check for availability of MDPROP
!
      MDPROP = 15
      lavail = .false.
      INQUIRE(FILE='MDPROP',EXIST=existing)
      IF (.NOT.existing) THEN
         write(*,*) '  ** Warning ***'
         write(*,*) '  ** Property integrals (MDPROP) do not exist.'
         write(*,*) '  ** No transition moments are calculated.'
         write(*,*) '  ** Continuing with energy-only calculations.'
         RETURN
      ELSE
        lavail=.true.
        write(*,*) '  ** MDPROP file found.'
      END IF
!     write(*,*)
!     write(*,*) 'Print level is set to ',polprp_printlev
!     write(*,*)
!
! allocate proper space for the fou kinds of property integrals.
!
      npaoo = 0
      npaov = 0
      npavo = 0
      npavv = 0
      do i=1,nrep
        npaoo = npaoo + MOO(i)
        npaov = npaov + MOV(i)
        npavo = npavo + MVO(i)
        npavv = npavv + MVV(i)
      enddo
!
! check dimensions for consistency
!
      if( (npaoo+npaov+npavo+npavv) .ne. (nsp*nsp)) then
        write(*,*) '  ** Consistency check failed with:'
        write(*,*) '  ** property array has length',nsp*nsp
        write(*,*) '  ** Calculated dims are:',
     &                npaoo+npaov+npavo+npavv
        STOP
!     else
!       write(*,*) '  ** Consistency check passed with:',nsp*nsp
!       write(*,*) '  ** Size of OO array:',npaoo
!       write(*,*) '  ** Size of OV array:',npaov
!       write(*,*) '  ** Size of VO array:',npavo
!       write(*,*) '  ** Size of VV array:',npavv
      endif
!
! allocate proper space for property integrals
!
      allocate(CBUF(2*NSP*NSP))   !Always complex, holding raw integrals
      allocate(CPOO(RCW*npaoo))
      allocate(CPOV(RCW*npaov))
      allocate(CPVO(RCW*npavo))
      allocate(CPVV(RCW*npavv))
!
! set list of property names required for the TMs
!
      naoplist(1) = 'XDIPLEN'
      naoplist(2) = 'YDIPLEN'
      naoplist(3) = 'ZDIPLEN'
!
! loop over the three properties and store the integrals in the
! corresponding files.
!

      Do iprop = 1,3

        nameofprop=naoplist(iprop)

        success = .false.

        CALL read_flush_prop(success,nameofprop,NSP,IREPSPI,
     &                       npaoo,npaov,npavo,npavv,
     &                       CBUF,CPOO,CPOV,CPVO,CPVV,
     &                       polprp_printlev)

        if(.not.success) then
          call insuff_prop_return(nameofprop)
          lavail = .false.
          deallocate(CPVV)
          deallocate(CPVO)
          deallocate(CPOV)
          deallocate(CPOO)
          deallocate(CBUF)
          return
        endif

      Enddo
!
! releasing heap memory
!
      deallocate(CPVV)
      deallocate(CPVO)
      deallocate(CPOV)
      deallocate(CPOO)
      deallocate(CBUF)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE insuff_prop_return(nameofprop)
      character*8  :: nameofprop
      write(*,*)
      write(*,*) '  *** WARNING ***'
      write(*,*)
      write(*,*) 'Property file MDPROP available, but required dipole'
      write(*,*) 'moment properties are not available there.'
      write(*,*)
      write(*,*) 'Please add the following lines to your'
      write(*,*) '**MOLTRA section, in order to create the required'
      write(*,*) 'operator names for property processing:'
      write(*,*)
      write(*,*) '**MOLTRA'
      write(*,*) '.ACTIVE'
      write(*,*) 'energy ...  ...   ...'
      write(*,*) '.PRPTRA'
      write(*,*) '*PRPTRA'
      write(*,*) '.OPERATOR'
      write(*,*) ' XDIPLEN'
      write(*,*) '.OPERATOR'
      write(*,*) ' YDIPLEN'
      write(*,*) '.OPERATOR'
      write(*,*) ' ZDIPLEN'
      write(*,*)
      write(*,*) 'Calculation continues without transition moments.'
      write(*,*)
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PROPFETCH(ioch,pph,
     &                     idoo,idov,idvo,idvv,
     &                     doox,dooy,dooz,
     &                     dovx,dovy,dovz,
     &                     dvox,dvoy,dvoz,
     &                     dvvx,dvvy,dvvz)
     &               

 
      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!     Reads all property matrices from disk
!
!---------------Calling variables--------------------------------------
!
      Integer                    :: ioch
      Complex*16,dimension(:)    :: pph
      Integer                    :: idoo,idov,idvo,idvv
      Real*8,dimension(:)        :: doox,dooy,dooz
      Real*8,dimension(:)        :: dovx,dovy,dovz
      Real*8,dimension(:)        :: dvox,dvoy,dvoz
      Real*8,dimension(:)        :: dvvx,dvvy,dvvz
!
!---------------Common blocks----------------------------------------
!
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      Character*1              :: co
      Character*8              :: filename
      Character*7              :: nameroot
      Integer                  :: icaoo,icaov,icavo,icavv,ixx,ipix,i
      Integer                  :: ist
      Complex*16               :: A0,A1,AI
      Real*8                   :: ur,ui
!
!---------------Executable code--------------------------------------
!
      A0=(0.0d0,0.0d0)
      A1=(1.0d0,0.0d0)
      AI=(0.0d0,1.0d0)
      nameroot = 'PROPINT'

!  ********** X

      co='X'
      ipix=1
      write(filename,'(A1,A7)') co,nameroot
      open(unit=ioch,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) icaoo,icaov,icavo,icavv
      if( (  (icaoo - idoo) +
     &       (icaov - idov) +
     &       (icavo - idvo) +
     &       (icavv - idvv) ).ne.0) then
        write(*,*) 'Dimension inconsistency in ',filename
        call quit('***')
      endif
      read(ioch) ur,ui
      read(ioch) (doox(ixx),ixx=1,rcw*idoo)
      read(ioch) (dovx(ixx),ixx=1,rcw*idov)
      read(ioch) (dvox(ixx),ixx=1,rcw*idvo)
      read(ioch) (dvvx(ixx),ixx=1,rcw*idvv)
      close(ioch)
      pph(ipix) = dcmplx(ur,ui)
      write(*,*) ' ** Phase for ',filename,'  ',pph(ipix)
      If(carith) then
          if(pph(ipix).ne.A0) then
            write(*,*) 'Phase inconsistency in ',filename
          endif
      endif
      write(*,*) ' ** Integrals from ',filename,' successfully read.'


!  ********** Y

      co='Y'
      ipix=2
      write(filename,'(A1,A7)') co,nameroot
      open(unit=ioch,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) icaoo,icaov,icavo,icavv
      if( (  (icaoo - idoo) +
     &       (icaov - idov) +
     &       (icavo - idvo) +
     &       (icavv - idvv) ).ne.0) then
        write(*,*) 'Dimension inconsistency in ',filename
        call quit('***')
      endif
      read(ioch) ur,ui
      read(ioch) (dooy(ixx),ixx=1,rcw*idoo)
      read(ioch) (dovy(ixx),ixx=1,rcw*idov)
      read(ioch) (dvoy(ixx),ixx=1,rcw*idvo)
      read(ioch) (dvvy(ixx),ixx=1,rcw*idvv)
      close(ioch)
      pph(ipix) = dcmplx(ur,ui)
      write(*,*) ' ** Phase for ',filename,'  ',pph(ipix)
      If(carith) then
          if(pph(ipix).ne.A0) then
            write(*,*) 'Phase inconsistency in ',filename
          endif
      endif
      write(*,*) ' ** Integrals from ',filename,' successfully read.'

!  ********** Z

      co='Z'
      ipix=3
      write(filename,'(A1,A7)') co,nameroot
      open(unit=ioch,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) icaoo,icaov,icavo,icavv
      if( (  (icaoo - idoo) +
     &       (icaov - idov) +
     &       (icavo - idvo) +
     &       (icavv - idvv) ).ne.0) then
        write(*,*) 'Dimension inconsistency in ',filename
        call quit('***')
      endif
      read(ioch) ur,ui
      read(ioch) (dooz(ixx),ixx=1,rcw*idoo)
      read(ioch) (dovz(ixx),ixx=1,rcw*idov)
      read(ioch) (dvoz(ixx),ixx=1,rcw*idvo)
      read(ioch) (dvvz(ixx),ixx=1,rcw*idvv)
      close(ioch)
      pph(ipix) = dcmplx(ur,ui)
      write(*,*) ' ** Phase for ',filename,'  ',pph(ipix)
      If(carith) then
          if(pph(ipix).ne.A0) then
            write(*,*) 'Phase inconsistency in ',filename
          endif
      endif
      write(*,*) ' ** Integrals from ',filename,' successfully read.'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE read_flush_prop(SUCC,NAMEA,NSP,IREPSPI,
     &                           dimpoo,dimpov,dimpvo,dimpvv,
     &                           PROP,AOO,AOV,AVO,AVV,
     &                           printlevel)
!
      use interface_to_mpi
      use qstack       !make stack system available

      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!
!     Reads property integral files for operator A.
!     The routine was originally devised by LV and modified/extended
!     by MP for the RELADC purposes.
!     This version reads from the file created by PRTRAN
!
!     Important background for the use of MO-transformed property integrals:
!
!     Due to the possibility of spinfree calculations the wave function can
!     still be complex leading to complex property integrals also in the case
!     of CARITH=false. CARITH does not indicate if a calculation is quaternionic
!     or not but tells arithmetic only for calculations after the MOLTRA step.
!     A Spinfree calculation therefore leads to complex MOPIs even if CARITH=
!     false! In this case we can get purely imaginary property integrals for the
!     x,y,a dipole moment. Since carith=false, the integrals end up in the AXX
!     arrays which causes AOV = -AVO **seemingly violating Hermiticity**
!
!     We therefore need a canary indicating if we have purely imaginary integrals
!     The real array is therefore to be interpreted as the imaginary part of the
!     property integral in the case CARITH=F. The APHASE array transports the corresponding
!     phase information (purely real or purely imaginary) to the calling routine.
!
!     This allows for sticking at real arrays in the case of CARITH=F and postmultiplying
!     the F_ai matrix elements by the imaginary unit. The F_ai's therefore have to be kept
!     complex also in the SPINFREE/NONREL cases.
!
!     ATT: the absolute spinor numbers obtained through the irepspi array are *not*
!     the same ones as listed in the list of active occ/virt orbitals at the beginning
!     of CCSETI output! The irepspi absolute indices serve as pointers into the 
!     PROP array where they correspond to the proper indexing.
!
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/param.inc"
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
#include "../relccsd/ccpar.inc"
#include "polprp_stacklines.h"
!
!---------------Calling variables--------------------------------------
!
      Logical                    SUCC
      Character*8                NAMEA
      INTEGER                    NSP
      INTEGER                    IREPSPI(NSP,MXREP,2)
      INTEGER                    dimpoo,dimpov,dimpvo,dimpvv
      REAL*8                     PROP(2,NSP,NSP)
      REAL*8                     AOO(*),AOV(*),AVO(*),AVV(*)
      integer                    printlevel
!
!---------------Local variables--------------------------------------
!
      Integer                  :: MDPROP, PROCPROP
      REAL*8,allocatable       :: DUMPROP(:,:,:)
      CHARACTER*32             :: ACHAR
      REAL*8                   :: APHASE(2)
      REAL*8                   :: AVOI,AOVI,AMAXR,AMAXI,ANORM
      LOGICAL                  :: PHASE
      REAL*8                   :: hdev
      COMPLEX*16               :: h1,h2
! --- flushing
      Character*1              :: co
      Character*8              :: filename
      Integer                  :: ixx,rawnum
      real*8, parameter        :: ZEROTHRS=1.0D-11
      integer,allocatable      :: acv(:,:)
      integer                  :: irow,icol,retcode
      integer                  :: loopcount
      Real*8, allocatable      :: GOV(:)    !generic order OV propints

      Integer                  :: qst_auxline
!
!---------------Executable code--------------------------------------
!
!     Read symmetry information (master node only) from MDPROP
!     The file is organized in a way that after a label identifying the
!     property a binary number stream without CRLF follows. If the label
!     is no match, we read til the next CRLF which corresponds to a skipping
!     of all the unwanted integrals.
!
      MDPROP = 15
      PROCPROP = MDPROP + 1

      IF(dimpov.ne.dimpvo) CALL QUIT('Inconsistent dimensions!')

      IF (MYPROC.EQ.MASTER) THEN

      allocate(DUMPROP(NSP,NSP,2))
      DUMPROP = 0.0d0
      OPEN (MDPROP,FILE='MDPROP',FORM='UNFORMATTED')
      INOP = 0
    1 READ (MDPROP,ERR=10,END=11) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.NAMEA) GOTO 1
      WRITE (6,1000) NAMEA
      READ (MDPROP,ERR=10,END=11) 
     &    ((DUMPROP(I,J,1),DUMPROP(I,J,2),I=1,NSP),J=1,NSP)
      PROP=0.0d0
      Do I=1,NSP
        Do J=1,NSP
          PROP(1,I,J) = DUMPROP(I,J,1)
          PROP(2,I,J) = DUMPROP(I,J,2)
        Enddo
      Enddo
      deallocate(DUMPROP)
      CLOSE (MDPROP,STATUS='KEEP')

      write (6,*)
      hdev = 0.0d0
      rawnum = 0
      do i=1,nsp
        do j=1,nsp
          if ( (dabs(PROP(1,i,j)).GT.ZEROTHRS).OR.
     &         (dabs(PROP(2,i,j)).GT.ZEROTHRS) ) THEN
            if(printlevel.ge.3) then
              write(6,1013) I,J,PROP(1,i,j),PROP(2,i,j)
            endif
            rawnum = rawnum + 1
          endif
          h1=dcmplx(prop(1,j,i),prop(2,j,i))
          h2=dcmplx(prop(1,i,j),prop(2,i,j))
          hdev = hdev + abs(h1 - dconjg(h2))
        enddo
      enddo
      write (6,*) ' *** Raw property integrals for ',NAMEA,':',
     & ' being above ',ZEROTHRS,' : ',rawnum
      write (6,*) ' *** Number of spinors (squared)',nsp,nsp*nsp
      write (6,*) ' *** Deviation from Hermiticity:',hdev
      write (6,*)

      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE
 1013 FORMAT (2I5,F18.10,2X,F18.10)

      ENDIF !Myproc.eq.master
!
!     Distribute the information read from MDPROP
!
!if defined (VAR_MPI)
!     IF (NMPROC .GT. 1) THEN
!        call interface_mpi_BCAST(INOP,1,MASTER,
!    &                  global_communicator)
!        call interface_mpi_BCAST(PROP,2*NSP*NSP,MASTER,
!    &                  global_communicator)
!     END IF
!endif
!
!     Error exit if the integrals could not be read
!
      IF (INOP.EQ.1) THEN
        WRITE (6,*) ' Error reading property ',NAMEA,' in MDPROP'
        WRITE (6,*) ' Error reading property integrals'
        SUCC = .false.
        return
      ENDIF
      IF (INOP.EQ.2) THEN
        WRITE (6,*) ' Property ',NAMEA,' not found in MDPROP'
        WRITE (6,*) ' premature end of file or integrals missing!'
        SUCC = .false.
        return
      ELSE
        WRITE (6,*) ' Property ',NAMEA,' available in MDPROP'
      ENDIF

      if(carith) then
        write(*,*) 'We have complex arithmetic. RCW=',RCW
      else
        write(*,*) 'We have real arithmetic. RCW=',RCW
      endif
!
!__________________________________________________________________
!|   perform array coverage test for OO/OV/VO/VV property arrays
!|   only if printlevel >= 2
!
      if(printlevel.ge.2) then
        i=0
        write(*,*) 'Index coverage test....'
        call idxcoverage(irepspi,nsp,i)
        if(i.ne.0) stop 'Incomplete array coverage!'
        write(*,*) '... passed.'
      endif
!|
!|_________________________________________________________________
!

      APHASE(1) = D1
      APHASE(2) = D0
      PHASE = .FALSE.

      loopcount = 1

!     if(printlevel.ge.2) then
!       write(*,*)
!       write(*,*) 'Symmetry ordering of PROP ints:'
!       write(*,*)
!     endif

    2 continue
      write(*,*)
      write(*,*) '-------------------------------------------------'
      write(*,*) 'Entering phase loop for the',loopcount,'. time'
      write(*,*) '-------------------------------------------------'
      write(*,*)
      loopcount = loopcount + 1
      AMAXI = D0
      AAMAX = 0
      IIMAX = 0
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
!           if(printlevel.ge.2)
!    &        write(*,*) AIREP,' (',repna(airep),') ',arep,
!    &            ' (',repna(arep),') ',irep,' (',repna(irep),')'
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  IF (CARITH) THEN
                      AVO(RCW*AI-1) = PROP(1,AA,II)
                      AVO(RCW*AI)   = PROP(2,AA,II)
                      AOV(RCW*AI-1) = PROP(1,II,AA)   !ov acquire same order as vo !
                      AOV(RCW*AI)   = PROP(2,II,AA)
                  ELSE
                      AVO(AI) = PROP(1,AA,II)*APHASE(1) -
     &                          PROP(2,AA,II)*APHASE(2)
                      AOV(AI) = PROP(1,II,AA)*APHASE(1) -
     &                          PROP(2,II,AA)*APHASE(2)
!       if(printlevel.ge.2)
!    &    write(*,*) '-avo- irepspi#:',aa,ii,'  ccseti#:',
!    &    io(nrep+1)+iv(arep)+a,  io(irep)+i,'  value:',
!    &    avo(ai)
!         write(*,*) '-avo/aov- ',avo(ai),aov(ai)

!         write(*,*) AI,' pai1,pai2,pia1,pia2:',
!    &     PROP(1,AA,II),PROP(2,AA,II),PROP(1,II,AA),PROP(2,II,AA),
!    &     AVO(AI),AOV(AI)

                      AVOI = PROP(1,AA,II)*APHASE(2) +
     &                       PROP(2,AA,II)*APHASE(1)
                      AOVI = PROP(1,II,AA)*APHASE(2) +
     &                       PROP(2,II,AA)*APHASE(1)
                      IF (ABS(AVOI).GT.ABS(AMAXI)) THEN
                         AMAXR = AVO(AI) 
                         AMAXI = AVOI
                         AAMAX = AA
                         IIMAX = II
                      ENDIF
                      IF (ABS(AOVI).GT.ABS(AMAXI)) THEN
                         AMAXR = AOV(AI) 
                         AMAXI = AOVI 
                         AAMAX = II
                         IIMAX = AA
                      ENDIF
                  ENDIF
                   AI = AI + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (.NOT.CARITH.AND.ABS(AMAXI).GT.ACCUR) THEN
         IF (PHASE) THEN
            WRITE (IW,1010) AA,II,AMAXI
            CALL QUIT('USE COMPLEX ARITHMETICS')
         ELSE
            APHASE(1) = AMAXR
            APHASE(2) = AMAXI 
            ANORM = SQRT(AMAXR*AMAXR+AMAXI*AMAXI)
            APHASE(1) = APHASE(1) / ANORM
            APHASE(2) = APHASE(2) / ANORM
            PHASE = .TRUE.
            WRITE (6,1020) NAMEA,APHASE
            WRITE (6,*) 'Phase correction applied, rescaling integrals.'
            WRITE (6,*) 'For this specific properties.'
!
! the following jump instruction leads to a second integral parsing
! and scaling with the proper phase.
!
            GOTO 2
         ENDIF
      ENDIF
C

      IF( (AI-1).ne.dimpov) Call Quit('AOV counting error!')

      write(6,*)
      write(6,'(A)') 'From the AOV integrals we determined a phase'
      write(6,'(A,2F10.5)') 'of ',APHASE(1),APHASE(2)
      write(6,'(A)') 'which is also applied to the AOO and AVV ints.'
      write(6,*)
!
!  generate AOO stream
!
!     write(*,*) '- - - - DOO STREAM:  - - - - - '
      IJ = 1
      DO IJREP = 1,NREP
!       write(*,*)
!       write(*,*) 'DESREP:',ijrep
!       write(*,*)
        DO JREP = 1, NREP
           IREP = MULTB (JREP,IJREP+NREP,2)
           IREP = MULTB (IREP,NREP+1,2) ! form gamma*_i
!         write(*,*) 'irep/jrep,no(irep),no(jrep)',irep,jrep,
!    &               no(irep),no(jrep)
           DO J = 1, NO(JREP)
              JJ = IREPSPI(J,JREP,1)
              DO I = 1, NO(IREP)
                 II = IREPSPI(I,IREP,1)
                 IF (CARITH) THEN
                    AOO(RCW*IJ-1) = PROP(1,II,JJ)
                    AOO(RCW*IJ)   = PROP(2,II,JJ)
                 ELSE
                    AOO(IJ) = PROP(1,II,JJ)*APHASE(1) -
     &                        PROP(2,II,JJ)*APHASE(2)
                 ENDIF

!      write(*,*) 'Integral',IJ,' irep/jrep',irep,jrep,' val:',
!    &      aoo(ij)
                 IJ = IJ + 1

              ENDDO  !I
           ENDDO   !J
        ENDDO  !JREP
      ENDDO  !IJREP
      IF( (IJ-1).ne.dimpoo) Call Quit('AOO counting error!')
!
!  generate AVV stream
!
      AB = 1
      DO ABREP = 1,NREP
        DO BREP = 1, NREP
           AREP = MULTB (BREP,ABREP+NREP,2)
           AREP = MULTB (AREP,NREP+1,2) ! form gamma*_a
           DO B = 1, NV(BREP)
              BB = IREPSPI(B,BREP,2)
              DO A = 1, NV(AREP)
                 AA = IREPSPI(A,AREP,2)
                 IF (CARITH) THEN
                    AVV(RCW*AB-1) = PROP(1,AA,BB)
                    AVV(RCW*AB)   = PROP(2,AA,BB)
                 ELSE
                    AVV(AB) = PROP(1,AA,BB)*APHASE(1) -
     &                        PROP(2,AA,BB)*APHASE(2)
                 ENDIF
                 AB = AB + 1
                 iavvcount = iavvcount + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      IF( (AB-1).ne.dimpvv) Call Quit('AVV counting error!')

      write(*,*)
      write(*,*) 'Total number of AOO ints:     ',dimpoo
      write(*,*) 'Total number of AOV/AVO ints: ',dimpov
      write(*,*) 'Total number of AVV ints:     ',dimpvv

      SUCC = .true.
!
! ********  flush property integrals to disk
!
      co = NAMEA(1:1)
      if(carith) then         ! indicator for aphase that we have carith
        aphase(1) = 0.0d0
        aphase(2) = 0.0d0
      endif
      write(filename,'(A1,A7)') co,'PROPINT'
      write(*,*) 'Flushing properties to ',filename
      open(unit=PROCPROP,file=filename,form='unformatted',
     &     access='sequential',status='unknown')
      write(PROCPROP) dimpoo,dimpov,dimpvo,dimpvv
      write(PROCPROP) aphase(1),aphase(2)
      write(PROCPROP) (aoo(ixx),ixx=1,RCW*dimpoo)
      write(PROCPROP) (aov(ixx),ixx=1,RCW*dimpov)
      write(PROCPROP) (avo(ixx),ixx=1,RCW*dimpvo)
      write(PROCPROP) (avv(ixx),ixx=1,RCW*dimpvv)
      close(PROCPROP)
!
! ********  create an ov stream in the generic OV order
!           and push it on the stack.
!
      allocate(GOV(dimpov*rcw))
      GOV = 0.0d0
      IA = 1
      DO IAREP = 1, NREP
         DO AREP = 1, NREP
            IREP = MULTB (AREP,IAREP+NREP,2)
            IREP = MULTB (IREP,NREP+1,2)
            DO A = 1, NV(AREP)
              AA = IREPSPI(A,AREP,2)
              DO I = 1, NO(IREP)
                 II = IREPSPI(I,IREP,1)
                 IF (CARITH) THEN
                   GOV(RCW*IA-1) = PROP(1,II,AA)
                   GOV(RCW*IA)   = PROP(2,II,AA)
                 ELSE
                   GOV(IA) = PROP(1,II,AA)*APHASE(1) -
     &                       PROP(2,II,AA)*APHASE(2)
                 ENDIF
                 IA = IA + 1
              ENDDO   !I
            ENDDO  !A
         ENDDO  !AREP
      ENDDO  !IAREP

      qst_auxline = AUX_STACKLINE
      if(qstack_push(qst_auxline,dimpov*rcw,GOV).ne.dimpov*rcw) 
     &   stop 'QPE'
      deallocate(GOV)

!
!********************************************
!***   testing the property streams     *****
!********************************************
!
!     write(*,*)
!     write(*,*) '----- Irep order in the property integrals'
!     write(*,*) '----- <A|I> == <V|O> == beta (bosonic)'
!     write(*,*)
!     write(*,*)
!     DO AIREP = 1, NREP
!        write(*,*) 'AIREP ==  *************',airep,' ***************'
!        write(*,*)
!        ipcount = 0
!        DO IREP = 1, NREP
!           AREP = MULTB (IREP,AIREP+NREP,2)
!           AREP = MULTB (AREP,NREP+1,2)
!           if(airep.ne.multb(arep,irep,2)) 
!    &         stop 'Property irep inconsistency!'
!           write(*,*) 'arep,irep',arep,irep
!           DO I = 1, NO(IREP)
!           DO A = 1, NV(AREP)
!             write(*,*) '       A,I:',A,I
!             ipcount = ipcount + 1
!           ENDDO
!           ENDDO
!        ENDDO
!        write (*,*) '<A|I> integrals in this irep:',ipcount
!        write (*,*)
!     ENDDO
!     write(*,*) '--------------------------------------'
!     write (*,*)


!
!********************************************
!***   printing the property integrals  *****
!********************************************
!
      if(printlevel.lt.6) return

      write(*,*) '---------------- AOO -----------------'
      IJ = 1
      DO IJREP = 1,NREP
        DO JREP = 1, NREP
           IREP = MULTB (JREP,IJREP+NREP,2)
           IREP = MULTB (IREP,NREP+1,2) ! form gamma*_i
           DO J = 1, NO(JREP)
              DO I = 1, NO(IREP)
                 irow = io(irep) + i
                 icol = io(jrep) + j
                 IF (CARITH) THEN
                    write(6,'(A,I6,X,2I4,X,2F16.10)') 'int #',ij,
     &                   irow,icol,AOO(2*IJ-1),AOO(2*IJ)
                 ELSE
                    write(6,'(A,I6,X,2I4,X,F16.10)') 'int #',ij,
     &                   irow,icol,AOO(IJ)
                 ENDIF
                 IJ = IJ + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO
!
!  AOV type
!
!     write(*,*) '---------------- AOV -----------------'
!     IA = 1
!     DO IAREP = 1, NREP
!        DO AREP = 1, NREP
!           IREP = MULTB (AREP,IAREP+NREP,2)
!           IREP = MULTB (IREP,NREP+1,2)
!           write(*,*) '   ----> IA,I,Arep:',iarep,irep,arep
!           DO A = 1, NV(AREP)
!             DO I = 1, NO(IREP)
!                irow = io(irep) + i
!                icol = io(nrep) + iv(arep) + a
!                 IF (CARITH) THEN
!                   write(6,'(A,I6,X,2I4,X,2F16.10)') 'j/a: ',IA,
!    &                    irow,icol,AOV(2*IA-1),AOV(2*IA)
!                 ELSE
!                   write(6,'(A,I6,X,2I4,X,F16.10)') 'j/a: ',IA,
!    &                    irow,icol,AOV(IA)
!                 ENDIF
!                 IA = IA + 1
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
!
!  AVO type
!
      write(*,*)
      write(*,*) '--------- AVO/AOV (same symmetry order!) --------'
      write(*,*) '--------- this makes eigenvector contraction'
      write(*,*) '--------- much easier!'
      write(*,*)
      AI = 1
      DO AIREP = 1, NREP
         DO IREP = 1, NREP
            AREP = MULTB (IREP,AIREP+NREP,2)
            AREP = MULTB (AREP,NREP+1,2)
!           write(*,*) '   ----> AI,A,Irep:',airep,arep,irep
            DO I = 1, NO(IREP)
               DO A = 1, NV(AREP)
                 irow = io(nrep) + iv(arep) + a
                 icol = io(irep) + i
                  IF (CARITH) THEN
                    write(6,'(A,I6,X,2I4,X,2F16.10,A,2F16.10)')
     &              'int # ',ai,irow,icol,AVO(2*AI-1),AVO(2*AI),
     &              '  ||  ',AOV(2*AI-1),AOV(2*AI)
                  ELSE
                    write(6,'(A,I6,X,2I4,X,F16.10,A,F16.10)') 'int # ',
     &              ai,irow,icol,AVO(AI),'  ||  ',AOV(AI)
                  ENDIF
                  AI = AI + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO

!
!  AVV type
!
      write(*,*)
      write(*,*) '------------- AVV ----------------'
      write(*,*)
      AB = 1
      DO ABREP = 1,NREP
        DO BREP = 1, NREP
           AREP = MULTB (BREP,ABREP+NREP,2)
           AREP = MULTB (AREP,NREP+1,2) ! form gamma*_a
           DO B = 1, NV(BREP)
              BB = IREPSPI(B,BREP,2)
              DO A = 1, NV(AREP)
                 AA = IREPSPI(A,AREP,2)
                 irow = io(nrep) + iv(arep) + a
                 icol = io(nrep) + iv(brep) + b
                 IF (CARITH) THEN
                    write(6,'(A,I6,X,2I4,X,2F16.10)') 'int # ',
     &              ab,irow,icol,AVV(2*AB-1),AVV(2*AB)
                 ELSE
                    write(6,'(A,I6,X,2I4,X,F16.10)') 'int # ',
     &              ab,irow,icol,AVV(AB)
                 ENDIF
                 AB = AB + 1
              ENDDO
           ENDDO
        ENDDO
      ENDDO

      RETURN

 1000 FORMAT (/' Read integrals of type ',A8,
     & A8,' storage info : ',A8)
 1010 FORMAT (/' Largest imaginary part of matrix element',2I5,F10.2)
 1020 FORMAT (/' Property ',A8,' scaled with phase factor ',2F10.4)
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE IDXCOVERAGE(IREPSPI,NSP,ok)

      IMPLICIT none
!     IMPLICIT INTEGER (A-Z)
!
! perform array coverage test for the property integral arrays
! where absolute Kramers-unique spinor numbers are related to relative
! symmetry-grouped arrays.
!
!----------------Common Blocks-----------------------------------------
!
#include "../relccsd/symm.inc"
!provides MXREP and MULTB
!
!----------------Calling variables---------------------------------------------
!
      INTEGER                       :: NSP,ok
      INTEGER                       :: IREPSPI(NSP,MXREP,2)
!
!----------------Local variables---------------------------------------------
!
      integer,allocatable           :: acv(:,:)
      integer                       :: irep,jrep,arep,brep
      integer                       :: airep,ijrep,abrep
      integer                       :: i,ii,j,jj,a,aa,b,bb
      integer                       :: rem
!
!----------------Executable code ---------------------------------------------
!
      write(*,*) '*********'
      write(*,*) '*** Array coverage test triggered.'
      write(*,*) '*********'

      ok = 0

      allocate(acv(nsp,nsp))
      acv = 1

      write(*,*) '*** testing VO/OV...'
      DO AIREP = 1, NREP
        DO IREP = 1, NREP
          AREP = MULTB (IREP,AIREP+NREP,2)
          AREP = MULTB (AREP,NREP+1,2)
!         write(*,*) 'AIREP,AREP,IREP',airep,arep,irep
            DO I = 1, NO(IREP)
               II = IREPSPI(I,IREP,1)
               DO A = 1, NV(AREP)
                  AA = IREPSPI(A,AREP,2)
                  if(AA.gt.nsp) then
                    stop 'Internal error IDXCoverage' 
                  endif
                  if(II.gt.nsp) then
                    stop 'Internal error IDXCoverage' 
                  endif
                  acv(AA,II)=0
                  acv(II,AA)=0
               ENDDO
            ENDDO
        ENDDO
      ENDDO

      write(*,*) '*** testing OO...'
      DO IJREP = 1,NREP
        DO JREP = 1, NREP
           IREP = MULTB (JREP,IJREP+NREP,2)
           IREP = MULTB (IREP,NREP+1,2) ! form gamma*_i
           DO J = 1, NO(JREP)
              JJ = IREPSPI(J,JREP,1)
              DO I = 1, NO(IREP)
                 II = IREPSPI(I,IREP,1)
                 acv(II,JJ)=0
              ENDDO
           ENDDO
        ENDDO
      ENDDO

      write(*,*) '*** testing VV...'
      DO ABREP = 1,NREP
        DO BREP = 1, NREP
           AREP = MULTB (BREP,ABREP+NREP,2)
           AREP = MULTB (AREP,NREP+1,2) ! form gamma*_a
           DO B = 1, NV(BREP)
              BB = IREPSPI(B,BREP,2)
              DO A = 1, NV(AREP)
                 AA = IREPSPI(A,AREP,2)
                 acv(AA,BB)=0
              ENDDO
           ENDDO
        ENDDO
      ENDDO

      rem=sum(acv)
      if(rem.ne.0) then
         write(*,*) ' --->>> ERROR <<<---'
         write(*,*) ' --->>> index coverage: ',rem,
     &       ' should be zero !'
         ok = -1
      endif
      deallocate(acv)
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
