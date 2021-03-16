!      Copyright (c) 2014 by the authors of DIRAC.
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetno */
      SUBROUTINE RGETNO(CREF,OCCNO,UNO,CMO,NATORB,MOISNO,save_1pdens,
     &                  IPRINT,WORK,KFREE,LFREE)
      use memory_allocator
C***********************************************************************
C
C     Calculate natural orbital occupation numbers and, 
C     if NATORB is true, transform orbital to natural orbitals.
C
C     Input:
C        CREF   - CI reference vector
C        CMO    - MO orbitals
C        NATORB - transform CMO to natural orbitals
C
C     Output:
C        OCCNO  - natural orbital occupation numbers
C        UNO    - transformation matrix
C        MOISNO - MO orbitals were already NO orbitals at entry.
C        CMO    - if NATORB, then NO orbitals, else unchanged
C
C     Written by J. Thyssen - Dec 21 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      DIMENSION CREF(*)
      DIMENSION OCCNO(*), UNO(NASHT,NASHT,*)
      DIMENSION CMO(*)
      DIMENSION WORK(*)
C
      LOGICAL NATORB, MOISNO, EX, FILEOPEN, save_1pdens
      LOGICAL FNDLAB, FND
#if defined MCSCF_DEBUG
      real*8,  allocatable :: onedens(:,:,:)
      real*8,  allocatable :: cref_x(:)
      integer, allocatable :: mjvec(:)
#endif
      integer, allocatable :: reorder_fermion(:)
      integer, allocatable :: reorder_list(:,:)
      real(8), parameter   :: dm2 = -2.0d0
C
      CALL QENTER('RGETNO')
      KFRSAV = KFREE

      IF (NATORB) THEN
         CALL QUIT('FATAL ERROR in RGETNO:'//
     &      ' Input option NATORB not implemented yet.')
      END IF
C
C     Calculate DV density matrix
C     ---------------------------
C
      CALL MEMGET2('REAL','DV',KDV,N2ASHXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','EIG',KEIG,NASHT,WORK,KFREE,LFREE)

#if defined MCSCF_DEBUG

      iread_from_mc_file = 0
      call alloc(cref_x,nzconfq)
      cref_x = 0
      luni = -1
      inquire(file='KRMCOLD',exist=ex,opened=fileopen,number=luni)
      if(ex.and.fileopen) then
        rewind(luni)
        call reakrmc(luni,'CREF    ',cref_x,nzconfq)
        iread_from_mc_file = 1
      else if( ex.and.(.not. fileopen ))then
        luni = 99
        call opnfil(luni,'KRMCOLD','OLD    ','RGETNO')
        rewind(luni)
        call reakrmc(luni,'CREF    ',cref_x,nzconfq)
        iread_from_mc_file = 1
        close(luni)
      else if( .not. ex )then
        call dcopy(nzconfq,cref,1,cref_x,1)
      end if
!     call wrtmatmn(cref_x,1,nzconfq,1,nzconfq,lupri)

      if( iread_from_mc_file .gt. 0)then
        write(lupri,'(/1X,A,I16)')
     &  '(RGETNO) NO calculation:'//
     &  ' cref read from KRMCOLD with # of dets =',nzconfq
      end if

      call rmakdv(work(kdv),cref_x,work,kfree,lfree)
      call dealloc(cref_x)
#else
      CALL RMAKDV(WORK(KDV),CREF,WORK,KFREE,LFREE)
#endif
!
!     save 1-particle density matrix on file
!     TODO: expand to full NORBT x NORBT x NZ matrix
!
      if(save_1pdens)then
        lunii = 98
        CALL OPNFIL(lunii,'MCDENS','UNKNOWN','STORE_')
        CALL NEWLAB('MC  DENS',lunii,LUPRI)
        CALL WRITT(lunii,NASHT*NASHT*NZ,WORK(KDV))
        CLOSE(Unit=lunii,Status='KEEP')
      end if
C
C     Calculate max off-diagonal elements of density matrix
C
      DVMAX = D0
      DO I = 2, NASHT
         DO J = 1, I - 1
            DVMAX = MAX(DVMAX, ABS(WORK(KDV + (J-1) + (I-1)*NASHT)))
         END DO
      END DO
      DO IZ = 2, NZ
         DO I = 1, NASHT
            DO J = 1, I
               DVMAX = MAX(DVMAX,
     &              ABS(WORK(KDV + (J-1) + (I-1)*NASHT+(IZ-1)*N2ASHX)))
            END DO
         END DO
      END DO
      IF (DVMAX .LE. THRZER) THEN
         MOISNO = .TRUE.
         IF (IPRINT .GE. 0) WRITE(LUPRI,'(1X,A)')
     &      '(RGETNO) The orbitals are already natural orbitals.'
         CALL DZERO(UNO,N2ASHXQ)
         CALL DUNIT(UNO,NASHT)
         DO I = 1, NASHT
            WORK(KEIG+I-1) = WORK(KDV+(I-1)+(I-1)*NASHT)
         END DO
      ELSE
         MOISNO = .FALSE.
      END IF
!
#if defined MCSCF_DEBUG
      call alloc(onedens,nasht,nasht,nz)
      onedens = 0
      call dcopy(n2ashxq,work(kdv),1,onedens,1)
      do i = 1, NFSYM
         write(lupri,'(/1X,A,I3)')
     &   '(RGETNO) one-electron density matrix for NO calculation'//
     &    '(active orb.) for symmetry ',I
         call prqmat(onedens(1+iash(i),1+iash(i),1),nash(i),nash(i),
     &               nasht,nasht,nz,ipqtoq(1,0),lupri)
      end do
      call dealloc(onedens)
#endif
C
C     Diagonalize density matrix
C     --------------------------

!     zero out coupling blocks between non-equivalent GA spaces
      if(.not.save_1pdens)then
        if(NATORB)then
          if(ngas_dc > 2) call gasno(work(kdv))
        end if
      end if

      IF (.NOT. MOISNO) THEN
         IF (NZ .EQ. 1) THEN
            IJOB = 1
            IORDER = 0
            IPACK = 0
            CALL RSJACO(NASHT,NASHT,NASHT,WORK(KDV),WORK(KEIG),
     &           IJOB,IORDER,IPACK,UNO)
         ELSE
            MATZ = 1 ! eigenvalues and eigenvectors
            IF ( NATOLCR ) MATZ = 0 ! only eigenvalues
            CALL QDIAG(NZ,NASHT,WORK(KDV),NASHT,NASHT,
     &           WORK(KEIG),MATZ,UNO,NASHT,NASHT,WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/A,I5)')
     &            ' *** ERROR in RGETNO ***: '//
     &            'QDIAG failed with error code',IERR
               CALL QUIT('*** ERROR in RGETNO ***')
            END IF
         END IF
C
C        FIXME: select positive phase of transformation vectors
C
      END IF
C
      JOCC = 0
      JACT = 0
      JACT_TOT = 0
      CALL MEMGET('REAL',KOCCNO2,NORBT  ,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KOCCNO2),NORBT)
*
*     FIXME: right order of active orbitals if NZ == 2 and 
*            NO occ. numbers are calculated with LUCIAREL
*
      IF( NATOLCR .and. NZ .eq. 2 ) THEN
        WRITE(LUPRI,'(/A/)') 
     & ' *** WARNING *** QDIAG sorts eigenvalues - NO occupation'//
     & ' might appear in wrong order.'
        DO I = 1, NFSYM
           DO J = 1, NISH(I)
              JOCC = JOCC + 1
              OCCNO(JOCC) = D2
           END DO
           JACT_TOT = JACT_TOT + NASH(I)
           DO J = 1, NASH(I)
              JOCC = JOCC + 1
              JACT = JACT + 1
              OCCNO(JOCC) = WORK(KEIG + JACT_TOT - JACT )
C             OCCNO(JOCC) = WORK(KEIG + JACT - 1)
           END DO
        END DO
      ELSE
        DO I = 1, NFSYM
           DO J = 1, NISH(I)
              JOCC = JOCC + 1
              OCCNO(JOCC) = D2
           END DO
           DO J = 1, NASH(I)
              JOCC = JOCC + 1
              JACT = JACT + 1
              OCCNO(JOCC) = WORK(KEIG + JACT - 1)
           END DO
        END DO
      END IF
C
      IF ( IPRINT .GE. 1 .or. NATOLCR ) THEN
         IF( NATOLCR )THEN
            WRITE(LUPRI,'(//A/A/A)')
     &' ***************************************************',
     &' *************** CI natural orbitals ***************',
     &' ***************************************************'
         WRITE(LUPRI,'(//A,I3,/A)')
     &        '   Natural orbital occupation numbers for state:', 
     &            IDENSLR_STATE,
     &        '   ------------------------------------------------'
         ELSE
            WRITE(LUPRI,'(//A/A/A)')
     &' ***************************************************',
     &' *************** MCSCF natural orbitals ************',
     &' ***************************************************'
         END IF
C        
C        save occupancies for later use 
         DO I = 1, NFSYM
           IOFF1 = IOCC(I) + 1
           IOFF2 = IORB(I)
           do jj = 1, NPSH(i)
             WORK(KOCCNO2+IOFF2+jj-1) = dm2
           end do
           IOFF2 = IORB(I) + NPSH(i) + 1
           IF(NOCC(I).gt.0)THEN
             CALL DCOPY(NOCC(I),OCCNO(IOFF1),1,WORK(KOCCNO2+IOFF2-1),1)
           END IF
         END DO
         DO I = 1, NFSYM
            NOCCI = NOCC(I)
            IF (NOCCI .EQ. 0) THEN
               WRITE(LUPRI,9001) I
            ELSE
               NISHI = NISH(I)
               NASHI = NASH(I)
               IOCCI = IOCC(I)
               IOCCA = IOCCI + NISHI
               OCCSUM = DSUM(NOCCI,OCCNO(IOCCI+1),1)
               IF (NASHI .GT. 0) THEN
                  WRITE(LUPRI,9002) I, (OCCNO(IOCCA+J),J=1,NASHI)
                  OCCSUMA = DSUM(NASHI,OCCNO(IOCCA+1),1)
               ELSE
                  OCCSUMA = 0.0D0
               END IF

               WRITE(LUPRI,9003) OCCSUM-OCCSUMA,OCCSUMA,OCCSUM

               IF (NASHI .GT. 0) THEN
                  CALL ORDER3(WORK(KFREE),OCCNO(IOCCI+1),
     &               NOCCI,NOCCI,0,-1)
                  WRITE(LUPRI,9004) I, (OCCNO(IOCCA+J),J=1,NASHI)
               END IF
            END IF
         END DO
      END IF
C
 9001 FORMAT(/'   Symmetry',I3,'.   No occupied orbitals.')
 9002 FORMAT(/'   Symmetry',I3,'.   Occupation of active orbitals:'
     &      //,(5F14.9))
 9003 FORMAT(/'   Total',F16.9,' +',F16.9,' =',F16.9)
 9004 FORMAT(/'   Symmetry',I3,'. ** REORDERED DESCENDING OCCUPATION **'
     &      //,(5F14.9))
C
      IF ( IPRINT .GE. 100 .AND. .NOT. MOISNO ) THEN
         DO I = 1, NFSYM
            WRITE(LUPRI,'(/A,I2)')
     &           ' (RGETNO) Transformation to NO '//
     &           '(active orb.) for symmetry',I
            CALL PRQMAT(UNO(1+IASH(I),1+IASH(I),1),NASH(I),NASH(I),
     &           NASHT,NASHT,NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
C
C     Transform to natural orbitals
C     -----------------------------
C
C
      IF (.NOT. NATORB) THEN
        IF ( IPRINT .GE. 1 )
     &    WRITE(LUPRI,'(/A)')
     &    ' (RGETNO) active orbitals are NOT transformed '//
     &    'to natural orbitals.'
      ELSE


C       save MO coefficients
        CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
        CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
        IF(.NOT. MOISNO) THEN
          IF(IPRINT.GE.1)
     &      WRITE(LUPRI,'(/1X,A)')
     &        '(RGETNO) active orbitals are transformed '//
     &        'to natural orbitals.'
C
C
          DO 10 I = 1, NFSYM
             IF (NASH(I) .EQ. 0) GO TO 10
C
             IF ( IPRINT .GE. 20 ) THEN
                WRITE(LUPRI,'(/1X,A,I3)')
     &               '(RGETNO) Before transformation: '//
     &               'MO coefficients, symmetry ',I
                CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
             END IF
C
C            Transform active orbitals:
C
C            CMO(NO) = CMO * UNO
C
             NBASI = NFBAS(I,0)
             NASHI = NASH(I)
             NORBI = NORB(I)
             JUNO = 1 + IASH(I)
             ICMOA = 1 + ICMOQ(I) + (NPSH(I)+NISH(I)) * NBASI
             CALL QGEMM(NBASI,NASHI,NASHI,D1,
     &            'N','N',IPQTOQ(1,0),
     &            WORK(KCMO+ICMOA-1),NBASI,NORBI,NZ,
     &            'N','N',IPQTOQ(1,0),
     &            UNO(JUNO,JUNO,1),NASHT,NASHT,NZ,
     &            D0,IPQTOQ(1,0),
     &            CMO(ICMOA),NBASI,NORBI,NZ)
C
             IF ( IPRINT .GE. 20 ) THEN
                WRITE(LUPRI,'(/1X,A,I3)')
     &               '(RGETNO) After transformation: '//
     &               'MO coefficients, symmetry ',I
                CALL PRQMAT(CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),
     &               NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
             END IF
C
 10       CONTINUE
          MOISNO = .TRUE.
        END IF
 
!       write NOs and occupancies to file with label MCCINATO and 
!       MCNATOCC, respectively.

        CALL MEMGET('INTE',KMJVEC,NORBT,WORK,KFREE,LFREE)
        CALL IZERO(WORK(KMJVEC),NORBT)
        LUNI            = -1
        mypath_write_no = -1

!       retrieve boson/mj-value information if applciable
        INQUIRE(FILE='KRMCSCF',EXIST=EX,OPENED=FILEOPEN,NUMBER=LUNI)
        IF(EX)THEN
          IF(FILEOPEN)THEN
            REWIND(LUNI)
            mypath_write_no = 1
          ELSE
            LUNI = 99
            CALL OPNFIL(LUNI,'KRMCSCF','UNKNOWN','RGETNO')
            REWIND(LUNI)
            mypath_write_no = 2
          END IF
          if(spinfr_krmc)then
            CALL IREAKRMC(LUNI,'IBEIG   ',WORK(KMJVEC),NORBT)
          else if(LINEAR)THEN
            CALL IREAKRMC(LUNI,'MJVEC   ',WORK(KMJVEC),NORBT)
          END IF
        ELSE IF(.NOT. EX)THEN
          LUNI = 99
          CALL OPNFIL(LUNI,'KRMCSCF','UNKNOWN','RGETNO')
          CALL NEWLAB('*RGETNO ',LUNI,LUPRI)
          REWIND(LUNI)
          mypath_write_no = 3
        END IF

!       save reordered NOs if requested
        if(save_reordered_nos)then
          WRITE(LUPRI,'(/A)') '  saving REORDERED NOs on KRMCSCF!'
          WRITE(LUPRI,'(A/)') '         ---------                '

          call alloc(reorder_fermion,2)
          call alloc(reorder_list,nesht,2)

          reorder_fermion(1:2) = nesh(1:2)

!         restore unsorted NO numbers in occno
          call dcopy(norbt,WORK(KOCCNO2),1,occno,1)
!         create list for reordering in descending order according to occupation numbers
!         ------------------------------------------------------------------------------
          do i = 1, nfsym

!           write(lupri,'(/a,i4,a,i4)') 
!    &                  '  reorder from active index... ',
!    &                     NISH(I)+1, ' ... to ... ', 
!    &                     NISH(I)+NASH(I)

!           a. initialize inactive, active and secondary shells to no reordering

            do jj = 1, nesh(i)
              reorder_list(jj,i) = jj
            end do

!           -------------
!           b. reorder active shells after no. isvrono(i)
!           -------------

            iash_ro = nish(i) + isvrono(i)    ! off-set in reorder_list to orbitals to reorder
            nash_ro = nash(i) - isvrono(i)    ! no. of orbitals to reorder
            jash_ro = iorb(i)+npsh(i)+iash_ro ! off-set to first to reorder in total orbital list

            myturn = 1
            do 
!             write(lupri,*) 'myturn,nash(i) ', myturn,nash(i)
!             exit condition
              if (myturn > nash_ro) exit

              mymax  = idamax(nash_ro,occno(jash_ro+1),1)
!             write(lupri,*) '  max index and value',
!    &        mymax, occno(iorb(i)+npsh(i)+nish(i)+mymax)

!             new max index
              reorder_list(iash_ro+myturn,i) = iash_ro + mymax
              occno(jash_ro+mymax) = 0.0d0

              myturn = myturn + 1
            end do

!           do jj = 1, nash(i)
!             write(lupri,*) ' active electronic shell # ',
!    &                         nish(i)+jj,
!    &                         ' goes to ',
!    &                         reorder_list(nish(i)+jj,i)
!           end do

          end do ! i = 1,nfsym

!         write(lupri,*) 'reorder list'
!         do i = 1, nfsym
!           write(lupri,*) ' fermion sym # ',i
!           do jj = 1, nesh(i)
!             write(lupri,*) ' electronic shell # ',jj,' goes to ',
!    &                         reorder_list(jj,i)
!           end do
!         end do

          CALL REORD(CMO,WORK(KOCCNO2),WORK(KMJVEC),
     &               reorder_list,reorder_fermion,nesht)

          call dealloc(reorder_list)
          call dealloc(reorder_fermion)
        end if

        REWIND(LUNI)
        FND = FNDLAB('SODLABEL',LUNI)
        BACKSPACE LUNI
        CALL NEWLAB('SODLABEL',LUNI,LUPRI)
        REWIND(LUNI)

        if(spinfr_krmc)then
          CALL IWRTKRMC(LUNI,'IBEIG   ',WORK(KMJVEC),NORBT)
        else if(linear)then
          CALL IWRTKRMC(LUNI,'MJVEC   ',WORK(KMJVEC),NORBT)
        end if
        CALL WRTKRMC(LUNI,'MCCINATO',CMO,NCMOTQ)
        CALL WRTKRMC(LUNI,'MCNATOCC',WORK(KOCCNO2),NORBT)

        if(.not.FILEOPEN) CLOSE(LUNI,STATUS='KEEP')

        WRITE(LUPRI,'(/A,I2)')
     &  '  Natural orbitals saved with label MCCINATO on KRMCSCF'//
     &  '  -  path ',mypath_write_no
C
C       copy back MO coefficients
        CALL DZERO(CMO,NCMOTQ)
        CALL DCOPY(NCMOTQ,WORK(KCMO),1,CMO,1)
      ENDIF
 100  CALL MEMREL('RGETNO.transform',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('RGETNO')
      RETURN
      END

!*******************************************************************************
      SUBROUTINE GASNO(DV)
!
!
! Zero blocks in one-electron density matrix DV
! connecting different GAS spaces (RAS1, RAS2, RAS3).
!
! stefan: heavily inspired by RASNO written by Hans Joergen Aa. Jensen in Dalton
!
! currently works only for GASCISD or similar, e.g. no cv-correlation with an 
! additional GAS space to take care of for zeroing... FIXME. 
!
!
      use memory_allocator
      implicit none

#include "priunit.h"
#include "dcborb.h"
#include "dcbopt.h"
#include "dgroup.h"
#include "dcbbas.h"

      integer , parameter     :: dp = 8
      real(dp), parameter     :: d0 = 0.0d0

      real(dp), intent(inout) :: dv(nasht,nasht,nz)

      integer,  allocatable   :: iactyp(:)
      integer                 :: i, j, iactj, iacti, ifsym, iz

      call alloc(iactyp,nasht)
      iactyp = 0

!     alter type only for GAS 1 !!!
      do ifsym = 1, nfsym
        do j = 1, ngsh(ifsym,1)
          iactyp(iash(ifsym)+j) = -1
        end do
      end do

      do ifsym = 1, nfsym
        do 200 j = 1, nash(ifsym)
          iactj = iactyp(j+iash(ifsym))
          do 100 i = 1, nash(ifsym)
            iacti = iactyp(i+iash(ifsym))
            if(iacti /= iactj)then
              write(lupri,*) ' zero DV element between orbital '//
     &        ' # ',j + iash(ifsym), ' and orbital # ',
     &              i + iash(ifsym), ' element: ',
     &        dv(i + iash(ifsym), j + iash(ifsym), 1)
              do iz = 1, nz
                dv(i + iash(ifsym), j + iash(ifsym), iz) = d0
                dv(j + iash(ifsym), i + iash(ifsym), iz) = d0
              end do
            end if
  100     continue
  200   continue
      end do

      call dealloc(iactyp)

      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmccno */
      SUBROUTINE KRMCCNO(KEYCNO,RTRACI,CREF,CMO,IPRINT,EIG,
     &                   WORK,KFREE,LFREE)
C***********************************************************************
C
C     Transform orbitals to natural and/or Fock type orbitals.
C
C     Input:
C        CREF - CI reference vector
C        CMO  - MO orbitals
C
C     Output:
C        CREF, CMO - transformed quantities.
C        EIG,  eigenvalues of blocked FD matrix (if TRTOFD true)
C              primarily for use in QCDHF
C
C     Written by J. Thyssen - Mar 4 2000
C     Last revision :
C
C                S. Knecht  - Jul 16 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
#include "thrzer.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbkrmc.h"
C
      DIMENSION CREF(*), CMO(*), EIG(*), WORK(*)
      CHARACTER KEYCNO*6
      LOGICAL RTRACI
C
      LOGICAL TRTOFD, TRTOCO, TRTONO
      DIMENSION J1MO(2), N1MO(2), J2MO(2), N2MO(2)
      integer iimax
C
      CALL QENTER('KRMCCNO')
      KFRSAV = KFREE
      
C
      IF (KEYCNO .EQ. 'ONLYNO' .OR. KEYCNO .EQ. 'NATORB') THEN
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,2A)')
     &           '(KRMCCNO) No transformation.'
            GOTO 2000
         END IF
         IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,2A)')
     &        '(KRMCCNO) Transformation of active ',
     &        'orbitals to natural orbitals.'
         TRTOFD = .FALSE.
         TRTOCO = .FALSE.
         TRTONO = .TRUE.
      ELSE IF (KEYCNO .EQ. 'ONLYFD' .OR. KEYCNO .EQ. 'ONLYFC' .OR.
     &         KEYCNO .EQ. 'ONLYCO' .OR. KEYCNO .EQ. 'CANORB') THEN
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,A/,A/,A/)')
     &           '(KRMCCNO) MCSCF iteration - transformation of  ',
     &           '          positronic/inactive/virtual orbitals ',
     &           '          to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .FALSE.
         ELSE
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,2A)')
     &           '(KRMCCNO) Transformation of all ',
     &           'orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .TRUE.
            TRTONO = .FALSE.
         END IF
      ELSE
         IF (.NOT. RTRACI) THEN
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,3A)')
     &           '(KRMCCNO) Transformation of           ',
     &           '          positronic/inactive/virtual ',
     &           '          orbitals to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .FALSE.
         ELSE
            IF (IPROPT .GE. 0) WRITE(LUPRI,'(/1X,4A)')
     &           '(KRMCCNO) Transformation of active ',
     &           'orbitals to natural orbitals and ',
     &           'positronic/inactive/virtual orbitals',
     &           'to Fock type orbitals.'
            TRTOFD = .TRUE.
            TRTOCO = .FALSE.
            TRTONO = .TRUE.
         END IF
      END IF
C
C
C     ***************************************
C     *** Transform to Fock type orbitals ***
C     ***************************************
C
C     1. Calculate one-electron density matrix
C     ----------------------------------------
C
C     Needed for natural orbitals and/or FV.
C
      CALL MEMGET2('REAL','DV',KDV,N2ASHXQ,WORK,KFREE,LFREE)
      IF (NASHT .GE. 1) THEN
         CALL RMAKDV(WORK(KDV),CREF,WORK,KFREE,LFREE)
      END IF
C
      IF ( .NOT. (TRTOFD .OR. TRTOCO) ) GOTO 1000
C
      CALL MEMGET2('REAL','EIGVC',KEIGVC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','FD',KFD,N2ORBXQ,WORK,KFREE,LFREE)
C     ... initialize
      CALL DZERO(WORK(KEIGVC),N2ORBXQ)
      CALL DZERO(WORK(KFD),   N2ORBXQ)
C
      IF (NASHT .GT. 0) THEN
         LFV = N2ORBXQ
      ELSE
         LFV = 0
      END IF
      CALL MEMGET2('REAL','FV',KFV,LFV,WORK,KFREE,LFREE)
C
C
C     2. Calculate Fock matrices
C     --------------------------
C
      KFC = KFD
      CALL RFCKMAT(WORK(KDV),CMO,EIX,EAX,WORK(KFC),WORK(KFV),
     &     WORK(KFREE),LFREE)
C
C     Calculate FD: FD = FC + FV
C
      IF (NASHT .GT. 0)
     &   CALL DAXPY(N2ORBXQ,1.0D0,WORK(KFV),1,WORK(KFD),1)
C
      IF (TRTOCO) THEN
C
C        Extract active-active part of FD for RTRACI.
C
         CALL RGETAC(WORK(KFD),WORK(KDV),IPROPT)
C        CALL RGETAC(FD,FDAC,IPRINT)
C
         IF (.NOT. TRTOFD) THEN
            CALL MEMREL('KRMCCNO.after RGETAC for FCAC',
     &                  WORK,1,KEIGVC,KFREE,LFREE)
            GO TO 1000
         END IF
      END IF
C
CSK   CALL QUIT('KRMCCNO: TRTOFD not implemented correctly yet')
C     hjaaj Feb 2004 TODO: because QDIAG sorts eigenvalues code
C     must be modified to only call QDIAG with a block at a time
C     (p-p block, i-i block, s-s block) instead of zeroing inter-shell
C     couplings and active-active block.
C
CSK   DONE! July 2008
C
      ipropt_save = ipropt
!     ipropt      = 10 ! debug

      IF (IPROPT .GE. 5) THEN
         DO I = 1, NFSYM
            WRITE(LUPRI,*) 'KRMCCNO: FOCK matrix for fsym = ',I
            CALL PRQMAT(WORK(KFD+I2ORBX(I,I)),NORB(I),NORB(I),
     &                  NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         ENDDO
      END IF
!
      CALL MEMGET('REAL',kibeig,norbt,WORK,KFREE,LFREE)

!     read symmetry informations for orbitals from file KRMCSCF
!     for the two cases that are currently implemented:
 
!     a. spinfree MCSCF
      if(spinfr_krmc)then
        rewind(lukrmc)
        call ireakrmc(lukrmc,'IBEIG   ',work(kibeig),norbt)
!     b. MCSCF in linear symmetry
      else if(opt_chckjz)then
        rewind(lukrmc)
        call ireakrmc(lukrmc,'MJVEC   ',work(kibeig),norbt)
!     c. initialize to 0
      else
        call izero(work(kibeig),norbt)
      end if
C
C     3. Diagonalize Fock matrix 
C     --------------------------
C

!     extract sub-blocks of the Fock matrix in terms of j_z values / boson
!     irreps for linear symmetry / spinfree calculations
      if(spinfr_krmc)then 
        max_loop = nbsym
      else if(opt_chckjz)then
        max_loop = iimax(norbt,work(kibeig),1)
      else
        max_loop = 1
      end if
C
C     3.1 i-i block
C     -------------
C
      IF ( NISHT .gt. 0 ) THEN
C
C        ... allocate for i-i block diagonalization of Fock-matrix
         CALL MEMGET('REAL',KFDI,   N2ISHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVCI,N2ISHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGI,  NISHT    ,WORK,KFREE,LFREE)
C        ... initialize
         CALL DZERO(WORK(KFDI),   N2ISHX*NZ)
         CALL DZERO(WORK(KEIGVCI),N2ISHX*NZ)
         CALL DZERO(WORK(KEIGI),  NISHT)
!
!        extract i-i part of WORK(KFD)

         do indx_loop = 1, max_loop

           if(spinfr_krmc)then
!            get boson irrep
             is_symm = indx_loop - 1
           else if(opt_chckjz)then
!            get mj-value
             is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
           else
             is_symm = 0
           end if
           CALL MATGAT_symm(WORK(KFD),NORBT,NORBT,WORK(KFDI),NISHT,
     &                      NISHT,IDXI2G,work(kibeig),is_symm,NZ)
         end do

         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*) 'KRMCCNO: inactive FOCK matrix (extracted)'
           CALL PRQMAT(WORK(KFDI),NISHT,NISHT,
     &                 NISHT,NISHT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
!
!        diagonalize
         CALL DIAG_piv2fock('I',WORK(KFDI),WORK(KEIGI),WORK(KEIGVCI),
     &                      IDXI2G,work(kibeig),WORK(KFREE),LFREE)

!        insert new i-i matrix block into transformation matrix
         do indx_loop = 1, max_loop

           if(spinfr_krmc)then
!            get boson irrep
             is_symm = indx_loop - 1
           else if(opt_chckjz)then
!            get mj-value
             is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
           else
             is_symm = 0
           end if
           CALL MATSCT_symm(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCI),
     &                      NISHT,NISHT,IDXI2G,work(kibeig),is_symm,NZ)
         end do

         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*)
     &     'KRMCCNO: transformation matrix (i-i)'
           CALL PRQMAT(WORK(KEIGVC),NORBT,NORBT,
     &                 NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
 
         CALL MEMREL('KRMCCNO.after i-i',WORK,1,KFDI,KFREE,LFREE)
      END IF

C
C     For the empty (virtual) orbitals:
C     calculate modified FD: FD' = FC + DKRMC_MVOFAC*FV = FD + (DKRMC_MVOFAC-1)*FV
C
      FAC_FV = DKRMC_MVOFAC - 1.0D0
      IF (NASHT .GT. 0 .AND. FAC_FV .NE. 0.0D0)
     &   CALL DAXPY(N2ORBXQ,FAC_FV,WORK(KFV),1,WORK(KFD),1)

C
C     3.2 p-p block
C     -------------
C
      IF ( NPSHT .gt. 0 ) THEN
C
C        ... allocate for p-p block diagonalization of Fock-matrix
         CALL MEMGET('REAL',KFDP,   N2PSHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVCP,N2PSHX*NZ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGP,  NPSHT    ,WORK,KFREE,LFREE)
C        ... initialize
         CALL DZERO(WORK(KFDP),   N2PSHX*NZ)
         CALL DZERO(WORK(KEIGVCP),N2PSHX*NZ)
         CALL DZERO(WORK(KEIGP),  NPSHT)
!
!        extract p-p part of WORK(KFD)

         do indx_loop = 1, max_loop

           if(spinfr_krmc)then
!            get boson irrep
             is_symm = indx_loop - 1
           else if(opt_chckjz)then
!            get mj-value
             is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
           else
             is_symm = 0
           end if
           CALL MATGAT_symm(WORK(KFD),NORBT,NORBT,WORK(KFDP),NPSHT,
     &                      NPSHT,IDXP2G,work(kibeig),is_symm,NZ)
         end do

         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*) 'KRMCCNO: positronic FOCK matrix (extracted)'
           CALL PRQMAT(WORK(KFDP),NPSHT,NPSHT,
     &                 NPSHT,NPSHT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
!
!        diagonalize
         CALL DIAG_piv2fock('P',WORK(KFDP),WORK(KEIGP),WORK(KEIGVCP),
     &                      IDXP2G,work(kibeig),WORK(KFREE),LFREE)

!        insert new p-p matrix block into transformation matrix
         do indx_loop = 1, max_loop

           if(spinfr_krmc)then
!            get boson irrep
             is_symm = indx_loop - 1
           else if(opt_chckjz)then
!            get mj-value
             is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
           else
             is_symm = 0
           end if
           CALL MATSCT_symm(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCP),
     &                      NPSHT,NPSHT,IDXP2G,work(kibeig),is_symm,NZ)
         end do

         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*)
     &     'KRMCCNO: transformation matrix (p-p + i-i)'
           CALL PRQMAT(WORK(KEIGVC),NORBT,NORBT,
     &                 NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF

         CALL MEMREL('KRMCCNO.after p-p',WORK,1,KFDP,KFREE,LFREE)
      END IF
!
!     3.3 'active' secondary orbitals: v-v block
!     -------------
!
      IF ( NVSHT .gt. 0 ) THEN
!
!        ... allocate for v-v block diagonalization of Fock-matrix
         CALL MEMGET('REAL',KFDVI,   N2VSHXQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVCVI,N2VSHXQ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KEIGVI,  NVSHT  ,WORK,KFREE,LFREE)
!        ... initialize
         CALL DZERO(WORK(KFDVI),   N2VSHXQ)
         CALL DZERO(WORK(KEIGVCVI),N2VSHXQ)
         CALL DZERO(WORK(KEIGVI),  NVSHT  )
!
!        extract v-v part of WORK(KFD)

         do indx_loop = 1, max_loop

           if(spinfr_krmc)then
!            get boson irrep
             is_symm = indx_loop - 1
           else if(opt_chckjz)then
!            get mj-value
             is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
           else
             is_symm = 0
           end if
           CALL MATGAT_symm(WORK(KFD),NORBT,NORBT,WORK(KFDVI),NVSHT,
     &                      NVSHT,IDXV2G,work(kibeig),is_symm,NZ)
         end do

         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*) 'KRMCCNO: virtual FOCK matrix (extracted)'
           CALL PRQMAT(WORK(KFDVI),NVSHT,NVSHT,
     &                 NVSHT,NVSHT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
!
!        diagonalize
         CALL DIAG_piv2fock('V',WORK(KFDVI),WORK(KEIGVI),WORK(KEIGVCVI),
     &                      IDXV2G,work(kibeig),WORK(KFREE),LFREE)

         IF (FAC_FV .NE. 0.0D0) THEN
         IF (IPROPT .GE. 1) THEN
            WRITE (LUPRI,'(/A,1P,D10.2)')
     &      ' Modified virtual orbitals for MCSCF from FC + a*FV, a =',
     &      DKRMC_MVOFAC
            CALL PREIGV(WORK(KEIGVI),NVSHT,0,WORK(KFREE),LFREE)
         END IF
         END IF
!
!        insert new v-v matrix block into transformation matrix
         do indx_loop = 1, max_loop

           if(spinfr_krmc)then
!            get boson irrep
             is_symm = indx_loop - 1
           else if(opt_chckjz)then
!            get mj-value
             is_symm = ((-1)**(indx_loop+1))*(2*(indx_loop-1)+1)
           else
             is_symm = 0
           end if
           CALL MATSCT_symm(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCVI),
     &                      NVSHT,NVSHT,IDXV2G,work(kibeig),is_symm,NZ)
         end do

         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*) 
     &     'KRMCCNO: transformation matrix (p-p + i-i + v-v)'
           CALL PRQMAT(WORK(KEIGVC),NORBT,NORBT,
     &                 NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
!
         CALL MEMREL('KRMCCNO.after v-v',WORK,1,KFDVI,KFREE,LFREE)
      END IF
C
C     3.4 deleted secondary orbitals: d-d block 
C     -------------
C
      IF( NDSHT .gt. 0 )THEN
C        ... allocate
         CALL MEMGET('REAL',KEIGVCDL,N2DSHXQ,WORK,KFREE,LFREE)
C
C        ... initialize and make unit matrix
C
         CALL DUNIT2(WORK(KEIGVCDL),NDSHT,NDSHT,NDSHT,NZ)
C
C        insert d-d matrix block into complete transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCDL),
     &               NDSHT,NDSHT,IDXD2G,NZ)
         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*) 
     &     'KRMCCNO: transformation matrix (p-p + i-i + v-v + d-d)'
           CALL PRQMAT(WORK(KEIGVC),NORBT,NORBT,
     &                 NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         CALL MEMREL('KRMCCNO.after d-d',WORK,1,KEIGVCDL,KFREE,LFREE)
      END IF
C
C     3.5 unit matrix for a-a block
C     -----------------------------
C
      IF (NASHT .GT. 0) THEN
C        ... allocate
         CALL MEMGET('REAL',KEIGVCA,N2ASHXQ,WORK,KFREE,LFREE)
C
C        ... initialize and make unit matrix
C
         CALL DUNIT2(WORK(KEIGVCA),NASHT,NASHT,NASHT,NZ)
C
C        insert a-a matrix block into complete transformation matrix
C
         CALL MATSCT(WORK(KEIGVC),NORBT,NORBT,WORK(KEIGVCA),
     &               NASHT,NASHT,IDXU2G,NZ)
C
         IF (IPROPT .GE. 5) THEN
           WRITE(LUPRI,*) 
     &  'KRMCCNO: transformation matrix (p-p + i-i + v-v + d-d + a-a)'
           CALL PRQMAT(WORK(KEIGVC),NORBT,NORBT,
     &                 NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
         CALL MEMREL('KRMCCNO.after a-a',WORK,1,KEIGVCA,KFREE,LFREE)
      END IF
C
      CALL MEMREL('KRMCCNO.after FD diag',WORK,1,KFD,KFREE,LFREE)
C
C     4. Transform orbitals
C     ---------------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
      DO I = 1, NFSYM
C
C        CMO(new) = CMO(old) * EIGVC
C
         NBASI = NFBAS(I,0)
         NORBI = NORB(I)
C
         IF( IPROPT .ge. 2 )THEN
           WRITE(LUPRI,'(A,I3)') ' KRMCCNO: '//
     &     'Coefficients before transformation for symmetry',I
           CALL PRQMAT(CMO(1+ICMOQ(I)),NBASI,NORBI,
     &          NBASI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
           WRITE(LUPRI,'(A,I3)')
     &     ' KRMCCNO: transformation matrix for symmetry',I
           CALL PRQMAT(WORK(KEIGVC+I2ORBX(I,I)),NORBI,NORBI,
     &          NORBT,NORBT,NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Transform all orbitals
C
         CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KCMO+ICMOQ(I)),NBASI,NORBI,NZ,
     &        'N','N',IPQTOQ(1,0),
     &        WORK(KEIGVC+I2ORBX(I,I)),NORBT,NORBT,NZ,
     &        D0,IPQTOQ(1,0),
     &        CMO(1+ICMOQ(I)),NBASI,NORBI,NZ)
C
C        Copy back active orbitals
C
         ICMOA = ICMOQ(I) + (NPSH(I)+NISH(I)) * NBASI
         CALL DCOPY(NBASI*NASH(I),WORK(KCMO+ICMOA),1,
     &        CMO(1+ICMOA),1)
C
         IF( NDSH(I) .gt. 0 )THEN
C
C          Copy back deleted orbitals from s-s block
C
           ICMOD = ICMOQ(I) + 
     &             (NPSH(I)+NISH(I)+NASH(I)+NVSH(I)) * NBASI
           CALL DCOPY(NBASI*NDSH(I),WORK(KCMO+ICMOD),1,
     &          CMO(1+ICMOD),1)
         END IF
C
         IF( IPROPT .ge. 2 )THEN
           WRITE(LUPRI,'(A,I3)') ' KRMCCNO: '//
     &     'MO coefficients after FD transformation for symmetry',I
           CALL PRQMAT(CMO(1+ICMOQ(I)),NBASI,NORBI,
     &          NBASI,NORBI,NZ,IPQTOQ(1,0),LUPRI)
         END IF
      END DO
C
C     Deallocate everything but the density matrix.
C
      CALL MEMREL('KRMCCNO.after FD rot',WORK,1,KEIGVC,KFREE,LFREE)
C
 1000 CONTINUE
C
      IF (.NOT. (TRTONO .OR. TRTOCO) .OR. NZCONF.EQ.0) GOTO 2000
C
C     **********************************************************
C     *** Transform to natural or canonical active  orbitals ***
C     **********************************************************
C
C
      CALL QUIT('*** ERROR in KRMCCNO ***' //
     &            'TRTOCO/TRTONO not implemeted yet!')
C
C
      CALL MEMGET('REAL',KEIG,NASHT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGVC,N2ASHXQ,WORK,KFREE,LFREE)
C
C     Find transformation matrix
C     --------------------------
C
      IF (MCTYPE .EQ. JGAS) THEN
C
C        Calculate pseudo-natural or ``pseudo''-Fock type
C        transformation matrix.
C
C        CALL RGASNO(WORK(KDV),WORK(KEIGVC),WORK(KEIG))
C
         call quit('*** ERROR in KRMCCNO *** ' //
     &        'NO not implemented for GAS')
C
      ELSE
C
C        CAS wave function.
C
C        Find transformation matrix.
C        (WORK(KDV) is either active part of FD or DV).
C
         call quit('*** ERROR in KRMCCNO *** ' //
     &        'NO not implemented for CAS yet')
C
C        hjaaj TODO: at least NORBT must be NASHT in RSJACO/QDIAG calls
C        check for other errors ... /feb 2004
         IF (NZ .EQ. 1) THEN
            IJOB = 1
            IORDER = 0
            IPACK = 0
            CALL RSJACO(NORBT,NORBT,NORBT,WORK(KDV),WORK(KEIG),
     &           IJOB,IORDER,IPACK,WORK(KEIGVC))
         ELSE
            MATZ = 1
            CALL QDIAG(NZ,NORBT,WORK(KDV),NORBT,NORBT,
     &           WORK(KEIG),MATZ,WORK(KEIGVC),NORBT,NORBT,
     &           WORK(KFREE),LFREE,IERR)
            IF (IERR .NE. 0) THEN
               WRITE(LUPRI,'(/,1X,2A,I4)')
     &              '*** ERROR in RGETNO ***: ',
     &              'QDIAG failed with error code ',IERR
               CALL QUIT('*** ERROR in RGETNO ***')
            END IF
         END IF
C
      END IF
C
C     Transform orbitals
C     ------------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL DCOPY(NCMOTQ,CMO,1,WORK(KCMO),1)
C
      DO I = 1, NFSYM
         IF (NASH(I) .NE. 0) THEN
C
C           CMO(NO) = CMO * UNO
C
            NBASI = NFBAS(I,0)
            NASHI = NASH(I)
            NORBI = NORB(I)
            ICMOA = 1 + ICMOQ(I) + (NPSH(I)+NISH(I)) * NBASI
            CALL QGEMM(NBASI,NASHI,NASHI,D1,
     &           'N','N',IPQTOQ(1,0),
     &           WORK(KCMO+ICMOA-1),NBASI,NORBI,NZ,
     &           'N','N',IPQTOQ(1,0),
     &           WORK(KEIGVC+I2ASHX(I,I)),NASHT,NASHT,NZ,
     &           D0,IPQTOQ(1,0),
     &           CMO(ICMOA),NBASI,NORBI,NZ)
C
         END IF
      END DO
C
C     Counter-rotate CI vector
C     ------------------------
C
C     CALL RTRACI(CREF,WORK(KEIGVC))
C
 2000 CONTINUE
      CALL QEXIT('KRMCCNO')
      ipropt = ipropt_save
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck RPROP */
      SUBROUTINE KRCI_PROP(CMO,WORK,LWORK)
C***********************************************************************
C
C     1. calculate property MO representation for list of (symmetric) 
C        one-electron (and later also two-electron) operators.
C
C     2. calculate properties , e.g., transition dipole moments. 
C
C     3. obtain final quantities , e.g., oscillator strength f. 
C
C     Written by S. Knecht - Aug 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "consts.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbopt.h"
#include "krciprop.h"
#include "dcbkrci.h"
C
      DIMENSION CMO(*), WORK(*)
      integer, allocatable :: iomval(:)
      logical selected
C
#include "memint.h"
      MZ = MIN(2,NZ)
C
C     introduce property task
C
      CALL TITLER('KR-CI property calculation ','*',199)
C
      IF (OPT_CIPROGRAM .ne. 'LUCIAREL') THEN
         WRITE(LUPRI,'(/A/A/2A)')
     &   'WARNING: Requested KR-CI property calculations skipped',
     &   'WARNING: because they are only implemented for LUCIAREL',
     &   'WARNING: and you have requested ',OPT_CIPROGRAM
         RETURN
      END IF

      WRITE(LUPRI,9000)
 9000 FORMAT(
     &     3X,'This is output from the DIRAC KR-CI property module',/,
     &     //,
     &     3X,'Written by: ',/,
     &     3X,'  Stefan Knecht and Hans Joergen Aa. Jensen ',//,
     &     3X,'With contributions from: ',/,
     &     3X,'  Malaya K. Nayak and Timo Fleig ')
      CALL TITLER(
     &'(P-odd, P,T-odd properties and magnetic hyperfine interaction)',
     & ' ',199)
      WRITE(LUPRI,9005)
 9005 FORMAT(
     &     1X,79('*'))
!
      WRITE(LUPRI,9006)
 9006 FORMAT(/3X,'Proper citations:',//,
     & 4X,'S.R. Knecht, PhD thesis, University of Duesseldorf, 2009.',/,
C    & 4X,'    Parallel Relativistic Multiconfiguration Methods:',/,
C    & 8X,'New Powerful Tools for Heavy-Element
C    & Electronic-Structure Studies,',
     & 4X,'T. Fleig, M.K. Nayak, Phys. Rev. A  88 (2013) 032514',/,
     & 4X,'T. Fleig, M.K. Nayak, J. Mol. Spectrosc.  300 (2014) 16',/,
     & 4X,'M. Denis et al., New J. Phys.  7 (2015) 043005',/,
     & 4X,'T. Fleig et al., Phys. Rev. A  93 (2016) 012505',/,
     &     //,
     & 1X,79('*'),/)

C
      WRITE (LUPRI,'(/A/A)')
     &   '  Symmetry   # of electronic states',
     &   '  ________   ______________________'
C               1              4
      NPROP_ROOTS_KRCI = 0
      DO I = 1, NKRCI_MAX_SYM
         NPROP_ROOTS_KRCI = NPROP_ROOTS_KRCI + NKRCI_CIROOTS(I)
         if (NPROP_ROOTS_KRCI.ne.0) then
           WRITE(LUPRI,'(I7,I15)') I, NKRCI_CIROOTS(I)
         end if
      END DO
      WRITE(LUPRI,'(A/A,I15/)')
     &   '  _________________________________',
     &   '  Total', NPROP_ROOTS_KRCI
C
C     eigenstate <--> symmetry relation
      IXROOTS     = 0
      IXROOTS_SYM = 0
      CALL IZERO(ISYMEIG_KRCI,NPROP_ROOTS_KRCI)
      DO I = 1, NKRCI_MAX_SYM
         IXROOTS_SYM = I
         DO J = 1, NKRCI_CIROOTS(I)
            IXROOTS = IXROOTS + 1
            ISYMEIG_KRCI(IXROOTS) = IXROOTS_SYM
         END DO
      END DO
C
C     KRCI_CVECS.x vector file extensions
      CALL DETXFLAB(XSYMFLAB,NKRCI_MAX_SYM)
C
C     ***************************************************
C     *** Memory allocation *****************************
C     ***************************************************
C
      LXPRPKRCI = MZ * NPROP_KRCI * NPROP_ROOTS_KRCI**2
      CALL MEMGET('REAL',KXPRPKRCI,LXPRPKRCI,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPRPCORE, NPROP_KRCI,WORK,KFREE,LFREE)
C
C     ... initialize
      CALL DZERO(WORK(KXPRPKRCI),LXPRPKRCI)
      CALL DZERO(WORK(KPRPCORE),NPROP_KRCI)
C
C     ***************************************************
C     *** operator(s) in MO representation **************
C     ***************************************************
C
C     active-active part in MOLFDIR format saved on 
C     file KRMC_FOCK
C
      CALL KRCIPRPMO(CMO,WORK(KPRPCORE),WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** calculate properties **************************
C     ***************************************************
C
      CALL KRCIGETPRP(WORK(KXPRPKRCI),WORK(KFREE),LFREE)
C
C     ***************************************************
C     *** get spectroscopic quantities ******************
C     ***************************************************
C
C     array of eigen- and omega values 
      CALL MEMGET('REAL',KXROOTS,NPROP_ROOTS_KRCI,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KXROOTS),NPROP_ROOTS_KRCI)
      allocate(iomval(nprop_roots_krci))
C
      CALL KRCIPRPANA(WORK(KXPRPKRCI),WORK(KXROOTS),WORK(KPRPCORE),
     &                iomval,WORK(KFREE),LFREE,MZ)
      CALL MEMREL('RPROP_KRCI',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     CI root selection based on the omega value
      if(state_select_krci)then
        iroot_cnt = 0
        do i = 1, nkrci_max_sym
          if(nkrci_selom(i).gt.0)then
            selected = .false.
!           do j = 1, nkrci_ciroots(i)
             iom_val_tmp = iabs(iomval(iroot_cnt+NOMEGASEL_statenr))
             if((iom_val_tmp.eq.nkrci_selom(i)).and.(.not.selected))then
                call krci_cp_sel_root(xsymflab(i),NOMEGASEL_statenr,
     &                                mz,iom_val_tmp)
                selected = .true.
             end if
!           end do
          end if
          iroot_cnt = iroot_cnt + nkrci_ciroots(i)
        end do
      end if
      deallocate(iomval)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck KRCIPRPMO */
      SUBROUTINE KRCIPRPMO(CMO,PRPCOREE,WORK,LWORK)
C***********************************************************************
C
C     Transform property integrals to MO basis.
C     Integrals are saved on file KRMC_FOCK.
C
C     Input:
C        CMO     - orbital coefficients.
C
C     Written by S. Knecht - Aug 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C write-protection for DUMMY.
#include "dummy.h"
#include "consts.h"
#include "maxorb.h"
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "krciprop.h"
#include "dcbxpr.h"
      DIMENSION CMO(*), PRPCOREE(*), WORK(*)
      CHARACTER*8 FILELAB
      LOGICAL LOPEN
#if defined MCSCF_DEBUG
      real*8, allocatable:: jz_momat(:)
#endif
C
#include "memint.h"
C
      IPROPT = 00
      MZ = MIN(2,NZ)
C
      CALL QENTER('KRCIPRPMO')
C
      INQUIRE(FILE='AOPROPER',OPENED=LOPEN)
      IF (.NOT. LOPEN) THEN
         OPEN (LU1INT,STATUS='OLD',FORM='UNFORMATTED',
     &        FILE='AOPROPER')
      ENDIF
#if defined MCSCF_DEBUG
      if( dojzexp) then
        allocate(jz_momat(n2ashxq))
        jz_momat = 0
        iprintjz = 0
      end if
#endif

C
C     ***************************************************
C     *** Get property matrices *************************
C     ***************************************************
C
      DO I = 1, NPROP_KRCI
C
         INDXPR  = LPROP_KRCI(I)
         ISYM    = IPRPSYM(INDXPR)
         ITRSYM  = IPRPTIM(INDXPR)
         IOPSY   = JBTOF(ISYM-1,1)
         FILELAB = PRPNAM(INDXPR)(1:8)
C
         CALL MEMGET('REAL',KRPRPMO,N2ORBXQ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KRPRPMO),N2ORBXQ)
         CALL PRPMAT(INDXPR,IOPSY,WORK(KRPRPMO),.TRUE.,WORK,CMO,
     &               DUMMY,ICMOQ,NORB,WORK,KFREE,LFREE,IPROPT)

         if(filelab(3:8) == 'dipole')then
!          scale with proper prefactor for dipole (length) integrals (see factor -1.0d0 in src/prp/pamprp.F - def_dipole)
           call dscal(NZ*NORBT**2,-1.0d0,WORK(KRPRPMO),1)
         end if
C
C        print full property matrix
C
         IF( IPROPT .ge. 5 ) THEN
            WRITE(LUPRI,'(/2A)')'   KRCIPRPMO: MO matrix of property: ',
     &                              FILELAB
            CALL PRQMAT(WORK(KRPRPMO),NORBT,NORBT,NORBT,
     &                  NORBT,NZ,IPQTOQ(1,ISYM-1),LUPRI )
         ENDIF
C
C        calculate PRPCORE(I) = 2.0D0 * SUM_(all inactive) P_ii
C
         IF( NISHT .gt. 0)THEN
            TEMP_PCORE  = SUM_DIAG_ELM('I',WORK(KRPRPMO),NORBT,NORBT)
            PRPCOREE(I) = D2 * TEMP_PCORE
         ELSE
            PRPCOREE(I) = D0
         END IF
         IF (IPROPT .GE. 5)
     &      WRITE(LUPRI,9000) FILELAB,PRPCOREE(I)
 9000 FORMAT(' (KRCIPRPMO): ',A,' property core value: ',F25.15/)
C
C        get active-active part
C
         CALL MEMGET('REAL',KRPRPMOAC,   N2ASHXQ,WORK,KFREE,LFREE)
CSK      IPROPT = 50
         CALL RGETAC(WORK(KRPRPMO),WORK(KRPRPMOAC),IPROPT)
C
C        transform to MOLFDIR format
         CALL MEMGET('REAL',KRPRPMOACM,4*N2ASHX*MZ,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KRPRPMOACM),  4*N2ASHX*MZ)
C        WRITE(LUPRI,'(/2A)')'   KRCIPRPMO: MO matrix of property: ',
C    &                           FILELAB
#if defined MCSCF_DEBUG
         if( dojzexp )then
           xmultfac = 1.0D0
           if( filelab .eq. 'z-angula') then
             xmultfac = -1.0D0
             call daxpy(n2ashxq,xmultfac,work(krprpmoac),1,jz_momat,1)
             iprintjz = iprintjz + 1
           else if ( filelab .eq. 'z-spin')then
             call daxpy(n2ashxq,xmultfac,work(krprpmoac),1,jz_momat,1)
             iprintjz = iprintjz + 2
           end if
           if( iprintjz .eq. 3) then
             write(lupri,*) ' property matrix (a-a part) of jz-operator'
             call prqmat(jz_momat,nasht,nasht,nasht,
     &                   nasht,nz,ipqtoq(1,isym-1),lupri )
           end if
         end if
         WRITE(LUPRI,'(A,3I5)')
     &   ' - KRCIPRPMO: Boson, Fermion, and Time reversal symmetry :',
     &   ISYM,IOPSY,ITRSYM
#endif
         CALL QFC2MFC(WORK(KRPRPMOAC),WORK(KRPRPMOACM),ISYM,
     &                ITRSYM,IPROPT)
C             QFC2MFC(DQFC,DMFC,ISYM,ITIM,IPRINT)
C
C        save on file KRMC_FOCK
C
         IF (FILELAB(4:6).EQ.'HYP'.OR.FILELAB(4:6).EQ.'ANA')THEN
           CALL FLIP_BLOCK(WORK(KRPRPMOACM),FILELAB,IPROPT)
         END IF
C
         CALL KRCI_PRPFILE(LUKRM3,FILELAB,WORK(KRPRPMOACM),
     &                     4*N2ASHX*MZ,0)
         CALL MEMREL('KRCIPRPMO',WORK,KWORK,KWORK,KFREE,LFREE)
      END DO
C
#if defined MCSCF_DEBUG
      if( dojzexp) then
        deallocate(jz_momat)
      end if
#endif
C
!     IPROPT = 00
      CLOSE(LU1INT,STATUS = 'KEEP')
      CALL QEXIT('KRCIPRPMO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck KRCIPRPANA */
      SUBROUTINE KRCIPRPANA(XPROPMAT,EVECS,EPROP_CORE,IOMVAL,
     &                      WORK,LWORK,MZ)
C***********************************************************************
C
C     Calculate requested spectroscopic quantities.
C
C     implemented: 1. oscillator strength for electronic transitions
C                  2. omega expectation values (for j_z = s_z + l_z)
C
C     Input: property matrix XPROPMAT
C
C     Written by S. Knecht - Nov 2008
C
C         T. Fleig - Mar 2016
C         Introduced .EEDM, .MHYP, .ENSPS flags programmed by MKN and TF
C
C         Malaya K. Nayak - September 2018
C         Introduced .NMQM flag programmed by MKN
C         Introduced a general property module for .OPERATOR technique
C
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
      PARAMETER ( XTHRESH = 1.0D-7 )
      PARAMETER ( DEGEN   = 1.0D-10 )
C
#include "dcborb.h"
#include "dcbopt.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "krciprop.h"
#include "dcbkrci.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "mxcent.h"
#include "dipole.h"
C#include "nuclei.h"
C#include "chrnos.h"
      DOUBLE COMPLEX XMOM
      REAL*8 OSCL, OSCLABS
      DIMENSION XPROPMAT(NPROP_ROOTS_KRCI,NPROP_ROOTS_KRCI,
     &                   NPROP_KRCI,2)
CMKN     &                   NPROP_KRCI,MZ)
      DIMENSION GTENSOR(NPROP_ROOTS_KRCI),GTENSORNORM(NPROP_ROOTS_KRCI),
     &          ENSPSMOM(NPROP_ROOTS_KRCI),HYFENFIN(NPROP_ROOTS_KRCI),
     &          HYFEDMFIN(NPROP_ROOTS_KRCI),EEFFAU(NPROP_ROOTS_KRCI)
      DIMENSION WORK(*), OSCL(3), EVECS(*), EPROP_CORE(*), OSCLABS(3)
      INTEGER IND_DIP(3), IND_PRP(12)
      CHARACTER XLABEL*16,OPTYPE*12
      LOGICAL PRINT_EV
!     dynamic memory allocation
      integer, allocatable :: indx(:)
      integer, intent(in)  :: mz
C
#include "memint.h"
C
      CALL QENTER('XKRCIPRPANA')
C
      IZERO = 0
      ZERO = 0.0d0
      IONE = 1
      ONE = 1.0d0
C
CTF   Initialize
      do IL = 1,NPROP_ROOTS_KRCI,1
        GTENSOR(IL) = ZERO
        GTENSORNORM(IL) = ZERO
        EEFFAU(IL) = ZERO
      end do
C
C     print the list of states if available otherwise skip it.
      PRINT_EV = .FALSE.
C
C     search and possibly get energies from file KRCI_CVECS.INFO
      CALL GET_KRCICVECS_E(EVECS,NKRCI_CIROOTS,XSYMFLAB,
     &                     NKRCI_MAX_SYM,PRINT_EV)
C
      IF( PRINT_EV )THEN
C
C       print list of eigenstates and energies
        WRITE(LUPRI,'(/,1X,A,F14.10,A,F17.9,A)')
     &   '( 1 au =',XTEV,' eV  / ',XTKAYS,' cm-1)'

        WRITE(LUPRI,'(/A/A)') 
     &   '  List of eigenstates in property calculation',
     &   '  ___________________________________________________'
            WRITE(LUPRI,'(A)') 
     &   '   # eigenstate    irrep           energy (a.u.)'
        DO IVEC=1, NPROP_ROOTS_KRCI
            WRITE(LUPRI,'(I11,7X,A4,5X,1F22.10)') 
     &            IVEC,XREPEIG(ISYMEIG_KRCI(IVEC)),EVECS(IVEC)
        END DO
        WRITE(LUPRI,'(A/)') 
     &   '  ___________________________________________________'
C
C       sort
        call alloc(INDX,NPROP_ROOTS_KRCI)
        INDX = 0
        CALL INDEXX(NPROP_ROOTS_KRCI,EVECS,INDX)
C
C       relative ordering
        WRITE(LUPRI,'(/A/A)') 
     &   '  List of eigenstates in property calculation (sorted)',
     & '  ___________________________________________________________'//
     &   '_________________________________________________________'
        WRITE (LUPRI,'(A)')
     &           '  Level | abs. eigenvalue     |'//
     &           ' rel. eigenvalue (a.u.) |'//
     &           ' rel. eigenvalue (eV) | rel. eigenvalue (cm-1) |'//
     &           ' irrep '
!
        ILEVEL=0
        EWF=EVECS(INDX(1))
        DO I=1,NPROP_ROOTS_KRCI
          E = EVECS(INDX(I))
          WRITE(LUPRI,1010) ILEVEL,E,E-EWF,XTEV*(E-EWF),XTKAYS*(E-EWF),
     &                      XREPEIG(ISYMEIG_KRCI(INDX(I)))
          ILEVEL=ILEVEL+1
        ENDDO
!
        WRITE(LUPRI,'(A/)') 
     & '  ___________________________________________________________'//
     &   '_____________________________________________'
C
C
      END IF
C
C     *****************************************************************
C     ************************** PROPERTIES ***************************
C     *****************************************************************
C     Thsi is a General Property Module Implemented for an arbitrary
C     one-electron operator (defined in .OPERATOR list).
C
C     Written by Malaya K. Nayak - September 28, 2018
C     *****************************************************************
C     ****************** GENERAL PROPERTY MODULE **********************
C     *****************************************************************
      if (dogenp_krci) then
C
C       print title
        write(lupri,*)
        print*,' NPROP_ROOTS_KRCI = ',NPROP_ROOTS_KRCI
        write(lupri,'(/28x,a/28x,a/28x,a)')
     & ' ********************************************************',
     & ' *********     General Property Module     **************',
     & ' ********************************************************'
C
      DO I = 1, NPROP_KRCI
C
         INDXPR  = LPROP_KRCI(I)
         XLABEL = PRPNAM(INDXPR)(1:8)
         IND_PRP(I) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
C
C       ----------------------------
C       print Property expectation value
C       ----------------------------
        write(lupri,'(/a,A6,a/a)')
     & '     < state |      ',XLABEL,'      | state >   '//
     & '    Property Expectation [a.u.]      Norm [a.u.]'
        write(lupri,'(a)')
     & '  __________________       _________________'//
     & '    ___________________________________________________'//
     & '_______'
C
       DO IVEC = 1,NPROP_ROOTS_KRCI
          JVEC  = IVEC
          INIST = IVEC
          JFLST = JVEC
          ISYMST_IN = ISYMEIG_KRCI(INIST)
          ISYMST_FL = ISYMEIG_KRCI(JFLST)
          EDMMOM = DCMPLX(XPROPMAT(INIST,JFLST,IND_PRP(I),1),
     &                    XPROPMAT(INIST,JFLST,IND_PRP(I),2))
          EDMNORM = SQRT(EDMMOM**2)
C
          write(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F20.12,5X,1F20.12)')
     & '     * state',INIST,'(',XREPEIG(ISYMST_IN),')       * state',
     &      JFLST,'(',XREPEIG(ISYMST_FL),')   ',
     &      EDMMOM,EDMNORM
       END DO
        write(lupri,'(a)')
     & '  __________________________________________'//
     & '    ___________________________________________________'//
     & '_______'
       END DO
C
      end if
C
C     ******************************************************************
C
C     1. (transition dipole moments) and oscillator strengths
C     =======================================================
C
      IF( DOOSCILLST )THEN
C
!       WRITE(lupri,*) ' content of XPROPMAT'
!       CALL WRTMATMN(XPROPMAT,1,NPROP_ROOTS_KRCI**2*Nprop_krci*mz,1,
!    &       NPROP_ROOTS_KRCI**2*nprop_krci*mz,lupri)

C       write title
        WRITE(LUPRI,'(/28X,A/28X,A/28X,A)') 
     & ' ********************************************************',
     & ' ********* printing transition dipole moment(s) *********',
     & ' ********************************************************'
C
C       calculate and print TDM - length gauge
C       ------------------------------------------------------
C
C       determine index for (x,y,z)
        XLABEL     = 'X dipole length'
        IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Y dipole length'
        IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Z dipole length'
        IND_DIP(3) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XREAL = 0.0D0
        XCOMP = 0.0D0
C
        WRITE(LUPRI,'(/,1X,A,F14.10,A)')
     &   ' ( 1 au =',debye,' Debye) TDM contributions with a norm'//
     &   ' smaller than 1.0d-07 Debye are discarded (not printed).'

C       ---------
C       print TDM
C       ---------
        WRITE(LUPRI,'(/A/A)')
     & '   <initial state   |x,y,z|   final state>  '//
     & '            TDM [a.u.] (x, y, z)           Norm [a.u.]'//
     & '  Norm [Debye]',
     & '  __________________       _________________'//
     & '    ________________________________________________________'//
     & '_______'
        DO 201 IVEC=1,NPROP_ROOTS_KRCI
         DO 101 JVEC= IVEC+1, NPROP_ROOTS_KRCI
            DEIG  = EVECS(IVEC) - EVECS(JVEC)
            IF( DEIG < ZERO )THEN
              INIST = IVEC
              JFLST = JVEC
            ELSE
              INIST = JVEC
              JFLST = IVEC
            END IF
            DO IPRP=1,3
C             length gauge
              IF( MZ .eq. 2 )
     &        XCOMP         = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),2)
              XREAL         = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),1)
              XMOM          = XREAL
              if(mz.eq.2)then
                xmom        = dcmplx(xreal,xcomp)
              end if
              OSCL(IPRP)    = XMOM
              OSCLABS(IPRP) = ABS(XMOM)
            END DO
            ISYMST_IN = ISYMEIG_KRCI(INIST)
            ISYMST_FL = ISYMEIG_KRCI(JFLST)

!           compute relative state number within a given symmetry irrep
            inist_rel = inist
            jflst_rel = jflst
            do i = 1, isymst_in-1
             inist_rel = inist_rel - NKRCI_CIROOTS(i)
            end do
            do i = 1, isymst_fl-1
             jflst_rel = jflst_rel - NKRCI_CIROOTS(i)
            end do
C
C           print only values above a given threshold - 1.0D-7
            IF( SQRT((OSCL(1)**2)+(OSCL(2)**2)+(OSCL(3)**2))*debye 
     &         < XTHRESH ) GOTO 101

       WRITE(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,3F12.8,1F14.8,1F14.8)')
     &  '   * state',INIST_rel,'(',XREPEIG(ISYMST_IN),
     &      ')         * state',
     &      JFLST_rel,'(',XREPEIG(ISYMST_FL),')   ',
     &      OSCL(1),OSCL(2),OSCL(3),
     &      SQRT((OSCL(1)**2)+(OSCL(2)**2)+(OSCL(3)**2)),
     &      SQRT((OSCL(1)**2)+(OSCL(2)**2)+(OSCL(3)**2))*debye
 101      CONTINUE
 201    CONTINUE
        WRITE(LUPRI,'(A)')
     & '  __________________________________________'//
     & '    ________________________________________________________'//
     & '_______'
C
#ifdef CODE_NOT_WORKING
C
C       determine index for (dx,dy,dz)
        XLABEL     = 'X dip vel'
        IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Y dip vel'
        IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Z dip vel'
        IND_DIP(3) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XREAL = 0.0D0
        XCOMP = 0.0D0
C
C       Check sum rule for oscillator strength
        SUMV = 0.0D0
        WRITE(LUPRI,'(/A/A)')
     & '    initial state             final state   '//
     & '      oscillator strength (velocity gauge)',
     & '  __________________       _________________'//
     & '      ____________________________________'
        DO 400 IVEC=1,NPROP_ROOTS_KRCI
          DO 300 JVEC= IVEC+1, NPROP_ROOTS_KRCI
            IF( EVECS(IVEC) .lt. EVECS(JVEC) )THEN
              INIST = IVEC
              JFLST = JVEC
            ELSE IF( EVECS(IVEC) .eq. EVECS(JVEC) )THEN
              GOTO 300
            ELSE IF( EVECS(IVEC) .gt. EVECS(JVEC) )THEN
              INIST = JVEC
              JFLST = IVEC
            END IF
            DEIG  = EVECS(JFLST) - EVECS(INIST)
            DO IPRP=1,3
C             velocity gauge
              XREAL = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),1)
              IF( MZ .eq. 2 ) 
     &        XCOMP = XPROPMAT(INIST,JFLST,IND_DIP(IPRP),2)
              XMOM = DCMPLX(XREAL,XCOMP)
C
C             f_(dx,dy,dz)  = 2.0D0*|<inist|(x,y,z)|jflst>|**2 * (E_fl - E_ini)
              OSCL(IPRP) = 2.0D0*ABS(XMOM)**2*(DEIG)
            END DO
            SUMV = SUMV + (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0
            ISYMST_IN = ISYMEIG_KRCI(INIST)
            ISYMST_FL = ISYMEIG_KRCI(JFLST)
        WRITE(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F10.6)')
     & '   * state',INIST,'(',XREPEIG(ISYMST_IN),')         * state',
     &       JFLST,'(',XREPEIG(ISYMST_FL),')                    ',
     &       (OSCL(1)+OSCL(2)+OSCL(3))/3.0D0 
        WRITE(LUPRI,'(A,1F16.8,A,1F16.8)')
     & '  E(au)',EVECS(INIST),'    E(au)',EVECS(JFLST)
 300      CONTINUE
 400    CONTINUE
        WRITE(LUPRI,*) 
     & ' -----------------------------------------------------------'//
     & '------------'
        WRITE(LUPRI,'(1A,1F10.6/)')
     &     '   Sum of oscillator strength : ',SUMV
        WRITE(LUPRI,'(1A,1F10.6/)')
     &     '   Total sum of oscillator strengths (length + velocity'//
     &     ' gauge)  : ',SUMV+SUML
#endif
      END IF

C     2. permanent dipole moments
C     =======================================================
      if(dodipmom_krci)then
C
C       print title
        write(lupri,'(/4x,a/4x,a/4x,a)') 
     & ' ********************************************************',
     & ' ********* permanent electric  dipole moment(s) *********',
     & ' ********************************************************'
C
C       calculate and print dipole moment - length gauge
C       ------------------------------------------------------
C
C       determine index for (x,y,z)
        XLABEL     = 'X dipole length'
        IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Y dipole length'
        IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XLABEL     = 'Z dipole length'
        IND_DIP(3) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
        XREAL = 0.0D0
        XCOMP = 0.0D0

C       Nuclear contribution to dipole moment
C       ========================================
C
        call dipnuc(0,.false.)

!
!       -------------------
!       print dipole moment
!       -------------------
!       write(lupri,'(/a/a)')
!    & '            < state |x,y,z| state>  '//
!    & '               dipole moment [a.u.] (x, y, z)     Norm [a.u.]'//
!    & '  Norm [Debye]',
!    & '  __________________       _________________'//
!    & '    ________________________________________________________'//
!    & '_______'
        do 202 ivec = 1, nprop_roots_krci

          jvec  = ivec
          inist = ivec
          jflst = jvec
          do iprp = 1, 3
!           length gauge
            if(mz.eq.2)
     &      xcomp         = xpropmat(inist,jflst,ind_dip(iprp),2)
            xreal         = xpropmat(inist,jflst,ind_dip(iprp),1)
            xmom          = xreal
            if(mz.eq.2)then
              xmom          = dcmplx(xreal,xcomp)
            end if
            oscl(iprp)    = xmom
            osclabs(iprp) = abs(xmom)
          end do
          isymst_in = isymeig_krci(inist)
          isymst_fl = isymeig_krci(jflst)
!
!         print only values above a given threshold - 1.0D-7
          if((osclabs(1)+osclabs(2)+osclabs(3)).lt.xthresh)then
!           print *, 'dipole moment too small for state #',inist
            call dzero(oscl,3)
          end if


!      write(lupri,'(a,i3,a,a4,a,i3,a,a4,a,3f12.8,1f14.8,1f14.8/)')
!    &  '   * state',inist,'(',xrepeig(isymst_in),')         * state',
!    &      jflst,'(',xrepeig(isymst_fl),')   ',
!    &      oscl(1),oscl(2),oscl(3),
!    &      sqrt((oscl(1)**2)+(oscl(2)**2)+(oscl(3)**2)),
!    &      sqrt((oscl(1)**2)+(oscl(2)**2)+(oscl(3)**2))*debye

       write(lupri,'(/a,i3,a,a4,a/)')
     &  '     * dipole moment for electronic state: ',
     &           inist,' (',xrepeig(isymst_in),')'

      CALL PRSYMB(LUPRI,'-',76,5)
      WRITE(LUPRI,'(10X,3(A12,13X))')
     &      'Electronic  ','Nuclear     ','Total       '
      WRITE(LUPRI,'(10X,3(A12,13X))')
     &      'contribution','contribution','contribution'
      CALL PRSYMB(LUPRI,'-',76,5)
      IOFF = ICHAR('w')
      ! print in debye
      DO I = 1,3
         WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &        CHAR(IOFF+I),(DEBYE*(oscl(I)+eprop_core(ind_dip(i)))),
     &                                    ' Debye',
     &                   (DEBYE*DIPMN(I)),' Debye',
     &        (DEBYE*(oscl(I)+eprop_core(ind_dip(i))+DIPMN(I))),' Debye'
      ENDDO
      CALL PRSYMB(LUPRI,'-',76,5)
      ! print in a.u.
      DO I = 1,3
         WRITE(LUPRI,'(A1,3(3X,F16.8,A))')
     &        CHAR(IOFF+I),(oscl(I)+eprop_core(ind_dip(i))),' a.u. ',
     &                   (DIPMN(I)),' a.u. ',
     &        (oscl(I)+eprop_core(ind_dip(i))+DIPMN(I)),' a.u.'
      ENDDO


 202    continue

      end if
!
!     3. omega expectation values
!     ===========================
!
      IF( DOJZEXP )THEN
C       ... allocate omega matrix
        CALL MEMGET('REAL',KOMEGAEXP,NZ*NPROP_ROOTS_KRCI**2,WORK,
     &              KFREE,LFREE)
C       s_z operator
        CALL MEMGET('REAL',Ksz_op,NZ*NPROP_ROOTS_KRCI**2,WORK,
     &              KFREE,LFREE)
C       l_z operator
        CALL MEMGET('REAL',Klz_op,NZ*NPROP_ROOTS_KRCI**2,WORK,
     &              KFREE,LFREE)
C       ... initialize
        CALL DZERO(WORK(KOMEGAEXP),NZ*NPROP_ROOTS_KRCI**2)
        CALL DZERO(WORK(Ksz_op)   ,NZ*NPROP_ROOTS_KRCI**2)
        CALL DZERO(WORK(Klz_op)   ,NZ*NPROP_ROOTS_KRCI**2)
C
C       get and diagonalize omega expectation value matrix
        CALL GET_OMEGA_EXP(WORK(KOMEGAEXP),WORK(Ksz_op),WORK(Klz_op),
     &                     NPROP_ROOTS_KRCI,NZ,
     &                     XPROPMAT,NPROP_KRCI,XREPEIG,
     &                     MAX_NKRCI_MAX_SYM,NKRCI_CIROOTS,
     &                     naelec,IOMVAL,WORK(KFREE),LFREE)
C       release memory        
        CALL MEMREL('XKRCIPRPANA',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
!
C
C     =======================================================
C     4. electron electric dipole moment
C        (effective electric field)
C     =======================================================
      if(doeedm_krci)then
C
C       print title
        write(lupri,*)
        write(lupri,'(/28x,a/28x,a/28x,a)')
     & ' ********************************************************',
     & ' ********* electron Electric Dipole Moment **************',
     & ' ********************************************************'
C
         XLABEL = 'P2-EDM'
         IND_EDM = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
C
C       ----------------------------
C       print eEDM expectation value
C       ----------------------------
        write(lupri,'(/a,A6,a/a)')
     & '  <state| ',XLABEL,'  2*ic gamma0 gamma5 p^2 | state >  '//
     & '  E_eff [a.u.]             E_eff [GV/cm]   '
        write(lupri,'(a)')
     & '  __________________       _________________'//
     & '    ___________________________________________________'//
     & '_______'
        do IVEC = 1,NPROP_ROOTS_KRCI,1
          JVEC  = IVEC
          INIST = IVEC
          JFLST = JVEC
          ISYMST_IN = ISYMEIG_KRCI(INIST)
          ISYMST_FL = ISYMEIG_KRCI(JFLST)
          EEFFAU(IVEC) = DCMPLX(XPROPMAT(INIST,JFLST,IND_EDM,1),
     &                    XPROPMAT(INIST,JFLST,IND_EDM,2))
CTF          EDMNORM = sqrt(EDMMOM**2)
          EEFFGVCM = ( EEFFAU(IVEC) * EFAUMKSA )
C
          write(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F20.12,5X,1F20.12)')
     & ' EEDM  state',INIST,'(',XREPEIG(ISYMST_IN),')         state',
     &      JFLST,'(',XREPEIG(ISYMST_FL),')   ',
     &      EEFFAU(IVEC),EEFFGVCM
        end do
        write(lupri,'(a)')
     & '  __________________________________________'//
     & '    ___________________________________________________'//
     & '_______'
C
      end if
C
C     =======================================================
C     5. Magnetic Hyperfine Structure Constants
C     Written by Malaya K. Nayak, Sept 13, 2013
C            and Timo Fleig, April 7-8, 2016
C     =======================================================
      if (domhyp_krci) then
C
C       print title
        write(lupri,*)
        write(lupri,'(/28x,a/28x,a/28x,a)') 
     & ' ********************************************************',
     & ' *********  Hyperfine Structure Constants  **************',
     & ' ********************************************************'
C
      DO I = 1, NPROP_KRCI
         INDXPR  = LPROP_KRCI(I)
         XLABEL = PRPNAM(INDXPR)(1:8)
         IF (XLABEL(4:6).EQ.'HYP') THEN
         IND_PRP(I) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
C
C       ----------------------------
C       print mHYP expectation value
C       ----------------------------
        write(lupri,'(/a,A6,a/a)')
     & '    < state | ',XLABEL,' [alpha x E] | state >  '//
     & '     Hyperfine Const. [a.u.]     Delta E_HypFin [a.u.]   '
        write(lupri,'(a)')
     & '  __________________       _________________'//
     & '    ___________________________________________________'//
     & '_______'
C
       DO IVEC = 1,NPROP_ROOTS_KRCI
          JVEC  = IVEC
          INIST = IVEC
          JFLST = JVEC
          ISYMST_IN = ISYMEIG_KRCI(INIST)
          ISYMST_FL = ISYMEIG_KRCI(JFLST)
         HFCONST = DCMPLX(XPROPMAT(INIST,JFLST,IND_PRP(I),1),
     &                    XPROPMAT(INIST,JFLST,IND_PRP(I),2))
CTF
C   prefactor = mu / ( 2c I mp )
         read(XLABEL(2:2),*) INUC
         PREFACNUC = VKRCI_NUCMAGMOM(INUC) /
     &               ( 2 * CVAL * VKRCI_NUCSPIN(INUC) * XFMP )
         if (IVEC.eq.IONE) print*,' PREFACNUC = ',PREFACNUC
C   The hyperfine energy splitting:
         HFENERGY = HFCONST * PREFACNUC
C   For a level energy, multiply by 1/2[F(F+1)-J(J+1)-I(I+1)]
C   and add to the electronic energy
C
          write(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F20.12,5X,1F20.12)')
     & ' MHYP  state',INIST,'(',XREPEIG(ISYMST_IN),')         state',
     &      JFLST,'(',XREPEIG(ISYMST_FL),')   ',
     &     HFCONST,HFENERGY
       END DO
        write(lupri,'(a)')
     & '  __________________________________________'//
     & '    ___________________________________________________'//
     & '_______'
         END IF
       END DO
C
      end if
C     =======================================================
C     6. electron-Nucleus scalar-pseudoscalar interaction
C     Written by Malaya K. Nayak, Feb 28, 2014
C     =======================================================
      if (doensps_krci) then
C
C       print title
        write(lupri,*)
        write(lupri,'(/28x,a/28x,a/28x,a)') 
     & ' ********************************************************',
     & ' *********  e-N S-PS interaction Parameter  *************',
     & ' ********************************************************'
C
      DO I = 1, NPROP_KRCI
         INDXPR = LPROP_KRCI(I)
         XLABEL = PRPNAM(INDXPR)(1:8)
         IF (XLABEL(4:6).EQ.'SPS') THEN
         IND_PRP(I) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
C
C       ----------------------------
C       print eNSPS expectation value
C       ----------------------------
        write(lupri,'(/a,A6,a/a)')
     & '  <state| ',XLABEL,' i[gamma0 gamma5]rho_N |state>'//
     & '    e-N S-PS Parameter [a.u.]        Norm [a.u.]'
        write(lupri,'(a)')
     & '  __________________       _________________'//
     & '    ___________________________________________________'//
     & '_______'
C
       DO IVEC  = 1,NPROP_ROOTS_KRCI
          JVEC  = IVEC
          INIST = IVEC
          JFLST = JVEC
          ISYMST_IN = ISYMEIG_KRCI(INIST)
          ISYMST_FL = ISYMEIG_KRCI(JFLST)
         ENSPSMOM(IVEC) = DCMPLX(XPROPMAT(INIST,JFLST,IND_PRP(I),1),
     &                    XPROPMAT(INIST,JFLST,IND_PRP(I),2))
         EDMNORM = SQRT(ENSPSMOM(IVEC)**2)
C
          write(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F20.12,5X,1F20.12)')
     & ' ENSPS state',INIST,'(',XREPEIG(ISYMST_IN),')         state',
     &      JFLST,'(',XREPEIG(ISYMST_FL),')   ',
     &     ENSPSMOM(IVEC),EDMNORM
       END DO
        write(lupri,'(a)')
     & '  __________________________________________'//
     & '    ___________________________________________________'//
     & '_______'
         END IF
       END DO
C
      end if
C     =======================================================
C     7. Nuclear Magnetic Quadrupole Moment Interaction 
C     Written by Malaya K. Nayak, Mar 23, 2015
C     =======================================================
      if (donmqm_krci) then
C
C       print title
        write(lupri,*)
        write(lupri,'(/28x,a/28x,a/28x,a)')
     & ' ********************************************************',
     & ' *********  Nuclear Magnetic Quadrupole Moment  *********',
     & ' ********************************************************'
C
      DO I = 1, NPROP_KRCI
         INDXPR = LPROP_KRCI(I)
         XLABEL = PRPNAM(INDXPR)(1:8)
         IF (XLABEL(4:6).EQ.'MQM') THEN
         IND_PRP(I) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROP_KRCI,XLABEL)
C
C       ----------------------------
C       print nMQM expectation value
C       ----------------------------
        write(lupri,'(/a,A6,a/a)')
     & '    < state | ',XLABEL,' [alpha X EFG] | state >  '//
     & '    Nuc. MQM Constant [a.u.]       Norm [a.u.]'
        write(lupri,'(a)')
     & '  __________________       _________________'//
     & '    ___________________________________________________'//
     & '_______'
C
       DO IVEC  = 1,NPROP_ROOTS_KRCI
          JVEC  = IVEC
          INIST = IVEC
          JFLST = JVEC
          ISYMST_IN = ISYMEIG_KRCI(INIST)
          ISYMST_FL = ISYMEIG_KRCI(JFLST)
          EDMMOM = DCMPLX(XPROPMAT(INIST,JFLST,IND_PRP(I),1),
     &                    XPROPMAT(INIST,JFLST,IND_PRP(I),2))
          EDMNORM = SQRT(EDMMOM**2)
C
          write(LUPRI,'(A,I3,A,A4,A,I3,A,A4,A,1F20.12,5X,1F20.12)')
     & ' NMQM  state',INIST,'(',XREPEIG(ISYMST_IN),')         state',
     &      JFLST,'(',XREPEIG(ISYMST_FL),')   ',
     &      EDMMOM,EDMNORM
       END DO
        write(lupri,'(a)')
     & '  __________________________________________'//
     & '    ___________________________________________________'//
     & '_______'
         END IF
       END DO
C
      end if
C     =======================================================
C
C     release dynamically allocated memory
      if(print_ev) call dealloc(INDX)
C 
      CALL QEXIT('XKRCIPRPANA')
 1010 FORMAT(1X,I5,2X,F22.10,F16.12,F26.12,F28.12,'      (',A4,')')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck FLIP_BLOCK */
      SUBROUTINE FLIP_BLOCK(PROP,ACHAR,IPRINT)
C
C     Flip the Off-diagonal blocks of property matrix to Diagonal blocks
C     for parallel components of Hyperfine and Anapole moment operators. 
C
C     Written by Malaya K. Nayak Sept 04 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dgroup.h"

      DIMENSION PROP(2*NASHT,2*NASHT,NZ)
      CHARACTER*8 ACHAR
C
      MORB = 2*NASHT
      IF ((ACHAR(4:6).EQ.'HYP'.OR.ACHAR(4:6).EQ.'ANA').AND.
     &    (ACHAR(1:1).NE.'Z')) THEN
      WRITE(6,*)'Diagonal and Off-diagonal blocks are exchanged 
     &           for ', ACHAR
C
      IF ((ACHAR(4:6).EQ.'HYP'.AND.ACHAR(1:1).EQ.'Y').OR.
     &    (ACHAR(4:6).EQ.'ANA'.AND.ACHAR(1:1).EQ.'X')) THEN
      DO I=1,MORB
      DO J=1,MORB
        IF (NZ.EQ.1) THEN
        PROP(I,J,1) =-PROP(I,J,1)
        ELSE
        PROP(I,J,1) =-PROP(I,J,1)
        PROP(I,J,2) =-PROP(I,J,2)
        END IF
      END DO
      END DO
      END IF
C      
      MBY2 = MORB/2
      DO I=1,MORB
      DO J=1,MORB
      TEMPR = 0.0D0
      TEMPI = 0.0D0
      IF (J .LE. MBY2) THEN
      K = I
      L = J+MBY2
        IF (NZ .EQ. 1) THEN
        TEMPR = PROP(I,J,1)
        PROP(I,J,1) = PROP(K,L,1)
        ELSE
        TEMPR = PROP(I,J,1)
        TEMPI = PROP(I,J,2)
        PROP(I,J,1) = PROP(K,L,1)
        PROP(I,J,2) = PROP(K,L,2)
        END IF
      IF ((ACHAR(4:6).EQ.'HYP'.AND.ACHAR(1:1).EQ.'Y').OR.
     &    (ACHAR(4:6).EQ.'ANA'.AND.ACHAR(1:1).EQ.'X')) THEN
        IF (NZ .EQ. 1) THEN
        PROP(K,L,1) = TEMPR 
        ELSE
        PROP(K,L,1) = TEMPR
        PROP(K,L,2) = TEMPI
        END IF
      ELSE
        IF (NZ .EQ. 1) THEN
        PROP(K,L,1) =-TEMPR
        ELSE
        PROP(K,L,1) = TEMPR
        PROP(K,L,2) =-TEMPI
        ENDIF
      END IF
      END IF
      END DO
      END DO
      IF ((ACHAR(4:6).EQ.'HYP'.AND.ACHAR(1:1).EQ.'X').OR.
     &    (ACHAR(4:6).EQ.'ANA'.AND.ACHAR(1:1).EQ.'Y')) THEN
      WRITE(6,*)'Sign of the lower diagonal block is changed 
     &           for ', ACHAR
      END IF
      END IF
      IF ( IPRINT .GE. 30 ) THEN
         CALL HEADER('New Output matrix:',-1)
         CALL PRQMAT(PROP,2*NASHT,2*NASHT,2*NASHT,2*NASHT,MIN(NZ,2),
     $        IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck GET_OMEGA_EXP */
      SUBROUTINE GET_OMEGA_EXP(OMEGAMAT,op_sz,op_lz,NOMEGAST,NZ,
     &                         XPROPMAT,NPROPX,
     &                         XREPOM,IALL_SYM,NKRCI_CIROOTS,naelec,
     &                         IOMVAL,WORK,LWORK)
C***********************************************************************
C
C     Diagonalize matrix containing j_z = s_z + l_z expectation values.
C
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "krciprop.h"
#include "dcbxpr.h"
      CHARACTER*4 XREPOM(*)
      DIMENSION OMEGAMAT(NOMEGAST,NOMEGAST,*), NKRCI_CIROOTS(*)
      DIMENSION op_sz(NOMEGAST,NOMEGAST,*), op_lz(NOMEGAST,NOMEGAST,*)
      DIMENSION XPROPMAT(NOMEGAST,NOMEGAST,NPROPX,*), WORK(*), IOMVAL(*)
      INTEGER IND_DIP(2)
      CHARACTER XLABEL*16, om_string*5, om_divst*2
#include "memint.h"
      CALL QENTER('GET_OMEGA_EXP')
      NXTEST = 0
C
C     determine index for sigma_z and l_z
      XLABEL     = 'z-spin'
      IND_DIP(1) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROPX,XLABEL)
      XLABEL     = 'z-angular mom'
      IND_DIP(2) = IXPRPINDEX(PRPNAM,LPROP_KRCI,NPROPX,XLABEL)
C
C     copy sz-matrix part and lz-matrix part to OMEGAMAT
      DO IZ = 1, NZ
        DO I = 1, NOMEGAST
C         ... extract sz and copy to OMEGAMAT
          CALL DCOPY(NOMEGAST,XPROPMAT(1,I,IND_DIP(1),IZ),1,
     &               OMEGAMAT(1,I,IZ),1)
          CALL DCOPY(NOMEGAST,XPROPMAT(1,I,IND_DIP(1),IZ),1,
     &               op_sz(1,I,IZ),1)
C         ... now extract and add lz 
          CALL DAXPY(NOMEGAST,1.0D0,XPROPMAT(1,I,IND_DIP(2),IZ),1,
     &               OMEGAMAT(1,I,IZ),1)
          CALL DCOPY(NOMEGAST,XPROPMAT(1,I,IND_DIP(2),IZ),1,
     &               op_lz(1,I,IZ),1)
        END DO
      END DO
C
C     debug print
      NXTEST = 0
      IF( NXTEST .gt. 0 )THEN
C       complete XPROPMAT
        WRITE(LUPRI,*) ' complete XPROPMAT'
        CALL WRTMATMN(XPROPMAT,1,LXPRPKRCI,1,LXPRPKRCI,LUPRI)
C       complete OMEGAMAT
        WRITE(LUPRI,*) ' complete OMEGAMAT'
        CALL WRTMATMN(OMEGAMAT,1,NZ*NOMEGAST**2,1,NZ*NOMEGAST**2,LUPRI)
      END IF
C
C     print omega expectation values
      WRITE(LUPRI,'(/A/,A)')
     & '    omega expectation values for calculated eigenstates',
     & '  ________________________________________________________'
      WRITE(LUPRI,'(/A/,A/)')
     & '    state #        symmetry      <omega>      <s_z>     <l_z>',
     & '  ----------     -----------   ------------  -------  --------'
C
C     determine: boson or fermion wave function...
      if(mod(naelec,2).eq.0)then
        om_multf = 1.0D0
        om_divst = '  '
      else
        om_multf = 2.0D0
        om_divst = '/2'
      end if
!     initialize array which will contain the final omega values
      do initvl = 1, nomegast
        iomval(initvl) = -99999
      end do
C
C     offset index
      IX_AC = 0
C
C     loop over all symmetries and block diagonalize
      DO 100 IX = 1, IALL_SYM
C
C       only symmetries with non-zero number of eigenstates
        IF( NKRCI_CIROOTS(IX) .le. 0 ) GOTO 100
C
        NXSTATES = NKRCI_CIROOTS(IX)
C       ... allocate for x-x block diagonalization of <OMEGA>-matrix
        CALL MEMGET('REAL',KOMX,   NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGVCX,NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGX,  NXSTATES      ,WORK,KFREE,LFREE)
C       s_z
        CALL MEMGET('REAL',KszX,   NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGszX,NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGsz, NXSTATES      ,WORK,KFREE,LFREE)
C       l_z
        CALL MEMGET('REAL',KlzX,   NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGlzX,NZ*NXSTATES**2,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIGlz, NXSTATES      ,WORK,KFREE,LFREE)
C       ... initialize
        CALL DZERO(WORK(KOMX),   NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGVCX),NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGX),  NXSTATES)
        CALL DZERO(WORK(KszX),   NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGszX),NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGsz),  NXSTATES)
        CALL DZERO(WORK(KlzX),   NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGlzX),NZ*NXSTATES**2)
        CALL DZERO(WORK(KEIGlz),  NXSTATES)
C
C       extract x-x part of <OMEGA>-matrix
        CALL MATGATX(OMEGAMAT,NOMEGAST,WORK(KOMX),NXSTATES,IX_AC,NZ)
C       extract x-x part of <s_z>-matrix
        CALL MATGATX(op_sz,NOMEGAST,WORK(KszX),NXSTATES,IX_AC,NZ)
C       extract x-x part of <l_z>-matrix
        CALL MATGATX(op_lz,NOMEGAST,WORK(KlzX),NXSTATES,IX_AC,NZ)
C
#ifdef MOD_DEBUG
        WRITE(LUPRI,*) ' OMEGAMAT (block) before diagonalization'
        CALL WRTMATMN(WORK(KOMX),1,NZ*NXSTATES**2,1,
     &                NZ*NXSTATES**2,LUPRI)
        WRITE(LUPRI,*) ' szMAT (block) before diagonalization'
        CALL WRTMATMN(WORK(KszX),1,NZ*NXSTATES**2,1,
     &                NZ*NXSTATES**2,LUPRI)
        WRITE(LUPRI,*) ' lzMAT (block) before diagonalization'
        CALL WRTMATMN(WORK(KlzX),1,NZ*NXSTATES**2,1,
     &                NZ*NXSTATES**2,LUPRI)
#endif
C
C       diagonalize
        CALL DIAGC('X',WORK(KOMX),WORK(KEIGX),WORK(KEIGVCX),
     &             NXSTATES,WORK(KFREE),LFREE)
        CALL DIAGC('X',WORK(KszX),WORK(KEIGsz),WORK(KEIGszX),
     &             NXSTATES,WORK(KFREE),LFREE)
        CALL DIAGC('X',WORK(KlzX),WORK(KEIGlz),WORK(KEIGlzX),
     &             NXSTATES,WORK(KFREE),LFREE)
C
        NXTEST = 0
        IF( NXTEST .gt. 0 )THEN
          WRITE(LUPRI,*) ' OMEGAMAT after diagonalization'
          CALL WRTMATMN(WORK(KOMX),1,NZ*NXSTATES**2,1,
     &                  NZ*NXSTATES**2,LUPRI)
          WRITE(LUPRI,*) ' OMEGAEIGV after diagonalization'
          CALL WRTMATMN(WORK(KEIGX),1,NXSTATES,1,NXSTATES,LUPRI)
          WRITE(LUPRI,*) ' OMEGAEIGVEC after diagonalization'
          CALL WRTMATMN(WORK(KEIGVCX),1,NZ*NXSTATES**2,1,
     &                  NZ*NXSTATES**2,LUPRI)
        END IF
        NXTEST = 0
C
C       in the complex case "NZ==2" the eigenvalues are sorted in
C       ascending order - find absolute largest element in each
C       eigenvector to identify the corresponding CI eigenstate.
        IF( NZ .ge. 2 )THEN
C
#ifdef CODE_NOT_WORKING
C         FIXME: this "reordering" code does not work properly ...
C
C         list of CI eigenstates
          CALL MEMGET('REAL',KEIGTX,  NXSTATES      ,WORK,KFREE,LFREE)
          CALL DZERO(WORK(KEIGTX),NXSTATES)
          CALL GET_XOMEGALIST(WORK(KEIGTX),WORK(KEIGVCX),NXSTATES,
     &                        XREPOM(IX),WORK(KEIGX))
#endif
          WRITE(LUPRI,'(/A/)') 
     &    ' ***** OMEGA VALUES ARE SORTED in ascending order *****'
          WRITE(LUPRI,'(A)') 
     &    ' To assign the states properly, rerun the property part'//
     &    ' with less number of roots'
          WRITE(LUPRI,'(A/)') 
     &    ' and do thus a stepwise identification... ;)'
C
        END IF
C
C       print omega values
        do i = 1, nxstates
!         om_state  =  work(keigx+i-1) * om_multf
          om_state  = (work(KEIGsz+i-1) + work(KEIGlz+i-1))* om_multf
          iom_state =  nint(om_state)
          if(iabs(iom_state).eq.0) iom_state = 0
          write(om_string,'(i3,a2)') iom_state,om_divst
          write(lupri,'(i11,8x,a4,10x,a5,4x,1f10.4,1f10.4)')
     &          i,xrepom(ix),om_string,
     &          work(KEIGsz+i-1),
     &          work(KEIGlz+i-1)
!    &          op_sz(IX_AC+i,IX_AC+i,1),
!    &          op_lz(IX_AC+i,IX_AC+i,1)
!         save omega value on array (for later use 
!         when CI root selection is enabled)
          iomval(ix_ac+i) = iom_state
        end do
#ifdef MOD_DEBUG
        do i = 1, nxstates
          write(lupri,'(/a)')' *** omega values in original format ***'
          write(lupri,'(i11,8X,a4,6x,1f16.10)')
     &          i,xrepom(ix),work(keigx+i-1)
        end do
#endif
C
C       release memory
        CALL MEMREL('KRMCCNO.after x-x',WORK,1,KOMX,KFREE,LFREE)
C
C       keep track of active states
        IX_AC = IX_AC + NXSTATES
 100  CONTINUE
C
C     ... done
      CALL QEXIT('GET_OMEGA_EXP')
      END
C----------------------------------------------------------------------
C     end of file krmccan.F
C----------------------------------------------------------------------
