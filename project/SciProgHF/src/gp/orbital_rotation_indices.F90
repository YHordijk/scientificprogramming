module orbital_rotation_indices

   implicit none

   public get_orbital_rotation_indices_pp
   public get_orbital_rotation_indices_pn
   public set_orbital_rotation_indices
   public setxop
   public update_orbital_rotation_indices

   private

   integer, allocatable, target :: orbital_rotation_indices_pp(:)
   integer, allocatable, target :: orbital_rotation_indices_pn(:)

   logical :: is_initialized = .false.

contains

   subroutine check_if_interface_is_initialized()
      if (.not. is_initialized) then
         print *, 'ERROR: you try to access orbital_rotation_indices'
         print *, '       but this interface is not initialized'
         stop 1
      end if
   end subroutine

   function get_orbital_rotation_indices_pp() result(p)
      integer, pointer :: p(:)
      call check_if_interface_is_initialized()
      p => orbital_rotation_indices_pp
   end function

   function get_orbital_rotation_indices_pn() result(p)
      integer, pointer :: p(:)
      call check_if_interface_is_initialized()
      p => orbital_rotation_indices_pn
   end function

   subroutine update_orbital_rotation_indices(in_nzhope, &
                                              in_nzxope, &
                                              in_nzxopp, &
                                              in_jzxope,  &
                                              in_jzxopp)
      integer, intent(in) :: in_nzhope
      integer, intent(in) :: in_nzxope
      integer, intent(in) :: in_nzxopp
      integer, intent(in) :: in_jzxope(*)
      integer, intent(in) :: in_jzxopp(*)

#include "dcbxrs.h"

      nzhope = in_nzhope
      nzxope = in_nzxope
      nzxopp = in_nzxopp

      if (nzhope > 0) then
         if (allocated(orbital_rotation_indices_pp)) deallocate(orbital_rotation_indices_pp)
         allocate(orbital_rotation_indices_pp(nzhope*2))
         orbital_rotation_indices_pp(1:nzhope*2) = in_jzxope(1:nzhope*2)
      end if

      if (nzxopp > 0) then
         if (allocated(orbital_rotation_indices_pn)) deallocate(orbital_rotation_indices_pn)
         allocate(orbital_rotation_indices_pn(nzxopp*2))
         orbital_rotation_indices_pn(1:nzxopp*2) = in_jzxopp(1:nzxopp*2)
      end if

   end subroutine

      SUBROUTINE SETXOP(INDSTR,SKIPEE,SKIPEP,GAS,IOPSY,NOROT,           &
                   NXOPE,NHOPE,NXOPP,                  &
                   IPRINT)
!***********************************************************************
!
!     Set orbital rotation index vectors
!
!     Orbitals are classified as:
!       inactive  : doubly occupied in all determinants
!       active    : singly occupied in some determinants
!       secondary : unoccupied in all determinants
!       positronic: positronic orbitals
!
!     Ordering is:
!       e-e rotations: (skipped if SKIPEE = T)
!          IA : inactive-active
!          IS : inactive-secondary
!          AS : active-secondary
!          AB : active-active
!       e-p rotations: (skipped if SKIPEP = T)
!          IP : inactive-positronic
!          AP : active-positronic
!         (SP : secondary-positronic   -  only for QED=T)
!
!     NXOPE is the number of non-redundant e-e rotations.
!     NHOPE is the total number of rotations.
!
!     Written by  Hans Joergen Aa. Jensen and T.Saue 27-Jan-1995
!     Last revision: July 21 1997 - tsaue
!
!***********************************************************************

#include "priunit.h"
!
! Used from common blocks:
!  DGROUP: NFSYM
!  DCBHAM: QED
!  DCBORB: IORB,NPSH,NFRO,NISH,NASH
!
#include "dgroup.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbbas.h"
!
      LOGICAL SKIPEP,SKIPEE,GAS
      CHARACTER INDSTR(3,2)*72
      integer :: NSTR(0:2,3,2)
      integer :: NOROT(*)
      integer, allocatable :: ibuf(:)
      integer              :: idummy(2)
      integer              :: nxope, nxopp, nhope
      integer              :: nbuf
      integer              :: iprint
      integer              :: iopsy
!
      CALL QENTER('SETXOP')
      NXOPE =  0
      NXOPP =  0
      NHOPE =  0
!
!     Allocate buffer areas
!     =====================
!
      NBUF = N2BBASX
      allocate(ibuf(mxfbas*3*nfsym))
!
!     ***********************************
!     *** Electron-electron rotations ***
!     ***********************************
!
      if (allocated(orbital_rotation_indices_pp)) deallocate(orbital_rotation_indices_pp)
      if (skipee) then
         allocate(orbital_rotation_indices_pp(2))
         orbital_rotation_indices_pp = 0
      else
!        ... Unpack strings
         CALL UNPSTR(INDSTR,3,IBUF,MXFBAS,NSTR,IPRINT)
         ! first get length (dryrun=.true.)
         CALL SETXEE(idummy,    IOPSY,NXOPE,NHOPE,ibuf,NSTR,GAS,IPRINT,NOROT,.true.)
         ! then fill array
         allocate(orbital_rotation_indices_pp(nhope*2))
         orbital_rotation_indices_pp = 0
         CALL SETXEE(orbital_rotation_indices_pp,IOPSY,NXOPE,NHOPE,ibuf,NSTR,GAS,IPRINT,NOROT,.false.)
      end if
!
!     ***********************************
!     *** Electron-positron rotations ***
!     ***********************************
!
      if (allocated(orbital_rotation_indices_pn)) deallocate(orbital_rotation_indices_pn)
      if (skipep) then
         allocate(orbital_rotation_indices_pn(2))
         orbital_rotation_indices_pn = 0
      else
!        ... Unpack strings
         CALL UNPSTR(INDSTR,3,IBUF,MXFBAS,NSTR,0)
         ! first get length (dryrun=.true.)
         CALL SETXEP(idummy,    IOPSY,NXOPP,QED,IBUF,NSTR,IPRINT,NOROT,.true.)
         ! then fill array
         allocate(orbital_rotation_indices_pn(nxopp*2))
         orbital_rotation_indices_pn = 0
         CALL SETXEP(orbital_rotation_indices_pn,IOPSY,NXOPP,QED,IBUF,NSTR,IPRINT,NOROT,.false.)
      end if

      deallocate(ibuf)
      is_initialized = .true.

      CALL QEXIT('SETXOP')

      end subroutine



      SUBROUTINE SETXEE(JXOPE,IOPSY,NXOPE,NHOPE,IND,NIND,GAS,IPRINT,    &
                   NOROT,dryrun)
!
!     Generate electron-electron rotations
!
!     Written by T.Saue July 18 1997
!
!***********************************************************************
#include "priunit.h"
#include "maxash.h"
!
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbibn.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbxrs.h"
!
      integer :: JXOPE(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      integer :: NOROT(*)
      integer :: nhope, nxope, iprint
      LOGICAL GAS
      LOGICAL dryrun
      integer :: isym, jsym
      integer :: isti, iendi, ista, ienda
      integer :: jsti, jendi, jsta, jenda
      integer :: jsts, jends
      integer :: ia, ja, js, in, j, iaa, jaa
      integer :: indi, inda, inds
      integer :: iopsy
!
      integer :: JACAC(MAXASH,MAXASH)
!
      CALL QENTER('SETXEE')
!
!
!     ************************************
!     *** Find non-redundant rotations ***
!     ************************************
!
!     No active-active rotations are included here,
!     if not CAS (told with GAS variable) these are added
!     in RASGWOP below.
!
      NXOPE = 0
      DO 10 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
!
         ISTI  = NIND(1,1,ISYM) + 1
         IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
         JSTS  = NIND(1,3,JSYM) + 1
         JENDS = NIND(1,3,JSYM) + NIND(2,3,JSYM)
!
         IF(NIND(1,1,ISYM).NE.0) THEN
            WRITE(LUPRI,'(A,I5)')                                       &
            '*** WARNING(SETXEE): Inactive positrons detected:',        &
            NIND(1,1,ISYM)
         ENDIF
         IF(NIND(2,1,ISYM).EQ.0) GOTO 300
!
!        For Sternheimer complement add inactive-inactive elements
!        =========================================================
         IF(STERNC) THEN
           JSTI  = NIND(1,1,JSYM) + 1
           JENDI = NIND(1,1,JSYM) + NIND(2,1,JSYM)
           DO IN = ISTI,IENDI
             INDI = IND(IN,1,ISYM)
             DO JA = JSTI,JENDI
               INDA = IND(JA,1,JSYM)
               IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDA) .EQ. 0) THEN
                  NXOPE          = NXOPE + 1
                  if (.not. dryrun) then
                     JXOPE(1,NXOPE) = INDI
                     JXOPE(2,NXOPE) = INDA
                  end if
               END IF
             ENDDO
           ENDDO
         ENDIF
!
!        Inactive - active elements
!        ==========================
!
 100     CONTINUE
         IF(NIND(1,2,JSYM).NE.0) THEN
            WRITE(LUPRI,'(A,I5)')                                       &
            '*** WARNING(SETXEE): Active positrons detected:',          &
            NIND(1,2,JSYM)
         ENDIF
!                 inactive               active
         IF (NIND(2,1,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 200
         DO IN = ISTI,IENDI
            INDI = IND(IN,1,ISYM)
            DO JA = JSTA,JENDA
               INDA = IND(JA,2,JSYM)
               IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDA) .EQ. 0) THEN
                  NXOPE          = NXOPE + 1
                  if (.not. dryrun) then
                     JXOPE(1,NXOPE) = INDI
                     JXOPE(2,NXOPE) = INDA
                  end if
               END IF
            ENDDO
         ENDDO
!
!        Inactive-secondary elements
!        ===========================
!
 200     CONTINUE
!                 inactive               secondary
         IF (NIND(2,1,ISYM).EQ.0.OR.NIND(2,3,JSYM).EQ.0) GO TO 300
         DO IN = ISTI,IENDI
            INDI = IND(IN,1,ISYM)
            DO JS = JSTS,JENDS
               INDS = IND(JS,3,JSYM)
               IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDS) .EQ. 0) THEN
                  NXOPE          = NXOPE + 1
                  if (.not. dryrun) then
                     JXOPE(1,NXOPE) = INDI
                     JXOPE(2,NXOPE) = INDS
                  end if
               END IF
            ENDDO
         ENDDO
!
!        Active-secondary elements
!        =========================
!
 300     CONTINUE
!                 active                 secondary
         IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,3,JSYM).EQ.0) GO TO 10
         DO IA = ISTA,IENDA
            INDA = IND(IA,2,ISYM)
            DO JS = JSTS,JENDS
               INDS = IND(JS,3,JSYM)
               IF (NOROT(INDA) .EQ. 0 .AND. NOROT(INDS) .EQ. 0) THEN
                  NXOPE          = NXOPE+1
                  if (.not. dryrun) then
                     JXOPE(1,NXOPE) = INDA
                     JXOPE(2,NXOPE) = INDS
                  end if
               END IF
            ENDDO
         ENDDO
 10   CONTINUE
!
!
!     **********************************************
!     *** Find redundant active-active rotations ***
!     **********************************************
!
!
      NHOPE = NXOPE
      CALL IZERO(JACAC,MAXASH**2)
!
      IF (GAS) THEN
!
!        Find non-redundant active-active rotations.
!
         CALL RGASWOP(JXOPE,IOPSY,NXOPE,IND,NIND,JACAC,IPRINT,NOROT,dryrun)
!
      END IF
!
      NHOPE = NXOPE
!
      DO 20 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
         IF (JSYM .LT. ISYM) GOTO 20
!
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
!
!        Active-active elements
!        ----------------------
!
         IF (ISYM .EQ. JSYM) THEN
            IF (NIND(2,2,ISYM).LE.1) GO TO 20
            DO IA = ISTA, IENDA
               IAA = IND(IA,2,ISYM)
               DO JA = IA + 1, IENDA
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.       &
                  NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NHOPE = NHOPE + 1
                     if (.not. dryrun) then
                        JXOPE(1,NHOPE) = IAA
                        JXOPE(2,NHOPE) = JAA
                     end if
                  END IF
               END DO
            END DO
         ELSE
            IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 20
            DO IA = ISTA,IENDA
               IAA = IND(IA,2,ISYM)
               DO JA = JSTA,JENDA
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.       &
                  NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NHOPE = NHOPE + 1
                     if (.not. dryrun) then
                        JXOPE(1,NHOPE) = IAA
                        JXOPE(2,NHOPE) = JAA
                     end if
                  END IF
               ENDDO
            ENDDO
         END IF
 20   CONTINUE
!
!     Print section
!     =============
!
      if (.not. dryrun) then
          IF(IPRINT.GE.2) THEN
            CALL HEADER('Output from SETXEE',-1)
            WRITE(LUPRI,'(A,I10)')                                          &
            '# of e-e rotations (non-redundant): ',NXOPE
            IF (IPRINT.GE.3) THEN
              WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,NXOPE)
            ENDIF
            WRITE(LUPRI,'(A,I10)') '# of extra e-e rotations (redundant): ',&
             NHOPE - NXOPE
            IF (IPRINT.GE.3) THEN
              WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=NXOPE+1,NHOPE)
            ENDIF
          ENDIF
      end if
!
      CALL QEXIT('SETXEE')
      end subroutine

      SUBROUTINE SETXEP(JXOPP,IOPSY,NXOPP,QED,IND,NIND,IPRINT,NOROT,dryrun)
!***********************************************************************
!
!     Generate electron-positron rotations
!
!     Written by T.Saue July 21 1997
!
!***********************************************************************
#include "priunit.h"
!
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
!
      LOGICAL QED, dryrun
      integer :: JXOPP(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      integer :: NOROT(*)
      integer :: nxopp, iprint
      integer :: indi, iopsy, isym, jsym, inds, is, iends, ists, inda
      integer :: isti, iendi, ista, ienda, in, jp, j, ia, indp
!
      CALL QENTER('SETXEP')
!
      NXOPP = 0
      INDI  = 0
      DO 10 ISYM = 1,NFSYM
        JSYM  = MOD(ISYM+IOPSY,2) + 1
!
        ISTI  = NIND(1,1,ISYM) + 1
        IENDI = NIND(1,1,ISYM) + NIND(2,1,ISYM)
        ISTA  = NIND(1,2,ISYM) + 1
        IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!
!       Inactive-positronic elements
!       ===========================
!
 100    CONTINUE
        IF (NIND(2,1,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 200
        DO IN = ISTI,IENDI
           INDI =  IND(IN,1,ISYM)
           DO JP = 1,NIND(1,3,JSYM)
              INDP = IND(JP,3,JSYM)
              IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDP) .EQ. 0) THEN
                 NXOPP          = NXOPP+1
                 if (.not. dryrun) then
                    JXOPP(1,NXOPP) = INDI
                    JXOPP(2,NXOPP) = INDP
                 end if
              END IF
           ENDDO
        ENDDO
!
!       Active-positronic elements
!       ==========================
!
 200    CONTINUE
        IF (NIND(2,2,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 300
        DO IA = ISTA,IENDA
           INDA = IND(IA,2,ISYM)
           DO JP = 1,NIND(1,3,JSYM)
              INDP = IND(JP,3,JSYM)
              IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDP) .EQ. 0) THEN
                 NXOPP          = NXOPP+1
                 if (.not. dryrun) then
                    JXOPP(1,NXOPP) = INDA
                    JXOPP(2,NXOPP) = INDP
                 end if
              END IF
           ENDDO
        ENDDO
!
!       Secondary-positronic elements
!       =============================
!
 300    CONTINUE
        IF (QED) THEN
           IF (NIND(1,3,ISYM).EQ.0.OR.NIND(1,3,JSYM).EQ.0) GO TO 10
           ISTS  = NIND(1,3,ISYM) + 1
           IENDS = NIND(1,3,ISYM) + NIND(2,3,ISYM)
           DO IS = ISTS,IENDS
              INDS = IND(IS,1,ISYM)
              DO JP = 1,NIND(1,3,JSYM)
                 INDP = IND(JP,3,JSYM)
                 IF (NOROT(INDI) .EQ. 0 .AND. NOROT(INDP) .EQ. 0) THEN
                    NXOPP          = NXOPP+1
                    if (.not. dryrun) then
                       JXOPP(1,NXOPP) = INDS
                       JXOPP(2,NXOPP) = INDP
                    end if
                 END IF
              ENDDO
           ENDDO
        ENDIF
!
 10   CONTINUE
!
!     Print section
!     =============
!
      if (.not. dryrun) then
         IF(IPRINT.GE.2) THEN
           CALL HEADER('Output from SETXEP',-1)
           WRITE(LUPRI,'(A,I10)') '# of e-p rotations: ',NXOPP
           IF(IPRINT.GE.3) THEN
             WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,NXOPP)
           ENDIF
         ENDIF
      end if
!
      CALL QEXIT('SETXEP')

   end subroutine

      SUBROUTINE UNPSTR(INDSTR,NTYP,IND,LIND,NIND,IPRINT)
!***********************************************************************
!
!     General routine for unpacking of orbital strings
!
!     Written by T.Saue July 19 1997
!
!***********************************************************************
#include "priunit.h"
!
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER INDSTR(NTYP,NFSYM)*72
      integer :: IND(LIND,NTYP,NFSYM),NIND(2,NTYP,NFSYM)
      integer :: iprint, ntyp, lind
      integer :: ifrp, ityp, nvec, i, ioff, joff
!
      DO IFRP = 1,NFSYM
        DO ITYP = 1,NTYP
          NVEC = 1
          CALL  NUMLST(INDSTR(ITYP,IFRP),IND(1,ITYP,IFRP),              &
                  NFBAS(IFRP,0),-NFBAS(IFRP,2),NFBAS(IFRP,1),           &
                  IFRP,NVEC)
          ! count number of p-type and of e-type orbitals
          CALL ORBCNT(IND(1,ITYP,IFRP),NVEC,                            &
                 NPSH(IFRP),NESH(IFRP),                                 &
                 NIND(1,ITYP,IFRP),NIND(2,ITYP,IFRP))
!
!         Convert to absolute indices of positronic orbitals.
!         (For ifrp=1: this means mapping -npsh(1):-1 to 1:npsh(1),
!          that is, -1 -> npsh(1), -2 -> npsh(1)-1, ..., -npsh(1) -> 1;
!          for ifrp=2 we add the total number of orbitals in ifrp=1)
!
          IOFF = NPSH(IFRP) + 1 + IORB(IFRP)
          DO I = 1,NIND(1,ITYP,IFRP)
            IND(I,ITYP,IFRP) = IOFF + IND(I,ITYP,IFRP)
          ENDDO
!
!         Convert to absolute indices of electronic orbitals.
!         (For ifrp=1: this means mapping 1:nesh(1) to npsh(1)+1:npsh(1)+nesh(1),
!          for ifrp=2 we add the total number of orbitals in ifrp=1)
!
          IOFF = NPSH(IFRP) + IORB(IFRP)
          JOFF = NIND(1,ITYP,IFRP)
          DO I = 1,NIND(2,ITYP,IFRP)
            IND(I+JOFF,ITYP,IFRP) = IOFF + IND(I+JOFF,ITYP,IFRP)
          ENDDO
        ENDDO
      ENDDO
!
!     Print section
!     =============
!
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Output from UNPSTR',-1)
        DO ITYP = 1,NTYP
          WRITE(LUPRI,'(A,I5)') '* Orbital string no.',ITYP
          DO IFRP = 1,NFSYM
            WRITE(LUPRI,'(3X,A,A3)') '* Fermion ircop ',FREP(IFRP)
            IF(NIND(1,ITYP,IFRP).EQ.0) THEN
              WRITE(LUPRI,'(5X,A)') 'No positronic orbitals.'
            ELSE
              WRITE(LUPRI,'(5X,A)') 'Positronic orbitals:'
              WRITE(LUPRI,'(8X,12I5)')                                  &
              (IND(I,ITYP,IFRP),I=1,NIND(1,ITYP,IFRP))
            ENDIF
            IF(NIND(2,ITYP,IFRP).EQ.0) THEN
              WRITE(LUPRI,'(5X,A)') 'No electronic orbitals.'
            ELSE
              JOFF = NIND(1,ITYP,IFRP)
              WRITE(LUPRI,'(5X,A)') 'Electronic orbitals:'
              WRITE(LUPRI,'(8X,12I5)')                                  &
              (IND(I+JOFF,ITYP,IFRP),I=1,NIND(2,ITYP,IFRP))
            ENDIF
          ENDDO
        ENDDO
      ENDIF

   end subroutine

      SUBROUTINE RGASWOP(JXOPE,IOPSY,NXOPE,IND,NIND,JACAC,IPRINT,NOROT,dryrun)
!***********************************************************************
!
!     Generate non-redundant active-active rotations.
!
!     Written by J. Thyssen - Jan 18 2001
!
!     Revised GAS setup, Timo Fleig, May 2002
!
!***********************************************************************
#include "priunit.h"
#include "maxash.h"
!
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbibn.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
!
      integer :: JXOPE(2,*),IND(MXFBAS,3,NFSYM),NIND(2,3,NFSYM)
      integer :: JACAC(MAXASH,MAXASH)
      integer :: NOROT(*)
!
      integer :: IGAS(MAXASH)
      integer :: iprint, nxope, iopsy
      LOGICAL GASCAS(MXGAS), dryrun
      integer :: nxopes, isym, i, j, ii, ngst, jsym
      integer :: mincas, maxcas
      integer :: ista, ienda, jsta, jenda
      integer :: ia, ja, iaa, jaa

!
      NXOPES = NXOPE
!
!     Set up help arrays
!     ------------------
!
!     IGAS point from active orbital to which GAS space it belongs.
!
      II = 0
      DO ISYM = 1, NFSYM
         DO I = 1, NGAS_DC
            DO J = 1, NGSH(ISYM,I)
               II = II + 1
               IGAS(II) = I
            END DO
         END DO
      END DO
!
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Output from RGASWOP',-1)
         WRITE(LUPRI,'(1X,A)')                                          &
         'Pointer from active orbitals to GAS spaces:'
         DO I = 1, NASHT
            WRITE(LUPRI,'(5X,A,I3,A,I3)')                               &
            'Orbital ',I,' is in GAS space ',IGAS(I)
         END DO
      END IF
!
!     Determine whether a GAS space is a CAS space:
!
      DO I = 1, NGAS_DC
!
!        Number of orbitals in this GAS shell
!
         NGST = 0
         DO ISYM = 1, NFSYM
            NGST = NGST + NGSH(ISYM,I)
         END DO
!
!        Minimum number of electrons in this GAS space:
!
         MINCAS = 2*NGST - (NGASSP(2,I) - NGASSP(1,I))
         IF (I.EQ.NGAS_DC.AND.NGAS_DC.GT.1) THEN
           MINCAS = NAELEC - NGASSP(2,I-1)
         ELSE IF (I.EQ.NGAS_DC.AND.NGAS_DC.EQ.1) then
           MINCAS = NAELEC
         END IF
!        print*,'GAS ',I,' MINCAS,MAXCAS ',MINCAS,MAXCAS
!
!        Maximum number of holes in this GAS space:
!
         MAXCAS = 2*NGST - MINCAS
!
!        This is a GAS space if the minimum number of electrons
!        after this space is less than MINCAS - and -
!        if the maximum number of holes in this GAS space is larger
!        than MAXCAS:
!
         IF ( (NGASSP(1,I) .LE. MINCAS .AND.                            &
         (NGASSP(2,I) - NGASSP(1,I)) .GE. MAXCAS) .OR.                  &
         I .EQ. NGAS_DC ) THEN
            GASCAS(I) = .TRUE.
         ELSE
            GASCAS(I) = .FALSE.
         END IF
      END DO
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(1X,A)')                                          &
         'GAS space is a CAS space?'
         DO I = 1, NGAS_DC
            WRITE(LUPRI,'(5X,A,I3,A,L1)')                               &
            'GAS space ',I,' : ',GASCAS(I)
         END DO
      END IF
!
!     Algorithm:
!     ----------
!
!     No inter-GAS-space active-active rotations
!     No active-active rotations into CAS spaces
!     No active-active rotations from a CAS space into the NEXT GAS space.
!
!     Active-active elements
!     ----------------------
!
      DO 20 ISYM = 1,NFSYM
         JSYM  = MOD(ISYM+IOPSY,2) + 1
         IF (JSYM .LT. ISYM) GOTO 20
!
         ISTA  = NIND(1,2,ISYM) + 1
         IENDA = NIND(1,2,ISYM) + NIND(2,2,ISYM)
!
         JSTA  = NIND(1,2,JSYM) + 1
         JENDA = NIND(1,2,JSYM) + NIND(2,2,JSYM)
         IF (ISYM .EQ. JSYM) THEN
            IF (NIND(2,2,ISYM).LE.1) GO TO 20
            DO IA = ISTA, IENDA
               DO JA = IA + 1, IENDA
                  IAA = IND(IA,2,ISYM)
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.       &
                  IGAS(IDXG2U(IAA)) .NE. IGAS(IDXG2U(JAA)) .AND.        &
                  .NOT. GASCAS(IGAS(IDXG2U(JAA))) .AND.                 &
                  .NOT. (                                               &
                     GASCAS(IGAS(IDXG2U(IAA))) .AND.                    &
                  IGAS(IDXG2U(JAA)) .LE. IGAS(IDXG2U(IAA)+1)) .AND.     &
                  NOROT(IAA) .EQ. 0 .AND. NOROT(JAA) .EQ. 0) THEN
                     NXOPE = NXOPE + 1
                     if (.not. dryrun) then
                        JXOPE(1,NXOPE) = IAA
                        JXOPE(2,NXOPE) = JAA
                        JACAC(IDXG2U(IAA),IDXG2U(JAA)) = 1
                        JACAC(IDXG2U(JAA),IDXG2U(IAA)) = 1
                     end if
                  END IF
               END DO
            END DO
         ELSE
            IF (NIND(2,2,ISYM).EQ.0.OR.NIND(2,2,JSYM).EQ.0) GO TO 20
            DO IA = ISTA,IENDA
               DO JA = JSTA,JENDA
                  IAA = IND(IA,2,ISYM)
                  JAA = IND(JA,2,JSYM)
                  IF (JACAC(IDXG2U(IAA),IDXG2U(JAA)) .EQ. 0 .AND.       &
                  IGAS(IDXG2U(IAA)) .NE. IGAS(IDXG2U(JAA)) .AND.        &
                  .NOT. GASCAS(IDXG2U(JAA)) .AND.                       &
                  .NOT. (                                               &
                     GASCAS(IDXG2U(IAA)) .AND.                          &
                  IGAS(IDXG2U(JAA)) .LE. IGAS(IDXG2U(IAA)+1))) THEN
                     NXOPE = NXOPE + 1
                     if (.not. dryrun) then
                        JXOPE(1,NXOPE) = IAA
                        JXOPE(2,NXOPE) = JAA
                        JACAC(IDXG2U(IAA),IDXG2U(JAA)) = 1
                        JACAC(IDXG2U(JAA),IDXG2U(IAA)) = 1
                     end if
                  END IF
               ENDDO
            ENDDO
         END IF
 20   CONTINUE
!
!
      if (.not. dryrun) then
         IF (IPRINT .GE. 10) THEN
            WRITE(LUPRI,'(1X,A,I3,A)')                                     &
            'The following ',NXOPE - NXOPES,                               &
            ' active-active rotations were considered non-redundant: '
            DO I = NXOPES + 1, NXOPE
               WRITE(LUPRI,'(5X,2I4)')                                     &
               JXOPE(1,I), JXOPE(2,I)
            END DO
         END IF
      end if

   end subroutine

   subroutine set_orbital_rotation_indices(                            &
                                           do_pp,                      &
                                           do_pn,                      &
                                           length_pp,                  &
                                           length_pn,                  &
                                           orbital_rotation_vector_pp, &
                                           orbital_rotation_vector_pn  &
                                          )

      logical :: do_pp
      logical :: do_pn
      integer :: length_pp
      integer :: length_pn
      integer :: orbital_rotation_vector_pp(*)
      integer :: orbital_rotation_vector_pn(*)

      integer :: i, j

      if (allocated(orbital_rotation_indices_pp)) deallocate(orbital_rotation_indices_pp)
      if (do_pp) then
         allocate(orbital_rotation_indices_pp(length_pp*2))
         orbital_rotation_indices_pp(1:length_pp*2) = orbital_rotation_vector_pp(1:length_pp*2)
      else
         allocate(orbital_rotation_indices_pp(2))
         orbital_rotation_indices_pp = 0
      end if

      if (allocated(orbital_rotation_indices_pn)) deallocate(orbital_rotation_indices_pn)
      if (do_pn) then
         allocate(orbital_rotation_indices_pn(length_pn*2))
         orbital_rotation_indices_pn(1:length_pn*2) = orbital_rotation_vector_pn(1:length_pn*2)
      else
         allocate(orbital_rotation_indices_pn(2))
         orbital_rotation_indices_pn = 0
      end if

      is_initialized = .true.

   end subroutine

end module
