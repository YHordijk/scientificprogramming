!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_COUSAT_R (PHPH_BLOCK,EPS,DESREP,LENCOUPLBL,
     &                             INTBUF,NBUFS,DOXEXT,WTHR,IOCH)
!
      IMPLICIT INTEGER (A-Z)
!     IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!  Computes (real) 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!
!---------------Calling variables--------------------------------------
!
      REAL*8                           :: PHPH_BLOCK(:)
      REAL*8                           :: EPS(:)
      INTEGER                          :: DESREP
      INTEGER                          :: LENCOUPLBL,INTBUF,NBUFS
      LOGICAL                          :: DOXEXT
      REAL*8                           :: WTHR
      INTEGER                          :: IOCH
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep,k
      INTEGER*8            :: ix8,n8  !needed for the ivovvt variable!
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,nmvod

      INTEGER              :: mxno,mxnv,ioff
      INTEGER, dimension (:,:,:), allocatable  :: oolt
      INTEGER, dimension (:,:,:), allocatable  :: vvlt

      REAL*8, dimension(:), allocatable :: vooo1
      REAL*8, dimension(:), allocatable :: vovv1
      REAL*8, dimension(:), allocatable :: bufcol,bufcolex
!
!  variables for writing out ADC matrix
!
      REAL*8,  dimension(:), allocatable :: wbufr
      INTEGER, dimension(:), allocatable :: wbufi, wbufj
      INTEGER              :: IBFP
      CHARACTER*6          :: NAMEX = 'ADCXPM'
!
!---------------Interface area --------------------------------------
!
      interface

        INTEGER FUNCTION XCOLDET(i1)
          INTEGER                         :: i1
        END FUNCTION

        SUBROUTINE XCOLMAKE(i1,i2,ra1,
     &                      i3,i4,i5,i6,
     &                      i7,i8,
     &                      ra2,ra3,ia1,ia2)
          INTEGER                         :: i1,i2
          REAL*8, dimension(:)            :: ra1
          INTEGER                         :: i3,i4,i5,i6
          INTEGER                         :: i7,i8
          REAL*8, dimension(:)            :: ra2,ra3
          INTEGER, dimension(:,:,:)       :: ia1,ia2
        END SUBROUTINE XCOLMAKE

        SUBROUTINE XCOLWRIT(i1,r1,ra1,ra2,ia1,ia2,i2,i3,i4,i5,i6)
          INTEGER                         :: i1
          REAL*8                          :: r1
          REAL*8, dimension(:)            :: ra1,ra2
          INTEGER, dimension(:)           :: ia1,ia2
          INTEGER                         :: i2,i3,i4,i5,i6
        END SUBROUTINE XCOLWRIT

        SUBROUTINE XSATMAKE(i1,i2,i3,i4,i5,r1,
     &              ra1,ia1,ia2,ra2,ia3,ia4,
     &              l1)
          INTEGER                   :: i1,i2,i3,i4,i5
          REAL*8                    :: r1
          REAL*8, dimension(:)      :: ra1
          INTEGER, dimension(:)     :: ia1,ia2
          REAL*8, dimension(:)      :: ra2
          INTEGER,dimension(:,:,:)  :: ia3,ia4
          LOGICAL                   :: l1
        END SUBROUTINE XSATMAKE

        SUBROUTINE XCOLFLSH(i1,i2,i3,i4,ra1,ia1,ia2)
          INTEGER                         :: i1,i2,i3,i4
          REAL*8, dimension(:)            :: ra1
          INTEGER, dimension(:)           :: ia1,ia2
        END SUBROUTINE XCOLFLSH

      end interface
!
!---------------Executable code--------------------------------------
!
!________________________________________________________
!|     allocate ADC MATRIX READ/WRITE BUFFERS and counters
!|
!|
      allocate(wbufr(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 1 could not be allocated.'
      allocate(wbufi(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 2 could not be allocated.'
      allocate(wbufj(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 3 could not be allocated.'

      wbufr = 0.0d0
      wbufi = 0
      wbufj = 0

      IBFP = 0   ! pointer to actual buffer level
      NBUFS = 0  ! counter for full buffers
!|
!|
!|_______________________________________________________
!________________________________________________________
!|     create OOT and VVT lookup tables
!|     tables are small, computation is too clumsy
!|
      mxno = 0; mxnv = 0
      do irep=1,nrep  
        if(no(irep).gt.mxno) mxno=no(irep)
        if(nv(irep).gt.mxnv) mxnv=nv(irep)
      enddo
      allocate(oolt(mxno,mxno,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #a in BUILD_X_COUSAT_R'
      allocate(vvlt(mxnv,mxnv,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #b in BUILD_X_COUSAT_R'
      OOLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NO(IREP)
          DO I=J+1,NO(IREP)
            OOLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
      VVLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NV(IREP)
          DO I=J+1,NV(IREP)
            VVLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
!|
!|
!|_______________________________________________________
!
!
!  Attention*** Some arrays that have to store *huge* numbers 
!  are defined as Integer*8 in symm.inc ! The generic sorters
!  operate with unspecified INTEGER declarations. Transferring such
!  an array as an argument to the SRT_xxx results in an erroneous
!  behaviour of the sorting routine.
!
! fetch <VO||OO> integrals needed for coupling block 
! 
      n1=ivooot(nrep+1)
      allocate(vooo1(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vooo1!'
      call getvooo(vooo1)
!
! fetch <VO||VV> integrals needed for coupling block 
! and test for completeness
! 
      n8=ivovvt(nrep+1)
      allocate(vovv1(n8),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vovv1!'
      do ix8 = 1,n8
        vovv1(ix8) = -9999.0d0
      enddo
      call rdvovv (vovv1)
      do ix8 = 1,n8
        if(vovv1(ix8).eq. -9999.0d0)
     &      call quit('BUILD_X_COUSAT_R error: gap in VOVV stream!')
      enddo
      write(*,*) 'BUILD_X_COUSAT_R: VOVV stream complete.'
!
! length of main block:
!
      NMVOD = MVO(desrep)
!
! determine actual length of coupling block
! and allocate column buffer
!
      lencouplbl =  XCOLDET(desrep)
      WRITE(IW,*) '   *)) Current symmetry:          ',desrep
      WRITE(IW,*) '   *)) Block length of main block:',NMVOD
      WRITE(IW,*) '   *)) Block length of coupling block:',lencouplbl
      WRITE(IW,'(A,I8,A)') '    *)) Matrix dimension: [',
     &        lencouplbl + NMVOD,' ]**2'

      ALLOCATE(BUFCOL(lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #3 in BUILD_X_COUSAT_R'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #4 in BUILD_X_COUSAT_R'
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2P2H/PH coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
      mjcol = 1
      do jrep = 1,nrep

        brep = multb(desrep+nrep,jrep,2)
        if(MULTB(brep,jrep,2).eq.desrep) then
!         write(*,*) 'BJ rep generation ok.'
        else
          stop 'BJ rep generation failed.'
        endif

        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)

          bufcol   = 0.0d0
          bufcolex = 0.0d0

          CALL XCOLMAKE(desrep,lencouplbl,bufcol,
     &                  jrep,brep,jfun,bfun,
     &                  mxno,mxnv,
     &                  vooo1,vovv1,oolt,vvlt)
!
! combine main block column with coupling block column
! istart is the start address for picking the right main block
! data to corresponding symetry DESREP.
!
          istart = JVOVO(desrep) + (mjcol-1)*nmvod
          bufcolex(1:nmvod)=phph_block(istart+1:istart+nmvod)
          bufcolex(nmvod+1:nmvod+lencouplbl)=
     &        bufcol(1:lencouplbl)
!
! write out full column.
!
          CALL XCOLWRIT(IOCH,wthr,bufcolex,wbufr,wbufi,wbufj,
     &                  IBFP,NBUFS,LENCOUPLBL+NMVOD,mjcol,INTBUF)
          mjcol = mjcol + 1

        enddo   !bfun
        enddo   !jfun

      enddo  ! jrep

      IF((mjcol-1).ne.mvo(desrep))
     &  stop 'error 3 in build_x_cousat_r!'
!________________________________________________________
!|
!|
!| from here the main and coupling block are written to
!| the external file. The SAT entries are also written to
!| the same file.
!| ATT**  we can release, and should release *all*
!| memory which we do not need anymore.
!| Due to the large size of the matrices we can not afford
!| a four-indexed lookup table as it was used in the two-particle
!| propagator!
!|
!|_______________________________________________________


      CALL XSATMAKE(IOCH,IBFP,NBUFS,INTBUF,DESREP,WTHR,
     &              WBUFR,WBUFI,WBUFJ,EPS,OOLT,VVLT,
     &              DOXEXT)
!
!  the SAT block in DESREP is now also on the file
!  and we can return in this symmetry and diagonalize...
!
!
!  flush buffer
!
      IF(IBFP.gt.0) THEN
        CALL XCOLFLSH(IOCH,IBFP,NBUFS,INTBUF,WBUFR,WBUFI,WBUFJ)
      ENDIF

      write(*,*) '--------------------------------------------'
      write(*,*) '---- Coupling and Satellite blocks finished.'
      write(*,*) '---- Number of required buffers:',NBUFS
      write(*,*) '--------------------------------------------'

      deallocate(bufcolex)
      deallocate(bufcol)
      deallocate(vovv1)
      deallocate(vooo1)
      deallocate(vvlt)
      deallocate(oolt)
!
! release storage for matrix writing
!
      deallocate(wbufj)
      deallocate(wbufi)
      deallocate(wbufr)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_COUSAT_C (PHPH_BLOCK,EPS,DESREP,LENCOUPLBL,
     &                             INTBUF,NBUFS,DOXEXT,WTHR,IOCH)
!
      IMPLICIT INTEGER (A-Z)
!     IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!  Computes (complex) 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!
!---------------Calling variables--------------------------------------
!
      REAL*8                           :: PHPH_BLOCK(:)
      REAL*8                           :: EPS(:)
      INTEGER                          :: DESREP
      INTEGER                          :: LENCOUPLBL,INTBUF,NBUFS
      LOGICAL                          :: DOXEXT
      REAL*8                           :: WTHR
      INTEGER                          :: IOCH
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep,k
      INTEGER*8            :: n8  !needed for the ivovvt variable!
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,nmvod

      INTEGER              :: mxno,mxnv,ioff
      INTEGER, dimension (:,:,:), allocatable  :: oolt
      INTEGER, dimension (:,:,:), allocatable  :: vvlt

      REAL*8, dimension(:), allocatable     :: vooo1
      REAL*8, dimension(:), allocatable     :: vovv1
      complex*16, dimension(:), allocatable :: bufcol,bufcolex
!
!  variables for writing out ADC matrix
!
      complex*16,dimension(:), allocatable   :: wbufr
      INTEGER, dimension(:), allocatable     :: wbufi, wbufj
      INTEGER              :: IBFP
      CHARACTER*6          :: NAMEX = 'ADCXPM'
!
!  auxiliary array for complex data type
!
      complex*16, allocatable, dimension(:)    :: phph_blc

      real*8,parameter        :: rzero=(0.0d0)
      complex*16,parameter    :: czero=(rzero,rzero)
!
!---------------Interface area --------------------------------------
!
      interface

        INTEGER FUNCTION XCOLDET(i1)
          INTEGER                         :: i1
        END FUNCTION

        SUBROUTINE XCOLMAKE_C(i1,i2,ca1,
     &                        i3,i4,i5,i6,
     &                        i7,i8,
     &                        ra2,ra3,ia1,ia2)
          INTEGER                         :: i1,i2
          complex*16, dimension(:)        :: ca1
          INTEGER                         :: i3,i4,i5,i6
          INTEGER                         :: i7,i8
          REAL*8, dimension(:)            :: ra2,ra3
          INTEGER, dimension(:,:,:)       :: ia1,ia2
        END SUBROUTINE XCOLMAKE_C

        SUBROUTINE XCOLWRIT_C(i1,r1,ca1,ca2,ia1,ia2,i2,i3,i4,i5,i6)
          INTEGER                         :: i1
          REAL*8                          :: r1
          complex*16, dimension(:)        :: ca1,ca2
          INTEGER, dimension(:)           :: ia1,ia2
          INTEGER                         :: i2,i3,i4,i5,i6
        END SUBROUTINE XCOLWRIT_C

        SUBROUTINE XSATMAKE_C (i1,i2,i3,i4,i5,r1,
     &              ca1,ia1,ia2,ra2,ia3,ia4,l1)
          INTEGER                     :: i1,i2,i3,i4,i5
          REAL*8                      :: r1
          complex*16, dimension(:)    :: ca1
          INTEGER, dimension(:)       :: ia1,ia2
          REAL*8, dimension(:)        :: ra2
          INTEGER,dimension(:,:,:)    :: ia3,ia4
          LOGICAL                     :: l1
        END SUBROUTINE XSATMAKE_C

        SUBROUTINE XCOLFLSH_C(i1,i2,i3,i4,ca1,ia1,ia2)
          INTEGER                         :: i1,i2,i3,i4
          complex*16, dimension(:)        :: ca1
          INTEGER, dimension(:)           :: ia1,ia2
        END SUBROUTINE XCOLFLSH_C

      end interface
!
!---------------Executable code--------------------------------------
!
      if(rcw.ne.2) then
        call quit('Internal error in BUILD_X_COUSAT_C !')
      endif
!___________________________________________________________
!|
!|  generate complex auxiliary array from original real one.

      n1 = jvovo(nrep+1)
      allocate(phph_blc(n1))
      phph_blc = czero
      do i=1,n1
        ix = 2*(i-1) + 1
        phph_blc(i)=dcmplx(phph_block(ix),phph_block(ix+1))
      enddo
!|
!|__________________________________________

!________________________________________________________
!|     allocate ADC MATRIX WRITE BUFFERS and counters
!|
!|
      allocate(wbufr(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 1 could not be allocated.'
      allocate(wbufi(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 2 could not be allocated.'
      allocate(wbufj(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 3 could not be allocated.'

      wbufr = czero
      wbufi = 0
      wbufj = 0
      IBFP = 0   ! pointer to actual buffer level
      NBUFS = 0  ! counter for full buffers
!|
!|
!|_______________________________________________________
!________________________________________________________
!|     create OOT and VVT lookup tables
!|     tables are small, computation is too clumsy
!|
      mxno = 0; mxnv = 0
      do irep=1,nrep  
        if(no(irep).gt.mxno) mxno=no(irep)
        if(nv(irep).gt.mxnv) mxnv=nv(irep)
      enddo
      allocate(oolt(mxno,mxno,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #a in BUILD_X_COUSAT_R'
      allocate(vvlt(mxnv,mxnv,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #b in BUILD_X_COUSAT_R'
      OOLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NO(IREP)
          DO I=J+1,NO(IREP)
            OOLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
      VVLT = 0
      DO IREP=1,NREP
        IOFF=1
        DO J=1,NV(IREP)
          DO I=J+1,NV(IREP)
            VVLT(I,J,IREP)=IOFF
            IOFF=IOFF+1
          ENDDO
        ENDDO
      ENDDO
!|
!|
!|_______________________________________________________
!
!
!  Attention*** Some arrays that have to store *huge* numbers 
!  Are defined as Integer*8 in symm.inc ! The generic sorters
!  operate with unspecified INTEGER declarations. Transferring such
!  an array as an argument to the SRT_xxx results in an erroneous
!  behaviour of the sorting routine.

!
! fetch <VO||OO> integrals needed for coupling block 
! 
      n1=ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vooo1!'
      call getvooo(vooo1)
!
! fetch <VO||VV> integrals needed for coupling block 
! 
      n8=ivovvt(nrep+1)
      allocate(vovv1(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vovv1!'
      call rdvovv (vovv1)
!
! length of main block:
!
      NMVOD = MVO(desrep)
!
! determine actual length of coupling block
! and allocate column buffer
!
      lencouplbl = XCOLDET(desrep)
      WRITE(IW,*) '   *)) Current symmetry:          ',desrep
      WRITE(IW,*) '   *)) Block length of main block:',NMVOD
      WRITE(IW,*) '   *)) Block length of coupling block:',lencouplbl
      WRITE(IW,'(A,I8,A)') '    *)) Matrix dimension: [',
     &        lencouplbl + NMVOD,' ]**2'

      ALLOCATE(BUFCOL(lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #3 in BUILD_X_COUSAT_R'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #4 in BUILD_X_COUSAT_R'
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2P2H/PH coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
      mjcol = 1
      do jrep = 1,nrep

        brep = multb(desrep+nrep,jrep,2)
        if(MULTB(brep,jrep,2).eq.desrep) then
!         write(*,*) 'BJ rep generation ok.'
        else
          stop 'BJ rep generation failed.'
        endif

        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)

          bufcol   = czero
          bufcolex = czero

          CALL XCOLMAKE_C(desrep,lencouplbl,bufcol,
     &                    jrep,brep,jfun,bfun,
     &                    mxno,mxnv,
     &                    vooo1,vovv1,oolt,vvlt)
!
! combine main block column with coupling block column
! istart is the start address for picking the right main block
! data to corresponding symetry DESREP.
!
          istart = JVOVO(desrep) + (mjcol-1)*nmvod
          bufcolex(1:nmvod)=phph_blc(istart+1:istart+nmvod)
          bufcolex(nmvod+1:nmvod+lencouplbl)=
     &        bufcol(1:lencouplbl)
!
! write out full column.
!
          CALL XCOLWRIT_C(IOCH,wthr,bufcolex,wbufr,wbufi,wbufj,
     &                    IBFP,NBUFS,LENCOUPLBL+NMVOD,mjcol,INTBUF)
          mjcol = mjcol + 1

        enddo
        enddo

      enddo

      IF((mjcol-1).ne.mvo(desrep))
     &  stop 'error 3 in BUILD_X_COUSAT_C!'
!________________________________________________________
!|
!|
!| from here the main and coupling block are written to
!| the external file. The SAT entries are also written to
!| the same file.
!| ATT**  we can release, and should release *all*
!| memory which we do not need anymore.
!| Due to the large size of the matrices we can not afford
!| a four-indexed lookup table as it was used in the two-particle
!| propagator!
!|
!|_______________________________________________________


      CALL XSATMAKE_C(IOCH,IBFP,NBUFS,INTBUF,DESREP,WTHR,
     &                WBUFR,WBUFI,WBUFJ,EPS,OOLT,VVLT,
     &                DOXEXT)
!
!  the SAT block in DESREP is now also on the file
!  and we can return in this symmetry and diagonalize...
!
!
!  flush buffer
!
      IF(IBFP.gt.0) THEN
        CALL XCOLFLSH_C(IOCH,IBFP,NBUFS,INTBUF,WBUFR,WBUFI,WBUFJ)
      ENDIF

      write(*,*) '--------------------------------------------'
      write(*,*) '---- Coupling and Satellite blocks finished.'
      write(*,*) '---- Number of required buffers:',NBUFS
      write(*,*) '--------------------------------------------'

      deallocate(bufcolex)
      deallocate(bufcol)
      deallocate(vovv1)
      deallocate(vooo1)
      deallocate(vvlt)
      deallocate(oolt)
!
! release storage for matrix writing
!
      deallocate(wbufj)
      deallocate(wbufi)
      deallocate(wbufr)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION XCOLDET(desrep)
!
      IMPLICIT INTEGER (A-Z)

      INTEGER                               :: XCOLDET
      INTEGER                               :: desrep
!
!---------------Description--------------------------------------------
!
!  determines length of excitation coupling (2p2h) block (V>V O>O: desrep)
!  Att: These are nontotally symmetric <cd||kl> streams where
!  c>d and k>l!
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Executable code--------------------------------------
!
      ICNT = 0

      DO KLREP = 1,NREP     ! klrep, cdrep and desrep are bosonic
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                  ICNT = ICNT + 1
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      XCOLDET = ICNT
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XKENEDIAG(desrep,eps)

      use qstack
!
      IMPLICIT INTEGER (A-Z)

      INTEGER                               :: desrep
      Real*8, dimension(:)                  :: eps
!
!---------------Description--------------------------------------------
!
!  routine calculates the 2p2h K_ab,ij diagonal value for strict ADC-2
!  and stores it on the stack. This considerably accelerates the
!  ADC-2s calculations since no 8-fold loop is to be executed in the
!  XSATMAKE routine.
!
!---------------Local variables--------------------------------------
!
      Real*8, allocatable,dimension (:)     :: ecdkl
      Real*8                                :: ec,ed,ek,el
      Integer                               :: icntrl
      Integer                               :: qst_eneline
!
!---------------interface region --------------------------------------
!
      interface

        INTEGER FUNCTION XCOLDET(i1)
          integer                      :: i1
        END FUNCTION

      end interface
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
#include "polprp_stacklines.h"
!
!---------------Executable code--------------------------------------
!
      qst_eneline = ENE_STACKLINE

      icntrl = xcoldet(desrep)

      allocate(ecdkl(icntrl))
      ecdkl = 0.0d0

      ICNT = 0
      DO KLREP = 1,NREP     ! klrep, cdrep and desrep are bosonic
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 30 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 30
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 40 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 40
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                  ICNT = ICNT + 1
                  coff=io(nrep+1) + iv(crep) + c
                  doff=io(nrep+1) + iv(drep) + d
                  koff=io(krep) + k
                  loff=io(lrep) + l
                  ec=eps(coff)
                  ed=eps(doff)
                  ek=eps(koff)
                  el=eps(loff)
                  ecdkl(icnt) = ec + ed - ek - el
                ENDDO
              ENDDO
 40         CONTINUE
         ENDDO
      ENDDO
 30   CONTINUE
      ENDDO

      if(icnt.ne.icntrl) stop 'Internal error XKENEDIAG!'

      if(qstack_push(qst_eneline,icnt,ecdkl).ne.icnt) stop 'QE'

      deallocate(ecdkl)

      write(*,*) 'ADC-2s sat energies written. Length check passed.'
 
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLMAKE (desrep,lencouplbl,bufcol,
     &                     jrep,brep,J,B,
     &                     mxno,mxnv,
     &                     voooz,vovvz,
     &                     oolt,vvlt)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Build the coupling block column for which we provide
! jrep,brep,J and B. BUFCOL was allocated by the caller
! and is filled in this subroutine.
!              *** real case ***
!
      INTEGER                          :: desrep,lencouplbl
      REAL*8, dimension(:)             :: bufcol
      INTEGER                          :: jrep,brep,J,B
      INTEGER                          :: mxno,mxnv
      REAL*8, dimension(:)             :: voooz,vovvz
      INTEGER, dimension(:,:,:)        :: oolt,vvlt
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Executable code--------------------------------------
!

!     write(*,*) 'Xcolmake:',desrep,lencouplbl,jrep,brep,j,b
!     write(*,*) 'Sequence of beta(cd) and beta(kl) for',desrep
!     write(*,*)
!     write(*,*) 'beta(cd)     beta(kl)'
      CDKL = 0
!
!  do the |>|> side...
!
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      if(MULTB(CDREP+NREP,KLREP+NREP,2).ne.desrep) then
         Stop 'xcolmake error #aa1'
      endif
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
!
!  do the <|<| side...
!
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!
! from here we have all ireps/spinor numbers.
! Now treat the four delta cases.
! CDKL is the local column index in the current irrep!
! The offsets into the integral fields are *calculated*
! since a lookup table is impossible.
!____________________________________________________________
!|
!|
!|  Case A:  -< dj||kl > delta_bc, we calc a VOOO-offset.
!|
      if(BREP.eq.CREP.and.B.eq.C) then
        DJREP = MULTB(DREP,JREP,1)
        IF(MULTB(DJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #44'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
!       write(*,*) 'Case A:',djrep,klrep,nvo(djrep),nvo(klrep)
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(DREP,JREP) + (J-1)*NV(DREP) + D
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1) 
     &     stop 'intern. error #66!'
        bufcol(cdkl) = bufcol(cdkl) - voooz(ioff)
      endif
!|
!|  Case B:  +< cj||kl > delta_bd, we calc a VOOO-offset.
!|
      if(BREP.eq.DREP.and.B.eq.D) then
        CJREP = MULTB(CREP,JREP,1)
        IF(MULTB(CJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #45'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
!       write(*,*) 'Case B:',cjrep,klrep,nvo(cjrep),nvo(klrep)
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(CREP,JREP) + (J-1)*NV(CREP) + C
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #67!'
        bufcol(cdkl) = bufcol(cdkl) + voooz(ioff)
      endif
!|
!|  Case C:  -< bl||cd >* delta_jk, we calc a VOVV offset
!|
      if(JREP.eq.KREP.and.J.eq.K) then
        BLREP = MULTB(BREP,LREP,1)
        IF(MULTB(BLREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #46'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
!       write(*,*) 'Case C:',blrep,cdrep,nvo(blrep),nvo(cdrep)
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,LREP) + (L-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #68!'
        bufcol(cdkl) = bufcol(cdkl) - vovvz(ioff)
      endif
!|
!|  Case D:  +< bk||cd >* delta_jl, we calc a VOVV offset
!|
      if(JREP.eq.LREP.and.J.eq.L) then
        BKREP = MULTB(BREP,KREP,1)
        IF(MULTB(BKREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #47'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
!       write(*,*) 'Case D:',bkrep,cdrep,nvo(bkrep),nvo(cdrep)
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,KREP) + (K-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #69!'
        bufcol(cdkl) = bufcol(cdkl) + vovvz(ioff)
      endif
!|
!|            END   Case A-D
!|___________________________________________________________
!
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(CDKL.ne.LENCOUPLBL) STOP 'Internal error #77!'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLMAKE_C (desrep,lencouplbl,bufcol,
     &                       jrep,brep,J,B,
     &                       mxno,mxnv,
     &                       voooz,vovvz,oolt,vvlt)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Build the coupling block column for which we provide
! jrep,brep,J and B. BUFCOL was allocated by the caller
! and is filled in this subroutine.
!              *** real case ***
! all arrays are transferred as ** assumed shape **
! the size in each dimension stems from the caller
!
      INTEGER                          :: desrep
      INTEGER                          :: lencouplbl
      complex*16, dimension(:)         :: bufcol
      INTEGER                          :: jrep,brep,J,B
      INTEGER                          :: mxno,mxnv
      REAL*8, dimension(:)             :: voooz
      REAL*8, dimension(:)             :: vovvz
      INTEGER, dimension(:,:,:)        :: oolt
      INTEGER, dimension(:,:,:)        :: vvlt
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Executable code--------------------------------------
!
!
! check consistency
!
      if(rcw.ne.2) then
        call quit('Internal error in XCOLMAKE_C !')
      endif

! start

      CDKL = 0
!
!  do the |>|> side...
!
      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
!
!  do the <|<| side...
!
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!
! from here we have all ireps and spinor numbers for filling the
! corresponding matrix element at the correct position.
! we have to abandon indentation here in order to fit some code
! on the line !
!
! treat the four delta cases.
! CDKL can be seen as the local column index here in this irrep!
! the offsets are *calculated* and not yielded via a lookup table
! again.
!____________________________________________________________
!|
!|
!|
!|  Case A:  < dj||kl > delta_bc, we calc a VOOO-offset.
!|
      if(BREP.eq.CREP.and.B.eq.C) then
        DJREP = MULTB(DREP,JREP,1)
        IF(MULTB(DJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #44'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(DREP,JREP) + (J-1)*NV(DREP) + D
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1) 
     &     stop 'intern. error #66!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) - dcmplx(voooz(iir),voooz(iir+1))
      endif
!|
!|  Case B:  < cj||kl > delta_bd, we calc a VOOO-offset.
!|
      if(BREP.eq.DREP.and.B.eq.D) then
        CJREP = MULTB(CREP,JREP,1)
        IF(MULTB(CJREP+NREP,KLREP+NREP,2).ne.1) STOP 'IE #45'
        IF(KREP.EQ.LREP) THEN
          IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
        ELSE
          IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
        ENDIF
        IOFF=IOFF * NVO(KLREP) + IVOOOT(KLREP)
        IOFF=IOFF + IIVO(CREP,JREP) + (J-1)*NV(CREP) + C
        if(ioff.gt.ivooot(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #67!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) + dcmplx(voooz(iir),voooz(iir+1))
      endif
!|
!|  Case C:  < bl||cd >* delta_jk, we calc a VOVV offset
!|
      if(JREP.eq.KREP.and.J.eq.K) then
        BLREP = MULTB(BREP,LREP,1)
        IF(MULTB(BLREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #46'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,LREP) + (L-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #68!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) - 
     &    dconjg(dcmplx(vovvz(iir),vovvz(iir+1)))
      endif
!|
!|  Case D:  < bk||cd >* delta_jl, we calc a VOVV offset
!|
      if(JREP.eq.LREP.and.J.eq.L) then
        BKREP = MULTB(BREP,KREP,1)
        IF(MULTB(BKREP+NREP,CDREP+NREP,2).ne.1) STOP 'IE #47'
        IF(CREP.EQ.DREP) THEN
          IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
        ELSE
          IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
        ENDIF
        IOFF=IOFF * NVO(CDREP) + IVOVVT(CDREP)
        IOFF=IOFF + IIVO(BREP,KREP) + (K-1)*NV(BREP) + B
        if(ioff.gt.ivovvt(nrep+1).or.ioff.lt.1)
     &     stop 'intern. error #69!'
        iir = (ioff - 1)*rcw + 1
        bufcol(cdkl) = bufcol(cdkl) + 
     &    dconjg(dcmplx(vovvz(iir),vovvz(iir+1)))
      endif
!|
!|
!|___________________________________________________________
!
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(CDKL.ne.LENCOUPLBL) STOP 'Internal error #77!'
!
!     write(iw,*) 'Xcolmake: current run finished.'
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLWRIT (IOCH,wthr,bufferin,buf,ioi,ioj,
     &                     IB,NBUFS,LENBUFFER,mjcol,INTBUF)
!
      IMPLICIT NONE

      integer                       ::  ioch
      real*8                        ::  wthr
      real*8, dimension(:)          ::  bufferin,buf
      integer, dimension(:)         ::  ioi,ioj
      integer                       ::  ib,nbufs,lenbuffer,mjcol,intbuf
!
!---------------Description--------------------------------------------
!
!  writes out a complete column of the ADC matrix held in bufferin
!  to the buckets provided by the caller. If the buckets are full,
!  they are emptied to the file in IOCH.
!
!---------------Local variables ---------------------------------------
!
      INTEGER            ::  i,jdummy,ixx
      REAL*8             ::  s
!
!---------------Execution ---------------------------------------------
!
      jdummy = 0
      DO i=mjcol,lenbuffer  !i hereby serves as the row index !
        s=bufferin(i)
        IF(abs(s).gt.wthr) THEN
          IB=IB+1
          BUF(IB)=S
          IOI(IB)=i
          IOJ(IB)=mjcol
          IF(IB.EQ.INTBUF) THEN
            NBUFS = NBUFS + 1
            WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                  (IOI(IXX),IXX=1,INTBUF),
     &                  (IOJ(IXX),IXX=1,INTBUF),
     &                   INTBUF,JDUMMY
            IB = 0
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLWRIT_C (ioch,wthr,bufferin,buf,ioi,ioj,
     &                       ib,nbufs,lenbuffer,mjcol,intbuf)
!
      IMPLICIT NONE

      integer                       ::  ioch
      real*8                        ::  wthr
      complex*16, dimension(:)      ::  bufferin,buf
      integer, dimension(:)         ::  ioi,ioj
      integer                       ::  ib,nbufs,lenbuffer,mjcol,intbuf
!
!---------------Description--------------------------------------------
!
!  writes out a complete column of the ADC matrix held in bufferin
!  to the buckets provided by the caller. If the buckets are full,
!  they are emptied to the file in IOCH.
!
!---------------Local variables ---------------------------------------
!
      INTEGER            ::  i,jdummy,ixx
      complex*16         ::  s
!
!---------------Execution ---------------------------------------------
!
      jdummy = 0
      DO i=mjcol,lenbuffer  !i hereby serves as the row index starting from mjcol 
        s=bufferin(i)
        IF(abs(s).gt.wthr) THEN
          IB=IB+1
          BUF(IB)=S
          IOI(IB)=i
          IOJ(IB)=mjcol
          IF(IB.EQ.INTBUF) THEN
            NBUFS = NBUFS + 1
            WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                  (IOI(IXX),IXX=1,INTBUF),
     &                  (IOJ(IXX),IXX=1,INTBUF),
     &                   INTBUF,JDUMMY
            IB = 0
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLFLSH (IOCH,IB,NBUFS,INTBUF,buf,ioi,ioj)
!
      IMPLICIT NONE

      integer                       :: IOCH,IB,NBUFS,INTBUF
      real*8, dimension(:)          :: buf
      integer, dimension(:)         :: ioi,ioj

      integer                       :: ixx,jdummy = 0

      NBUFS = NBUFS + 1
      WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &            (IOI(IXX),IXX=1,INTBUF),
     &            (IOJ(IXX),IXX=1,INTBUF),
     &            IB,JDUMMY
      IB = 0

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XCOLFLSH_C (IOCH,IB,NBUFS,INTBUF,buf,ioi,ioj)
!
      IMPLICIT NONE

      integer                       :: IOCH,IB,NBUFS,INTBUF
      complex*16, dimension(:)      :: buf
      integer, dimension(:)         :: ioi,ioj

      integer                       :: ixx,jdummy = 0

      NBUFS = NBUFS + 1
      WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &            (IOI(IXX),IXX=1,INTBUF),
     &            (IOJ(IXX),IXX=1,INTBUF),
     &            IB,JDUMMY
      IB = 0

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSATMAKE(ioch,ibfp,nbufs,intbuf,desrep,wthr,
     &                    buf,ioi,ioj,eps,oolt,vvlt,dox)

      use qstack
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!  generates the ADC-2x SAT block.
!  the column counter MJCOL carries the actual value at this stage!
!
!
!---------------Calling variables--------------------------------------
!
      INTEGER                   :: ioch,ibfp,nbufs,intbuf
      REAL*8                    :: wthr
      INTEGER                   :: mjcol,desrep
      REAL*8, dimension(:)      :: buf
      INTEGER, dimension(:)     :: ioi,ioj
      REAL*8, dimension(:)      :: eps
      INTEGER,dimension(:,:,:)  :: oolt, vvlt
      LOGICAL                   :: dox  !tells the routine ADC2-x
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "polprp_stacklines.h"
!
!---------------Local variables--------------------------------------
!
      REAL*8, allocatable, dimension(:)    :: vvvvz, ooooz, vovoz
      INTEGER                              :: n1,ialloc
      INTEGER*8                            :: i8,n8
      REAL*8                               :: EA,EB,EI,EJ
      INTEGER                              :: JDUMMY

      real*8                               :: S

      Integer                              :: kenelen
      REAL*8, allocatable, dimension(:)    :: eabij

      Integer                              :: qst_eneline
!
!---------------Interface area --------------------------------------
!
      interface

        SUBROUTINE mpi_master_vvvv_complete(ra1,i1)
          real*8, dimension(:)            :: ra1
          integer*8                       :: i1
        END SUBROUTINE

        INTEGER FUNCTION XCOLDET(i1)
          integer                      :: i1
        END FUNCTION

      end interface
!
!---------------Start execution--------------------------------------
!
      call pst('Constructing real satellite part+')

      qst_eneline = ENE_STACKLINE
!
! fetch diagonal energies of the satellite block from stack
! (always done since it was pushed by XKENEDIAG)
!
      kenelen = xcoldet(desrep)
      allocate(eabij(kenelen))
      eabij = 0.0d0
      if(qstack_popf(qst_eneline,eabij).ne.kenelen) stop 'QE'

!______________________________________________________
!|
!|  if we do a restricted calculation only we just add
!|  the energy diagonal to the write buffer and we are done.
!|
      IF(.not.DOX) THEN

        write(*,*) 'for restricted ADC-2'
        Do iene = 1,kenelen

          S = eabij(iene)

          IF(ABS(S).gt.WTHR) THEN
            IBFP=IBFP+1
            BUF(IBFP)=S
            IOI(IBFP)=MVO(DESREP) + iene
            IOJ(IBFP)=MVO(DESREP) + iene
            IF(IBFP.EQ.INTBUF) THEN
              NBUFS = NBUFS + 1
              WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                    (IOI(IXX),IXX=1,INTBUF),
     &                    (IOJ(IXX),IXX=1,INTBUF),
     &                     INTBUF,JDUMMY
              IBFP = 0
            ENDIF
          ENDIF
        Enddo ! iene

        write(*,*) 'energy diagonal written.'
        deallocate(eabij)
        return
      ENDIF  ! .not.DOX
!|
!|_____________________________________________________

!_______________  extended ADC 2 from here  ___________
!|
!|  Allocate required arrays for extended mode.
!|  VVVV, OOOO and VOVO integrals are not needed in strict mode.
!|
!|  The if(dox) can go away because from here dox is set.
!|

      write(*,*) 'for extended ADC-2'
      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VVVV buffer could not be allocated.'

#if defined (VAR_MPI)
      do i8=1,n8
        vvvvz(i8) = -9999.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      call mpi_master_vvvv_complete(vvvvz,n8)  !.... and completes the stream.
      write(*,*) 'VVVV check...'
      do i8=1,n8
        if(vvvvz(i8).eq. -9999.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) '... complete.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif

      n1 = ioooott(nrep+1)
      allocate(ooooz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'OOOO buffer could not be allocated.'

      n1 = ivovo(nrep+1)
      allocate(vovoz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVO buffer could not be allocated.'
 
      CALL GETOOOO(ooooz)
      CALL GETVOVO(vovoz)
!
! SAT elements are calculated only for the lower triangle which
! means that ABIJ >= CDKL !
! no lookup table (too big!)
!
      CDKL=0    !column counter

      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the cdkl loop.
!  we abandon indentation here and start the abij loop.
!  This part of code unfortunately scales as N**8 but we have
!  no remedy at the moment.

      ABIJ=0     !row counter

      DO IJREP = 1,NREP
      ABREP = MULTB(DESREP+NREP,IJREP+NREP,2)
      DO 44 JREP = 1, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 44
      DO J = 1, NO(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, NO(IREP)
            DO 33 BREP = 1, NREP
              AREP = MULTB(BREP,ABREP+NREP,2)
              IF (AREP.LT.BREP) GOTO 33
              DO B = 1, NV(BREP)
                AMIN = 1
                IF (AREP.EQ.BREP) AMIN = B + 1
                DO A = AMIN, NV(AREP)
                   ABIJ = ABIJ + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the abij loop.
!  Again, we abandon indentation here and start the delta cases
!
      IF(ABIJ.GE.CDKL) THEN
!
! we are on or below the diagonal and start calculation
!
! otherwise: DO NOTHING!
!
!  a: check diagonal
!
       S = 0.0d0

       if(arep.eq.crep.and.brep.eq.drep.and.
     &    irep.eq.krep.and.jrep.eq.lrep) then
         if(a.eq.c.and.b.eq.d.and.i.eq.k.and.j.eq.l) then
           aoff=io(nrep+1) + iv(arep) + a
           boff=io(nrep+1) + iv(brep) + b
           ioff=io(irep) + i
           joff=io(jrep) + j
           ea=eps(aoff); eb=eps(boff); ei=eps(ioff); ej=eps(joff)
           S = S + (ea + eb - ei - ej)
         endif
       endif
!
!  b: check VVVV contributions  <AB||CD> delta_ik delta_jl
!
       if(irep.eq.krep.and.jrep.eq.lrep) then
         if(i.eq.k.and.j.eq.l) then
           ZCDREP=MULTB(CREP,DREP,1)
           IF(ZCDREP.ne.CDREP) STOP 'internal error 1 XSATMAKE'
           IF(CREP.EQ.DREP) THEN
             IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
           ELSE
             IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
           ENDIF
           IOFF=IOFF * NVVT(ZCDREP) + IVVVVTT(ZCDREP)
           IF(AREP.EQ.BREP) THEN
             IOFF=IOFF + IIVVT(AREP,BREP) + VVLT(A,B,AREP)
           ELSE
             IOFF=IOFF + IIVVT(AREP,BREP) +
     &         (B-1)*NV(AREP) + A
           ENDIF
           if(ioff.gt.ivvvvtt(nrep+1)) stop 'vvvv offset error !'
           S = S + vvvvz(IOFF)
         endif
       endif
!
!  c: check OOOO contributions  <IJ||KL>* delta_ac delta_bd
!
       if(arep.eq.crep.and.brep.eq.drep) then
         if(a.eq.c.and.b.eq.d) then
           ZKLREP=MULTB(KREP,LREP,1)
           IF(ZKLREP.ne.KLREP) STOP 'internal error 2 XSATMAKE'
           IF(KREP.EQ.LREP) THEN
             IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
           ELSE
             IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
           ENDIF
           IOFF=IOFF * NOOT(ZKLREP) + IOOOOTT(ZKLREP)
           IF(IREP.EQ.JREP) THEN
             IOFF=IOFF + IIOOT(IREP,JREP) + OOLT(I,J,IREP)
           ELSE
             IOFF=IOFF + IIOOT(IREP,JREP) +
     &         (J-1)*NO(IREP) + I
           ENDIF
           if(ioff.gt.ioooott(nrep+1)) stop 'oooo offset error !'
           S = S + ooooz(IOFF)
         endif
       endif
!
! Block 1: check the four delta_bd  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.drep.and.b.eq.d)then

! 1a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 1b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 1c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 1d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 

       endif  ! delta_bd
!
! Block 2: check the four delta_ac  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(arep.eq.crep.and.a.eq.c)then
! 2a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 2b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 2c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 2d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 

       endif  !delta_ac
!
! Block 3: check the four delta_bc  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.crep.and.b.eq.c)then
! 3a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 3b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 3c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 3d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 

       endif  !delta_bc
!
! Block 4: check the four delta_ad  contributions
!
       if(arep.eq.drep.and.a.eq.d)then
! 4a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 4b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S + vovoz(IOFF)
         endif 
! 4c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 
! 4d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           S = S - vovoz(IOFF)
         endif 

       endif  !delta_ad
!
!  write out corresponding matrix entry with row and column number
!  note: local row and column counters both start at value 1.
!  row and column counter has to be offset by size of main block
!
       IF(ABS(S).gt.WTHR) THEN
         IBFP=IBFP+1
         BUF(IBFP)=S
         IOI(IBFP)=ABIJ + MVO(DESREP)
         IOJ(IBFP)=CDKL + MVO(DESREP)
         IF(IBFP.EQ.INTBUF) THEN
           NBUFS = NBUFS + 1
           WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                  INTBUF,JDUMMY
           IBFP = 0
         ENDIF
       ENDIF

      ENDIF     !(ABIJ.GE.CDKL)
!  
!
!_______________________________________________________________
                ENDDO
              ENDDO
 33         CONTINUE
         ENDDO
      ENDDO
 44   CONTINUE
      ENDDO
!_______________________________________________________________
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO  !KLREP

! we deallocate only in extended mode
! otherwise the arrays are not allocated.

      deallocate(vovoz)
      deallocate(ooooz)
      deallocate(vvvvz)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XSATMAKE_C (ioch,ibfp,nbufs,intbuf,desrep,wthr,
     &                       buf,ioi,ioj,eps,oolt,vvlt,dox)

      use qstack
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!  generates the complex ADC-2x SAT block.
!  the column counter MJCOL carries the actual value at this stage!
!
!
!---------------Calling variables--------------------------------------
!
      INTEGER                    :: ioch,ibfp,nbufs,intbuf
      REAL*8                     :: wthr
      INTEGER                    :: mjcol,desrep
      complex*16, dimension(:)   :: buf
      INTEGER, dimension(:)      :: ioi,ioj
      REAL*8, dimension(:)       :: eps
      INTEGER,dimension(:,:,:)   :: oolt, vvlt
      LOGICAL                    :: dox  !tells the routine ADC2-x
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "polprp_stacklines.h"
!
!---------------Local variables--------------------------------------
!
      REAL*8, allocatable, dimension(:)    :: vvvvz, ooooz, vovoz
      INTEGER                              :: n1,ialloc
      INTEGER*8                            :: i8,n8
      REAL*8                               :: EA,EB,EI,EJ
      INTEGER                              :: JDUMMY

      complex*16                           :: S

      Integer                              :: kenelen
      REAL*8, allocatable, dimension(:)    :: eabij

      Integer                              :: qst_eneline
!
!---------------Interface area --------------------------------------
!
      INTERFACE

        SUBROUTINE mpi_master_vvvv_complete(ra1,i1)
          real*8, dimension(:)            :: ra1
          integer*8                       :: i1
        END SUBROUTINE

        INTEGER FUNCTION XCOLDET(i1)
          integer                      :: i1
        END FUNCTION

      END INTERFACE
!
!---------------Start execution--------------------------------------
!
      call pst('Constructing complex satellite part+')

      qst_eneline = ENE_STACKLINE
!
! fetch diagonal energies of the satellite block from stack
! (always done since it was pushed by XKENEDIAG)
!
      kenelen = xcoldet(desrep)
      allocate(eabij(kenelen))
      eabij = 0.0d0
      if(qstack_popf(qst_eneline,eabij).ne.kenelen) stop 'QE'
!______________________________________________________
!|
!|  if we do a restricted calculation only we just add
!|  the energy diagonal to the write buffer and we are done.
!|
      IF(.not.DOX) THEN

        write(*,*) 'for restricted ADC-2'
        Do iene = 1,kenelen

          S = dcmplx(eabij(iene),0.0d0)

          IF(ABS(S).gt.WTHR) THEN
            IBFP=IBFP+1
            BUF(IBFP)=S
            IOI(IBFP)=MVO(DESREP) + iene
            IOJ(IBFP)=MVO(DESREP) + iene
            IF(IBFP.EQ.INTBUF) THEN
              NBUFS = NBUFS + 1
              WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                    (IOI(IXX),IXX=1,INTBUF),
     &                    (IOJ(IXX),IXX=1,INTBUF),
     &                     INTBUF,JDUMMY
              IBFP = 0
            ENDIF
          ENDIF
        Enddo ! iene

        write(*,*) 'energy diagonal written.'
        deallocate(eabij)
        return
      ENDIF  ! .not.DOX
!|
!|_____________________________________________________

!_______________  extended ADC 2 from here  ___________
!|
!|  we allocate these arrays only in extended mode.
!|  VVVV, OOOO and VOVO integrals are not needed in strict mode.
!|
!|  The if(dox) can go away because from here dox is set.
!|
      write(*,*) 'for extended ADC-2'
      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VVVV buffer could not be allocated.'

#if defined (VAR_MPI)
      do i8=1,n8*rcw
        vvvvz(i8) = -9999.0d0
      enddo
      call rdvvvv(vvvvz)  ! master reads its own chunks
      call mpi_master_vvvv_complete(vvvvz,n8)  !.... and completes the stream.
      write(*,*) 'VVVV check...'
      do i8=1,n8*rcw
        if(vvvvz(i8).eq. -9999.0d0)
     &    call quit('Internal error: gap in the VVVV stream!')
      enddo
      write(*,*) '... complete.'
#else
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
#endif

      n1 = ioooott(nrep+1)*rcw
      allocate(ooooz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'OOOO buffer could not be allocated.'

      n1 = ivovo(nrep+1)*rcw
      allocate(vovoz(n1),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVO buffer could not be allocated.'

      CALL GETOOOO(ooooz)
      CALL GETVOVO(vovoz)

!
! SAT elements are calculated only for the lower triangle which
! means that ABIJ >= CDKL !
! no lookup table (too big!)
!
      CDKL=0    !column counter

      DO KLREP = 1,NREP     ! KLREP (bosonic!) loops through
      CDREP = MULTB(DESREP+NREP,KLREP+NREP,2)
      DO 10 LREP = 1, NREP
      KREP = MULTB(LREP,KLREP+NREP,2)
      IF (KREP.LT.LREP) GOTO 10
      DO L = 1, NO(LREP)
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, NO(KREP)
            DO 20 DREP = 1, NREP
              CREP = MULTB(DREP,CDREP+NREP,2)
              IF (CREP.LT.DREP) GOTO 20
              DO D = 1, NV(DREP)
                CMIN = 1
                IF (CREP.EQ.DREP) CMIN = D + 1
                DO C = CMIN, NV(CREP)
                   CDKL = CDKL + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the cdkl loop.
!  we abandon indentation here and start the abij loop.
!  This part of code unfortunately scales as N**8 but we have
!  no remedy at the moment.

      ABIJ=0     !row counter

      DO IJREP = 1,NREP
      ABREP = MULTB(DESREP+NREP,IJREP+NREP,2)
      DO 44 JREP = 1, NREP
      IREP = MULTB(JREP,IJREP+NREP,2)
      IF (IREP.LT.JREP) GOTO 44
      DO J = 1, NO(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, NO(IREP)
            DO 33 BREP = 1, NREP
              AREP = MULTB(BREP,ABREP+NREP,2)
              IF (AREP.LT.BREP) GOTO 33
              DO B = 1, NV(BREP)
                AMIN = 1
                IF (AREP.EQ.BREP) AMIN = B + 1
                DO A = AMIN, NV(AREP)
                   ABIJ = ABIJ + 1
!_______________________________________________________________
!
!  here we are at the innermost position of the abij loop.
!  Again, we abandon indentation here and start the delta cases
!
      IF(ABIJ.GE.CDKL) THEN
!
! we are on or below the diagonal and start selection
!
! otherwise: DO NOTHING!
!
!  a: check diagonal
!
       S = (0.0d0,0.0d0)

       if(arep.eq.crep.and.brep.eq.drep.and.
     &    irep.eq.krep.and.jrep.eq.lrep) then
         if(a.eq.c.and.b.eq.d.and.i.eq.k.and.j.eq.l) then
           aoff=io(nrep+1) + iv(arep) + a
           boff=io(nrep+1) + iv(brep) + b
           ioff=io(irep) + i
           joff=io(jrep) + j
           ea=eps(aoff); eb=eps(boff); ei=eps(ioff); ej=eps(joff)
           S = S + dcmplx((ea + eb - ei - ej),0.0d0)
         endif
       endif
!
!  b: check VVVV contributions  <AB||CD> delta_ik delta_jl
!
       if(irep.eq.krep.and.jrep.eq.lrep) then
         if(i.eq.k.and.j.eq.l) then
           ZCDREP=MULTB(CREP,DREP,1)
           IF(ZCDREP.ne.CDREP) STOP 'internal error 1 XSATMAKE'
           IF(CREP.EQ.DREP) THEN
             IOFF=IIVVT(CREP,DREP) + VVLT(C,D,CREP) - 1
           ELSE
             IOFF=IIVVT(CREP,DREP) + (D-1)*NV(CREP) + C - 1
           ENDIF
           IOFF=IOFF * NVVT(ZCDREP) + IVVVVTT(ZCDREP)
           IF(AREP.EQ.BREP) THEN
             IOFF=IOFF + IIVVT(AREP,BREP) + VVLT(A,B,AREP)
           ELSE
             IOFF=IOFF + IIVVT(AREP,BREP) +
     &         (B-1)*NV(AREP) + A
           ENDIF
           if(ioff.gt.ivvvvtt(nrep+1)) stop 'vvvv offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vvvvz(iir),vvvvz(iir+1))
         endif
       endif
!
!  c: check OOOO contributions  <IJ||KL>* delta_ac delta_bd
!
       if(arep.eq.crep.and.brep.eq.drep) then
         if(a.eq.c.and.b.eq.d) then
           ZKLREP=MULTB(KREP,LREP,1)
           IF(ZKLREP.ne.KLREP) STOP 'internal error 2 XSATMAKE'
           IF(KREP.EQ.LREP) THEN
             IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1
           ELSE
             IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1
           ENDIF
           IOFF=IOFF * NOOT(ZKLREP) + IOOOOTT(ZKLREP)
           IF(IREP.EQ.JREP) THEN
             IOFF=IOFF + IIOOT(IREP,JREP) + OOLT(I,J,IREP)
           ELSE
             IOFF=IOFF + IIOOT(IREP,JREP) +
     &         (J-1)*NO(IREP) + I
           ENDIF
           if(ioff.gt.ioooott(nrep+1)) stop 'oooo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dconjg(dcmplx(ooooz(iir),ooooz(iir+1)))
         endif
       endif
!
! Block 1: check the four delta_bd  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.drep.and.b.eq.d)then

! 1a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 1b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 1c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 1d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  ! delta_bd
!
! Block 2: check the four delta_ac  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(arep.eq.crep.and.a.eq.c)then
! 2a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 2b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 2c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 2d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  !delta_ac
!
! Block 3: check the four delta_bc  contributions (see written
!      implementation notes, whose convention is followed !)
!
       if(brep.eq.crep.and.b.eq.c)then
! 3a
         if(jrep.eq.lrep.and.j.eq.l) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 3b
         if(irep.eq.krep.and.i.eq.k) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 3c
         if(jrep.eq.krep.and.j.eq.k) then
           DIREP=MULTB(DREP,IREP,1)
           IOFF=IIVO(DREP,IREP) + (I-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DIREP) + IVOVO(DIREP)
           IOFF=IOFF + IIVO(AREP,LREP) +
     &         (L-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 3d
         if(irep.eq.lrep.and.i.eq.l) then
           DJREP=MULTB(DREP,JREP,1)
           IOFF=IIVO(DREP,JREP) + (J-1)*NV(DREP) + D - 1
           IOFF=IOFF*NVO(DJREP) + IVOVO(DJREP)
           IOFF=IOFF + IIVO(AREP,KREP) +
     &         (K-1)*NV(AREP) + A
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  !delta_bc
!
! Block 4: check the four delta_ad  contributions
!
       if(arep.eq.drep.and.a.eq.d)then
! 4a
         if(jrep.eq.lrep.and.j.eq.l) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 4b
         if(irep.eq.krep.and.i.eq.k) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S + dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 4c
         if(jrep.eq.krep.and.j.eq.k) then
           CIREP=MULTB(CREP,IREP,1)
           IOFF=IIVO(CREP,IREP) + (I-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CIREP) + IVOVO(CIREP)
           IOFF=IOFF + IIVO(BREP,LREP) +
     &         (L-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 
! 4d
         if(irep.eq.lrep.and.i.eq.l) then
           CJREP=MULTB(CREP,JREP,1)
           IOFF=IIVO(CREP,JREP) + (J-1)*NV(CREP) + C - 1
           IOFF=IOFF*NVO(CJREP) + IVOVO(CJREP)
           IOFF=IOFF + IIVO(BREP,KREP) +
     &         (K-1)*NV(BREP) + B
           if(ioff.gt.ivovo(nrep+1)) stop 'vovo offset error !'
           iir = (ioff - 1)*rcw + 1
           S = S - dcmplx(vovoz(iir),vovoz(iir+1))
         endif 

       endif  !delta_ad
!
!  write out corresponding matrix entry with row and column number
!  note: local row and column counters both start at value 1.
!  row and column counter has to be offset by size of main block
!
       IF(ABS(S).gt.WTHR) THEN
         IBFP=IBFP+1
         BUF(IBFP)=S
         IOI(IBFP)=ABIJ + MVO(DESREP)
         IOJ(IBFP)=CDKL + MVO(DESREP)
         IF(IBFP.EQ.INTBUF) THEN
           NBUFS = NBUFS + 1
           WRITE(IOCH) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                  INTBUF,JDUMMY
           IBFP = 0
         ENDIF
       ENDIF

      ENDIF     !(ABIJ.GE.CDKL)
!  
!
!_______________________________________________________________
                ENDDO
              ENDDO
 33         CONTINUE
         ENDDO
      ENDDO
 44   CONTINUE
      ENDDO
!_______________________________________________________________
                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO !KLREP

! we deallocate only in extended mode
! otherwise the arrays are not allocated.

      deallocate(vovoz)
      deallocate(ooooz)
      deallocate(vvvvz)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE print_adc_matrix(iw,iobase,intbuf,nbufs)
!
!     For small matrices only: print the whole matrix.

      IMPLICIT NONE

!-------------- formal parameters

      integer        ::  iw
      integer        ::  iobase
      integer        ::  intbuf
      integer        ::  nbufs

!-------------- common variables

!-------------- local variables

      integer                             :: i,ixx,nact,jdummy
      real*8, allocatable, dimension (:)  :: buf
      integer, allocatable, dimension (:) :: ioi,ioj
      integer                             :: icol,irow
      real*8                              :: a

!-------------- execution

      call PST('Printing the complete ADC matrix+')
      write(iw,*)
      write(iw,*)
!
!  allocate read/write buffers
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))

      do i=1,nbufs
        read(iobase,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                       (ioi(ixx),ixx=1,intbuf),
     &                       (ioj(ixx),ixx=1,intbuf),
     &                       nact,jdummy
          do ixx = 1,nact
            irow = ioi(ixx)
            icol = ioj(ixx)
            a    = buf(ixx)
            write(iw,'(2I8,E22.14)') irow,icol,a
          enddo
      enddo
      write(iw,*)
      write(iw,*) '---------------------------------------------'

      deallocate(ioi,ioj,buf)
      return

 889  call quit('Error in reading ADC matrix')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
