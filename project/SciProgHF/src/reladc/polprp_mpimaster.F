#if defined (VAR_MPI)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_master_vvvv_complete(vvvv,idi)
!
      use interface_to_mpi
!
      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!  routine for gathering the VVVV integrals from the nodes.
!
!       ****** executed by MASTER ******
!
!  Attention: Also here we follow the CCSETI integral convention that
!  the same real array is used for real and complex integrals but is twice
!  the size in case of complex arithmetic.
!
!  Attention! the transferred vovv array *must not* be set to zero because
!  it is incrementally filled. On entry it already contains partial contributions from
!  the master node!
!
!---------------Calling variables--------------------------------------
!
      real*8, dimension(:)     :: vvvv
      integer*8                :: idi
!
!---------------Common Blocks--------------------------------------
!
#include  "polprp_servercodes.h"
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER                           :: irep,node
      INTEGER                           :: i
      INTEGER                           :: intlength,localoffset
      REAL*8,allocatable,dimension(:)   :: buf

      write(*,*) ' ---- Master: entering mpi_master_vvvv_complete.'
      write(*,*) ' ---- Master: number of ireps:',nrep

      do irep=1,nrep    ! after each irep, slaves wait for a new bcast!

        call interface_mpi_BCAST(SERVER_VVVV,1,MASTER,
     &                           global_communicator)
!
!  fetch contributions to irep from each node
!  in the length/offset numbers the RCW is already accounted for!
!
        do node=MASTER+1,NMPROC-1
          call interface_mpi_send(irep,1,node,MSGN,global_communicator)
          call interface_mpi_recv(intlength,1,node,MSGN,
     &                            global_communicator)
          call interface_mpi_recv(localoffset,1,node,MSGN,
     &                            global_communicator)
!         write(*,*) '   ++ master obtained irep/intlength/offset',
!    &        irep,intlength,localoffset,' from',node
          if(intlength.eq.0) cycle
!
!  only allocate if nonzero length was transferred
!
          allocate(buf(intlength))
          call interface_mpi_recv(buf,intlength,node,MSGN,
     &                            global_communicator)
          do i=1,intlength
            vvvv(localoffset + i) = buf(i)
          enddo
          deallocate(buf)
!         write(*,*) ' ---- Master: finished fetching from node',node
        enddo  !node
!       write(*,*) ' ---- Master: finished fetching from irep',irep
      enddo  ! irep

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_master_vovv_complete(vovv,idi)
!
      use interface_to_mpi
!
      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!  routine for gathering the VOVV integrals from the nodes.
!
!       ****** executed by MASTER ******
!
!  Attention: Also here we follow the CCSETI integral convention that
!  the same real array is used for real and complex integrals but is twice
!  the size in case of complex arithmetic.

!  Attention! the transferred vovv array *must not* be set to zero because
!  it is incrementally filled. On entry it already contains partial contributions from
!  the master node!
!
!---------------Calling variables--------------------------------------
!
      real*8, dimension(:)       :: vovv(:)
      integer*8                  :: idi
!
!---------------Common Blocks--------------------------------------
!
#include  "polprp_servercodes.h"
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER                           :: irep,node
      INTEGER                           :: i
      INTEGER                           :: intlength,localoffset
      REAL*8,allocatable,dimension(:)   :: buf

      write(*,*) ' ---- Master: entering mpi_master_vovv_complete.'
      write(*,*) ' ---- Master: number of ireps:',nrep

      do irep=1,nrep
        call interface_mpi_BCAST(SERVER_VOVV,1,MASTER,
     &                           global_communicator)
!
!  fetch contributions to irep from each node
!  in the length/offset numbers the RCW is already accounted for!
!
        do node=MASTER+1,NMPROC-1
          call interface_mpi_send(irep,1,node,MSGN,global_communicator)
          call interface_mpi_recv(intlength,1,node,MSGN,
     &                            global_communicator)
          call interface_mpi_recv(localoffset,1,node,MSGN,
     &                            global_communicator)
!         write(*,*) '   ++ master obtained length/offs',intlength,
!    &                localoffset,' from',node,irep
          if(intlength.eq.0) cycle
!
! only allocate if intlength greater than zero.
!
          allocate(buf(intlength))
          call interface_mpi_recv(buf,intlength,node,MSGN,
     &                            global_communicator)
          do i=1,intlength
            vovv(localoffset + i) = buf(i)
          enddo
          deallocate(buf)
!         write(*,*) ' ---- Master: finished fetching from node',node
        enddo  !node
      enddo  ! irep
!
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_bxcs_r(PHPH_BLOCK,
     &                      dimphph,
     &                      EPS,
     &                      DESREP,
     &                      LENCOUPLBL,
     &                      INTBUF,
     &                      NBUFS,
     &                      DOXEXT,
     &                      WTHR,
     &                      IOCH,
     &                      IW)
!
      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!    ************** executed by MASTER when parallel run is done
!
!  Computes *real* 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!
!---------------Calling variables--------------------------------------
!
      real*8, dimension(:)             :: PHPH_BLOCK
      integer                          :: dimphph
      real*8, dimension(*)             :: EPS
      integer                          :: DESREP
      integer                          :: LENCOUPLBL
      integer                          :: INTBUF
      integer                          :: NBUFS
      logical                          :: DOXEXT
      real*8                           :: WTHR
      integer                          :: IOCH
      integer                          :: IW
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep,i,j,k
      INTEGER*8            :: i8,n8  !needed for the ivovvt variable!
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,nmvod,istart

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

        SUBROUTINE mpi_master_vovv_complete(ra1,i8)
          REAL*8, dimension(:)             :: ra1
          INTEGER*8                        :: i8
        END SUBROUTINE

        INTEGER FUNCTION  XCOLDET(i1)
          INTEGER                         :: i1
        END FUNCTION XCOLDET

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
          REAL*8, dimension(*)      :: ra2
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
      wbufr=0.0d0; wbufi = 0; wbufj = 0
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
!  Attention*** Some arrays that have to store *huge* numbers 
!  Are defined as Integer*8 in symm.inc ! The generic sorters
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
! It is important to note that the VOVV integrals
! are distributed over the nodes. Unfortunately we need them
! in a heavy-duty random access manner for the construction of
! the ADC matrix. Therefore we complete the integrals on one
! node and avoid communicating a huge load of individual numbers
! and index computations. (less computer time, more memory needs)
!
! convention: integers denote the *number* of integrals irrespective of
! their actual length (real/complex)
! 
      n8=ivovvt(nrep+1)
      allocate(vovv1(n8),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('allocation error of vovv1 array')
      endif   
      do i8=1,n8
        vovv1(i8) = -9999.0d0
      enddo
!
! rdvovv now only reads integral portions on master node !
! and we need to complete the remaining integrals.
!
      call rdvovv (vovv1)
      call mpi_master_vovv_complete(vovv1,n8)
      do i8=1,n8
        if(vovv1(i8).eq.-9999.0d0)
     &  call quit('Internal error: illegal gaps in the vovv stream.')
      enddo
      write(*,*) 'Master: <VO||VV> stream ready.'
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
      if(ialloc.ne.0) STOP 'alloc error #355'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #356'
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2H2P/HP coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
      mjcol = 1
      do jrep = 1,nrep

        brep = multb(desrep+nrep,jrep,2)

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
        enddo
        enddo
      enddo

      IF((mjcol-1).ne.mvo(desrep))
     &  stop 'error 3 in build_x_cousat_r!'
!
!  free space no longer in use.
!
      deallocate(vovv1)
      deallocate(vooo1)
      deallocate(bufcolex)
      deallocate(bufcol)
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
      write(iw,*) '----------- after xsatmake ----------------------'
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
      write(iw,*) '----------- after xcolflsh ----------------------'

      write(*,*) '--------------------------------------------'
      write(*,*) '---- Coupling and Satellite blocks finished.'
      write(*,*) '---- Number of required buffers:',NBUFS
      write(*,*) '--------------------------------------------'
!
! release all intermediate storage
!
      deallocate(vvlt)
      deallocate(oolt)
      deallocate(wbufj)
      deallocate(wbufi)
      deallocate(wbufr)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_bxcs_c(PHPH_BLOCK,
     &                      dimphph,
     &                      EPS,
     &                      DESREP,
     &                      LENCOUPLBL,
     &                      INTBUF,
     &                      NBUFS,
     &                      DOXEXT,
     &                      WTHR,
     &                      IOCH,
     &                      IW)
!
      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!    ************** executed by MASTER when parallel run is done
!
!  Computes *complex* 2H2P/HP and SAT block of the exc. ADC Matrix
!  In the corresponding symmetry DESREP.
!
!---------------Calling variables--------------------------------------
!
      real*8, dimension(:)             :: PHPH_BLOCK
      integer                          :: dimphph
      real*8, dimension(*)             :: EPS
      integer                          :: DESREP
      integer                          :: LENCOUPLBL
      integer                          :: INTBUF
      integer                          :: NBUFS
      logical                          :: DOXEXT
      real*8                           :: WTHR
      integer                          :: IOCH
      integer                          :: IW
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      INTEGER              :: n1,irep,i,j,k,ix
      INTEGER*8            :: i8,n8  !needed for the ivovvt variable!
      INTEGER              :: jrep,brep,jfun,bfun
      INTEGER              :: ialloc,mjcol,nmvod,istart

      INTEGER              :: mxno,mxnv,ioff
      INTEGER, dimension (:,:,:), allocatable  :: oolt
      INTEGER, dimension (:,:,:), allocatable  :: vvlt

      REAL*8, dimension(:), allocatable     :: vooo1
      REAL*8, dimension(:), allocatable     :: vovv1
      complex*16, dimension(:), allocatable :: bufcol,bufcolex
!
!  variables for writing out ADC matrix
!
      complex*16,  dimension(:), allocatable :: wbufr
      INTEGER, dimension(:), allocatable     :: wbufi, wbufj
      INTEGER              :: IBFP
      CHARACTER*6          :: NAMEX = 'ADCXPM'
!
!  auxiliary complex array
!
      complex*16,  dimension(:), allocatable :: phph_blc
!
!  parameters
!
      real*8,parameter        :: rzero=0.0d0
      complex*16,parameter    :: czero=(rzero,rzero)
!
!---------------Interface area --------------------------------------
!
      interface

        SUBROUTINE mpi_master_vovv_complete(ra1,i8)
          REAL*8, dimension(:)             :: ra1
          INTEGER*8                        :: i8
        END SUBROUTINE

        INTEGER FUNCTION XCOLDET(i1)
          INTEGER                         :: i1
        END FUNCTION XCOLDET

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
          REAL*8, dimension(*)        :: ra2
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
!________________________________________________________
!|     allocate ADC MATRIX WRITE BUFFERS and counters
!|
!|
      allocate(wbufr(intbuf),stat=ialloc)      !already complex!
      if(ialloc.ne.0) STOP 'data buffer 1 could not be allocated.'
      allocate(wbufi(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 2 could not be allocated.'
      allocate(wbufj(intbuf),stat=ialloc)
      if(ialloc.ne.0) STOP 'data buffer 3 could not be allocated.'
      wbufr = czero
      wbufi = 0 
      wbufj = 0
      IBFP  = 0  ! pointer to actual buffer level
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
      if(ialloc.ne.0) STOP 'alloc error #a in mpi_bxcs_c'
      allocate(vvlt(mxnv,mxnv,nrep),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #b in mpi_bxcs_c'
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

!___________________________________________________________
!|
!|  transform phph block stored as a real*rcw array into a 
!|  true complex array for further processing.

      allocate(phph_blc(dimphph))
      phph_blc = czero
      do i=1,dimphph
        ix = 2*(i-1) + 1
        phph_blc(i)=dcmplx(phph_block(ix),phph_block(ix+1))
      enddo
!|  end transformation
!|__________________________________________
!
    
! fetch <VO||OO> integrals needed for coupling block 
! 
      n1=ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error vooo1!'
      call getvooo(vooo1)
!
! fetch <VO||VV> integrals needed for coupling block 
! It is important to note that the VOVV integrals
! are distributed over the nodes. Unfortunately we need them
! in a heavy-duty random access manner for the construction of
! the ADC matrix. Therefore we complete the integrals on one
! node and avoid communicating a huge load of individual numbers
! and index computations. (less computer time, more memory needs)
! 
      n8=ivovvt(nrep+1)
      allocate(vovv1(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) then
        call quit('allocation error of vovv1 array')
      endif   
      do i8=1,n8*rcw
        vovv1(i8) = -9999.0d0
      enddo
!
! rdvovv now only reads integral portions on master node !
! and we need to complete the remaining integrals.
!
      call rdvovv (vovv1)
      call mpi_master_vovv_complete(vovv1,n8)
      do i8=1,n8*rcw
        if(vovv1(i8).eq.-9999.0d0)
     &  call quit('Internal error: illegal gaps in the vovv stream.')
      enddo
      write(*,*) 'Master: <VO||VV> stream ready.'
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
      if(ialloc.ne.0) STOP 'alloc error #357'
      ALLOCATE(BUFCOLEX(NMVOD + lencouplbl),stat=ialloc)
      if(ialloc.ne.0) STOP 'alloc error #358'
!
!  according to the SRT16 convention, 
!  the columns also have DESREP symmetry. 
!  We loop over the columns in the same manner and construct
!  the columns of the 2H2P/HP coupling block. 
!  after construction we write out the column together with the
!  components of the main block.
!
      mjcol = 1
      do jrep = 1,nrep

        brep = multb(desrep+nrep,jrep,2)

        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)

          bufcol   = (0.0d0,0.0d0)
          bufcolex = (0.0d0,0.0d0)

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
     &  stop 'error 3 in mpi_bxcs_c!'
!
!  free space no longer in use.
!
      deallocate(vovv1)
      deallocate(vooo1)
      deallocate(bufcolex)
      deallocate(bufcol)
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
      write(iw,*) '----------- after xsatmake_c ----------------------'
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
      write(iw,*) '----------- after xcolflsh_c ----------------------'

      write(*,*) '--------------------------------------------'
      write(*,*) '---- Coupling and Satellite blocks finished.'
      write(*,*) '---- Number of required buffers:',NBUFS
      write(*,*) '--------------------------------------------'

!
! release all intermediate storage
!
      deallocate(vvlt)
      deallocate(oolt)
      deallocate(wbufj)
      deallocate(wbufi)
      deallocate(wbufr)
      deallocate(phph_blc)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_master_mat_dist(iobase,nmbase,
     &                               isym,intbuf,nbufs,ladc,iw)
!
!     Master Matrix distribution routine. 
!
!     In addition, master also extracts diagonal from the total ADC
!     matrix and writes it to a file to be retrieved later on. 
!
      use interface_to_mpi

      IMPLICIT NONE

!-------------- formal parameters

      integer        ::  iobase
      character*6    ::  nmbase
      integer        ::  isym
      integer        ::  intbuf
      integer        ::  nbufs
      integer        ::  ladc
      integer        ::  iw

!-------------- common variables

#include "../relccsd/ccpar.inc"
#include "../relccsd/complex.inc"
#include "polprp_servercodes.h"

!-------------- local variables

      character*10                            :: locfn
      integer,allocatable,dimension(:)        :: nbufs_per_node
      integer                                 :: iobase1,io_diag
      integer                                 :: i,ip,iq,ixx,nact,jdummy
      real*8, allocatable, dimension (:)      :: buf
      real*8, allocatable, dimension (:)      :: bufr, bufi
      complex*16, allocatable, dimension (:)  :: cbuf
      integer, allocatable, dimension (:)     :: ioi,ioj
      integer                                 :: icol,irow
      integer                                 :: curslave,nbufscur
      real*8, allocatable, dimension (:)      :: adcdiagonal
      complex*16,allocatable, dimension (:)   :: cadcdiagonal

      logical                                 :: isthere

      real*8,parameter        :: rzero=0.0d0
      complex*16,parameter    :: czero=(rzero,rzero)

!-------------- execution

      IF(MYPROC.NE.MASTER) THEN
        write(*,*) ' ***** Internal inconsistency *****'
        write(*,*) ' ***** in master node number  *****'
        call quit('error')
      ENDIF

      if(carith) then
        write(*,*)  'Master: Entering (complex) matrix distribution'
      else
        write(*,*)  'Master: Entering (real) matrix distribution'
      endif
      write(*,*) ' +++ Total # of nodes/loc node: ',NMPROC,MYPROC
!
!  if number of buffers is less than the number of nodes this is a clear
!  case for a serial run. The INTBUF buffer is dimensioned such that the
!  serial matrix is small enough to go through in a short period of time
!
      write(iw,*) '+++++++ BUFSIZE:',intbuf
      write(iw,*) '+++++++ NBUFS:',nbufs
      write(iw,*) '+++++++ NPROCS:',NMPROC
      write(iw,*) '+++++++ LADC:',ladc
      if(nbufs.lt.NMPROC) then
        write(iw,*) '+++++++ WARNING. Matrix too small for'
        write(iw,*) '+++++++ parallel run!. Please do a serial'
        write(iw,*) '+++++++ calculation.'
        !call quit('+++ Aborting. Matrix too small for parcal.')
      endif
!
!  allocate read/write buffers (master and slaves)
!
      if(carith) then
        allocate(cbuf(intbuf))
        allocate(bufr(intbuf))
        allocate(bufi(intbuf))
      else
        allocate(buf(intbuf))
      endif
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))

! calculate distribution of batches

      allocate(nbufs_per_node(NMPROC))
      ip = nbufs/NMPROC
      iq = mod(nbufs,NMPROC)
      do i=1,NMPROC
        nbufs_per_node(i) = ip
      enddo
      if(iq.gt.0) then
         nbufs_per_node(NMPROC) = nbufs_per_node(NMPROC) + iq
      endif
      write(iw,*) 'Buffer distribution for each node:'
      do i=1,NMPROC
        write(iw,*) 'Node: ',i-1,' number of batches: ',
     &          nbufs_per_node(i)
      enddo
      write(iw,*)
      
! generate local file names and open corresponding file.

      IF (MYPROC .LT. 10) THEN
         WRITE (locfn,'(A6,A1,I1)') nmbase,'.',MYPROC
      ELSEIF (MYPROC .LT. 100) THEN
         WRITE (locfn,'(A6,A1,I2)') nmbase,'.',MYPROC
      ELSEIF (MYPROC .LT. 1000) THEN
         WRITE (locfn,'(A6,A1,I3)') nmbase,'.',MYPROC
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND mpi_xmat_dist routine!")
      ENDIF
      write(*,*) 'Master opens ',locfn,' for writing ADC matrix.'

      inquire(file=locfn,exist=isthere)
      if(isthere) then
        open(unit=iobase,file=locfn,status='old')
        close(iobase,status='delete')
      endif
      open(unit=iobase,file=locfn, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='new')

      jdummy = 0
!
!  allocate required buffers
!
      if(carith) then
        allocate(cadcdiagonal(ladc))
        cadcdiagonal = czero
      else
        allocate(adcdiagonal(ladc))
        adcdiagonal = rzero
      endif

! open complete ADC matrix

      iobase1 = iobase + 1
      open(unit=iobase1,file=nmbase, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='old')
      rewind(iobase1)

! open file for adc diagonal

      io_diag = iobase + 2
      open(unit=io_diag,file='ADCDGTMP', access='SEQUENTIAL',
     &     form='UNFORMATTED',status='unknown')

! write master's own chunk, at the beginning the number of batches

      write(iobase,ERR=889) nbufs_per_node(MASTER+1)
      do i=1,nbufs_per_node(MASTER+1)

        if(carith) then
          read(iobase1,ERR=889) (cbuf(ixx),ixx=1,intbuf),
     &                          (ioi(ixx),ixx=1,intbuf),
     &                          (ioj(ixx),ixx=1,intbuf),
     &                          nact,jdummy
        else
          read(iobase1,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                          (ioi(ixx),ixx=1,intbuf),
     &                          (ioj(ixx),ixx=1,intbuf),
     &                          nact,jdummy
        endif

! search for diagonal entries

        do ixx = 1,nact
          irow = ioi(ixx)
          icol = ioj(ixx)
          if(irow.eq.icol) then
             if(carith) then
               cadcdiagonal(irow) = cbuf(ixx)
             else
               adcdiagonal(irow) = buf(ixx)
             endif
          endif
        enddo

        if(carith) then
          write(iobase,ERR=889) (cbuf(ixx),ixx=1,intbuf),
     &                          (ioi(ixx),ixx=1,intbuf),
     &                          (ioj(ixx),ixx=1,intbuf),
     &                          nact,jdummy
        else
          write(iobase,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                          (ioi(ixx),ixx=1,intbuf),
     &                          (ioj(ixx),ixx=1,intbuf),
     &                          nact,jdummy
        endif
      enddo  !nbufs_per_node
      flush(iobase)
      close(iobase)   !Master's local ADC matrix written.
      write(iw,*) '------MASTER HAS CLOSED LOCAL ADCFILE-------'

! now spread slave chunks

      call interface_mpi_BCAST(SERVER_MATDIS,1,MASTER,
     &                         global_communicator)
      do curslave=MASTER+2,NMPROC

        write(iw,*) 'Master sends',nbufs_per_node(curslave),
     &     ' buffers to slave',curslave-1
        call interface_mpi_send(nbufs_per_node(curslave),1,curslave-1,
     &                          MSGN,global_communicator)

        do i=1,nbufs_per_node(curslave)

          if(carith) then
            read(iobase1,ERR=889) (cbuf(ixx),ixx=1,intbuf),
     &                            (ioi(ixx),ixx=1,intbuf),
     &                            (ioj(ixx),ixx=1,intbuf),
     &                            nact,jdummy
          else
            read(iobase1,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                            (ioi(ixx),ixx=1,intbuf),
     &                            (ioj(ixx),ixx=1,intbuf),
     &                            nact,jdummy
          endif

! search for diagonal entries

          do ixx = 1,nact
            irow = ioi(ixx)
            icol = ioj(ixx)
            if(irow.eq.icol) then
              if(carith) then
                cadcdiagonal(irow) = cbuf(ixx)
              else
                adcdiagonal(irow) = buf(ixx)
              endif
            endif
          enddo

! subroutine interface_mpi_send_i0(x,ndim,receiver,message_id,communicator)
! nodes count from zero, this is the blocking send. Execution will not
! continue until slave has completely received the data chunk

          if(carith) then
            bufr = real(cbuf)
            bufi = aimag(cbuf)
            call interface_mpi_send(bufr,intbuf,curslave-1,
     &                              MSGN,global_communicator)
            call interface_mpi_send(bufi,intbuf,curslave-1,
     &                              MSGN,global_communicator)
          else
            call interface_mpi_send(buf,intbuf,curslave-1,
     &                              MSGN,global_communicator)
          endif
          call interface_mpi_send(ioi,intbuf,curslave-1,
     &                            MSGN,global_communicator)
          call interface_mpi_send(ioj,intbuf,curslave-1,
     &                            MSGN,global_communicator)
          call interface_mpi_send(nact,1,curslave-1,
     &                            MSGN,global_communicator)
!         write(iw,*) ' +++ buffer',i,nact,' sent to slave',curslave-1

        enddo  ! loop over nbufs_per_node

      enddo    ! curslave loop

      if(carith) then
        write(io_diag,ERR=888) (cadcdiagonal(ixx),ixx=1,ladc)
      else
        write(io_diag,ERR=888) (adcdiagonal(ixx),ixx=1,ladc)
      endif

      close(io_diag)

      close(iobase1,status='DELETE')   !remove large Master ADC matrix

      deallocate(nbufs_per_node)
      if(carith) then
        deallocate(cadcdiagonal)
      else
        deallocate(adcdiagonal)
      endif
!
! after distribution we sync all nodes.
!
      deallocate(ioj)
      deallocate(ioi)

      if(carith) then
        deallocate(cbuf)
        deallocate(bufr)
        deallocate(bufi)
      else
        deallocate(buf)
      endif

      return 

 888  call quit('problem R/W ADC diagonal in mpi_master_mat_dist !')
 889  call quit('problem R/W ADC matrix in mpi_master_xmat_dist !')

      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#endif
