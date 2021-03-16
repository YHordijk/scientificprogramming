!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WCONDAT_X(LUN,FILENAME,IRECL,DESREP,RECCOUNT)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Writes explicit configuration information for the excited
!     states in symmetry DESREP.
!     Important: The absolute spinor numbers calculated here correspond
!     to the list in the CCSETI output and not the irepspi numbering!
!     This makes retrieval of orbitals easier for the user
!
!---------------Calling variables--------------------------------------
!
      INTEGER LUN
      CHARACTER*6 FILENAME
      INTEGER IRECL
      INTEGER DESREP
      INTEGER RECCOUNT
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      CHARACTER*70 FIELD
      INTEGER bl_conf
      CHARACTER*4 bl_irna
      
!
!---------------Executable code--------------------------------------
!
! open config file
!
      IRECL = 70
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
!
! First do the P--->H part. Column counter is the decisive
! IRREP.
!
      BL_CONF = 0
      BL_IRNA = '    '
      RECCOUNT = 1

!____________ START WITH H ---> P __________________________
!|
!|
!|

      do jrep = 1,nrep   !jrep,brep are fermionic, desrep is bosonic
        brep = multb(desrep+nrep,jrep,2)
        do jfun = 1,NO(jrep)
        do bfun = 1,NV(brep)
          jabs=IO(jrep) + jfun
          babs=IO(nrep+1) + IV(brep) + bfun

          WRITE(FIELD,333)
     &         babs,brep,REPNA(brep),
     &         jabs,jrep,REPNA(jrep),
     &         BL_CONF,BL_CONF,BL_IRNA,
     &         BL_CONF,BL_CONF,BL_IRNA 
 333  FORMAT(2(I3,' (',I2,',',A4,')  <---- ',I3,' (',I2,',',A4,') '))

          WRITE(LUN,REC=RECCOUNT) FIELD
          RECCOUNT = RECCOUNT + 1

        enddo
        enddo
      enddo
      if( (reccount-1).ne.MVO(desrep))
     &   stop 'Inconsistent rec counter in WCONDAT_X!'
      write(*,*)
!|
!|__________________________________________________________
!
!
!____________ NOW DO  2H ---> 2P ___________________________
!|
!|
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

          kabs =IO(krep) + K
          labs =IO(lrep) + L
          ccabs=IO(nrep+1) + IV(crep) + C
          ddabs=IO(nrep+1) + IV(drep) + D

          WRITE(FIELD,333)
     &         ccabs,crep,REPNA(crep),
     &         kabs,krep,REPNA(krep),
     &         ddabs,drep,REPNA(drep),
     &         labs,lrep,REPNA(lrep)

          WRITE(LUN,REC=RECCOUNT) FIELD
          RECCOUNT = RECCOUNT + 1

                ENDDO
              ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
!|
!|__________________________________________________________

      CLOSE(LUN)

      RECCOUNT = RECCOUNT - 1
 
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVVOO_PP (EPS,T2)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Divide VVOO integrals by denominators ea + eb - ej - ei
!     If we have complex numbers, the program accounts for this by
!     treating real/imaginary part in the *real* array.
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)                 ::  EPS
      REAL*8,dimension(:)                 ::  T2
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 FAC,FAC1,FAC2,FAC3
!
!---------------Executable code--------------------------------------
!
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = - EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = - EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1
                  IF (CARITH) THEN
                     T2(2*ABIJ-1) = T2(2*ABIJ-1)/FAC
                     T2(2*ABIJ)   = T2(2*ABIJ)/FAC
                  ELSE
                     T2(ABIJ) = T2(ABIJ)/FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(ABIJ.ne.IVVOOTT(nrep+1)) THEN
        WRITE(*,*) 'ABIJ/IVVOOTT:',ABIJ,IVVOOTT(NREP+1)
        STOP 'COUNTING ERROR IN DENOMVVOO_PP!'
      ENDIF
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EARR_VVOO (EPS,T2)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     This routine creates an energy denominator array in the
!     VVOO bra-ket storage <|<|,|>|> in order to be sorted with SRT16
!     Normally the VVOO integrals are tridiagonally stored. In this case
!     the full ranges are covered because SRT 16 does not handle tridiagonal
!     storage cases!
!     The energies are all real and for the complex case we set the imaginary
!     part to zero!
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)                 ::  EPS
      REAL*8,dimension(:)                 ::  T2
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 FAC,FAC1,FAC2,FAC3
!
!---------------Executable code--------------------------------------
!
!     call PST('Creating ABIJ energy denominators+')
!     write(*,*) 'jvovo:',jvovo(nrep+1)
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      if(multb(irp,jrp,1).ne.ijrp) stop 'earr error'
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = - EPS(JJ)
         DO I = 1, NO(IRP)
            II = IOFF + I
            FAC2 = - EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
      if(multb(arp,brp,1).ne.ijrp) stop 'earr error'
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               DO A = 1, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1
                  IF(CARITH) THEN
                    T2(2*ABIJ - 1) = 1.0/FAC  !real part
                    T2(2*ABIJ) = 0.0  ! imaginary part
                  ELSE
                    T2(ABIJ) = 1.0/FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(ABIJ.ne.JVOVO(nrep+1)) THEN
        WRITE(*,*) 'ABIJ/JVOVO:',ABIJ,JVOVO(nrep+1)
        STOP 'COUNTING ERROR IN EARR_VVOO!'
      ENDIF
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVOVO(EPS,T)
 
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
! Divides the <a|i><b|j> integral stream by the corresponding
! denominators e_a - e_i + e_b - e_j.
! The stream has JVOVO organization.
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)              :: EPS
      REAL*8,dimension(:)              :: T
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
#include "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8                           :: EA,EI,EB,EJ
      REAL*8                           :: ESM
!
!---------------Executable code--------------------------------------
!
      AIBJ = 0

      DO BJREP = 1, NREP   !bosonic irep
        DO 10 JREP = 1, NREP   !fermionic irep
          JOFF = IO(JREP)
          BREP = MULTB(JREP,BJREP+NREP,2)
          BOFF = IO(NREP+1) + IV(BREP)
          DO J = 1, NO(JREP)
            JJ = JOFF + J
            EJ = eps(JJ)
            DO B = 1, NV(BREP)
               BB = BOFF + B
               EB = eps(BB)
               DO 20 IREP = 1, NREP
                 IOFF = IO(IREP)
                 AREP = MULTB(IREP,BJREP+NREP,2)
                 AOFF = IO(NREP+1) + IV(AREP)
                 DO I = 1, NO(IREP)
                   II = IOFF + I
                   EI = eps(II)
                   DO A = 1, NV(AREP)
                     AA = AOFF + A
                     EA = eps(AA)

                     AIBJ = AIBJ + 1
                     IF (CARITH) THEN
                       ESM = EA+EB-EI-EJ
                       T(2*AIBJ-1) = T(2*AIBJ-1)/ESM
                       T(2*AIBJ) = T(2*AIBJ)/ESM
                     ELSE
                       T(AIBJ) = T(AIBJ)/(EA+EB-EI-EJ)
                     ENDIF

                   ENDDO !A
                 ENDDO !I
 20           ENDDO  !IREP
            ENDDO  !B
          ENDDO  !J
 10     ENDDO !JREP
      ENDDO  !BJREP

      IF((AIBJ).ne.JVOVO(nrep+1)) Then
        write(*,*) 'Consistency error in DENOMVOVO.'
        write(*,*) 'Counter:',AIBJ
        write(*,*) 'Required value:',JVOVO(nrep+1)
        stop
      ENDIF

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMVOTS(eps,aak)
 
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Divides the totally symmetric A_ak array by eps_a - eps_k
!
!---------------Calling variables--------------------------------------
!
      REAL*8,dimension(:)              :: EPS
      REAL*8,dimension(:)              :: AAK
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8                   :: epsa,epsk,ediff
C
C---------------Executable code--------------------------------------
C
      AK = 0
      DO KREP = 1, NREP
        AREP = KREP    ! array is totally symmetric
        KOFF = IO(KREP)
        AOFF = IO(NREP+1) + IV(AREP)
        DO K = 1,NO(KREP)
          DO A = 1,NV(AREP)
            epsa = EPS(AOFF + A)
            epsk = EPS(KOFF + K)
            ediff = epsa - epsk
            AK = AK+1
            AAK(AK) = AAK(AK)/ediff
!           write(*,*) 'DENOMVOTS: ak,aak(ak):',ak,aak(ak)
          ENDDO
        ENDDO
      ENDDO
      IF(AK.ne.MVO(1)) STOP 'Inconsistency in DENOMVOTS!'

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE rcaxpy(n,carith,phase,as,ad)
      Implicit none
!
!---------------- description ----------------
!
! this routine xaxpys from as to ad depending on phase
! ad is always complex, as is either real or complex (following carith)
!
!---------------- formal parameters ----------------
!
      integer                      :: n
      logical                      :: carith
      complex*16                   :: phase
      real*8,dimension(:)          :: as
      complex*16,dimension(:)      :: ad
!
!---------------- local variables ----------------
!
      integer                      :: rcw
      integer                      :: i
!
!---------------- executable code ----------------
!
      rcw = 1
      if(carith) rcw = 2

      if(.not.carith) then       ! ** as = real case + phase **
        do i=1,n
          ad(i) = ad(i) + phase*dcmplx(as(i),0.0d0)
        enddo
      else                       ! ** as = complex case, no phase **
        do i=1,n
          ad(i) = ad(i) + dcmplx(as(2*i-1),as(2*i))
        enddo
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCA_NEGATE(n,arr)
!
!  this routine changes sign of a real or a complex array
!  depending on the CARITH flag. In case of complex:
!  z -->  -z  not  z --> z*  !
!
      implicit none
      integer                      :: n
      real*8, dimension(:)         :: arr
!
! it is assumed that n counts the elements and not the array length!
!

#include  "../relccsd/complex.inc"

      integer                      :: i

      DO I=1,N
        IF(CARITH) THEN
          arr(2*i-1) = -1.0D0*arr(2*i-1)
          arr(2*i)   = -1.0D0*arr(2*i)
        ELSE
          arr(i) = -1.0D0*arr(i)
        ENDIF
      ENDDO
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE print_input_info(iw,nrep)
      use polprp_cfg     !these data are visible on master node only!

      implicit none

      integer               :: iw         ! I/O file handle for proper screen output
      integer               :: nrep       ! number of total symmetries in the calculation
!
      integer               :: i,n,icount

#include "../relccsd/complex.inc"

      write(iw,*) 'POLPRP control variables:'
      write(iw,*) '-------------------------'
      write(iw,*)
      write(iw,*) 'Extended ADC-2 calculation: ',polprp_doextended
      write(iw,*) 'Compute transition moments : ',polprp_dotrmo
      write(iw,*) 'Complex arithmetic: ',carith
      write(iw,'(A,F10.5)') 
     &   ' Threshold for matrix elements: ',polprp_writethr
      write(iw,*) 'Davidson diagonalization:',polprp_dodiag
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + polprp_statesym(i)
      enddo
      if(n.eq.0) then
        write(iw,*) 'All final state symmetries calculated.'
        do i=1,nrep
          polprp_statesym(i) = i
        enddo
      else
        write(iw,*) 'User has selected the following final'
        write(iw,*) 'state symmetries:'
        icount=1
        do i=1,nrep
          if(polprp_statesym(i).ne.0) then
            write(iw,*) '   No: ',icount,
     &                  '  symmetry: ',polprp_statesym(i)
            icount = icount + 1
            if(polprp_statesym(i).gt.nrep) then
              write(iw,*) 'This state symmetry is higher than'
              write(iw,*) 'the number of available symmetries.'
              write(iw,*) 'Setting it to the maximum value.'
              polprp_statesym(i) = nrep
            endif
          endif
        enddo
      endif
      write(iw,*)
!
! report diagonalization parameters
!
      write(iw,*) 'Davidson control variables:'
      write(iw,*) '---------------------------'
      write(iw,*)

      if(polprp_dodiag) then
        WRITE(IW,*) 'Number of eigenstates to be calculated (roots)',
     &     polprp_davroots
        WRITE(IW,*) 'Maximum subspace size: ',polprp_davmaxsp
        WRITE(IW,*) 'Number of Davidson iterations: ',
     &     polprp_davmaxit
        WRITE(IW,'(A,E20.4)') ' Requested energy conv.: ',polprp_davconv
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRECONTRACT_PAR(iounit,iw,eps)
!
!
!  this routine is called by master only. 
!  If we have a parallel run the nodes do their share of the contraction
!        <ab||cd> * <cd||ji>/(ec + ed - ej - ei)
!  and send the partial result back to master where it is combined and
!  stored on disk. Allocation of large VVVV buffers is hereby avoided.
!
!  In serial run it can be assumed that the VVVV integrals fit into a
!  single buffer, otherwise parallel execution is advised.
!

#if defined (VAR_MPI)
      use interface_to_mpi
#endif

      implicit none

      integer                  :: iounit         ! I/O file handle for writing to file
      integer                  :: iw             ! file handle for writing to screen
      real*8, dimension(:)     :: eps
!
! local variables
!
      real*8, allocatable, dimension(:)        :: vvvvz
      real*8, allocatable, dimension(:)        :: vvoo1
      real*8, allocatable, dimension(:)        :: vvoop
      real*8, allocatable, dimension(:)        :: vvoos
      integer                                  :: n1,ixx,inode
      integer                                  :: ialloc
      integer                                  :: signal
      integer*8                                :: n8,i8
      real*8                                   :: R0
      COMPLEX*16                               :: A0,A1
      character*9                              :: preconfn='VVVVXVVOO'
      logical                                  :: isthere


#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/ccpar.inc"
#include  "polprp_servercodes.h"
!
!   function interfaces
!
      interface

        SUBROUTINE denomvvoo_pp(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE do_partial_precon(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

      end interface
!
! execution
!
      if(carith.and.rcw.eq.2) then
        write(iw,*) ' -> Performing (complex) precontraction'
      else if(.not.carith.and.rcw.eq.1) then
        write(iw,*) ' -> Performing (real) precontraction'
      else
        write(iw,*) ' Control error in precontract_new!'
        call quit('error')
      endif
      R0=0.0d0
      A0=(0.0D00,0.0D00)
      A1=(1.0D00,0.0D00)

! allocate buffers and get integrals needed in serial and parallel

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoop(n1*rcw))

      vvoop = R0

      call getvvoo(vvoo1)
      call denomvvoo_pp(eps,vvoo1)
!____________________________________________________________________
!|
!|  Do parallel execution with partial contractions for each VVVV batch
!|
#if defined (VAR_MPI)
      write(iw,*) ' -> in parallel mode.'
      call interface_mpi_bcast(SERVER_PRECON,1,MASTER,
     &                         global_communicator)

!  sending  <cd||ji>/(ec + ed - ej - ei) to the slaves

      do inode=master+1,nmproc-1
        write(*,*) 'Sending vvoo stream to node',inode
        call interface_mpi_send(vvoo1,n1*rcw,inode,
     &                          MSGN,global_communicator)
      enddo

! also perform local partial contraction

      vvoop = 0.0d0
      call do_partial_precon(vvoo1,vvoop)

! collect the partial streams from the slaves

      allocate(vvoos(n1*rcw))
      do inode=master+1,nmproc-1
        vvoos = R0
        write(*,*) 'Fetching vvoo stream from node',inode
        call interface_mpi_recv(vvoos,n1*rcw,inode,MSGN,
     &                          global_communicator)
        vvoop = vvoop + vvoos
      enddo
      deallocate(vvoos)

#else
!____________________________________________________________________
!|
!|  Do serial execution using complete VVVV batch
!|
      write(iw,*) ' -> in serial mode.'
      n8 = ivvvvtt(nrep+1)
      allocate(vvvvz(n8*rcw),stat=ialloc)
      if(ialloc.ne.0) then
        STOP 'VVVV buffer could not be allocated.'
      else
        write(*,*) n8,' words for VVVV buffer were allocated.'
      endif
      call rdvvvv(vvvvz)  ! master reads complete stream in serial mode
      CALL CNTRCT('N','N',NVVT,NOOT,NVVT,A1,vvvvz,vvoo1,A0,vvoop,NREP)
      deallocate(vvvvz)
#endif
!|  
!|    from here the internal state is identical in ser/par:
!|    contracted integrals sit in vvoop array and are written to disk
!|___________________________________________________________________
!
      write(*,*) 'Precontracted VVOO stream chksum:',
     &    dot_product(vvoop,vvoop)
      Inquire(file=preconfn,exist=isthere)
      If(isthere) then
         write(*,*) 'Deleting previous precontraction file.'
         open(unit=iounit, file=preconfn, status='old')
         close(iounit, status='delete')
      Endif
!
!  open corresponding file and write out contracted integrals.
!
      open(unit=iounit,file=preconfn, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
      write(iounit) (vvoop(ixx), ixx=1,n1*rcw)
      close(iounit)
      write(*,*) 'File ',preconfn,' generated.'
!
!  free buffers used in ser/par mode.
!
      deallocate(vvoop)
      deallocate(vvoo1)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine do_partial_precon(vvoo_in,vvoo_ou)
      implicit none
!
      real*8,allocatable,dimension (:)    :: vvoo_in, vvoo_ou

#include  "polprp_servercodes.h"
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
#include  "../relccsd/complex.inc"
!
! local variables
! 
      integer                             :: irp,istart,mint,nbatch
      integer                             :: n2,nint,m,n,k,localstart
      integer                             :: off1,off2,off3
      real*8,allocatable,dimension(:)     :: bufr
      logical                             :: done,deb
      complex*16,parameter                :: cone=(1.0d0,0.0d0)
!
! executable code
! 
      deb = .false.
      If(SERVER_TRACE.ne.0) deb = .true.

      if(deb) then
        n2 = idist(1,1,irp)
        write(*,*) 'ioff start on node',myproc,n2
        n2 = idist(2,1,irp)
        write(*,*) 'ioff end on node',myproc,n2
        n2 = idist(2,1,irp) - idist(1,1,irp)
        write(*,*) '# of NVVT batches on node',myproc,n2
        n2 = idist(3,1,irp)
        write(*,*) 'absolute start add. of this batch: ',myproc,n2
        n2 = idist(4,1,irp)
        write(*,*) 'absolute end add. of this batch: ',myproc,n2
        n2 = idist(5,1,irp)
        write(*,*) 'relative start add. of this batch: ',myproc,n2
        n2 = idist(6,1,irp)
        write(*,*) 'relative end add. of this batch: ',myproc,n2
      endif

      Do irp = 1,nrep
        m = nvvt(irp)  !maximum number of batches in this irrep
        n = noot(irp)
        IF(m*n.GT.0) THEN
          istart = 0
          nbatch = idist(2,1,irp) - idist(1,1,irp) !current batches on this node
          n2 = nbatch * nvvt(irp) ! number of VVVV integrals on this node

          allocate(bufr(n2*rcw))
          CALL GETVVVV (irp,istart,nint,done,bufr,m)
          if(nint.ne.nbatch.or.(.not.done)) then
            write(*,*) 'Serious batch distribution error!'
            write(*,*) 'Contact programmer of the module.'
            call quit()
          else
            if(deb) write(*,*) 'Node',myproc,' local VVVV checksum:',
     &        dot_product(bufr,bufr)
          endif
          localstart = idist(1,1,irp)
          off1 = 1
          off2 = (ivvoott(irp) + localstart)*rcw + 1
          off3 = ivvoott(irp)*rcw + 1
          CALL XGEMM ('N','N',m,n,nbatch,cone,bufr(off1),m,
     &              vvoo_in(off2),m,cone,vvoo_ou(off3),m)

          deallocate(bufr)
          if(deb) write(*,*) 'Node',myproc,'$$ loc. precon checksum:',
     &       dot_product(vvoo_ou,vvoo_ou)
        ENDIF
      Enddo  !irp

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine compare_evc_r(io1,desrep,nroots,ladc,
     &                         davibase,fullbase)
      implicit none
!
! this routine evaluates eigenvectors for the available Davidson roots
! with respect to the full diagonalization results. It determines phases and
! absolute deviations
!
      integer        :: io1,desrep,nroots,ladc
      character*8    :: davibase
      character*7    :: fullbase 
!
! local variables
! 
      character*11                        :: daviname
      character*10                        :: fullname
      logical                             :: fullisthere
      integer                             :: io2
      integer                             :: i,j,k
      integer                             :: k1,k2
      integer                             :: nmain
      real*8                              :: lambda1,lambda2,dummy
      real*8,parameter                    :: autoev = 27.2113957D0
      real*8,allocatable,dimension (:)    :: e1,e2
      real*8                              :: devmod,devabs,devphs
!
! executable code
! 
      Call PST('Comparing D/F eigenvectors+')

      IF(desrep.GT.9) THEN
        WRITE(daviname,'(A8,A1,I2)') davibase,'.',desrep
        WRITE(fullname,'(A7,A1,I2)') fullbase,'.',desrep
      ELSE
        WRITE(daviname,'(A8,A2,I1)') davibase,'.0',desrep
        WRITE(fullname,'(A7,A2,I1)') fullbase,'.0',desrep
      ENDIF
!
! do we have an eigenvector file from fulldia ?
!
      fullisthere=.false.
      INQUIRE(FILE=fullname,EXIST=fullisthere)
      if(.not.fullisthere) then
        write(*,*) 'The eigenvector file from a full diagonalization'
        write(*,*) 'is not available. No comparison is made. Returning.'
        return
      endif
      write(*,*) 'Comparing eigenvectors for',nroots,' roots'
      write(*,*)
!
! comparing eigenvectors
!
      io2 = io1 + 1
      open(io1,file=daviname,access='sequential',
     &     status='unknown',form='unformatted')
      open(io2,file=fullname,access='sequential',
     &     status='unknown',form='unformatted')

      allocate(e1(ladc))
      allocate(e2(ladc))
      e1 = 0.0d0
      e2 = 0.0d0
!
! three integer reads from the Davidson file
!
      read(io1) k
      read(io1) nmain
      read(io1) k
      write(*,*) 'Main space has dimension',nmain
      write(*,*)

      write(*,'(18X,2A7,4X,3A8)') 'root1','root2',' devmod',
     &                          'devabs','devphs'
      do k=1,nroots
        read(io1) k1
        read(io2) k2
        if(.not.(k1.eq.k2 .and. k1.eq.k)) then
          write(*,*) 'Inconsistency in root numbering!'
          write(*,*) 'This is serious, continuation questionable!'
          write(*,*) 'k,k1,k2:',k,k1,k2
          exit  !leave this do loop
        endif
        read(io1) lambda1
        read(io2) lambda2
        read(io1) dummy
        read(io2) dummy
        read(io1) (e1(i),i=1,ladc)
        read(io2) (e2(i),i=1,ladc)
!
! compare vectors
!
        devmod = 0.0d0
        devabs = 0.0d0
        devphs = 0.0d0
        do i=1,ladc
          devmod = devmod + (abs(e1(i)) - abs(e2(i)))
          devabs = devabs + (e1(i) - e2(i))
          devphs = devphs + (e1(i) + e2(i))
        enddo
        write(*,'(A,I3,2X,2F12.6,2X,3F14.8)') 'Comparing root',k,
     &         lambda1,lambda2,devmod,devabs,devphs
        write(*,*)
        do i=1,min(10,nmain)
          write(*,*) '       ',e1(i),e2(i)
        enddo
        write(*,*) '..................'
        write(*,*)
      enddo
!
! release ressources
!
      deallocate(e1)
      deallocate(e2)

      close(io1)
      close(io2)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      function qmasymh(e,lde)
!
!     computes deviation from Hermiticity of a quadratic matrix
!
      implicit none
      real*8                       :: qmasymh
      complex*16, dimension(:,:)   :: e
      integer                      :: lde

! local variables

      real*8                       :: asym
      integer                      :: i,j
!
!---------------executable code--------------------------------------
!
      asym=0.0D0

      if(lde.eq.1) then
       qmasymh = asym
       return
      endif

      do i=2,lde
        do j=1,i-1
          asym = asym + abs(dconjg(e(i,j)) - e(j,i))
        enddo
      enddo

      qmasymh = asym

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATXMAT(n,a,b,c)
 
      IMPLICIT none
!
!---------------Description--------------------------------------------
!
! This routine performs an *elementwise* multiplication as
!                           c = a * b
! a,b,c can be real or complex. The usefulness of this routine is that
! the reference always occurs through real arrays as in the X... routines.
!
!---------------Calling variables--------------------------------------
!
      Integer                          :: n
      REAL*8,dimension(:)              :: a,b,c
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/complex.inc"
!
!---------------Local Variables--------------------------------------
!
      Integer         :: i
      Complex*16      :: s1,s2,s3
!
!---------------Executable --------------------------------------
!
      If(carith) then
        Do i=1,n
          s1=dcmplx(a(2*i-1),a(2*i))
          s2=dcmplx(b(2*i-1),b(2*i))
          s3 = s1 * s2
          c(2*i-1) = dble(s3)
          c(2*i)   = aimag(s3)
        Enddo
      Else
        Do i=1,n
          c(i) = a(i) * b(i)
        Enddo
      Endif

      Return
      End
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1T3_NTS (NREP,MULTB,DOINV,NPAIR,KFIE,LFIE,NTRIPL,
     &                       OFF,OFF1,BUF1,BUF2,DESREP,ICEIL)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     BUF(IJ,K>L:KLREP) <--> BUF(IJ K,L:LREP)
!     where the product of all four ireps can be non-totally symmetric
!     If (DOINV) inverse sort.
!
!                 <IJ| x |KL> = desrep
!
!---------------Calling variables--------------------------------------
!
      INTEGER                  NREP
      INTEGER                  MULTB(64,64,2)
      LOGICAL                  DOINV
      INTEGER                  NPAIR(NREP),KFIE(NREP)
      INTEGER                  LFIE(NREP),NTRIPL(NREP)
      INTEGER                  OFF(NREP+1),OFF1(32,32)
      REAL*8                   BUF1(*),BUF2(*)
      INTEGER                  DESREP,ICEIL
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
       integer summer,sumkl
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      summer = 0
      IJKL = 1
      DO KLREP = 1, NREP

       IJREP = MULTB (DESREP+NREP,KLREP+NREP,2)
       IF(MULTB(ijrep+nrep,klrep+nrep,2).ne.desrep) 
     &        stop 'Symm error 1 in SRT_1T3_NTS'

       N = NPAIR(IJREP)
       DO 10 LREP = 1, NREP
        KREP = MULTB(LREP,KLREP+NREP,2)
        IF (KREP.LT.LREP) GOTO 10

        IJKREP = MULTB (IJREP+NREP,KREP,2)
        IJLREP = MULTB (IJREP+NREP,LREP,2)
        if(multb(ijkrep,lrep,1).ne.desrep) stop 'Symm error 2'
        if(multb(ijlrep,krep,1).ne.desrep) stop 'Symm error 3'

        LOFF = (OFF(LREP) + OFF1(IJREP,KREP)) * RCW
        KOFF = (OFF(KREP) + OFF1(IJREP,LREP)) * RCW

        DO L = 1, LFIE(LREP)
         IJL  = (L-1) * N * RCW + 1
         KMIN = 1
         IF (KREP.EQ.LREP) KMIN = L + 1
         DO K = KMIN, KFIE(KREP)
            IJK  = (K-1) * N * RCW + 1
            IJKL1 = LOFF + (L-1) * NTRIPL(IJKREP) * RCW + IJK 
            IJLK1 = KOFF + (K-1) * NTRIPL(IJLREP) * RCW + IJL 
            if(IJKL1+N.gt.iceil*rcw+1) write(*,*) '*** serious alert A'
            if(IJLK1+N.gt.iceil*rcw+1) write(*,*) '*** serious alert B'
            IF (DOINV) THEN
               CALL XCOPY (N,BUF1(IJKL1),1,BUF2(IJKL),1)
               CALL XAXPY (N,-A1,BUF1(IJLK1),1,BUF2(IJKL),1)
            ELSE
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJKL1),1)
               CALL XCOPY (N,BUF1(IJKL),1,BUF2(IJLK1),1)
               CALL XSCAL (N,-A1,BUF2(IJLK1),1)
            ENDIF
            IJKL = IJKL + N * RCW
         ENDDO
        ENDDO
 10    CONTINUE
      ENDDO
!
      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRT1T2_NTS (NREP,MULTB,DOINV,NPAIR1,IFIE,JFIE,
     &                       NPAIR2,OFF,OFF2,BUF1,BUF2,DESREP)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     BUF(I>J,KL:KLREP) <--> BUF(I,J KL:JKLREP)
!     If (DOINV) inverse sort.
!
!---------------Last modified------------------------------------------
!
!    by MP for nontotally symmetric cases.
!
!---------------Calling variables--------------------------------------
!
      REAL*8 BUF1(*),BUF2(*)
      INTEGER NREP,DESREP,MULTB(64,64,2)
      INTEGER NPAIR1(NREP),IFIE(NREP),JFIE(NREP),NPAIR2(NREP)
      INTEGER OFF(NREP+1),OFF2(32,32)
      LOGICAL DOINV
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/param.inc"
!
!---------------Executable code--------------------------------------
!
      IF (.NOT.DOINV) THEN
         CALL XCOPY(OFF(NREP+1),A0,0,BUF2,1)
      ENDIF
      IJKL = 0
      DO KLREP = 1, NREP

! new for nontotally symmetric cases:
       IJREP = MULTB (DESREP+NREP,KLREP+NREP,2)
       IF(MULTB(ijrep+nrep,klrep+nrep,2).ne.desrep) stop 'Error!'

       M = NPAIR1(IJREP)
       N = NPAIR2(KLREP)
       IJ = 1
       DO 10 JREP = 1, NREP
        IREP = MULTB(JREP,IJREP+NREP,2)
        IF (IREP.LT.JREP) GOTO 10

! new for nontotally symmetric cases:

        JKLREP = MULTB(JREP,KLREP+NREP,2)
        IKLREP = MULTB(IREP,KLREP+NREP,2)

        JKLOFF = (OFF(JKLREP)+OFF2(JREP,KLREP)*IFIE(IREP)) * RCW
        IKLOFF = (OFF(IKLREP)+OFF2(IREP,KLREP)*JFIE(JREP)) * RCW
        M2 = IFIE(IREP) * JFIE(JREP)
        DO J = 1, JFIE(JREP)
         IMIN = 1
         IF (IREP.EQ.JREP) IMIN = J + 1
         DO I = IMIN, IFIE(IREP)
            IJKL1 = JKLOFF + ((J-1) * IFIE(IREP) + I - 1) * RCW + 1
            JIKL1 = IKLOFF + ((I-1) * JFIE(JREP) + J - 1) * RCW + 1
            IF (DOINV) THEN
               CALL XCOPY(N,BUF1(IJKL1),M2,BUF2(IJKL+IJ),M)
               CALL XAXPY(N,-A1,BUF1(JIKL1),M2,BUF2(IJKL+IJ),M)
            ELSE
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(IJKL1),M2)
               CALL XCOPY(N,BUF1(IJKL+IJ),M,BUF2(JIKL1),M2)
               CALL XSCAL(N,-A1,BUF2(JIKL1),M2)
            ENDIF
            IJ = IJ + RCW
         ENDDO
        ENDDO
 10    CONTINUE
       IJKL = IJKL + M * N * RCW
      ENDDO

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Function get_file_unit()
!
! returns a unit number for file I/O that is not in use yet
!
      Implicit none

#define LUMIN 70
#define LUMAX 99

      Integer                   :: get_file_unit
      integer                   :: lu,iostat
      logical                   :: istaken
!
      istaken = .true.
      do lu = LUMIN, LUMAX
         inquire (unit=lu, opened=istaken, iostat=iostat)
         if (iostat.ne.0) cycle
         if (.not.istaken) exit
      end do
!
      if(istaken) then
        write(*,*) ' *** Attention! Reached LUMAX=',LUMAX
        write(*,*) ' *** No free file handle available. Terminating.'
        call quit('I/O problem')
      else
!       write(*,*) 'Obtained file handle',lu
      endif
      get_file_unit = lu
      return
      end Function
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DISORT(n,arr,brr)
!
! adapted from numerical recipes in fortran 77
! to simultaneously shuffle the index array brr
! with the sorted array arr
!
      implicit none
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      INTEGER brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,atemp
      INTEGER b,btemp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        atemp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=atemp
        btemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=btemp
        if(arr(l+1).gt.arr(ir))then
          atemp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l).gt.arr(ir))then
          atemp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=atemp
          btemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=btemp
        endif
        if(arr(l+1).gt.arr(l))then
          atemp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=atemp
          btemp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=btemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        atemp=arr(i)
        arr(i)=arr(j)
        arr(j)=atemp
        btemp=brr(i)
        brr(i)=brr(j)
        brr(j)=btemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          stop 'increase NSTACK in DISORT'
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      double complex function pp_zdotc(n,zx,zy)

      integer                           :: n
      double complex, dimension(:)      :: zx,zy
!
!   local variables
!
      double complex                    :: a
      integer                           :: i
      double complex,parameter          :: czero = (0.0d0,0.0d0)
!
!   executable code
!
      a = czero

      do i=1,n
        a = a + dconjg(zx(i)) * zy(i)
      enddo

      pp_zdotc = a

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      real*8 function rzvecnorm(n,z)
!
! this function calculates the (real) norm of a complex vector of length n.
!
      integer                     :: n
      complex*16, dimension (:)   :: z

      integer          :: i
      complex*16       :: sum,x

      sum = (0.0d0,0.0d0)
      do i=1,n
        x = z(i)
        sum = sum + dconjg(x) * x
      enddo
      if(abs(aimag(sum)).gt.1.0E-10) then
        write(*,*) 'Too large imaginary part in rzvecnorm!'
        call quit('numerical problem')
      endif
      rzvecnorm = sqrt(real(sum))

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xmvmul_sp_r(iobase,intbuf,ladc,PO,PN)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase,intbuf,ladc
      real*8,dimension(:)            :: PO,PN
!
!------------------------ local variables -------------------------
!
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A
      integer                                 :: nbufsloc
      integer                                 :: ioresult

!
!---------------Description--------------------------------------------
!
!  Perform the real symmetric multiplication of the partial or complete
!  ADC matrix with PO and store the result in the vector PN
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  nbufsloc is an indicator if the ADC matrix is complete or partial
!  In serial mode it is complete and has -1 at the beginning. Here we need the
!  iostat mechanism to determine the file end
!
!  In parallel mode nbufsloc contains the actual number of buffers. This is known due
!  to the matrix distribution procedure.
!
!  By this we avoid an awkward handling of the number of buffers as arguments, or common blocks...
!
      pn = 0.0d0

      rewind(iobase)
      read(iobase) nbufsloc

      If(nbufsloc.lt.0) then  !any negative number indicates serial
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

      Do

        READ(iobase,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                               (IOI(IXX),IXX=1,INTBUF),
     &                               (IOJ(IXX),IXX=1,INTBUF),
     &                               NACT,JDUMMY
        if (ioresult > 0) then
           call quit('ADC matrix read error (real). Terminating')
        endif
        if (ioresult < 0) then  ! EOF detected
           exit    !leave infinite do loop
        endif

! useful data in buffer

        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO

      Enddo
!|
!|   
!|__________________________________________________________
      else
!___________________________________________________________
!|
!|   parallel case, exact buffer number is known
!|
        DO I = 1,nbufsloc
          READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                         (IOI(IXX),IXX=1,INTBUF),
     &                         (IOJ(IXX),IXX=1,INTBUF),
     &                          NACT,JDUMMY
          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

            PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
          ENDDO
        ENDDO
!|
!|   
!|__________________________________________________________

      endif

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in xmvmul_sp_r')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE xmvmul_sp_c(iobase,intbuf,ladc,PO,PN)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase,intbuf,ladc
      complex*16,dimension(:)        :: PO,PN
!
!------------------------ local variables -------------------------
!
      complex*16, allocatable, dimension(:)   :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      complex*16                              :: A
      integer                                 :: nbufsloc
      integer                                 :: ioresult

!
!---------------Description--------------------------------------------
!
!  Perform the hermitian multiplication of the external partial 
!  ADC matrix with PO and store the result in the vector PN
!
!---------------Executable code--------------------------------------
!
!  for the matrix elements we note that only the lower triangular
!  part *AND* the diagonal of the ADC-matrix is stored. 
!  the if takes care of not double counting the diagonal.

      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
!
!  nbufsloc is an indicator if the ADC matrix is complete or partial
!  In serial mode it is complete and has -1 at the beginning. Here we need the
!  iostat mechanism to determine the file end
!
!  In parallel mode nbufsloc contains the actual number of buffers. This is known due
!  to the matrix distribution procedure.
!
!  By this we avoid an awkward handling of the number of buffers as arguments, or common blocks...
!
      pn = (0.0d0,0.0d0)

      rewind(iobase)
      read(iobase) nbufsloc

      If(nbufsloc.lt.0) then  !any negative number indicates serial
!___________________________________________________________
!|
!|   serial case, exact buffer number is not known
!|   we read stream blocks until end of file is reached.

        Do

          READ(iobase,iostat=ioresult) (BUF(IXX),IXX=1,INTBUF),
     &                                 (IOI(IXX),IXX=1,INTBUF),
     &                                 (IOJ(IXX),IXX=1,INTBUF),
     &                                 NACT,JDUMMY
          if (ioresult > 0) then
             call quit('ADC matrix read error (complex). Terminating')
          endif
          if (ioresult < 0) then  ! EOF detected
             exit    !leave infinite do loop
          endif

! useful data in buffer

          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

             PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW)*dconjg(A)
          ENDDO

        Enddo
!|
!|   
!|__________________________________________________________
      else
!___________________________________________________________
!|
!|   parallel case, exact buffer number is known
!|
        DO I = 1,nbufsloc
          READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                         (IOI(IXX),IXX=1,INTBUF),
     &                         (IOJ(IXX),IXX=1,INTBUF),
     &                          NACT,JDUMMY
          DO K = 1, NACT
            IROW = IOI(K)
            ICOL = IOJ(K)
            A    = BUF(K)

            PN(IROW) = PN(IROW) + PO(ICOL) * A
            IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW)*dconjg(A)
          ENDDO
        ENDDO
!|
!|   
!|__________________________________________________________

      endif

      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 889  CALL QUIT('Matrix read error in xmvmul_sp_c')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDIAG_INI_R(iobase,intbuf,nbufs,ladc)

!  Fetch the diagonal entries of the sparse matrix A, needed for the
!  Davidson iterations.
!  This routine is only called in the serial case and obtains NBUFS.
!  It generates a local file handle and writes out diagonal.

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase
      integer                        :: intbuf,nbufs,ladc
!
!------------------------ local variables -------------------------
!
      Character*8,parameter                   :: diagfile='ADCDGTMP'
      integer                                 :: serialindicator
      REAL*8,  allocatable, dimension(:)      :: BUF
      INTEGER, allocatable, dimension(:)      :: IOI,IOJ
      REAL*8,  allocatable, dimension(:)      :: adiag
      INTEGER                                 :: I,K,IXX
      INTEGER                                 :: NACT,JDUMMY,IROW,ICOL
      REAL*8                                  :: A
      INTEGER                                 :: io_out

!
!--------------- Interfaces --------------------------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

      end interface
!
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      allocate(adiag(ladc))

!  .. read elements

      REWIND(iobase)
      read(iobase) serialindicator   !discarded, but should be less than zero!

      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          IF(IROW.EQ.ICOL) adiag(irow) = A

        ENDDO
      ENDDO

      io_out = get_file_unit()
      open(unit=io_out,file=diagfile, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
      write(io_out,ERR=887) (adiag(ixx),ixx=1,ladc)
      close(io_out)

      deallocate(adiag)
      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 887  CALL QUIT('Matrix write error in XDIAG_INI_R')
 889  CALL QUIT('Matrix read error in XDIAG_INI_R')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XDIAG_INI_C(iobase,intbuf,nbufs,ladc)

! complex extension of XDIAG_INI_R

      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      integer                        :: iobase
      integer                        :: intbuf,nbufs,ladc
!
!------------------------ local variables -------------------------
!
      Character*8,parameter                     :: diagfile='ADCDGTMP'
      integer                                   :: serialindicator
      complex*16, allocatable, dimension(:)     :: BUF
      INTEGER, allocatable, dimension(:)        :: IOI,IOJ
      complex*16, allocatable, dimension(:)     :: adiag
      INTEGER                                   :: I,K,IXX
      INTEGER                                   :: NACT,JDUMMY,IROW,ICOL
      complex*16                                :: A
      INTEGER                                   :: io_out
!
!--------------- Interfaces --------------------------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

      end interface
!
!---------------Executable code--------------------------------------
!
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      allocate(adiag(ladc))

c  .. read elements

      REWIND(iobase)
      read(iobase) serialindicator   !discarded, but should be less than zero!

      DO I = 1,NBUFS
        READ(iobase,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          IF(IROW.EQ.ICOL) adiag(irow) = A

        ENDDO
      ENDDO

      io_out = get_file_unit()
      open(unit=io_out,file=diagfile, access='SEQUENTIAL',
     &     form='UNFORMATTED',status='NEW')
       write(io_out,ERR=887) (adiag(ixx),ixx=1,ladc)
       close(io_out)

      deallocate(adiag)
      deallocate(ioj)
      deallocate(ioi)
      deallocate(buf)

      RETURN

 887  CALL QUIT('Matrix write error in XDIAG_INI_C')
 889  CALL QUIT('Matrix read error in XDIAG_INI_C')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ortho_warn()

      implicit none
!
!---------------Description--------------------------------------------
!
! prints out a descriptive warning message for the user
!

      write(*,*)
      write(*,*)
      write(*,*) '*****************************************************'
      write(*,*) '*'
      write(*,*) '*           W A R N I N G'
      write(*,*) '*'
      write(*,*) '*  Orthogonality of ADC eigenvectors starts to'
      write(*,*) '*  disappear!'
      write(*,*) '*'
      write(*,*) '*  This is due to the fact that the current Krylov'
      write(*,*) '*  space dimension becomes insufficient in order to'
      write(*,*) '*  provide enough information for the vector search.'
      write(*,*) '*'
      write(*,*) '*  We stop the iterations here. Otherwise'
      write(*,*) '*  the already converged vectors lose significance.'
      write(*,*) '*'
      write(*,*) '*  You should restart the calculation with an'
      write(*,*) '*  increased value for DVMAXSP.'
      write(*,*) '*'
      write(*,*) '*  We set number of macroiterations to maximum value.'
      write(*,*) '*  You will see this in the corresponding output.'
      write(*,*) '*'
      write(*,*) '*****************************************************'
      write(*,*)
      write(*,*)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
