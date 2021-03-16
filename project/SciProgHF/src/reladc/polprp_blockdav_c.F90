!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine mpi_block_davidson_c(iobase, &
                                      intbuf, &
                                      ladc, &
                                      nmain, &
                                      desrep, &
                                      nmbase)
!
!  real Davidson diagonalizer for serial and parallel runs.
!
#if defined (VAR_MPI)
      use interface_to_mpi
#endif
      use qstack    ! activate stack system
      use polprp_cfg   ! make davidson control variables accessible.

      Implicit none
!
!  formal parameters of the subroutine
!
      integer                               :: iobase
      integer                               :: intbuf
      integer                               :: ladc
      integer                               :: nmain
      integer                               :: desrep
      character*6                           :: nmbase
!
!  control variables for Davidson
!
      integer                               :: nroots
      integer                               :: maxdavsp
      integer                               :: maxdavit
      real*8                                :: convthr

#if defined (VAR_MPI)

#include  "../relccsd/ccpar.inc"
#include  "polprp_servercodes.h"

#endif

#include  "polprp_stacklines.h"
!
!  local variables
!
      real*8,parameter                       :: rzero = 0.0d0
      real*8,parameter                       :: prethr = 0.001d0
      real*8                                 :: prenorm1,prenorm2
      real*8                                 :: time_totl_s,time_totl_e
      real*8                                 :: time_root_s,time_root_e
      real*8                                 :: time_kryl_s,time_kryl_e
      real*8                                 :: time_resi_s,time_resi_e
      real*8                                 :: time_cur
      real*8                                 :: ortho_cur
      real*8,allocatable,dimension(:)        :: lambda
      real*8,allocatable,dimension(:)        :: resnorm
      real*8,allocatable,dimension(:)        :: rwork
      real*8,allocatable,dimension(:)        :: mpiaux
      real*8,allocatable,dimension(:)        :: mpires1,mpires2
      real*8,allocatable,dimension(:)        :: bwc   ! aux array for qstack interoperability
                                             
      complex*16,parameter                   :: czero = (rzero,rzero)
      complex*16                             :: hdx
      complex*16,allocatable,dimension(:)    :: adiag
      complex*16,allocatable,dimension(:)    :: bv,bw,bx
      complex*16,allocatable,dimension(:,:)  :: hij,hijcal
      complex*16,allocatable,dimension(:)    :: prc
      complex*16,allocatable,dimension(:)    :: pro
      complex*16,allocatable,dimension(:,:)  :: res
      complex*16,allocatable,dimension(:)    :: work
                                             
      integer                                :: i,j,ix
      integer                                :: krydim,krysav
      integer                                :: maciter
      integer                                :: io_bold
      integer                                :: lwork,info
      integer                                :: kroot
      integer                                :: ida_iostat
      integer                                :: iallocstat
      Integer                                :: qst_hboline
                                             
      logical                                :: allconverged
      logical                                :: oviol
      logical,allocatable,dimension(:)       :: converged
                                             
      character*10                           :: fileadc
      character*60                           :: listout
!
!  interfaces for subroutines
!
      interface

        SUBROUTINE init_krylov_c(i1,i2,i3,i4,i5,l1)
          integer                          :: i1,i2,i3,i4,i5
          logical                          :: l1
        END SUBROUTINE

        SUBROUTINE xmvmul_sp_c(i1,i2,i3,ca1,ca2)
          integer                          :: i1,i2,i3
          complex*16, dimension(:)         :: ca1,ca2
        END SUBROUTINE

        SUBROUTINE calc_resvec_c(i1,i2,i3,i4,ca1,ra1,ca2)
          integer                      :: i1,i2,i3,i4
          complex*16, dimension(:,:)   :: ca1
          real*8, dimension(:)         :: ra1
          complex*16, dimension(:,:)   :: ca2
        END SUBROUTINE

        FUNCTION rzvecnorm(i1,ca1)
          real*8                           :: rzvecnorm
          integer                          :: i1
          complex*16, dimension(:)         :: ca1
        END FUNCTION

        SUBROUTINE calc_precon_c(i1,ca1,r1,ca2,ca3)
          integer                      :: i1
          complex*16, dimension(:)     :: ca1
          real*8                       :: r1
          complex*16, dimension(:)     :: ca2
          complex*16, dimension(:)     :: ca3
        END SUBROUTINE

        SUBROUTINE calc_preort_c(i1,i2,i3,ca1,ca2)
          integer                      :: i1,i2,i3
          complex*16, dimension(:)     :: ca1
          complex*16, dimension(:)     :: ca2
        END SUBROUTINE

        SUBROUTINE write_davevecs_c(i1,i2,i3,i4,i5,i6,ra1,ra2,ca1,r1)
          integer                      :: i1,i2,i3,i4,i5,i6
          real*8, dimension(:)         :: ra1
          real*8, dimension(:)         :: ra2
          complex*16, dimension(:,:)   :: ca1
          real*8                       :: r1
        END SUBROUTINE

      end interface
!
! say hello and mark start time
!
      Call PST('Entering complex (S/P) Davidson diagonalizer+')
      call cpu_time(time_totl_s)
!         
!  assign stack line number in order to avoid direct integers
!  as actual parameters for the qstack calls (I8!)
!       
      qst_hboline = HBO_STACKLINE
!
! initialize davidson control variables from POLPRP input
!
      nroots   =   polprp_davroots   !number of requested Davidson roots
      maxdavsp =   polprp_davmaxsp   !maximum dimension of Davidson space
      maxdavit =   polprp_davmaxit   !maximum number of Dav. iterations
      convthr  =   polprp_davconv    !convergence of Davidson eigenvectors
!
! the orthogonalization violation flag is false at the beginning
!
      oviol = .false.
!
! allocate and initialize arrays according to calculation dimensions
!
! complex arrays:

      allocate(adiag(ladc))
         adiag = czero
      allocate(bv(ladc))
         bv = czero
      allocate(bw(ladc))
         bw = czero
      allocate(bx(ladc))
         bx = czero
      allocate(hij(maxdavsp,maxdavsp))
         hij = czero
      allocate(hijcal(maxdavsp,maxdavsp))
         hijcal = czero
      allocate(prc(ladc))
         prc = czero
      allocate(pro(ladc))
         pro = czero
      allocate(res(ladc,nroots))
         res = czero

! real arrays:

      allocate(lambda(maxdavsp))
         lambda = rzero
      allocate(resnorm(nroots))
         resnorm = rzero
      lwork = 3*maxdavsp
      allocate(work(lwork))
         work = czero
      allocate(rwork(lwork))
         rwork = rzero
      allocate(mpiaux(ladc))
         mpiaux = rzero
      allocate(mpires1(ladc))
         mpires1 = rzero
      allocate(mpires2(ladc))
         mpires2 = rzero
      allocate(bwc(2*ladc))
         bwc = rzero
      allocate(converged(nroots))
         converged = .false.

!*******************************************************
! print essential information to the user
!*******************************************************

      WRITE(*,*)
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*) ' ------ Davidson parameters -----------'
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*)
#if defined (VAR_MPI)
      WRITE (*,*) '*** parallel run on',NMPROC,' nodes.'
#else
      WRITE (*,*) '*** serial run'
#endif
      WRITE(*,*)
      WRITE(*,*) ' ------->> symmetry:              ',desrep
      WRITE(*,*) 'ADC matrix dimension (# of rows): ',ladc
      WRITE(*,*) 'Size of main space:               ',nmain
      WRITE(*,*)
      WRITE(*,*) 'Number of Davidson roots:         ',nroots
      WRITE(*,*) 'Maximum Davidson subspace:        ',maxdavsp
      WRITE(*,*) 'Number of Macro iterations:       ',maxdavit
      WRITE(*,'(A,E14.5)') ' Eigenvalue convergence up to:   ',convthr
      WRITE(*,*)
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*) ' ------ Progress of calculation:'
      WRITE(*,*) ' --------------------------------------'
      WRITE(*,*)


!*******************************************************
! Read ADC matrix diagonal (is complex but imag part is zero)
!*******************************************************

      open(unit=iobase,file='ADCDGTMP', access='SEQUENTIAL', &
           form='UNFORMATTED',status='unknown')
      read(iobase) (adiag(i),i=1,ladc)
      close(iobase,status='DELETE')
      write(*,*) 'Diagonal elements read'

!*******************************************************
! Open main ADC matrix file. Remains open during cycles
!*******************************************************

#if defined (VAR_MPI)
      write (fileadc,'(A6,A1,I1)') nmbase,'.',MASTER
#else
      write (fileadc,'(A6)') nmbase
#endif
      open(iobase,FILE=fileadc,FORM='UNFORMATTED',STATUS='UNKNOWN')
      write(*,*) 'ADC file ',fileadc,' opened.'

      maciter = 1
      io_bold = iobase + 1

!*******************************************************
!  ENTRY POINT FOR MACROITERATIONS
!*******************************************************

 9999 CONTINUE
      write(*,*) 'Starting macroiteration #',maciter
      hij = czero

!*******************************************************
! Open direct access file for current Krylov space
!*******************************************************

      open(io_bold,file='DAVXBO',access='direct', &
           recl=ladc*16,status='unknown')
      write(*,*) 'Empty direct access file for complex ', &
                 'Krylov space opened.'

!*******************************************************
! Initialize Krylov space with (restart) vectors.
!*******************************************************

      call init_krylov_c(io_bold,desrep,nroots,maxdavsp,ladc, &
                         polprp_davreort)
      krydim = nroots
      krysav = 0

      WRITE(*,*)
      WRITE(*,'(3X,A,4X,A,7X,A,7X,A,5X,A,A)') 'root','dim',  &
                 'exc.energy','error','status','  time/root (s)'
      WRITE(*,'(70A1)') ('-',i=1,70)
      allconverged = .false.
!_______________________________________________________________
!|
!| Microiteration Davidson loop until all roots converged or the
!| Maximum Davidson space is reached.
!|
!|
      DO WHILE (allconverged.eqv..false.)

!*******************************************************
! Step X: Form subspace Hamiltonian, <b_I|H|b_J>
!         This is the only place where the large matrix is 
!         accessed --> parallelization.
!*******************************************************

      call cpu_time(time_kryl_s)


#if defined (VAR_MPI)
      do j=krysav+1,krydim
        call interface_mpi_bcast(SERVER_MATMUL,1,MASTER, &
                                 global_communicator)
        call interface_mpi_bcast(ladc,1,MASTER,global_communicator)
        read(io_bold,rec=j,iostat=ida_iostat) bv    

        mpiaux = real(bv)  ! communicate real part
        call interface_mpi_bcast(mpiaux,ladc,MASTER,global_communicator)
        mpiaux = aimag(bv)  ! communicate imaginary part
        call interface_mpi_bcast(mpiaux,ladc,MASTER,global_communicator)

        call xmvmul_sp_c(iobase,intbuf,ladc,bv,bx)  !master's own share --> bx

        mpires1 = rzero    ! this is the array that accepts the reduced contributions
        mpiaux = real(bx)
        call interface_mpi_reduce(mpiaux,mpires1,ladc, &
                       op_MPI_SUM,MASTER,global_communicator)
        mpires2 = rzero    ! this is the array that accepts the reduced contributions
        mpiaux = aimag(bx)
        call interface_mpi_reduce(mpiaux,mpires2,ladc, &
                       op_MPI_SUM,MASTER,global_communicator)

        do ix=1,ladc
          bwc(2*ix-1) = mpires1(ix)   !real part
          bwc(2*ix)   = mpires2(ix)   !imag part
        enddo

        if(qstack_push(qst_hboline,2*ladc,bwc).ne.2*ladc) stop 'QE push'
      enddo
#else
      do j=krysav+1,krydim
        read(io_bold,rec=j,iostat=ida_iostat) bv    
        call xmvmul_sp_c(iobase,intbuf,ladc,bv,bw)
        do ix=1,ladc
          bwc(2*ix-1) = real(bw(ix))   !real part
          bwc(2*ix)   = aimag(bw(ix))   !imag part
        enddo
        if(qstack_push(qst_hboline,2*ladc,bwc).ne.2*ladc) stop 'QE push'
      enddo
#endif


! form incremental contributions to subspace projected Hamiltonian
! call dot_product routine that directly accesses C stack.
! also written in C (necessarily) and correct krysav

      do i=1,krydim
        read(io_bold,rec=i) bv
        do ix=1,ladc
          bwc(2*ix-1) = real(bv(ix))    !real part
          bwc(2*ix)   = aimag(bv(ix))   !imag part
        enddo
        do j=krysav+1,krydim
          if(qstack_directx_c(qst_hboline,ladc,bwc,j,hdx).ne.ladc) &
             stop 'QE x_c'
          hij(i,j) = hdx
          hij(j,i) = dconjg(hdx)
        enddo
      enddo

      call cpu_time(time_kryl_e)
      time_cur = time_kryl_e - time_kryl_s
      write(*,'(A,F10.2)') 'Time for constructing Krylov projection', &
                           time_cur

      krysav = krydim

!*******************************************************
! Step X+Y: Diagonalize current subspace Hamiltonian
! hijcal is copy, contains eigenvectors after DSYEV call
!*******************************************************

      hijcal = hij
      call zheev('V','L',krydim,hijcal,maxdavsp,lambda, &
                 work,lwork,rwork,info)
      if(info.ne.0) then
        write(*,*) 'ZHEEV problem in complex Davidson, code:',info
        call quit()
      endif

!*******************************************************
! Step X+Y: calculate nroots residual vectors (kets!)
! Current dimension of Krylov space!
! resvec_c accesses C stack!
!*******************************************************

      call cpu_time(time_resi_s)
      call calc_resvec_c(io_bold,ladc,nroots,krydim, &
                         hijcal,lambda,res)
      call cpu_time(time_resi_e)
      time_cur = time_resi_e - time_resi_s
      write(*,'(A,F10.2)') 'Time for constructing residual vectors', &
                           time_cur
      write(*,*)

!*******************************************************
! Step X+Y: Enter loop over the roots.
! If convergence is not reached determine root-specific
! preconditioner, possibly enlarge space
!*******************************************************

      do kroot = 1,nroots

        call cpu_time(time_root_s)

        resnorm(kroot)=rzvecnorm(ladc,res(:,kroot))

        write(listout,'(A3,I3,A3,I5,6X,F10.7,6X,F10.7)') '   ',  &
            kroot,'   ',krydim,lambda(kroot),resnorm(kroot)


        if(resnorm(kroot).le.convthr) then

          converged(kroot) = .true.
          listout(47:48) =  ' *'
          listout(52:54) = '   '
          call cpu_time(time_root_e)
          time_cur = time_root_e - time_root_s
          write(*,'(A60,4X,F10.2)') listout,time_cur
          cycle

        else

          converged(kroot) = .false.
          listout(47:48) =  ' o'

          call calc_precon_c(ladc,adiag,lambda(kroot),res(:,kroot),prc)
          prenorm1=rzvecnorm(ladc,prc)

          call calc_preort_c(io_bold,ladc,krydim,prc,pro)
          prenorm2=rzvecnorm(ladc,pro)

          if( (prenorm2/prenorm1).ge.prethr) then

            oviol = .false.    !  preconditioner has large enough norm
            pro = pro/dcmplx(prenorm2,rzero)
            krydim = krydim + 1
            write(io_bold,rec=krydim) pro  !extend trial space
            listout(52:54) = ' XT'

          else

            oviol = .true.    !  preconditioner norm is critical
            pro = pro/dcmplx(prenorm2,rzero)
            krydim = krydim + 1
            write(io_bold,rec=krydim) pro  !extend trial space
            listout(52:54) = ' NX'

          endif
          call cpu_time(time_root_e)
          time_cur = time_root_e - time_root_s
          write(*,'(A60,4X,F10.2)') listout,time_cur

        endif

      enddo ! kroot
      write(*,*)
!
!  check if all roots converged.
!
      allconverged = .true.
      do kroot=1,nroots
        if(converged(kroot).eqv..false.) then
          allconverged = .false.
        endif
      enddo
!
!  check if orthogonality problems pertain after last root
!
      if(oviol) then
        write(*,*) ' ****            Warning                  ****'
        write(*,*) ' **** Orthogonality of preconditioner is weak.'
        write(*,*) ' **** Check convergence of final results.'
        write(*,*) ' ****                                     ****'
        write(*,*)
      endif
!
!  check if dimensionality extension is still in allowed range.
!
      if( (krydim + nroots).gt.maxdavsp) then
        write(*,*) 'Maximum Davidson space reached. Exiting'
        exit
      endif


      END DO   !Microiteration
!|
!|
!|
!|       end of microiteration loop
!|
!|______________________________________________________


!#############################################################
!############ Check break conditions for macroiterations #####
!############ and clear up things                        #####
!#############################################################

      write(*,*) 'Macroiteration',maciter,' ended with'
      write(*,*) 'Krylov space dimension:',krydim

! clear stack (essential! large memory demands)

      i=qstack_drop(qst_hboline)

! write out eigenvectors or restart vectors from B_old file.
! If not all roots converged we collapse subspace again to
! krydim = nroots.

      if(allconverged) then
        write(*,*) 'All roots converged. Writing eigenvectors'
      else
        write(*,*) 'Nonconverged roots remain. Collapsing subspace.'
      endif
!     call write_davevecs_c(io_bold,desrep,nroots,nmain,ladc, &
!                           krydim,lambda,resnorm,hijcal,ortho_cur)
      call write_davevecs_c(io_bold,desrep,nroots,nmain,ladc, &
                            krysav,lambda,resnorm,hijcal,ortho_cur)

! in both cases we have to get rid off the large B old file

      close(io_bold,status='DELETE')

! check if Krylov space becomes insufficent. This is indicated by 
! deteriorating orthogonality of the (long) eigenvectors.

      if(abs(ortho_cur).gt.1.0d-10) then
        call ortho_warn()
        maciter = maxdavit  !forcing maximum macro iterations, leads to loop break
      endif

! now check if we should enter the next macroiteration.
! If yes we start over

      if(.not.allconverged) then
        maciter = maciter + 1
        if(maciter.le.maxdavit) then
          write(*,*) 'Entering macroiteration',maciter
          Goto 9999
        else
          write(*,*) 'User limit of macroiterations reached.'
          write(*,*) 'Check convergence of remaining roots and'
          write(*,*) 'be careful with further eigenvector usage.'
        endif
      else
          write(*,*) 'Overall Davidson convergence reached.'
          write(*,*) 'Macroiterations used:',maciter
      endif
!
!  deallocate complex arrays
!
      deallocate(adiag)
      deallocate(bv)
      deallocate(bw)
      deallocate(bx)
      deallocate(hij)
      deallocate(hijcal)
      deallocate(prc)
      deallocate(pro)
      deallocate(res)

!  deallocate real arrays

      deallocate(lambda)
      deallocate(work)
      deallocate(rwork)
      deallocate(mpiaux)
      deallocate(mpires1)
      deallocate(mpires2)
      deallocate(bwc)
      deallocate(converged)

!  close and delete ADC matrix file

      close(iobase,status='delete')    !delete ADC matrix

      call cpu_time(time_totl_e)

      write(*,*)
      write(*,'(70A1)') ('-',i=1,70)
      write(*,*) 'Davidson finished.'
      time_cur = time_totl_e - time_totl_s
      write(*,'(A,F10.2)') 'Total time spent in Davidson: ',time_cur
      write(*,'(70A1)') ('-',i=1,70)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE init_krylov_c(io_bold,desrep,nroots,maxdavsp,ladc,dao)

      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!   This routine checks the existence of some restart files and tries to fit
!   found restart vectors into the vector array. The order of importance is:
!
!   ADCEVECS.XX
!   DAVSTART.XX
!   Unity matrix
!
!------------------------ calling variables -------------------------
!
      integer                    :: io_bold,desrep,nroots,maxdavsp,ladc
      logical                    :: dao
!
!------------------------ local variables -------------------------
!
      logical                                :: isthere
      integer                                :: fileavail,io
      integer                                :: i,j,k,ix,nsiz
      integer                                :: ida_iostat
      Character*11                           :: filename_d,filename_e
      Real*8                                 :: ortho,a1
      complex*16,allocatable,dimension(:)        :: qm,ec
      complex*16,allocatable,dimension(:,:)      :: bstart
      integer                                :: xnroots,xnmain,xladc
      integer                                :: unroots,uladc
      integer                                :: jdummy
      Real*8                                 :: rdummy
      character*8                            :: hp_name,ev_name
      Real*8,parameter                       :: rzero=0.0d0
      complex*16                             :: cax
      complex*16,parameter                   :: czero=(rzero,rzero)
      complex*16,parameter                   :: cone=(1.0d0,rzero)
!
!------------------------ interface -------------------------
!
      Interface

        function get_file_unit()
           Integer    :: get_file_unit
        end function

        function pp_zdotc(i1,ca1,ca2)
          double complex                   :: pp_zdotc
          INTEGER                          :: i1
          double complex, dimension(:)     :: ca1,ca2
        end function

      End interface
!
!------------------------ executable code -------------------------
!
      fileavail = 0
      hp_name = 'DAVSTART'
      ev_name = 'ADCEVECS'
!
! check availability of DAVSTART.XX file
!
      IF(desrep.GT.9) THEN
        WRITE(filename_d,'(A8,A1,I2)') hp_name,'.',desrep
      ELSE
        WRITE(filename_d,'(A8,A2,I1)') hp_name,'.0',desrep
      ENDIF
      inquire(file=filename_d,exist=isthere)
      if(isthere) fileavail = fileavail + 1
!
! check availability of ADCEVECS.XX file
!
      IF(desrep.GT.9) THEN
        WRITE(filename_e,'(A8,A1,I2)') ev_name,'.',desrep
      ELSE
        WRITE(filename_e,'(A8,A2,I1)') ev_name,'.0',desrep
      ENDIF
      inquire(file=filename_e,exist=isthere)
      if(isthere) fileavail = fileavail + 1
!
!  read corresponding file or pack in unity matrix
!
      io = get_file_unit()

      Select case(fileavail)

!------------------------------------------------------------

      case(2)                 !  ADCEVECS is present

      write(*,*) 'Using ',filename_e,' as restart file.'
      OPEN(io,file=filename_e,access='sequential', &
           status='unknown',form='unformatted')
      READ(io) xnroots
      READ(io) xnmain
      READ(io) xladc
      write(*,*) 'There are',nroots,' vectors on the file'
      write(*,*) 'Vectors have length',xladc

      unroots = nroots
      if(nroots.gt.xnroots) then
        write(*,*) 'Only',xnroots,' available and used.'
        unroots = xnroots
      endif

      uladc = ladc
      if(ladc.ne.xladc) then
        write(*,*) 'Vector length does not match.'
        write(*,*) 'It will be adjusted but convergence can suffer!'
        uladc = min(ladc,xladc)
      endif

      allocate(qm(ladc))
      allocate(ec(xladc))
      Do k=1,unroots
        read(io) jdummy
        read(io) rdummy
        read(io) rdummy
        read(io) (ec(ix),ix=1,xladc)
        qm = czero
        qm(1:uladc) = ec(1:uladc)  !rest is zero if xladc .lt. ladc
        write(io_bold,rec=k,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call  &
           quit('Error writing direct access file for bold')
      Enddo
      write(*,*) 'Wrote',unroots,' start vectors into initial block.'

! filling the remaining start vectors with diagonal one

      Do k=unroots+1,nroots
        write(*,*) 'Writing defective vector #',k
        qm = czero
        qm(k) = cone
        write(io_bold,rec=k,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call  &
           quit('Error writing direct access file for bold')
      Enddo
      deallocate(ec)
      deallocate(qm)

      if(dao) write(*,*) 'Reorthogonalizing ADCEVECS.'
      if(dao) call reortho_ext_c(io_bold,ladc,nroots)

      Close(io)

!------------------------------------------------------------

      case(1)                 !  DAVSTART is present

      write(*,*) 'Using ',filename_d,' as restart file.'
      OPEN(io,file=filename_d,access='sequential',  &
           status='unknown',form='unformatted')
      READ(io) nsiz
      WRITE(*,*) 'Found',nsiz,' Start vectors.'

      If(2*nroots.gt.maxdavsp) then
         Call Quit('Increase Davidson space!')
      Endif

      If(nroots.gt.nsiz) then
        Write(*,*) 'You request more states than there',  &
                   ' are in the Main space.'
        Write(*,*) 'Adapting requested roots',nroots,' to',nsiz
        nroots = nsiz
      Endif

      allocate(bstart(nsiz,nroots))
      bstart = czero

      do i=1,nroots
        Read(io) (bstart(j,i),j=1,nsiz)
      enddo
      write(*,*) nroots,' start vectors read.'
      write(*,*) 'Length of start vectors:',nsiz

! check orthogonality of the start vectors

      ortho = 0.0d0
      Do i=1,nroots
      Do j=1,nroots
        cax = pp_zdotc(nsiz,bstart(:,i),bstart(:,j))
        ortho = ortho + abs(cax)
      Enddo
      Enddo
      ortho = ortho - dble(nroots)
      WRITE(*,*) 'Orthogonality:',ortho
!
! now make start vectors the initial Davidson vectors for the first 
! macro iteration
! 
      allocate(qm(ladc))
      Do i=1,nroots
        qm(1:nsiz) = bstart(1:nsiz,i)
        qm(nsiz+1:ladc) = 0.0d0
        write(io_bold,rec=i,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call  &
           quit('Error writing direct access file for bold')
      Enddo

      deallocate(qm)
      deallocate(bstart)

      Close(io)

!------------------------------------------------------------

      case(0)                 !  no file. Unity matrix

      write(*,*) 'No restart file found. Using unity matrix.'

      allocate(qm(ladc))
      qm = czero
      Do i=1,nroots
        qm(i) = cone
        write(io_bold,rec=i,iostat=ida_iostat) qm
        if(ida_iostat.ne.0) call  &
           quit('Error writing direct access file for bold')
        qm(i) = czero
      Enddo
      deallocate(qm)



!------------------------------------------------------------

      End select

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine reortho_ext_c(io_bold,ladc,nroots)
!
      IMPLICIT NONE
!
!------------------------ calling variables -------------------------
!
      INTEGER                        :: io_bold,ladc,nroots
!
!------------------------ local variables -------------------------
!
      complex*16, allocatable, dimension(:)   :: x,ax
      complex*16                              :: ca5
      real*8                                  :: rn
      INTEGER                                 :: i,k,ida_iostat
!
!------------------------ interfaces -------------------------------
!
      interface

        FUNCTION pp_zdotc(i1,ca1,ca2)
          double complex                   :: pp_zdotc
          INTEGER                          :: i1
          double complex, dimension(:)     :: ca1,ca2
        END FUNCTION

        FUNCTION rzvecnorm(i1,ca1)
          real*8                           :: rzvecnorm
          integer                          :: i1
          double complex, dimension(:)     :: ca1
        END FUNCTION

      end interface
!
!---------------Executable ----------------------------------------
!
! reorthogonalize nroots vectors stored columnwise in the array v
!
      allocate(x(ladc))
      allocate(ax(ladc))

      Do k=1,nroots
        read(io_bold,rec=k,iostat=ida_iostat) x
        Do i=1,k-1
          read(io_bold,rec=i,iostat=ida_iostat) ax
          ca5 = pp_zdotc(ladc,x(:),ax(:))
          x = x - ca5*ax(:)
        Enddo
        rn = rzvecnorm(ladc,x)
        ax = x/dcmplx(rn,0.0d0)
        write(io_bold,rec=k,iostat=ida_iostat) ax
      Enddo

      deallocate(x)
      deallocate(ax)
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine calc_resvec_c(iu,ladc,nroots,krydim, &
                               cji,lambda,res)
!
!  compute residual vectors res, filled within this routine.
!  Addresses the corresponding H*b_old vector from stack.
!
      use qstack

      implicit none
!
! formal arguments
!
      integer                               :: iu
      integer                               :: ladc
      integer                               :: nroots
      integer                               :: krydim
      complex*16,dimension(:,:)             :: cji
      real*8,dimension(:)                   :: lambda
      complex*16,dimension(:,:)             :: res
!
! local variables
!
      integer                               :: ix,j,kroot
      Integer                               :: qst_hboline

      real*8,parameter                      :: rzero = 0.0d0
      real*8,allocatable,dimension(:)       :: hrc

      complex*16,parameter                  :: czero = (rzero,rzero)
      complex*16,allocatable,dimension(:)   :: aux
      complex*16,allocatable,dimension(:)   :: bv
      complex*16,allocatable,dimension(:)   :: hbo

#include  "polprp_stacklines.h"
!
! executable code
!
      qst_hboline = HBO_STACKLINE
!
! form sum_j (C_ji*(H*b_j - lambda_i*b_j))
!
      allocate(hrc(2*ladc))
      hrc = rzero

      allocate(aux(ladc))
      allocate(bv(ladc))
      allocate(hbo(ladc))
!
! first loop over Krylov space (minimize disk reads!)
! Inner loop: over nroots
!
      res = czero
      Do j=1,krydim
        read(iu,rec=j) bv
        if(qstack_peekn(qst_hboline,hrc,j).ne.2*ladc) stop 'QE'

! reconstruct complex vector from real stack line

        do ix=1,ladc
          hbo(ix) = dcmplx(hrc(2*ix-1),hrc(2*ix))
        enddo
        Do kroot=1,nroots
          aux(:) = hbo - dcmplx(lambda(kroot),rzero)*bv(:)
          aux(:) = aux(:)*cji(j,kroot)
          res(:,kroot) = res(:,kroot) + aux(:)
        Enddo  !kroot
      Enddo  !j

      deallocate(hrc)
      deallocate(aux)
      deallocate(bv)
      deallocate(hbo)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine calc_precon_c(ladc,adiag,lambdai,res,prc)
!
! compute preconditioner
!
      implicit none
!
! formal arguments
!
      integer                               :: ladc
      complex*16,dimension(:)               :: adiag
      real*8                                :: lambdai
      complex*16,dimension(:)               :: res
      complex*16,dimension(:)               :: prc
!
! local variables
!
      integer                               :: j
      real*8                                :: a1
!
! executable code
!

      Do j=1,ladc
        a1=lambdai - real(adiag(j))
        if(dabs(a1).lt.1.0e-03) then
          prc(j) = (100.0d0,0.0d0) * res(j)
        else
          prc(j)=res(j)/dcmplx(a1,0.0d0)
        endif
      Enddo

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine calc_preort_c(iu,ladc,krydim,prc,pro)
!
! orthogonalize preconditioner to current subspace
! prc goes in, unchanged
! pro goes out, assigned
!
      implicit none
!
! formal arguments
!
      integer                               :: iu
      integer                               :: ladc
      integer                               :: krydim
      complex*16,dimension(:)               :: prc
      complex*16,dimension(:)               :: pro

      Interface

        function pp_zdotc(i1,ca1,ca2)
          double complex                   :: pp_zdotc
          INTEGER                          :: i1
          double complex, dimension(:)     :: ca1,ca2
        end function

      End interface
!
! local variables
!
      integer                               :: j
      complex*16                            :: a1
      complex*16,allocatable,dimension(:)   :: aux
!
! executable code
!
      allocate(aux(ladc))


      pro = prc
      Do j = 1,krydim
        read(iu,rec=j) aux
        a1=pp_zdotc(ladc,aux,prc)
        pro(:) = pro(:) - a1*aux(:)
      Enddo

      deallocate(aux)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE write_davevecs_c(io_bold,desrep,nroots,nmain,ladc, &
                                  krydim,lambda,resnorm,cji,ortho)

      implicit none
!
!---------------Description--------------------------------------------
!
! eigenvectors (complex) are written out.
!
!------------------------ calling variables -------------------------
!
      integer                         :: io_bold
      integer                         :: desrep,nroots
      integer                         :: nmain,ladc,krydim
      real*8,dimension(:)             :: lambda
      real*8,dimension(:)             :: resnorm
      complex*16,dimension(:,:)       :: cji
      real*8                          :: ortho
!
!------------------------ local variables -------------------------
!
      integer                                :: i,j,k
      integer                                :: iu,iocnf
      integer                                :: irecl,istat
      integer,allocatable,dimension(:)       :: ievec
      real*8                                 :: a1,enorm
      real*8,parameter                       :: autoev = 27.2113957D0
      real*8,parameter                       :: rzero=(0.0d0)
      real*8,allocatable,dimension(:)        :: sx
      complex*16                             :: cax
      complex*16,parameter                   :: czero=(rzero,rzero)
      complex*16,allocatable,dimension(:)    :: aux
      complex*16,allocatable,dimension(:)    :: bv
      complex*16,allocatable,dimension(:,:)  :: evc
      character*6,parameter                  :: filecnf='XCONFG'
      character*11                           :: davevcname
      character*70                           :: field
      logical                                :: isthere

!
!------------------------ local variables -------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

        FUNCTION rzvecnorm(i1,ca1)
          real*8                           :: rzvecnorm
          INTEGER                          :: i1
          complex*16, dimension(:)         :: ca1
        END FUNCTION

        FUNCTION pp_zdotc(i1,ca1,ca2)
          double complex                   :: pp_zdotc
          INTEGER                          :: i1
          double complex, dimension(:)     :: ca1,ca2
        END FUNCTION

      end interface
!
!------------------------ executable code -------------------------
!
      allocate(ievec(ladc))
      ievec = 0
      allocate(sx(ladc))
      sx = rzero
      allocate(aux(ladc))
      aux = czero
      allocate(bv(ladc))
      bv = czero
      allocate(evc(ladc,nroots))
      evc = czero

! *** calculate eigenvector approximations from trial space
! *** this is done with a temporary array in order to minimize
! *** read operations from direct access file

      Do j = 1,krydim
        read(io_bold,rec=j) aux
        Do k = 1,nroots
          evc(:,k) = evc(:,k) + cji(j,k)*aux
        Enddo
      Enddo

! *** open file

      iu = get_file_unit()
!     write(*,*) 'Unit for writing eigenvectors:',iu
      if(desrep.GT.9) THEN
        write(davevcname,'(A8,A1,I2)') 'ADCEVECS','.',desrep
      else
        write(davevcname,'(A8,A2,I1)') 'ADCEVECS','.0',desrep
      endif
      write(*,*) 'Writing out eigenvectors to file ',davevcname
      write(*,*)
      open(iu,file=davevcname,access='sequential', &
           status='unknown',form='unformatted')

! *** write header information

      write(iu) nroots
      write(iu) nmain
      write(iu) ladc

! *** loop over the roots

      Do k = 1,nroots

        write(iu) k
        write(iu) lambda(k)

! *** compute pole strength of kth eigenvector and write it out

        a1 = 0.0d0
        do j=1,nmain
          a1 = a1 + real(dconjg(evc(j,k))*evc(j,k))
        enddo
        write(iu) a1
        write(iu) (evc(j,k),j=1,ladc)

      Enddo   !k. 
      close(iu)

!   We check orthonormality of all eigenvectors

      ortho = 0.0d0
      do i=1,nroots
        do j=i,nroots
          ortho = ortho + abs(pp_zdotc(ladc,evc(:,i),evc(:,j)))
        enddo
      enddo
      ortho = ortho - dble(nroots)
      write(*,*) 'Orthogonality of written vectors:',ortho
      write(*,*)

!__________________________________________________________
!|
!|
!|                We perform state analysis
!|
!|

!   first check validity of configuration file, otherwise no reasonable 
!   analysis possible. get record length of configuration file and open it

      Inquire(file=filecnf,exist=isthere,recl=irecl,iostat=istat)
      if(.not.isthere) then
        write(*,*) '******  warning  *******'
        write(*,*) 'No configuration data available. State'
        write(*,*) 'analysis not possible.'
        write(*,*) 'ADCSTATE.XX not created.'
        write(*,*) '******  warning  *******'
        deallocate(ievec)
        deallocate(aux)
        deallocate(bv)
        deallocate(evc)
        return
      endif

      irecl = 70
      iocnf = get_file_unit()
      OPEN(iocnf,file=filecnf,ACCESS='DIRECT',RECL=IRECL, &
           STATUS='UNKNOWN')
!     write(*,*) 'Config file opened with',iocnf

      iu = get_file_unit()
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') 'ADCSTATE','.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') 'ADCSTATE','.0',desrep
      ENDIF
      OPEN(iu,file=davevcname,access='sequential', &
           status='unknown',form='formatted')
!     write(*,*) 'State analysis file opened with',iu

! now loop over states and collect most contributing configurations

      Do k = 1,nroots
        a1 = rzero
        do j=1,nmain
          a1 = a1 + real(dconjg(evc(j,k))*evc(j,k))
        enddo

        Write(iu,'(A14,X,I4)') 'Excited state:',k
        Write(iu,'(A19)') '-----------------'
        Write(iu,'(A12,F10.6,A4,F10.6,A8,F10.6,A5,F10.6,A5,I4)')  &
         'Exc. energy:', &
         lambda(k),' au ',lambda(k)*autoev, &
         ' eV  PS:',a1,' Err:',resnorm(k),' Sym:',desrep
        Write(iu,'(A28,F15.6)') 'Norm^2 of double exc. contr:', &
               1.0d0 - a1
        Write(iu,'(A)') ' '
!
! determine 20 most relevant configurations
! sort the corresponding eigenvector with respect to coeffs.
!
        do j=1,ladc
          ievec(j)=j
        enddo
        sx=-abs(evc(:,k))
        Call DISORT(ladc,sx,ievec)  !sorts from smallest to largest
        Do j=1,min(20,ladc)
          Read(iocnf,rec=ievec(j)) field
          cax = evc(ievec(j),k)
          Write(iu,'(A70,2F15.6,3X,F15.6)') field,cax, &
                  real(dconjg(cax)*cax)
        Enddo
        Write(iu,'(A17)') '-----------------'
        Write(iu,'(A)') ' '

      Enddo ! nroots
!|
!|
!|_______________________________________________
!

! cleanup. Close files, deallocate arrays

      close(iu)
      close(iocnf)
      deallocate(ievec)
      deallocate(sx)
      deallocate(aux)
      deallocate(bv)
      deallocate(evc)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
