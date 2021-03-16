#if defined (VAR_MPI)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE mpi_task_server(iobase,intbuf)

      use interface_to_mpi

      IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
! This is a server module that executes all slave node tasks required in
! the parallel POLPRP run. The services are controlled by a corresponding
! broadcast code. The codes are contained in the polprp_servercodes.h
! file. A termination signal triggers exit from the server loop.
!
!--------------- Includes--------------------------------------
!
#include  "polprp_servercodes.h"
#include  "../relccsd/symm.inc"
#include  "../relccsd/ccpar.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      Integer                           :: iobase,intbuf

      Logical                           :: deb = .false.

      Integer                           :: taskid = 0   ! defines the server's task
      REAL*8,allocatable,dimension(:)   :: buf
      logical                           :: done
      INTEGER                           :: irp,istart,nint,mint
      INTEGER                           :: ialloc,locoff,istat
      INTEGER                           :: efflen,n1,n2
      INTEGER                           :: m,n,k
!
! variables for partial precontraction
!
      real*8,allocatable,dimension(:)       :: vvoo_in,vvoo_ou
!
! variables for matrix distribution
!
      Integer                               :: i,ixx,jdummy
      REAL*8,allocatable,dimension(:)       :: bufr,bufi
      Complex*16,allocatable,dimension(:)   :: cbuf
      Integer,allocatable,dimension(:)      :: ioi,ioj
      Integer                               :: nbufscur,nact
      Character*10                          :: locfn
      Character*6                           :: nmbase='PP_ELE'
      Logical                               :: isthere
!
! variables for matrix multiplication
!
      Complex*16,allocatable,dimension(:)   :: cveco,cvecn
      REAL*8,allocatable,dimension(:)       :: veco,vecn
      REAL*8,allocatable,dimension(:)       :: vecr,veci
      REAL*8,allocatable,dimension(:)       :: dumarr
      Integer                               :: ladc
      Integer                               :: irow,icol
      Integer                               :: nbufsloc
      Complex*16                            :: a
      Real*8                                :: b

      complex*16, parameter  :: czero=(0.0d0,0.0d0)
      complex*16, parameter  :: cone =(1.0d0,0.0d0)
      real*8, parameter      :: rzero=0.0d0

     
!
!---------------Interface area --------------------------------------
!
      interface

        SUBROUTINE do_partial_precon(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

      end interface
!
!---------------Executable code--------------------------------------
!
!----------------------------------------------------------------
!        determine if user wants to trace MPI communication
!----------------------------------------------------------------

      If(SERVER_TRACE.ne.0) deb = .true.

!----------------------------------------------------------------
      DO       ! ---------- INFINITE SERVER MAIN LOOP -----------
!----------------------------------------------------------------

        call interface_mpi_BCAST(taskid,1,MASTER,global_communicator)
        if(deb) write(*,*) 'Node',MYPROC,' has received taskid',taskid

        select case(taskid)


!----------------------------------------------------------------
!-----  TASK: say hello
!----------------------------------------------------------------


          case(SERVER_HELLO)
          write(*,*) 'Node',MYPROC,' confirms existence.'


!----------------------------------------------------------------
!-----  TASK: terminate server loop
!----------------------------------------------------------------


          case(SERVER_DONE)
          write(*,*) 'Node',MYPROC,' terminates service.'
          exit   ! terminates infinite do loop


!----------------------------------------------------------------
!-----  TASK: perform local share of precontraction
!----------------------------------------------------------------


          case(SERVER_PRECON)
          write(*,*) 'Node',MYPROC,' starts partial precontraction.'

          n1 = ivvoott(nrep+1)
          allocate(vvoo_in(n1*rcw))
          allocate(vvoo_ou(n1*rcw))
          vvoo_in = 0.0d0
          vvoo_ou = 0.0d0
          call interface_mpi_recv(vvoo_in,n1*rcw,MASTER,MSGN,
     &                            global_communicator)
!
!             do the partial (local) precontraction
!
          call do_partial_precon(vvoo_in,vvoo_ou)

          call interface_mpi_send(vvoo_ou,n1*rcw,MASTER,MSGN,
     &                            global_communicator)
          if(deb) write(*,*) 'Node',myproc,' has sent output ints.'
          deallocate(vvoo_in)
          deallocate(vvoo_ou)


!----------------------------------------------------------------
!-----  TASK: provide local share of VVVV integrals in irep IRP
!----------------------------------------------------------------


          case(SERVER_VVVV)
          call interface_mpi_recv(irp,1,MASTER,MSGN,global_communicator)
          if(deb) write(*,*) 'Node',MYPROC,' obtained irp:',irp
          if(irp.lt.1.or.irp.gt.nrep) then
            write(*,*)  '**** Error on node',MYPROC,':'
            write(*,*)  '**** irep number out of bounds,'
            call quit('IREP error')
          endif
          efflen = NVVT(irp) * NVVT(irp) !efflen is max number of communicated integrals
                                         !in this irep
          if(efflen.gt.0) then
            allocate(buf(efflen*RCW),stat = ialloc)
            if(ialloc.ne.0)
     &          call quit('Allocation problem in mpi_task_server!')
            buf = 0.0d0

            istart = 0

            if(deb) write(*,*) 'Istart before:',istart
            MINT = NVVT(IRP)  !maximum number of possible batches (last indices count)
            CALL GETVVVV (irp,istart,nint,done,buf,mint)

            if(deb) write(*,*) 'Recovered ',nint,' batches from',
     &        mint,' requested ones on node',MYPROC,' istart now',
     &        istart
!
! NINT contains the actual number of batches on *this* node.
!
            if (.not.done) then
               write(*,*) 'irrep',irp,' unfinished on node',myproc
               call quit('Error reading VVVV ints in mpi_task_server.')
            endif
            if(deb) write(*,*) '-----> max length:',efflen,
     &            ' actual len:',NVVT(irp) * NINT

            efflen = NVVT(irp) * NINT
          endif

          call interface_mpi_send(efflen*rcw,1,MASTER,MSGN,
     &                            global_communicator)
          locoff = IDIST(3,1,irp)
          call interface_mpi_send(locoff*rcw,1,MASTER,MSGN,
     &                            global_communicator)

          if(efflen.gt.0) then
            call interface_mpi_send(buf,efflen*rcw,MASTER,MSGN,
     &                            global_communicator)
          endif
          if (allocated(buf))  deallocate(buf)


!----------------------------------------------------------------
!-----  TASK: provide local share of VOVV integrals in irep IRP
!----------------------------------------------------------------


          case(SERVER_VOVV)
          call interface_mpi_recv(irp,1,MASTER,MSGN,global_communicator)
          if(deb) write(*,*) 'Node',MYPROC,' obtained irp:',irp
          if(irp.lt.1.or.irp.gt.nrep) then
            write(*,*)  '**** Error on node',MYPROC,':'
            write(*,*)  '**** irep number out of bounds,'
            call quit('IREP error')
          endif
          efflen = NVO(irp) * NVVT(irp)  !efflen is max number of communicated integrals
                                       !in this irep
          if(efflen.gt.0) then
            allocate(buf(efflen*RCW),stat = ialloc)
            if(ialloc.ne.0)
     &        call quit('Allocation problem in mpi_task_server!')
            buf = 0.0d0

            istart = 0

            MINT = NVVT(IRP)  !maximum number of possible batches
            CALL GETVOVV (irp,istart,nint,done,buf,mint)

            if(deb) write(*,*) 'Recovered ',nint,' batches on',MYPROC
!
! NINT contains the actual number of batches on *this* node.
!
            if (.not.done) then
              write(*,*) 'irrep',irp,' unfinished on node',myproc
              call quit('Error reading VOVV ints in mpi_task_server')
            endif

            efflen = NVO(irp) * NINT
          endif
!
          call interface_mpi_send(efflen*rcw,1,MASTER,MSGN,
     &                            global_communicator)
          locoff = IDIST(3,2,irp)
          call interface_mpi_send(locoff*rcw,1,MASTER,MSGN,
     &                            global_communicator)

          if(efflen.gt.0) then
            call interface_mpi_send(buf,efflen*rcw,MASTER,MSGN,
     &                              global_communicator)
          endif

          if (allocated(buf))  deallocate(buf)


!----------------------------------------------------------------
!-----  TASK: store fraction of ADC matrix on local node
!----------------------------------------------------------------


          case(SERVER_MATDIS)
!
! create local file
! file must remain after construction. However,
! if a partial ADC matrix exists at this point it must be deleted
! because in each symmetry this matrix is constructed exactly once.
!
          IF (MYPROC .LT. 10) THEN
            WRITE (locfn,'(A6,A1,I1)') nmbase,'.',MYPROC
          ELSEIF (MYPROC .LT. 100) THEN
            WRITE (locfn,'(A6,A1,I2)') nmbase,'.',MYPROC
          ELSEIF (MYPROC .LT. 1000) THEN
            WRITE (locfn,'(A6,A1,I3)') nmbase,'.',MYPROC
           ELSE
          CALL QUIT("NMPROC.GT.1000! EXTEND mpi_task_server routine!")
          ENDIF
          if(deb) write(*,*) 'Node',myproc,' opened ',
     &              locfn,' for writing.'
          inquire(file=locfn,exist=isthere)
          if(isthere) then
            open(unit=iobase,file=locfn, access='SEQUENTIAL',
     &           form='UNFORMATTED',status='old')
            close(iobase,status='delete')
          endif
          open(unit=iobase,file=locfn, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='new')
!
! create local buffers, freed afterwards
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

          jdummy = 0
!
! enter loop for receiving chunks of ADC matrix
!
          call interface_mpi_recv(nbufscur,1,MASTER,MSGN,
     &                           global_communicator)
          if(deb) write(*,*) 'Slave',MYPROC,' will receive',nbufscur,
     &        ' batches from master.'
          write(iobase,ERR=889) nbufscur
          do i=1,nbufscur
            if(carith) then
              call interface_mpi_recv(bufr,intbuf,MASTER,MSGN,
     &                                global_communicator)
              call interface_mpi_recv(bufi,intbuf,MASTER,MSGN,
     &                                global_communicator)
              do ixx=1,intbuf
                cbuf(ixx) = dcmplx(bufr(ixx),bufi(ixx))
              enddo
            else
              call interface_mpi_recv(buf,intbuf,MASTER,MSGN,
     &                                global_communicator)
            endif
            call interface_mpi_recv(ioi,intbuf,MASTER,MSGN,
     &                              global_communicator)
            call interface_mpi_recv(ioj,intbuf,MASTER,MSGN,
     &                              global_communicator)
            call interface_mpi_recv(nact,1,MASTER,MSGN,
     &                              global_communicator)
            if(carith) then
              write(iobase,ERR=889) (cbuf(ixx),ixx=1,intbuf),
     &                              (ioi(ixx),ixx=1,intbuf),
     &                              (ioj(ixx),ixx=1,intbuf),
     &                              nact,jdummy
            else
              write(iobase,ERR=889) (buf(ixx),ixx=1,intbuf),
     &                              (ioi(ixx),ixx=1,intbuf),
     &                              (ioj(ixx),ixx=1,intbuf),
     &                              nact,jdummy
            endif
          enddo
          flush(iobase)
          close(iobase)
          if(deb) write(*,*) 'Node',myproc,' has closed local ADC file.'
!
! deallocating local buffers
!
          if(carith) then
            deallocate(cbuf,bufr,bufi)
          else
            deallocate(buf)
          endif
          deallocate(ioi,ioj)


!----------------------------------------------------------------
!-----  TASK: perform partial matrix multiplication (r+c)
!----------------------------------------------------------------


          case(SERVER_MATMUL)
          call interface_mpi_BCAST(ladc,1,MASTER,global_communicator)
          if(deb) write(*,*) 'Node',myproc,' has obtained dim',ladc
!
! open local matrix file
!
          if (myproc .lt. 10) then
             write (locfn,'(a6,a1,i1)') nmbase,'.',MYPROC
          elseif (myproc .lt. 100) then
             write (locfn,'(a6,a1,i2)') nmbase,'.',MYPROC
          elseif (myproc .lt. 1000) then
             write (locfn,'(a6,a1,i3)') nmbase,'.',MYPROC
          else
             call quit("nmproc.gt.1000!")
          ENDIF

          open(unit=iobase,file=locfn, access='SEQUENTIAL',
     &         form='UNFORMATTED',status='old')
          if(deb) write(*,*) 'Local file ',locfn,' opened for reading.'
!
! complex case
!
          if(carith) then

            allocate(cveco(ladc))
            allocate(cvecn(ladc))
            allocate(vecr(ladc))
            allocate(veci(ladc))
            allocate(dumarr(ladc))
            allocate(ioi(intbuf))
            allocate(ioj(intbuf))
            allocate(cbuf(intbuf))
            cveco = czero
            cvecn = czero
            vecr = rzero
            veci = rzero
            dumarr = rzero
            ioi = 0
            ioj = 0
            cbuf = czero

            call interface_mpi_BCAST(vecr,ladc,MASTER,
     &                               global_communicator)
            call interface_mpi_BCAST(veci,ladc,MASTER,
     &                               global_communicator)
            do i=1,ladc
              cveco(i) = dcmplx(vecr(i),veci(i))
            enddo

            cvecn = czero
            rewind(iobase)
            read(iobase) nbufsloc
            do i = 1,nbufsloc
              read(iobase,err=889) (cbuf(ixx),ixx=1,intbuf),
     &                             (ioi(ixx),ixx=1,intbuf),
     &                             (ioj(ixx),ixx=1,intbuf),
     &                             nact,jdummy
              do k = 1, nact
                irow = ioi(k)
                icol = ioj(k)
                a    = cbuf(k)
                cvecn(irow) = cvecn(irow) + cveco(icol) * a
                if(irow.ne.icol) cvecn(icol) = cvecn(icol) +
     &                           cveco(irow)*dconjg(a)
              enddo
            enddo
            vecr = real(cvecn)
            veci = aimag(cvecn)
            call interface_mpi_reduce(vecr,dumarr,ladc,
     &                 op_MPI_SUM,MASTER,global_communicator)
            call interface_mpi_reduce(veci,dumarr,ladc,
     &                 op_MPI_SUM,MASTER,global_communicator)

            deallocate(cveco)
            deallocate(cvecn)
            deallocate(vecr)
            deallocate(veci)
            deallocate(dumarr)
            deallocate(ioi)
            deallocate(ioj)
            deallocate(cbuf)

          else
!
! real case
!
            allocate(veco(ladc))
            allocate(vecn(ladc))
            allocate(dumarr(ladc))
            allocate(ioi(intbuf))
            allocate(ioj(intbuf))
            allocate(buf(intbuf))

            veco = 0.0d0
            vecn = 0.0d0
            dumarr = 0.0d0
            ioi = 0
            ioj = 0
            buf = 0.0d0

            call interface_mpi_BCAST(veco,ladc,MASTER,
     &                               global_communicator)
            vecn = 0.0d0
            rewind(iobase)
            read(iobase) nbufsloc
            do i = 1,nbufsloc
              read(iobase,err=889) (buf(ixx),ixx=1,intbuf),
     &                             (ioi(ixx),ixx=1,intbuf),
     &                             (ioj(ixx),ixx=1,intbuf),
     &                             nact,jdummy
              do k = 1, nact
                irow = ioi(k)
                icol = ioj(k)
                b    = buf(k)
                vecn(irow) = vecn(irow) + veco(icol) * b
                if(irow.ne.icol) vecn(icol) = vecn(icol) + veco(irow)*b
              enddo
            enddo
            call interface_mpi_reduce(vecn,dumarr,ladc,
     &                 op_MPI_SUM,MASTER,global_communicator)


            deallocate(veco)
            deallocate(vecn)
            deallocate(dumarr)
            deallocate(ioi)
            deallocate(ioj)
            deallocate(buf)

          endif
          close(iobase)
          if(deb) write(*,*) 'Node',MYPROC,' closed local file ',locfn


        end select

!----------------------------------------------------------------
      ENDDO    ! ----------- END OF SERVER MAIN LOOP ------------
!----------------------------------------------------------------


      return
!
!  error exit for read/write routine.
!
 889  call quit('problem on node with reading/writing partial matrix')
      END SUBROUTINE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#endif
