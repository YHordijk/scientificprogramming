!      Copyright (c) 2019 by the authors of DIRAC.
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine LVECANALY_R(io1,io2,io3,desrep,ladc,nocc, &
                             irecl,filecnf,ionizlevel)
!
! purpose: reads in and analyzes long eigenvectors of 
! single/double ionization calculations 
! depending on the value of ionizlevel.
! Since eigenvectors can be *very* long we restrict the
! printout to eigenvector components larger than 0.001.
! Smaller ones are not considered to be relevant.
!
! Important: The old way to get rid of linearly dependent long eigenvectors 
! as done in RADCVECANALYZE_F90
! is questionable. Either we have true physical degeneracy or artificial degeneracy
! due to overiteration or low symmetry intensity spreading. In all these cases a
! Lanczos reiteration will lead to artificial long eigenvectors with high pole strength
! but these are spurious! We will therefore introduce a canonical orthogonalization
! procedure for these vectors in the corresponding subspace. Remember: the spaces of
! a possible linear dependence are in general very small so we do not produce numerical
! overhead.
!
! Since possible huge memory chunks are freed at this point (incore Lanczos) we can
! read in the long eigenvectors in memory (one is in general only interested in a very
! small subset of eigenvectors compared to the size of the ADC matrix).
!
      use adc_mat
      use adc_fano_exchange
      use memory_allocator

      implicit none
!
! transferred variables
!
      integer, intent(in) :: io1,io2,io3,desrep,ladc,nocc,irecl
      integer, intent(in) :: ionizlevel
      character*6, intent(in) :: filecnf
!
! local variables
!
      CHARACTER*10                :: LEVECFN,ANALOUT
      CHARACTER(len=irecl)        :: field
      CHARACTER*11                :: blankf='           '
      CHARACTER*11                :: ostr='    OCC    '
      CHARACTER*11                :: vstr='    VIR    '
      logical                     :: isthere
      integer                     :: neigenv,ix,ladcx,ispcount,nv,idum
      integer, allocatable        :: indarr(:), indarr2(:)
      integer, allocatable        :: irel(:)
      integer                     :: iallocstat,koff,loff,lwork,info
      integer                     :: i,j,k,l
      integer                     :: ispoff,nindep,icount
      real*8, allocatable         :: evec(:,:),eval(:),smat(:,:),work(:)
      real*8, allocatable         :: sevl(:),evli(:,:)
      real*8                      :: enorm,polestr,polethr,diffthr,seval
      real*8                      :: autoev  = 27.2113957d0
      real*8                      :: coefthr = 0.01d0
      real*8                      :: woutthr = 0.0001d0
!
! ------------------------ execution starts ----------------------
!
      call pst('Analyzer/Purger for long ADC real eigenvectors+')
      SELECT CASE (ionizlevel)
        CASE (1)
          write(*,*) 'Analyzing single ionizations.'   
        CASE (2)
          write(*,*) 'Analyzing double ionizations.'   
        CASE (3)
          write(*,*) '** Attention! Excitations now only accessible'
          write(*,*) '** in the new POLPRP module (see manual)'
          call quit('** Use POLPRP module')
        CASE DEFAULT
          call quit('** internal error ** ionizlevel selector')
      END SELECT
      write(*,*)
      write(*,*) 'Calling arguments:'
      write(*,*) 'Symmetry of final state:',DESREP
      write(*,*) 'Length of ADC matrix:   ',LADC
      write(*,*) 'Length of main space:   ',NOCC
      write(*,*) 'Record length of config file:',IRECL
      write(*,*) 'Generic config file  name: ',FILECNF
!     write(*,*) 'File handles: ',io1,io2,io3
!
! clear variable length character
!
      DO i=1,len(field)
        field(i:i)=' '
      ENDDO

! constructing eigenvector and analysis output file name

      IF(DESREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',DESREP
        SELECT CASE (ionizlevel)
          CASE(1)
            WRITE(ANALOUT,'(A8,I2)') 'SEVCANL.',DESREP
          CASE(2)
            WRITE(ANALOUT,'(A8,I2)') 'DEVCANL.',DESREP
          CASE(3)
            WRITE(ANALOUT,'(A8,I2)') 'XEVCANL.',DESREP
        END SELECT
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',DESREP
        SELECT CASE (ionizlevel)
          CASE(1)
            WRITE(ANALOUT,'(A9,I1)') 'SEVCANL.0',DESREP
          CASE(2)
            WRITE(ANALOUT,'(A9,I1)') 'DEVCANL.0',DESREP
          CASE(3)
            WRITE(ANALOUT,'(A9,I1)') 'XEVCANL.0',DESREP
        END SELECT
      ENDIF
      write(*,*) 'Fetching eigenvectors from: ',LEVECFN
      write(*,*) 'Writing analysis results to: ',ANALOUT

! check the existence of eigenvector and configuration file

      INQUIRE(file=levecfn,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Eigenvector file with ADC vectors not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif
      INQUIRE(file=filecnf,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Configuration file not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif

! opening eigenvector, configuration and output file

      OPEN(io1, FILE=LEVECFN, FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='old')
      OPEN(io2,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL, &
          STATUS='UNKNOWN')
      OPEN(io3, FILE=ANALOUT, FORM='FORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')


! start analysis

      read(io1) ladcx,neigenv
      write(*,*) neigenv,' eigenvectors are stored in ',levecfn
      if(ladcx.ne.ladc) then
        call quit(' ** internal error ** Length mismatch!')
      endif


      allocate(evec(ladcx,neigenv),stat=iallocstat)
      if(iallocstat.ne.0) then
        write(*,*) 'Not enough dynamic memory to store'
        write(*,*) 'eigenvectors! We have to skip the analysis.'
        return
      endif
      allocate(eval(neigenv))
      allocate(indarr(neigenv))

      polethr = 0.001d0; indarr=1

! read in the (reiterated) eigenvectors from file

      DO j=1,neigenv
        read(io1) eval(j)
        read(io1) (evec(ix,j),ix=1,LADCX)
        enorm=dot_product(evec(:,j),evec(:,j))
!       evec=evec/sqrt(enorm)
        polestr=dot_product(evec(1:nocc,j),evec(1:nocc,j))
        write(*,'(A,I5,A,3F16.8)') 'eigenvector',j,': E/PS/NRM:', &
              eval(j)*autoev,polestr,enorm
      ENDDO

! determine subspace dimensions according to eigenvalue degeneracy.
! generate corresponding index array.

      diffthr = 1.0D-06
      ispcount=1
      IF(neigenv.gt.1) THEN
        DO J=1,neigenv-1
          IF(dabs(eval(j+1) - eval(j)).LT.diffthr) THEN
            indarr(ispcount) = indarr(ispcount) + 1
          ELSE
            ispcount = ispcount + 1           
          ENDIF
        ENDDO
        write(*,*) 'Number of (pseudo-)subspaces in this symmetry:',&
                   ispcount
        ispoff=0
        DO J=1,ispcount
          write(*,*) ' (pseudo-)space',J,'energy',&
                     eval(ispoff+indarr(J))*autoev, &
                      ' dim:',indarr(J)
          ispoff = ispoff + indarr(J)
        ENDDO
      ENDIF

! In case of a fanoadc run, allocate the arrays for the initial state
! and the energies and polestrengths

      IF(reladc_md_isfano) THEN
        nr_in_evecs = 0
        CALL alloc(in_evecs,ladcx,ispcount)
        CALL alloc(in_energy,2,ispcount)
      END IF

! now loop over all subspaces and form the corresponding set of
! linearly independent vectors

      ispoff = 0
      Do 888 J=1,ispcount
! grab energy eigenvalue in this subspace
        seval = eval(1+ispoff)
        write(*,*) 'Treating space',J,seval*autoev
        nv = indarr(J)
        If (nv.gt.1) Then
! create overlap matrix
          nv = indarr(J); lwork=10*nv
          allocate(smat(nv,nv),work(lwork),sevl(nv))
          smat=0.0d0; work=0.0d0; sevl=0.0d0
          Do k=1,nv
          Do l=1,nv
            koff=ispoff+k
            loff=ispoff+l
            smat(k,l)=dot_product(evec(:,koff),evec(:,loff))
          Enddo
          Enddo
! diagonalize overlap matrix
          CALL DSYEV('V','U',nv,smat,nv,sevl,work,lwork,info)
! smat now contains eigenvectors, eigenvalues are ordered
          write(*,*)' Eigenvalues/vectors of OVL matrix:'
          Do l=1,nv
            write(*,*) l,sevl(l)
            write(*,*) smat(:,l)
          Enddo
! determine indices of relevant eigenvectors
          allocate(irel(nv));irel=0
          nindep=0
          do l=1,nv
            if(dabs(sevl(l)).gt.0.1d0) then
              nindep=nindep+1
              irel(nindep)=l
              write(*,*) '****** irel index:',l
            endif
          enddo
          write(*,*) 'Found',nindep,' linearly independent eigenvecors.'
! now form the subset of linearly independent eigenvectors
          allocate(evli(ladc,nindep))
          evli=0.0d0
          do k=1,nindep
            do l=1,nv
              evli(:,k) = evli(:,k) + smat(l,irel(k))*evec(:,ispoff+l)
            enddo
          enddo
! display norm of new eigenvectors
          do k=1,nindep
            write(*,*) '  Renormalizing new eigenvector:',k
            enorm=dot_product(evli(:,k),evli(:,k))
            evli(:,k) = evli(:,k)/sqrt(enorm)
            polestr=dot_product(evli(1:nocc,k),evli(1:nocc,k))
            write(*,*) '      PS of new eigenvector(s):',polestr
          enddo
          ispoff=ispoff+nv
          deallocate(irel,sevl,work,smat)
        Else 
! we have subspace dimension 1 and transfer vector for unified writing.
          write(*,*) 'No action taken, space is one-dimensional.'
          ispoff=ispoff+1
          allocate(evli(ladc,1))
          evli(:,1)=evec(:,ispoff)
          enorm=dot_product(evli(:,1),evli(:,1))
          write(*,*) '   Norm of eigenvector:',enorm
          polestr=dot_product(evli(1:nocc,1),evli(1:nocc,1))
          write(*,*) '      PS of eigenvector:',polestr
          nindep=1
        Endif

        write(*,*) 'Dim. check (nv,nindep,ispoff):',nv,nindep,ispoff
!
! start writing of vectors
! branch according to ionizlevel
!
!______________________________________________________
!|
!|

        DO 177 K=1,NINDEP
!
! it can happen that even a set of linearly independent vectors still does not
! acquire polestrength. We do not write these vectors into the ouput file!
!
        polestr=dot_product(evli(1:nocc,k),evli(1:nocc,k))
        if(polestr.lt.woutthr) then
          write(*,*) 'We skip',seval*autoev, &
                     ' due to vanishing pole strength.'
          cycle
        endif
        write(*,*) 'writing eigenvector for',seval,seval*autoev
        write(io3,'(A,F16.10,A,F16.10)') 'Eigenvector ', &
              seval*autoev,' eV,  Pole strength: ',polestr

        IF(ionizlevel.eq.1) then
!
          write(io3,'(A8,2X,3A11,A20)') &
            'No.',ostr,ostr,vstr,' Coeff(real)'
          write(io3,'(A)') ' '
          icount = 1
          DO ix=1,ladcx
            read(io2,REC=ix) field
            if(dabs(evli(ix,K)).gt.coefthr) then 
              IF(ix.gt.nocc) then
                write(io3,'(I8,2X,A33,F20.12)') icount,field,evli(ix,K)
              ELSE
                write(io3,'(I8,2X,3A11,F20.12)') &
                     icount,field(12:22),blankf,blankf,evli(ix,K)
              ENDIF
              icount = icount + 1
            endif
          ENDDO

          IF(reladc_md_isfano) THEN
            nr_in_evecs = nr_in_evecs + 1
            in_evecs(1:ladcx,nr_in_evecs) = evli(1:ladcx,k)
            in_energy(1,nr_in_evecs)      = seval
            in_energy(2,nr_in_evecs)      = polestr
          END IF
!
        ELSE IF(ionizlevel.eq.2) then
!
          write(io3,'(A8,2X,4A11,A20)') &
            'No.',ostr,ostr,ostr,vstr,'Coeff(real)'
          write(io3,'(A)') ' '
          icount = 1
          DO ix=1,ladcx
            read(io2,REC=ix) field
            if(dabs(evli(ix,K)).gt.coefthr) then
              IF(ix.gt.nocc) then
                write(io3,'(I8,2X,A44,F20.12)') icount,field,evli(ix,K)
              ELSE
                write(io3,'(I8,2X,A22,2A11,F20.12)') &
                     icount,field,blankf,blankf,evli(ix,K)
              ENDIF
              icount = icount + 1
            endif
          ENDDO
!
        ELSE   ! obsolete....
!
          call quit('** internal error, no excitations available here.')
!         write(io3,'(A8,2X,4A14,A20)') &
!           'Cnf.',ostr,vstr,ostr,vstr,'Coeff(real)'
!         write(io3,'(A)') ' '
!         icount = 1
!         DO ix=1,ladcx
!           read(io2,REC=ix) field
!           if(dabs(evli(ix,K)).gt.coefthr) then
!             IF(ix.gt.nocc) then
!               write(io3,'(I8,2X,A56,F20.12)') ix,field,evli(ix,K)
!             ELSE
!               write(io3,'(I8,2X,A28,2A14,F20.12)') &
!                    ix,field,blankf,blankf,evli(ix,K)
!             ENDIF
!             icount = icount + 1
!           endif
!         ENDDO
        ENDIF
        write(io3,'(A)') ' '
        write(io3,'(A)') '--------------------------------------------'
        write(io3,'(A)') ' '
        write(*,'(A)') ' '
        write(*,'(A)') '--------------------------------------------'
        write(*,'(A)') ' '

 177    CONTINUE
!|
!|
!|_____________________________________________________

        deallocate(evli)

 888  Continue


! release memory

      deallocate(indarr)
      deallocate(eval)
      deallocate(evec)

      close(io1)
      close(io2)
      close(io3)
      
      end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine LVECANALY_C(io1,io2,io3,desrep,ladc,nocc, &
                             irecl,filecnf,ionizlevel)
!
!  complex counterpart of LVECANALY_R
!  for comments, see this routine.
!
      use adc_mat
      use adc_fano_exchange
      use memory_allocator
!
      implicit none
!
! transferred variables
!
      integer, intent(in) :: io1,io2,io3,desrep,ladc,nocc,irecl
      integer, intent(in) :: ionizlevel
      character*6, intent(in) :: filecnf
!
!  define function interfaces called by this subroutine
!  Since the caller is in Fortran90 we need these interfaces.
!  Omitting them inevitably leads to an error.
!
      interface

         function dot_productz(a,b)
           complex*16  dot_productz
           complex*16, dimension(:), intent(in) :: a,b
         end function dot_productz

         function vnorm_cmplx(a)
           real*8 vnorm_cmplx
           complex*16, dimension(:), intent(in) :: a
         end function vnorm_cmplx

      end interface
!
! local variables
!
      CHARACTER*10                :: LEVECFN,ANALOUT
      CHARACTER*44                :: field
      CHARACTER*11                :: blankf='           '
      CHARACTER*11                :: ostr='    OCC    '
      CHARACTER*11                :: vstr='    VIR    '
      logical                     :: isthere
      integer                     :: neigenv,j,ix,ladcx,ispcount,nv,idum
      integer, allocatable        :: indarr(:), indarr2(:)
      integer, allocatable        :: irel(:)
      integer                     :: iallocstat,k,l,koff,loff,lwork,info
      integer                     :: ispoff,nindep,icount
      complex*16, allocatable     :: evec(:,:)    !eigenvectors are now complex!
      complex*16, allocatable     :: evli(:,:)    !aux array for final complex eigenvectors
      complex*16, allocatable     :: smat(:,:),work(:)   !complex ovl, complex work array
      real*8, allocatable         :: rwork(:)      ! real work array for ZHEEV
      real*8, allocatable         :: eval(:)      !eigenvalues of read-in vectors
      real*8, allocatable         :: sevl(:)   !  *real* eigenvalues of OVL matrix
      real*8                      :: enorm,polestr,polethr,diffthr,seval
      real*8                      :: autoev  = 27.2113957d0
      real*8                      :: coefthr = 0.01d0
      real*8                      :: woutthr = 0.0001d0
      complex*16                  :: A0
      real*8                      :: R0
!
! ------------------------ execution starts ----------------------
!
      call pst('Analyzer/Purger for long ADC complex eigenvectors+')
      SELECT CASE (ionizlevel)
        CASE (1)
          write(*,*) 'Analyzing single ionizations.'   
        CASE (2)
          write(*,*) 'Analyzing double ionizations.'   
        CASE DEFAULT
          call quit('** internal error ** ionizlevel selector')
      END SELECT
      write(*,*)
      write(*,*) 'Calling arguments:'
      write(*,*) 'Symmetry of final state:',DESREP
      write(*,*) 'Length of ADC matrix:   ',LADC
      write(*,*) 'Length of main space:   ',NOCC
      write(*,*) 'Record length of config file:',IRECL
      write(*,*) 'Generic config file  name: ',FILECNF
!     write(*,*) 'File handles: ',io1,io2,io3

! constructing eigenvector and analysis output file name

      IF(DESREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',DESREP
        IF(ionizlevel.eq.1) THEN
          WRITE(ANALOUT,'(A8,I2)') 'SEVCANL.',DESREP
        ELSE
          WRITE(ANALOUT,'(A8,I2)') 'DEVCANL.',DESREP
        ENDIF
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',DESREP
        IF(ionizlevel.eq.1) THEN
          WRITE(ANALOUT,'(A9,I1)') 'SEVCANL.0',DESREP
        ELSE
          WRITE(ANALOUT,'(A9,I1)') 'DEVCANL.0',DESREP
        ENDIF
      ENDIF
      write(*,*) 'Fetching eigenvectors from: ',LEVECFN
      write(*,*) 'Writing analysis results to: ',ANALOUT

! check the existence of eigenvector and configuration file

      INQUIRE(file=levecfn,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Eigenvector file with ADC vectors not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif
      INQUIRE(file=filecnf,exist=isthere)
      if(.not.isthere) then
        write(*,*) '******* ATTENTION ********'
        write(*,*) 'Configuration file not present'
        write(*,*) 'for symmetry',desrep
        write(*,*) 'Skipping analysis'
        return
      endif

! opening eigenvector, configuration and output file

      OPEN(io1, FILE=LEVECFN, FORM='UNFORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='old')
      OPEN(io2,FILE=FILECNF,ACCESS='DIRECT',RECL=IRECL, &
          STATUS='UNKNOWN')
      OPEN(io3, FILE=ANALOUT, FORM='FORMATTED', &
          ACCESS='SEQUENTIAL', STATUS='UNKNOWN')


! start analysis

      A0=(0.0D0,0.0D0); R0=0.0d0
      read(io1) ladcx,neigenv
      write(*,*) neigenv,' eigenvectors are stored in ',levecfn
      if(ladcx.ne.ladc) then
        call quit(' ** internal error ** Length mismatch!')
      endif


      allocate(evec(ladcx,neigenv),stat=iallocstat)
      if(iallocstat.ne.0) then
        write(*,*) 'Not enough dynamic memory to store'
        write(*,*) 'eigenvectors! We have to skip the analysis.'
        return
      endif
      allocate(eval(neigenv))
      allocate(indarr(neigenv))

      polethr = 0.001d0; indarr=1

! read in the (reiterated) eigenvectors from file

      DO j=1,neigenv
        read(io1) eval(j)
        read(io1) (evec(ix,j),ix=1,LADCX)
!       write(*,*) 'waaahhhh:', &
!         size(evec(:,j)),size(evec)
        enorm=vnorm_cmplx(evec(:,j))
!       evec=evec/sqrt(enorm)
        polestr=vnorm_cmplx(evec(1:nocc,j))
        write(*,'(A,I5,A,3F16.8)') 'eigenvector',j,': E/PS/NRM:', &
              eval(j)*autoev,polestr,enorm
      ENDDO

! determine subspace dimensions according to eigenvalue degeneracy.
! generate corresponding index array.

      diffthr = 1.0D-06
      ispcount=1
      IF(neigenv.gt.1) THEN
        DO J=1,neigenv-1
          IF(dabs(eval(j+1) - eval(j)).LT.diffthr) THEN
            indarr(ispcount) = indarr(ispcount) + 1
          ELSE
            ispcount = ispcount + 1           
          ENDIF
        ENDDO
        write(*,*) 'Number of (pseudo-)subspaces in this symmetry:',&
                   ispcount
        ispoff=0
        DO J=1,ispcount
          write(*,*) ' (pseudo-)space',J,'energy',&
                     eval(ispoff+indarr(J))*autoev, &
                      ' dim:',indarr(J)
          ispoff = ispoff + indarr(J)
        ENDDO
      ENDIF

! In case of a fanoadc run, allocate the arrays for the initial state
! and the energies and polestrengths

      IF(reladc_md_isfano) THEN
        nr_in_evecs = 0
        CALL alloc(cin_evecs,ladcx,ispcount)
        CALL alloc(in_energy,2,ispcount)
      END IF

! now loop over all subspaces and form the corresponding set of
! linearly independent vectors

      ispoff = 0
      Do 888 J=1,ispcount
! grab energy eigenvalue in this subspace
        seval = eval(1+ispoff)
        write(*,*) 'Treating space',J,seval*autoev
        nv = indarr(J)
        If (nv.gt.1) Then
! create overlap matrix. Can have complex entries but is Hermitian!
          nv = indarr(J); lwork=10*nv
          allocate(smat(nv,nv),work(lwork),rwork(lwork),sevl(nv))
          smat=A0; work=A0; rwork=R0; sevl=R0
          Do k=1,nv
          Do l=1,nv
            koff=ispoff+k
            loff=ispoff+l
            smat(k,l)=dot_productz(evec(:,koff),evec(:,loff))
          Enddo
          Enddo
! diagonalize hermitian overlap matrix
          CALL ZHEEV('V','U',nv,smat,nv,sevl,work,lwork,rwork,info)
! smat now contains eigenvectors, eigenvalues are ordered
          write(*,*)' Eigenvalues/vectors of OVL matrix:'
          Do l=1,nv
            write(*,*) l,sevl(l)
            Do k=1,nv
              write(*,'(2F20.10)') smat(k,l)
            Enddo
            write(*,*)
          Enddo
! determine indices of relevant eigenvectors
          allocate(irel(nv));irel=0
          nindep=0
          do l=1,nv
            if(dabs(sevl(l)).gt.0.1d0) then
              nindep=nindep+1
              irel(nindep)=l
              write(*,*) '****** irel index:',l
            endif
          enddo
          write(*,*) 'Found',nindep,' linearly independent eigenvecors.'
! now form the subset of linearly independent eigenvectors
          allocate(evli(ladc,nindep))
          evli=A0
          do k=1,nindep
            do l=1,nv
              evli(:,k) = evli(:,k) + smat(l,irel(k))*evec(:,ispoff+l)
            enddo
          enddo
! display norm of new eigenvectors
          do k=1,nindep
            write(*,*) '  Renormalizing new eigenvector:',k
            enorm=vnorm_cmplx(evli(:,k))
            evli(:,k) = evli(:,k)/cmplx(sqrt(enorm),0.0d0)
            polestr=vnorm_cmplx(evli(1:nocc,k))
            write(*,*) '      PS of new eigenvectors:',polestr
          enddo
          ispoff=ispoff+nv
          deallocate(irel,sevl,work,rwork,smat)
        Else 
! we have subspace dimension 1 and transfer vector for unified writing.
          write(*,*) 'No action taken, space is one-dimensional.'
          ispoff=ispoff+1
          allocate(evli(ladc,1))
          evli(:,1)=evec(:,ispoff)
          enorm=vnorm_cmplx(evli(:,1))
          write(*,*) '   Norm of eigenvector:',enorm
          polestr=vnorm_cmplx(evli(1:nocc,1))
          write(*,*) '      PS of eigenvector:',polestr
          nindep=1
        Endif

        write(*,*) 'Dim. check (nv,nindep,ispoff):',nv,nindep,ispoff
!
! start writing of vectors
! branch according to ionizlevel
!
!______________________________________________________
!|
!|
        !IF(reladc_md_isfano) THEN                                      
        !  nr_in_evecs = nindep                                         
        !  CALL alloc(cin_evecs,ladcx,nr_in_evecs)                       
        !END IF    

        DO 177 K=1,NINDEP

        polestr=vnorm_cmplx(evli(1:nocc,k))
        if(polestr.lt.woutthr) then
          write(*,*) 'We skip',seval*autoev, &
                     ' due to vanishing pole strength.'
          cycle
        endif
        write(*,*) 'writing eigenvector for',seval,seval*autoev
        write(io3,'(A,F16.10,A,F16.10)') 'Eigenvector ', &
              seval*autoev,' eV,  Pole strength: ',polestr

        IF(ionizlevel.eq.1) then
          write(io3,'(A8,2X,3A11,2A20)') &
            'No.',ostr,ostr,vstr,' Coeff(real)',' Coeff(imag)'
          write(io3,'(A)') ' '
          icount = 1
          DO ix=1,ladcx
            read(io2,REC=ix) field
            if(abs(evli(ix,K)).gt.coefthr) then 
              IF(ix.gt.nocc) then
                write(io3,'(I8,2X,A33,2F20.12)') icount,field,evli(ix,K)
              ELSE
                write(io3,'(I8,2X,3A11,2F20.12)') &
                     icount,field(12:22),blankf,blankf,evli(ix,K)
              ENDIF
              icount = icount + 1
            endif
          ENDDO

          IF(reladc_md_isfano) THEN
            nr_in_evecs = nr_in_evecs + 1
            cin_evecs(1:ladcx,k) = evli(1:ladcx,k)
            in_energy(1,nr_in_evecs)      = seval
            in_energy(2,nr_in_evecs)      = polestr
          END IF

        ELSE
          write(io3,'(A8,2X,4A11,2A20)') &
            'No.',ostr,ostr,ostr,vstr,' Coeff(real)',' Coeff(imag)'
          write(io3,'(A)') ' '
          icount = 1
          DO ix=1,ladcx
            read(io2,REC=ix) field
            if(abs(evli(ix,K)).gt.coefthr) then
              IF(ix.gt.nocc) then
                write(io3,'(I8,2X,A44,2F20.12)') icount,field,evli(ix,K)
              ELSE
                write(io3,'(I8,2X,A22,2A11,2F20.12)') &
                     icount,field,blankf,blankf,evli(ix,K)
              ENDIF
              icount = icount + 1
            endif
          ENDDO
        ENDIF
        write(io3,'(A)') ' '
        write(io3,'(A)') '--------------------------------------------'
        write(io3,'(A)') ' '
        write(*,'(A)') ' '
        write(*,'(A)') '--------------------------------------------'
        write(*,'(A)') ' '

 177    CONTINUE
!|
!|
!|_____________________________________________________

        deallocate(evli)

 888  Continue


! release memory

      deallocate(indarr)
      deallocate(eval)
      deallocate(evec)

      close(io1)
      close(io2)
      close(io3)
      
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine INITINCORE_R(ioin,ioout,nbufs,nincw,intbuf,  &
                   buf,ioi,ioj,bufc,ioic,iojc,ivalid)
!
!  initialize incore Lanczos. Read in as much of the matrix as possible
!
      implicit none
      integer, intent(in)      :: ioin,ioout  ! file handles for accessing ADC matrix (is open!)
      integer, intent(inout)   :: nbufs       ! number of buffers the ADC matrix has at the moment!
      integer, intent(in)      :: nincw       ! number of words available for incore storage
      integer, intent(in)      :: intbuf      ! buffer size
      real*8                   :: buf(intbuf)
      integer                  :: ioi(intbuf),ioj(intbuf)
      real*8                   :: bufc(nincw)
      integer                  :: ioic(nincw),iojc(nincw)
      integer, intent(inout)   :: ivalid      !ivalid contains the number of valid elements in the huge array.
!
!  local variables
!
      integer      :: nbufinc,nbufooc   !number of buffers in/out of core
      integer      :: i,j,k,irec,itest,ixx,nact,jdummy
      character*6  :: newfile='ADCXPK'
!
!  execution starts
!
      call pst('Initializing Incore(R) Lanczos. Reading matrix...+')
      write(*,*) 'ADC matrix has',nbufs,' buffers of length',intbuf
      write(*,*) 'Available words for incore storage:',nincw
      nbufinc = nincw/intbuf
      if(nbufinc.gt.nbufs) nbufinc = nbufs
      write(*,*) 'Number of buffers held in core:',nbufinc
      nbufooc = nbufs - nbufinc
      write(*,*) 'Number of buffers remaining external:',nbufooc
!
!  read nbufinc buffer
!
      ivalid=0
      REWIND(ioin)
      DO IREC = 1,nbufinc
        READ(ioin,iostat=itest) (BUF(IXX),IXX=1,INTBUF),  &
                                (IOI(IXX),IXX=1,INTBUF),  &
                                (IOJ(IXX),IXX=1,INTBUF),  &
                                NACT,JDUMMY
        if(itest.ne.0)  &
          call quit('Matrix read error (1) in INITINCORE')
        DO K = 1, nact
          ivalid = ivalid + 1
          ioic(ivalid) = IOI(K)
          iojc(ivalid) = IOJ(K)
          bufc(ivalid) = BUF(K)
        ENDDO
      ENDDO
      write(*,*) nbufinc,' records read.'
      write(*,*) 'Number of valid elements:',ivalid
!
!  write the remaining records in a new packed ADC matrix file
!  and transfer control to this new file with reduced number of buffers
!  If ADC matrix is completely incore the file disappears from disk
!  and is never accessed for the rest of the program anymore.
!
      OPEN(unit=ioout,file=newfile,form='UNFORMATTED',  &
           status='UNKNOWN')
      DO IREC = nbufinc+1,nbufs
        READ(ioin,iostat=itest) (BUF(IXX),IXX=1,INTBUF),  &
                                (IOI(IXX),IXX=1,INTBUF),  &
                                (IOJ(IXX),IXX=1,INTBUF),  &
                                NACT,JDUMMY
        if(itest.ne.0)  &
          call quit('Matrix read error (2) in INITINCORE')
        WRITE(ioout,iostat=itest) (BUF(IXX),IXX=1,INTBUF), &
                                  (IOI(IXX),IXX=1,INTBUF), &
                                  (IOJ(IXX),IXX=1,INTBUF), &
                                   NACT,JDUMMY
        if(itest.ne.0)  &
          call quit('Matrix write error in INITINCORE')
      ENDDO

!
!  remaining ADC file is written. The old one is deleted and
!  buffer sizes are adjusted
!
      CLOSE(IOOUT)
      CLOSE(IOIN,status='delete')
      OPEN(unit=IOIN,file=newfile,form='UNFORMATTED',  &
           status='UNKNOWN')
      NBUFS = nbufooc
      if(nbufs.eq.0) then
        write(*,*) 'ADC matrix is completely in core.'
      else
        write(*,*) 'Adjusted number of buffers to',nbufs
      endif

      return
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine INITINCORE_C(ioin,ioout,nbufs,nincw,intbuf,  &
                   buf,ioi,ioj,bufc,ioic,iojc,ivalid)
!
!  initialize incore Lanczos. Read in as much of the matrix as possible
!
      implicit none
      integer, intent(in)      :: ioin,ioout  ! file handles for accessing ADC matrix (is open!)
      integer, intent(inout)   :: nbufs       ! number of buffers the ADC matrix has at the moment!
      integer, intent(in)      :: nincw       ! number of words available for incore storage
      integer, intent(in)      :: intbuf      ! buffer size
      complex*16               :: buf(intbuf)
      integer                  :: ioi(intbuf),ioj(intbuf)
      complex*16               :: bufc(nincw)
      integer                  :: ioic(nincw),iojc(nincw)
      integer, intent(inout)   :: ivalid      !ivalid contains the number of valid elements in the huge array.
!
!  local variables
!
      integer      :: nbufinc,nbufooc   !number of buffers in/out of core
      integer      :: i,j,k,irec,itest,ixx,nact,jdummy
      character*6  :: newfile='ADCXPK'
!
!  execution starts
!
      call pst('Initializing Incore(C) Lanczos. Reading matrix...+')
      write(*,*) 'ADC matrix has',nbufs,' buffers of length',intbuf
      write(*,*) 'Available words for incore storage:',nincw
      nbufinc = nincw/intbuf
      if(nbufinc.gt.nbufs) nbufinc = nbufs
      write(*,*) 'Number of buffers held in core:',nbufinc
      nbufooc = nbufs - nbufinc
      write(*,*) 'Number of buffers remaining external:',nbufooc
!
!  read nbufinc buffer
!
      ivalid=0
      REWIND(ioin)
      DO IREC = 1,nbufinc
        READ(ioin,iostat=itest) (BUF(IXX),IXX=1,INTBUF),  &
                                (IOI(IXX),IXX=1,INTBUF),  &
                                (IOJ(IXX),IXX=1,INTBUF),  &
                                NACT,JDUMMY
        if(itest.ne.0)  &
          call quit('Matrix read error (1) in INITINCORE')
        DO K = 1, nact
          ivalid = ivalid + 1
          ioic(ivalid) = IOI(K)
          iojc(ivalid) = IOJ(K)
          bufc(ivalid) = BUF(K)
        ENDDO
      ENDDO
      write(*,*) nbufinc,' records read.'
      write(*,*) 'Number of valid elements:',ivalid
!
!  write the remaining records in a new packed ADC matrix file
!  and transfer control to this new file with reduced number of buffers
!  If ADC matrix is completely incore the file disappears from disk
!  and is never accessed for the rest of the program anymore.
!
      OPEN(unit=ioout,file=newfile,form='UNFORMATTED',  &
           status='UNKNOWN')
      DO IREC = nbufinc+1,nbufs
        READ(ioin,iostat=itest) (BUF(IXX),IXX=1,INTBUF),  &
                                (IOI(IXX),IXX=1,INTBUF),  &
                                (IOJ(IXX),IXX=1,INTBUF),  &
                                NACT,JDUMMY
        if(itest.ne.0)  &
          call quit('Matrix read error (2) in INITINCORE')
        WRITE(ioout,iostat=itest) (BUF(IXX),IXX=1,INTBUF), &
                                  (IOI(IXX),IXX=1,INTBUF), &
                                  (IOJ(IXX),IXX=1,INTBUF), &
                                   NACT,JDUMMY
        if(itest.ne.0)  &
          call quit('Matrix write error in INITINCORE')
      ENDDO

!
!  remaining ADC file is written. The old one is deleted and
!  buffer sizes are adjusted
!
      CLOSE(IOOUT)
      CLOSE(IOIN,status='delete')
      OPEN(unit=IOIN,file=newfile,form='UNFORMATTED',  &
           status='UNKNOWN')
      NBUFS = nbufooc
      if(nbufs.eq.0) then
        write(*,*) 'ADC matrix is completely in core.'
      else
        write(*,*) 'Adjusted number of buffers to',nbufs
      endif

      return
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine TWOHPOP_R(lupri,io_lan,ibeta,nrep,multb,no,noot)
!
!  Purpose: Prepare the propagator contributions to the reduced density matrix.
!
!  the Delta_mn matrix is then written to a file 'DADCPOP'
!  in the format
!              integer dimension
!              real(:,:) dmn   :  real part, vector 1
!              real(:,:) dmn   :  imaginary part, vector 1
!                       ....
!              real(:,:) dmn   :  real part, vector N
!              real(:,:) dmn   :  imaginary part, vector N
!
      implicit none
      integer, intent(in)  :: lupri    ! file handle for printing to stdout
      integer, intent(in)  :: io_lan   ! file handle start
      integer, intent(in)  :: ibeta    ! symmetry of the two-hole final state 
      integer, intent(in)  :: nrep     ! number of ireps in the overall run
      integer, intent(in)  :: multb(64,64,2)     ! fixed dimension multiplication table
      integer, intent(in)  :: no(nrep)    ! number of occupied orbitals in each irep
      integer, intent(in)  :: noot(nrep)  ! number of OO-combinations in each irep
!
!  local variables
!
      real*8,  allocatable  :: dmn(:,:), xmna(:,:), xmn2(:,:),xmnc(:,:)
      real*8,  allocatable  :: xmnd(:,:),thvecs(:,:),evec(:)
      real*8,  allocatable  :: cevl(:), polestr(:)
      real*8,  allocatable  :: dmnext(:,:,:)
      real*8                :: polethr = 0.001d0
      real*8                :: autoev = 27.2113957D0
      real*8                :: A0=0.0d0, A1=1.0d0
      real*8                :: evltri,polestrtri
      integer, allocatable  :: no_offs(:)
      integer               :: i,j,ix,ivcount,iwwcount
      integer               :: n1,n2
      integer               :: nocc,ipoint,no_sum
      integer               :: irp,jrp,imin,irow,jcol
      logical               :: isthere
      character*7           :: dmn_filename='DADCPOP'

!  say hello and check basic information

      call PST('Performing 2H population analysis (real vectors)+')
      nocc=noot(ibeta)
      write(lupri,*) 'Final state symmetry:',ibeta
      write(lupri,*) 'Number of 2h states in this symm:',nocc

!__________________________________________________________________
!|
!|  read the file of Lanczos eigenvectors and extract the 2h/2h block
!|        ----->   in this symmetry !!  <-----
!|  from it. 
!|  Fill the corresponding arrays with the read vectors and determine
!|  the polestrengths. The raw Lanc eigenvector goes to *evec* and
!|  the 2h/2h block goes to *thvecs*

      Inquire(file='TMATEVC',exist=isthere)
      IF(.not.isthere) THEN
        write(lupri,*) '**************'
        write(lupri,*) '**** Attention'
        write(lupri,*) '**************'
        write(lupri,*) 'TMATEVC file not present.'
        write(lupri,*) 'No population analysis possible'
        return
      ENDIF

      OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(IO_LAN)
      READ(IO_LAN) N1
      READ(IO_LAN) N2
      WRITE(lupri,*) 'Found',N1,' Lanczos eigenvectors of length', &
                   N2,' on TMATEVC'
!     WRITE(lupri,*) nocc,' of them are 2H states!'

      allocate(thvecs(nocc,N1), &
               evec(n2),          &
               cevl(N1),        &
               polestr(N1))

      thvecs = 0.0d0; evec = 0.0d0; cevl = 0.0d0; polestr = 0.0d0
!
! extract relevant vectors from TMATEVC, as long as spurious vectors are
! still sitting there (will be changed in the future)
!
! only for the relevant vectors we keep the eigenvalue, the eigenvector and
! the pole strength in the corresponding arrays.
!
      iwwcount = 0
      DO I=1,N1   ! now check all vectors sitting in TMATEVC for suitable 2h analysis
        READ(IO_LAN) evltri
        READ(IO_LAN) (evec(ix),ix=1,n2)
        polestrtri=dot_product(evec(1:nocc),evec(1:nocc))
        if(polestrtri.gt.polethr) then
              write(lupri,'(A,I5,2(A,F16.6),A)')    &
                  'Taking vector ',i,             &
                  ' of energy',evltri*autoev,     &
                  ' and polestr.',polestrtri,     &
                  ' for population analysis'
              iwwcount = iwwcount + 1
              thvecs(:,iwwcount) = evec(1:nocc)    ! copy corresponding hole components
              cevl(iwwcount)     = evltri          ! copy corresponding eigenvalue (in a.u.)
              polestr(iwwcount)  = polestrtri      ! copy corresponding eigenvalue (in a.u.)
        endif 
      ENDDO
      write(lupri,*) iwwcount, &
        ' vectors found from',N1,' vectors suitable for 2H analysis'
      CLOSE(IO_LAN,STATUS='DELETE')
!|
!|
!|
!|
!|_______ 2h/2h block in memory together with the polestrengths ______


!  prepare occupied orbital offset table

      allocate(no_offs(nrep+1))
      no_offs = 0; no_sum = 0
      DO i=1,nrep
        no_offs(i+1) = no_offs(i) + no(i)
        no_sum = no_sum + no(i)
      ENDDO
!     write(lupri,*) 'offset table:',no_offs
!     write(lupri,*) 'total number of occupied orbitals:',no_sum

!  allocate space for Delta_mn matrix and clear matrix

      allocate(dmn(no_sum,no_sum))
      allocate(dmnext(no_sum,no_sum,2))
      allocate(xmna(no_sum,no_sum))
      allocate(xmnc(no_sum,no_sum))
      allocate(xmnd(no_sum,no_sum))
      allocate(xmn2(no_sum,no_sum))

!  open DADCPOP file for sequential access and write 
!  number of matrices and corresponding matrix dimension
!
!  ATT We check if file exists already and append the next symmetry
!  to it. This liberates us from writing symmetry-specific files altogether.
 
      inquire(file=dmn_filename,exist=isthere)

      select case(isthere)

      case (.FALSE.)      ! Delta_mn file was not created
      
        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='new',         &
             access='sequential',  &
             form='unformatted')

      case (.TRUE.)       ! Delta_mn file is there.

        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='old',         &
             access='sequential',  &
             form='unformatted',   &
             position='append')

      end select

!  write relevant data in this symmetry to the DNM file
!  symmetry, # of analyzed vectors, sum of occupied orbitals

      write(IO_LAN) ibeta, iwwcount, no_sum

!__________________________________________________________________
!|
!|  Loop over all vectors in the 2h/2h block possessing high enough 
!|  pole strength

      DO ivcount=1,iwwcount 

        xmna = 0.0d0;  xmn2 = 0.0d0; xmnc = 0.0d0; xmnd = 0.0d0

        write(lupri,*) '   ** doing analysis for vector:',ivcount
        write(lupri,*) '   ** eigenvalue: ',cevl(ivcount)*autoev

!
!  construct x_mn matrix
!
        ipoint = 0
        DO 14 JRP = 1, NREP
          IRP = MULTB(JRP,ibeta+NREP,2)
          IF (IRP.LT.JRP) GOTO 14
          DO J = 1, NO(JRP)
            IMIN = 1
            IF (IRP.EQ.JRP) IMIN = J + 1
            DO I = IMIN, NO(IRP)
              ipoint = ipoint + 1   ! get next eigenvector component
              irow   = no_offs(irp) + I
              jcol   = no_offs(jrp) + J
              xmna(irow,jcol) = thvecs(ipoint,ivcount)
!             write(*,'(A,F10.6,2I5)') '   read in components:', &
!       thvecs(ipoint,ivcount),irow,jcol
!             write(*,'(A,4I4,A,2I4,F10.6)') & 
!                   ' irp,jrp,i,j',irp,jrp,i,j, &
!                   ' row/col/entry:',irow,jcol,xmna(irow,jcol)
            ENDDO
          ENDDO 
  14    CONTINUE

        IF( (ipoint).NE.NOOT(ibeta)) THEN
          WRITE(lupri,*) 'ipoint,noot',ipoint,noot(ibeta)
          CALL QUIT('Error in TWOHPOP for 2h configurations!')
        ENDIF

!       write(lupri,*) 'Contributions from xmna:'
!       do i=1,no_sum
!         write(lupri,'(40F10.5)') (xmna(i,j),j=1,no_sum)
!       enddo

!  contract x_ml and x*_nl. XGEMM can be used because all irrelevant entries are
!  filled with zeros and will not contribute to the matrix product.
!  Second call: contract x_lm and x_ln. Transpose the first matrix.

! first part: contract x_ml . x_ln* over l<m,n
! result goes into xmnc
! 
        xmn2 = transpose(xmna)  ! form x_ln from x_nl. 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmna,no_sum, &
                    xmn2,no_sum,A0,xmnc,no_sum)

!       write(lupri,*) 'Contributions from xmnc (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmnc(i,j)),j=1,no_sum)
!       enddo
!
! second part: form x_lm from x_ml (in xmna) and then
! form x_ln* result goes into xmnd
!
        xmn2 = transpose(xmna) 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmn2,no_sum, &
                    xmna,no_sum,A0,xmnd,no_sum)
        xmnd = transpose(xmnd)
        xmn2 = xmnc + xmnd

!       write(lupri,*) 'Contributions from part 1+2 (real only):'
!       do i=1,no_sum
!         write(lupri,'(40F10.5)') (xmn2(i,j),j=1,no_sum)
!       enddo
!       stop

!  now form Delta_mn matrix. Part 1 is just one on the diagonal,
!  part 2 is just minus the pole strength on the dagonal.
        
        dmn = 0.0d0
        DO i=1,no_sum
          dmn(i,i) = 1.0d0 - polestr(ivcount)
!         write(*,*) 'real polestrength:',polestr(ivcount)
        ENDDO
        dmn = dmn + xmn2

        write(lupri,*) 'Delta_mn:'
        do i=1,no_sum
          write(lupri,'(40F10.5)') (dmn(i,j),j=1,no_sum)
        enddo

!  write delta_mn matrix to file in the format requested by dirac conversion
!  routine.

        dmnext=0.0d0
        dmnext(:,:,1)=dmn(:,:)
        dmn=0.0d0
        dmnext(:,:,2)=dmn(:,:)
        write(io_lan) dmnext

      ENDDO   !ivcount

!  perform all deallocations in reverse order of the allocations
!  and close file

!     write(lupri,*) 'entering deallocation section'

      deallocate(xmn2)
      deallocate(xmnd)
      deallocate(xmnc)
      deallocate(xmna)
      deallocate(dmnext)
      deallocate(dmn)
      deallocate(no_offs)
      deallocate(polestr)
      deallocate(cevl)
      deallocate(evec)
      deallocate(thvecs)

!     write(lupri,*) 'deallocation finished'

      close(io_lan)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine TWOHPOP_C(lupri,io_lan,ibeta,nrep,multb,no,noot)
!
!  Purpose: Prepare the propagator contributions to the reduced density matrix.
!
!  the Delta_mn matrix is then written to a file 'DADCPOP'
!  in the format
!              integer dimension
!              real(:,:) dmn   :  real part, vector 1
!              real(:,:) dmn   :  imaginary part, vector 1
!                       ....
!              real(:,:) dmn   :  real part, vector N
!              real(:,:) dmn   :  imaginary part, vector N
!
      implicit none
      integer, intent(in)  :: lupri    ! file handle for printing to stdout
      integer, intent(in)  :: io_lan   ! file handle start
      integer, intent(in)  :: ibeta    ! symmetry of the two-hole final state 
      integer, intent(in)  :: nrep     ! number of ireps in the overall run
      integer, intent(in)  :: multb(64,64,2)     ! fixed dimension multiplication table
      integer, intent(in)  :: no(nrep)    ! number of occupied orbitals in each irep
      integer, intent(in)  :: noot(nrep)  ! number of OO-combinations in each irep
!
!  function interface
!
      interface

         function vnorm_cmplx(a)
           real*8 vnorm_cmplx
           complex*16, dimension(:), intent(in) :: a
         end function vnorm_cmplx

      end interface
!
!  local variables
!
      complex*16, allocatable   :: dmn(:,:),        &
                                   xmna(:,:),       &
                                   xmn2(:,:),       &
                                   xmnc(:,:),       &
                                   xmnd(:,:)
      complex*16, allocatable   :: thvecs(:,:), evec(:)
      real*8, allocatable       :: cevl(:), polestr(:), dmnext(:,:,:)
      real*8                    :: polethr = 0.001d0
      real*8                    :: autoev = 27.2113957D0
      complex*16                :: A0=(0.0d0,0.0d0), A1=(1.0d0,0.0d0)
      integer, allocatable      :: no_offs(:)
      integer                   :: i,j,ix,ivcount,iwwcount
      integer                   :: n1,n2
      integer                   :: nocc,ipoint,no_sum
      integer                   :: irp,jrp,imin,irow,jcol
      logical                   :: isthere
      character*7               :: dmn_filename='DADCPOP'

!  say hello and check basic information

      call PST('Performing 2H population analysis (complex vectors)+')
      nocc=noot(ibeta)
      write(lupri,*) 'Final state symmetry:',ibeta
      write(lupri,*) 'Number of 2h states in this symm:',nocc

!__________________________________________________________________
!|
!|  read the file of Lanczos eigenvectors and extract the 2h/2h block
!|        ----->   in this symmetry !!  <-----
!|  from it. 
!|  Fill the corresponding arrays with the read vectors and determine
!|  the polestrengths. The raw Lanc eigenvector goes to *evec* and
!|  the 2h/2h block goes to *thvecs*

      Inquire(file='TMATEVC',exist=isthere)
      IF(.not.isthere) THEN
        write(lupri,*) '**************'
        write(lupri,*) '**** Attention'
        write(lupri,*) '**************'
        write(lupri,*) 'TMATEVC file not present.'
        write(lupri,*) 'No population analysis possible'
        return
      ENDIF

      OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='OLD')
      REWIND(IO_LAN)
      READ(IO_LAN) N1
      READ(IO_LAN) N2
      WRITE(lupri,*) 'Found',N1,' Lanczos eigenvectors of length', &
                   N2,' on TMATEVC'
!     WRITE(lupri,*) nocc,' of them are 2H states!'

      allocate(thvecs(nocc,nocc), &
               evec(n2),          &
               cevl(nocc),        &
               polestr(nocc))

      thvecs = A0; evec = A0; cevl = 0.0d0; polestr = 0.0d0

      iwwcount = 0
      DO I=1,nocc
        READ(IO_LAN) cevl(i)
        READ(IO_LAN) (evec(ix),ix=1,n2)
!       write(*,*) 'cevl,evec:',cevl,'\\',evec
!       pause
        polestr(i)=vnorm_cmplx(evec)
        if(polestr(i).gt.polethr) then
              write(lupri,'(A,I5,A,2F16.6,A)')           &
                  'Taking vector ',i,                    &
                  ' with', cevl(i)*autoev,polestr(i),    &
                  ' for population analysis'
              iwwcount = iwwcount + 1
              thvecs(:,i) = evec(1:nocc)   ! copy corresponding components
!             write(*,*) 'corr ev:',evec(1:nocc)
        endif 
      ENDDO
      write(lupri,*) iwwcount, &
        ' vectors have enough pole strength for 2H analysis'
      CLOSE(IO_LAN,STATUS='DELETE')
!|
!|
!|
!|
!|_______ 2h/2h block in memory together with the polestrengths ______


!  prepare occupied orbital offset table

      allocate(no_offs(nrep+1))
      no_offs = 0; no_sum = 0
      DO i=1,nrep
        no_offs(i+1) = no_offs(i) + no(i)
        no_sum = no_sum + no(i)
      ENDDO
!     write(lupri,*) 'offset table:',no_offs
!     write(lupri,*) 'total number of occupied orbitals:',no_sum

!  allocate space for Delta_mn matrix and clear matrix

      allocate(dmn(no_sum,no_sum))
      allocate(dmnext(no_sum,no_sum,2))
      allocate(xmna(no_sum,no_sum))
      allocate(xmnc(no_sum,no_sum))
      allocate(xmnd(no_sum,no_sum))
      allocate(xmn2(no_sum,no_sum))

!  open DADCPOP file for sequential access and write 
!  number of matrices and corresponding matrix dimension
!
!  ATT We check if file exists already and append the next symmetry
!  to it. This liberates us from writing symmetry-specific files altogether.
 
      inquire(file=dmn_filename,exist=isthere)

      select case(isthere)

      case (.FALSE.)      ! Delta_mn file was not created
      
        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='new',         &
             access='sequential',  &
             form='unformatted')

      case (.TRUE.)       ! Delta_mn file is there.

        open(unit=IO_LAN,          & 
             file=dmn_filename,    &
             status='old',         &
             access='sequential',  &
             form='unformatted',   &
             position='append')

      end select

!  write relevant data in this symmetry to the DNM file
!  symmetry, # of analyzed vectors, matrix dimension

      write(IO_LAN) ibeta, iwwcount, no_sum

!__________________________________________________________________
!|
!|  Loop over all vectors in the 2h/2h block possessing high enough 
!|  pole strength

      DO ivcount=1,nocc 

        IF(polestr(ivcount).lt.polethr) cycle ! too small, no POPANA

        xmna = A0; xmn2 = A0; xmnc = A0; xmnd = A0

        write(lupri,*) '   ** doing analysis for vector:',ivcount
        write(lupri,*) '   ** eigenvalue: ',cevl(ivcount)*autoev

!
!  found a vector, construct x_mn matrix
!
        ipoint = 0
        DO 14 JRP = 1, NREP
          IRP = MULTB(JRP,ibeta+NREP,2)
          IF (IRP.LT.JRP) GOTO 14
          DO J = 1, NO(JRP)
            IMIN = 1
            IF (IRP.EQ.JRP) IMIN = J + 1
            DO I = IMIN, NO(IRP)
              ipoint = ipoint + 1   ! get next eigenvector component
              irow   = no_offs(irp) + I
              jcol   = no_offs(jrp) + J
              xmna(irow,jcol) = thvecs(ipoint,ivcount)
!             write(*,*) '   read in components:',thvecs(ipoint,ivcount)
!             write(*,*) 'check: irow,icol:',irow,jcol ! irow must always be greater than jcol
!             write(*,'(A,4I4,A,2I4,F10.6)') & 
!                   ' irp,jrp,i,j',irp,jrp,i,j, &
!                   ' row/col/entry:',irow,jcol,xmna(irow,jcol)
            ENDDO
          ENDDO 
  14    CONTINUE

        IF( (ipoint).NE.NOOT(ibeta)) THEN
          WRITE(lupri,*) 'ipoint,noot',ipoint,noot(ibeta)
          CALL QUIT('Error in TWOHPOP for 2h configurations!')
        ENDIF

!       write(lupri,*) 'Contributions from xmna (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmna(i,j)),j=1,no_sum)
!       enddo
!       write(lupri,*) 'Contributions from xmna (imag):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (AIMAG(xmna(i,j)),j=1,no_sum)
!       enddo

!  contract x_ml and x*_nl. XGEMM can be used because all irrelevant entries are
!  filled with zeros and will not contribute to the matrix product.
!  Second call: contract x_lm and x_ln. Transpose the first matrix.

! first part: contract x_ml . x_ln* over l<m,n
! result goes into xmnc

        xmn2 = transpose(xmna)  ! form x_ln from x_nl. 
        xmn2 = dconjg(xmn2)     ! form x_ln* from x_nl. 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmna,no_sum, &
                    xmn2,no_sum,A0,xmnc,no_sum)

!       write(lupri,*) 'Contributions from xmnc (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmnc(i,j)),j=1,no_sum)
!       enddo
!
! second part: form x_lm from x_ml (in xmna) and then
! form x_ln* result goes into xmnd
!
        xmn2 = dconjg(transpose(xmna))  ! form x_ln* from x_ml. 
        CALL XGEMM ('N','N',no_sum,no_sum,no_sum,A1,xmn2,no_sum, &
                    xmna,no_sum,A0,xmnd,no_sum)
        xmnd = transpose(xmnd)
        xmn2 = xmnc + xmnd
!
!
!       write(lupri,*) 'Contributions from part1+2 (real):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (REAL(xmn2(i,j)),j=1,no_sum)
!       enddo
!       write(lupri,*) 'Contributions from part1+2 (imag):'
!       do i=1,no_sum
!         write(lupri,'(40F8.4)') (AIMAG(xmn2(i,j)),j=1,no_sum)
!       enddo
!       stop


!  now form Delta_mn matrix. Part 1 is just one on the diagonal,
!  part 2 is just minus the pole strength on the dagonal.
        
        dmn = 0.0d0
        DO i=1,no_sum
          dmn(i,i) = A1 - DCMPLX(polestr(ivcount),0.0d0)
        ENDDO
        dmn = dmn + xmn2

        write(lupri,*) 'Delta_mn (real part):'
        do i=1,no_sum
          write(lupri,'(40F10.5)') (REAL(dmn(i,j)),j=1,no_sum)
        enddo
        write(lupri,*) 'Delta_mn (imag part):'
        do i=1,no_sum
          write(lupri,'(40F10.5)') (AIMAG(dmn(i,j)),j=1,no_sum)
        enddo
        write(lupri,*)
        write(lupri,'(A)') '------------------------------------------'
        write(lupri,*)

!  write delta_mn matrix to file in the format requested by dirac conversion
!  routine. dmnext can be a real array since real part is first layer and 
!  imaginary part is second layer written out.

        dmnext=0.0d0
        dmnext(:,:,1)=real(dmn(:,:))
        dmnext(:,:,2)=aimag(dmn(:,:))
        write(io_lan) dmnext

      ENDDO   !ivcount

!  perform all deallocations in reverse order of the allocations
!  and close file

!     write(lupri,*) 'entering deallocation section'

      deallocate(xmn2)
      deallocate(xmnd)
      deallocate(xmnc)
      deallocate(xmna)
      deallocate(dmnext)
      deallocate(dmn)
      deallocate(no_offs)
      deallocate(polestr)
      deallocate(cevl)
      deallocate(evec)
      deallocate(thvecs)

!     write(lupri,*) 'deallocation finished'

      close(io_lan)

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine DNMMATCHK(lupri,io_lan)
!
      implicit none
      integer :: lupri    ! file handle for printing to stdout
      integer :: io_lan   ! file handle start
      logical :: errflg   ! error flag
!
!  local variables
!
      real*8,  allocatable  :: dmnext(:,:,:)
      integer               :: i,j,ivcount,ibeta,iwwcount,no_sum
      character*7           :: dmn_filename='DADCPOP'

      call PST('Performing DMN matrix check+')

!  open DNM matrix and read in everything

      open(unit=io_lan,             &
           file=dmn_filename,       &
           status='old',            &
           access='sequential',     &
           form='unformatted')

      errflg=.false.
 12   read(io_lan,end=800,err=900) ibeta,iwwcount,no_sum
!     write(lupri,*) 'Record:',ibeta,iwwcount,no_sum
      allocate(dmnext(no_sum,no_sum,2))
      dmnext = 0.0d0
      do ivcount=1,iwwcount
        read(io_lan) dmnext
      enddo

!     write(lupri,*) 'Delta_mn (real):'
!     do i=1,no_sum
!       write(lupri,'(40F10.5)') (dmnext(i,j,1),j=1,no_sum)
!     enddo
!     write(lupri,*) 'Delta_mn (complex):'
!     do i=1,no_sum
!       write(lupri,'(40F10.5)') (dmnext(i,j,2),j=1,no_sum)
!     enddo

      write(lupri,*) 'Symmetry',ibeta,' passed.'
      deallocate(dmnext)
      goto 12

 800  return

 900  errflg=.false.
      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION DOT_PRODUCTZ(a,b)
      complex*16 dot_productz
!
!  This function computes
!        Dot_productz = conjg(a) * b
!  and can have complex results for a.ne.b
!
      complex*16, Dimension(:), intent(in)    ::  a,b
!
!  local variables
!
      complex*16 :: ca1
      integer    :: i,n1,n2,n

      n1 = size(a)
      n2 = size(b)
      if(n1.ne.n2) STOP 'Vectors have different length in dot_productz!'
      n = n1
!     write(*,*) 'Vector length in dot_productz:',n

      ca1=(0.0d0,0.0d0)
      DO i=1,n
        ca1 = ca1 + conjg(a(i))*b(i)
      ENDDO

      write(*,*) 'Total real part:',dble(ca1)
      write(*,*) 'Total imag part:',imag(ca1)

      dot_productz = ca1
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION vnorm_cmplx(a)
      real*8 vnorm_cmplx
!
!  This function computes
!        vnorm_cmplx = conjg(a) * a
!  as a real value.
!
      complex*16, Dimension(:), intent(in)    ::  a
!
!  local variables
!
      complex*16 :: ca1
      real*8     :: rsum
      integer    :: i,n


      rsum=0.0d0

      n = size(a)
!     write(*,*) 'Vector length in vnorm_cmplx:',n

      DO i=1,n
        ca1 = conjg(a(i))*a(i)
        rsum = rsum + dble(ca1)
      ENDDO

      vnorm_cmplx = rsum
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
