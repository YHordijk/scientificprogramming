C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_X_MAIN (daibj,eps,dimphph)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Computes PH/PH block (main block) of the Excitation ADC Matrix
!     (polarization propagator). 
!     Works for real and complex case.
!
!---------------Calling variables--------------------------------------
!
!
      REAL*8,dimension(:)         :: daibj
      REAL*8,dimension(:)         :: eps
      INTEGER                     :: dimphph
!
!---------------Common Blocks, symmetry arrays from CCSETI --------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
!
!---------------Local variables--------------------------------------
!
      interface
        SUBROUTINE denomvvoo_pp(ra1,ra2)
          real*8, dimension(:)  :: ra1,ra2
        END SUBROUTINE
      end interface
!
!---------------Local variables--------------------------------------
!
      real*8, allocatable        :: sa_vvoott(:)
      real*8, allocatable        :: sb_vvoott(:)
      real*8, allocatable        :: buf1(:)
      real*8, allocatable        :: buf2(:)
      real*8, allocatable        :: buf3(:)
      real*8, allocatable        :: buf4(:)
      real*8, allocatable        :: buf5(:)
      real*8, allocatable        :: bab(:)
      real*8, allocatable        :: dij(:)
      real*8, allocatable        :: vovo1(:)

      integer                    :: ialloc

      character*31               :: ierror=
     &                              'Internal error in BUILD_X_MAIN!'
!
!---------------Executable code--------------------------------------
!
      IF(CARITH) THEN
        CALL PST('Building complex PH/PH block (all symmetries)+')
      ELSE
        CALL PST('Building real PH/PH block (all symmetries)+')
      ENDIF
      if(carith.and.(rcw.ne.2)) then
        call quit(ierror)
      endif
      if(.not.carith.and.(rcw.ne.1)) then
        call quit(ierror)
      endif
!
!******************* preparations
!
!  often needed v_ac,kl array which is V_ac,kl (VVOO)
!  divided by the energy denominator
!  e_a + e_c - e_k - e_l and the original V_ac,kl
!  we therefore calculate it once and store it in this subroutine.
!  the variables icX count the number of integrals. For the actual
!  storage requirements one has to multiply with RCW!
!
      ic1 = ivvoott(nrep+1)
      allocate(sa_vvoott(ic1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'allocation problem 1 in BUILD_X_MAIN'
      sa_vvoott = 0.0d0
      allocate(sb_vvoott(ic1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'allocation problem 2 in BUILD_X_MAIN'
      sb_vvoott = 0.0d0

      call getvvoo(sa_vvoott)
      call xcopy(ic1,sa_vvoott,1,sb_vvoott,1)
      call denomvvoo_pp(eps,sa_vvoott)
!
! we have v_ac,kl in sa_vvoott and
! V_ac,kl         in sb_vvoott
!
!  allocate remaining arrays
!  BAB: since the delta_ij enforces gamma_i = gamma_j only the totally
!  symmetric entries in BAB will survive!
!
      ic2 = NFVV
      allocate(bab(ic2*rcw))
      bab = 0.0d0

      ic3 = max(lvvoot(nrep+1),kvvoot(nrep+1))
      allocate(buf1(ic3*rcw))  ! array has to accomodate A,BIJ and ABI,J storage
      allocate(buf2(ic3*rcw))
      buf1 = 0.0d0; buf2 = 0.0d0

      ic4 = NFOO
      allocate(dij(ic4*rcw)); 
      dij = 0.0d0

      ic5 = max(j2vovo(nrep+1),jvovo(nrep+1))
      allocate(buf3(ic5*rcw))  ! all arrays have to hold VOVO integrals
      allocate(buf4(ic5*rcw))  ! in MVO*MVO instead of NVO*NVO manner!
      allocate(buf5(ic5*rcw))
      buf3 = 0.0d0; buf4 = 0.0d0; buf5 = 0.0d0

      ic6 = ivovo(nrep+1)      ! size of the generic VOVO array in generic BKC
      allocate(vovo1(ic6*rcw))  ! in MVO*MVO instead of NVO*NVO manner!
      vovo1 = 0.0d0
!
!
!_____________________________________________________
!|
!| Calculate B_ab intermediate   TERM C(2)A
!| result in bab  <a||b> BKC
!
!  sort v_ac,kl  -->  v_a,ckl    (buf1)
!  sort V*_bc,kl  -->  V*_b,ckl  (buf2), c.c. will bw accounted for
!                                in the corresponding contrct call!
!  v_a,ckl  *  V*_b,ckl   -->  bab   <||>
!  'hermitize'  the result in bab (add the transposed cc)
!  **att** the array bab is packed as nv1*nv1(beta1),...
!
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             sa_vvoott,buf1)
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             sb_vvoott,buf2)
      call cntrct ('N','C',nv,nv,nvoot,a1,buf1,buf2,a1,bab,nrep)
      do krep=1,nrep
        ioff=jjvv(krep,krep)*rcw + 1
        if(carith) then
          call cmatsym(bab(ioff),nv(krep),nv(krep))
        else
          call matsym(bab(ioff),nv(krep),nv(krep))
        endif
      enddo
!|
!|____________________________________________________
!
!_____________________________________________________
!|
!| Calculate D_ij intermediate   TERM C(2)B
!| result in dij  |i><j|
!
!  sort v_cd,ki   -->  v_cdk,i  <c|<d||k>,|i>  (buf1)
!  V_cd,kj        -->  V*_cd,kj
!  sort V*_cd,kj   -->  V*_cdk,j  |c>|d><k|,<j|  (buf2)
!  v_cdk,i  *  V*_cdk,j  -->  Dij   -->  'hermitize'
!  Dij   has |i><j| 
!
      buf1 = 0.0d0; buf2 = 0.0d0
      call srt1t3 (nrep,multb,lfa,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &             sa_vvoott,buf1)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call srt1t3 (nrep,multb,lfa,nvvt,no,no,nvvot,kvvoot,kkvvot,
     &             sb_vvoott,buf2)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call cntrct ('T','N',no,no,nvvot,a1,buf1,buf2,a0,dij,nrep)
      do krep=1,nrep
        ioff=jjoo(krep,krep)*rcw + 1
        if(carith) then
          call cmatsym(dij(ioff),no(krep),no(krep))
        else
          call matsym(dij(ioff),no(krep),no(krep))
        endif
      enddo
!|
!|____________________________________________________
!
! ****************** why this conjugation ???? ***************
      if(carith) call conjuga(ic4,dij,1)
!
!_____________________________________________________
!|
!|  Calculate TERM C(2)C, part[A] herm.    |
!  conjugation required, result in buf3   |
!
!  sort v_ac,ik    -->  v_ai,ck    (buf3)
!  V*_bc,jk = V*_cb,kj
!  sort V*_cb,kj   -->  V*_ck,bj   (buf4)
!
      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,
     &              mvo,j2vovo,jjvo,jjvo,sa_vvoott,buf3)
      if(carith) call conjuga(ic1,sb_vvoott,1)
      call srt1tt4 (nrep,multb,lfa,lfa,nv,nv,no,no,
     &              mvo,j2vovo,jjvo,jjvo,sb_vvoott,buf4)
      if(carith) call conjuga(ic1,sb_vvoott,1)
!
!      CALL CNTRCT('N','T',MVO,MVO,MVO,A1,buf3,buf4,A0,buf5,NREP)
!
!
!  here again special contraction of v_ai,ck  X  V*_ck,bj
!  buf3  *  buf4. Result has BKC:  <a||i>,|b><j| as the storage mode.
!
      ioff1 = 1
      ioff3 = 1
      do irep = 1, nrep
         jrep = multb(irep+nrep,1+nrep,2)
         if(multb(jrep+nrep,irep+nrep,1).ne.1)
     &      stop 'serious symmetriy violation'
         m = mvo(irep)
         n = mvo(irep)
         k = mvo(jrep)
         ioff2 = j2vovo(jrep) * rcw + 1
         call xgemm ('N','N',m,n,k,a1,buf3(ioff1),m,buf4(ioff2),n,
     &               a0,buf5(ioff3),m)
         ioff1 = ioff1 + m * k * rcw
         ioff3 = ioff3 + m * n * rcw
      enddo
!
! This is now the first half contribution to C(2)C. Part[B] is the hermitian
! conjugated one. Multiply result by -0.5 according to the equations.
!
      do krep=1,nrep
        ioff=jvovo(krep)*rcw + 1
        if(carith) then
          call cmatsym(buf5(ioff),mvo(krep),mvo(krep))
        else
          call matsym(buf5(ioff),mvo(krep),mvo(krep))
        endif
      enddo
      call xscal(jvovo(nrep+1),-ap5,buf5,1)
!
!  buf5 now contains C(2)C in <a||i>,|b><j| storage mode (Y_ai,bj).
!  Do intermediate resort of Y_ai,bj (buf5) ---->  Y_aj,bi  (buf3) with *generic* VOVO BKC.
!  in SRT16 the buffer order is inverted for inverse sorting.
!
      buf3 = 0.0d0
      CALL SRT16 (nrep,multb,lfa,ltr,nv,no,nv,no,mvo,jvovo,jjvo,jjvo,
     &            buf3,buf5)
!
!  ----------------------------------
!  Then we can add C(2)C and C(1)   |
!  due to identical BKC structure   |
!  ----------------------------------
!
      call getvovo(vovo1)
      call xscal(ic6,-A1,vovo1,1)
      call xaxpy (ic6,A1,buf3,1,vovo1,1)
!
!--------------------------------------------------------
!  'weave' into the C(0), C(2A) and C(2B) contributions |
!  as diagonal energies, 0.5*Bab and 0.5*Dij            |
!--------------------------------------------------------
!
      if(carith) then
        call bweaverc(vovo1,eps,bab,dij,ic6,ic2,ic4)
      else
        call bweaverr(vovo1,eps,bab,dij,ic6,ic2,ic4)
      endif
!
!--------------------------------------------------------
!  vovo1 now contains all the C(2) ai,bj contributions  |
!  we sort it back now to the PHPH storage convention.  |
!  This is now a forward sort <<>>  -->  <>><           |
!--------------------------------------------------------
!
      CALL SRT16(NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,
     &           JJVO,JJVO,vovo1,daibj)

      deallocate(vovo1)
      deallocate(buf5)
      deallocate(buf4)
      deallocate(buf3)
      deallocate(dij)
      deallocate(buf2)
      deallocate(buf1)
      deallocate(bab)
      deallocate(sb_vvoott)
      deallocate(sa_vvoott)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BWEAVERR(bphph,eps,bab,dij,ic6,ic2,ic4)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Crawl through the VOVO structure and distribute the
!     C(0), C(2A) and C(2B) contributions to the phph block
!
!---------------Calling variables--------------------------------------
!
      integer, intent(in)         :: ic6, ic2, ic4
      real*8, intent(inout)       :: bphph(ic6)
      real*8, intent(in)          :: eps(*)
      real*8, intent(in)          :: bab(ic2)
      real*8, intent(in)          :: dij(ic4)
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8                           ::  esum
!
!---------------Executable code--------------------------------------
!
!  the loop structure follows the generic VO||VO storage with <|<|  |>|>
!
      AIBJ = 0
      DO BJRP = 1, NREP
      DO 10 JRP = 1, NREP
       BRP = MULTB(JRP,BJRP+NREP,2)
       DO J = 1, NO(JRP)
       DO B = 1, NV(BRP)
         DO 20 IRP = 1, NREP
          ARP = MULTB(IRP,BJRP+NREP,2)
          DO I = 1, NO(IRP)
          DO A = 1, NV(ARP)
            AIBJ = AIBJ + 1  ! count elements
!
!  from here we know arep,irep,brep,jrep, a,i,b,j
!
!
!  treat the diagonal element which is C(0)=(e_a-e_i)*d_ab*d_ij
!
            IF(ARP.eq.BRP.and.IRP.eq.JRP.and.a.eq.b.and.i.eq.j) THEN
              AOFF = IO(NREP+1) + IV(ARP) + A
              IOFF = IO(IRP) + I
              esum = EPS(AOFF) - EPS(IOFF)
              bphph(aibj) = bphph(aibj) + esum
            ENDIF
!
!  treat the C(2A) contribution: C(2A) = +0.5*Delta_ij*B_ab
!
            IF(IRP.eq.JRP) THEN
              IF(I.eq.J) THEN
                if(ARP.NE.BRP) call quit('bweaverr_error C(2A)')
                aboff = IVV(BRP) + (A-1)*NV(BRP) + B
                bphph(aibj) = bphph(aibj) + 0.5d0*bab(aboff)
              ENDIF
            ENDIF
!
!  treat the C(2B) contribution: C(2B) = +0.5*Delta_ab*D_ij
!
            IF(ARP.eq.BRP) THEN
              IF(A.eq.B) THEN
                if(IRP.NE.JRP) call quit('bweaverr_error C(2B)')
                ijoff = IOO(JRP) + (I-1)*NO(JRP) + J
                bphph(aibj) = bphph(aibj) + 0.5d0*dij(ijoff)
              ENDIF
            ENDIF
!
!  all contributions "woven into"
!
            ENDDO
            ENDDO
 20      CONTINUE
       ENDDO
       ENDDO
 10   CONTINUE
      ENDDO  !BJRP
!
      if(aibj.ne.ivovo(nrep+1)) then
        call quit('VOVO counting error in BWEAVERR!')
      endif

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BWEAVERC(bphph,eps,bab,dij,ic6,ic2,ic4)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Count through the VOVO integrals and distribute the
!     C(0), C(2A) and C(2B) contributions to the phph block
!
!---------------Calling variables--------------------------------------
!
      integer, intent(in)         :: ic6, ic2, ic4
      complex*16, intent(inout)   :: bphph(ic6)
      real*8, intent(in)          :: eps(*)
      complex*16, intent(in)      :: bab(ic2)
      complex*16, intent(in)      :: dij(ic4)
      complex*16,parameter        :: coha=(0.5d0,0.0d0)
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8                           ::  esum
!
!---------------Executable code--------------------------------------
!
!  the loop structure follows the generic VO||VO storage with <|<|  |>|>
!
      AIBJ = 0
      DO BJRP = 1, NREP
      DO 10 JRP = 1, NREP
       BRP = MULTB(JRP,BJRP+NREP,2)
       DO J = 1, NO(JRP)
       DO B = 1, NV(BRP)
         DO 20 IRP = 1, NREP
          ARP = MULTB(IRP,BJRP+NREP,2)
          DO I = 1, NO(IRP)
          DO A = 1, NV(ARP)
            AIBJ = AIBJ + 1  ! count elements
!
!  from here we know arep,irep,brep,jrep, a,i,b,j
!
!
!  treat the diagonal element which is C(0)=(e_a-e_i)*d_ab*d_ij
!
            IF(ARP.eq.BRP.and.IRP.eq.JRP.and.a.eq.b.and.i.eq.j) THEN
              AOFF = IO(NREP+1) + IV(ARP) + A
              IOFF = IO(IRP) + I
              esum = EPS(AOFF) - EPS(IOFF)
              bphph(aibj) = bphph(aibj) + dcmplx(esum,0.0d0)
            ENDIF
!
!  treat the C(2A) contribution: C(2A) = +0.5*Delta_ij*B_ab
!
            IF(IRP.eq.JRP) THEN
              IF(I.eq.J) THEN
                if(ARP.NE.BRP) call quit('bweaverc_error C(2A)')
                aboff = IVV(BRP) + (A-1)*NV(BRP) + B
                bphph(aibj) = bphph(aibj) + coha*bab(aboff)
              ENDIF
            ENDIF
!
!  treat the C(2B) contribution: C(2B) = +0.5*Delta_ab*D_ij
!
            IF(ARP.eq.BRP) THEN
              IF(A.eq.B) THEN
                if(IRP.NE.JRP) call quit('bweaverc_error C(2B)')
                ijoff = IOO(JRP) + (I-1)*NO(JRP) + J
                bphph(aibj) = bphph(aibj) + coha*dij(ijoff)
              ENDIF
            ENDIF
!
!  all contributions "woven into"
!
            ENDDO
            ENDDO
 20      CONTINUE
       ENDDO
       ENDDO
 10   CONTINUE
      ENDDO  !BJRP
!
      if(aibj.ne.ivovo(nrep+1)) then
        call quit('VOVO counting error in BWEAVERC!')
      endif

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SOLVE_PHPH_R(daibj,n,iounit)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Diagonalizes the phph block in each symmetry separately.
!     For the singles amplitudes a local configuration analysis is
!     provided for convenience.
!
!---------------Calling variables--------------------------------------
!
      real*8,dimension(:)         :: daibj
      integer                     :: n,iounit
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      real*8, allocatable         :: amat(:,:),xmat(:)
      integer                     :: desrep,jrep,brep,trep
      integer                     :: ioff,nsiz,ialloc,icount
      integer                     :: nx,lfun,kfun,jabs,babs
      integer                     :: resh_array(2)
      real*8, allocatable         :: evl(:),work(:)
      integer                     :: lwork,info
      real*8,parameter            :: AUTOEV = 27.2113957D0
      real*8                      :: x2
      character*26,allocatable    :: cfield(:)
      character*11                :: davstartname
      character*8                 :: dsroot='DAVSTART'
!
!---------------Interface region-------------------------------------
!
      interface

        function qmasymr(ra1,i1)
          real*8                 :: qmasymr
          real*8,dimension(*)    :: ra1
          integer                :: i1
        end function qmasymr

      end interface
!
!---------------Executable code--------------------------------------
!
      call PST('Generating (real) Dav. start vectors (all symm.)+')
!
!  create local single configuration descriptor array
!  remember: first index is virtual^*, second is occupied
!  this corresponds to the storage mode!
!
      nx = jvoxx(nrep+1)
      allocate(cfield(nx))

      icount = 0
      do desrep=1,nrep   !desrep is a bosonic loop
        do jrep = 1,nrep   !jrep,brep are fermionic
          brep = multb(desrep+nrep,jrep,2)
          if(multb(brep,jrep,2).ne.desrep) then
            write(*,*) 'PH/PH ordering problem!'
            stop 'solve_phph_r'
          endif
          do jfun = 1,NO(jrep)
          do bfun = 1,NV(brep)
            jabs=IO(jrep) + jfun
            babs=IO(nrep+1) + IV(brep) + bfun
            icount = icount + 1
            write(cfield(icount),777) jabs,repna(jrep),
     &                                babs,repna(brep)
          enddo
          enddo
        enddo
      enddo
 777  Format(I3,' (',A4,') --->',1X,I3,' (',A4,')')
      if(icount.ne.jvoxx(nrep+1)) then
         write(*,*) icount,jvoxx(nrep+1)
         stop 'conf array error!'
      endif
!     write(*,*) 'HP configuration counting successful'
!
!  diagonalize phph submatrices
!  and write corresponding orthogonal eigenvectors
!  as starting vectors for the Davidson algorithm.
!
!___________________________________________________________
!|
!|
!|
      do desrep = 1,nrep

        ioff = JVOVO(desrep)   ! calculate total offset in corresponding irrep
        nsiz = MVO(desrep)
        if(nsiz.eq.0) then
          write(*,'(1X,A27,6X,I3,3X,A4)') 'No states in symmetry',
     &          desrep,repna(desrep)
          cycle
        endif
        write(*,*)
!
!  open corresponding file
!
        IF(desrep.GT.9) THEN
          WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
        ELSE
          WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
        ENDIF
        OPEN(iounit,file=davstartname,access='sequential',
     &       status='unknown',form='unformatted')

        allocate(amat(nsiz,nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 1 in SOLVE_PHPH_R'
        allocate(xmat(nsiz*nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 2 in SOLVE_PHPH_R'
        amat = 0.0d0; xmat = 0.0d0
!
!  we convert the linear MVO*MVO array to quadratic array
!  hereby we need xmat, compiler does not understand it otherwise
!
        resh_array(1)=nsiz; resh_array(2)=nsiz
        xmat(1:nsiz*nsiz) = daibj(ioff+1:ioff + nsiz*nsiz)
        amat = reshape(xmat,resh_array)
!
!  check asymmetry of obtained PHPH matrix
!
        x2 = dabs(qmasymr(amat,nsiz))
        if(x2.gt.1.0E-10) then
          write(*,*) 'Warning: Asymmetry of PH/PH matrix',
     &         ' in symmetry',desrep,' amounts to',x2
        endif
!
!  amat array has ph/ph structure and is transferred to LAPACK
!
        lwork=5*nsiz
        allocate(evl(nsiz));allocate(work(lwork))
        CALL DSYEV('V','L',nsiz,amat,nsiz,evl,work,lwork,info)
        if(info.ne.0) then
          write(*,*) 'DSYEV problem!  Info code:',info
          stop
        endif
!
!  write out eigenvalues and largest eigenvectors for the current symmetry
!  on screen.
!
!       write(*,*)
        write(*,'(1X,A27,6X,I3,3X,A4,A)') 'Vectors available for symm.',
     &            desrep,repna(desrep),'  (max. 10 are printed)'
        write(*,*)
        do i=1,min(nsiz,10)
          write(*,'(I5,2F22.12,2X,A4,X,A2)') i,evl(i),evl(i)*autoev,
     &             repna(desrep),'@X'
!         write(*,'(5X,A)') '-----------------------------------------'
!         do j=1,nsiz
!           if(dabs(amat(j,i)).gt.0.1d0) then
!             ioff=jvoxx(desrep)
!             write(*,*) j,':  ',cfield(ioff+j),amat(j,i)
!           endif
!         enddo
!         write(*,*)
        enddo
!
!  write out start vectors to specific file name.
!
        Write(iounit) nsiz
        Do i=1,nsiz
          Write(iounit) (amat(j,i),j=1,nsiz)
        Enddo
        Close(iounit)
        
!
!  deallocate all auxiliary arrays and go to next symmetry
!
        deallocate(work)
        deallocate(evl)
        deallocate(xmat)
        deallocate(amat)

      enddo  !desrep
      deallocate(cfield)
      write(*,*) '...success.'
!|
!|
!|
!|__________________________________________________________


      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SOLVE_PHPH_C(daibj,nele,iounit)
!
      IMPLICIT INTEGER (A-Z)
!
!---------------Description--------------------------------------------
!
!     Diagonalizes the hermitian PHPH block in each symmetry separately.
!     For the singles amplitudes a local coniguration analysis is
!     provided for convenience.
!     Due to the more restrictive F90 compiler with function prototypes
!     we can not reinterpret the real array as a complex one. We introduce
!     an *intermediate complex* array. The main block size is small enough
!     not to run into trouble.
!     nele is the number of elements in the array and *not* the physical
!     extension of the real array. The real array has an allocation of
!     2*nele !!
!
!---------------Calling variables--------------------------------------
!
      real*8,dimension(:)         :: daibj
      integer                     :: nele,iounit
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
!
!---------------Local variables--------------------------------------
!
      complex*16, allocatable, dimension(:,:)      :: amat
      complex*16, allocatable, dimension(:)        :: xmat,daibjc
      integer                             :: desrep,jrep,brep,trep
      integer                             :: ioff,nsiz,ialloc,icount
      integer                             :: nx,lfun,kfun,jabs,babs
      integer                             :: resh_array(2)
      integer                             :: i,ix
      real*8, allocatable, dimension(:)        :: evl,rwork
      complex*16, allocatable, dimension(:)    :: work
      integer                             :: lwork,info
      real*8,parameter                    :: AUTOEV = 27.2113957D0
      real*8                              :: x2
      character*26,allocatable, dimension(:)   :: cfield
      character*11                        :: davstartname
      character*8                         :: dsroot='DAVSTART'

      real*8,parameter        :: rzero=(0.0d0)
      complex*16,parameter    :: czero=(rzero,rzero)
!
!---------------Interface region-------------------------------------
!
      interface

        function qmasymh(ca1,i1)
          real*8                       :: qmasymh
          complex*16,dimension(:,:)    :: ca1
          integer                      :: i1
        end function qmasymh

        subroutine ZHEEV(c1,c2,i1,ca1,i2,ra1,ca2,i3,ra2,i4)
           character*1                         :: c1,c2
           integer                             :: i1,i2,i3,i4
           double complex, dimension(i2,*)     :: ca1
           double precision, dimension(*)      :: ra1,ra2
           double complex, dimension(*)        :: ca2
        end subroutine

      end interface
!
!---------------Executable code--------------------------------------
!
      call PST('Generating (complex) Dav. start vectors (all symm.)+')

      if(rcw.ne.2) then
        write(*,*) 'Variable RCW is not 2 despite complex mode!'
        call quit('*** internal consistency violated ***')
      endif
!
!  create local single configuration descriptor array
!  remember: first index is virtual^*, second is occupied
!  this corresponds to the storage mode!
!
      nx = jvoxx(nrep+1)
      allocate(cfield(nx))  !character field
!___________________________________________
!|
!|  configuration generator
!|
      icount = 0
      do desrep=1,nrep   !desrep is a bosonic loop
        do jrep = 1,nrep   !jrep,brep are fermionic
          brep = multb(desrep+nrep,jrep,2)
          if(multb(brep,jrep,2).ne.desrep) then
            write(*,*) 'PH/PH ordering problem!'
            stop 'SOLVE_PHPH_C'
          endif
          do jfun = 1,NO(jrep)
          do bfun = 1,NV(brep)
            jabs=IO(jrep) + jfun
            babs=IO(nrep+1) + IV(brep) + bfun
            icount = icount + 1
            write(cfield(icount),777) jabs,repna(jrep),
     &                                babs,repna(brep)
          enddo
          enddo
        enddo
      enddo
 777  Format(I3,' (',A4,') --->',1X,I3,' (',A4,')')
      if(icount.ne.jvoxx(nrep+1)) then
         write(*,*) icount,jvoxx(nrep+1)
         stop 'conf array error!'
      endif
!|
!|  end configuration generator
!|__________________________________________
!
!  diagonalize PHPH submatrices
!  and write corresponding complex orthogonal eigenvectors
!  as starting vectors for the Davidson algorithm.
!
!___________________________________________________________
!|
!|  generate complex intermediate array from original real one.

      allocate(daibjc(nele))
      daibjc = czero
      do i=1,nele
        ix = 2*(i-1) + 1
        daibjc(i)=dcmplx(daibj(ix),daibj(ix+1))
      enddo
!|
!|__________________________________________
!
!___________________________________________________________
!|
!| form quadratic matrix in specific symmetry
!| from total linear array and diagonalize
!|
      do desrep = 1,nrep

        ioff = JVOVO(desrep)   ! calculate total offset in corresponding irrep
        nsiz = MVO(desrep)
        if(nsiz.eq.0) then
          write(*,'(1X,A27,6X,I3,3X,A4)') 'No states in symmetry',
     &          desrep,repna(desrep)
          cycle
        endif
!
!  open corresponding file
!
        IF(desrep.GT.9) THEN
          WRITE(davstartname,'(A8,A1,I2)') dsroot,'.',desrep
        ELSE
          WRITE(davstartname,'(A8,A2,I1)') dsroot,'.0',desrep
        ENDIF
        OPEN(iounit,file=davstartname,access='sequential',
     &       status='unknown',form='unformatted')

        allocate(amat(nsiz,nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 1 in SOLVE_PHPH_C'
        allocate(xmat(nsiz*nsiz),stat=ialloc)
        if(ialloc.ne.0) STOP 'Allocation problem 2 in SOLVE_PHPH_C'
        amat = czero
        xmat = czero
!
!  we convert the linear MVO*MVO array to quadratic array
!  hereby we need xmat, compiler does not understand it otherwise
!
        resh_array(1)=nsiz; resh_array(2)=nsiz
        xmat(1:nsiz*nsiz) = daibjc(ioff+1:ioff + nsiz*nsiz)
        amat = reshape(xmat,resh_array)
!
!  check deviation from hermiticity of obtained PHPH matrix
!
        x2 = dabs(qmasymh(amat,nsiz))
        if(x2.gt.1.0E-10) then
          write(*,*) 'Warning: Deviation from hermiticity of PH/PH',
     &       ' matrix in symmetry',desrep,' amounts to',x2
        else
          write(*,*) 'Hermiticity of PH/PH block:',x2
        endif
        write(*,*)
!
!  amat array has ph/ph structure and is transferred to LAPACK
!
        lwork=5*nsiz
        allocate(evl(nsiz))
        allocate(work(lwork))
        allocate(rwork(lwork))
        CALL ZHEEV('V','L',nsiz,amat,nsiz,evl,work,lwork,rwork,info)
        if(info.ne.0) then
          write(*,*) 'ZHEEV problem in SOLVE_PHPH_C! code:',info
          call quit('numerical problem')
        endif
!
!  write out eigenvalues and largest eigenvectors for the current symmetry
!  on screen.
!
        write(*,'(1X,A27,6X,I3,3X,A4,A)') 'Vectors available for symm.',
     &            desrep,repna(desrep),'  (max. 10 are printed)'
        write(*,*)
        do i=1,min(nsiz,10)
          write(*,'(I5,2F22.12,2X,A4,X,A2)') i,evl(i),evl(i)*autoev,
     &             repna(desrep),'@X'
        enddo
!
!  write out start vectors to specific file name.
!
        Write(iounit) nsiz
        Do i=1,nsiz
          Write(iounit) (amat(j,i),j=1,nsiz)
        Enddo
        Close(iounit)
        
!
!  deallocate all auxiliary arrays and go to next symmetry
!
        deallocate(work)
        deallocate(rwork)
        deallocate(evl)
        deallocate(xmat)
        deallocate(amat)

      enddo
      deallocate(cfield)
      deallocate(daibjc)
      write(*,*) '...success.'
!|
!|
!|
!|__________________________________________________________

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
