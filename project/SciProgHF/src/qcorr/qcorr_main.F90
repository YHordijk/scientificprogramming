!     ***********************************************************************
      subroutine Q_correction(C,NBLOCK,IBLOCK,IBLTP,LUC,                &
     &                        LUSCR,ICISTR,IRC,WORK,relQC)
!     incoming variables:
!                         C      - scratch space for ci vector
!                         NBLOCK - # of ttsblocks
!                         IBLOCK - array containing information about each ttsblock
!                         IBLTP  - block types
!                         LUC    - internal file pointer for CI vector file
!                         LUSCR  - file pointer for scratch vector file
!                         ICISTR - currently set to 2 by default
!                         IRC:   - algebra control value: 1 (real) or 2 (complex)
      use qcorr_cfg
      use qcorr_interface

      implicit none

!     input
      integer, intent(in)   :: IBLOCK(8,NBLOCK),IBLTP(*)
      real(8), intent(inout):: C(*), WORK(*)
      integer, intent(in)   :: nblock, luc, luscr, irc, icistr
      logical, intent(in)   :: relQC
!     scratch
      integer, allocatable  :: astring(:,:)
      integer, allocatable  :: bstring(:,:)
      integer, allocatable  :: l_ref(:,:)
      integer, allocatable  :: ref_order(:,:)
      real(8), allocatable  :: ref_wavefunction_coeff(:)
      real(8), allocatable  :: e_corr(:), e_tot(:)
      integer               :: nstates_ci
      real(8)               :: q_corr, edci
      real(8)               :: Edc_1, Erdc_1, Epc_1, Emc_1,Epcp_1,Emc_N1
      real(8)               :: Edc_2, Erdc_2, Epc_2, Emc_2,Epcp_2,Emc_N2
      real(8)               :: reference_energy, reference_energy_tmp
      integer               :: iactp, ibctp, nbl_c, nael, nbel
      integer               :: ntest, lblk, iboff, ionem
      integer               :: irir, imk2_r, i, j
      integer               :: la, lb, max_lref_len
      integer               :: idum, lucoef
      real(8)               :: dum
      logical               :: tobe = .false.

      call QENTER('QCORR')

!     initialize
      call init_qcorr(relQC,.true.)

!     print flag
      ntest = 000
      ntest = max(ntest,print_qcorr)

      allocate(ref_wavefunction_coeff(nref_qcorr))
      ref_wavefunction_coeff = 0.0d0 ! initialize all to 0.0d0

      allocate(l_ref(nref_qcorr,2))
      allocate(nash_qcorr(nfsym_q))
      nash_qcorr = 0

      ! read reference configurations from file (either written by MCSCF or by user input for HF reference)
      open(1234,file='refvec.luci',status='old',form='unformatted',     &
      access='sequential',action="read",position='rewind')

      read(1234) nref_qcorr,nash_qcorr(1:nfsym_q)
      read(1234) nref_e_qcorr
      allocate(reference_energy_qcorr(nref_e_qcorr))
      reference_energy_qcorr = 0
      read(1234) reference_energy_qcorr(1:nref_e_qcorr)
!     read # reference configurations
      allocate(ref_wavefunction_coeff_qcorr(nref_qcorr))
      allocate(ref_wavefunction_qcorr(nref_qcorr,2))
      max_lref_len = 0
      do i = 1, nref_qcorr
        read(1234) l_ref(i,1),l_ref(i,2),                               &
                         ref_wavefunction_qcorr(i,1)(1:l_ref(i,1)),     &
                         ref_wavefunction_qcorr(i,2)(1:l_ref(i,2)),     &
     &                   ref_wavefunction_coeff_qcorr(i)
!#define DEBUG_QCORR
#ifdef DEBUG_QCORR
         print *, ref_wavefunction_qcorr(i,1)(1:l_ref(i,1))
         print *, ref_wavefunction_qcorr(i,2)(1:l_ref(i,2))
         print *, ref_wavefunction_coeff_qcorr(i)
#endif
          max_lref_len = max(max_lref_len,l_ref(i,1),l_ref(i,2)) 
      
      end do
      
      allocate(ref_order(max_lref_len,2*nref_qcorr))
      ref_order = 0
      !> sort numbers from 1 ... x in reference strings such that we do not miss anything in linear symmetry (reordered spinors wrt mj)
      do j = 1, nref_qcorr
        do i = 1, l_ref(j,1)/4
          read(ref_wavefunction_qcorr(j,1)(1+(i-1)*4:4*i),'(i4)')       &
     &         ref_order(i,j)
        end do
        call sort(ref_order(1,j),l_ref(j,1)/4)
        do i = 1, l_ref(j,2)/4
          read(ref_wavefunction_qcorr(j,2)(1+(i-1)*4:4*i),'(i4)')       &
     &         ref_order(i,j+nref_qcorr)
        end do
        call sort(ref_order(1,j+nref_qcorr),l_ref(j,2)/4)
        !> write sorted strings
        ref_wavefunction_qcorr(j,1)(1:l_ref(j,1)) = ' '
        ref_wavefunction_qcorr(j,2)(1:l_ref(j,2)) = ' '
        write(ref_wavefunction_qcorr(j,1),"(60i4)")                     &
     &  (ref_order(i,j),i=1,l_ref(j,1)/4)
        write(ref_wavefunction_qcorr(j,2),"(60i4)")                     &
     &  (ref_order(i,j+nref_qcorr),i=1,l_ref(j,2)/4)
!#define DEBUG_QCORR
#ifdef DEBUG_QCORR
         print *, 'sorted reference ...'
         print *, ref_wavefunction_qcorr(j,1)
         print *, ref_wavefunction_qcorr(j,2)
         print *, ref_wavefunction_coeff_qcorr(j)
#endif

      end do

      deallocate(ref_order)


      close(1234,status="keep")

      if(.not.set_nash_q_qcorr)nash_q(1:nfsym_q)=nash_qcorr(1:nfsym_q)

      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call REWINE(LUC,-1)
      call itods(-1,1,-1,LUSCR)

!     loop over real and imaginary parts
      do IRIR = 1,IRC

        if (IRIR.EQ.2) then
!         skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if

        do IMK2_R = 1,NMS2VAL_q,1

          if(relQC)then
            IACTP = IST_FOR_DT_q(1,IMK2_R)
            IBCTP = IST_FOR_DT_q(2,IMK2_R)
          else
            IACTP = 1
            IBCTP = 2
          end if
!  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2_q(IMK2_R-1)
          end if
!
          if(relQC)then 
            NBL_C = NBLK_MS2_q(IMK2_R)
          else
            NBL_C = NBLOCK
          end if

          call REWINE(LUSCR,-1)
          call COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          call ITODS(-1,1,-1,LUSCR)

          NAEL = NELEC_q(IACTP)
          NBEL = NELEC_q(IBCTP)

          allocate(astring(NAEL,MXNSTR_q))
          allocate(bstring(NBEL,MXNSTR_q))
!
          call get_Q_correction(C,LUSCR,NBL_C,                          &
     &                          WORK(KNSTSO_q(IACTP)),                  &
     &                          WORK(KNSTSO2_q(IBCTP)),                 &
     &                          IBLOCK(1,IBOFF),                        &
     &                          NAEL,NBEL,astring, bstring,             &
     &                          IBLTP,NSMST_q,                          &
     &                          IACTP,IBCTP,ICISTR,ntest,               &
     &                          ref_wavefunction_qcorr,                 &
     &                          ref_wavefunction_coeff,nref_qcorr,      &
     &                          nash_qcorr,nash_q,nfsym_q,              &
     &                          NORB_q,relQC, wunit_q,l_ref,cvorb_qcorr)
!         NAEL/NBEL may be zero and therefore no allocation of string arrays...
          if(allocated(bstring)) deallocate(bstring)
          if(allocated(astring)) deallocate(astring)
        end do
      end do

      if(ntest > 2)then
        write(wunit_q,'(A)') '                                       '
        write(wunit_q,'(A)') '  Reference and CI coefficients       :' 
        write(wunit_q,'(A)') '  ====================================='

        do i = 1, nref_qcorr

          write(wunit_q,'(a,i4)')"  reference wave function (wf) #  ",i
          write(wunit_q,'(a,es13.6)')                                   &
     &    "     CI coefficient of the wf:  ",ref_wavefunction_coeff(i)
          write(wunit_q,'(a,es13.6)')                                   &
     &    "  MC/HF coefficient of the wf:  ",                           &
     &                                 ref_wavefunction_coeff_qcorr(i)
        end do
      end if

      ! read data for +Q correction
      open(file="energies.CI",unit=10,status="old",                     &
           form="unformatted",access="sequential")

      read(10) nstates_ci,reference_energy_tmp

      allocate(e_corr(nstates_ci))
      allocate(e_tot(nstates_ci))

      read(10) (e_tot(i),i=1,nstates_ci)
      close(10,status="keep")

      write(wunit_q,"(/A,/A,/A,/A,/A,/A)")                              &
     &                  "  The acronyms for each Q correction below"//  &
     &                  "  refer to: ",                                 &
     & "  Edc:  Davidson correction              [Eq. (34)]",           &
     & "  Erdc: renormalized Davidson correction [Eq. (36)]",           &
!fixme "  Epc:  Pople correction                 [Eq. (38)]",           &
     & "  Epc`: modified Pople correction        [Eq. (39)]",           &
     & "  Emc : Meissner correction              [Eq. (41)]"

      do i = 1, nstates_ci

        reference_energy = reference_energy_qcorr(i)
        e_corr(i)        = e_tot(i) - reference_energy
        
        call E_Q(e_corr(i),ref_wavefunction_coeff,nref_qcorr,           &
     &           ref_wavefunction_coeff_qcorr,NACTEL_q,                 &
     &           Edc_1, Erdc_1, Epc_1, Emc_1, Epcp_1, Emc_N1,           &
     &           Edc_2, Erdc_2, Epc_2, Emc_2, Epcp_2, Emc_N2,           &
     &           ntest, wunit_q)

        write(wunit_q,"(/A,i3,a)")'  Q corrections for state ',i,       &
     &  ' with c0^2 according to Eq. (43)'
      ! output of the q-correction       
        if(ntest > 2)then
          write(wunit_q,"( A,1f16.8)")'  reference   energy: ',         &
     &                                   reference_energy
          write(wunit_q,"( A,1f16.8)")'  correlation energy: ',e_corr(i)
          write(wunit_q,"( A,1f16.8)")'  original CI energy: ',e_tot(i)

        end if

        !> MC/HF wave function
        if(reference_energy == 0.0d0)then
          call quit('qcorr module: no reference energy given!')
        end if 

        CALL PRSYMB(wunit_q,'-',120,2)
        WRITE(wunit_q,'(2X,a12,2x,4(A16,1X),5x,a16)') 'model',          &
     &        'Edc ','Erdc ','Epc` ', 'Emc ','Emc (Nact - 2e)'
        WRITE(wunit_q,'(3X,a12,4x,5(f16.6,1X))')                        &
     &  'CI+Q energy ',e_tot(i)+Edc_1,e_tot(i)+Erdc_1,e_tot(i)+Epcp_1,  &
     &                 e_tot(i) + Emc_1, Emc_N1 + e_tot(i)
        CALL PRSYMB(wunit_q,'-',120,2)

!       now using eq. (44) for c0^2: c0^2 == <ref|ci>^2 == ( sum\limits_{p=1}^{nref} c_p^{(0)} * c_p )^2 
!       where c_p^{(0)} are the coefficients in the normalized 0th-order reference (MCSCF) wave function
        write(wunit_q,"(/A,i3,a)")'  Q corrections for state ',i,       &
     &  ' with c0^2 according to Eq. (44)'

!       write(wunit_q,"(/A)")  " Q corrections according to Eq. (44)"//
!    &                   " where c0^2 == <ref|ci>^2 == "//
!    &                   "( sum\limits_{p=1}^{nref} c_p^{(0)} * c_p )^2"
!       write(wunit_q,"(A/)") " Note: c_p^{(0)} are the reference"//
!    &                  " coefficients in the"//
!    &                  " normalized 0th-order reference (MCSCF) wave"//
!    &                  " function!"
        CALL PRSYMB(wunit_q,'-',120,2)
        WRITE(wunit_q,'(2X,a12,2x,4(A16,1X),5x,a16)') 'model',          &
     &        'Edc ','Erdc ','Epc` ', 'Emc ','Emc (Nact - 2e)'
        WRITE(wunit_q,'(3X,a12,4x,5(f16.6,1X))')                        &
     &  'CI+Q energy ',e_tot(i)+Edc_2,e_tot(i)+Erdc_2,e_tot(i)+Epcp_2,  &
     &                 e_tot(i) + Emc_2, Emc_N2 + e_tot(i)
        CALL PRSYMB(wunit_q,'-',120,2)
 
      end do


      deallocate(e_corr)
      deallocate(e_tot)

      deallocate(nash_qcorr)
      deallocate(l_ref)
      deallocate(reference_energy_qcorr)
      deallocate(ref_wavefunction_qcorr)
      deallocate(ref_wavefunction_coeff_qcorr)
      deallocate(ref_wavefunction_coeff)

!     finalize
      call exit_qcorr(relQC)

      call QEXIT('QCORR')

      end subroutine Q_correction
!     ***********************************************************************

      subroutine get_Q_correction(C,LUC,NBL_C,                          &
     &                            NSSOA,                                &
     &                            NSSOB,                                &
     &                            IBLOCK,                               &
     &                            NAEL,NBEL,IASTR,IBSTR,                &
     &                            IBLTP,NSMST,                          &
     &                            IASPGPTP,IBSPGPTP,ICISTR,IPRNT,       &
     &                            ref_wavefunction,                     &
     &                            ref_wavefunction_coeff,nref,nashref,  &
     &                            nash,nfsym,                           &
     &                            nrorb,relQC, wunit_q, lref,cvorb)
!     ***********************************************************************

      use mospinor_info
      implicit none
      real(8), intent(inout)          :: C(*)
      integer, intent(inout)          :: IASTR(NAEL,*), IBSTR(NBEL,*)
      integer, intent(in)             :: lref(nref,2),nashref(nfsym)
      integer, intent(in)             :: luc, nbl_c, nsmst, nfsym
      integer, intent(in)             :: nael, nbel, icistr, iprnt
      integer, intent(in)             :: NSSOA(NSMST,*), NSSOB(NSMST,*)
      integer, intent(in)             :: IBLTP(*), nash(nfsym)
      integer, intent(in)             :: IBLOCK(8,NBL_C)
      integer, intent(in)             :: IASPGPTP, IBSPGPTP
      integer, intent(in)             :: nrorb
      integer, intent(in)             :: nref, wunit_q
      integer, intent(in)             :: cvorb(2)
      character(len=240), intent(in)  :: ref_wavefunction(nref,2)
      logical, intent(in)             :: relQC

!     output
      real(8) , intent(inout)         :: ref_wavefunction_coeff(nref)
     
!     scratch
      character (len=240), allocatable:: ref_wavefunction_tmp(:,:)
      integer,             allocatable:: istr(:,:)
      integer,             allocatable:: ref_order(:,:)
      integer                         :: iel, i, j
      integer                         :: ia, iabas, ib, ibbas
      integer                         :: iampack, imzero, idet, idum
      integer                         :: ilena, ilena2, ilenb, ilenb2
      integer                         :: iasm, iatp, ibsm, ibtp
      integer                         :: iref_tmp, irestr_ana, jblock
      integer                         :: minia, nastr1, nbstr1, nia, nib
      integer                         :: off_a, off_b
      integer                         :: nal, nbt, threshold
      logical                         :: is_not_ref

      allocate(ref_wavefunction_tmp(nref,2))
      ref_wavefunction_tmp(:,1) = ' '
      ref_wavefunction_tmp(:,2) = ' '
      allocate(istr(max(nael,nbel),2))

      if( ICISTR .ge. 2 ) call REWINE(LUC,-1)

      threshold = cvorb(1) + cvorb(2)
      do i = 1, size(nashref)
        threshold = threshold + nashref(i)
      end do

      allocate(ref_order(threshold,2))

      if(iprnt > 100)                                                   &
      print *, 'nfsym, nash(1:nfsym), threshold',                       &
                nfsym, nash(1:nfsym), threshold

      IDET     = 0
      do JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)

!        Obtain alpha/beta strings of sym IASM/IBSM and type IATP/IBTP
         if(relQC)then 
#ifdef PRG_DIRAC
           call GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,          &
     &                                NASTR1,IASTR)
           call GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,          &
     &                                NBSTR1,IBSTR)
#endif
         else

           IDUM = 0
           call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,        &
     &                            nrorb,0,IDUM,IDUM)
           call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,        &
     &                            nrorb,0,IDUM,IDUM)
         end if

         if(IBLTP(IASM).eq.2) then
            IRESTR_ANA = 1
         else
            IRESTR_ANA = 0
         end if
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IMZERO = 0
         if(ICISTR.ge.2) then
!          read in a Type-Type-symmetry block
           call IFRMDS(IDET,1,-1,LUC)
           call FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
           IDET = 0
         end if
         if(IMZERO.ne.1) then
            IBBAS = 1
            IABAS = 1
            do IB = IBBAS,IBBAS+NIB-1                                  ! loop over beta strings
               if(IRESTR_ANA.eq.1.and.IATP.eq.IBTP) then
                  MINIA = IB - IBBAS + IABAS
               else
                  MINIA = IABAS
               end if
               do IA = MINIA,IABAS+NIA-1                              ! loop over alpha strings 
                 IDET = IDET + 1
                   
                is_not_ref = .false.

                !> temporary hack for CAS core orbitals in MRCI, for example orbitals that were not part of the CAS reference
!                  TODO: make a permanent solution for non-rel CI
                if(relQC)then

                   nal   = 0
                   nbt   = 0
                   istr  = 0

                   do j = 1,nael

                     if(iastr(j,ia) > threshold)then
                       is_not_ref = .true.
                       exit
                     end if
                     off_a = imosp_luci2dirac1(ireots (iastr(j,ia)))

                     if(iprnt > 200)then
                       print *,                                         &
                       'imosp_luci2dirac1(ireots (iastr(iel,ia)))',     &
                        imosp_luci2dirac1(ireots (iastr(j,ia))) 
                     end if

                     ! check if actual count in gerade/ungerade
                     if(off_a > nash(1))then ! ungerade, cvorb(2)

!                      check for non-CAS "core" orbital
                       if(off_a-nash(1) > cvorb(2))then
                         nal         = nal + 1
                         istr(nal,1) = off_a-cvorb(2)-nash(1)+nashref(1)
                         if(iprnt>200) print *, 'set istr..',istr(nal,1)
                       end if
                     else ! gerade, cvorb(1)

!                      check for non-CAS "core" orbital
                       if(off_a > cvorb(1))then
                         nal         = nal + 1
                         istr(nal,1) = off_a-cvorb(1)
                         if(iprnt>200) print *, 'set istr..',istr(nal,1)
                       end if
                     end if
                   end do

                   do j = 1,nbel

                     if(ibstr(j,ib) > threshold)then
                       is_not_ref = .true.
                       exit
                     end if

                     off_b = imosp_luci2dirac2(ireots (ibstr(j,ib)))
                     if(iprnt > 200)then
                       print *,                                         &
                       'imosp_luci2dirac2(ireots (ibstr(iel,ib)))',     &
                        imosp_luci2dirac2(ireots (ibstr(j,ib))) 
                     end if

                     ! check if actual count in gerade/ungerade
                     if(off_b > nash(1))then ! ungerade, cvorb(2)

!                      check for non-CAS "core" orbital
                       if(off_b-nash(1) > cvorb(2))then
                         nbt         = nbt + 1
                         istr(nbt,2) = off_b-cvorb(2)-nash(1)+nashref(1)
                         if(iprnt>200) print *, 'set istr..',istr(nbt,2)
                       end if
                     else ! gerade, cvorb(1)

!                      check for non-CAS "core" orbital
                       if(off_b > cvorb(1))then
                         nbt         = nbt + 1
                         istr(nbt,2) = off_b-cvorb(1)
                         if(iprnt>200) print *, 'set istr..',istr(nbt,2)
                       end if
                     end if
                   end do

                   if(iprnt>100)then
                     print *, 'a-string: ',istr(1:nal,1)
                     print *, 'b-string: ',istr(1:nbt,2)
                   end if

                   write(ref_wavefunction_tmp(1,1),"(60i4)")            &
     &             (istr(iel,1),iel=1,nal)
                   write(ref_wavefunction_tmp(1,2),"(60i4)")            &
     &             (istr(iel,2),iel=1,nbt)

                 else ! non-rel CI
                   if((cvorb(1)+cvorb(2)) > 0)                          &
               call quit('no core orbitals yet in non-rel CI+Q! FIXME!')
                   write(ref_wavefunction_tmp(1,1),"(60i4)")            &
     &             (iastr(iel,ia),iel=1,nael)                         
                   write(ref_wavefunction_tmp(1,2),"(60i4)")            &
     &             (ibstr(iel,ib),iel=1,nbel)
                 end if ! relQC switch

                 if(is_not_ref) cycle ! cycle if current determinant is for sure not included in the reference space

                 ilena   = len(trim(ref_wavefunction_tmp(1,1)))      
                 ilenb   = len(trim(ref_wavefunction_tmp(1,2)))      


                 !> sort numbers from 1 ... x such that we do not miss anything in linear symmetry (reordered spinors wrt mj)
                 ref_order = 0
                 do i = 1, ilena/4
                   read(ref_wavefunction_tmp(1,1)(1+(i-1)*4:4*i),'(i4)')&
     &             ref_order(i,1)
                 end do
                 call sort(ref_order(1,1),ilena/4)
                 do i = 1, ilenb/4
                   read(ref_wavefunction_tmp(1,2)(1+(i-1)*4:4*i),'(i4)')&
     &             ref_order(i,2)
                 end do
                 call sort(ref_order(1,2),ilenb/4)

                 !> write sorted strings
                 ref_wavefunction_tmp(1,1)(1:ilena) = ' '
                 ref_wavefunction_tmp(1,2)(1:ilenb) = ' '
                 write(ref_wavefunction_tmp(1,1),"(60i4)")              &
     &           (ref_order(i,1),i=1,ilena/4)
                 write(ref_wavefunction_tmp(1,2),"(60i4)")              &
     &           (ref_order(i,2),i=1,ilenb/4)
                 
                 if(iprnt > 100)then
                   write(wunit_q,*) ' comparison strings     ',         &
     &                                ilena,ilenb
                   write(wunit_q,*)ref_wavefunction_tmp(1,1)(1:ilena)
                   write(wunit_q,*)ref_wavefunction_tmp(1,2)(1:ilenb)
                 end if
                                                                     
                 iref_tmp = 0

                 do 

                   iref_tmp = iref_tmp + 1                           ! exit if iref_tmp > nref
                   if(iref_tmp > nref) exit

                   if(iprnt > 100)then
                     write(wunit_q,*) ' reference strings: for ',       &
     &                                  iref_tmp
                     write(wunit_q,*)                                   &
     &               ref_wavefunction(iref_tmp,1)(1:lref(iref_tmp,1))
                     write(wunit_q,*)                                   &
     &               ref_wavefunction(iref_tmp,2)(1:lref(iref_tmp,2))
                   end if
                                                                     
                   if(ref_wavefunction_tmp(1,1)(1:ilena) ==             &! compare alpha strings
     &                ref_wavefunction(iref_tmp,1)(1:lref(iref_tmp,1))  &
     &             .and.                                                &
     &                ref_wavefunction_tmp(1,2)(1:ilenb) ==             &! compare beta strings
     &                ref_wavefunction(iref_tmp,2)(1:lref(iref_tmp,2))) &
     &             then
!                    print *, 'total match found , coefficient ',c(idet)
                     ref_wavefunction_coeff(iref_tmp) = C(IDET)         ! if true, save coeff
                     exit
                   end if
                 end do
               end do
            end do
         end if
      end do

      deallocate(ref_order)
      deallocate(ref_wavefunction_tmp)
      deallocate(istr)

      end subroutine get_Q_correction
!     ***********************************************************************

      subroutine E_Q(e_corr,ref_wavefunction_coeff,nref,                &
     &               ref_wavefunction_coeff_i,N,                        &
     &               Edc_1, Edrc_1, Epc_1, Emc_1, Epcp_1,Emc_1N,        &
     &               Edc_2, Edrc_2, Epc_2, Emc_2, Epcp_2,Emc_2N,        &
     &               iprnt, wunit_q)
!     ***********************************************************************
!
!
! Davidson Correction
!
! Equations from: Chemical Reviews 2012, 112 page 216-218
!
! Edc  = eq. 34 p. 116
! Edrc = eq. 36 p. 116
! Epc  = eq. 38 p. 116
! Emc  = eq. 41 p. 116
!
! c01  = eq. 43 p. 117
! c02  = eq. 44 p. 117
!
!**********************************************************************
!
      implicit none

      real*8,  intent(in)  :: e_corr
      real*8,  intent(in)  :: ref_wavefunction_coeff(nref)
      real*8,  intent(in)  :: ref_wavefunction_coeff_i(nref)
      real*8,  intent(out) :: Edc_1, Edrc_1, Epc_1, Emc_1, Epcp_1,Emc_1N
      real*8,  intent(out) :: Edc_2, Edrc_2, Epc_2, Emc_2, Epcp_2,Emc_2N
      real*8,  intent(in)  :: N
      integer, intent(in)  :: nref
      integer, intent(in)  :: iprnt, wunit_q

      real*8  :: c01             ! c01    = sum ( cp ) **2
      real*8  :: c02             ! c02    = (sum ( cref * cp )) **2
      real*8  :: c02_2  = 0      ! c02_2  = c02 * c02, entspricht c01
      real*8  :: teta_1 = 0      ! teta_1 = arccos( sqrt(c01) )
      real*8  :: teta_2 = 0      ! teta_2 = arccos( sqrt(c02_2))
      real*8  :: NN     = 0.0d0
      integer :: i      = 0

      c01 = 0.0d0      
      c02 = 0.0d0      

      NN = N - 2.0d0

      do i = 1, nref
      c01 = c01 + ref_wavefunction_coeff(i) ** 2
      end do

      teta_1 = acos(sqrt(c01))

      do i = 1, nref
      c02 = c02 + abs(ref_wavefunction_coeff(i)*                        &
     &                ref_wavefunction_coeff_i(i))
      end do

      c02_2  = c02**2
      teta_2 = acos(sqrt(c02_2))

      if(iprnt > 10)then
        write(wunit_q,'(/a )')  " coefficients entering the q equations"
        write(wunit_q,'( a )')  " -------------------------------------"
        write(wunit_q,*) " c01       :",c01
        write(wunit_q,*) " sqrt (c01):",sqrt(c01)
        write(wunit_q,*) " c02       :",c02
        write(wunit_q,*) " c02_2     :",c02_2
        write(wunit_q,'( a/)')  " -------------------------------------"
      end if
     
!     print *, 'tan(2*teta_1)',tan(2*teta_1)
!     print *, 'alter...     ',(sqrt(1-c01)/sqrt(c01))/(1-((1-c01)/c01))

      Edc_1  = (1 - c01) * e_corr            
      Edrc_1 = ((1-c01) / c01) * e_corr
      Epc_1  = e_corr*(sqrt(N**2 + 2*N * (tan(2*teta_1))**2) - N) /     &
     &         (2*(1/cos(2*teta_1) -1))
!     print *, 'blubb ...'
      Epcp_1  = e_corr*(sqrt(N**2 + 2*N * (tan(2*teta_1))**2) - N) /    &
     &         (2*((1+((1-c01)/c01))/(2*((1-c01)/c01))-1))
!     print *, 'blubb ...',Epcp_1
      Epcp_1  = e_corr*(sqrt(N**2 + 2*N * (tan(2*teta_1))**2) - N) /    &
     &         (2*((1/(2*c01-1))-1))
!     print *, 'blubb 2...',Epcp_1
      Epcp_1  = (1 - (2/N)) * ((1-c01)/c01) * e_corr

      if(N > 3.0d0) Emc_1 = ((N - 2) * (N - 3) / (N * (N - 1))) *       &
     &                      ((1 - c01) / c01) * e_corr  
      if(NN > 3.0d0) Emc_1N = ((NN - 2) * (NN - 3) / (NN * (NN - 1))) * &
     &                      ((1 - c01) / c01) * e_corr  

      Edc_2  = (1 - c02_2) * e_corr
      Edrc_2 = ((1-c02_2) / c02_2) * e_corr
      Epc_2  = (sqrt(N **2 + 2*N * (tan(2*teta_2))**2) - N) /           &
     &         (2*(1/cos(2*teta_2)) - 1) * e_corr
      Epcp_2  = (1 - (2/N)) * ((1-c02_2)/c02_2) * e_corr

      if(N > 3.0d0) Emc_2   = ((N - 2) * (N - 3) / (N * (N - 1))) *     &
     &                       ((1 - c02_2) / c02_2) * e_corr  
      if(NN > 3.0d0) Emc_2N = ((NN - 2) * (NN - 3) / (NN * (NN - 1))) * &
     &                      ((1 - c02_2) / c02_2) * e_corr  
      
      end subroutine E_Q
!     ***********************************************************************

      subroutine put_refvec_for_qcorr(C,NBLOCK,IBLOCK,IBLTP,LUC,        &
     &                                LUSCR,ICISTR,IRC,WORK,relQC,nref, &
     &                                reference_energy)
!     incoming variables:
!                         C      - scratch space for ci vector
!                         NBLOCK - # of ttsblocks
!                         IBLOCK - array containing information about each ttsblock
!                         IBLTP  - block types
!                         LUC    - internal file pointer for CI vector file
!                         LUSCR  - file pointer for scratch vector file
!                         ICISTR - currently set to 2 by default
!                         IRC:   - algebra control value: 1 (real) or 2 (complex)
      use qcorr_cfg
      use qcorr_interface

      implicit none

!     input
      integer, intent(in)   :: IBLOCK(8,NBLOCK),IBLTP(*)
      real(8), intent(inout):: C(*), WORK(*)
      real(8), intent(in)   :: reference_energy
      integer, intent(in)   :: nblock, luc, luscr, irc, icistr
      integer, intent(in)   :: nref
      logical, intent(in)   :: relQC
!     scratch
      integer, allocatable  :: astring(:,:)
      integer, allocatable  :: bstring(:,:)
      real(8), allocatable  :: e_corr(:), e_tot(:)
      integer               :: nstates_ci
      integer               :: nref_e
      real(8)               :: q_corr, edci
      integer               :: iactp, ibctp, nbl_c, nael, nbel
      integer               :: ntest, lblk, iboff, ionem
      integer               :: irir, imk2_r, i

      call QENTER('REFVEC')

!     initialize
      call init_qcorr(relQC,.false.)

!     print flag
      ntest = 000
      ntest = max(ntest,print_qcorr)

      open(1234,file='refvec.luci',status='replace',form='unformatted', &
     &access='sequential',action="readwrite",position='rewind')

      nref_e = 1

      write(1234) nref,nash_q(1:nfsym_q)
      write(1234) nref_e
      write(1234) reference_energy
!     write(1234) nref

      print *, 'printing all reference confs, total: ',nref, 'ref e:',  &
     &          reference_energy,'active orbitals: ',nash_q(1:nfsym_q)

      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call REWINE(LUC,-1)
      call itods(-1,1,-1,LUSCR)

!     loop over real and imaginary parts
      do IRIR = 1,IRC

        if (IRIR.EQ.2) then
!         skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if

        do IMK2_R = 1,NMS2VAL_q,1

          if(relQC)then
            IACTP = IST_FOR_DT_q(1,IMK2_R)
            IBCTP = IST_FOR_DT_q(2,IMK2_R)
          else
            IACTP = 1
            IBCTP = 2
          end if
!  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2_q(IMK2_R-1)
          end if
!
          if(relQC)then 
            NBL_C = NBLK_MS2_q(IMK2_R)
          else
            NBL_C = NBLOCK
          end if

          call REWINE(LUSCR,-1)
          call COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          call ITODS(-1,1,-1,LUSCR)

          NAEL = NELEC_q(IACTP)
          NBEL = NELEC_q(IBCTP)

          allocate(astring(NAEL,MXNSTR_q))
          allocate(bstring(NBEL,MXNSTR_q))
!
          call put_refvec_for_qcorrs(C,LUSCR,NBL_C,                     &
     &                               WORK(KNSTSO_q(IACTP)),             &
     &                               WORK(KNSTSO2_q(IBCTP)),            &
     &                               IBLOCK(1,IBOFF),                   &
     &                               NAEL,NBEL,astring, bstring,        &
     &                               IBLTP,NSMST_q,                     &
     &                               IACTP,IBCTP,ICISTR,ntest,          &
     &                               NORB_q,relQC,wunit_q,nref)
!         NAEL/NBEL may be zero and therefore no allocation of string arrays...
          if(allocated(bstring)) deallocate(bstring)
          if(allocated(astring)) deallocate(astring)
        end do
      end do

      close(1234,status="keep")

!     finalize
      call exit_qcorr(relQC)

      call QEXIT('REFVEC')

      end subroutine put_refvec_for_qcorr
!     ***********************************************************************

      subroutine put_refvec_for_qcorrs(C,LUC,NBL_C,                     &
     &                                 NSSOA,                           &
     &                                 NSSOB,                           &
     &                                 IBLOCK,                          &
     &                                 NAEL,NBEL,IASTR,IBSTR,           &
     &                                 IBLTP,NSMST,                     &
     &                                 IASPGPTP,IBSPGPTP,ICISTR,IPRNT,  &
     &                                 nrorb,relQC,wunit_q,nref)
!     ***********************************************************************

      use mospinor_info
      implicit none

      real(8), intent(inout)          :: C(*)
      integer, intent(inout)          :: IASTR(NAEL,*), IBSTR(NBEL,*)
      integer, intent(in)             :: nref
      integer, intent(in)             :: luc, nbl_c, nsmst
      integer, intent(in)             :: nael, nbel, icistr, iprnt
      integer, intent(in)             :: NSSOA(NSMST,*), NSSOB(NSMST,*)
      integer, intent(in)             :: IBLTP(*)
      integer, intent(in)             :: IBLOCK(8,NBL_C)
      integer, intent(in)             :: IASPGPTP, IBSPGPTP
      integer, intent(in)             :: nrorb
      integer, intent(in)             :: wunit_q
      logical, intent(in)             :: relQC

!     scratch
      character (len=240), allocatable:: ref_wavefunction_tmp(:,:)
      integer                         :: iel
      integer                         :: ia, iabas, ib, ibbas
      integer                         :: iampack, imzero, idet, idum
      integer                         :: ilena, ilena2, ilenb, ilenb2
      integer                         :: iasm, iatp, ibsm, ibtp
      integer                         :: iref_tmp, irestr_ana, jblock
      integer                         :: minia, nastr1, nbstr1, nia, nib

      allocate(ref_wavefunction_tmp(2,2))

      if( ICISTR .ge. 2 ) call REWINE(LUC,-1)
      
      IDET     = 0
      do JBLOCK = 1, NBL_C
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)

!        Obtain alpha/beta strings of sym IASM/IBSM and type IATP/IBTP
         if(relQC)then 
#ifdef PRG_DIRAC
           call GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,          &
     &                                NASTR1,IASTR)
           call GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,          &
     &                                NBSTR1,IBSTR)
#endif
         else

           IDUM = 0
           call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,        &
     &                            nrorb,0,IDUM,IDUM)
           call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,        &
     &                            nrorb,0,IDUM,IDUM)
         end if

         if(IBLTP(IASM).eq.2) then
            IRESTR_ANA = 1
         else
            IRESTR_ANA = 0
         end if
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IMZERO = 0
         if(ICISTR.ge.2) then
!          read in a Type-Type-symmetry block
           call IFRMDS(IDET,1,-1,LUC)
           call FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
           IDET = 0
         end if
         if(IMZERO.ne.1) then
            IBBAS = 1
            IABAS = 1
            do IB = IBBAS,IBBAS+NIB-1                                  ! loop over beta strings
               if(IRESTR_ANA.eq.1.and.IATP.eq.IBTP) then
                  MINIA = IB - IBBAS + IABAS
               else
                  MINIA = IABAS
               end if
               do  IA = MINIA,IABAS+NIA-1                              ! loop over alpha strings 
                    IDET = IDET + 1
                   
                    if(iprnt > 100)then
                      write(wunit_q,*) 'idet, c',idet,c(idet),nael,nbel
                      write(wunit_q,"(60i4)")                           &
     &                (iastr(iel,ia),iel=1,nael)
                      write(wunit_q,"(60i4)")                           &
     &                (ibstr(iel,ib),iel=1,nbel)
                    end if

!                   write(ref_wavefunction_tmp(1,1),"(60i4)")           &
!    &              (iastr(iel,ia),iel=1,nael)                         
!                   write(ref_wavefunction_tmp(1,2),"(60i4)")           &
!    &              (ibstr(iel,ib),iel=1,nbel)
                    write(ref_wavefunction_tmp(1,1),"(60i4)")           &
     &             (imosp_luci2dirac1(ireots(iastr(iel,ia))),iel=1,nael)
                    write(ref_wavefunction_tmp(1,2),"(60i4)")           &
     &             (imosp_luci2dirac2(ireots(ibstr(iel,ib))),iel=1,nbel)

                    ilena    = len(trim(ref_wavefunction_tmp(1,1)))
                    ilenb    = len(trim(ref_wavefunction_tmp(1,2)))
                                                                        
                    write(1234) ilena,ilenb,                            &
     &                          ref_wavefunction_tmp(1,1)(1:ilena),     &
     &                          ref_wavefunction_tmp(1,2)(1:ilenb),     &
                                C(IDET)
!                   print *, ref_wavefunction_tmp(1,1)(1:ilena)
!                   print *, ref_wavefunction_tmp(1,2)(1:ilenb)
!                   print *, C(IDET)
   
               end do
            end do
         end if
      end do

      deallocate(ref_wavefunction_tmp)

      end subroutine put_refvec_for_qcorrs
!     ***********************************************************************

!     !> routines for sorting... found on http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
!     !> slightly modified to make it compile - stknecht may 2014

      INTEGER FUNCTION  FindMinimum(x, Start, finish)
         IMPLICIT  NONE
         INTEGER, INTENT(IN)                   :: Start, finish
         INTEGER, DIMENSION(1:finish), INTENT(IN) :: x
         INTEGER                               :: Minimum
         INTEGER                               :: Location
         INTEGER                               :: i
    
         Minimum  = x(Start)! assume the first is the min
         Location = Start! record its position
         DO i = Start+1, finish ! start with next elements
            IF (x(i) < Minimum) THEN!   if x(i) less than the min?
               Minimum  = x(i)!      Yes, a new minimum found
               Location = i                !      record its position
            END IF
         END DO
         FindMinimum = Location        ! return the position
      END FUNCTION  FindMinimum

      SUBROUTINE  Swap(a, b)
         IMPLICIT  NONE
         INTEGER, INTENT(INOUT) :: a, b
         INTEGER                :: Temp

         Temp = a
         a    = b
         b    = Temp
      END SUBROUTINE  Swap

      SUBROUTINE  Sort(x, is_size)
        IMPLICIT  NONE
        INTEGER, INTENT(IN)                       :: is_size
        INTEGER, DIMENSION(1:is_size), INTENT(INOUT) :: x
        INTEGER                                   :: i
        INTEGER                                   :: Location
        INTEGER                                   :: FindMinimum

        DO i = 1, is_size-1! except for the last
           Location = FindMinimum(x, i, is_size)! find min from this to last
           CALL  Swap(x(i), x(Location))! swap this and the minimum
        END DO
      END SUBROUTINE  Sort
