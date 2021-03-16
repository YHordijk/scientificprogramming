!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_TRAMO(desrep,ioch,eps)

      use qstack
 
!     IMPLICIT INTEGER (A-Z)
      IMPLICIT none
!
!---------------Description--------------------------------------------
!
!     Calculates transition moments from the transformed
!     x,y,z dipole moment integrals and the Davidson eigenvectors
!
!     Important: Also for real two-electron integrals (high symmetry and/or 
!     NONREL,SPINFREE Hamiltonians) the property integrals are in complex form.
!     These cases will be handled by the corresponding phase factors stored at the
!     beginning of the property integral file.
!
!---------------Calling variables--------------------------------------
!
      INTEGER                          :: desrep,ioch
      REAL*8,dimension(:)              :: eps
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/complex.inc"
#include "../relccsd/symm.inc"
#include  "polprp_stacklines.h"
!
!---------------Local variables--------------------------------------
!
! a whole bunch of required integer variables for
! counting and offsetting

      INTEGER               ::  i, k, m, n
      INTEGER               ::  arep, brep, crep, abirep, bijrep
      INTEGER               ::  idaoo, idaov, idavo, idavv
      INTEGER               ::  ioff1, ioff2, ioff3
      INTEGER               ::  irep, irep_b, jrep, krep
      INTEGER               ::  off, off1, off2, off3

! filename of property integrals

      Character*8,dimension(3)               :: trafile

! constants

      Complex*16                             :: A0,A1,AI

! index arrays for the nontotally symmetric (NTS) property arrays

      INTEGER,allocatable,dimension(:)       :: JPOOXX,JPOVXX
      INTEGER,allocatable,dimension(:)       :: JPVOXX,JPVVXX

! offset and length arrays for NTS two-electron contributions
                                            
      INTEGER,allocatable,dimension(:)       :: LVVONO
      INTEGER,allocatable,dimension(:)       :: LVNVOO

! phase information for the property x,y,z arrays in case of real arithmetic

      Complex*16,allocatable,dimension(:)    :: pphase
      Complex*16                             :: pphasex
      Complex*16                             :: pphasey
      Complex*16                             :: pphasez

! intermediate and final F_ai, F_aibj matrices (as vectors)

      REAL*8,allocatable,dimension(:)        :: faipx,faipy,faipz
      Complex*16,allocatable,dimension(:)    :: faitx,faity,faitz
      Complex*16,allocatable,dimension(:)    :: fabijtx,fabijty,fabijtz

! dipole moment integrals in MO basis for all three components

      REAL*8,allocatable,dimension(:)        :: doox,dooy,dooz
      REAL*8,allocatable,dimension(:)        :: dovx,dovy,dovz
      REAL*8,allocatable,dimension(:)        :: dvox,dvoy,dvoz
      REAL*8,allocatable,dimension(:)        :: dvvx,dvvy,dvvz

! OV dipole moment integrals holding the generic order OV property ints

      REAL*8,allocatable,dimension(:)        :: dovxg,dovyg,dovzg
                                            
! auxiliary arrays for contractions

      REAL*8,allocatable,dimension(:)        :: aac,aki,aak,abj
                                            
! integral storage arrays

      REAL*8,allocatable,dimension(:)        :: vvoo1
      REAL*8,allocatable,dimension(:)        :: vvoo2
      REAL*8,allocatable,dimension(:)        :: vvoo3
      REAL*8,allocatable,dimension(:)        :: vvoo4
      REAL*8,allocatable,dimension(:)        :: vvoop
                                            
      REAL*8,allocatable,dimension(:)        :: vooo1,vooo2
      REAL*8,allocatable,dimension(:)        :: vovo1,vovo2,vovo3
      REAL*8,allocatable,dimension(:)        :: vovv1,vovv2
                                            
      REAL*8,allocatable,dimension(:)        :: oooo

      REAL*8,allocatable,dimension(:)        :: vabijx
      REAL*8,allocatable,dimension(:)        :: vabijy
      REAL*8,allocatable,dimension(:)        :: vabijz
                                            
! energy denominator array

      REAL*8,allocatable,dimension(:)        :: evvoo

! some parameters and bookkeeping variables

      Logical                                :: LFA = .false.
      Logical                                :: LTR = .true.
      Logical                                :: istaken = .false.

      character*9                            :: preconfn='VVVVXVVOO'
      character*40                           :: iorder
      integer                                :: incsat
      real*8                                 :: x1,x2,x3,x4,xcntrl
      Integer                                :: lp1,lp2
      integer                                :: nabixj,naxbij,nabxij

      Integer                                :: n1,n2,n3,n4,n1e
      Integer                                :: ialloc,iat
      
      Integer*8                              :: n8

      Integer                                :: qst_prpline
      Integer                                :: qst_auxline
!
!---------------interface region --------------------------------------
!
      interface

        SUBROUTINE propfetch(i1,ca1,
     &                          i2,i3,i4,i5,
     &                          ra1,ra2,ra3,
     &                          ra4,ra5,ra6,
     &                          ra7,ra8,ra9,
     &                          ra10,ra11,ra12)
          integer                      :: i1
          complex*16, dimension(:)     :: ca1
          integer                      :: i2,i3,i4,i5
          real*8, dimension(:)         :: ra1,ra2,ra3
          real*8, dimension(:)         :: ra4,ra5,ra6
          real*8, dimension(:)         :: ra7,ra8,ra9
          real*8, dimension(:)         :: ra10,ra11,ra12
        END SUBROUTINE

        INTEGER FUNCTION XCOLDET(i1)
          integer                      :: i1
        END FUNCTION

        SUBROUTINE denomvots(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE denomvovo(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE EARR_VVOO(ra1,ra2)
          real*8, dimension(:)         :: ra1,ra2
        END SUBROUTINE

        SUBROUTINE MATXMAT(i1,ra1,ra2,ra3)
          integer                      :: i1
          real*8, dimension(:)         :: ra1,ra2,ra3
        END SUBROUTINE

        SUBROUTINE rcaxpy(i1,l1,c1,ra1,ca1)
          integer                      :: i1
          logical                      :: l1
          complex*16                   :: c1
          real*8, dimension(:)         :: ra1
          complex*16, dimension(:)     :: ca1
        END SUBROUTINE

        SUBROUTINE rca_negate(i1,ra1)
          integer                      :: i1
          real*8, dimension(:)         :: ra1
        END SUBROUTINE

        SUBROUTINE TM_GEN(i1,i2,i3,ch1,
     &                    l1,
     &                    ca1,ca2,ca3,
     &                    ca4,ca5,ca6)
          integer                       :: i1,i2,i3
          character*30                  :: ch1
          logical                       :: l1
          complex*16, dimension(:)      :: ca1,ca2,ca3
          complex*16, dimension(:)      :: ca4,ca5,ca6
        END SUBROUTINE

        SUBROUTINE mpi_master_vovv_complete(ra1,i1)
          real*8, dimension(:)            :: ra1
          integer*8                       :: i1
        END SUBROUTINE mpi_master_vovv_complete

        Function zdotc (i1,ca1,i2,ca2,i3)
         complex*16                        :: zdotc
         Integer                           :: i1,i2,i3
         Real*8, dimension(*)              :: ca1,ca2
        End Function

      end interface
!
!---------------Executable code--------------------------------------
!
! say hello
!
      if(carith) then
        call pst('Transition moments (complex run)+')
      else
        call pst('Transition moments (real run)+')
      endif

!****************************************
! important remark:
!         Note: Even in SPINFREE/NONREL cases the TMs'
!         can be complex due to the quaternionic structure'
!         of the wave function.'
!****************************************
!     
!  assign stack line numbers in order to avoid direct integers
!  as actual parameters for the qstack calls.
!
      qst_prpline = PRP_STACKLINE
      qst_auxline = AUX_STACKLINE
!
! clear stack line for pushing intermediates
! if stack line is cleared, no action is taken.
! return value i is irrelevant.
!
      i = qstack_drop(qst_prpline)
!
! initialization of constant variables
!
      trafile = (/'XPROPINT','YPROPINT','ZPROPINT'/)
      do i=1,30
        iorder(i:i)=' '
      enddo
      A0=(0.0d0,0.0d0)
      A1=(1.0d0,0.0d0)
      AI=(0.0d0,1.0d0)
!
! Generate index arrays for OO,OV,VO,VV property arrays
!
      n1 = nrep + 1
      allocate(JPOOXX(n1))
      allocate(JPVOXX(n1))
      allocate(JPOVXX(n1))
      allocate(JPVVXX(n1))

      JPOOXX=0
      JPVOXX=0
      JPOVXX=0
      JPVVXX=0

      DO IREP = 1, NREP
        JPOOXX(IREP+1) = JPOOXX(IREP) + MOO(IREP)
        JPOVXX(IREP+1) = JPOVXX(IREP) + MOV(IREP)
        JPVOXX(IREP+1) = JPVOXX(IREP) + MVO(IREP)
        JPVVXX(IREP+1) = JPVVXX(IREP) + MVV(IREP)
      ENDDO

!*****************************************
! complex array for the phases.
!*****************************************

      allocate(pphase(3))
      pphase = A0
!
!***************************************************
! ****  Transition element vectors. Attention
! ****  Even in a real calculation they will be kept
! ****  complex!
! ****  the final transition moment arrays are always complex.
! ****  the intermediate transition moment arrays follow CARITH.
!***************************************************
!
      n1 = mvo(desrep)
      lp1 = n1*rcw   ! actual number of reals in faip(x,y,z) array (r/c)

      allocate(faitx(n1))
      allocate(faity(n1))
      allocate(faitz(n1))
      faitx = A0
      faity = A0
      faitz = A0
      
      allocate(faipx(rcw*n1))
      allocate(faipy(rcw*n1))
      allocate(faipz(rcw*n1))
      faipx = 0.0d0
      faipy = 0.0d0
      faipz = 0.0d0
!
!***************************************************
! ****  Arrays holding the property integrals. 
! ****  They are either
! ****  complex if CARITH=T or real, if CARITH=F.
! ****  However, in the latter case a phase information
! ****  indicating purely imaginary property integrals
! ****  has to be additionally stored and processed.
!***************************************************
!
      idaoo = JPOOXX(NREP+1)   !  dimensions
      idaov = JPOVXX(NREP+1)
      idavo = JPVOXX(NREP+1)
      idavv = JPVVXX(NREP+1)

      allocate(doox(idaoo*rcw),dooy(idaoo*rcw),dooz(idaoo*rcw))
      allocate(dovx(idaov*rcw), dovy(idaov*rcw), dovz(idaov*rcw))
      allocate(dovxg(idaov*rcw),dovyg(idaov*rcw),dovzg(idaov*rcw))
      allocate(dvox(idavo*rcw),dvoy(idavo*rcw),dvoz(idavo*rcw))
      allocate(dvvx(idavv*rcw),dvvy(idavv*rcw),dvvz(idavv*rcw))

      doox  = 0.0d0
      dooy  = 0.0d0
      dooz  = 0.0d0
      dovx  = 0.0d0
      dovy  = 0.0d0
      dovz  = 0.0d0
      dvox  = 0.0d0
      dvoy  = 0.0d0
      dvoz  = 0.0d0
      dvvx  = 0.0d0
      dvvy  = 0.0d0
      dvvz  = 0.0d0
      dovxg = 0.0d0
      dovyg = 0.0d0
      dovzg = 0.0d0

!*********************************************************************
!  VVOO-type integral arrays, they always remain allocated since they
!  hold integrals repeatedly needed in nearly every step.
!*********************************************************************

      n1 = ivvoott(nrep+1)
      allocate(vvoo1(n1*rcw))
      allocate(vvoo2(n1*rcw))
      vvoo1 = 0.0d0
      vvoo2 = 0.0d0

      n1 = max0(lvvoot(nrep+1),jvovo(nrep+1),
     &          ivvoott(nrep+1),ivvoot(nrep+1),
     &          JVOOV(nrep+1))

      allocate(vvoo3(n1*rcw))
      allocate(vvoo4(n1*rcw))
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
!
!  auxiliary arrays
!
      n1 = MVV(1)
      allocate(aac(n1*rcw))
      n1 = MOO(1)
      allocate(aki(n1*rcw))
      n1 = MVO(1)
      allocate(aak(n1*rcw))
      n1 = JPVOXX(NREP+1) 
      allocate(abj(n1*rcw))
      aac = 0.0d0
      aki = 0.0d0
      aak = 0.0d0
      abj = 0.0d0
!______________________________________________________
!|
!|
!|
!|
!|
!|       START CALCULATIONS
!|
!|
!|
!|
!|_____________________________________________________
!
!  Initialize VVOO integral arrays
!  v_abij ---> vvoo1  /  v*_cbjk ---> vvoo2
!
      n1 = ivvoott(nrep+1)
      call getvvoo(vvoo1)
      call xcopy(n1,vvoo1,1,vvoo2,1)
      if(carith) call conjuga(n1,vvoo2,1)
      call denomvvoo(eps,vvoo1,vvoo1)
      call denomvvoo(eps,vvoo2,vvoo2)
!
!  get the property integrals and phases from corresponding files.
!
      CALL propfetch(ioch,pphase,    !fetch properties, we transfer dimensions!
     &               idaoo,idaov,idavo,idavv,
     &               doox,dooy,dooz,
     &               dovx,dovy,dovz,
     &               dvox,dvoy,dvoz,
     &               dvvx,dvvy,dvvz)
      write(*,*) 'Phase values from propfetch:'
      write(*,*) 'Phasex:',pphase(1)
      write(*,*) 'Phasey:',pphase(2)
      write(*,*) 'Phasez:',pphase(3)

!     if(carith) then
!       write(*,*) 'control sumc 1:',zdotc(idaoo,doox,1,doox,1)
!       write(*,*) 'control sumc 1:',zdotc(idaoo,dooy,1,dooy,1)
!       write(*,*) 'control sumc 1:',zdotc(idaoo,dooz,1,dooz,1)
!       write(*,*) 'control sumc 2:',zdotc(idaov,dovx,1,dovx,1)
!       write(*,*) 'control sumc 2:',zdotc(idaov,dovy,1,dovy,1)
!       write(*,*) 'control sumc 2:',zdotc(idaov,dovz,1,dovz,1)
!       write(*,*) 'control sumc 3:',zdotc(idavo,dvox,1,dvox,1)
!       write(*,*) 'control sumc 3:',zdotc(idavo,dvoy,1,dvoy,1)
!       write(*,*) 'control sumc 3:',zdotc(idavo,dvoz,1,dvoz,1)
!       write(*,*) 'control sumc 4:',zdotc(idavv,dvvx,1,dvvx,1)
!       write(*,*) 'control sumc 4:',zdotc(idavv,dvvy,1,dvvy,1)
!       write(*,*) 'control sumc 4:',zdotc(idavv,dvvz,1,dvvz,1)
!     else
!       write(*,*) 'control sumr 1:',dot_product(doox,doox)
!       write(*,*) 'control sumr 1:',dot_product(dooy,dooy)
!       write(*,*) 'control sumr 1:',dot_product(dooz,dooz)
!       write(*,*) 'control sumr 2:',dot_product(dovx,dovx)
!       write(*,*) 'control sumr 2:',dot_product(dovy,dovy)
!       write(*,*) 'control sumr 2:',dot_product(dovz,dovz)
!       write(*,*) 'control sumr 3:',dot_product(dvox,dvox)
!       write(*,*) 'control sumr 3:',dot_product(dvoy,dvoy)
!       write(*,*) 'control sumr 3:',dot_product(dvoz,dvoz)
!       write(*,*) 'control sumr 4:',dot_product(dvvx,dvvx)
!       write(*,*) 'control sumr 4:',dot_product(dvvy,dvvy)
!       write(*,*) 'control sumr 4:',dot_product(dvvz,dvvz)
!     endif

! canary values for catching errors. For complex arithmetic these
! phase variables must not be accessed or used.

      pphasex = -999999.0d0  
      pphasey = -999999.0d0
      pphasez = -999999.0d0
      
      If(carith.eqv..false.) then
        if(aimag(pphase(1)).ne.0) then
          write(*,*) 'Real mode, but adjusting dvox'
          pphasex = AI
        else
!         write(*,*) 'Unchanged dvox'
          pphasex = A1
        endif
        if(aimag(pphase(2)).ne.0) then
          write(*,*) 'Real mode, but adjusting dvoy'
          pphasey = AI
        else
!         write(*,*) 'Unchanged dvoy'
          pphasey = A1
        endif
        if(aimag(pphase(3)).ne.0) then
          write(*,*) 'Real mode, but adjusting dvoz'
          pphasez = AI
        else
!         write(*,*) 'Unchanged dvoz'
          pphasez = A1
        endif
      Endif   !carith
!
!  storage of property contributions is F_ai. Therefore the
!  First contribution is D_ai and not D_ia
!
!------------------------
!------  F_ai(0)
!------------------------
!
      n1 = MVO(desrep)       !number of integrals   (OV = VO!)
      off = JPVOXX(desrep)*rcw + 1    !offset for prop. integrals (OV = VO!)

! even if the dimensions are equivalent we copy the dov* integrals since
! they bear the correct phases in case of imaginary property integrals.

      call xcopy(n1,dvox(off),1,faipx,1)
      call xcopy(n1,dvoy(off),1,faipy,1)
      call xcopy(n1,dvoz(off),1,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th order                     '
      incsat = 0
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 
!
!------------------------
!------  F_ai(1)
!------------------------
!
! sort v_ab,ij --> v_ai,bj
! the minus sign in the contraction is accounted for in the 
! XGEMV routines!
! We contract over D_jb = D*_bj !!
!
      vvoo3 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MVO(desrep)
      CALL XGEMV ('N',M,K,-A1,vvoo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,-A1,vvoo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,-A1,vvoo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)
!
! add result on top of the final arrays
!
      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
!
!  calculate corresponding 0th + 1st order TMs
!
      iorder='0th + 1st order               '
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!
!  enter calculation for second-order contributions
!
!------------------------
!------  F_ai(2)A
!------------------------
!
!  sort  v_ab,jk  -->  v_a,bjk    (vvoo3)
!  sort v*_cb,jk  -->  v*_c,bjk   (vvoo4)
!         v_a,bjk  *  v*_c,bjk   -->  aux1
!
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      call srt1t2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo2,vvoo4)
      call cntrct ('N','T',nv,nv,nvoot,a1,vvoo3,vvoo4,A0,aac,nrep)
!
!  contract A_ac and D_ci(x,y,z) over c in special loop because
!  D_ci is not totally symmetric and the order of Y_c in D_ci does not
!  follow the one in A_ac.
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         crep = multb(desrep+nrep,irep,2)
         arep = crep
         if(MULTB(arep,irep,2).ne.desrep) 
     &         stop 'Symmetry problem in F_ai(2)A'
         m = NV(arep)
         n = NO(irep)
         k = NV(crep)
         ioff1 = JJVV(arep,crep)*rcw + 1
         ioff2 = (JPVOXX(desrep) + JJVO(crep,irep))*rcw + 1
         ioff3 = JJVO(arep,irep)*rcw + 1
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvoy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aac(ioff1),m,dvoz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo
      n1 = MVO(desrep)


      CALL XSCAL(n1,-A1/2.0d0,faipx,1)
      CALL XSCAL(n1,-A1/2.0d0,faipy,1)
      CALL XSCAL(n1,-A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

!     write(*,*) '------  F_ai(2)A finished'

      iorder='0th + 1st + 2ndA'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!
!------------------------
!------  F_ai(2)B
!------------------------
!
!  sort  v*_bc,jk (vvoo2)  -->  v*_bcj,k  (vvoo3)  (1T3)
!  sort  v_bc,ji  (vvoo1)  -->  v_bcj,i   (vvoo4)  (1T3)
!
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo2,vvoo3)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo1,vvoo4)
      CALL CNTRCT('T','N',NO,NO,NVVOT,A1,vvoo3,vvoo4,A0,aki,NREP)
!
!  contract D_ak(x,y,z) and A_ki over k in special loop because
!  D_ak is not totally symmetric!
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         krep = irep
         if(MULTB(arep,irep,2).ne.desrep) 
     &         stop 'Symmetry problem in F_ai(2)B'
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = (JPVOXX(desrep) + JJVO(arep,krep))*rcw + 1
         ioff2 = JJOO(krep,irep)*rcw + 1
         ioff3 = JJVO(arep,irep)*rcw + 1
         call xgemm ('N','N',m,n,k,a1,dvox(ioff1),m,aki(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvoy(ioff1),m,aki(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvoz(ioff1),m,aki(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      CALL XSCAL(n1,-A1/2.0d0,faipx,1)
      CALL XSCAL(n1,-A1/2.0d0,faipy,1)
      CALL XSCAL(n1,-A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

!     write(*,*) '------  F_ai(2)B finished'

      iorder='0th + 1st + 2ndAB'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!
!------------------------
!------  F_ai(2)C
!------------------------
!
!
!  sort  v_ab,ij  (vvoo1)  -->  v_ai,bj    (vvoo3)  (1TT4)
!  v*_cb,kj = v*_bc,jk !
!  sort  v*_bc,jk (vvoo2)  -->  v*_bj,ck   (vvoo4)  (1TT4)

!  additionally, allocate proper array
!
      vvoo3 = 0.0d0
      vvoo4 = 0.0d0
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo1,vvoo3)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &                 MVO,JVOVO,JJVO,JJVO,vvoo2,vvoo4)
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      vovo1 = 0.0d0
      if(ialloc.ne.0)  call quit('vovo1 allocation problem (tramo)')
!
!  do contraction of v_ai,bj  *  v*_bj,ck  =  A_ai,ck
!
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,vvoo4(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  do contraction of A_ai,ck  *  D_ck (determine proper symmetry!)
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      M = MVO(desrep)
      K = MVO(desrep)
      OFF1 = jvovo(desrep)*rcw + 1
      OFF2 = JPVOXX(desrep)*rcw + 1
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvox(OFF2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvoy(OFF2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vovo1(OFF1),M,dvoz(OFF2),1,A0,faipz,1)
      deallocate(vovo1)
!
!  scale result by +0.5d0 and add to final values
!
      n1 = MVO(desrep)
      CALL XSCAL(n1,A1/2.0d0,faipx,1)
      CALL XSCAL(n1,A1/2.0d0,faipy,1)
      CALL XSCAL(n1,A1/2.0d0,faipz,1)

      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)
      
!     write(*,*) '------  F_ai(2)C finished'

      iorder='0th + 1st + 2ndABC'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 
!
!
!--------------------------------------
!------  F_ai(2,1), a DOO contribution
!--------------------------------------
!
!  In parallel case we need to collect all VOVV integrals!
!  only temporarily allocate space for the VOVV integrals.
!
      iat = 0
      n1 = ivovvt(nrep+1)
      allocate(vovv1(n1*rcw),stat=ialloc); iat = iat + ialloc
      n1 = lvovvt(nrep+1)
      allocate(vovv2(n1*rcw),stat=ialloc); iat = iat + ialloc
      if(ialloc.ne.0) STOP 'VOVV allocation problem in TRAMO!'
!
!  sort V_aj,bc  --->  Va,jbc   (1S2)    vovv1 --> vovv2
!  sort v_bc,kj  --->  v_jbc,k  (SRT22)  vvoo1 --> vvoo3
!  then contract Va,jbc  *  v_jbc,k  ---> A_ak (tot symm)
!  divide A_ak by (eps_a  -  eps_k), and multiply by -1.0
!  and contract A_ak with D_ki(x,y,z) (special loop!)
!
#if defined (VAR_MPI)
      call rdvovv(vovv1)
      n8 = ivovvt(nrep+1)
      call mpi_master_vovv_complete(vovv1,n8*rcw)  !fetch integrals from slaves
#else
      CALL RDVOVV(vovv1); 
#endif
      vovv2 = 0.0d0
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             vovv1,vovv2)
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            vvoo1,vvoo3)
      CALL CNTRCT ('N','N',NV,NO,NOVVT,A1,vovv2,vvoo3,A0,aak,NREP)
      CALL DENOMVOTS (EPS,aak)
      CALL XSCAL(MVO(1),-A1,aak,1)
      deallocate(vovv1,vovv2)
!
!  contract A_ak/(e_a - e_k) and D_ki(x,y,z) over k in special loop because
!  D_ki is not totally symmetric and the order of Y_k in D_ki does not
!  follow the one in A_ak.
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         krep = multb(desrep+nrep,irep,2)
         arep = krep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JJVO(arep,krep)*rcw + 1
         ioff2 = (JPOOXX(desrep) + JJOO(krep,irep))*rcw + 1
         ioff3 = JJVO(arep,irep)*rcw + 1
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,doox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,1'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!     write(*,*) '------  F_ai(2,1) finished'
      
!
!--------------------------------------
!------  F_ai(2,2), a DOO contribution
!--------------------------------------
!
      iat = 0
      n1 = ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc); iat = iat + ialloc
      n1 = lovoot(nrep+1)
      allocate(vooo2(n1*rcw),stat=ialloc); iat = iat + ialloc
      if(ialloc.ne.0) STOP 'VOOO allocation problem in TRAMO!'
!
!  sort v_ab,jk  --->  v_a,bjk   (still in vvoo1 ---> vvoo3)
!  sort V*_bl,jk --->  V*_bjk,l  vooo1 ---> vooo2
!  contract  v_a,bjk  *  V*_bjk,l  ---> A_al (aak)
!  divide A_ak by (eps_a  -  eps_k), and multiply by -1.0
!  and contract A_ak with D_ki(x,y,z)
!
      call SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      CALL GETVOOO(vooo1);vooo2 = 0.0d0
      IF(CARITH) CALL CONJUGA(ivooot(nrep+1),vooo1,1)
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo1,vooo2)
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,vvoo3,vooo2,A0,aak,NREP)
      CALL DENOMVOTS (EPS,aak)
      CALL XSCAL(MVO(1),-A1,aak,1)
      deallocate(vooo1,vooo2)
!
!  special loop as above for contracting A_al and D_li (identical to
!  A_ak and D_ki)
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      do irep = 1, nrep
         krep = multb(desrep+nrep,irep,2)
         arep = krep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NO(krep)
         ioff1 = JJVO(arep,krep)*rcw + 1
         ioff2 = (JPOOXX(desrep) + JJOO(krep,irep))*rcw + 1
         ioff3 = JJVO(arep,irep)*rcw + 1
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,doox(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooy(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,aak(ioff1),m,dooz(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,2'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 
      
!     write(*,*) '------  F_ai(2,2) finished'
!
!-------------------------------------
!------  F_ai(2,3), a DVV contribution
!-------------------------------------
      n1 = ivooot(nrep+1)
      allocate(vooo1(n1*rcw),stat=ialloc)
      n1 = lovoot(nrep+1)
      allocate(vooo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOOO allocation problem in TRAMO!'
!
!  sort v_bc,jk  --->  v_b,cjk   (still in vvoo1 ---> vvoo3)
!  sort V*_ci,jk --->  V*_cjk,i  vooo1 ---> vooo2
!  contract  v_b,cjk  *  V*_cjk,i  ---> A_bi (aak)
!  divide A_bi by (eps_b  -  eps_i)
!  and contract D_ab(x,y,z) with A_bi
!
      call SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
      CALL GETVOOO(vooo1)
      IF(CARITH) CALL CONJUGA(ivooot(nrep+1),vooo1,1)
      vooo2 = 0.0d0
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           vooo1,vooo2)
      CALL CNTRCT ('N','N',NV,NO,NVOOT,A1,vvoo3,vooo2,A0,aak,NREP)
      deallocate(vooo1,vooo2)
      CALL DENOMVOTS (EPS,aak) !A_bi now in aak.
!
! do special contraction of D_ab and A_bi over b
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         brep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(brep)
         ioff1 = (JPVVXX(desrep) + JJVV(arep,brep))*rcw + 1
         ioff2 = JJVO(brep,irep)*rcw + 1
         ioff3 = JJVO(arep,irep)*rcw + 1
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,aak(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,aak(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,aak(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,3'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!     write(*,*) '------  F_ai(2,3) finished'
!
!-------------------------------------
!------  F_ai(2,4), a DVV contribution
!-------------------------------------
!
!  only temporarily allocate space for the VOVV integrals.
!
      n1 = ivovvt(nrep+1)
      allocate(vovv1(n1*rcw),stat=ialloc)
      n1 = lvovvt(nrep+1)
      allocate(vovv2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0) STOP 'VOVV allocation problem in TRAMO!'
!
!  sort v_bc,ij  --->  v_jbc,i  (SRT22)  vvoo1 --> vvoo3
!  sort V_dj,bc  --->  Vd,jbc   (1S2)    vovv1 --> vovv2
!  then contract Vd,jbc  *  v_jbc,i  ---> A_di (tot symm)
!  divide A_di by (eps_d  -  eps_i)
!  and contract D_ad(x,y,z) with Adi
!
      CALL SRT22 (NREP,MULTB,LFA,NVVT,NO,NO,NOVVT,LOOVVT,LLOVVT,
     &            vvoo1,vvoo3)
#if defined (VAR_MPI)
      call rdvovv(vovv1)
      n8 = ivovvt(nrep+1)
      call mpi_master_vovv_complete(vovv1,n8*rcw)  !fetch integrals from slaves
#else
      CALL RDVOVV(vovv1)
#endif
      vovv2 = 0.0d0
      CALL SRT1S2 (NREP,MULTB,LFA,NVO,NV,NO,NVVT,LVOVVT,LLOVVT,
     &             vovv1,vovv2)
      CALL CNTRCT ('N','N',NV,NO,NOVVT,A1,vovv2,vvoo3,A0,aak,NREP)
      deallocate(vovv1, vovv2)
      CALL DENOMVOTS (EPS,aak)
!
! do special contraction of D_ad and A_di over d! (loop from above)
! is the same as of D_ab and A_bi over b
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0

      do irep = 1, nrep
         arep = multb(desrep+nrep,irep,2)
         brep = irep
         if(MULTB(arep,irep,2).ne.desrep) stop 'munch:-('
         m = NV(arep)
         n = NO(irep)
         k = NV(brep)
         ioff1 = (JPVVXX(desrep) + JJVV(arep,brep))*rcw + 1
         ioff2 = JJVO(brep,irep)*rcw + 1
         ioff3 = JJVO(arep,irep)*rcw + 1
         call xgemm ('N','N',m,n,k,a1,dvvx(ioff1),m,aak(ioff2),k,
     &               a0,faipx(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvy(ioff1),m,aak(ioff2),k,
     &               a0,faipy(ioff3),m)
         call xgemm ('N','N',m,n,k,a1,dvvz(ioff1),m,aak(ioff2),k,
     &               a0,faipz(ioff3),m)
      enddo

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,4'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!     write(*,*) '------  F_ai(2,4) finished'
!
!-------------------------------------
!------  F_ai(2,5), a DOV contribution
!-------------------------------------
!
!  read Vaj,bi ---> Vai,bj (SRT16) (vovo1 ---> vovo2)
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo2(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo3(n1*rcw),stat=ialloc)
        if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      CALL getvovo(vovo1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
!
!  v_b>c,k>j  --->  -v_b>c,j>k  (RCA_NEGATE) vvoo1 --> -vvoo1
!  -v_b>c,j>k  --->  -v_bj,ck   (SRT1TT4) -vvoo1 --> vovo1
!  negate vvoo1 again in order to restore original values.

      n1 = ivvoott(nrep+1) 
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo1)
      CALL RCA_NEGATE(n1,vvoo1)

!  contract Vai,bj (vovo2) and -v_bj,ck (vovo1) to Aai,ck (vovo3)
!  Aai,ck  attains  <||>,<||>  BKC  structure

      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = JVOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo2(OFF1),M,vovo1(OFF2),K,
     &               A0,vovo3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO

! ************** ATT !  no backsort, instead denomvovo!
! divide A_ai,ck array (JVOVO structure) by the corresponding
! energy denominators e_a - e_i + e_c - e_k

      CALL DENOMVOVO(eps,vovo3)

!  contract A_ai,kc (vovo3)  and Dkc over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MVO(desrep)

      if(carith) call conjuga(M,dvox(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      if(carith) call conjuga(M,dvox(off2),1)

      if(carith) call conjuga(M,dvoy(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      if(carith) call conjuga(M,dvoy(off2),1)

      if(carith) call conjuga(M,dvoz(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)
      if(carith) call conjuga(M,dvoz(off2),1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      deallocate(vovo1,vovo2,vovo3)

      iorder='0th + 1st + 2ndABC + 2,5'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!     write(*,*) '------  F_ai(2,5) finished'
!
!-------------------------------------
!------  F_ai(2,6): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  sort the -v_ab,ji (vvoo1) ---> -v_ai,bj (vovo3)
!  and the V*_bk,cj ---> V*_bj,ck (SRT 16)
!  then contract -v_ai,bj and V*_bj,ck to  A_ai,ck
!       VO,VO  x  VO,VO  --->  VO,VO <||>,<||>
!
      n1 = ivvoott(nrep+1) 
      CALL RCA_NEGATE(n1,vvoo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo3)
      CALL RCA_NEGATE(n1,vvoo1)

      CALL getvovo(vovo1)
      n1 = ivovo(nrep+1) 
      if(carith) call conjuga(n1,vovo1,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = J2VOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo3(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo1(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  A_ai,ck now in vovo1
!  divide by energy denominators
!
      CALL DENOMVOVO(eps,vovo1)

!  contract B_ai,ck (vovo1) and Dkc = D*ck over corresponding symmetry
!  and add to total result.

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MVO(desrep)

      if(carith) call conjuga(M,dvox(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      if(carith) call conjuga(M,dvox(off2),1)

      if(carith) call conjuga(M,dvoy(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      if(carith) call conjuga(M,dvoy(off2),1)

      if(carith) call conjuga(M,dvoz(off2),1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovz(off2),1,A0,faipz,1)
      if(carith) call conjuga(M,dvoz(off2),1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      deallocate(vovo1,vovo2,vovo3)

      iorder='0th + 1st + 2ndABC + 2,6'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

!     write(*,*) '------  F_ai(2,6) finished'
!
!-------------------------------------
!------  F_ai(2,7): a DOV contribution
!-------------------------------------
!
!
! For the final main block contributions we need the OV
! property integrals in generic order and not in the VO order!
! we fetch them from the system stack (pushed by propertyhandler)
! in a FIFO manner! Since the push was done only once we have to
! repush the data back one stack
!
      if(qstack_popf(qst_auxline,dovxg).ne.idaov*rcw) stop 'Qpferr'
      if(qstack_popf(qst_auxline,dovyg).ne.idaov*rcw) stop 'Qpferr'
      if(qstack_popf(qst_auxline,dovzg).ne.idaov*rcw) stop 'Qpferr'
      if(qstack_push(qst_auxline,idaov*rcw,dovxg).ne.idaov*rcw)
     &   stop 'Qpferr'
      if(qstack_push(qst_auxline,idaov*rcw,dovyg).ne.idaov*rcw)
     &   stop 'Qpferr'
      if(qstack_push(qst_auxline,idaov*rcw,dovzg).ne.idaov*rcw)
     &   stop 'Qpferr'
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo2(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
      n1e = max(n1,ivvoo(nrep+1))
      allocate(evvoo(n1e * rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('evvoo allocation problem (tramo)')
      vovo1=0.0d0
      vovo2=0.0d0
      vovo3=0.0d0
      evvoo=0.0d0
!
!  sort the V_aj,ck ---> A_ak,cj (SRT 16), vovo1
!  sort the v_cb,ji (vvoo1) ---> B_cj,bi (1TT4), vovo2
!  contract A_ak,cj  x  B_cj,bi  to C_ak,bi , vovo3
!
      CALL getvovo(vovo3)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo3,vovo1)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = JVOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo1(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
!  backsort the C_ak,bi array (vovo3) to D_ab,ki (vovo2) via 1SS4 inverse
!  BKC from <||>,<||>   ---> <|<|,|>|>
!
      CALL SRT1SS4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,JVOVO,JJVO,JJVO,
     &              vovo3,vovo2)
!
!  create an energy denominator array in VVOO structure but with the
!  full range in a,b,i,j and multiply ti elementwise with the 
!  D_ab,ki array in vovo2. Result goes to vovo3
!
      if(IVVOO(nrep+1).ne.JVOVO(nrep+1)) then
        write(*,*) 'Array lenghts are mismatching for energy denoms!'
        stop
      endif
      n1 = JVOVO(nrep+1)
      call EARR_VVOO(eps,evvoo)
      Call MATXMAT(n1,evvoo,vovo2,vovo3) 
!
!  vovo3 has now the correct BKC to be used in a SRT 16 step from
!  c_ab,ki  <|<|,|>|>   --->  c_ai,bk with <||>,|><|  BKC
!
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,JVOOV,JJVO,JJOV,
     &            vovo3,vovo1)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOOV(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MOV(desrep)

      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovxg(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovyg(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovzg(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,7'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

      deallocate(vovo1,vovo2,vovo3)

!     write(*,*) '------  F_ai(2,7) finished'
!
!-------------------------------------
!------  F_ai(2,8): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo1(n1*rcw),stat=ialloc)
      allocate(vovo2(n1*rcw),stat=ialloc)
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem in F_ai(2,8)')
!
!  sort the v_ac,kj (vvoo1) ---> A_ak,cj (1TT4) (vovo1)
!  and the V*_ci,bj ---> B_cj,bk (SRT 16)  (vovo2)
!  contract A_ak,cj  and   B_cj,bi to C_ak,bi  (vovo3)
!
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo1,vovo1)
      CALL getvovo(vovo3)
      if(carith) CALL conjuga(ivovo(nrep+1),vovo3,1)
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NO,NV,NO,MVO,JVOVO,JJVO,JJVO,
     &            vovo3,vovo2)
      OFF1 = 1
      OFF3 = 1
      DO IREP = 1, NREP
         JREP = MULTB(IREP+NREP,1+NREP,2)
         M = MVO(IREP)
         N = MVO(IREP)
         K = MVO(JREP)
         OFF2 = JVOVO(JREP) * RCW + 1
         CALL XGEMM ('N','N',M,N,K,A1,vovo1(OFF1),M,vovo2(OFF2),K,
     &               A0,vovo3(OFF3),M)
         OFF1 = OFF1 + M * K * RCW
         OFF3 = OFF3 + M * N * RCW
      ENDDO
!
! C_ak,bi is now in vovo3, and we can proceed exactly in the same way
! as for F_ai(2,7)!
!
!
!  backsort the C_ak,bi array (vovo3) to D_ab,ki (vovo2) via 1SS4 inverse
!  BKC from <||>,<||>   ---> <|<|,|>|>
!
      CALL SRT1SS4 (NREP,MULTB,LFA,LTR,NV,NV,NO,NO,MVO,JVOVO,JJVO,JJVO,
     &              vovo3,vovo2)
!
!  use the energy denominator array from above and do the same
!  multiplications.
!
      if(IVVOO(nrep+1).ne.JVOVO(nrep+1)) then
        write(*,*) 'Array lenghts are mismatching for energy denoms!'
        stop
      endif
      n1 = JVOVO(nrep+1)
      Call MATXMAT(n1,evvoo,vovo2,vovo3)
!
!  vovo3 has now the correct BKC to be used in a SRT 16 step from
!  c_ab,ki  <|<|,|>|>   --->  c_ai,bk with <||>,|><|  BKC
!
      CALL SRT16 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,MVO,JVOOV,JJVO,JJOV,
     &            vovo3,vovo1)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOOV(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MOV(desrep)

      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovxg(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovyg(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vovo1(off1),M,
     &            dovzg(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,8'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

      deallocate(vovo1,vovo2,vovo3)

!     write(*,*) '------  F_ai(2,8) finished'
!
!-------------------------------------
!------  F_ai(2,9): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
! contraction of V_ab,cd and v_cd,ij is done in precontraction section
! of main part. We read in the VVOO precontracted integrals.
!
      n1 = ivvoott(nrep+1)
      allocate(vvoop(n1*rcw))
      open(ioch,file=preconfn,form='unformatted',
     &     access='sequential',status='unknown')
      read(ioch) (vvoop(K),K=1,n1*RCW)
      close(ioch)
      CALL DENOMVVOO (eps,vvoop,vvoop)
!
! The ~A_ab,ij ints are now in vvoop
! resort ~A_ab,ij to Bai,bj (vovo3) with <||>,<||> BKC.
! and contract with OV propints but in VO order.
! 
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoop,vovo3)

      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MVO(desrep)

      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)
      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,9'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

      deallocate(vvoop)
      deallocate(vovo3)

!     write(*,*) '------  F_ai(2,9) finished'
!
!-------------------------------------
!-----  F_ai(2,10): a DOV contribution
!-------------------------------------
!
      n1 = max(ivovo(nrep+1),jvovo(nrep+1))
      allocate(vovo3(n1*rcw),stat=ialloc)
      if(ialloc.ne.0)  call quit('vovo allocation problem (tramo)')
!
!  contract v_ab,kl and V_kl,ij to A_ab,ij (vvoo3)
!  from here we go on as in F_ai(2,9)
!
      n1=ioooott(nrep+1)
      allocate(oooo(n1*rcw))
      call getoooo(oooo)
      vvoo3 = 0.0d0
      CALL CNTRCT('N','N',NVVT,NOOT,NOOT,A1,vvoo1,oooo,A0,vvoo3,NREP)
!
!  next, divide A_ab,ij by the corresponding energy denominators
!  and resort it to B_ai,bj (vovo3)
!
      CALL DENOMVVOO(eps,vvoo3,vvoo3)
      CALL SRT1TT4 (NREP,MULTB,LFA,LFA,NV,NV,NO,NO,
     &              MVO,JVOVO,JJVO,JJVO,vvoo3,vovo3)
!
!  finally contract the B_ai,bj with the D_jb propints but in VO order!
!
      faipx = 0.0d0; faipy = 0.0d0; faipz = 0.0d0
      off1 = JVOVO(desrep)*rcw + 1
      off2 = JPVOXX(desrep)*rcw + 1
      M=MVO(desrep)
      K=MVO(desrep)

      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovx(off2),1,A0,faipx,1)

      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovy(off2),1,A0,faipy,1)

      CALL XGEMV ('N',M,K,A1,vovo3(off1),M,
     &            dovz(off2),1,A0,faipz,1)

      n1 = MVO(desrep)
      call rcaxpy(n1,carith,pphasex,faipx,faitx)
      call rcaxpy(n1,carith,pphasey,faipy,faity)
      call rcaxpy(n1,carith,pphasez,faipz,faitz)

      iorder='0th + 1st + 2ndABC + 2,10'
      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LFA,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)
      if(qstack_push(qst_prpline,lp1,faipx).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipy).ne.lp1) stop 'QE' 
      if(qstack_push(qst_prpline,lp1,faipz).ne.lp1) stop 'QE' 

      deallocate(oooo)
      deallocate(vovo3)

!     write(*,*) '------  F_ai(2,10) finished'
!
!  ----------------------------------------------------------
!  ----------------------------------------------------------
!  ----------------------------------------------------------
!  ----------------------------------------------------------
!  -------     2h2p second order contributions to F  --------
!  ----------------------------------------------------------
!  ----------------------------------------------------------
!  ----------------------------------------------------------
!  ----------------------------------------------------------
!
      call pst('calculating the 2h2p 2nd order contributions+')

!______________________________________________________________
!|
!|
!|  Generate index arrays for NTS two-particle contribs
!|
!
      allocate(LVVONO(NREP+1))
      allocate(LVNVOO(NREP+1))

      LVVONO  = 0   !offset array for nontotally symmetric VVO,O array
      LVNVOO  = 0   !offset array for nontotally symmetric V,VOO array

      nabxij = 0  ! total number of nts ab,ij  integrals
      nabixj = 0  ! total number of nts abi,j  integrals
      naxbij = 0  ! total number of nts a,bij  integrals
!|
!| fill the AB|IJ, ABI|J and A|BIJ offset arrays for NTS handling.
!|
!| In the abi,j and a,bij cases both ireps are fermionic whereas in the
!| ab,ij case both ireps are bosonic
!|
      DO JREP = 1, NREP

        IREP = MULTB (DESREP+NREP,JREP,2)
         if(multb(irep,jrep,2).ne.desrep) 
     &      stop '2p2h fermion symm error!'
        IREP_B = MULTB(DESREP+NREP,JREP+NREP,2)
         if(multb(irep_b+nrep,jrep+nrep,2).ne.desrep) 
     &      stop '2p2h bosonic symm error!'

        nabxij = nabxij + nvvt(irep_b)*noot(jrep)

        LVVONO(JREP+1) = LVVONO(JREP) + nvvot(irep)*no(jrep)
        nabixj = nabixj + nvvot(irep)*no(jrep)

        LVNVOO(JREP+1) = LVNVOO(JREP) + nv(irep)*nvoot(jrep)
        naxbij = naxbij + nv(irep)*nvoot(jrep)

      ENDDO !JREP
      write(*,*) 'size of F_abij(x,y,z) buffer:        ',nabxij
      write(*,*) 'size of vabij(x,y,z) buffer (part A):',nabixj
      write(*,*) 'size of vabij(x,y,z) buffer (part B):',naxbij
!|
!|______________________________________________________________
!
!
!  allocate NTS ab,ij buffers
!
      IF(nabxij.ne.xcoldet(desrep)) 
     &    stop 'fabij length does not match sat space length!'

      allocate(fabijtx(nabxij))
      allocate(fabijty(nabxij))
      allocate(fabijtz(nabxij))
      allocate(vvoop(nabxij*rcw))    ! buffer for nts vvoo stream
      fabijtx = A0
      fabijty = A0
      fabijtz = A0
      vvoop = 0.0d0

      lp2 = nabxij*rcw   ! actual length in reals of F_abij array (r/c)
!
!  allocate NTS abi,j buffers
!
      allocate(vabijx(nabixj*rcw))
      allocate(vabijy(nabixj*rcw))
      allocate(vabijz(nabixj*rcw))
      vabijx = 0.0d0
      vabijy = 0.0d0
      vabijz = 0.0d0
!
!---------------------------------------------------
!-----  F_ab,ij(2nd order PART A)  -----------------
!-----  see master sheet -TM 20-   -----------------
!---------------------------------------------------
!
!_________________________________________________________________
!|
!|   test forward antisymmetrization with nontotally symmetric stream
!|   and !!new!! SRT1T3_NTS routine.
 
      write(*,*) 'Performing antisymmetry test of SRT1T3_NTS...'

      if(carith) then
        do i=1,nabxij
          vvoop(2*i-1) = dble(i)
          vvoop(2*i)   = dble(i)
        enddo
      else
        do i=1,nabxij
          vvoop(i) = dble(i)
        enddo
      endif

      vabijx = 0.0d0
      CALL SRT1T3_NTS (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,
     &                LVVONO,KKVVOT,vvoop,vabijx,desrep,nabixj)
      xcntrl = 0.0d0

      if(carith) then
        do i=1,lvvono(nrep+1)
          xcntrl = xcntrl + vabijx(2*i-1) + vabijx(2*i)
        enddo
      else
        do i=1,lvvono(nrep+1)
          xcntrl = xcntrl + vabijx(i)
        enddo
      endif
      if(abs(xcntrl).gt.1.0E-16) stop '** FAILED **'
      write(*,*) '...passed.'
      vabijx = 0.0d0
!|
!|___________________________________________________________________
!
!  ... start production calculations
!
!  sort v_ab,ik to A_abi,k  (vvoo3)
!
      vvoo3 = 0.0d0
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             vvoo1,vvoo3)
!
!  contract A_abi,k (VVOO3) with D_kj   --->  B_abi,j (vabijx,y,z)
!
      OFF3 = 1
      DO JREP = 1, NREP    ! loop over rightmost irep of D_kj

        ABIREP = MULTB (DESREP+NREP,JREP,2)
        if(multb(abirep,jrep,2).ne.desrep)
     &       stop 'Part A fermionic symm. error'

        OFF1 = KVVOOT(ABIREP)*rcw + 1 
        OFF2 = (JPOOXX(desrep) + JJOO(ABIREP,JREP))*rcw + 1 

        M = NVVOT(ABIREP)
        N = NO(JREP)
        K = NO(ABIREP)

        CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,doox(OFF2),K,
     &              A0,vabijx(OFF3),M)
        CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,dooy(OFF2),K,
     &              A0,vabijy(OFF3),M)
        CALL XGEMM ('N','N',M,N,K,A1,vvoo3(OFF1),M,dooz(OFF2),K,
     &              A0,vabijz(OFF3),M)

        OFF3 = OFF3 + M * N * RCW

      ENDDO  !JREP
!
!  do antisymmetrization via nontotally symmetric sorter.
!
!  x,y,z results are subsequently stored in vvoop, therefore we have to push
!  immediately!
!
      n2 = LVVONO(nrep+1)
      vvoop = 0.0d0
      CALL SRT1T3_NTS (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,
     &                LVVONO,KKVVOT,vabijx,vvoop,desrep,n2)
      CALL rca_negate(nabxij,vvoop)
      call rcaxpy(nabxij,carith,pphasex,vvoop,fabijtx)
      if(qstack_push(qst_prpline,lp2,vvoop).ne.lp2) stop 'QE' 

      vvoop = 0.0d0
      CALL SRT1T3_NTS (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,
     &                LVVONO,KKVVOT,vabijy,vvoop,desrep,n2)
      CALL rca_negate(nabxij,vvoop)
      call rcaxpy(nabxij,carith,pphasey,vvoop,fabijty)
      if(qstack_push(qst_prpline,lp2,vvoop).ne.lp2) stop 'QE' 

      vvoop = 0.0d0
      CALL SRT1T3_NTS (NREP,MULTB,LTR,NVVT,NO,NO,NVVOT,
     &                LVVONO,KKVVOT,vabijz,vvoop,desrep,n2)
      CALL rca_negate(nabxij,vvoop)
      call rcaxpy(nabxij,carith,pphasez,vvoop,fabijtz)
      if(qstack_push(qst_prpline,lp2,vvoop).ne.lp2) stop 'QE' 

      iorder='0th-2nd order ph + 2p2h (part A)'
      incsat = nabxij

      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LTR,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      deallocate(vabijx)
      deallocate(vabijy)
      deallocate(vabijz)

!     write(*,*) 'F_ab,ij(2,part A) contributions finished.'
!
!---------------------------------------------------
!-----  F_ab,ij(2nd order PART B)  -----------------
!-----  see master sheet -TM 21-   -----------------
!---------------------------------------------------
!
!  allocate NTS a,bij buffers
!
      allocate(vabijx(naxbij*rcw))
      allocate(vabijy(naxbij*rcw))
      allocate(vabijz(naxbij*rcw))
      vabijx = 0.0d0
      vabijy = 0.0d0
      vabijz = 0.0d0
!_________________________________________________________________
!|
!|   test forward antisymmetrization with nontotally symmetric stream
!|   and !!new!! SRT1T2_NTS routine.
 
      write(*,*) 'Performing antisymmetry test of SRT1T2_NTS...'

      if(carith) then
        do i=1,nabxij
          vvoop(2*i-1) = dble(i)
          vvoop(2*i)   = dble(i)
        enddo
      else
        do i=1,nabxij
          vvoop(i) = dble(i)
        enddo
      endif

      vabijx = 0.0d0
      CALL SRT1T2_NTS (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVNVOO,LLVOOT,
     &                vvoop,vabijx,desrep)
      xcntrl = 0.0d0

      if(carith) then
        do i=1,lvnvoo(nrep+1)
          xcntrl = xcntrl + vabijx(2*i-1) + vabijx(2*i)
        enddo
      else
        do i=1,lvnvoo(nrep+1)
          xcntrl = xcntrl + vabijx(i)
        enddo
      endif
      if(abs(xcntrl).gt.1.0E-16) stop '** FAILED **'
      write(*,*) '...passed.'
      vabijx = 0.0d0
!|
!|___________________________________________________________________
!
!  ... start production calculations
!
!  sort v_cb,ij to A_c,bij (vvoo3)
!
      vvoo3 = 0.0d0
      CALL SRT1T2 (NREP,MULTB,LFA,NVVT,NV,NV,NOOT,LVVOOT,LLVOOT,
     &             vvoo1,vvoo3)
!
!  contract Dac with v_c,bij to B_a,bij (vabijx,y,z)
!
      OFF3 = 1
      DO BIJREP = 1, NREP   !loop over righmost irep of v_c,bij
         
         AREP = MULTB(DESREP+NREP,BIJREP,2)
         if(multb(arep,bijrep,2).ne.desrep) stop 'Part B symm. error'

         OFF1 = (JPVVXX(DESREP) + JJVV(AREP,BIJREP))*rcw + 1
         OFF2 = LVVOOT(BIJREP)*rcw + 1

         M = NV(AREP)
         N = NVOOT(BIJREP)
         K = NV(BIJREP)

         CALL XGEMM ('N','N',M,N,K,A1,dvvx(off1),M,vvoo3(off2),K,
     &               A0,vabijx(off3),m)
         CALL XGEMM ('N','N',M,N,K,A1,dvvy(off1),M,vvoo3(off2),K,
     &               A0,vabijy(off3),m)
         CALL XGEMM ('N','N',M,N,K,A1,dvvz(off1),M,vvoo3(off2),K,
     &               A0,vabijz(off3),m)

        OFF3 = OFF3 + M * N * RCW

      enddo  ! BIJREP
!
!  do antisymmetrization via nontotally symmetric sorter.
!
!  and, as above, push immediately.
!
      vvoop=0
      CALL SRT1T2_NTS (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,
     &                 LVNVOO,LLVOOT,vabijx,vvoop,desrep)
      call rcaxpy(nabxij,carith,pphasex,vvoop,fabijtx)
      if(qstack_push(qst_prpline,lp2,vvoop).ne.lp2) stop 'QE' 

      vvoop=0
      CALL SRT1T2_NTS (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,
     &                 LVNVOO,LLVOOT,vabijy,vvoop,desrep)
      call rcaxpy(nabxij,carith,pphasey,vvoop,fabijty)
      if(qstack_push(qst_prpline,lp2,vvoop).ne.lp2) stop 'QE' 

      vvoop=0
      CALL SRT1T2_NTS (NREP,MULTB,LTR,NVVT,NV,NV,NOOT,
     &                 LVNVOO,LLVOOT,vabijz,vvoop,desrep)
      call rcaxpy(nabxij,carith,pphasez,vvoop,fabijtz)
      if(qstack_push(qst_prpline,lp2,vvoop).ne.lp2) stop 'QE' 

      iorder='0th-2nd order ph + 2p2h (part A+B)'
      incsat = nabxij

      CALL TM_GEN(ioch,desrep,incsat,iorder,
     &            LTR,
     &            faitx,faity,faitz,
     &            fabijtx,fabijty,fabijtz)

      deallocate(vabijx)
      deallocate(vabijy)
      deallocate(vabijz)

      write(*,*) 'F_ab,ij(2,part B) contributions finished.'

      deallocate(fabijtx,fabijty,fabijtz)
      deallocate(vvoop)
!______________________________________________________
!|
!|
!|
!|
!|
!|  END CALCULATIONS, Deallocate all remaining arrays
!|
!|
!|
!|
!|_____________________________________________________
!
!
      deallocate(abj,aak,aki,aac)

      deallocate(vvoo3,vvoo4)
      deallocate(vvoo1,vvoo2)

      deallocate(doox,dooy,dooz)
      deallocate(dovx,dovy,dovz)
      deallocate(dvox,dvoy,dvoz)
      deallocate(dvvx,dvvy,dvvz)

      deallocate(dovxg,dovyg,dovzg)

      deallocate(faipx,faipy,faipz)
      deallocate(faitx,faity,faitz)

      deallocate(pphase)

      deallocate(JPOOXX,JPVOXX,JPOVXX,JPVVXX)
!
!  test the new tm out routine. It is maximally standalone.
!  and uses the F-vectors from the stack.
!
      call TM_OUT(desrep,pphasex,pphasey,pphasez)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TM_GEN(ioch,desrep,isatctl,iorder,
     &                  do2h2p,
     &                  faitx,faity,faitz,
     &                  fabijtx,fabijty,fabijtz)

      implicit none
!
!  this routine contracts the Davidson eigenvector with the 
!  F_ai/F_abij matrix elements to form the x,y,z transition moments.
!
!  The routine can be called for each order of perturbation theory in the
!  F_ai matrix elements. It just needs access to the Davidson eigenvectors.
!
!  calling parameters
!
      integer                      :: ioch,desrep,isatctl
      character*40                 :: iorder
      logical                      :: do2h2p
      complex*16, dimension(:)     :: faitx,faity,faitz
      complex*16, dimension(:)     :: fabijtx,fabijty,fabijtz

!
!  includes
!
#include  "../relccsd/complex.inc"
!
!  locals
!
      character*8, parameter                  :: basename='ADCEVECS'
      character*11                            :: davevcname
      logical                                 :: isthere

      real*8, allocatable, dimension(:)       :: rcoef
      complex*16, allocatable, dimension(:)   :: ccoef
      complex*16,allocatable,dimension(:)     :: tmx,tmy,tmz

      integer                                 :: ioresult
      integer                                 :: nroots,nmain,ladc
      real*8,allocatable,dimension(:)         :: evl,xnr
      real*8                                  :: ps
      real*8, parameter                       :: D0=0.0d0
      complex*16, parameter                   :: A0=(D0,D0)
      integer                                 :: i,k,ist,iln,numofstacks
      complex*16                              :: ax,ay,az
      real*8, parameter                       :: autoev = 27.2113957D0
      integer                                 :: n2h2p
      character*7,dimension(2)                :: numtype
!
!  executable
!
      if(carith.and.RCW.ne.2) STOP 'Carith error'
      if(.not.carith.and.RCW.ne.1) STOP 'Carith error'
      numtype(1) = 'real'
      numtype(2) = 'complex'
      IF(desrep.GT.9) THEN
        WRITE(davevcname,'(A8,A1,I2)') basename,'.',desrep
      ELSE
        WRITE(davevcname,'(A8,A2,I1)') basename,'.0',desrep
      ENDIF
      INQUIRE(file=davevcname,exist=isthere)
      if(.not.isthere) then
        write(*,*) '***************************************'
        write(*,*) '*****  Warning:                   *****'
        write(*,*) '***************************************'
        write(*,*) 'No eigenvectors found in symmetry',desrep,'!'
        write(*,*) 'No TMs accessible!'
        write(*,*)
        return
      endif
      OPEN(ioch,file=davevcname,access='sequential',
     &     status='unknown',form='unformatted')
      REWIND(ioch)

!_________________________  perform the contractions
!|
!| If we are in complex arithmetic, the eigenvector obtained
!| from Davidson is to ne complex conjugated before contraction
!|

      READ(ioch,iostat=ioresult) nroots
      READ(ioch,iostat=ioresult) nmain
      READ(ioch,iostat=ioresult) ladc
      n2h2p = ladc-nmain

!     write(*,*) 'Found',nroots,' states in the eigenvector file.'
!     write(*,*) 'Length of ADC matrix (from file):',ladc
!     write(*,*) 'Length of main space (from file):',nmain
!     write(*,*) 'Length of sat space:',n2h2p

      if(do2h2p.and.isatctl.ne.n2h2p) then
        call quit('Internal error for coupling block length')
      endif
!
! allocate eigenvalue, transition moment and coefficient arrays
!
      allocate(evl(nroots))
      evl = D0

      allocate(tmx(nroots))
      allocate(tmy(nroots))
      allocate(tmz(nroots))
      tmx = A0
      tmy = A0
      tmz = A0
      if(carith) then
        allocate(ccoef(ladc))
      else
        allocate(rcoef(ladc))
      endif

      allocate(xnr(nroots))
      xnr = D0

      Do k=1,nroots
        read(ioch) ist
        read(ioch) evl(k)
        read(ioch) ps
!       write(*,*)
!       write(*,'(A,X,I3,X,A,F10.4,A,X,F10.4,A,F10.4)') 
!    &    'State',ist,':exc. energy',evl(k),'(au)',evl(k)*autoev,
!    &    ' (eV) and PS',ps

        if(carith) then
          read(ioch)(ccoef(i),i=1,ladc)
          do i=1,ladc
            xnr(k)=xnr(k) + 
     &         (real(ccoef(i))**2.0d0 + aimag(ccoef(i))**2.0d0)
          enddo
!
!  do c.c.
!
          ccoef = dconjg(ccoef)

        else
          read(ioch)(rcoef(i),i=1,ladc)
          do i=1,ladc
            xnr(k)=xnr(k) + rcoef(i)**2.0d0
          enddo
        endif
        xnr(k) = dsqrt(xnr(k))
!
! the do2h2p flag controls the level of contraction.
!
        ax=A0
        ay=A0
        az=A0
        if(carith) then
          do i=1,nmain
            ax = ax + ccoef(i) * faitx(i)
            ay = ay + ccoef(i) * faity(i)
            az = az + ccoef(i) * faitz(i)
          enddo
          if(do2h2p) then
            do i=1,n2h2p
              ax = ax + ccoef(i+nmain) * fabijtx(i)
              ay = ay + ccoef(i+nmain) * fabijty(i)
              az = az + ccoef(i+nmain) * fabijtz(i)
            enddo
          endif
        else
          do i=1,nmain
            ax = ax + dcmplx(rcoef(i),D0)*faitx(i)
            ay = ay + dcmplx(rcoef(i),D0)*faity(i)
            az = az + dcmplx(rcoef(i),D0)*faitz(i)
          enddo
          if(do2h2p) then
            do i=1,n2h2p
              ax = ax + dcmplx(rcoef(i+nmain),D0)*fabijtx(i)
              ay = ay + dcmplx(rcoef(i+nmain),D0)*fabijty(i)
              az = az + dcmplx(rcoef(i+nmain),D0)*fabijtz(i)
            enddo
          endif
        endif
        tmx(k) = ax
        tmy(k) = ay
        tmz(k) = az

      Enddo   ! nroots

      CLOSE(ioch)
!|________________________ end reading section

!_________________________  print output
!|
!|
      write(*,*)
      write(*,'(A35,A40)') '     *******************   Level:  ',
     &   iorder
      write(*,*)
      write(*,'(A,4X,A,2X,A,3X,A,9X,A,10X,A,10X,A)')
     &       'State','Energy (au)','Energy(eV)','Norm',
     &       'TM(x)_r        TM(x)_i',
     &       'TM(y)_r        TM(y)_i',
     &       'TM(z)_r        TM(z)_i'
      write(*,'(3A)') '----------------------------------------------',
     &  '----------------------------------------------',
     &  '----------------------------------------------'

      do k=1,nroots
        write(*,'(I4,F14.7,F14.7,F7.2,A,2F15.9,A,2F15.9,A,2F15.9)') k,
     &   evl(k),evl(k)*autoev,xnr(k),'  : ',
     &   tmx(k),' | ',tmy(k),' | ',tmz(k)
      enddo
      write(*,'(3A)') '----------------------------------------------',
     &  '----------------------------------------------',
     &  '----------------------------------------------'
!|
!|________________________ end print output

      if(carith) then
        deallocate(ccoef)
      else
        deallocate(rcoef)
      endif
      deallocate(tmx)
      deallocate(tmy)
      deallocate(tmz)

      deallocate(evl)
      deallocate(xnr)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TM_OUT(desrep,phasex,phasey,phasez)

      use qstack
 
      implicit none
!
!  this is a final output routine summarizing the results for the obtained
!  transition moments in symmetry desrep. 
!  Hereby the F_ai and F_ab,ij intermediates in x,y,z
!  are fetched from the stack and the TM contributions are listed for each state. 
!  The eigenvectors
!  in contrast to TM_GEN therefore have to be read only once. This is the method
!  of choice if it comes to long eigenvectors.
!
!  calling parameters
!
      integer           :: desrep
      complex*16        :: phasex,phasey,phasez

! 
! remark: also in relativistic real symmetries the phases of the property
!         integrals can be complex and we have to account for that.
!         Resulting TMs therefore are generally defined as complex numbers
!         in the final output.
!
!  includes
!
#include  "../relccsd/complex.inc"
#include  "polprp_stacklines.h"
!
!  defines
!
#define FAI_CONTRIBS 15
#define FABIJ_CONTRIBS 2
!
!---------------interface region --------------------------------------
!
      interface

        Function get_file_unit()
          integer                      :: get_file_unit
        END Function

        Subroutine writetmln(l1,ch1,i1,c1,c2,c3)
          logical                      :: l1
          Character(len=11)            :: ch1
          integer                      :: i1
          complex*16                   :: c1,c2,c3
        END Subroutine

      end interface
!
!------------------------- locals -----------------------------------
!
      real*8, parameter                       :: autoev = 27.2113957D0
      character*8, parameter                  :: evcnm_root='ADCEVECS'
      character*11                            :: evc_fname
      character*8                             :: out_fname='ADCTRMOM'
      logical                                 :: isthere

      integer                                 :: io1,io2
      integer                                 :: nroots,nmain,ladc,n2p2h

      complex*16, allocatable, dimension(:)   :: ccoef
      real*8, allocatable, dimension(:)       :: rcoef

      integer                                 :: ist !root counter
      real*8                                  :: evl ! kth eigenvalue
      real*8                                  :: ps  ! kth pole strength
      real*8                                  :: xnr ! kth norm
  
      integer                                 :: i,k,iroot,ioresult
      complex*16                              :: ax,ay,az  !local contribs per root
      complex*16                              :: tmx,tmy,tmz !global contribs per root
      complex*16                              :: fosc
      real*8                                  :: foscr

      real*8, allocatable,dimension(:,:)      :: faipx,faipy,faipz
      real*8, allocatable,dimension(:,:)      :: fabijpx,fabijpy,fabijpz
      real*8, parameter                       :: D0=0.0d0
      complex*16, parameter                   :: A0=(D0,D0)
      character(len=11),dimension(18)         :: iorder

      Integer                                 :: qst_prpline
      Integer                                 :: aux_prpline
!
!  executable
!
      call pst('Calculating transition moments (ADCTRMOM file)+')
      iorder(1)  = '0th        '
      iorder(2)  = '1st        '
      iorder(3)  = '2ndA       '
      iorder(4)  = '2ndB       '
      iorder(5)  = '2ndC       '
      iorder(6)  = '2nd 2,1    '
      iorder(7)  = '2nd 2,2    '
      iorder(8)  = '2nd 2,3    '
      iorder(9)  = '2nd 2,4    '
      iorder(10) = '2nd 2,5    '
      iorder(11) = '2nd 2,6    '
      iorder(12) = '2nd 2,7    '
      iorder(13) = '2nd 2,8    '
      iorder(14) = '2nd 2,9    '
      iorder(15) = '2nd 2,10   '
      iorder(16) = '2nd 2p2h A '
      iorder(17) = '2nd 2p2h B '
      iorder(18) = 'total      '
!
!  assign number of stacklines in order to avoid direct integer assignments
!
      qst_prpline = PRP_STACKLINE
      aux_prpline = AUX_STACKLINE
!
! generate symmetry-specific file name for eigenvectors
!
      IF(desrep.GT.9) THEN
        WRITE(evc_fname,'(A8,A1,I2)') evcnm_root,'.',desrep
      ELSE
        WRITE(evc_fname,'(A8,A2,I1)') evcnm_root,'.0',desrep
      ENDIF
      INQUIRE(file=evc_fname,exist=isthere)
      if(.not.isthere) then
        write(*,*) '***************************************'
        write(*,*) '*****  Warning:                   *****'
        write(*,*) '***************************************'
        write(*,*) 'No eigenvectors found in symmetry',desrep,'!'
        write(*,*) 'No TMs accessible!'
        write(*,*)
        return
      endif
!
! open eigenvector file
!
      io1 = get_file_unit()
      OPEN(io1,file=evc_fname,access='sequential',
     &     status='unknown',form='unformatted')
      rewind(io1)
      write(*,*) '-- opened eigenvector file: ',evc_fname
!
! open master output file and append data if it exists
!
      io2 = get_file_unit()
      INQUIRE(file=out_fname,exist=isthere)
      If(isthere) then  !append
        OPEN(io2,file=out_fname,access='sequential',
     &       position='append',status='old',form='formatted')
        write(*,*) '-- opened output file: ',out_fname
      Else
        OPEN(io2,file=out_fname,access='sequential',
     &       status='new',form='formatted')
        write(*,*) '-- created output file: ',out_fname
      Endif

!_________________________  perform the contractions
!|
!|
!|  header information on eigenvector file:
!|
      READ(io1,iostat=ioresult) nroots
      READ(io1,iostat=ioresult) nmain
      READ(io1,iostat=ioresult) ladc
      n2p2h = ladc-nmain

      write(*,*)
      write(*,*) 'Treating symmetry:',desrep
      write(*,*) 'Found',nroots,' states in the eigenvector file.'
      write(*,*) 'Length of ADC matrix (from file):',ladc
      write(*,*) 'Length of main space (from file):',nmain
      write(*,*) 'Length of sat space:',n2p2h
      write(*,*)
!
! allocate required coefficient arrays (not possible earlier)
!
      if(carith) then
        allocate(ccoef(ladc))
      else
        allocate(rcoef(ladc))
      endif
!
! allocate F_ai and F_abij arrays (not possible earlier)
! these arrays are declared as real. In the complex case they acquire
! twice the size. In relativistic real cases, however, possibly
! complex property phases exist and they have to be accounted for 
! in the contractions properly which means that contraction results are
! always treated as complex numbers.
!
      allocate(faipx(nmain*rcw,FAI_CONTRIBS))
      allocate(faipy(nmain*rcw,FAI_CONTRIBS))
      allocate(faipz(nmain*rcw,FAI_CONTRIBS))

      allocate(fabijpx(n2p2h*rcw,FABIJ_CONTRIBS))
      allocate(fabijpy(n2p2h*rcw,FABIJ_CONTRIBS))
      allocate(fabijpz(n2p2h*rcw,FABIJ_CONTRIBS))
!
! pull contribs from stack, use popf so we get a FIFO type
! and we can stick to natural ordering.
! Stackline is PRP_STACKLINE.
!
      Do k=1,FAI_CONTRIBS
        if(qstack_popf(qst_prpline,faipx(:,k)).ne.nmain*rcw) stop 'QE'
        if(qstack_popf(qst_prpline,faipy(:,k)).ne.nmain*rcw) stop 'QE'
        if(qstack_popf(qst_prpline,faipz(:,k)).ne.nmain*rcw) stop 'QE'
      EndDo
      Do k=1,FABIJ_CONTRIBS
        if(qstack_popf(qst_prpline,fabijpx(:,k)).ne.n2p2h*rcw)
     &        stop 'QE'
        if(qstack_popf(qst_prpline,fabijpy(:,k)).ne.n2p2h*rcw)
     &        stop 'QE'
        if(qstack_popf(qst_prpline,fabijpz(:,k)).ne.n2p2h*rcw)
     &        stop 'QE'
      EndDo

! ****************** start writing to ADCTRMOM *******************

      write(io2,'(A)') ' '
      write(io2,'(A)') '*********************************************'
      write(io2,'(A,I12)') '**  Final State Symmetry: ',desrep
      write(io2,'(A)') '*********************************************'
      write(io2,'(A)') ' '

!*******************************************
!**** MASTER LOOP **************************
!*******************************************

      Do iroot=1,nroots
        read(io1) ist
        read(io1) evl
        read(io1) ps

        write(io2,'(A)') ' '
        write(io2,'(A,X,I3,X,A,X,I3,X,A,F12.6,A,X,F12.6,A,F12.6,A)') 
     &    'State',ist,'symm',desrep,'exc. energy',evl,' (au)',
     &    evl*autoev,' (eV)   PS',ps,'   @E'
        write(io2,'(A,15X,A,10X,A,10X,A)')
     &       'Order', 'TM(x)_r        TM(x)_i',
     &       'TM(y)_r        TM(y)_i','TM(z)_r        TM(z)_i'

!  read real or complex coefficients, form norm
!  In case of complex algebra we have to complex conjugate the
!  eigenvectors before contraction.

        if(carith) then
          read(io1)(ccoef(i),i=1,ladc)
          do i=1,ladc
            xnr=xnr + 
     &         (real(ccoef(i))**2.0d0 + aimag(ccoef(i))**2.0d0)
          enddo
          ccoef = dconjg(ccoef)
        else
          read(io1)(rcoef(i),i=1,ladc)
          do i=1,ladc
            xnr=xnr + rcoef(i)**2.0d0
          enddo
        endif
        xnr = dsqrt(xnr)
!
!****   CONTRIBUTION LOOP F_ai   ****
!  
        tmx=A0
        tmy=A0
        tmz=A0
        Do k = 1,FAI_CONTRIBS

         ax=A0
         ay=A0
         az=A0
         if(.not.carith) then
           do i=1,nmain
            ax = ax+dcmplx(rcoef(i),D0)*phasex*dcmplx(faipx(i,k),0.0d0)
            ay = ay+dcmplx(rcoef(i),D0)*phasey*dcmplx(faipy(i,k),0.0d0)
            az = az+dcmplx(rcoef(i),D0)*phasez*dcmplx(faipz(i,k),0.0d0)
           enddo
         else
           do i=1,nmain
            ax = ax + ccoef(i) * dcmplx(faipx(2*i-1,k),faipx(2*i,k))
            ay = ay + ccoef(i) * dcmplx(faipy(2*i-1,k),faipy(2*i,k))
            az = az + ccoef(i) * dcmplx(faipz(2*i-1,k),faipz(2*i,k))
           enddo
         endif
         tmx = tmx + ax
         tmy = tmy + ay
         tmz = tmz + az

         Call writetmln(.false.,iorder(k),io2,tmx,tmy,tmz)

        Enddo   !k for P/H contribs

! --------------- now add 2p2h contributions -------------------

        Do k = 1,FABIJ_CONTRIBS
         ax=A0
         ay=A0
         az=A0
         if(.not.carith) then
           do i=1,n2p2h
         ax=ax+dcmplx(rcoef(i+nmain),D0)*phasex*
     &          dcmplx(fabijpx(i,k),0.0d0)
         ay=ay+dcmplx(rcoef(i+nmain),D0)*phasey*
     &          dcmplx(fabijpy(i,k),0.0d0)
         az=az+dcmplx(rcoef(i+nmain),D0)*phasez*
     &          dcmplx(fabijpz(i,k),0.0d0)
           enddo
         else
           do i=1,n2p2h
         ax = ax+ccoef(i+nmain)*
     &          dcmplx(fabijpx(2*i-1,k),fabijpx(2*i,k))
         ay = ay+ccoef(i+nmain)*
     &          dcmplx(fabijpy(2*i-1,k),fabijpy(2*i,k))
         az = az+ccoef(i+nmain)
     &          *dcmplx(fabijpz(2*i-1,k),fabijpz(2*i,k))
           enddo
         endif

         tmx = tmx + ax
         tmy = tmy + ay
         tmz = tmz + az

         Call writetmln(.false.,iorder(k+FAI_CONTRIBS),
     &                  io2,tmx,tmy,tmz)

        Enddo   !k for 2P2H contribs

        Call writetmln(.true.,iorder(18),io2,tmx,tmy,tmz)
        fosc = tmx*dconjg(tmx) + tmy*dconjg(tmy) + tmz*dconjg(tmz)
        if(dabs(dimag(fosc)).gt.1.0E-14)
     %   write(io2,'(A)') '*** Oscill. strength acquires imag contribs!'
        foscr = (2.0d0/3.0d0)*evl*real(fosc)
        write(io2,'(A,F16.8,A)') '   Osc. strength:     -->>  ',
     &           foscr,'  a.u.      @O'
         write(io2,'(3A)') '-----------------------------------------',
     &  '-----------------------------------------',
     &  '--------------------------'
         write(io2,'(A)') ' '
         write(io2,'(A)') ' '
         write(io2,'(A)') ' '

         write(*,*) 'Root #',iroot,' written to ',out_fname

      Enddo   ! iroot
      CLOSE(io1)
      CLOSE(io2)
!|
!|________________________ end print output

      if(carith) then
        deallocate(ccoef)
      else
        deallocate(rcoef)
      endif

      deallocate(faipx)
      deallocate(faipy)
      deallocate(faipz)
      deallocate(fabijpx)
      deallocate(fabijpy)
      deallocate(fabijpz)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE writetmln(last,iorder,io,tx,ty,tz)
!
! this routine writes out a contribution line to the transition moments
! If last line then it sums up and also calculates the oscillator strength
! in x,y,z for the state under consideration.

      logical                            :: last
      character*11                       :: iorder
      Integer                            :: io
      Complex*16                         :: tx,ty,tz

      if(.not.last) then
        write(io,'(A14,2F15.9,A,2F15.9,A,2F15.9)')  iorder,
     &        tx,' | ',ty,' | ',tz
      else
         write(io,'(3A)') '-----------------------------------------',
     &  '-----------------------------------------',
     &  '--------------------------'
        write(io,'(A14,2F15.9,A,2F15.9,A,2F15.9)')  iorder,
     &        tx,' | ',ty,' | ',tz
      endif

      return
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

