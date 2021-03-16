CMI AMFI include files:
C=====================================================================================
C in AMFI directory: 
C   alphac.h  datapow.h  dofuc.h  ired.h  para.h  param.h  Regge.h  scfarr.h  
C
C in DIRAC include directory:
C     dummy.h implicit.h pi.h priunit.h 

      Subroutine amfi(IN,LUPROP,BREIT,FINITE,EXP_FIN,
     *                 CVEL_AMFI,EXAMF,IC_AMFIX,IACH,WRK,LFREE) 
CBS
CBS IN:  Input file, to be replaced by direct reading from DALTON arrays..
CBS LUPROP: Unit for writing the atomic integrals 
CBS BREIT: FLAG to switch to Breit-Pauli (Douglas-Kroll is the default)
CBS FINITE: flag whether to use a finite nucleus or not ...
CB  EXP_FIN: the finite nucleus exponent (if required)
CBS WRK, KFREE,LFREE standard work-array parameters in DALTON
CBS
CBS
CMI  IC_AMFIX ... specify what SO to include:
CMI  IACH     ... charge (int) of the nucleus - for modification of the occupation numbers in relscf !
CMI  CVEL     ... speed of light (in au) entering from DIRAC
CMI              (is in alphac.h )
CMI
CMI=====================================================================
#include "implicit.h"
#include "priunit.h"
#include "alphac.h"
c###########################################################################
c
c          A M F I 
c
c    Atomic Mean-Field Spin-Orbit Integral Program   
c
c Integral-code to generate the one- and two-electron spin-orbit integrals 
c in the no-pair approximation for an atom.
c 
c basis set is built by atomic functions of the form:
c 
c     f(r,Omega)= r**l Y_(lm) (Omega) 
c
c Allthough the code is created with a lot of care and love for 
c the details, the author doesn't give any warranty for it's 
c correctness. 
c
c B.Schimmelpfennig  Fysikum/Stockholm Summer 1996 
c
c If you use this code, please honour the authors work 
c by citing this work properly. 
c
c The author would like to thank the Deutsche Forschungsgemeinschaft 
c for financing this project by a Forschungsstipendium.
c
c
c   The spatial integrals are expected to be used with a spin part 
c   expressed in Paulis spin-matrices rather than with the Spin-operator
c   itself. So if a factor of two is somehow missing, check whether the 
c   same form of the operator is used.
c
c
c   WARNING !!!   WARNING !!   WARNING !!  WARNING !!   WARNING !!
c 
c   when writing spin-same-orbit and spin-other-oribt with sigma_i:
c 
c   For the spin-other-orbit-integrals particle 1 and 2 are exchanged
c   on the arrays carteXOO,carteYOO,carteZOO!!!!!!!!!
c 
c   The reason is to use most of the same-orbit part again and to 
c   have the same symmetry for the integrals on the arrays.
c
c
c   if the spin-other-orbit-part is used in the formulation with 
c   sigma_j, the particles are of cause not interchanged.
c
c
c
c   (i|HSO_mean|j) = (ij) + 1/2 * sum_M  occ(M) {
c                   2(ij|MM)_same - (iM|jM)_same -2(iM|jM)_other
c                   + (jM|iM)_same +2(jM|iM)_other } 
c
c   in the subroutines some signs are changed  to reorder indices 
c   in the integrals to (iM|jM) or (Mi|Mj) accoding to the way they 
c   were calculated before. 
c
c
c
c   one-particle integrals (really one-particle or mean-field) 
c   are written to files in CONTANDMULT. Look there for information on 
c   the format of files.  
c
c
c  BUGS:  There is still a strange sign-error in the two-electron-integrals 
c  if one applies straight-forward the formulae of the documentation.
c  This problem has been solved by the the cheater...
c
c  Everybody is welcome to find the problem in the formulas ........
c
c  First reasonable results on Thallium (SD with frozen 5D) 14.10.96
c
c
c
c
c
c  Connection to MOLCAS: 
c  How wonderful, they normalize the functions exactly as I do, which 
c  means they use the correct linear combinations.  
c
c  Exponents and coefficients are expected in the MOLCAS-Format      
c  first exponents 
c  coefficients afterwards   
c
c                                           8.5.97   
c                                                           
c  New version for DALTON canibalized from the MOLCAS version september 2000
c                                                           
c###########################################################################
#include "para.h"
      logical keep    ! parameter to decide about keeping angular 
cbs                     ! integrals in memory 
      logical keepcart    ! parameter to decide about keeping cartesian
cbs                         ! integrals in memory 
      logical makemean   ! parameter to decide about generating a meanfield   
      logical bonn       ! if bonn is set, Bonn-approach for spin-other orbit
      logical breit      ! if breit is set, BREIT-PAULI only                 
      logical SAMEORB    ! parameter for same-orbit only 
      logical AIMP       ! parameter to delete CORE for AIMP     
      logical oneonly    ! parameter to use only oneelectron integrals    
      logical FINITE

      logical ADDSO1  ! Add also one-electron SO integrals
      LOGICAL DOEXIT ! Exit AMFi when charge in readbas is 0
      LOGICAL EXAMF ! if .true., leave AMFI after readbas
      character*4  symmetry  
#include "datapow.h"
      common ipowxyz(3,-Lmax:Lmax,0:Lmax)
      dimension WRK(LFREE)

#include "amfi_if.h"
c     logical amfi_verbose
c     amfi_verbose = .false.

CMI   ... when we have CVEL, do calculate ALPHAC, ALPHA2
      CVEL = CVEL_AMFI
      ALPHAC = 1.0D0/CVEL
      ALPHA2 = ALPHAC*ALPHAC
      speed2=CVEL*CVEL
      speed4=speed2*speed2
c##########################################################################
cbs  #####################################################################
cbs         version with all angular integrals in memory 
c         keep=.true.
cbs  #####################################################################
cbs         version without  all angular integrals in memory 
          keep=.false. 
cbs  #####################################################################
cbs         version without  all cartesian integrals in memory 
          keepcart=.false. 
cbs  #####################################################################
cbs         version with all cartesian integrals in memory 
c         keepcart=.true.
cbs  #####################################################################
cbs   initialize tables with double facultatives... 
      call inidf
cbs   move some powers of x,y,z to the right place   BEGIN 
cbs   check if Lpowmax is high enough..
      if (Lpowmax.lt.Lmax) then 
      CALL QUIT('AMFI: increase lpowmax and edit ixyzpow')
      endif 
      jrun=1
      do irun=0,Lmax
      do Mval=-irun,irun
      ipowxyz(1,Mval,irun)=ixyzpow(jrun)
      ipowxyz(2,Mval,irun)=ixyzpow(jrun+1)
      ipowxyz(3,Mval,irun)=ixyzpow(jrun+2)
      jrun=jrun+3
      enddo
      enddo
cbs   move some powers of x,y,z to the right place   END   
      if (FINITE) then 
        ifinite=1
      else
        ifinite=0
      endif 
cbs   read the input 
      call readbas(Lhigh,makemean,bonn,breit,
     *symmetry,AIMP,ncont4,numballcart,IN,
     *ifinite,EXP_FIN,DOEXIT) 

CMI   IF (EXAMF.and.amfi_verbose) THEN
      IF(EXAMF)THEN
        WRITE(LUPRI,*)
     &  ' skip explicit AMFI - reading AMFI integrals from file'//
     &  ' AOPROPER_MNF.xxx!'
        RETURN
      ENDIF

CMI    .... when AOPROPER file exit, leave the AMFI !

CMI    set the oneonly option !!
       IC_AMFI = IC_AMFIX
CMI    if (amfi_verbose) then
       if (IPR_AMFI.GE.1) then
          WRITE(LUPRI,*) 'AMFI: IC_AMFI=',IC_AMFI
          WRITE(LUPRI,*) 'AMFI: charge of atom IACH=',IACH
          WRITE(LUPRI,'(A,F18.7)') 'AMFI: the speed of light:',CVEL
          if(ifinite >  0) 
     &    write(lupri,*) ' using finite nucleus in AMFI'
       endif
       IF (DOEXIT) THEN
        IF (IC_AMFI.EQ.2.OR.IC_AMFI.EQ.3.OR.
     &      IC_AMFI.EQ.4.OR.IC_AMFI.EQ.5)  THEN
C          if (amfi_verbose) 
           WRITE(LUPRI,*) '...zero charge, so no mean-field orbitals !'
        ENDIF
        IC_AMFI = 1
C       if (amfi_verbose)
           WRITE(LUPRI,*) 'The charge read in MNF.INP is ZERO!'//
     &     ' IC_AMFI set to 1 !'
       ENDIF
       IF (IC_AMFI.EQ.1) THEN
        oneonly=.true.
        ADDSO1=.TRUE.
C       WRITE(LUPRI,*) 'AMFI: ONLY SO1!' 
        IF (DOEXIT) 
     &     WRITE(LUPRI,*) 'AMFI: charge=0! Zero integrals!'
       ELSE IF (IC_AMFI.EQ.2) THEN
        oneonly=.false.
        sameorb=.true.   
        ADDSO1=.FALSE.
C       WRITE(LUPRI,*) 'AMFI: ONLY MFSSO2 !' 
       ELSE IF (IC_AMFI.EQ.3) THEN
        oneonly=.false.
        sameorb=.false.   
        ADDSO1=.FALSE.
C       WRITE(LUPRI,*) 'AMFI: ONLY MFSO2 !' 
       ELSE IF (IC_AMFI.EQ.4) THEN
        oneonly=.false.
        sameorb=.true.   
        ADDSO1=.TRUE.
C       WRITE(LUPRI,*) 'AMFI: SO1 + MFSSO2 !' 
       ELSE IF (IC_AMFI.EQ.5) THEN
        oneonly=.false.
        sameorb=.false.   
        ADDSO1=.TRUE.
C       WRITE(LUPRI,*) 'AMFI: SO1 + MFSO2 !' 
       ELSE
         write(lupri,'(a,i3,a)') ' *** error in AMFI: 
     & IC_AMFI is out of range (1 <= IC_AMFI <= 5):',IC_AMFI,' ***' 
         CALL QUIT('AMFI: IC_AMFI out of range')
       ENDIF    

CMI    if (amfi_verbose) then
       if (IPR_AMFI.GE.1) then
          if (IC_AMFI.eq.1) then
             WRITE(LUPRI,*) ': ONLY SO1!'
          else if (IC_AMFI.le.3) then
             WRITE(LUPRI,*) ': ONLY MFSSO2 !'
          else if (IC_AMFI.le.5) then
             WRITE(LUPRI,*) ': SO1 + MFSO2 !'
          endif
       endif
cbs
 123  if (ifinite.eq.2) call finsub  
cbs
! Lhigh is the highest l-value in the basis set

      call genpowers(Lhigh) !generate powers of exponents and overlaps 
cbs   start generating modified contraction coefficients
cbs   generate starting adresses of contraction coefficients  on 
cbs   contrarray 
      call genstar(Lhigh)  
cbs   generate ovlp of normalized primitives 
      call genovlp(Lhigh) 
      do lrun=0,Lhigh
cbs   cont(L) arranges all the contraction coefficients for a given L-value
cbs   and renormalizes them 
      call cont(lrun,breit,ifinite)
      enddo 

      if (makemean.and.(.not.oneonly).and.ifinite.le.1) 
     *call getAOsx(Lhigh,IACH,WRK,LFREE)

cbs                        
cbs        beginning the angular part  
      if (.not.oneonly) then  
CBS   write(6,*) '***************************************************' 
CMI   if (amfi_verbose) write(LUPRI,*)
      if (IPR_AMFI.GE.2) write(LUPRI,*)
     &      '********   beginning the 2e-part ******************' 
CBS   write(6,*) '***************************************************' 
cbs                        
cbs  ##################################################################################### 
cbs  ##################################################################################### 
cbs  ##################################################################################### 
cbs
cbs    
      call angular(Lhigh,keep,keepcart,makemean,bonn,breit,
     *sameorb,ifinite,WRK,LFREE) ! subroutine for angular part 
      endif 
      if (ifinite.eq.1) then ! redo everything for finite core
CBS   write(6,*) 'once more the two-electron integrals'
      ifinite=2
      goto 123
      endif 
cbs ########################################################################################
cbs ########################################################################################
cbs ########################################################################################
CBS   write(6,*) '***************************************************' 
      if (IPR_AMFI.GE.2) write(LUPRI,*)
     &      '*******   beginning the 1-electron-part  **********' 
CBS   write(6,*) '***************************************************' 

cbs    the one-electron spin-orbit integrals   
      call gen1overR3(Lhigh)   ! generates the 1/r**3 integrals  for normalized functions 
      call contandmult(Lhigh,makemean,AIMP,oneonly,
     &     ADDSO1,numballcart,LUPROP,
     *     ifinite,WRK,LFREE) ! multiplies radial integrals with l,m-dependent
cbs                             factors and contraction coefficients 
CBS   write(6,*) '***************************************************' 
CMI   if (amfi_verbose) write(LUPRI,*)
      if (IPR_AMFI.GE.2) write(LUPRI,*)
     &         '*******   end of  the 1-electron-part    **********' 
CBS   write(6,*) '***************************************************' 
cbs ########################################################################################
cbs ########################################################################################
cbs ########################################################################################

      Return
      end

   
      subroutine finsub
cbs
cbs   subroutine to set up parameters for finite nucleus. The s-functions are replaced 
cbs   by just one exponent which models the nucleus.
cbs
#include "implicit.h"
#include "para.h"
#include "param.h"
      common /nucleus/ charge,Exp_finite    
      noccorb(0)=1  
      do l=1,lmax_occ
      noccorb(l)=0              
      enddo 
      occup(1,0)=-charge 
      nprimit_keep=nprimit(0)
      ncontrac_keep=ncontrac(0)
      nprimit(0)=1
      ncontrac(0)=1
      exponents(1,0)=0.5d0*Exp_finite       
      return 
      end  
      
      
      subroutine angular(Lhigh,keep,keepcart,makemean,bonn,
     *breit,sameorb,ifinite,WRK,LFREE)
c
cbs   COMBINES THE RADIAL INTEGRALS WITH THE ANGULAR FACTORS 
c
cbs   if keep=.true. then 
cbs   all the integrals will be kept in memory. 
cbs   Perhaps, there will be the option to make the 
cbs   transformation to the cartesian basis-sets 
cbs   everytime, they are required. 
cbs   Therefore, the integrals are kept in memory and 
cbs   can be further transformed, whenever required.   
cbs   in order not to waste to much memory, the atomic  
cbs   integrals are thrown away after each l,l,l,l-block 
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "param.h"
      logical keep,keepcart,icheck,mcheckxy,mcheckz,makemean,bonn,   
     *breiT,sameorb,cleaner,NFINI

#include "amfi_if.h"
c     logical amfi_verbose
cbs   NFINI means not finite nucleus  
      dimension l2block(0:Lmax,0:Lmax,0:Lmax,0:Lmax)
      dimension WRK(LFREE)
cbs #####################################################################
cbs   some preparation of factors needed later on..                     #
cbs ######################################################################
      ipnt(i,j)=(max(i,j)*max(i,j)-max(i,j))/2+min(i,j) 
      roottwo=dsqrt(2d0) 
c     amfi_verbose = .false.
cbs   calculate some prefactors that will be needed quite often      
      call prefac(Lmax,preroots,clebsch) 
        if (ifinite.ne.2) then 
cbs     clean array for one electron integrals
        iprod=MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)*Lmax     
        call dzero(onecartX,iprod) 
        call dzero(onecartY,iprod) 
        call dzero(onecartZ,iprod) 
        NFINI=.true.   
        else
        NFINI=.false.   
        endif 
cbs   generate an array with sign for (even/odd) m-values
      isignM(0)=1
      do I=2,Lmax,2
      isignM(I)=1
      isignM(-I)=1
      enddo
      do I=1,Lmax,2
      isignM(I)=-1
      isignM(-I)=-1
      enddo
cbs #####################################################################
cbs   prefactors preXZ und preY include the factors 1/root(2)
cbs   for the +/- linear combinations of spherical harmonics 
cbs #####################################################################
      do M4=-Lmax,Lmax
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,m2,m3,m4)=0.25d0
      enddo
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,m2,m3,0)=preXZ(m1,m2,m3,0)*roottwo  
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,m2,0,m3)=preXZ(m1,m2,0,m3)*roottwo  
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,0,m2,m3)=preXZ(m1,0,m2,m3)*roottwo  
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(0,m1,m2,m3)=preXZ(0,m1,m2,m3)*roottwo  
      enddo
      enddo
      enddo
      do M4=-Lmax,Lmax
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preY(m1,m2,m3,m4)=preXZ(m1,m2,m3,m4) 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   additional (-) signs from the (-i) factors  in the 
cbs   (-) linear combinations   (see tosigX(Y,Z).f)
cbs #####################################################################
cbs   + - - -   =>   minus 
      do M4=-Lmax,-1    
      do M3=-Lmax,-1    
         do M2=-Lmax,-1   
         do M1= 0,Lmax
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
cbs   - + - -   =>   minus 
         do M2= 0,Lmax
         do M1=-Lmax,-1   
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
      enddo
      enddo
      do M2= 0,Lmax   
      do M1= 0,Lmax
cbs   + + + -   =>   minus 
         do M4=-Lmax,-1    
         do M3= 0,Lmax    
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
cbs   + + - +   =>   minus 
         do M4= 0,Lmax    
         do M3=-Lmax,-1    
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
      enddo
      enddo
cbs   + +  - -  >   - 
      do M4=-Lmax,-1  
      do M3=-Lmax,-1  
      do M2=0,Lmax
      do M1=0,Lmax
	      preY(m1,m2,m3,m4)=-preY(m1,m2,m3,m4) 
      enddo
      enddo
      enddo
      enddo
cbs   - -  + +  >   - 
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=-Lmax,-1  
      do M1=-Lmax,-1  
	      preY(m1,m2,m3,m4)=-preY(m1,m2,m3,m4) 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   some quick decision for interaction    
cbs #####################################################################
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
	      icheck=mcheckxy(m1,m2,m3,m4)
	      icheckxy(m1,m2,m3,m4)=icheck 
	      icheckxy(m1,m2,m3,-m4)=icheck 
	      icheckxy(m1,m2,-m3,m4)=icheck 
	      icheckxy(m1,-m2,m3,m4)=icheck 
	      icheckxy(-m1,m2,m3,m4)=icheck 
	      icheckxy(m1,m2,-m3,-m4)=icheck 
	      icheckxy(m1,-m2,m3,-m4)=icheck 
	      icheckxy(m1,-m2,-m3,m4)=icheck 
	      icheckxy(m1,-m2,-m3,-m4)=icheck 
	      icheckxy(-m1,m2,m3,-m4)=icheck 
	      icheckxy(-m1,m2,-m3,m4)=icheck 
	      icheckxy(-m1,m2,-m3,-m4)=icheck 
	      icheckxy(-m1,-m2,m3,m4)=icheck 
	      icheckxy(-m1,-m2,m3,-m4)=icheck 
	      icheckxy(-m1,-m2,-m3,m4)=icheck 
	      icheckxy(-m1,-m2,-m3,-m4)=icheck 
      enddo 
      enddo 
      enddo 
      enddo 
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
	      icheck=mcheckz(m1,m2,m3,m4)
	      icheckz(m1,m2,m3,m4)=icheck 
	      icheckz(m1,m2,m3,-m4)=icheck 
	      icheckz(m1,m2,-m3,m4)=icheck 
	      icheckz(m1,m2,-m3,-m4)=icheck 
	      icheckz(m1,-m2,m3,m4)=icheck 
	      icheckz(m1,-m2,m3,-m4)=icheck 
	      icheckz(m1,-m2,-m3,m4)=icheck 
	      icheckz(m1,-m2,-m3,-m4)=icheck 
	      icheckz(-m1,m2,m3,m4)=icheck 
	      icheckz(-m1,m2,m3,-m4)=icheck 
	      icheckz(-m1,m2,-m3,m4)=icheck 
	      icheckz(-m1,m2,-m3,-m4)=icheck 
	      icheckz(-m1,-m2,m3,m4)=icheck 
	      icheckz(-m1,-m2,m3,-m4)=icheck 
	      icheckz(-m1,-m2,-m3,m4)=icheck 
	      icheckz(-m1,-m2,-m3,-m4)=icheck 
      enddo 
      enddo 
      enddo 
      enddo 
cbs #####################################################################
cbs   there are at most 16 possible combinations of signs ( 2**4) 
cbs #####################################################################
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
      do irun=1,16
      interxyz(irun,m1,m2,m3,m4)=0         
      enddo
      enddo
      enddo
      enddo
      enddo
cbs   the following M values are the ones from the cartesian 
cbs   linear combinations. interxyz gives the sign sequence 
cbs   for interacting spherical functions, starting with 
cbs   type 1 (++++) and ending with type 16 (-++-)   
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
      if (icheckxy(m1,m2,m3,m4).or.icheckz(m1,m2,m3,m4)) then 
          irun=0     
          if (iabs(m1+m2-m3-m4).le.1) then 
          irun=irun+1
          interxyz(irun,m1,m2,m3,m4)=1          ! + + + + 
	          if (m1.gt.0.and.m2.gt.0.and.
     *            m3.gt.0.and.m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=2  ! - - - - 
	          endif 
          endif 
          if (iabs(m1+m2-m3+m4).le.1) then 
	          if (m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=3  ! + + + - 
	          endif 
	          if (m1.gt.0.and.m2.gt.0.and.
     *            m3.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=4  ! - - - + 
	          endif 
          endif 
          if (iabs(m1+m2+m3-m4).le.1) then 
	          if (m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=5  ! + + - + 
	          endif 
	          if (m1.gt.0.and.m2.gt.0.and.
     *            m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=6  ! - - + - 
	          endif 
          endif 
          if (iabs(m1-m2-m3-m4).le.1) then 
	          if (m2.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=7  ! + - + + 
	          endif 
	          if (m1.gt.0.and.m3.gt.0.and.
     *            m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=8  ! - + - - 
	          endif 
          endif 
          if (iabs(-m1+m2-m3-m4).le.1) then 
	          if (m1.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=9  ! - + + + 
	          endif 
	          if (m2.gt.0.and.m3.gt.0.and.
     *            m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=10 ! + - - - 
	          endif 
          endif 
          if (iabs(m1+m2+m3+m4).le.1) then 
	          if (m3.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=11 ! + + - - 
	          endif 
	          if (m1.gt.0.and.m2.gt.0) then  
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=12 ! - - + +   
	          endif 
          endif 
          if (iabs(m1-m2-m3+m4).le.1) then 
	          if (m2.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=13 ! + - + - 
	          endif 
	          if (m1.gt.0.and.m3.gt.0) then  
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=14 ! - + - +   
	          endif 
          endif 
          if (iabs(m1-m2+m3-m4).le.1) then 
	          if (m2.gt.0.and.m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=15 ! + - - + 
	          endif 
	          if (m1.gt.0.and.m4.gt.0) then  
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=16 ! - + + -   
	          endif 
          endif 
      endif 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   isgnprod gives the sign due to powers (-1)**M  this are again 
cbs   angular m-values 
cbs #####################################################################
      do M4=-Lmax,Lmax
      if (M4.gt.0) then 
      inter4=isignM(M4) 
      else 
      inter4=1
      endif 
      do M3=-Lmax,Lmax
      if (M3.gt.0) then 
      inter3=inter4*isignM(M3) 
      else
      inter3=inter4 
      endif 
      do M2=-Lmax,Lmax
      if (M2.gt.0) then 
      inter2=inter3*isignM(M2) 
      else
      inter2=inter3 
      endif 
      do M1=-Lmax,Lmax
      if (M1.gt.0) then 
      isgnprod(m1,m2,m3,m4)=inter2*isignM(M1)
      else 
      isgnprod(m1,m2,m3,m4)=inter2
      endif   
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   some preparation of factors needed later on..  finished           # 
cbs #####################################################################
c
c
c
cbs   counter for total number of cartesian integrals                   !  set some counters 
      numbcart=0                                                        !
cbs   same orbit integrals integrals  on carteXSO carteYSO and carteSO                            
cbs   other orbit integrals  on carteXOO carteYOO and carteOO                            
      iangfirst=0 ! first block of angular integrals                   
cbs #####################################################################
cbs   loop over all (l,l,l,l) blocks generated in the radial part       # 
cbs #####################################################################
      do lrun4=0,Lmax
      do lrun3=0,Lmax
      do lrun2=0,Lmax
      do lrun1=0,Lmax
      l2block(lrun1,lrun2,lrun3,lrun4)=0
      enddo
      enddo
      enddo
      enddo
cbs   loop over all possible < l1 l2, l3 l4 > blocks
CMI   if (amfi_verbose) write(LUPRI,'(A)') '   L1   L2   L3   L4' 
      if (IPR_AMFI.GE.2) write(LUPRI,'(A)') '   L1   L2   L3   L4' 
      do l1=0,Lhigh   ! improving is probably possible...
      do l2=0,Lhigh
      do l3=0,l1
      do l4=0,l2
cbs   check parity
      if (mod(l1+l2+l3+l4,2).eq.0) then
cbs   check that Lleft and Lright do not always differ by more than one
cbs   a difference of two means two spin flips and is therefore not allowed
      Lleftmax=l1+l2
      Lrightmax=l3+l4
      Lleftmin=iabs(l1-l2)
      Lrightmin=iabs(l3-l4)
      if ((Lrightmin-Lleftmax.le.1.and.Lrightmax-Lleftmin.gt.-1).or.
     *(Lleftmin-Lrightmax.le.1.and.Lleftmax-Lrightmin.gt.-1)) then
cbs   additional check for mean-field
      if ((l1.eq.l3.and.l2.eq.l4).or.(l1.eq.l2.and.l3.eq.l4)) then
      if (l1+l3.ne.0) then
CMI      if (amfi_verbose) write(LUPRI,'(4I5)') l1,l2,l3,l4
         if (IPR_AMFI.GE.2) write(LUPRI,'(4I5)') l1,l2,l3,l4
CBS   now I determine the size of the angular integral arrays
        jblock=0
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3           
        m4=m1+m2-m3+1
        if (iabs(m4).le.l4) then
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then   
        jblock=jblock+1
        endif 
        endif 
        enddo 
        enddo 
        enddo 
        do m1=  0,l1
        do m2=-l2,l2
        do m3=-l3,l3 
        m4=m1+m2-m3
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then         
        if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then !  all m eqal 0 make no sense....
        if (iabs(m4).le.l4)  then                                               
        jblock=jblock+1
        endif 
        endif 
        endif 
        enddo 
        enddo 
        enddo 
CBS   done !!                                                     
cbs     number of contracted integrals for each block 
        ncont=ncontrac(l1)*ncontrac(l2)*
     *  ncontrac(l3)*ncontrac(l4) 
      mxangint=jblock*ncont
cbs   determine the size icont4 for the radial integrals 
      call gencoulDIM(l1,l2,l3,l4,makemean,bonn,breit,
     *sameorb,icont4) 
      IANGSO = 1
      iangOO=iangSO+mxangint
      icartSO=iangOO+mxangint   
      icartOO=icartSO+ncont
      iconSO=icartOO+ncont
      iconOO=iconSO+icont4   
      KLAST = ICONOO + ICONT4
CMI   IF (KLAST .GT. LFREE) CALL STOPIT('AMFI  ','ANGULAR',KLAST,LFREE)
      IF (KLAST .GT. LFREE) THEN 
        WRITE(LUPRI,*) 
        WRITE(LUPRI,*) 'angular: KLAST,LFREE:',KLAST,LFREE
        CALL QUIT(
     & 'angular: work space exceeded;allocate more memory in DIRAC !') 
      ENDIF
      LLEFT = LFREE - KLAST + 1
      call gencoul(l1,l2,l3,l4,makemean,bonn,breit,
     *sameorb,WRK(iconSO),WRK(iconOO),icont4,
     *WRK(KLAST),LLEFT) ! generates and transforms integrals
        l2block(l1,l2,l3,l4)=1  ! can be used for getting the
cbs   local counter for integral adresses  
        mblock=0 ! counter of (m,m,m,m)-blocks for (l1,l2,l3,l4)   
cbs     if keep is set to false, the angular integrals are 
cbs     thrown away after each block of l-values 
cbs     which means integrals start at address 0 
        if (.not.keep) iangfirst=0
        locstar=iangfirst ! local starting adress counter 
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
        do m4=-l4,l4
        mcombina(1,m1,m2,m3,m4)=0  ! will hold type of integrals (1,2,3) 
        mcombina(2,m1,m2,m3,m4)=0  ! will hold number of block 
        enddo 
        enddo 
        enddo 
        enddo 
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
cbs     m4 is more or less fixed by m1-3 
c####################################################################################
c####################################################################################
c########## the L- -type block to be combined with sigma+ ###########################
c####################################################################################
c####################################################################################
        m4=m1+m2-m3+1
        if (iabs(m4).le.l4) then !the  L- -block to be combined with sigma+
cbs     not all m-combinations are needed for the mean-field 
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then 
        mcombina(1,m1,m2,m3,m4)=1    
        mblock=mblock+1
        if (locstar+ncont.gt.mxangint) then 
        write(LUPRI,*)'not enough space allocated for angular integrals'
        write(LUPRI,*) 'increase mxangint to at least ',
     *  locstar+ncont        
        CALL QUIT('Out of dimensional bounds in AMFI')
        endif  
cbs mkangLmin = make_angular_integrals_for_L- type operator 
cbs really generates  the angular prefactors and combines them with 
cbs the radial integrals
        call mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       WRK(iangSO+locstar),
     *       WRK(iangOO+locstar),
     *       Lfirst(1),Llast(1),Lblocks(1),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       WRK(iconSO+Lstarter(1)-1),
     *       WRK(iconSO+Lstarter(2)-1),
     *       WRK(iconSO+Lstarter(3)-1),
     *       WRK(iconSO+Lstarter(4)-1),
     *       WRK(iconOO+Lstarter(1)-1),
     *       WRK(iconOO+Lstarter(2)-1),
     *       WRK(iconOO+Lstarter(3)-1),
     *       WRK(iconOO+Lstarter(4)-1),
     *       preroots,clebsch,scratch4,bonn,breit,
     *       sameorb) 
        locstar=locstar+ncont ! increase starting address 
        mcombina(2,m1,m2,m3,m4)=mblock  ! set the block number 
c####################################################################################
c####################################################################################
c########## the L+ -type block to be combined with sigma- ###########################
c####################################################################################
c####################################################################################
c
c   these integrals are obtained by changing the signs of the m-values.
c   As the integrals are the same, the pointer points to the same integrals...
c
c
        mcombina(1,-m1,-m2,-m3,-m4)=3
        mcombina(2,-m1,-m2,-m3,-m4)=mblock  
        endif 
        Endif 
        enddo 
        enddo 
        enddo 
c####################################################################################
c####################################################################################
c########## the L0 -type block to be combined with sigma0 ###########################
c####################################################################################
c####################################################################################
        do m1=  0,l1
        do m2=-l2,l2
        do m3=-l3,l3
cbs     m4 is more or less fixed by m1-3 
        m4=m1+m2-m3 ! the L0-block to be combined with sigma0 
cbs     not all m-combinations are needed for the mean-field 
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then 
c       
        if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then !  all m eqal 0 make no sense....
        if (iabs(m4).le.l4)  then    
        mcombina(1,m1,m2,m3,m4)=2
        mblock=mblock+1
        if (locstar+ncont.gt.mxangint) then 
        write(LUPRI,*)'not enough space allocated for angular integrals'
        write(LUPRI,*) 'increase mxangint to at least ',
     *  locstar+ncont        
        CALL QUIT('Out of dimensional bounds in AMFI')
        endif  
        call mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       WRK(iangSO+locstar),
     *       WRK(iangOO+locstar),
     *       Lfirst(1),Llast(1),Lblocks(1),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       WRK(iconSO+Lstarter(1)-1),
     *       WRK(iconSO+Lstarter(2)-1),
     *       WRK(iconSO+Lstarter(3)-1),
     *       WRK(iconSO+Lstarter(4)-1),
     *       WRK(iconOO+Lstarter(1)-1),
     *       WRK(iconOO+Lstarter(2)-1),
     *       WRK(iconOO+Lstarter(3)-1),
     *       WRK(iconOO+Lstarter(4)-1),
     *       preroots,clebsch,scratch4,bonn,breit,
     *       sameorb) 
        locstar=locstar+ncont
        mcombina(2,m1,m2,m3,m4)=mblock  
        endif 
        endif 
        endif 
        enddo
        enddo
        enddo
cbs  ##################################################################################
cbs  ##################################################################################
cbs     transformation to l,m dependent integrals is finished 
cbs  ##################################################################################
c
c
c
c
cbs  ##################################################################################
cbs     begin transformation to cartesian integrals 
cbs  ##################################################################################
cbs  ##################################################################################
cbs     check out, which combinations of m-values will 
cbs     contribute to cartesian integrals    
        do m1=-l1,l1       !    
        do m2=-l2,l2       ! these indices now run over the real harmonics      
        do m3=-l3,l3       !
        do m4=-l4,l4       !
        mcombcart(1,m1,m2,m3,m4)=0     ! will hold the type  x=1 y=2 z=3 
        mcombcart(2,m1,m2,m3,m4)=0     ! will hold the block number
        enddo 
        enddo 
        enddo 
        enddo 
        mblockx=0
        mblocky=0
        mblockz=0
        do m3=-l3,l3            
        do m4=-l4,l4       
cbs     if the l-values are the same : triangular matrix over m-values is sufficient       
        if (l1.eq.l3) then    
        m1upper=m3
        else
        m1upper=l1
        endif 
        if (makemean) m1upper=l1
cbs     if the l-values are the same : triangular matrix over m-values is sufficient       
        if (l2.eq.l4) then 
        m2upper=m4
        else
        m2upper=l2 
        endif 
        if (makemean) m2upper=l2
        do m1=-l1,m1upper   
        If (l1.eq.l3.and.m1.eq.m3) then ! clean real zeros by symmetry to be exactly zero 
cbs     this a problem of the spin-other-orbit integrals, as they are by formula 
cbs     not antisymmetric in the indices for particle 1. 
        cleaner=.true.
        else
        cleaner=.false.
        endif  
        do m2=-l2,m2upper   
cbs     not all m-combinations are needed for the mean-field 
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.m2.eq.m4).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.(m1.eq.m2.or.m3.eq.m4))) then 
C
        indx=ipowxyz(1,m1,l1)+ipowxyz(1,m2,l2)+
     *  ipowxyz(1,m3,l3)+ipowxyz(1,m4,l4)
        indy=ipowxyz(2,m1,l1)+ipowxyz(2,m2,l2)+
     *  ipowxyz(2,m3,l3)+ipowxyz(2,m4,l4)
        indz=ipowxyz(3,m1,l1)+ipowxyz(3,m2,l2)+
     *  ipowxyz(3,m3,l3)+ipowxyz(3,m4,l4)
        indx=mod(indx,2)
        indy=mod(indy,2)
        indz=mod(indz,2)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA X      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (indx.eq.0.and.indy.eq.1.and.indz.eq.1.and.      
     *  icheckxy(m1,m2,m3,m4)) then  ! Y*Z ->  transforms like  L_x (B1)
cbs     integrals for sigma_x 
        mblockx=mblockx+1
        mcombcart(1,m1,m2,m3,m4)=1
        mcombcart(2,m1,m2,m3,m4)=mblockx                               
        call tosigX(m1,m2,m3,m4,WRK(iangSO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartSO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breiT)) 
     *  call tosigX(m1,m2,m3,m4,WRK(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartOO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner)
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.and.
     *       iabs(m1+m3).eq.1.and.l1.ne.0) then 
             call two2mean13(WRK(icartSO),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartx(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.iabs(m4+m3).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartSO),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if(NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.iabs(m4+m3).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartSO),       
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34b(WRK(icartSO),    
     *       WRK(icartOO),    
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.iabs(m1+m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12a(WRK(icartSO),   
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12a(WRK(icartSO),     
     *       WRK(icartOO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.and.iabs(m1+m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12b(WRK(icartSO),    
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12b(WRK(icartSO),    
     *       WRK(icartOO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             endif 
             endif 
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Y      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.0.and.indz.eq.1.and.  
     *  icheckxy(m1,m2,m3,m4)) then  ! X*Z transforms like L_y  (B2) 
cbs     integrals for sigma_y 
        mblocky=mblocky+1
        mcombcart(1,m1,m2,m3,m4)=2
        mcombcart(2,m1,m2,m3,m4)=mblocky                               
        call tosigY(m1,m2,m3,m4,WRK(iangSO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartSO),preY,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breit)) 
     *  call tosigY(m1,m2,m3,m4,WRK(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartOO),preY,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.iabs(m3-m1).eq.1.and.l1.ne.0) then 
             call two2mean13(WRK(icartSO),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartY(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.iabs(m3-m4).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartOO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.iabs(m3-m4).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartOO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.iabs(m1-m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartOO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.anD.Iabs(m1-m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartOO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             endif 
             endif 
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Z      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.1.and.indz.eq.0.and.  
     *  icheckz(m1,m2,m3,m4)) then ! X*Y transforms like L_z  (A2) 
cbs     integrals for sigma_z 
        mblockz=mblockz+1
        mcombcart(1,m1,m2,m3,m4)=3
        mcombcart(2,m1,m2,m3,m4)=mblockz                               
        call tosigZ(m1,m2,m3,m4,WRK(iangSO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartSO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breit)) 
     *  call tosigZ(m1,m2,m3,m4,WRK(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartOO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.m1.eq.-m3.and.l1.ne.0) then 
             call two2mean13(WRK(icartSO),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartz(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.m3.eq.-m4) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartOO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.m3.eq.-m4) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.m1.eq.-m2) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.and.m1.eq.-m2) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartSO),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             endif 
             endif 
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif 
        endif    
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        endif ! for check of significance for meanfield.  
        enddo 
        enddo 
        enddo 
        enddo 
        numbcart=numbcart+(mblockx+mblocky+mblockz)*ncont
cbs   just controlling if x and y integrals have the same number of blocks 
      if (mblockx.ne.mblocky) then 
      write(LUPRI,*) 
     *'numbers of integrals for sigma_x and sigma_y not equal!'    
      write(LUPRI,'(A12,4I3,2(A3,I5))') 
     *'l1,l2,l3,l4 ',l1,l2,l3,l4,' X:',mblockx,' Y:',mblocky  
      write(LUPRI,*) ' check the ipowxyz-array'
      CALL QUIT('Problems with IPOWXYA array in AMFI')
      endif   
cbs   start adresses for the next <ll|ll> block of integrals 
      endif
      endif
      endif
      endif
      enddo
      enddo
      enddo
      enddo
      return 
      end  
      subroutine buildcoul(l1,l2,l3,l4,! angular momenta of primitives
     *incl1,incl3, ! shifts for different radial integrals
     *Lrun, ! L-value for coulomb integrals 
     *prmints,
     *nprim1,nprim2,nprim3,nprim4,  ! number of primitives
     *expo1,expo2,expo3,expo4, ! arrays with the exponents
     *power13,
     *power24,
     *quotpow1,quotpow2
     *)
cbs ##################################################################
c
cbs  purpose: builds up the coulomb integrals 
cbs  inbetween primitives and multiplies 
cbs  with extra factors to correct the 
cbs  normalization                  
c
cbs ##################################################################
#include "implicit.h"
#include "para.h"
#include "param.h"
#include "dofuc.h"
#include "pi.h"
      dimension expo1(nprim1),
     *expo2(nprim2),
     *expo3(nprim3),
     *expo4(nprim4), ! the exponents
     *prmints(nprim1,nprim2,nprim3,nprim4), ! scratch array for integrals over primitives
     *power13(MxprimL,MxprimL),
     *power24(MxprimL,MxprimL),
     *quotpow1(nprim1,nprim2,nprim3,nprim4),
     *quotpow2(nprim1,nprim2,nprim3,nprim4),
     *fraclist1(0:Lmax+3),fraclist2(0:Lmax+3),fact(MxprimL),
     *frac(MxprimL),cfunctx1(MxprimL),cfunctx2(MxprimL)
      root8ovpi=dsqrt(8d0/pi)    
cbs ##################################################################
cbs        prepare indices for coulint
cbs ##################################################################
      n1=l1+incl1+1
      n2=l2+1
      n3=l3+incl3+1
      n4=l4+1
      n13=n1+n3
      n24=n2+n4
      index1=N13-Lrun-1
      index2=n24+Lrun
      index3=N24-Lrun-1
      index4=n13+Lrun
      do krun=0,(index1-1)/2
      fraclist1(krun)=dffrac(krun+krun+index2-1,krun+krun)*
     *dffrac(1,index2-1)
      enddo
      do krun=0,(index3-1)/2
      fraclist2(krun)=dffrac(krun+krun+index4-1,krun+krun)*
     *dffrac(1,index4-1)
      enddo
cbs ##################################################################
cbs   common factors including double factorials 
cbs ##################################################################
      doff1=dffrac(index1-1,n13-1)*dffrac(n24+Lrun-1,n24-1)
      doff2=dffrac(index3-1,n24-1)*dffrac(n13+Lrun-1,n13-1)
      if (index1.eq.1) then 
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
   	        pow24inv=doff1/power24(irun4,irun2)
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                prmints(irun1,irun2,irun3,irun4)=    
     *          quotpow1(irun1,irun2,irun3,irun4)*          
     *          dsqrt(0.5d0*(expo1(irun1)+expo3(irun3)))*
     *          power13(irun3,irun1)*pow24inv
      		enddo 
      		enddo 
      		enddo 
      		enddo 
      else 
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
                alpha24inv=1d0/(expo2(irun2)+expo4(irun4))
                pow24inv=doff1/power24(irun4,irun2)
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                a1324= alpha24inv*(expo1(irun1)+expo3(irun3))
                   Cfunctx1(irun1)=fraclist1(0)
                   frac(irun1)=a1324/(1d0+a1324)
                   fact(irun1)=frac(irun1) 
                enddo
*vocl    loop,repeat(Lmax+3)
                   do k=1,(index1-1)/2
                   do irun1=1,limit1  
                   Cfunctx1(irun1)=Cfunctx1(irun1)+fraclist1(k)
     *            *fact(irun1)
                   enddo
                   do irun1=1,limit1  
                   fact(irun1)=fact(irun1)*frac(irun1)
                   enddo
                   enddo
                do irun1=1,limit1  
                alpha13=0.5d0*(expo1(irun1)+expo3(irun3))
                prmints(irun1,irun2,irun3,irun4)=    
     *          quotpow1(irun1,irun2,irun3,irun4)*       
     *          dsqrt(alpha13)*power13(irun3,irun1)*pow24inv*
     *          Cfunctx1(irun1) 
		enddo
		enddo
		enddo
		enddo
      endif   
      if (index3.eq.1) then 
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
		pow24=doff2*power24(irun4,irun2)*
     *          dsqrt(0.5d0*(expo2(irun2)+expo4(irun4)))
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                prmints(irun1,irun2,irun3,irun4)=    
     *          prmints(irun1,irun2,irun3,irun4)+    
     *          pow24*quotpow2(irun1,irun2,irun3,irun4)/
     *          power13(irun3,irun1)
                enddo 
                enddo 
                enddo 
                enddo 
      else
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
                alpha24=expo2(irun2)+expo4(irun4)
		pow24=doff2*power24(irun4,irun2)*
     *          dsqrt(0.5d0*alpha24)                                  
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                a2413= alpha24/(expo1(irun1)+expo3(irun3))
                   Cfunctx2(irun1)=fraclist2(0)
                   frac(irun1)=a2413/(1d0+a2413)
                   fact(irun1)=frac(irun1)
        	enddo 
*vocl    loop,repeat(Lmax+3)
                   do k=1,(index3-1)/2
                   do irun1=1,limit1 
                   Cfunctx2(irun1)=Cfunctx2(irun1)+
     *             fraclist2(k)*fact(irun1)
                   enddo
                   do irun1=1,limit1 
                   fact(irun1)=fact(irun1)*frac(irun1)
                   enddo
                   enddo
                do irun1=1,limit1 
                prmints(irun1,irun2,irun3,irun4)=    
     *          prmints(irun1,irun2,irun3,irun4)+    
     *          quotpow2(irun1,irun2,irun3,irun4)*
     *          Cfunctx2(irun1)*
     *          pow24/power13(irun3,irun1)
                enddo 
                enddo 
                enddo 
                enddo 
      endif 
cbs   make some mirroring for identical l-values
cbs   for the case that l1=l3 
      if (l1.eq.l3) then 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=irun3+1,nprim1
      prmints(irun1,irun2,irun3,irun4)=    
     *prmints(irun3,irun2,irun1,irun4)
      enddo 
      enddo 
      enddo 
      enddo 
      endif  
cbs   for the case that l2=l4 
      if (l2.eq.l4) then 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=irun4+1,nprim2
      do irun1=1,nprim1
      prmints(irun1,irun2,irun3,irun4)=    
     *prmints(irun1,irun4,irun3,irun2)
      enddo 
      enddo 
      enddo 
      enddo 
      endif 
cbs   some factors which are the same for all cases 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      prmints(irun1,irun2,irun3,irun4)=    
     *prmints(irun1,irun2,irun3,irun4)*
     *coulovlp(irun4,irun2,0,0,l4,l2)*
     *coulovlp(irun3,irun1,incl3,incl1,l3,l1)*
     *root8ovpi
      enddo 
      enddo 
      enddo 
      enddo 
cbs   
cbs  look for additional factors, as the 
cbs  coulomb integrals are calculated 
cbs  for normalized functions with that 
cbs  specific l 
cbs  
cbs  if l was increased by one, the factor is
cbs  0.5*dsqrt((2l+3)/(exponent))
cbs  if l was decreased by one, the factor is
cbs  2d0*dsqrt(exponent/(2l+1))
cbs  
cbs
cbs   check for first function 
cbs
cbs  
      if (incl1.eq.1) then 
      fact1=0.5d0*dsqrt(dfloat(l1+l1+3))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1/dsqrt(expo1(irun1))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      elseif (incl1.eq.-1) then 
      fact1=2d0/dsqrt(dfloat(l1+l1+1))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1*dsqrt(expo1(irun1))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      endif
cbs  
cbs
cbs   check for third function 
cbs
cbs  
      if (incl3.eq.1) then 
      fact1=0.5d0*dsqrt(dfloat(l3+l3+3))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1/dsqrt(expo3(irun3))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      elseif (incl3.eq.-1) then 
      fact1=2d0/dsqrt(dfloat(l3+l3+1))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1*dsqrt(expo3(irun3))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      endif
      return 
      end 
      subroutine cartoneX(L,Lmax,onecontr,ncontrac,
     *MxcontL,onecartX)
#include "implicit.h"
      dimension onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *onecartX(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))  
cbs   arranges the cartesian one-elctron-integrals for X  on a 
cbs   quadratic matrix 
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   - + Integrals    m || mprime     mprime=m+1
      do Mprime=2,L
      M=mprime-1 
      iaddr=ipnt(Mprime+L+1,-M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartX(icont,jcont,iaddr)=
     *onecartX(icont,jcont,iaddr)
     *-0.25d0*(
     *onecontr(icont,jcont,Mprime,1)+
     *onecontr(icont,jcont,-Mprime,3))
      enddo
      enddo 
      enddo 
cbs   - + Integrals    m || mprime     mprime=m-1
      do Mprime=1,L-1
      M=mprime+1 
      iaddr=ipnt(Mprime+L+1,-M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartX(icont,jcont,iaddr)=
     *onecartX(icont,jcont,iaddr)
     *-0.25d0*(
     *onecontr(icont,jcont,Mprime,3)+
     *onecontr(icont,jcont,-Mprime,1))
      enddo
      enddo 
      enddo 
cbs   -1 || 0 integrals 
      pre=dsqrt(0.125d0) 
      iaddr=ipnt(L,L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartX(icont,jcont,iaddr)=
     *onecartX(icont,jcont,iaddr)
     *-pre* (onecontr(icont,jcont,0,3)+
     *onecontr(icont,jcont,0,1) )
      enddo
      enddo 
      return
      end 

      subroutine cartoneY(L,Lmax,onecontr,ncontrac,
     *MxcontL,onecartY)
#include "implicit.h"
      dimension onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *onecartY(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))  
cbs   arranges the cartesian one-electron integrals for Y  
cbs   on a quadratic matrix 
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   + + Integrals    m || mprime     mprime=m+1
      do Mprime=2,L
      M=mprime-1 
      iaddr=ipnt(Mprime+L+1,M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartY(icont,jcont,iaddr)=
     *onecartY(icont,jcont,iaddr) 
     *-0.25d0*(
     *onecontr(icont,jcont,Mprime,1)+
     *onecontr(icont,jcont,-Mprime,3))
      enddo
      enddo 
      enddo 
cbs   - - Integrals    m || mprime     mprime=m-1
      do Mprime=1,L-1
      M=mprime+1 
      iaddr=ipnt(-Mprime+L+1,-M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartY(icont,jcont,iaddr)=
     *onecartY(icont,jcont,iaddr) 
     *+0.25d0*(
     *onecontr(icont,jcont,Mprime,3)+
     *onecontr(icont,jcont,-Mprime,1))
      enddo
      enddo 
      enddo 
cbs    0 || 1 integrals 
      pre=-dsqrt(0.125d0)      
      iaddr=ipnt(L+1,L+2)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartY(icont,jcont,iaddr)=
     *onecartY(icont,jcont,iaddr) 
     *+pre*    
     *(onecontr(icont,jcont,1,1)+
     *onecontr(icont,jcont,-1,3)) 
      enddo 
      enddo 
      return
      end 
      subroutine cartoneZ(L,Lmax,onecontr,ncontrac,
     *MxcontL,onecartZ)
#include "implicit.h"
      dimension onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))  
cbs   arranges the cartesian one-electron integrals for Z  
cbs   on a quadratic matrix 
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   - + Integrals    m || mprime     mprime=m
      do Mprime=1,L
      iaddr=ipnt(Mprime+L+1,-mprime+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartZ(icont,jcont,iaddr)=
     *onecartZ(icont,jcont,iaddr)+
     *0.5d0*(
     *onecontr(icont,jcont,Mprime,2)-
     *onecontr(icont,jcont,-Mprime,2))
      enddo
      enddo 
      enddo 
      return
      end

 
      subroutine chngcont(coeffs,coeffst1,coeffst1a,coeffst2,
     *coeffst2a,ncont,nprims,evec,
     *type1,type2,work,work2,work3,MxprimL,
     *rootOVLP,OVLPinv,exponents)              
c###############################################################################
cbs   purpose: makes out of old contraction coefficients(in normalized functions)
cbs   new coefficients including the kinematical factors
cbs   using the diagonal matrices on type1 and type2 (see subroutine kinemat)
cbs   coeffst1a and coeffst2a additionally include the exponents alpha 
cbs   (that is why ....a). So the exponents in the integrals are moved 
cbs   to the contraction coefficients and not in some way into the primitive 
cbs   integrals. 
cbs
cbs   the different cases for contracted integrals differ later on in the
cbs   choice of different sets of contraction coefficients. 
cbs
c###############################################################################
#include "implicit.h"
      dimension coeffs(nprims,ncont),    ! original contraction coefficients   
     *coeffst1(nprims,ncont),            ! A * contraction coefficients
     *coeffst1a(nprims,ncont),           ! A * alpha*contraction coefficients
     *coeffst2a(nprims,ncont),           ! c*A/(E+m) * contraction coefficients
     *coeffst2(nprims,ncont),            ! c*A/(E+m) * alpha *contraction coefficients    
     *evec(nprims,nprims),
     *work(nprims,nprims) ,
     *work2(nprims,nprims) ,
     *work3(nprims,nprims) ,
     *rootOVLP(MxprimL,*),
     *OVLPinv(MxprimL,*),
     *type1(*),type2(*),
     *exponents(*) 
cbs   
cbs   first new coefficients for type1 (A) 
cbs   generate a transformation matrix on work
cbs   
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      work2(I,J)=0d0
      work3(I,J)=0d0
      enddo
      enddo
cbs   build up the transformation matrix 
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(I,J)=work(I,J)+evec(I,K)*type1(K)*evec(J,K)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work2(I,J)=work2(I,J)+work(I,K)*rootOVLP(K,J)    
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work3(I,J)=work3(I,J)+rootOVLP(I,K)*work2(K,J)    
      enddo
      enddo
      enddo
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(J,I)=work(J,I)+OVLPinv(I,K)*work3(K,J)       
      enddo
      enddo
      enddo
      do K=1,ncont
      do I=1,nprims
      coeffst1(I,K)=0d0
      enddo
      enddo
cbs   now transform the vectors 
      do K=1,ncont
      do J=1,nprims
      do I=1,nprims
      coeffst1(I,K)=coeffst1(I,K)+work(J,I)*coeffs(J,K)
      enddo 
      enddo 
      enddo 
cbs  
cbs   now with exponent   
cbs    
      do K=1,ncont
      do I=1,nprims
      coeffst1a(I,K)=exponents(I)*coeffst1(I,K) 
      enddo
      enddo
cbs   
cbs   and now the same for the other type  A/(E+m) 
cbs   
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      work2(I,J)=0d0
      work3(I,J)=0d0
      enddo
      enddo
cbs   build up the transformation matrix 
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(I,J)=work(I,J)+evec(I,K)*type2(K)*evec(J,K)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work2(I,J)=work2(I,J)+work(I,K)*rootOVLP(K,J)    
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work3(I,J)=work3(I,J)+rootOVLP(I,K)*work2(K,J)    
      enddo
      enddo
      enddo
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(J,I)=work(J,I)+OVLPinv(I,K)*work3(K,J)       
      enddo
      enddo
      enddo
      do K=1,ncont
      do I=1,nprims
      coeffst2(I,K)=0d0
      enddo
      enddo
cbs   now transform the vectors 
      do K=1,ncont
      do J=1,nprims
      do I=1,nprims
      coeffst2(I,K)=coeffst2(I,K)+work(J,I)*coeffs(J,K)
      enddo 
      enddo 
      enddo 
cbs  
cbs   now with exponent   
cbs    
      do K=1,ncont
      do I=1,nprims
      coeffst2a(I,K)=exponents(I)*coeffst2(I,K) 
      enddo
      enddo
      return 
      end

 
      subroutine cont(L,breit,ifinite)
cbs###########################################################################
cbs   cont prepares all required contraction coefficients for functions 
cbs   with angular momentum L
cbs###########################################################################
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension tkintria((MxprimL*MxprimL+MxprimL)/2)
      logical breit,breit_finite
      breit_finite=.true.
cbs   transcon transfers and normalizes contracted functions
cbs   ore more precizely the coefficients  
CMI ... modified according to 'amfi_wizard'  
      call transcon(cntscrtch(1,1,L),MxprimL,
     *MxcontL,normovlp(1,1,L),
     *contrarray(iaddori(L)),nprimit(L),ncontrac(L),ovlpcont(1,1,L))
cbs   gentkin generates the matrix of kinetic energy TKIN 
      call gentkin(L,TKIN,nprimit(L),exponents(1,L),rootOVLPinv(1,1,L))
cbs   kindiag diagonalizes TKIN 
cbs   for finite nucleus 
      if (ifinite.eq.2.and.L.eq.0) then  
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit_finite)  
      else 
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit)
      endif 
cbs   kinemat generates kinematic factors in 
cbs   the basis of eigenvectors   
      call kinemat(L,nprimit(L),eval,type1,type2,Energy)
      incr=nprimit(L)*ncontrac(L)
cbs   chngcont= changecont generates the contraction coeffs
cbs   including kinematic factors and even exponents as factors
      call chngcont(
     *contrarray(iaddori(L)),                  
     *contrarray(iaddtyp1(L)),                  
     *contrarray(iaddtyp2(L)),                  
     *contrarray(iaddtyp3(L)),                  
     *contrarray(iaddtyp4(L)),                  
     *ncontrac(L),nprimit(L),evec,
     *type1,type2,scratch4,scratch4(nprimit(L)*nprimit(L)+1),
     *scratch4(2*nprimit(L)*nprimit(L)+1),MxprimL,
     *rootOVLP(1,1,L),OVLPinv(1,1,L),
     *exponents(1,L))
      return
      end
