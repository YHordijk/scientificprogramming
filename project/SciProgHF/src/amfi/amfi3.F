      Subroutine contandmult(Lhigh,makemean,AIMP,oneonly,
     &  ADDSO1,numballcart,
     *  LUPROP,ifinite,WRK,LWRK)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "param.h"
#include "ired.h"
#include "alphac.h"
      logical ADDSO1
      logical makemean,AIMP,oneonly  
      character*8 xa,ya,za
      dimension xa(4),ya(4),za(4)
      DIMENSION WRK(LWRK)
      common /nucleus/ charge,Exp_Finite
      double precision normasHERMIT(-Lmax:Lmax,0:Lmax)
Chjaaj jan07:  sss_normasHERMIT is a "trick" so the data
C     statement later also works for gfortran compiler;
C     it would not accept the variable loop limits in the ml loop.
      double precision sss_normasHERMIT(25)
      equivalence (normasHERMIT, sss_normasHERMIT)
CMI   ... addition 
      integer function indx_LML

      data sss_normasHERMIT
     &    /1.0d0, 1.0d0,1.0d0,1.0d0,
     &     1.0d0,1.0d0,3.46410162d0,1.0d0,2.0d0,
     &     4.8989795d0,1.0d0,6.3245553d0,-2.5819889d0,6.3245553d0,
     &     2.0d0,-1.6329932d0,
     &   3.4641016d0,4.89897949d0,9.16515139d0,4.3204938d0,
     &  -3.4156503d0,4.3204938d0,18.330303d0,-1.6329932d0,-6.9282032d0/ 
      IPNT(I,J)=(J*J-J)/2+I 
cbs   get back the real number of functions for the finite nucleus   
      if (ifinite.eq.2) ncontrac(0)=ncontrac_keep
c###############################################################################
cbs   subroutine to contract radial one-electron integrals 
cbs   and multiply them with angular factors 
c###############################################################################
      xa(1)='********'
      ya(1)='********'
      za(1)='********'
      xa(2)='        '
      ya(2)='        '
      Za(2)='        '
      xa(3)='ANTISYMM'
      ya(3)='ANTISYMM'
      Za(3)='ANTISYMM'
      xa(4)='X1MNF-SO'
      ya(4)='Y1MNF-SO'
      ZA(4)='Z1MNF-SO'
c
cbs   clean the arrays for cartesian integrals
C
      length3=(numbalLcart*numbalLcart+numbalLcart)/2
      iloca=length3
      IOCAX = 1
      iocay=iocax+iloca
      iocaz=iocay+iloca
      iocax2=iocaz+iloca
      iocay2=iocax2+iloca
      iocaz2=iocay2+iloca

      KLAST = IOCAZ2 + ILOCA
CMI   IF (KLAST .GT. LWRK) CALL STOPIT('AMFI  ','CAMUL',KLAST,LFREE)
      IF (KLAST .GT. LWRK) CALL QUIT('AMFI CAMUL')
      call dzero(WRK(iocax),6*length3)
c
c
c
c
cbs   one-electron-integrals:
cbs   1. index: number of first contracted function
cbs   2. index: number of second contracted function
cbs   3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
cbs   4. index: L-value
cbs    onecartX(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
cbs    onecartY(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
cbs    onecartZ(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax)
c
C===================================================================
cbs   generate one-electron integrals for all L greater/equal 1
C===================================================================
      if (ifinite.eq.2) charge=0d0 ! nuclear integrals are modelled for finite nucleus somewhere else

CMI ... based on BS advise to remove 1 el.SO integrals..
CMI ... this is only temporary solution for now...
       IF (.NOT.ADDSO1) charge=0d0
  
      do L=1,Lhigh   
      	call contone(L,oneoverr3(1,L),onecontr(1,1,-Lmax,1,L),
     *  Lmax,contrarray(iaddtyp3(L)),nprimit(L),ncontrac(L),
     *	MxcontL,dummyone,
     *	onecartx(1,1,1,L),onecartY(1,1,1,L),onecartZ(1,1,1,L),
     *	charge,oneonly)
      Enddo 
c
cbs   ***********************************************************************
cbs   now move all integrals to one big arrays for X,Y,Z
cbs   ***********************************************************************
      do Lrun=1,Lhigh  !loop over L-values (integrals are diagonal in L)
      mrun=0
      do Msec=-Lrun,Lrun    ! cartesian M-values  (Mfirst,Msec) with 
      do Mfirst=-Lrun,Msec  ! Mfirst <= Msec (actually '=' does never appear
c                             as there is no L-component  in Ag
C
c
cbs   determine  if L_X L_Y or L_Z 
        ipowx=ipowxyz(1,mfirst,Lrun)+ipowxyz(1,msec,Lrun)
        ipowy=ipowxyz(2,mfirst,Lrun)+ipowxyz(2,msec,Lrun)
        ipowz=ipowxyz(3,mfirst,Lrun)+ipowxyz(3,msec,Lrun)
c
      	mrun=mrun+1
cbs     now determine the irreducable representations 
      	iredfirst=iredLM(Mfirst,Lrun)
      	iredsec=iredLM(Msec,Lrun)
cbs     check out which IR is the lower one. 
      	if (iredfirst.le.iredsec) then 
cbs     calculate shift to get to the beginning of the block 
      	   iredired=shiftIRIR((iredsec*iredsec-iredsec)/2+iredfirst)
     *	   +incrlm(Mfirst,Lrun)*itotalperIR(iredsec)+
     *	   incrLM(Msec,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then  
      		do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
      		do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
CBS                print *, 'iocax',iocax,iredired,icartsec
      		WRK(iocax+iredired+(icartsec-1))=
     *		WRK(iocax+iredired+(icartsec-1)) 
     *		+onecartx(icartfirst,icartsec,mrun,Lrun)
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredsec)
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then  
      		do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
      		do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
      		WRK(iocay+iredired+(icartsec-1))=
     *		WRK(iocay+iredired+(icartsec-1)) 
     *	        +onecarty(icartfirst,icartsec,mrun,Lrun)
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredsec)
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then  
      		do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
      		do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
      		WRK(iocaz+iredired+(icartsec-1))=
     *		WRK(iocaz+iredired+(icartsec-1)) 
     *	        +onecartz(icartfirst,icartsec,mrun,Lrun)
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredsec)
      		enddo
        endif 
      	elseif (iredfirst.gt.iredsec) then 
cbs     In this case, indices are exchanged with respect to former 
cbs     symmetry of blocks. Therefore, there will be a minus sign 
c
cbs     calculate shift to get to the beginning of the block 
      		iredired=shiftIRIR((iredfirst*iredfirst-iredfirst)/2+
     *		iredsec)+
     *		incrLM(Msec,Lrun)*itotalperIR(iredfirst)+  
     *		incrLM(Mfirst,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then  
      		do icartsec=1,ncontrac(Lrun) !loop over functions second index
      		do icartfirst=1,ncontrac(Lrun) !loop over functions first index
      		WRK(iocax+iredired+(icartfirst-1))=
     *		WRK(iocax+iredired+(icartfirst-1)) 
     *         -onecartx(icartsec,icartfirst,mrun,Lrun)
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredfirst)
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then  
      		do icartsec=1,ncontrac(Lrun) !loop over functions second index
      		do icartfirst=1,ncontrac(Lrun) !loop over functions first index
      		WRK(iocay+iredired+(icartfirst-1))=
     *		WRK(iocay+iredired+(icartfirst-1)) 
     *	       -onecarty(icartsec,icartfirst,mrun,Lrun)
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredfirst)
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then  
      		do icartsec=1,ncontrac(Lrun) !loop over functions second index
      		do icartfirst=1,ncontrac(Lrun) !loop over functions first index
      		WRK(iocaz+iredired+(icartfirst-1))=
     *		WRK(iocaz+iredired+(icartfirst-1)) 
     *	       -onecartz(icartsec,icartfirst,mrun,Lrun)
      		enddo
      		iredired=iredired+itotalperIR(iredfirst)
      		enddo
        endif 
      endif 
      enddo
      enddo
      enddo
C
C
cbs   copy integrals on arrays with no symmetry blocking at all 
cbs   which means huge triangular matrices 
      irun=0
      do norb2=1,numballcarT  
      ired2=iredoffunctnew(norb2)
      norbsh2=norb2-shiftIRED(ired2) 
      do norb1=1,norb2    
      ired1=iredoffunctnew(norb1)
      norbsh1=noRb1-shiftIRED(ired1) 
      irun=irun+1 
      iredirEd=shiftIRIR((ired2*ired2-ired2)/2+
     *		ired1)
      if (ired1.ne.ired2) then 
        WRK(iocax2+irun-1)=WRK(iocax-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2)) 
        WRK(iocay2+irun-1)=WRK(iocay-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2)) 
        WRK(iocaz2+irun-1)=WRK(iocaz-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2)) 
      else 
       WRK(iocax2+irun-1)=WRK(iocax-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)                 
       WRK(iocay2+irun-1)=WRK(iocay-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)                 
       WRK(iocaz2+irun-1)=WRK(iocaz-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)                 
      endif 
      Enddo 
      enddo

CMI   ... prepare the factor 1/4c^2
      fac = 0.25D0*ALPHA2 

CMI ... here it is : Waw !
c     write a hermit-like file   b.s. 4.10.96   
CBS   write(6,*) 'number of orbitals ',numbalLcarT
CBS   write(6,*) 'length of triangular matrix ', length3

CMI           write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
CMI  *        i=1,numbofsym),
CMI  *        numballcart,(Loffunction(I),I=1,numballcart), 
CMI  *        (Moffunction(I),I=1,numballcart),
CMI  *        Lhigh,(ncontrac(I),I=0,Lhigh)          

              write(LUPROP)  xa
              write(LUPROP)  numballcart
C        ... write only type 5HERLIT) and index 
          do I=1,numballcart

           IXL  = Loffunction(I)
           IXML = Moffunction(I)
           IXINLML = indx_LML(Loffunction,Moffunction,I)

           IPTYP=0
        IF (IXL.EQ.0.AND.IXML.EQ.0) THEN
            IPTYP=1
        ELSE IF (IXL.EQ.1.AND.IXML.EQ.1) THEN
            IPTYP=2 ! 2px
        ELSE IF (IXL.EQ.1.AND.IXML.EQ.-1) THEN
            IPTYP=3 ! 2py
        ELSE IF (IXL.EQ.1.AND.IXML.EQ. 0) THEN
            IPTYP=4 ! 2pz
        ELSE IF (IXL.EQ.2.AND.IXML.EQ.-2) THEN
            IPTYP=5 ! 3d2-
        ELSE IF (IXL.EQ.2.AND.IXML.EQ.-1) THEN
            IPTYP=6 ! 3d1-
        ELSE IF (IXL.EQ.2.AND.IXML.EQ. 0) THEN
            IPTYP=7 ! 3d0
        ELSE IF (IXL.EQ.2.AND.IXML.EQ.+1) THEN
            IPTYP=8 ! 3d1+
        ELSE IF (IXL.EQ.2.AND.IXML.EQ.+2) THEN
            IPTYP=9 ! 3d2+
        ELSE IF (IXL.EQ.3.AND.IXML.EQ.-3) THEN
            IPTYP=10 ! 4f3-
        ELSE IF (IXL.EQ.3.AND.IXML.EQ.-2) THEN
            IPTYP=11 ! 4f2-
        ELSE IF (IXL.EQ.3.AND.IXML.EQ.-1) THEN
            IPTYP=12 ! 4f1-
        ELSE IF (IXL.EQ.3.AND.IXML.EQ. 0) THEN
            IPTYP=13 ! 4f0
        ELSE IF (IXL.EQ.3.AND.IXML.EQ.+1) THEN
            IPTYP=14 ! 4f1+
        ELSE IF (IXL.EQ.3.AND.IXML.EQ.+2) THEN
            IPTYP=15 ! 4f2+
        ELSE IF (IXL.EQ.3.AND.IXML.EQ.+3) THEN
            IPTYP=16 ! 4f3+
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.-4) THEN
            IPTYP=17 ! 5g4-
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.-3) THEN
            IPTYP=18 ! 5g3-
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.-2) THEN
            IPTYP=19 ! 5g2-
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.-1) THEN
            IPTYP=20 ! 5g1-
        ELSE IF (IXL.EQ.4.AND.IXML.EQ. 0) THEN
            IPTYP=21 ! 5g0
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.+1) THEN
            IPTYP=22 ! 5g1+
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.+2) THEN
            IPTYP=23 ! 5g2+
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.+3) THEN
            IPTYP=24 ! 5g3+
        ELSE IF (IXL.EQ.4.AND.IXML.EQ.+4) THEN
            IPTYP=25 ! 5g4+
        ELSE
         write(lupri,*) 'contandmult: I/numballcart/Mxcart=',I,'/',
     &   numballcart,'/',Mxcart,
     &   'IXL=',IXL,' IXML=',IXML,' IXINLML=',IXINLML 
         CALL QUIT('AMFI:IPTYP not defined or another error !!!')
        ENDIF
           write(LUPROP) IPTYP,IXINLML
          enddo
              write(LUPROP) Lhigh,(ncontrac(I),I=0,Lhigh)
CMO           write(LUPROP)  xa
              write(LUPROP) (WRK(iocax2+irun)*fac,irun=0,length3-1) 
              write(LUPROP)  Ya  
              write(LUPROP) (WRK(iocay2+irun)*fac,irun=0,length3-1) 
CMI           .... I don't know why...
C             write(LUPROP)  Za    
              write(LUPROP) (WRK(iocaz2+irun)*fac,irun=0,length3-1) 

C            write(lupri,*) (WRK(iocax2+irun)*fac,irun=0,length3-1)
cbs   
cbs   that is it!!
cbs   
      return
      end


      subroutine contcasaOO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (+2)  with alpha1*alpha3
cbs   other-orbit term 
cbs   use averaged integrals by interchanging kinematic factors   
cbs   this is case a in the documentation 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)), !A *alpha                                       
     *contrarray(iaddtyp3(l2)), !A/E+m                                         
     *contrarray(iaddtyp4(l3)), !A/E+m *alpha     
     *contrarray(iaddtyp1(l4)), !A                                       
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)), 
     *contrarray(iaddtyp3(l2)), 
     *contrarray(iaddtyp2(l3)), 
     *contrarray(iaddtyp1(l4)), 
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp4(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)), 
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp2(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
      return 
      end
      subroutine contcasaSO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (+2)  with alpha1*alpha3
cbs   same orbit term 
cbs   this is case a in the documentation 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
c     write(6,*) 'scratch1 ',(scratch1(I),I=1,ilength)          
c     write(6,*) 'contraction coeff'
c     write(6,*) (contrarray(iaddtyp4(l1)+I),I=0,nprim(1)-1)
c     write(6,*) (contrarray(iaddtyp1(l2)+I),I=0,nprim(2)-1)
c     write(6,*) (contrarray(iaddtyp4(l3)+I),I=0,nprim(3)-1)
c     write(6,*) (contrarray(iaddtyp1(l4)+I),I=0,nprim(4)-1)
      call contract(
     *contrarray(iaddtyp4(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
c     write(6,*) 'nstart ',nstart 
      do irun=1,nprod
      cont4SO(nstart+irun-1)=scratch1(irun)
      enddo 
      return 
      end
      subroutine contcasb1OO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (0)  with alpha1
cbs   this is one of the cases b in the documentation 
cbs   use averaged integrals by interchanging kinematic factors 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
C
C
c
cbs   copy primitive integrals to scratch1
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),    
     *contrarray(iaddtyp3(l2)),  
     *contrarray(iaddtyp3(l3)),  
     *contrarray(iaddtyp1(l4)), 
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)),
     *contrarray(iaddtyp3(l2)),
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return 
      end


      subroutine contcasb1SO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (0)  with alpha1
cbs   this is one of the cases b in the documentation 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
cbs   copy primitive integrals to scratch1
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      return 
      end
      subroutine contcasb2OO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (0)  with alpha3
cbs   this is one of the cases b in the documentation 
cbs   use averaged integrals by interchanging kinematic factors 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp2(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp2(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return 
      end


      subroutine contcasb2SO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (0)  with alpha3
cbs   this is one of the cases b in the documentation 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      return 
      end


      SUBroutine contcascOO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (-2)  with factor 1 
cbs   this is case c in the documentation 
cbs   use averaged integrals by interchanging kinematic factors 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return                           
      end
      subroutine contcascSO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (-2)  with factor 1 
cbs   this is case c in the documentation 
#include "implicit.h"
#include "para.h"
#include "param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      return 
      end


      subroutine contone(L,oneoverr3,onecontr,Lmax,
     *contcoeff,nprim,ncont,MxcontL,dummy,
     *onecartx,onecartY,onecartZ,charge,oneonly)

cbs   contracts one-electron integrals and multiplies with l,m-dependent 
cbs   factors for L-,L0,L+ 

#include "implicit.h"
      dimension oneoverR3(*),
     *onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *contcoeff(nprim,ncont),dummy(ncont,ncont),
     *onecartx(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),  
     *onecarty(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),  
     *onecartz(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1))   
      logical oneonly  
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   first of all cleaning dummy and onecontr
      do jrun=1,ncont
      do irun=1,ncont
      dummy(irun,jrun)=0d0   
      enddo 
      enddo 

      if (oneonly) then
      iprod=MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)
      call dzero(onecartx,iprod)
      call dzero(onecarty,iprod)
      call dzero(onecartz,iprod)
      endif 

      iprod=3*(Lmax+lmax+1)*MxcontL*MxcontL
      call dzero(onecontr,iprod)
cbs   contract onto dummy 
      do icont2=1,ncont
      do icont1=1,ncont 
      do iprim2=1,nprim
      do iprim1=1,nprim
      dummy(icont1,icont2)=dummy(icont1,icont2)+
     *contcoeff(iprim1,icont1)*contcoeff(iprim2,icont2)*
     *oneoverR3(ipnt(iprim1,iprim2))
      enddo
      enddo
      enddo 
      enddo 
      do icont2=1,ncont
      do icont1=1,ncont 
      dummy(icont1,icont2)=dummy(icont1,icont2)*charge 
      enddo
      enddo
cbs   start to add l,m dependent factors    
      do M=-L,L 
      factormin=dsqrt(dfloat(L*L-M*M+L+M))
      factor0=dfloat(M)
      factorplus=dsqrt(dfloat(L*L-M*M+L-M))
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,1)=dummy(jrun,irun)*factormin  ! L-minus 
      enddo
      enddo
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,2)=dummy(jrun,irun)*factor0    ! L-0     
      enddo
      enddo
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,3)=dummy(jrun,irun)*factorplus ! L-plus  
      enddo
      enddo
      enddo

cbs   make the final cartesian integrals 
      call cartoneX(L,Lmax,onecontr,ncont,
     *MxcontL,onecartX(1,1,1))
      call cartoneY(L,Lmax,onecontr,ncont,
     *MxcontL,onecartY(1,1,1))
      call cartoneZ(L,Lmax,onecontr,ncont,
     *MxcontL,onecartZ(1,1,1))
      return 
      end


      subroutine contract(
     *coeffs1, !(nprim(1),ncont(1)) modified contraction coefficients
     *coeffs2, !(nprim(2),ncont(2)) modified contraction coefficients
     *coeffs3, !(nprim(3),ncont(3)) modified contraction coefficients
     *coeffs4, !(nprim(4),ncont(4)) modified contraction coefficients
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
cbs  array one contains at the beginning the uncontracted integrals 
     *arr1,  ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
     *arr2   ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
     *)
#include "implicit.h"
      dimension coeffs1(*),coeffs2(*),coeffs3(*),coeffs4(*),   
     *arr1(*),arr2(*),ncont(4),nprim(4),nolds(4),nnew(4)
C
cbs   makes four indextransformations in a row....
cbs   try to find out, which indices should be transformed first...
c
      ratio1=dfloat(nprim(1))/dfloat(ncont(1))
      ratio2=dfloat(nprim(2))/dfloat(ncont(2))
      ratio3=dfloat(nprim(3))/dfloat(ncont(3))
      ratio4=dfloat(nprim(4))/dfloat(ncont(4))
      do IBM=1,4
      nolds(IBM)=nprim(IBM)  
      nnew(IBM)=nprim(IBM)  
      enddo
cbs   determine first, second,third and last index
cbs   determine the first
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ifirst=1
      ratio1=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio2) then 
      ifirst=2
      ratio2=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio3) then 
      ifirst=3
      ratio3=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio4) then 
      ifirst=4
      ratio4=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      endif
      nolds(ifirst)=nnew(ifirst)
cbs   determine the second   
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      isec=1
      ratio1=0
      nnew(isec)=ncont(isec)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio2) then 
      isec=2
      ratio2=0
      nnew(isec)=ncont(isec)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio3) then 
      isec=3
      ratio3=0
      nnew(isec)=ncont(isec)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio4) then 
      isec=4
      ratio4=0
      nnew(isec)=ncont(isec)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      endif
      nolds(isec)=nnew(isec)
cbs   determine the third    
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ithird=1
      ratio1=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio2) then 
      ithird=2
      ratio2=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio3) then 
      ithird=3
      ratio3=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio4) then 
      ithird=4
      ratio4=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      endif
      nolds(ithird)=nnew(ithird)
cbs   determine the last    
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ifourth=1
      ratio1=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio2) then 
      ifourth=2
      ratio2=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio3) then 
      ifourth=3
      ratio3=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio4) then 
      ifourth=4
      ratio4=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      endif
cbs   contracted integrals are now on 
cbs   arr1(ncont1,ncont2,ncont3,ncont4)   
      return 
      end


      double precision function  couple3J(
     *l1,     ! integer  l1
     *l2,     ! integer  l2
     *l3,     ! integer  l3
     *m1,     ! integer  m1
     *m2,     ! integer  m2
     *m3)     ! integer  m3
cbs this routine calculates the coupling of three angular momenta to  zero   
cbs   
cbs   
cbs   Int dOmega i^(l1+l2+l3) Y^l1_m1 (Omega) Y^l2_m2 (Omega) Y^l3_m3 (Omega) =
cbs   sqrt( (2l1+1)(2l2+1)(2l2+3)/ 4Pi)  * 3J(l1,l2,l3,0,0,0) *
cbs   3J(l1,l2,l3,m1,m2,m3)
cbs   
cbs   
#include "implicit.h"
#include "pi.h"
      double precision inv4pi
cbs   (4*PI)**-1 
      inv4pi=0.25d0/pi                                   
cbs   initialize couple3J-coefficient
      couple3J=0d0
cbs   quick check 
      if (m1+m2+m3.ne.0) return 
cbs   double all values for regge3j   
      l1d=l1+l1
      l2d=l2+l2
      l3d=l3+l3
      m1d=m1+m1
      m2d=m2+m2
      m3d=m3+m3
      fac1=dsqrt(dfloat(l1d+1)*dfloat(l2d+1)*dfloat(l3d+1)*inv4pi)
      fac2=regge3j(l1d,l2d,l3d,0,0,0)
      fac3=regge3j(l1d,l2d,l3d,m1d,m2d,m3d)
      couple3J=fac1*fac2*fac3
      return 
      end 

      subroutine daxpint(from,to,fact,ndim1,ndim2,ndim3,ndim4)
#include "implicit.h"
cbs   subroutine similar to daxpy with interchange of two indices 
cbs   change from physicists notation to chemists notaion 
cbs   to(i,j,k,l)=to(i,j,k,l)+fact*from(i,k,j,l) 
      dimension from(ndim1,ndim2,ndim3,ndim4),
     *to(ndim1,ndim3,ndim2,ndim4)
      if (fact.eq.0d0) return 
      do irun4=1,ndim4
      do irun3=1,ndim3
      do irun2=1,ndim2
      do irun1=1,ndim1
      to(irun1,irun3,irun2,irun4)=to(irun1,irun3,irun2,irun4)+
     *fact*from(irun1,irun2,irun3,irun4)  
      enddo
      enddo
      enddo
      enddo
      return   
      end 


      subroutine gen1overR3(Lhigh)
#include "implicit.h"
cbs   generates the radial integrals  for the one electron spin orbit integrals   
cbs   taken the 1/r**3 formula from the documentation and included additional 
cbs   factors for normalization   
#include "para.h"
#include "param.h"
#include "dofuc.h"
#include "pi.h"
      do L=1,Lhigh    
      icount=0
      do iprim2=1,nprimit(L)
      alpha2=exponents(iprim2,L)
      do iprim1=1,iprim2    
      alpha1=exponents(iprim1,L)
      icount=icount+1
      oneoverR3(icount,L)=dsqrt(2d0/pi)*
     *(df(L+L-2)*2**(L+3)*
     *(alpha1*alpha2)**(0.25d0*
     *(L+L+3)))/((alpha1+alpha2)**L*df(L+L+1))
      enddo
      enddo
      enddo 
      return 
      end 
  


      subroutine gencoul(l1,l2,l3,l4,makemean,
     *bonn,breit,sameorb,cont4SO,cont4OO,icont4,
     *WRK,LFREE)
#include "implicit.h"      
cbs   SUBROUTINE to generate all required radial 
cbs   integrals for the four angular momenta l1-l4
#include "priunit.h"
#include "para.h"
#include "param.h"
      logical makemean,bonn,breit,sameorb 
      dimension cont4SO(*),cont4OO(*),WRK(LFREE)

#include "amfi_if.h"
c     logical amfi3_verbose
c     amfi3_verbose = .false.

      max1=1  !starting values for limits of precalculated 
c             ! powers of function Cfunct(X)
      max2=1
cbs   first of all, this routine determines, for which L 
cbs   values the radial integrals have to be solved
cbs   initialize the number of blocks for the different 
cbs   l-combinations  
cbs   no (ss|ss) contributions 
      if (l1.eq.0.and.l2.eq.0.and.l3.eq.0.and.l4.eq.0) return  ! no integrals for <ss|ss> 
      if (makemean) then 
      		nblock=1  ! sp sp are the first, so the first block
      		Lstarter(1)=1
      else 
      CALL QUIT('only mean-field with this version')
      endif 
cbs   keep track of L-values for later purposes   
      Lvalues(1)=l1
      Lvalues(2)=l2
      Lvalues(3)=l3
      Lvalues(4)=l4
cbs   now nanz is given the new value
      nanz=ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
      nprimprod=nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)
      IQUOT1 = 1
      iquot2=iquot1+nprimprod
      iquotp1=iquot2+nprimprod
      iquotp2=iquotp1+nprimprod
      iprim=iquotp2+nprimprod
      iscr1=iprim+nprimprod
      iscr2=iscr1+nprimprod
      KLAST = ISCR2 + NPRIMPROD
CMI   IF (KLAST .GT. LFREE) CALL STOPIT('AMFI  ','GENCOU',KLAST,LFREE)
      IF (KLAST .GT. LFREE) CALL QUIT('AMFI-GENCOU')
c
      call initfrac(nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),
     *WRK(iquot1),WRK(iquot2),exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4))
cbs   prepare the powers needed for cfunctx
c
c
c     There are seven different CASES of integrals following 
c       (   A  --  C) 
c
c     The structure is the same for all cases, therefore comments can be found only on case A
c
c
c
cbs   ###########################################################################################################
cbs   the (+2) cases          CASE A
cbs   ##########################################################################################################
      incl1=1  !  Those increments define the case 
      incl3=1
cbs   determine the possible L-values for the integrals by checking for triangular equation 
c
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
c
cbs   returns first and last L-values (Lanf,Lend), for which 
cbs   radial integrals have to be calculated 
      if(Lend-Lanf.ge.0) then 
cbs   if there are blocks 
      	Lblocks(1)=(Lend-Lanf)/2+1 ! L increases in steps of 2, 
cbs                                       due to parity conservation 
      	Lfirst(1)=Lanf
      	Llast(1)=Lend 
      else 
      	Lblocks(1)=0 
      endif 
      if (Lblocks(1).gt.0) then    ! integrals have to be calculated 
cbs### check, whether integrals fit on array ################
      if  (Lstarter(1)+nanz*Lblocks(1).gt.icont4) then 
         write(LUPRI,*) 'end at: ',Lstarter(1)+nanz*Lblocks(1) 
         CALL QUIT('increase icont4 in amfi.F')
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(1) ! gives the address, where to write the contracted integrals 
cbs   ipow1 and ipow2 are the the numbers of powers in the prefactor
cbs   of the function Cfunct 
cbs   now loop over possible L-values 
      do Lrun= Lfirst(1),Llast(1),2
      		ipow1=2+(l2+l4+Lrun)/2
      		ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
cbs   those powers have to be generated... 
      call getpow(ipow1,WRK(iquot1),WRK(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
cbs   those powers have to be generated... 
      call getpow(ipow2,WRK(iquot2),WRK(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
c     in buildcoul the radial integrals are calculated 
      		call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *		Lrun,WRK(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *		nprimit(l4), 
     *		exponents(1,l1),exponents(1,l2),
     *		exponents(1,l3),exponents(1,l4),
     *		powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *		WRK(iquotp1),WRK(iquotp2))
cbs   in the contcas_ routines the integrals are contracted, including exponents as prefactors... 
      		if (bonn.or.breit.or.sameorb) then      
      		call contcasASO(l1,l2,l3,l4,istart,WRK(iprim),
     *           WRK(iscr1),WRK(iscr2),cont4SO)
                else 
      		call contcasASO(l1,l2,l3,l4,istart,WRK(iprim),
     *           WRK(iscr1),WRK(iscr2),cont4SO)
                call contcasAOO(l1,l2,l3,l4,istart,WRK(iprim),
     *           WRK(iscr1),WRK(iscr2),cont4OO)
                endif 
      		istart=istart+nanz  ! start-address for the next block of contracted integrals 
      enddo 
      endif   
cbs   ##########################################################################################################
cbs   the (0) cases         CASE  B
cbs   ##########################################################################################################
      incl1=0
      incl3=0
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(2)=(Lend-Lanf)/2+1
      Lfirst(2)=Lanf
      Llast(2)=Lend 
      Lblocks(3)=(Lend-Lanf)/2+1
      Lfirst(3)=Lanf
      Llast(3)=Lend 
      else 
      Lblocks(2)=0 
      Lblocks(3)=0 
      endif 
      Lstarter(2)=Lstarter(1)+
     *nanz*Lblocks(1)
      Lstarter(3)=Lstarter(2)+
     *nanz*Lblocks(2)
cbs   primitive integrals are the same for type 2 and 3  !!!!!
      if (Lblocks(2).gt.0) then    
cbs### check, whether integrals fit on array ################
      if  (Lstarter(2)+2*nanz*Lblocks(2).gt.icont4) then 
      write(LUPRI,*) 'end at: ',Lstarter(2)+2*nanz*Lblocks(2) 
      CALL QUIT('increase icont4 in amfi.F')
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(2)
      istart2=Lstarter(3)
      do Lrun= Lfirst(2),Llast(2),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      call getpow(ipow1,WRK(iquot1),WRK(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call getpow(ipow2,WRK(iquot2),WRK(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,WRK(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4), 
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),   
     *WRK(iquotp1),WRK(iquotp2))
      if (bonn.or.breit.or.sameorb) then      
      call contcasB1SO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      else 
      call contcasB1SO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      Call contcasB1OO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4OO)
      Call contcasB2OO(l1,l2,l3,l4,istart2,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4OO)
      endif   
      istart=istart+nanz
      istart2=istart2+nanz
      enddo 
      endif   
cbs   ##########################################################################################################
cbs   the (-2) cases      CASE C
cbs   ##########################################################################################################
      if (l1.eq.0.or.l3.eq.0) then 
      Lblocks(4)=0
      else 
      incl1=-1
      incl3=-1
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(4)=(Lend-Lanf)/2+1
      Lfirst(4)=Lanf
      Llast(4)=Lend 
      else 
      Lblocks(4)=0 
      endif 
      endif
      Lstarter(4)=Lstarter(3)+
     *nanz*Lblocks(3)
      if (Lblocks(4).gt.0) then    
cbs### check, whether integrals fit on array ################
      if  (Lstarter(4)+nanz*Lblocks(4).gt.icont4) then 
      write(LUPRI,*) 'end at: ',Lstarter(4)+nanz*Lblocks(4) 
      CALL QUIT('increase icont4 in amfi.F')
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(4)
      do Lrun= Lfirst(4),Llast(4),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      call getpow(ipow1,WRK(iquot1),WRK(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call getpow(ipow2,WRK(iquot2),WRK(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,WRK(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4), 
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *WRK(iquotp1),WRK(iquotp2))
      if (bonn.or.breit.or.sameorb) then      
      call contcasCSO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      else 
      call contcasCSO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      call contcasCOO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4OO)
      endif 
      istart=istart+nanz
      enddo 
      endif   
      return 
      end  
      subroutine gencoulDIM(l1,l2,l3,l4,makemean,
     *bonn,breit,sameorb,icont4)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "param.h"
cbs   SUBROUTINE to calculate the dimemsion of the radial integral 
cbs   arrays. BASICALLY GENCOUL WITHOUT EXPLICIT INTEGRAL CALCULATION 
cbs   integrals for the four angular momenta l1-l4
      logical makemean,bonn,breit,sameorb 
      max1=1  !starting values for limits of precalculated 
c             ! powers of function Cfunct(X)
      max2=1
c
      incont4=0
c
cbs   first of all, this routine determines, for which L 
cbs   values the radial integrals have to be solved
cbs   initialize the number of blocks for the different 
cbs   l-combinations  
cbs   no (ss|ss) contributions 
      if (l1.eq.0.and.l2.eq.0.and.l3.eq.0.and.l4.eq.0) return  ! no integrals for <ss|ss> 
      if (makemean) then 
      		nblock=1  ! sp sp are the first, so the first block
      		Lstarter(1)=1
      else 
      CALL QUIT('only mean-field with this version')
      endif 
cbs   keep track of L-values for later purposes   
      Lvalues(1)=l1
      Lvalues(2)=l2
      Lvalues(3)=l3
      Lvalues(4)=l4
cbs   now nanz is given the new value
      nanz=ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
      nprimprod=nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)
c
cbs   prepare the powers needed for cfunctx
c
c
c     There are seven different CASES of integrals following 
c       (   A  --  C) 
c
c     The structure is the same for all cases, therefore comments can be found only on case A
c
c
c
cbs   ###########################################################################################################
cbs   the (+2) cases          CASE A
cbs   ##########################################################################################################
      incl1=1  !  Those increments define the case 
      incl3=1
cbs   determine the possible L-values for the integrals by checking for triangular equation 
c
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
c
cbs   returns first and last L-values (Lanf,Lend), for which 
cbs   radial integrals have to be calculated 
      if(Lend-Lanf.ge.0) then 
cbs   if there are blocks 
      	Lblocks(1)=(Lend-Lanf)/2+1 ! L increases in steps of 2, 
cbs                                       due to parity conservation 
      	Lfirst(1)=Lanf
      	Llast(1)=Lend 
      else 
      	Lblocks(1)=0 
      endif 
cbs   ##########################################################################################################
cbs   the (0) cases         CASE  B
cbs   ##########################################################################################################
      incl1=0
      incl3=0
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(2)=(Lend-Lanf)/2+1
      Lfirst(2)=Lanf
      Llast(2)=Lend 
      Lblocks(3)=(Lend-Lanf)/2+1
      Lfirst(3)=Lanf
      Llast(3)=Lend 
      else 
      Lblocks(2)=0 
      Lblocks(3)=0 
      endif 
      Lstarter(2)=Lstarter(1)+
     *nanz*Lblocks(1)
      Lstarter(3)=Lstarter(2)+
     *nanz*Lblocks(2)
cbs   ##########################################################################################################
cbs   the (-2) cases      CASE C
cbs   ##########################################################################################################
      if (l1.eq.0.or.l3.eq.0) then 
      Lblocks(4)=0
      else 
      incl1=-1
      incl3=-1
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(4)=(Lend-Lanf)/2+1
      Lfirst(4)=Lanf
      Llast(4)=Lend 
      else 
      Lblocks(4)=0 
      endif 
      endif
      Lstarter(4)=Lstarter(3)+
     *nanz*Lblocks(3)
c
CBS   now the hole purpose of this routine 
c
      icont4=Lstarter(4)+nanz*Lblocks(4)
      return 
      end  
      subroutine genovlp(Lhigh)
#include "implicit.h"
#include "para.h"
#include "param.h"
cbs   generates overlap of normalized  primitives. 
      dimension evecinv(MxprimL,MxprimL)
      do L=0,Lhigh 
      	do Jrun=1,nprimit(L)
      	do Irun=1,nprimit(L)
        normovlp(Irun,Jrun,L)=coulovlp(irun,jrun,0,0,
     *  L,L)   
      	enddo 
      	enddo 
cbs   invert the matrix, not very elegant, but sufficient
      ipnt=0
      do jrun=1,nprimit(L)
      do irun=1,jrun
      ipnt=ipnt+1
      scratchinv(ipnt)=normovlp(irun,jrun,L)
      enddo
      enddo
      do Jrun=1,nprimit(L)
      do Irun=1,nprimit(L)
      evecinv(Irun,Jrun)=0d0                                
      enddo 
      enddo 
      do Jrun=1,nprimit(L)
      evecinv(jrun,jrun)=1d0
      enddo 
      call jacobi(scratchinv,evecinv,nprimit(L),MxprimL) 
      do irun=1,nprimit(L)
      eval(irun)=dsqrt(scratchinv((irun*irun+irun)/2))
      enddo
cbs   ensure normalization of the vectors.
      do IRUN=1,nprimit(L)
      fact=0d0
      do JRUN=1,nprimit(L)
      fact=fact+evecinv(JRUN,IRUN)*evecinv(JRUN,IRUN)
      enddo
      fact=1d0/dsqrt(fact)
      do JRUN=1,nprimit(L)
      evecinv(JRUN,IRUN)=fact*evecinv(JRUN,IRUN)
      enddo
      enddo  
cbs   now generate rootOVLP           
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLP(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLP(irun,jrun,L)=rootOVLP(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
cbs   now generate rootOVLPinv           
      do irun=1,nprimit(L)
      eval(irun)=1d0/eval(irun)                            
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLPinv(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLPinv(irun,jrun,L)=rootOVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
cbs   now generate OVLPinv           
      do irun=1,nprimit(L)
      eval(irun)=eval(irun)*eval(irun)
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      OVLPinv(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      OVLPinv(irun,jrun,L)=OVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
      enddo 
      return 
      end   
      subroutine genpowers(Lhigh)   
#include "implicit.h"
#include "para.h"
#include "param.h"
#include "dofuc.h"
cbs   set some often used powers of exponents 
      do L2=0,Lhigh 
      do L1=0,L2    
      do irun1=1,nprimit(L1)
      do irun2=1,nprimit(L2)
      powexp(irun1,irun2,L1,L2,0)=1d0
      enddo
      enddo
      enddo
      enddo
      do L2=0,Lhigh 
      do L1=0,L2    
      do Lrun=1,(L1+L2+5)
      do irun2=1,nprimit(L2)
      do irun1=1,nprimit(L1)
      fact=dsqrt(0.5d0*(exponents(irun1,L1)+exponents(irun2,L2)))
      powexp(irun1,irun2,L1,L2,Lrun)= powexp(irun1,irun2,L1,L2,Lrun-1)*
     *fact
      enddo
      enddo
      enddo
      enddo
      enddo
cbs   generate coulovlp = overlap for normalized functions, but sometimes
cbs   with shifted l-values
      do l2=0,lhigh
      do incl2=-1,1
         if (l2+incl2.ge.0) then  ! do not lower l for s-functions 
         n2=l2+incl2+1
         df2=1d0/dsqrt(df(n2+n2-1))
         do l1=0,l2
         do incl1=-1,1
         if (l1+incl1.ge.0) then ! do not lower l for s-functions 
         n1=l1+incl1+1
         df1=1d0/dsqrt(df(n1+n1-1))
         df12=df(n1+n2-1)
         do iprim2=1,nprimit(l2)
         fact2=dsqrt(powexp(iprim2,iprim2,l2,l2,n2+n2+1))
         factor=fact2*df1*df2*df12
         do iprim1=1,nprimit(l1)
         fact1=dsqrt(powexp(iprim1,iprim1,l1,l1,n1+n1+1))
         coulovlp(iprim1,iprim2,incl1,incl2,l1,l2)= 
     *   fact1*factor/powexp(iprim1,iprim2,l1,l2,n1+n2+1) 
         enddo
         enddo
         endif
         enddo
         enddo
         endif
      enddo
      enddo
      return 
      end 
 
    
      subroutine genstar(Lhigh) 
#include "implicit.h"
#include "para.h"
#include "param.h"
cbs   purpose: generate start adresses of contraction coeffs on 
cbs   contrarray for the different L-Blocks 
      istart=1
      do L=0,Lhigh
      inc=nprimit(L)*ncontrac(L)
      iaddori(L)=istart
      istart=istart+inc  
      iaddtyp1(L)=istart 
      istart=istart+inc  
      iaddtyp2(L)=istart 
      istart=istart+inc  
      iaddtyp3(L)=istart 
      istart=istart+inc  
      iaddtyp4(L)=istart 
      istart=istart+inc  
      enddo 
      return 
      end  
      subroutine gentkin(L,TKIN,nprims,exponents,rootOVLPinv)   
#include "implicit.h"
#include "para.h"
cbs   subroutine to generate the kinetic energy 
      dimension TKIN(nprims,nprims),exponents(*),
     *dummy(MxprimL,MxprimL),dummy2(MxprimL,MxprimL),
     *rootOVLPinv(MxprimL,MxprimL)
cbs   one triangular part of the matrix
      do irun2=1,nprims  
      do irun1=1,irun2   
      	dummy(irun1,irun2)=
     *	Tkinet(l,exponents(irun1),
     *	exponents(irun2))
      enddo     
      enddo     
cbs   copy to the other triangular part.... 
      do irun2=1,nprims-1 
      do irun1=irun2+1,nprims      
      	dummy(irun1,irun2)=dummy(irun2,irun1)                       
      enddo     
      enddo     
cbs   now transform by rootovlp*dummy*rootovlp 
      do jrun=1,nprims
      do irun=1,nprims
      	TKIN(irun,jrun)=0d0
      	dummy2(irun,jrun)=0d0
      enddo
      enddo
      do irun=1,nprims
      do jrun=1,nprims
      do krun=1,nprims
      	dummy2(irun,jrun)=dummy2(irun,jrun)+
     *	dummy(irun,krun)*rootovlpinv(krun,jrun)
      enddo
      enddo
      enddo
      do irun=1,nprims
      do jrun=1,nprims
      do krun=1,nprims
      	Tkin(irun,jrun)=Tkin(irun,jrun)+
     *	dummy2(krun,jrun)*rootovlpinv(irun,krun)
      enddo
      enddo
      enddo
      return 
      end  

CMI ..................................
      subroutine getAOsx(lhigh,IACH,WORK,LWORK)                  
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "para.h"
#include "param.h"
cbs   get expansions of atomic orbitals in contracted functions 
      character*12    occtext,occread
      character*18  textnorbmf,textnorbmf2    
      logical EX
      DIMENSION WORK(LWORK)

#include "amfi_if.h"
c     logical amfi3_verbose
c     amfi3_verbose = .true. 

      occtext='OCCUPATION: ' 
      textnorbmf='Number of orbitals'
      Inquire(File='AO-expansion',exist=EX) 
      if (.not.EX)  then
CBS   write(6,*) 'get occupations from DATA-block'
      call getAOs2(lhigh,IACH,WORK,LWORK)
      return
      endif 

      LUAOEX = -1
      CALL GPOPEN(LUAOEX,'AO-expansion','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

CMI   if (amfi3_verbose) write(LUPRI,*) 'Orbitals for mean-field' 
      if (IPR_AMFI.GE.2) write(LUPRI,*) 'Orbitals for mean-field' 

      do lrun=0,lhigh 

CMI      if (amfi3_verbose) write(LUPRI,'(A3,I3)') 'L= ',lrun 
         if (IPR_AMFI.GE.2) write(LUPRI,'(A3,I3)') 'L= ',lrun 

         read(LUAOEX,'(A18,I3)') textnorbmf2,noccorb(lrun)
         if (textnorbmf.ne.textnorbmf2) 
     *      CALL QUIT('wrong keyword for number of orbitals in getAOsx')

CMI      if (amfi3_verbose) 
         if (IPR_AMFI.GE.2) 
     &      write(LUPRI,*) 'number of orbitals ',noccorb(lrun)

         do iorbital=1,noccorb(lrun)     
            read(LUAOEX,'(A12,F5.3)')  occread,occup(iorbital,lrun)

CMI         if (amfi3_verbose) 
            if (IPR_AMFI.GE.2) 
     &         write(LUPRI,'(A,F8.4)') occtext,occup(iorbital,lrun)
            if (occread.ne.occtext) 
     &         CALL QUIT('error reading AOs')
            read(LUAOEX,*) 
     &         (AOcoeffs(icont,iorbital,lrun),icont=1,ncontrac(lrun))

CMI         if ( amfi3_verbose ) then 
            if ( IPR_AMFI.GE.2 ) then 
               write(LUPRI,'(8F10.4)') (AOcoeffs(icont,iorbital,lrun),
     &            icont=1,ncontrac(lrun))
               write(LUPRI,*) ' '
            endif
            read(LUAOEX,*) 
         enddo
      enddo
      call gpclose(LUAOEX,'KEEP')
      return 
      end 


      subroutine getAOs2(lhigh,IACH,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "param.h"
#include "scfarr.h"
cbs   get expansions of atomic orbitals in contracted functions
      common /nucleus/ charge,Exp_finite
      character*12    occtext
      integer closedshells(0:LMAX),openshells(0:LMAX)
      logical doscf, EX
      dimension xorb(MxprimL,MxContL,0:lmax)
      dimension ovlpprimcont(MxcontL,MxprimL),
     *dummytrafo(MxcontL,MxprimL)
      real(8), allocatable :: zeta(:)
      DIMENSION WORK(LWORK)

#include "amfi_if.h"
c     logical amfi3_verbose
      data doscf / .true. /
#include "memint.h"

c     amfi3_verbose = .false.

      null=0
      ione=1   

      call getocc_ao(int(charge),IACH,closedshells,openshells)

      occtext='OCCUPATION: '
      do lrun=0,lhigh
      do irun=1,MxcontL
      do jrun=1,MxcontL
      AOcoeffs(jrun,irun,lrun)=0d0
      enddo
      enddo
      enddo

!radovan: why is this logical here?
      IF (doscf) THEN
CMI ... run the relscf.x stuff
c
c
c     getting occupied orbitals from the atomic SCF   
c
c
CBS   PREPARE INPUT FOR SCF-CODE
      open(unit=17,file='TOSCF',form='formatted',
     *status='unknown')
      rewind(17)
      write(17,*) charge
      write(17,'(4I3)') (nprimit(i),i=0,Lmax_occ)    
      write(17,'(4I3)') (closedshells(i),i=0,Lmax_occ)
      write(17,'(4I3)') (openshells(i),i=0,Lmax_occ)

CMI    .... write out important info....
      if (IPR_AMFI.GE.0) then
         write(lupri,'(1X,A)')
     &   '****  Written to the file TOSCF for "relscf" ****'
         write(lupri,'(1X,A,F6.3)')
     &   '       charge: ',charge
         write(lupri,'(A,4I3)')
     &   '      nprimit: ',(nprimit(i),i=0,Lmax_occ)  
         write(lupri,'(A,4I3)')
     &   '   closed sh.: ',(closedshells(i),i=0,Lmax_occ)
         write(lupri,'(A,4I3)')
     &   '     open sh.: ',(openshells(i),i=0,Lmax_occ)
      endif

CMI   if(breit) then
CMI   write(17,'(I5)') null   
CMI   else
      write(17,'(I5)') ione    
CMI   endif

!     Pack exponents in a one dimensional array as they will be read that way
!     avoids bug with gnu8 compiler in directly writing in a nested loop
      allocate(zeta(sum(nprimit)))
      ii = 0
      do l=0,Lmax_occ
         do irun=1,nprimit(l)
            ii = ii + 1
            zeta(ii) = exponents(irun,l)
         end do
      end do
      write(17,'(4(E20.14,1X))') zeta
      deallocate(zeta)

      rewind(17)
      close(17,status='keep')
C
C     set file labels for file RELSCF_COEF
      CALL GETLAB_RELSCF(INT(CHARGE))

CMI ... run internal atomic SCF code which reads TOSCF input file
      CALL INT_AT_SCF(WORK,LWORK)
      IF(ONRELSCF) THEN
C     ... that is all - SK - Jan 09
        WRITE(LUPRI,*)' *** INFO from getocc_ao: quit as you wanted ***'
        CALL QUIT(' *** getocc_ao: quit as you wanted *** ')
      END IF
cbs   set the ocupation numbers 
      do lrun=0,lhigh
      do i=1,closedshells(lrun)
      occup(i,lrun)=2.0
      enddo
      noccorb(lrun)=closedshells(lrun)
      if (openshells(lrun).gt.0) then
      i=closedshells(lrun)+1
      occup(i,lrun)=1d0*openshells(lrun)/dfloat(lrun+lrun+1)
      noccorb(lrun)=i
      endif
      enddo
cbs    
cbs   read the MO coefficients   
CMI ======================================================
CMI   if (amfi3_verbose) then
C     IPR_AMFI_SAVE = IPR_AMFI
C     IPR_AMFI = 2
      if (IPR_AMFI.GE.2) then
         write(LUPRI,*) 'AMFI: atomic orbitals to be read from'//
     &                  ' RELSCF_COEF' 
         write(LUPRI,*) 'lhigh = ',lhigh   
      endif
C
C     check if 'RELSCF_COEF' exists.
      LURELSCFCOEF = 7
      INQUIRE ( FILE = 'RELSCF_COEF', EXIST = EX )
      IF(.not. EX) THEN
        WRITE(LUPRI,*) ' *** ERROR in getAOs2: coefficient file not'//
     &                 ' present - missing RELSCF_COEF! ***'
        CALL QUIT('* ERROR in getAOs2: coefficient file not present *')
      END IF
      CALL OPNFIL(LURELSCFCOEF,'RELSCF_COEF','OLD','GETAOs')
      REWIND LURELSCFCOEF
C
C
      CALL IZERO(IRELSCF_INFO,MXRELSCF_INFO)

      CALL MEMGET('REAL',KRELCOEFFX,LDIMC+2,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KRELCOEFFX),LDIMC+2)
C
C     read info array
      CALL IREAKRMC(LURELSCFCOEF,RELSCFLABI,IRELSCF_INFO,
     &              MXRELSCF_INFO)
      IPLACEVLEN = 0
C     WRITE(LUPRI,*) ' IRELSCF_INFO(1)',IRELSCF_INFO(1)
      IPLACEVLEN = (IRELSCF_INFO(1) * 2) + 2
C     WRITE(LUPRI,*) ' IPLACEVLEN = ',IPLACEVLEN
      ISIZEVLEN  = IRELSCF_INFO(IPLACEVLEN)
C
C     read coefficients from file RELSCF_COEF to temporary array
      CALL REAAMFI(LURELSCFCOEF,RELSCFLABX,WORK(KRELCOEFFX),ISIZEVLEN)
C     WRITE(LUPRI,*) ' coeff matrix from file: len = ',ISIZEVLEN
C     CALL WRTMATMN(WORK(KRELCOEFFX),1,ISIZEVLEN,1,ISIZEVLEN,LUPRI)

      CLOSE(UNIT=LURELSCFCOEF,STATUS='KEEP')
C
C     distribute to xorb
      lrelscfvec = 1
      do lrun=0,lhigh
CMI      if (amfi3_verbose)
         if (IPR_AMFI.GE.2)
     &      write(LUPRI,*) 'lrun = ',lrun,
     &         ' nprimit=',nprimit(lrun),' nocc=',noccorb(lrun)

         if (noccorb(lrun).gt.0) then
C           consistency check
            if (IRELSCF_INFO((1+(2*lrun+1))).ne.nprimit(lrun))
     &         stop 'fuck 1'
            if (IRELSCF_INFO((1+(2*lrun+2))).ne.noccorb(lrun))
     &         stop 'fuck 2'

            do irun=1,nprimit(lrun)
               do jrun=1,noccorb(lrun)
                 xorb(irun,jrun,lrun) = WORK(KRELCOEFFX+lrelscfvec-1)
                 lrelscfvec = lrelscfvec + 1
               enddo
CMI            if (amfi3_verbose)
               if (IPR_AMFI.GE.2)
     &            WRITE(LUPRI,*)
     &                   (xorb(irun,jrun,lrun),jrun=1,noccorb(lrun))
            enddo
         endif
      enddo
C     ... release memmory
      CALL MEMREL('getAOs2',WORK,KWORK,KWORK,KFREE,LFREE)
      CLOSE(UNIT=LURELSCFCOEF,STATUS='KEEP')

CMI   if (amfi3_verbose) write(LUPRI,*) '-- orbitals read from fort.7' 
      if (IPR_AMFI.GE.2) write(LUPRI,*)
     &  '-- orbitals read from RELSCF_COEF, closed'

C     IPR_AMFI = IPR_AMFI_SAVE

      do lrun=0,lhigh
      if (noccorb(lrun).gt.0) then
cbs     prepare projection and normalization of SCF-vectors
cbs     on contracted basis set
        do irun=1,ncontrac(lrun)
        do jrun=1,ncontrac(lrun)
        ovlpcontinv(jrun,irun,lrun)=ovlpcont(jrun,irun,lrun)
        enddo
        enddo
CMI ... hmm.. a possible source of dammage...
        call invers(ovlpcontinv(1,1,lrun),MxcontL,ncontrac(lrun))

        do jrun=1,nprimit(lrun)
        do irun=1,ncontrac(lrun)
        ovlpprimcont(irun,jrun)=0d0
        dummytrafo(irun,jrun)=0d0
        enddo
        enddo
        do irun=1,ncontrac(lrun)
        do jrun=1,nprimit(lrun)
        do krun=1,nprimit(lrun)
        dummytrafo(irun,jrun)=dummytrafo(irun,jrun)+
     *  contrarray(iaddori(lrun)-1+(irun-1)*nprimit(lrun)+
     *  krun)*normovlp(krun,jrun,lrun)
        enddo
        enddo
        enddo
        do irun=1,ncontrac(lrun)
        do jrun=1,ncontrac(lrun)
        do krun=1,nprimit(lrun)
        ovlpprimcont(irun,krun)=ovlpprimcont(irun,krun)+
     *  dummytrafo(jrun,krun)*ovlpcontinv(irun,jrun,lrun)
        enddo
        enddo
        enddo
        do iorb=1,noccorb(lrun)
        do irun=1,ncontrac(lrun)
        xsum=0d0
        do krun=1,nprimit(lrun)
        xsum=xsum+xorb(krun,iorb,lrun)*ovlpprimcont(irun,krun)
        enddo
        AOcoeffs(irun,iorb,lrun)=xsum
        enddo
cbs   normalize orbitals
      dnorm=getnorm(ncontrac(lrun),AOcoeffs(1,iorb,lrun),
     *ovlpcont(1,1,lrun),MxcontL)

CMI   if (amfi3_verbose)
      if (IPR_AMFI.GE.2)
     &   write(LUPRI,*) 'squared overlap of orbital ',iorb,' L-value ',
     *lrun,dnorm

      do icont=1,ncontrac(lrun)
      AOcoeffs(icont,iorb,lrun)=
     *AOcoeffs(icont,iorb,lrun)/dsqrt(dnorm)
      enddo
      enddo
      endif    
      Enddo
C

      ELSE
CMI ... this part of the code is frozen 
CMI    as we do not utilize the contractions coefficients as
CMI    mean-field orbitals
CBS   write(6,*) 'Orbitals for mean-field'
      do lrun=0,lhigh
CBS   write(6,'(A3,I3)') 'L= ',lrun
      do i=1,closedshells(lrun)
      occup(i,lrun)=2.0
      AOcoeffs(i,i,lrun)=1d0
      enddo
      noccorb(lrun)=closedshells(lrun)
      if (openshells(lrun).gt.0) then
      i=closedshells(lrun)+1
      occup(i,lrun)=1d0*openshells(lrun)/dfloat(lrun+lrun+1)
      AOcoeffs(i,i,lrun)=1d0
      noccorb(lrun)=i
      endif
      enddo

      ENDIF
      return
      end
cbs
      subroutine getocc_ao(icharge,IACH,iclosed,iopen)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
      parameter (ichargemax=118)
      dimension iclocc(0:Lmax_occ,0:ichargemax)
      dimension iopocc(0:Lmax_occ,0:ichargemax)
      character*30 occtxt(0:ichargemax) 
      character*32 txt  
#include "amfi_if.h"
c     logical amfi3_verbose
CMI   data txt/'SO-integrals are calculated for '/
      data txt/'    ...electronic occupation of '/

      dimension iclosed(0:LMAX),iopen(0:LMAX)
      data (occtxt(i),i=0,ichargemax) /
     *'dummy atom (no integrals)     ',  ! 0
CMI  *' H: no mean-field             ',
     *' H: 1s^1                      ','He: 1s^2                      ', ! 1-2
     *'Li: [He]2s^1                  ','Be: [He]2s^2                  ', ! 3-4
     *' B: [He]2s^2 2p^1             ',' C: [He]2s^2 2p^2             ', ! 5-6
     *' N: [He]2s^2 2p^3             ',' O: [He]2s^2 2p^4             ', ! 7-8
     *' F: [He]2s^2 2p^5             ','Ne: [He]2s^2 2p^6             ', ! 9-10
     *'Na: [Ne]3s^1                  ','Mg: [Ne]3s^2                  ', ! 11-12
     *'Al: [Ne]3s^2 3p^1             ','Si: [Ne]3s^2 3p^2             ',
     *' P: [Ne]3s^2 3p^3             ',' S: [Ne]3s^2 3p^4             ',
     *'Cl: [Ne]3s^2 3p^5             ','Ar: [Ne]3s^2 3p^6             ',
     *' K: [Ar]4s^1                  ','Ca: [Ar]4s^2                  ',
     *'Sc: [Ar]4s^2 3d^1             ','Ti: [Ar]4s^2 3d^2             ',
     *' V: [Ar]4s^2 3d^3             ','Cr: [Ar]4s^2 3d^4             ',
     *'Mn: [Ar]4s^2 3d^5             ','Fe: [Ar]4s^2 3d^6             ',
     *'Co: [Ar]4s^2 3d^7             ','Ni: [Ar]4s^2 3d^8             ',
     *'Cu: [Ar]4s^1 3d^10            ','Zn: [Ar]4s^2 3d^10            ',
     *'Ga: [Ar]4s^2 3d^10 4p^1       ','Ge: [Ar]4s^2 3d^10 4p^2       ',
     *'As: [Ar]4s^2 3d^10 4p^3       ','Se: [Ar]4s^2 3d^10 4p^4       ',
     *'Br: [Ar]4s^2 3d^10 4p^5       ','Kr: [Ar]4s^2 3d^10 4p^6       ',
     *'Rb: [Kr]5s^1                  ','Sr: [Kr]5s^2                  ',
     *' Y: [Kr]5s^2 4d^1             ','Zr: [Kr]5s^2 4d^2             ',
     *'Nb: [Kr]5s^2 4d^3             ',
     *'Mo: [Kr]5s^2 4d^4             ',
     *'Tc: [Kr]5s^2 4d^5             ',
     *'Ru: [Kr]5s^2 4d^6             ',
     *'Rh: [Kr]5s^2 4d^7             ',
     *'Pd: [Kr]5s^2 4d^8             ',
     *'Ag: [Kr]5s^1 4d^10            ',
     *'Cd: [Kr]5s^2 4d^10            ',
     *'In: [Kr]5s^2 4d^10 5p^1       ',
     *'Sn: [Kr]5s^2 4d^10 5p^2       ',
     *'Sb: [Kr]5s^2 4d^10 5p^3       ',
     *'Te: [Kr]5s^2 4d^10 5p^4       ',
     *' I: [Kr]5s^2 4d^10 5p^5       ',
     *'Xe: [Kr]5s^2 4d^10 5p^6       ',
     *'Cs: [Xe]6s^1                  ',
     *'Ba: [Xe]6s^2                  ',
     *'La: [Xe]6s^2 5d^1             ',
     *'Ce: [Xe]6s^2 4f^2             ',
     *'Pr: [Xe]6s^2 4f^3             ',
     *'Nd: [Xe]6s^2 4f^4             ',
     *'Pm: [Xe]6s^2 4f^5             ',
     *'Sm: [Xe]6s^2 4f^6             ',
     *'Eu: [Xe]6s^2 4f^7             ',
     *'Gd: [Xe]6s^2 4f^8             ',
     *'Tb: [Xe]6s^2 4f^9             ',
     *'Dy: [Xe]6s^2 4f^10            ',
     *'Ho: [Xe]6s^2 4f^11            ',
     *'Er: [Xe]6s^2 4f^12            ',
     *'Tm: [Xe]6s^2 4f^13            ',
     *'Yb: [Xe]6s^2 4f^14            ',
     *'Lu: [Xe+4f^14]6s^2 5d^1       ',
     *'Hf: [Xe+4f^14]6s^2 5d^2       ',
     *'Ta: [Xe+4f^14]6s^2 5d^3       ',
     *' W: [Xe+4f^14]6s^2 5d^4       ',
     *'Re: [Xe+4f^14]6s^2 5d^5       ',
     *'Os: [Xe+4f^14]6s^2 5d^6       ',
     *'Ir: [Xe+4f^14]6s^2 5d^7       ',
     *'Pt: [Xe+4f^14]6s^1 5d^9       ', 
     *'Au: [Xe+4f^14]6s^1 5d^10      ',
     *'Hg: [Xe+4f^14]6s^2 5d^10      ',
     *'Tl: [Xe+4f^14+5d^10]6s^2 6p^1 ',
     *'Pb: [Xe+4f^14+5d^10]6s^2 6p^2 ',
     *'Bi: [Xe+4f^14+5d^10]6s^2 6p^3 ',
     *'Po: [Xe+4f^14+5d^10]6s^2 6p^4 ',
     *'At: [Xe+4f^14+5d^10]6s^2 6p^5 ',
     *'Rn: [Xe+4f^14+5d^10]6s^2 6p^6 ',
     *'Fr: [Rn]7s^1                  ',
     *'Ra: [Rn]7s^2                  ',
     *'Ac: [Rn]7s^2 6d^1             ',
     *'Th: [Rn]7s^2 6d^2             ',
     *'Pa: [Rn]7s^2 6d^1 5f^2        ',
     *' U: [Rn]7s^2 6d^1 5f^3        ',
     *'Np: [Rn]7s^2 6d^1 5f^4        ',
     *'Pu: [Rn]7s^2 6d^0 5f^6        ',
     *'Am: [Rn]7s^2 6d^0 5f^7        ',
     *'Cm: [Rn]7s^2 6d^0 5f^8        ',  ! 96
     *'Bk: [Rn]7s^2 6d^0 5f^9        ',  ! 97
     *'Cf: [Rn]7s^2 6d^0 5f^10       ',  ! 98
     *'Es: [Rn]7s^2 6d^0 5f^11       ',  ! 99
     *'Fm: [Rn]7s^2 6d^0 5f^12       ',  ! 100
     *'Md: [Rn]7s^2 6d^0 5f^13       ',  ! 101
     *'No: [Rn]7s^2 6d^0 5f^14       ',  ! 102
     *'Lr: [Rn]7s^2 7p^1 6d^0 5f^14  ',  ! 103
     *'Rf: [Rn]7s^2 7p^0 6d^2 5f^14  ',  ! 104
     *'Db: [Rn]7s^2 7p^0 6d^3 5f^14  ',  ! 105
     *'Sg: [Rn]7s^2 7p^0 6d^4 5f^14  ',  ! 106
     *'Bh: [Rn]7s^2 7p^0 6d^5 5f^14  ',  ! 107
     *'Hs: [Rn]7s^2 7p^0 6d^6 5f^14  ',  ! 108
     *'Mt: [Rn]7s^2 7p^0 6d^7 5f^14  ',  ! 109
     *'Ds: [Rn]7s^1 7p^0 6d^9 5f^14  ',  ! 110
     *'Rg: [Rn]7s^1 7p^0 6d^10 5f^14 ',  ! 111
     *'Cn: [Rn]7s^2 7p^0 6d^10 5f^14 ',  ! 112
     *'Nh: [Rn]7s^2 7p^1 6d^10 5f^14 ',  ! 113
     *'Fl: [Rn]7s^2 7p^2 6d^10 5f^14 ',  ! 114
     *'Mc: [Rn]7s^2 7p^3 6d^10 5f^14 ',  ! 115
     *'Lv: [Rn]7s^2 7p^4 6d^10 5f^14 ',  ! 116
     *'Ts: [Rn]7s^2 7p^5 6d^10 5f^14 ',  ! 117
     *'Og: [Rn]7s^2 7p^6 6d^10 5f^14 '/  ! 118
      data ((iclocc(i,j),i=0,LMAX_occ),j=0,ichargemax) /
     & 0 , 0, 0, 0,  0 , 0, 0, 0,  1 , 0, 0, 0,  ! 0-2
     & 1 , 0, 0, 0,  2 , 0, 0, 0,  2 , 0, 0, 0,  ! 3-5      
     & 2 , 0, 0, 0,  2 , 0, 0, 0,  2 , 0, 0, 0,  ! 6-8       
     & 2 , 0, 0, 0,  2 , 1, 0, 0,  2 , 1, 0, 0,  ! 9-11    
c
     & 3 , 1, 0, 0,  3 , 1, 0, 0,  3 , 1, 0, 0,  ! 12-14     
     & 3 , 1, 0, 0,  3 , 1, 0, 0,  3 , 1, 0, 0,  ! 15-17
     & 3 , 2, 0, 0,  3 , 2, 0, 0,  4 , 2, 0, 0,  ! 18-20  
c 
     & 4 , 2, 0, 0,  4 , 2, 0, 0,  4 , 2, 0, 0,  ! 21-23    
     & 4 , 2, 0, 0,  4 , 2, 0, 0,  4 , 2, 0, 0,  ! 24-26    
     & 4 , 2, 0, 0,  4 , 2, 0, 0,  3 , 2, 1, 0,  ! 27-29    
     & 4 , 2, 1, 0,       !30   
c 
     & 4 , 2, 1, 0,       !31 
     & 4 , 2, 1, 0,       !32 
     & 4 , 2, 1, 0,       !33  
     & 4 , 2, 1, 0,       !34 
     & 4 , 2, 1, 0,       !35 
     & 4 , 3, 1, 0,       !36 
     & 4 , 3, 1, 0,       !37 
     & 5 , 3, 1, 0,       !38 
     & 5 , 3, 1, 0,       !39 
     & 5 , 3, 1, 0,       !40 
c 
     & 5 , 3, 1, 0,       !41 
     & 5 , 3, 1, 0,       !42 
     & 5 , 3, 1, 0,       !43 
     & 5 , 3, 1, 0,       !44 
     & 5 , 3, 1, 0,       !45 
     & 5 , 3, 1, 0,       !46 
     & 4 , 3, 2, 0,       !47 
     & 5 , 3, 2, 0,       !48 
c 
     & 5 , 3, 2, 0,       !49  
     & 5 , 3, 2, 0,       !50  
     & 5 , 3, 2, 0,       !51  
     & 5 , 3, 2, 0,       !52  
     & 5 , 3, 2, 0,       !53  
     & 5 , 4, 2, 0,       !54  
     & 5 , 4, 2, 0,       !55  
     & 6 , 4, 2, 0,       !56 
     & 6 , 4, 2, 0,       !57 
     & 6 , 4, 2, 0,       !58 
     & 6 , 4, 2, 0,       !59 
     & 6 , 4, 2, 0,       !60 
c  
     & 6 , 4, 2, 0,       !61  
     & 6 , 4, 2, 0,       !62  
     & 6 , 4, 2, 0,       !63  
     & 6 , 4, 2, 0,       !64  
     & 6 , 4, 2, 0,       !65  
     & 6 , 4, 2, 0,       !66  
     & 6 , 4, 2, 0,       !67  
     & 6 , 4, 2, 0,       !68  
     & 6 , 4, 2, 0,       !69  
     & 6 , 4, 2, 1,       !70  
c   
     & 6 , 4, 2, 1,       !71 
     & 6 , 4, 2, 1,       !72 
     & 6 , 4, 2, 1,       !73 
     & 6 , 4, 2, 1,       !74 
     & 6 , 4, 2, 1,       !75 
     & 6 , 4, 2, 1,       !76 
     & 6 , 4, 2, 1,       !77 
     & 5 , 4, 2, 1,       !78 
     & 5 , 4, 3, 1,       !79 
     & 6 , 4, 3, 1,       !80 
c
     & 6 , 4, 3, 1,       !81  
     & 6 , 4, 3, 1,       !82  
     & 6 , 4, 3, 1,       !83  
     & 6 , 4, 3, 1,       !84  
     & 6 , 4, 3, 1,       !85  
     & 6 , 5, 3, 1,       !86  
     & 6 , 5, 3, 1,       !87  
     & 7 , 5, 3, 1,       !88 
     & 7 , 5, 3, 1,       !89 
     & 7 , 5, 3, 1,       !90 
c    
     & 7 , 5, 3, 1,       !91  
     & 7 , 5, 3, 1,       !92  
     & 7 , 5, 3, 1,       !93  
     & 7 , 5, 3, 1,       !94  
     & 7 , 5, 3, 1,       !95 
     & 7 , 5, 3, 1,       !96 
     & 7 , 5, 3, 1,       !97 
     & 7 , 5, 3, 1,       !98 
     & 7 , 5, 3, 1,       !99 
     & 7 , 5, 3, 1,       !100
     & 7 , 5, 3, 1,       !101
     & 7 , 5, 3, 2,       !102
     & 7 , 5, 3, 2,       !103
     & 7 , 5, 3, 2,       !104
     & 7 , 5, 3, 2,       !105
     & 7 , 5, 3, 2,       !106
     & 7 , 5, 3, 2,       !107
     & 7 , 5, 3, 2,       !108
     & 7 , 5, 3, 2,       !109
     & 6 , 5, 3, 2,       !110
     & 6 , 5, 4, 2,       !111, cs+os
     & 7 , 5, 4, 2,       !112, cs
     & 7 , 5, 4, 2,       !113, cs+os
     & 7 , 5, 4, 2,       !114, cs+os
     & 7 , 5, 4, 2,       !115, cs+os
     & 7 , 5, 4, 2,       !116, cs+os
     & 7 , 5, 4, 2,       !117, cs+os
     & 7 , 6, 4, 2/       !118, cs
cbs
      data ((iopocc(i,j),i=0,LMAX_occ),j=0,ichargemax) /
     & 0 , 0, 0, 0,   1 , 0, 0, 0,   0 , 0, 0, 0, ! 0-2
     & 1 , 0, 0, 0,   0 , 0, 0, 0,   0 , 1, 0, 0, ! 3-5
     & 0 , 2, 0, 0,   0 , 3, 0, 0,   0 , 4, 0, 0, ! 6-8
     & 0 , 5, 0, 0,   0 , 0, 0, 0,   1 , 0, 0, 0, ! 9-11
     & 0 , 0, 0, 0,   0 , 1, 0, 0,   0 , 2, 0, 0, !12-14 
     & 0 , 3, 0, 0,   0 , 4, 0, 0,   0 , 5, 0, 0, !15-17  
     & 0 , 0, 0, 0,   1 , 0, 0, 0,   0 , 0, 0, 0, !18-20   
     & 0 , 0, 1, 0,   0 , 0, 2, 0,   0 , 0, 3, 0, !21-23
     & 0 , 0, 4, 0,   0 , 0, 5, 0,   0 , 0, 6, 0, ! 24-26   
     & 0 , 0, 7, 0,   0 , 0, 8, 0,   1 , 0, 0, 0, ! 27-29
     & 0 , 0, 0, 0,    ! 30  
c
     & 0 , 1, 0, 0,    ! 31  
     & 0 , 2, 0, 0,    ! 32  
     & 0 , 3, 0, 0,    ! 33  
     & 0 , 4, 0, 0,    ! 34   
     & 0 , 5, 0, 0,    ! 35  
     & 0 , 0, 0, 0,    ! 36   
     & 1 , 0, 0, 0,    ! 37   
     & 0 , 0, 0, 0,    ! 38  
     & 0 , 0, 1, 0,    ! 39  
     & 0 , 0, 2, 0,    ! 40  
c
     & 0 , 0, 3, 0,    ! 41 
     & 0 , 0, 4, 0,    ! 42  
     & 0 , 0, 5, 0,    ! 43  
     & 0 , 0, 6, 0,    ! 44  
     & 0 , 0, 7, 0,    ! 45   
     & 0 , 0, 8, 0,    ! 46   
     & 1 , 0, 0, 0,    ! 47  
     & 0 , 0, 0, 0,    ! 48  
     & 0 , 1, 0, 0,    ! 49  
     & 0 , 2, 0, 0,    ! 50  
c
     & 0 , 3, 0, 0,    ! 51  
     & 0 , 4, 0, 0,    ! 52  
     & 0 , 5, 0, 0,    ! 53  
     & 0 , 0, 0, 0,    ! 54  
     & 1 , 0, 0, 0,    ! 55  
     & 0 , 0, 0, 0,    ! 56     
     & 0 , 0, 1, 0,    ! 57  
     & 0 , 0, 0, 2,    ! 58  
     & 0 , 0, 0, 3,    ! 59  
     & 0 , 0, 0, 4,    ! 60  
c
     & 0 , 0, 0, 5,    ! 61  
     & 0 , 0, 0, 6,    ! 62  
     & 0 , 0, 0, 7,    ! 63  
     & 0 , 0, 0, 8,    ! 64  
     & 0 , 0, 0, 9,    ! 65  
     & 0 , 0, 0, 10,    ! 66   
     & 0 , 0, 0, 11,    ! 67  
     & 0 , 0, 0, 12,    ! 68  
     & 0 , 0, 0, 13,    ! 69  
     & 0 , 0, 0,  0,    ! 70  
c
     & 0 , 0, 1, 0,    ! 71  
     & 0 , 0, 2, 0,    ! 72  
     & 0 , 0, 3, 0,    ! 73  
     & 0 , 0, 4, 0,    ! 74  
     & 0 , 0, 5, 0,    ! 75  
     & 0 , 0, 6, 0,    ! 76  
     & 0 , 0, 7, 0,    ! 77   
     & 1 , 0, 9, 0,    ! 78  
     & 1 , 0, 0, 0,    ! 79   
     & 0 , 0, 0, 0,    ! 80   
c
     & 0 , 1, 0, 0,    ! 81  
     & 0 , 2, 0, 0,    ! 82  
     & 0 , 3, 0, 0,    ! 83  
     & 0 , 4, 0, 0,    ! 84  
     & 0 , 5, 0, 0,    ! 85  
     & 0 , 0, 0, 0,    ! 86  
     & 1 , 0, 0, 0,    ! 87  
     & 0 , 0, 0, 0,    ! 88   
     & 0 , 0, 1, 0,    ! 89   
     & 0 , 0, 2, 0,    ! 90  
c
     & 0 , 0, 1, 2,    ! 91  
     & 0 , 0, 1, 3,    ! 92  
     & 0 , 0, 1, 4,    ! 93  
     & 0 , 0, 0, 6,    ! 94  
     & 0 , 0, 0, 7,    ! 95  
     & 0 , 0, 0, 8,    ! 96   
     & 0 , 0, 0, 9,    ! 97   
     & 0 , 0, 0,10,    ! 98   
     & 0 , 0, 0,11,    ! 99   
     & 0 , 0, 0,12,    ! 100  
     & 0 , 0, 0,13,    ! 101  
     & 0 , 0, 0, 0,    ! 102  
     & 0 , 1, 0, 0,    ! 103  
     & 0 , 0, 2, 0,    ! 104  
     & 0 , 0, 3, 0,    ! 105  
     & 0 , 0, 4, 0,    ! 106  
     & 0 , 0, 5, 0,    ! 107  
     & 0 , 0, 6, 0,    ! 108  
     & 0 , 0, 7, 0,    ! 109  
     & 1 , 0, 9, 0,    ! 110  
     & 1 , 0, 0, 0,    ! 111  
     & 0 , 0, 0, 0,    ! 112  
     & 0 , 1, 0, 0,    ! 113  
     & 0 , 2, 0, 0,    ! 114  
     & 0 , 3, 0, 0,    ! 115  
     & 0 , 4, 0, 0,    ! 116  
     & 0 , 5, 0, 0,    ! 117  
     & 0 , 0, 0, 0/    ! 118  
cbs
c     amfi3_verbose = .false.
      if (icharge.gt.ichargemax) then
         CALL QUIT('occupations not implemented')
      endif

CMI   if (amfi3_verbose) then
CMI   if (IPR_AMFI.GE.1) then
      if (IPR_AMFI.GE.0) then
CMI      write(LUPRI,'(A32,A30)') txt,occtxt(icharge)
         write(LUPRI,'(2x,a,i3)') 'charge on the calculated atom:',IACH
         write(LUPRI,'(2x,a,i4)') 
     &    'Mean-field summation for electrons #:',icharge-IACH
         write(LUPRI,'(A32,A30)') txt,occtxt(icharge-IACH)
      endif

      do irun=0,min(lmax,lmax_occ)
CMI      iclosed(irun)=iclocc(irun,icharge)
         iclosed(irun)=iclocc(irun,icharge-IACH)
CMI      iopen(irun)=iopocc(irun,icharge)
         iopen(irun)=iopocc(irun,icharge-IACH)
      end do
      do irun=min(lmax,lmax_occ)+1,lmax
         iclosed(irun)=0
         iopen(irun)=0
      end do
      return
      end
      subroutine getlab_relscf(icharge)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
      parameter (ichargemax=118)
      character*3 occtxt(0:ichargemax) 
      character*4 textin, textrl
#include "amfi_if.h"

      data (occtxt(i),i=0,ichargemax) /
     *'dum',  ! 0
CMI  *' H: no mean-field             ',
     *' H ','He ', ! 1-2
     *'Li ','Be ', ! 3-4
     *' B ',' C ', ! 5-6
     *' N ',' O ', ! 7-8
     *' F ','Ne ', ! 9-10
     *'Na ','Mg ', ! 11-12
     *'Al ','Si ',
     *' P ',' S ',
     *'Cl ','Ar ',
     *' K ','Ca ',
     *'Sc ','Ti ',
     *' V ','Cr ',
     *'Mn ','Fe ',
     *'Co ','Ni ',
     *'Cu ','Zn ',
     *'Ga ','Ge ',
     *'As ','Se ',
     *'Br ','Kr ',
     *'Rb ','Sr ',
     *' Y ','Zr ',
     *'Nb ',
     *'Mo ',
     *'Tc ',
     *'Ru ',
     *'Rh ',
     *'Pd ',
     *'Ag ',
     *'Cd ',
     *'In ',
     *'Sn ',
     *'Sb ',
     *'Te ',
     *' I ',
     *'Xe ',
     *'Cs ',
     *'Ba ',
     *'La ',
     *'Ce ',
     *'Pr ',
     *'Nd ',
     *'Pm ',
     *'Sm ',
     *'Eu ',
     *'Gd ',
     *'Tb ',
     *'Dy ',
     *'Ho ',
     *'Er ',
     *'Tm ',
     *'Yb ',
     *'Lu ',
     *'Hf ',
     *'Ta ',
     *' W ',
     *'Re ',
     *'Os ',
     *'Ir ',
     *'Pt ', 
     *'Au ',
     *'Hg ',
     *'Tl ',
     *'Pb ',
     *'Bi ',
     *'Po ',
     *'At ',
     *'Rn ',
     *'Fr ',
     *'Ra ',
     *'Ac ',
     *'Th ',
     *'Pa ',
     *' U ',
     *'Np ',
     *'Pu ',
     *'Am ',
     *'Cm ',  ! 96
     *'Bk ',  ! 97
     *'Cf ',  ! 98
     *'Es ',  ! 99
     *'Fm ',  ! 100
     *'Md ',  ! 101
     *'No ',  ! 102
     *'Lr ',  ! 103
     *'Rf ',  ! 104
     *'Db ',  ! 105
     *'Sg ',  ! 106
     *'Bh ',  ! 107
     *'Hs ',  ! 108
     *'Mt ',  ! 109
     *'Ds ',  ! 110
     *'Rg ',  ! 111
     *'Cn' ,  ! 112
     *'Nh' ,  ! 113
     *'Fl' ,  ! 114
     *'Mc' ,  ! 115
     *'Lv' ,  ! 116
     *'Ts' ,  ! 117
     *'Og' /  ! 118
cbs
      if (icharge.gt.ichargemax) then
         CALL QUIT('occupations not implemented')
      endif
      textin = 'inte'
      textrl = 'real'
C
C     set file labels for file RELSCF_COEF
      WRITE(RELSCFLABI,'(A3,A1,A4)') occtxt(icharge),'.',textin
      WRITE(RELSCFLABX,'(A3,A1,A4)') occtxt(icharge),'.',textrl

      end

      double precision function  getCG(
     *j1,     ! integer  2*j1
     *j2,     ! integer  2*j2
     *j3,     ! integer  2*j3
     *m1,     ! integer  2*m1
     *m2,     ! integer  2*m2
     *m3)     ! integer  2*m2
cbs this routine calculates the Clebsch-Gordon-coefficients
cbs by actually calculating the 3j-symbol 
cbs  ---                 ---
cbs  |  j1   j2    |   j3   |         j1+m1+j2-m2     
cbs  |             |        |  =  (-)                 sqrt (2  j3+1) *
cbs  |  m1   m2    |   m3   |
cbs  ---                 ---
cbs
cbs                             ---             ---
cbs                             |  j1   j2   j3   |     
cbs                             |                 |     
cbs                             |  m1   m2  -m3   |
cbs                              ---            ---
#include "implicit.h"
cbs   initialize CG-coefficient
      getCG=0d0
cbs   quick check 
      if (m1+m2.ne.m3) return 
      if (j1.lt.0.or.j2.lt.0.or.j3.lt.0) return   
cbs   check the correct sign    beginning  
      idummy=(j1+j2+m1-m2)/2
      if (mod(idummy,2).eq.0) then 
      isign=1
      else
      isign=-1
      endif 
cbs   check the correct sign    end          
      fac1=dsqrt(dfloat(j3+1))
      fac2=regge3j(j1,j2,j3,m1,m2,-m3)
      getCG=isign*fac1*fac2
      return 
      end 


      subroutine getLIMIT(l1,l2,l3,l4,Lanf,Lend)
#include "implicit.h"
#include "priunit.h"
cbs   get the minimum and maximum L-values 
cbs   of the the coulomb-potential to interact 
cbs   with l1-l4
      lower1=iabs(l1-l3) 
      lower2=iabs(l2-l4)         
      lupper1=l1+l3 
      lupper2=l2+l4
      Lanf=max(lower1,lower2)
      Lend=min(lupper1,lupper2) 
cbs     check for parity   
      lsum=Lanf+l1+l3
      if (mod(lsum,2).eq.1) Lanf=Lanf+1
      lsum=Lend+l1+l3
      if (mod(lsum,2).eq.1) Lend=Lend-1   
cbs   check the other parity 
      lsum=Lanf+l2+l4
      if (mod(lsum,2).eq.1) then 
      write(LUPRI,*) ' error in getLIMIT: '
      write(LUPRI,*) ' parity inconsistency for '
      write(LUPRI,*) 'l1,l2,l3,l4= ',l1,l2,l3,l4
      CALL QUIT('Inconsistency error in getLIMIT')              
      endif  
      return 
      end

 
      subroutine getpow(max,quot,quotpow,
     *nprim1,nprim2,nprim3,nprim4)
cbs   generates some powers of for the prefactors of cfunct(X)
cbs   look out for details there and in initfrac   
#include "implicit.h"
#include "para.h"
      dimension quotpow(nprim1,nprim2,
     *nprim3,nprim4),
     *quot(nprim1,nprim2,nprim3,nprim4) 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4)=
     *dsqrt(quot(irun1,irun2,irun3,irun4))
      enddo
      enddo
      enddo
      enddo
      if (max.eq.1) return 
cbs  
      do irun=2,max     
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4)=
     *quotpow(irun1,irun2,irun3,irun4)*
     *quot(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      enddo
      return 
      end

 
      subroutine inidf
cbs   initializes the df on common block  with double facultatives
#include "implicit.h"
#include "para.h"
#include "param.h"
#include "dofuc.h"
      df(0)=1.d0
      df(1)=1.d0
      do irun=2,ndfmx
      df(irun)=dfloat(irun)*df(irun-2)
      enddo
      do jbm=0,ndfmx-1  
      do ibm=jbm,ndfmx   
      dffrac(ibm,jbm)=df(ibm)/df(jbm) 
      enddo 
      enddo 
      do jbm=1,ndfmx  
      do ibm=0,jbm-1       
      dffrac(ibm,jbm)=1d0/dffrac(jbm,ibm)   
      enddo 
      enddo 
      return 
      end
      subroutine initfrac(nprimit1,nprimit2,
     *nprimit3,nprimit4,
     *quot1,quot2,expo1,expo2,
     *expo3,expo4)
cbs   initialize some arrays with factors  needed for cfunct(x) 
#include "implicit.h"
      dimension expo1(*),expo2(*),expo3(*),expo4(*), 
     *quot1(nprimit1,nprimit2,nprimit3,nprimit4),
     *quot2(nprimit1,nprimit2,nprimit3,nprimit4)
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
         sum24=expo2(irun2)+expo4(irun4)
              do irun1=1,nprimit1 
      		quot1(irun1,irun2,irun3,irun4)= 
     *		1d0/(1d0+(expo1(irun1)+expo3(irun3))/
     *		sum24)                             
              enddo
      enddo
      enddo
      enddo
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
      	sum24=expo2(irun2)+expo4(irun4)
      		do irun1=1,nprimit1 
      		quot2(irun1,irun2,irun3,irun4)= 
     *		1d0/(1d0+sum24/
     *		(expo1(irun1)+expo3(irun3)))
      		enddo
      enddo
      enddo
      enddo
      return
      end

      double precision function getnorm(norb,
     *AOcoeffs,ovlpcont,idim)
      implicit double precision(a-h,o-z)
      dimension AOcoeffs(*),ovlpcont(idim,*)
      getnorm=0d0
      do irun2=1,norb
      do irun1=1,norb
      getnorm=getnorm+AOcoeffs(irun2)*AOcoeffs(irun1)*
     *ovlpcont(irun1,irun2)
      enddo
      enddo
      return
      end



C********************************************************************
      SUBROUTINE INVERS(A,NMAX,N)
C********************************************************************
C     THIS SUBROUTINE IS TAKEN FROM ANPASS PROGRAM BY PETER TAYLOR
C********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     THIS SUBROUTINE INVERTS THE MATRIX A OF
C     REAL DIMENSION N x N AND NOMINAL DIMENSION NMAX x NMAX.
C     THE INVERTED MATRIX IS RETURNED AS A.
C
      DIMENSION A(NMAX,NMAX)
      D=1.D0
      DO 80 K=1,N
        BIGA=A(K,K)
   48   DO 55 I=1,N
          IF(I.NE.K)A(I,K)=-A(I,K)/BIGA
   55   CONTINUE
        DO 65 I=1,N
          HOLD=A(I,K)
          DO 62 J=1,N
            IF((I.NE.K).AND.(J.NE.K)) A(I,J)=HOLD*A(K,J)+A(I,J)
   62     CONTINUE
   65   CONTINUE
        DO 75 J=1,N
          IF (J.NE.K) A(K,J)=A(K,J)/BIGA
   75   CONTINUE
        D=D*BIGA
      A(K,K)=1.D0/BIGA
   80 CONTINUE
      RETURN
      END
