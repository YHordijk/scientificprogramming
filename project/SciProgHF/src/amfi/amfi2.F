      subroutine initired
#include "implicit.h"
cbs   initialize all information for ireducible representations 
cbs   later on, it might be useful to have a switch for 
cbs    changing to other orders of IREDs like e.g. in TURBOMOLE
c
c
c   HOW2ADD another symmetry:
c
c   1. add it in readbas.f to be accepted. Add the number of IRs
c
c   2. copy one of the symmetry-blocks in this subroutine and 
c      edit the multiplication-table for the group
c
c   3. assign the right IRs to L_X, L_Y and L_Z 
c
c   that is  all. Good luck!!!
c
#include "priunit.h"
#include "para.h"
#include "ired.h"
      character*3 symmetry 
      symmetry='D2H'  ! MOLCAS-Version
      if (symmetry.eq.'D2H') then 
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
      mult(5,1)=5
      mult(6,1)=6
      mult(7,1)=7
      mult(8,1)=8
c  
      mult(3,2)=4
      mult(4,2)=3
      mult(5,2)=6
      mult(6,2)=5
      mult(7,2)=8
      mult(8,2)=7
c  
      mult(4,3)=2
      mult(5,3)=7
      mult(6,3)=8
      mult(7,3)=5
      mult(8,3)=6
c  
      mult(5,4)=8
      mult(6,4)=7
      mult(7,4)=6
      mult(8,4)=5
c  
      mult(6,5)=2
      mult(7,5)=3
      mult(8,5)=4
c  
      mult(7,6)=4
      mult(8,6)=3
c  
      mult(8,7)=2
c  
C      
      do ired=1,8
      mult(ired,ired)=1
      enddo
      do irun=2,8
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
CBS   write(6,*) 
CBS   write(6,*) 
CBS  *'multiplicitation table (atkins,child and phillips)'
CBS   write(6,*) 
CBS   do ired=1,8
CBS   write(6,'(8I5)') (mult(jred,ired),jred=1,8) 
CBS   write(6,*) 
CBS   enddo
      
c     
      IRLX=4
      IRLY=3
      IRLZ=2
cbs   assume same order of ireds as Atkins Child and Phillips use..
cbs   would lead to an order with 1 to 1, 2 to 2 ...
cbs   however, this is the molecule/ seward order.   
      iredorder(1)=1
      iredorder(2)=4
      iredorder(3)=6
      iredorder(4)=7
      iredorder(5)=8
      iredorder(6)=5
      iredorder(7)=3
      iredorder(8)=2
      do ired=1,8
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(5)
      ipow2ired(0,0,1)=iredorder(6)
      ipow2ired(0,1,0)=iredorder(7)
      ipow2ired(1,0,0)=iredorder(8)
c     write(6,*) 'interacting IRs '
      do ired=1,8
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'C2V') then 
cbs   1. A1 2. A2 3. B1 4. B2
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c  
      mult(3,2)=4
      mult(4,2)=3
c  
      mult(4,3)=2
C      
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,4
      write(LUPRI,'(4I5)') (mult(jred,ired),jred=1,4) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=4
      IRLY=3
      IRLZ=2
cbs   this is the molecule/ seward order.   
      iredorder(1)=1
      iredorder(2)=4
      iredorder(3)=2
      iredorder(4)=3
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(1)
      ipow2ired(0,1,0)=iredorder(4)
      ipow2ired(1,0,0)=iredorder(3)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'D2 ') then 
cbs   1. A1 2. B1 3. B2 4. B3
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c  
      mult(3,2)=4
      mult(4,2)=3
      mult(4,3)=2
C      
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,4
      write(LUPRI,'(4I5)') (mult(jred,ired),jred=1,4) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=4
      IRLY=3
      IRLZ=2
      iredorder(1)=1
      iredorder(2)=2
      iredorder(3)=3
      iredorder(4)=4
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(1)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(3)
      ipow2ired(1,0,0)=iredorder(4)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'C2H') then 
cbs   assume 1.Ag 2.Au 3.Bg 4.Bu 
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c  
      mult(3,2)=4
      mult(4,2)=3
c  
      mult(4,3)=2
C      
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,4
      write(LUPRI,'(4I5)') (mult(jred,ired),jred=1,4) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=3
      IRLY=3
      IRLZ=1
      iredorder(1)=1
      iredorder(2)=2
      iredorder(3)=3
      iredorder(4)=4
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(1)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(3)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(4)
      ipow2ired(1,0,0)=iredorder(4)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'CS ') then 
      write(LUPRI,*) 'CS in initired '
cbs   assume 1.A' 2.A'                       
      mult(2,1)=2
C      
      do ired=1,2
      mult(ired,ired)=1
      enddo
      do irun=2,2
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,2
      write(LUPRI,'(2I5)') (mult(jred,ired),jred=1,2) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=2
      IRLY=2
      IRLZ=1
      iredorder(1)=1
      iredorder(2)=2
      do ired=1,2
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(1)
      ipow2ired(1,0,1)=iredorder(2)
      ipow2ired(0,1,1)=iredorder(2)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(1)
      ipow2ired(1,0,0)=iredorder(1)
c     write(6,*) 'interacting IRs '
      do ired=1,2
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      endif   
      return 
      end 
      subroutine kindiag(TKIN,TKINTRIA,ndim,evec,eval,breit)
#include "implicit.h"
cbs   determines eigenvectors and -values of TKIN  
      dimension tkin(ndim,ndim),
     *TKINTRIA((ndim*ndim+ndim)/2),eval(ndim),evec(ndim,ndim)
      logical breit
cbs   move symmetric matrix to triangular matrix 
      itria=1
      do irun2=1,ndim
      do irun1=1,irun2 
      TKINTRIA(itria)=TKIN(irun1,irun2)
      itria=itria+1 
      enddo
      enddo
      do irun2=1,ndim
      do irun1=1,ndim
      evec(irun1,irun2)=0d0
      enddo
      enddo
      do irun1=1,ndim
      evec(irun1,irun1)=1d0
      enddo
cbs   now diagonalize  
            CALL jacobi(TKINTRIA,evec,ndim,ndim)    
cbs   get the eigenvalues   
      do irun=1,ndim
      eval(irun)=TKINTRIA((irun*irun+irun)/2)
      enddo
      if (breit) then
      do irun=1,ndim
      eval(irun)=0d0 
      enddo
      endif 
cbs   ensure normalization of the vectors. 
      do IRUN=1,ndim
      fact=0d0
      do JRUN=1,ndim 
      fact=fact+evec(JRUN,IRUN)*evec(JRUN,IRUN) 
      enddo
      fact=1d0/dsqrt(fact)
      do JRUN=1,ndim
      evec(JRUN,IRUN)=fact*evec(JRUN,IRUN)
      enddo 
      enddo
      return   
      end   
      Subroutine kinemat(L,ndim,evtkin,type1,type2,Energy)
#include "implicit.h"
#include "alphac.h"
cbs   at least it's identical with Odd's valuE
CMI .... these values are defined in alphac.h, and are calculated at the start of AMFI !
CMI   parameter (speed2=CVEL*CVEL) 
CMI   parameter (speed4=speed2*speed2) 

cbs   this routine generates the kinematic A-factors=dsqrt((E+mc^2)/(2E))  
cbs   (type1) and   c*A/(E+mc^2) (type2)
cbs   The c in the second kinematic factor comes from Jan Almloef and 
cbs   Odd Gropen in Rev in Comp.Chem. 8(1996)
      dimension evtkin(*),type1(*),type2(*),Energy(*)  
c     E= dsqrt(p**2 c**2 + m**2 c**4) 
c     p**2= 2*m*TKIN    
c     with m = 1 
      do Irun=1,ndim
      if (evtkin(Irun).lt.0) CALL QUIT('strange kinetic energy ')
      Energy(Irun)=(evtkin(Irun)+evtkin(Irun))*speed2+speed4 
      enddo
      do Irun=1,ndim
      Energy(Irun)=dsqrt(energy(irun))
      enddo
      do Irun=1,ndim
!     dsqrt((E+mc^2)/(2E)):
      type1(Irun)=dsqrt(0.5d0*(1d0+speed2/Energy(Irun)))
      enddo
!      c*A/(E+mc^2) 
      do Irun=1,ndim
      type2(Irun)=CVEL*type1(Irun)/(Energy(Irun)+speed2)
      enddo
              do Irun=1,ndim
              type2(Irun)=2*CVEL*type2(Irun)
              enddo                        
      return 
      end
      Double precision function LMdepang(
     *L,M,l1,l2,l3,l4,m1,m2,m3,m4,cheater)
cbs   l1-l4 and m1-m4 are already shifted !!
cbs   purpose: calculates the angular part of the   
cbs   coulomb-type integrals. See documentation for details...
cbs   LMdepang= LM dependent angular factors 
cbs   cheater included for a correcting signs, as there were some 
cbs   signs (only signs!!!!) missing when compared to HERMIT  
cbs                                        B.S.  08.10.96 
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
      LMdepang=0d0
cbs   some quick checks
      if (L.lt.abs(M)) return 
      if (l1.lt.abs(m1)) return 
      if (l2.lt.abs(m2)) return 
      if (l3.lt.abs(m3)) return 
      if (l4.lt.abs(m4)) return 
cbs   prefactor
      fact1=4d0*pi/dfloat(L+L+1)
cbs   determining the sign
      isum=-l3-l1-l4-l2+2*(M+m3+m4)   !???? I am not sure 
      if (mod(isum,4).eq.0) then 
      isign=1
      elseif (iabs(mod(isum,4)).eq.2) then 
      isign=-1
      else 
      write(LUPRI,*) 'L,l1,l2,l3,l4,M,m1,m2,m3,m4'
      write(LUPRI,'(10I3)') L,l1,l2,l3,l4,M,m1,m2,m3,m4
      write(LUPRI,*) 'isum= ',isum,' mod = ',mod(isum,4)
      CALL QUIT('error in lmdepang')
      endif
      fact2=couple3J(L,l3,l1,-M,m3,-m1) 
      fact3=couple3J(L,l4,l2,M,m4,-m2)
C     write(6,*) 'fact2,fact3 ',fact2,fact3
      LMdepang=cheater*dfloat(isign)*fact1*fact2*fact3
      return 
      end 
      logical function mcheckxy(m1,m2,m3,m4)
      integer m1,m2,m3,m4,int12a,int12b,
     *int34a,int34b
cbs   makes a check, if there is an interaction inbetween cartesian functions 
cbs   with m-values m1-m4 
      mcheckxy=.true.   
      int12a=m1+m2
      int12b=-m1+m2
      int34a=m3+m4
      int34b=-m3+m4
cbs   lots of checks 
      if (iabs(int12a+int34a).eq.1) return
      if (iabs(int12a-int34a).eq.1) return
      if (iabs(int12b+int34b).eq.1) return
      if (iabs(int12b-int34b).eq.1) return
      if (iabs(int12a+int34b).eq.1) return
      if (iabs(int12a-int34b).eq.1) return
      if (iabs(int12b+int34a).eq.1) return
      if (iabs(int12b-int34a).eq.1) return
      mcheckxy=.false.
      return 
      end 
      logical function mcheckz(m1,m2,m3,m4)
cbs   makes a check, if there is an interaction inbetween cartesian functions 
cbs   with m-values m1-m4 
      integer m1,m2,m3,m4,int12a,int12b,
     *int34a,int34b
      mcheckz=.true.   
      int12a=m1+m2
      int12b=-m1+m2
      int34a=m3+m4
      int34b=-m3+m4
cbs   lots of checks 
      if (iabs(int12a+int34a).eq.0) return
      if (iabs(int12a-int34a).eq.0) return
      if (iabs(int12b+int34b).eq.0) return
      if (iabs(int12b-int34b).eq.0) return
      if (iabs(int12a+int34b).eq.0) return
      if (iabs(int12a-int34b).eq.0) return
      if (iabs(int12b+int34a).eq.0) return
      if (iabs(int12b-int34a).eq.0) return
      mcheckz=.false.
      return 
      end 
      subroutine mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *angintSO,angintOO,
     *Lfirst,Llast,Lblocks,
     *ncont1,ncont2,ncont3,
     *ncont4,
     *caseaSO,caseb1SO,caseb2SO,casecSO,
     *caseaOO,caseb1OO,caseb2OO,casecOO,
     *preroots,clebsch,dummy,bonn,breit,
     *sameorb)
#include "implicit.h"
cbs   subroutine for combining radial integrals with angular 
cbs   factors for the block with l1,l2,l3,l4,m1,m2,m3m,m4  
cbs   this routine mkangL0 = make angular factors for the L0-part 
cbs   includes both, spin-same and spin-other-orbit parts. 
      double precision LMdepang
      dimension 
     *angintSO(ncont1,ncont2,ncont3,ncont4),
     *angintOO(ncont1,ncont2,ncont3,ncont4),
     *Lfirst(*),Llast(*),Lblocks(*),
cbs   all the arrays with the radial integrals for 
cbs   this combination of l-values   
     *caseaSO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   integrals with alpha1*alpha3
     *caseb1SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha1
     *caseb2SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha3
     *casecSO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  integrals with factor 1
     *caseaOO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   integrals with alpha1*alpha3
     *caseb1OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha1
     *caseb2OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha3
     *casecOO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  integrals with factor 1
     *preroots(2,0:Lmax),                    ! some prefactors: dsqrt( (l(+1))/(2l+1))
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)         ! some clebsch gordans, that appear regulary
      dimension dummy(0:*)                                    
      logical bonn,breiT,sameorb 
c     write(6,*) 'begin mkangL0 ',
c    *l1,l2,l3,l4,m1,m2,m3,m4
cbs  
      ncontall=ncont1*ncont2*ncont3*ncont4
cbs   cheater introduced to correct signs, because they were different from HERMIT 
      if (mod(l1+l2+l3+l4,4).eq.2) then                  
      cheater=1d0
      else 
      cheater=-1d0
      endif 
cbs   cleaning up 
      if (bonn.or.breit.or.sameorb) then 
      call dzero(angintSO,ncontall) 
      else 
      call dzero(angintSO,ncontall) 
      call dzero(angintOO,ncontall) 
      endif  
cbs  starting with the same-orbit-contributions 
cbs  first term: ###########################################################################
      factor=-preroots(2,l1)*preroots(2,l3)*
     *clebsch(1,2,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1-1,m2,m3-1,m4,cheater) 
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then 
         Call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),CaseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   second term: ###########################################################################
      factor=-preroots(1,l1)*preroots(2,l3)*
     *clebsch(1,1,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1-1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,AngintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0) then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   third term: ###########################################################################
      factor=-preroots(2,l1)*preroots(1,l3)*
     *clebsch(1,2,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1-1,m2,
     *m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),CaseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l3)*factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   fourth term: ###########################################################################
      factor=-preroots(1,l1)*preroots(1,l3)*
     *clebsch(1,1,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,AngintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l3)*factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l1)*factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
	      if (Lfirst(4).lt.Kfirst) then 
	      do L=Lfirst(4),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(4)
	      endif 
	      if (Llast(4).gt.Klast) then 
	      do L=Klast,Llast(4),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(4)   
	      endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0)  then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   (4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs  fifth term: ###########################################################################
      factor=preroots(2,l1)*preroots(2,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(3,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
      else  
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
cbs   sixth  term: ###########################################################################
      factor=preroots(1,l1)*preroots(2,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(3,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0)  then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   seventh term: ###########################################################################
      factor=preroots(2,l1)*preroots(1,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(3,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         Call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         Call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         Call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   eigth term: ###########################################################################
      factor=preroots(1,l1)*preroots(1,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(3,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0) then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l3)*factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
	      if (Lfirst(4).lt.Kfirst) then 
	      do L=Lfirst(4),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(4)
	      endif 
	      if (Llast(4).gt.Klast) then 
	      do L=Klast,Llast(4),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
	      enddo 
	      Klast=Llast(4)   
	      endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   (4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
      return   
      end   
      subroutine mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,             
     *angintSO,angintOO,
     *Lfirst,Llast,Lblocks,
     *ncont1,ncont2,ncont3,
     *ncont4,
     *caseaSO,caseb1SO,caseb2SO,casecSO,
     *caseaOO,caseb1OO,caseb2OO,casecOO,
     *preroots,clebsch,dummy,bonn,breit,
     *sameorb)
#include "implicit.h"
cbs   subroutine for combining radial intgrls with angular 
cbs   factors for the block with l1,l2,l3,l4,m1,m2,m3m,m4  
cbs   this routine mkangLmin = make angular factors for the L- -part 
cbs   includes both, spin-same and spin-other-orbit parts. 
      double precision LMdepang
      dimension 
     *angintSO(ncont1,ncont2,ncont3,ncont4),
     *angintOO(ncont1,ncont2,ncont3,ncont4),
     *Lfirst(*),Llast(*),Lblocks(*),
cbs   all the arrays with the radial intgrls for 
cbs   this combination of l-values   
     *caseaSO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   intgrls with alpha1*alpha3
     *caseb1SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha1
     *caseb2SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha3
     *casecSO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  intgrls with factor 1          
     *caseaOO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   intgrls with alpha1*alpha3
     *caseb1OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha1
     *caseb2OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha3
     *casecOO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  intgrls with factor 1          
     *preroots(2,0:Lmax),                    ! some prefactors: dsqrt( (l(+1))/(2l+1)) 
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)         ! some clebsch gordans, that appear regulary 
      dimension dummy(0:*)                                   
      logical bonn,breiT,sameorb 
      root2=dsqrt(2.0d0) 
      root2inv=1d0/root2            
c     write(6,*) 'begin mkangL- ',
c    *l1,l2,l3,l4,m1,m2,m3,m4
cbs  
      ncontall=ncont1*ncont2*ncont3*ncont4
cbs   cheater introduced to correct signs, because they were different from HERMIT
      if (mod(l1+l2+l3+l4,4).eq.2) then
      cheater=1d0
      else
      cheater=-1d0
      endiF
cbs   cleaning up 
      if (bonn.or.breit.or.sameorb) then 
      call dzero(angintSO,ncontall)
      else 
      call dzero(angintSO,ncontall)
      call dzero(angintOO,ncontall)
      endif  
cbs  starting with the same-orbit-contributions 
cbs  first term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(2,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(2,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      Endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   second term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(2,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(2,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater) 
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
      if (Lfirst(3).lt.Kfirst) then 
      do L=Lfirst(3),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater) 
      enddo 
      Kfirst=Lfirst(3)
      endif 
      if (Llast(3).gt.Klast) then 
      do L=Klast,Llast(3),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Klast=Llast(3)
      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   third term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(1,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(2,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
      if (Lfirst(2).lt.Kfirst) then 
      do L=Lfirst(2),Kfirst,2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,
     *                 m3,m4,Cheater)
      enddo 
      Kfirst=Lfirst(2)
      endif 
      if (Llast(2).gt.Klast) then 
      do L=Klast,Llast(2),2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Klast=Llast(2)
      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else   
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   fourth term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(1,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(2,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
      if (Lfirst(2).lt.Kfirst) then 
      do L=Lfirst(2),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Kfirst=Lfirst(2)
      endif 
      if (Llast(2).gt.Klast) then 
      do L=Klast,Llast(2),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Klast=Llast(2)
      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif  
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
      if (Lfirst(3).lt.Kfirst) then 
      do L=Lfirst(3),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Kfirst=Lfirst(3)
      endif 
      if (Llast(3).gt.Klast) then 
      do L=Klast,Llast(3),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Klast=Llast(3)
      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
      if (Lfirst(4).lt.Kfirst) then 
      do L=Lfirst(4),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Kfirst=Lfirst(4)
      endif 
      if (Llast(4).gt.Klast) then 
      do L=Klast,Llast(4),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      enddo 
      Klast=Llast(4)   
      endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else  
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   (4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs  fifth term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(2,l3)*
     *clebsch(2,2,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
cbs   sixth  term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(2,l3)*
     *clebsch(2,1,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,4*
     *   factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
      if (Lfirst(3).lt.Kfirst) then 
      do L=Lfirst(3),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Kfirst=Lfirst(3)
      endif 
      if (Llast(3).gt.Klast) then 
      do L=Klast,Llast(3),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Klast=Llast(3)
      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   seventh term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(1,l3)*
     *clebsch(2,2,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
      if (Lfirst(2).lt.Kfirst) then 
      do L=Lfirst(2),Kfirst,2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Kfirst=Lfirst(2)
      endif 
      if (Llast(2).gt.Klast) then 
      do L=Klast,Llast(2),2
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Klast=Llast(2)
      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   eigth term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(1,l3)*
     *clebsch(2,1,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,4*
     *   factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
      if (Lfirst(2).lt.Kfirst) then 
      do L=Lfirst(2),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Kfirst=Lfirst(2)
      endif 
      if (Llast(2).gt.Klast) then 
      do L=Klast,Llast(2),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Klast=Llast(2)
      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
        call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *  caseb1SO(1,Lrun),1,angintSO,1)
      else 
        call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *  caseb1SO(1,Lrun),1,angintSO,1)
        call daxpy(ncontall,-(2+4*l3)*
     *factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
      if (Lfirst(3).lt.Kfirst) then 
      do L=Lfirst(3),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Kfirst=Lfirst(3)
      endif 
      if (Llast(3).gt.Klast) then 
      do L=Klast,Llast(3),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Klast=Llast(3)
      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
      if (Lfirst(4).lt.Kfirst) then 
      do L=Lfirst(4),Kfirst,2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Kfirst=Lfirst(4)
      endif 
      if (Llast(4).gt.Klast) then 
      do L=Klast,Llast(4),2
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      enddo 
      Klast=Llast(4)   
      endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else   
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   (4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
      return   
      end
   
      subroutine prefac(Lmax,preroots,clebsch)
#include "implicit.h"
      dimension preroots(2,0:Lmax),
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)
cbs   the roots appearing in front of all 
cbs   the contributions 
c     write(6,*) 'begin of prefac'
      do L=0,Lmax
      fact=1d0/dsqrt(dfloat(L+L+1))
      preroots(1,L)=dsqrt(dfloat(L))*fact 
      preroots(2,L)=dsqrt(dfloat(L+1))*fact 
      enddo
cbs   there are Clebsch-Gordon-Coefficients 
cbs   which always appear:
cbs 
cbs   -----                       ------
cbs  |                                 |
cbs  |  l +/- 1     1        |      l  |
cbs  |                       |         |
cbs  |                       |         |   
cbs  |  m+/-1,0   -1,1,0     |      m  |   
cbs  |                       |         |   
cbs  |                                 |   
cbs   -----                       -----
cbs 
cbs 
cbs  array clebsch (3,2,-Lmax:Lmax,0:Lmax)
cbs  first index    1:  m-1 
cbs                 2:  m 
cbs                 3:  m+1
cbs  second index   1:  l-1
cbs                 2:  l+1
cbs  third index        m 
cbs  fourth index       l 
cbs  
c     write(6,*),'start to generate CGs' 
      do L=0,Lmax
      L2=L+L
      do M=-L,L
c     write(6,*) 'L,M: ',L,M 
      M2=M+M
cbs   getCG calculates CG-coeffecients. In order to avoid fractions,
cbs   e.g. for spins, arguments are doubled values...
      clebsch(1,1,M,L)=
     *getCG(L2-2,2,L2,M2-2,2,M2)
      clebsch(2,1,M,L)=
     *getCG(L2-2,2,L2,M2,0,M2)
      clebsch(3,1,M,L)=
     *getCG(L2-2,2,L2,M2+2,-2,M2)
      clebsch(1,2,M,L)=
     *getCG(L2+2,2,L2,M2-2,2,M2)
      clebsch(2,2,M,L)=
     *getCG(L2+2,2,L2,M2,0,M2)
      clebsch(3,2,M,L)=
     *getCG(L2+2,2,L2,M2+2,-2,M2)
      enddo
      enddo 
      return 
      end
