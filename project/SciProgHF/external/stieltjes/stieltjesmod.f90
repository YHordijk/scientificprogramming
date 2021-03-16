MODULE StieltjesMod

!  USE General
!  USE Configs, ONLY : e_init

  IMPLICIT NONE

  REAL*16,PRIVATE      :: offset = 0.5 ! the distance from lowest energy to 0
  REAL(16),PRIVATE     :: overmax = 10. ! 1/overmax -> allowing 10% overlap for polynomials
  REAL(16),PRIVATE    :: conv_thresh0 = 5.e-2
  REAL(16),PRIVATE     :: conv_fac = 1. !> 1 higher moments,<1 lower moments preferred, handle by SIdetect
  INTEGER,PRIVATE      :: conv_max = 10 ! allow more deviation from orthogonality

CONTAINS

!  SUBROUTINE init_stieltjes(ioffset,iovermax,iconv_thresh0,iconv_fac &
!                ,iconv_max)
!    REAL(idk),INTENT(in) :: ioffset,iovermax,iconv_thresh0,iconv_fac
!    INTEGER,INTENT(in) :: iconv_max
!    !***************************************************************
!    offset = real(ioffset,kind=16)
!    overmax = real(iovermax,kind=16)
!    conv_thresh0 = iconv_thresh0
!    conv_fac = iconv_fac
!    conv_max = iconv_max
!  END SUBROUTINE init_stieltjes

  SUBROUTINE stieltjes(num,e8,g8,e_init,outfile,gamma0,unitout,printflag)
! Vitali Averbukh (2003). Comments to: vitali@tc.pci.uni-heidelberg.de
! This routine receives the sequence of the energy-gamma pairs and transforms it 
! into a shorter sequence of the Stieltjes energy-gamma pairs using the algorithm 
! of Muller-Plathe & Dierksen [Electronic Structure of Atoms, Molecules and Solids, 
! Proceedings of the II Escola Brasileira de Estruture Eletronica (World Scientific, 
! Singapore, 1990), p.1-29]. The tutorial is available on the net: 
!                              http://citeseer.nj.nec.com/465113.html
! Sequences of more than three points and up to NP points are allowed. 
! The recursive relations for the coefficients of the orthogonal polynomials
! are implemented in quadruple precision. The maximal polynomial order is determined 
! by the magnitude of the polynomial overlap approacing the polynomial norm.
! If only low-order (MAXORD < 5) approximation is available, the program produces the 
! interpolated Gamma(E) function. Otherwise, a sequence of Gamma(E) functions for the 
! approximation orders from 5 to MAXORD so that the convergence can be seen.
! The tridiagonal coefficient matrix is diagonalized by the REAL*16 version of the TQL2 
! routine (source available).
! The cumulative gamma is differentiated numerically and the resulting set of points is 
! interpolated using monotonicity-preserving piecewise cubic Hermite interpolant 
! (NAG: E01BEF,E01BFF). This interpolation scheme ensures that the Gamma(E) 
! is non-negative. The gamma at the ENERGY(1) is calculated by successively 
! higher order Stieltjes approximation. By default, the average of the three highest-order 
! results is taken as the final one. If convergence is detected at some lower order, the 
! converged result is preferred.
!
! The following parameters are transferred from the main program:
!*
! CONV_MAX - the maximal number of convergence search loops
!*
! OFFSET -  the energy shift (should not affect the final result)
!*
! OVERMAX -  the maximal permitted norm-to-overlap ratio of the 
!            adjacent orthogonal polynomials
!*
! CONV_THRESH0 - initial convergence threshold
!*
! CONV_FAC - convergence threshold factor making the convergence requirement
!            more strict for lower Stieltjes orders

      use stieltjes_ql_diag

      real*8 e8(:),g8(:)
      REAL(8), INTENT(IN)       :: e_init
      CHARACTER(11), INTENT(IN) :: outfile

      integer num,nmax,npol
      integer np,maxord,ierr,imax,min,max,unitout,printflag
      integer i,j,ifail,converge,iconv

! the maximal number of the energy points is NP
! the maximal order of polynomial in reduced moment problem is NMAX
      parameter (nmax=100)

      real*16 e_point(size(e8)),g_point(size(e8)) 
      real*16 e_new(nmax),g_new(nmax)
      real*16 bcoef(0:nmax),acoef(nmax),qpol(0:nmax,size(e8))
      real*16 diag(nmax),offdiag(nmax),abvec(nmax,nmax)
      real*16 e_min,e_max,bprod,asum,qnorm,qoverlap

      real*8 e0(nmax),gamma(nmax),der(nmax)
      real*8 energy(1),gammas(1),gamma0,gamma_h(nmax)
      real*8 dif12,dif23,dif13,gmean,difmean
      real*8 dfm1
      real*8 conv_thresh

!DEC$ IF DEFINED (DBGSTIELTJES)
      integer unitdbg
!DEC$ ENDIF

      np = size(e8)

! store the Stieltjes input
!DEC$ IF DEFINED (DBGSTIELTJES)
      unitdbg=unitout+3000
      write(unitdbg,*) 'num, unitout, printflag'
      write(unitdbg,*) num,unitout,printflag
      write(unitdbg,*) 'gamma0'
      write(unitdbg,*) gamma0
      write(unitdbg,*) 'offset,overmax'
      write(unitdbg,*) offset,overmax
      write(unitdbg,*) 'conv_thresh0,conv_fac'
      write(unitdbg,*) conv_thresh0,conv_fac
      write(unitdbg,*) 'conv_max'
      write(unitdbg,*) conv_max
      write(unitdbg,*) 'np = ',np
      write(unitdbg,*)
      write(unitdbg,*) 'energy points'
      do i = 1,num
        write(unitdbg,*) e8(i),g8(i)
      end do
!DEC$ ENDIF

! check the number of the input energy points 
      if (num.gt.size(e8)) then
         write(6,'("***WARNING*** Stieltjes: too many energy points")')
         write(6,'("          NUM (",i6,") > ",i6)') num,np
         gamma0=0.d0
         return
      end if
      if (num.le.3) then
         write(6,'("***ERROR*** Stieltjes: not enough energy points")')
         write(6,'("          NUM = ",i3)') num
         gamma0=0.d0
         return
      end if

! transform the energy and gamma values into quadruple precision
      if (printflag.eq.2) write(6,*) 'Input energy points'
      do i=1,num
         e_point(i)=real(e8(i),kind=16)
         g_point(i)=real(g8(i),kind=16)
         if (printflag.eq.2) write(6,*) e8(i),g8(i),i
      end do
      if (printflag.eq.2) write(6,*)

! define the minimal energy and the maximal energies
      e_min=e_point(1)
      e_max=e_point(1)
      do i=2,num
         if( e_min.gt.e_point(i)) e_min=e_point(i)
         if( e_max.lt.e_point(i)) e_max=e_point(i)
      end do

! define the energy of interest: gamma(energy) = Siegert gamma
      energy(1) = e_init

! check whether the defined energy belongs to the interval covered by E_POINT
      if (energy(1).lt.e_min.or.energy(1).gt.e_max) then
         gmean=0.d0
         do i=1,num
            gmean=gmean+g8(i)
         end do
         gmean=gmean/dfloat(num)
         write(6,'("***ERROR*** Stieltjes:")')
         write(6,'(" the required energy is out of range")')
         write(6,'("          ",e16.10,"  is not in [",e16.10,",",e16.10,"]" &
           )') energy(1),e_min,e_max
         write(6,'(" the number of energy points is ",i6)') num
         write(6,'(" the mean gamma in the covered range is ",es16.10)') gmean
         write(6,*)
         gamma0=0.d0
         return
      end if

! shifting the energies to avoid the small denumerator problem
      if (offset.gt.0.0q0) then
        do i=1,num
          e_point(i)=e_point(i)-e_min+offset
        end do
      end if

! initiate the recursive computation of the a,b coefficients and the orthogonal 
! polynomials according to (3.3.20-23) of Mueller-Plathe & Dierksen (1990)
       bcoef(0)=0.q0
       acoef(1)=0.q0
       do i=1,num
          bcoef(0)=bcoef(0)+g_point(i)
          acoef(1)=acoef(1)+g_point(i)/e_point(i)
       end do
       acoef(1)=acoef(1)/bcoef(0)

       do i=1,num
          qpol(0,i)=1.q0
          qpol(1,i)=1.q0/e_point(i)-acoef(1)
       end do

       bcoef(1)=0.q0
       acoef(2)=0.q0
       do i=1,num
          bcoef(1)=bcoef(1)+qpol(1,i)*g_point(i)/e_point(i)
          acoef(2)=acoef(2)+qpol(1,i)*g_point(i)/(e_point(i)**2)
       end do
       bcoef(1)= bcoef(1)/bcoef(0)
       acoef(2)=acoef(2)/(bcoef(0)*bcoef(1))-acoef(1)

! calculate the higher-order coefficients and polynomials recursively
! up to the (NUM-1)th order (total of NUM polynomials)
! if NUM > NMAX, we calculate only NMAX polynomials
       if (num.le.nmax) then
          npol=num
       else
          npol=nmax
       end if

       asum=acoef(1)
       do i=3,npol

          asum=asum+acoef(i-1)

          do j=1,num
             qpol(i-1,j)=(1.q0/e_point(j)-acoef(i-1))*qpol(i-2,j)   &
                  -bcoef(i-2)*qpol(i-3,j)
          end do

          bprod=bcoef(0)
          do j=1,i-2
             bprod=bprod*bcoef(j)
          end do

          bcoef(i-1)=0.q0
          do j=1,num
             bcoef(i-1)=bcoef(i-1)+qpol(i-1,j)*g_point(j)   &
                  /(e_point(j)**(i-1))
          end do
          bcoef(i-1)=bcoef(i-1)/bprod

          bprod=bprod*bcoef(i-1)

          acoef(i)=0.q0
          do j=1,num
             acoef(i)=acoef(i)+qpol(i-1,j)*g_point(j)/(e_point(j)**i)
          end do
          acoef(i)=acoef(i)/bprod-asum

       end do

! calculate the NUM-th (NPOL-th) order polynomial just for the orthogonality check 
       do j=1,npol
          qpol(npol,j)=(1.q0/e_point(j)-acoef(npol))*qpol(npol-1,j) &
               -bcoef(npol-1)*qpol(npol-2,j)
       end do

! check the orthogonality of the polynomials to define the maximal approximation order 
! if the orthogonality is preserved for all orders, MAXORD is set to NPOL
       maxord=npol
       qnorm=bcoef(0)
       if (printflag.ne.0) write(6,*) 'Stieltjes: orthogonality check'
       do i=1,npol
          qnorm=0.q0
          qoverlap=0.q0
          do j=1,num 
             qnorm=qnorm+qpol(i,j)**2*g_point(j)
             qoverlap=qoverlap+qpol(i,j)*qpol(i-1,j)*g_point(j)
          end do
          if (abs(qoverlap).lt.1.q-50) qoverlap=1.q-50
          !if (qabs(qoverlap).lt.1.q-50) qoverlap=1.q-50
          if (printflag.ne.0) write(6,*) i,qoverlap,qnorm
          !if (qnorm/qabs(qoverlap).le.overmax) then
          if (qnorm/abs(qoverlap).le.overmax) then
! MAXORD=I-1 is appropriate since the polynomial failing 
! the orthogonality check should not be used
             maxord=i-1
             go to 10
          end if
       end do

 10    continue
       if (printflag.ne.0) write(6,*)

! look how many Stieltjes orders are available
       if (maxord.lt.5) then
          min=maxord
          max=maxord
          print*, '***WARNING*** Stieltjes:' 
          print*, ' only very low-order approximation is available'
          print*, ' MAXORD=',maxord
       else
          min=5
          max=maxord
       end if

! seek for convergence if enough Stieltjes orders are available
       converge=0
       if (maxord.ge.7) converge=1

! perform the gamma calculation using the successive approximations 
! N=5,...,MAXORD (if MAXORD > 5)
       open(150,file=outfile)
       do 20 imax=min,max
       
       write(150,*) imax

! fill the coefficients matrix
       do i=1,imax
          diag(i)=acoef(i)
       end do
       do i=2,imax
          offdiag(i)=-sqrt(bcoef(i-1))
          !offdiag(i)=-qsqrt(bcoef(i-1))
!          WRITE(*,*) offdiag(i)
       end do

! diagonalize the coefficients matrix
! initialize the arrays
       do i=1,nmax
          do j=1,nmax
             abvec(i,j)=0.q0
          end do
          abvec(i,i)=1.q0
       end do
       call tql2(nmax,imax,diag,offdiag,abvec,ierr)
       if (ierr.ne.0) then
          print*, '***WARNING*** Stieltjes:'
          print*, ' the eigenvalue no. ',ierr,' failed to converge'
       end if
! fill the Stieltjes energy and gamma arrays
! note that the eigenvalues are inverse energies and are given in ascending order 
       do i=1,imax
          e_new(i)=1.q0/diag(imax+1-i)
          g_new(i)=bcoef(0)*abvec(1,imax+1-i)**2
       end do

! calculate the gamma's by simple numerical differentiation at the middle 
! point of each [E_NEW(I),E_NEW(I+1)] interval
       do i=1,imax-1
          e0(i)=(e_new(i)+e_new(i+1))/2.d0
          if (offset.gt.0.q0) then
            e0(i) = e0(i)+e_min-offset
          end if
          gamma(i)=5.d-1*(g_new(i+1)+g_new(i))/(e_new(i+1)-e_new(i))
          print*, e0(i),gamma(i),imax
          write(150,*) e0(i),gamma(i)
       end do
       print*, ' ' 

! check whether the required energy is inside the interval covered by the E0 grid,
! the Stieltjes orders for which this is not the case are excluded from the 
! convergence search
      if (energy(1).lt.e0(1)) then 
        write(6,'("***WARNING*** Stieltjes:")')
        write(6,'("          ",es16.10," < ",es16.10," (the first grid point)")') &
          energy(1),e0(1)
        write(6,'(" large inaccuracy expected")')
        gammas(1)=5.d-1*g_new(1)/e_new(1)
        print*, imax,gammas(1)
        gamma_h(imax)=gammas(1)
        min=min+1
        go to 20
      end if
      if (energy(1).gt.e0(imax-1)) then 
        write(6,'("***WARNING*** Stieltjes:")')
        write(6,'("          ",es16.10," < ",es16.10," (the last grid point)")')  &
          energy(1),e0(imax-1)
        gammas(1)=0.d0
        print*, imax,gammas(1)
        gamma_h(imax)=gammas(1)
        min=min+1
        go to 20
      end if

! and interpolate the result using monotonicity-preserving piecewise 
! cubic Hermite interpolant 
!       call e01bef (imax-1,e0,gamma,der,ifail)
!       call e01bff (imax-1,e0,gamma,der,1,energy,gammas,ifail)
       call herm_spline_init(imax-1,e0,gamma,der)
       call herm_spline(imax-1,e0,gamma,der,1,energy,gammas,ifail)
       if (ifail.eq.6) then
         write(6,*) "*** Warning, an error in the spline interpolation."
       else if (ifail.lt.0) then
         write(6,*) "*** Warning, the requested value had to be"
         write(6,*) "    extrapolated in the cubic spline."
       end if
       gamma_h(imax)=gammas(1)
       print*, imax,gammas(1)

! print the gamma as a function of the Stieltjes order to a separate file
       if (unitout.ne.0) write(unitout,*) imax,gammas(1)

 20   continue

      close(150)

      print*, ' '

! if the energy of interest is not covered in too many Stieltjes orders, no 
! convergence check will be initiated
      if (min+2.ge.max) converge=0

! if not enough gamma values are available, the highest-order approximation is taken
      if (converge.eq.0) then
         gamma0=gamma_h(max)

! if enough gamma values are available to perform the convergence check, 
! the final gamma is taken to be an average of the three highest-order approximations
! unless convergence is detected; the converged result is preferred
      else if (converge.eq.1) then
         gamma0=(gamma_h(max-2)+gamma_h(max-1)+gamma_h(max))/3.d0
         conv_thresh=conv_thresh0 
         iconv=1

 40      continue

         do imax=max,min+2,-1

            dif12=dabs(gamma_h(imax)-gamma_h(imax-1))
            dif23=dabs(gamma_h(imax-1)-gamma_h(imax-2))
            dif13=dabs(gamma_h(imax)-gamma_h(imax-2))
            gmean=(gamma_h(imax-2)+gamma_h(imax-1)+gamma_h(imax))/3.d0
            dfm1=(dif12+dif23+dif13)/3.d0
            difmean=(dif12+dif23+dif13)/3.d0

! prefer the convergence at higher orders over the one at lower orders
! to this end, the convergence criterion becomes more strict as we go to lower orders
            if (difmean.le.gmean*conv_thresh*conv_fac**(max-imax)) then
              write(6,'(" Convergence detected at ",i3,"th order")') imax
              write(6,'(" Convergence threshold = ",f5.1,"%")') &
              conv_thresh*1.d2
              gamma0=gmean
              go to 30
            end if

         end do

! if we came here, the convergence was not found, so we increase the threshold
! maximal number of CONV_MAX convergence search loops is allowed
         iconv=iconv+1

         if (iconv.eq.conv_max) go to 30
         conv_thresh=conv_thresh*1.2d0

         go to 40 

      end if

 30   continue

      print*, 'The final result: Gamma =',gamma0
      print*, ' '
      if (unitout.ne.0) write(unitout,*) ' '
  END SUBROUTINE stieltjes

END MODULE StieltjesMod
