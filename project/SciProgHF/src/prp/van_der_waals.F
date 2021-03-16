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

      BLOCK DATA gauss_legendre_initialization 
#include "gauss_legendre.h"
      DATA FREQ_6PT/
     &      0.010483552568D0, 0.061182644938D0, 0.184410387584D0,
     &      0.488041922035D0, 1.471005382177D0, 8.584876110674D0/     
      DATA FREQ_8PT/
     &      0.006077184458D0, 0.033951797700D0, 0.093305311706D0,
     &      0.206998847637D0, 0.434785029129D0, 0.964575310389D0,
     &      2.650816925655D0, 14.809489593503D0/
      DATA FREQ_10PT/
     &      0.003965760957D0, 0.021704887199D0, 0.057268418214D0, 
     &      0.118586527107D0, 0.222250327683D0, 0.404904878427D0, 
     &      0.758939503461D0, 1.571546810728D0, 4.146531570245D0, 
     &      22.694257409977D0/
      DATA FREQ_12PT/
     &      0.002791642925D0, 0.015106644820D0, 0.039001691306D0,
     &      0.077996102406D0, 0.138650521237D0, 0.233222703362D0,
     &      0.385897250578D0, 0.649114040084D0, 1.153903813439D0,
     &      2.307592234713D0, 5.957643213833D0, 32.239080141755D0/
      DATA FREQ_14PT/
     &      0.002071636175D0, 0.011133139668D0, 0.028371042124D0,
     &      0.055599195782D0, 0.095974618028D0, 0.154851318305D0,
     &      0.241489391568D0, 0.372687178576D0, 0.581202672246D0,
     &      0.937747936369D0, 1.618728449827D0, 3.172248647280D0,
     &      8.083972956642D0, 43.443921818061D0/   
      DATA FREQ_16PT/
     &      0.001598330153D0, 0.008550707934D0, 0.021606970998D0, 
     &      0.041801579820D0, 0.070856548562D0, 0.111517898281D0,
     &      0.168163496921D0, 0.247938945545D0, 0.362992589977D0,
     &      0.535193437621D0, 0.807045338795D0, 1.270171943543D0,
     &      2.153028674657D0, 4.165322386321D0, 10.525444290040D0,
     &      56.308766862655D0/ 
      DATA ABSCISSAS_6PT/
     &     -0.932469514203D0, -0.661209386466D0, -0.238619186083D0,
     &      0.238619186083D0,  0.661209386466D0,  0.932469514203D0/
      DATA ABSCISSAS_8PT/
     &     -0.960289856498D0, -0.796666477414D0, -0.525532409916D0,
     &     -0.183434642496D0,  0.183434642496D0,  0.525532409916D0,
     &      0.796666477414D0,  0.960289856498D0/ 
      DATA ABSCISSAS_10PT/
     &     -0.973906528517D0, -0.865063366689D0, -0.679409568299D0,
     &     -0.433395394129D0, -0.148874338982D0,  0.148874338982D0,
     &      0.433395394129D0,  0.679409568299D0,  0.865063366689D0,
     &      0.973906528517D0/ 
      DATA ABSCISSAS_12PT/
     &     -0.981560634247D0, -0.90411725637D0,  -0.769902674194D0, 
     &     -0.587317954287D0, -0.367831498998D0, -0.125233408511D0,
     &      0.125233408511D0,  0.367831498998D0,  0.587317954287D0,
     &      0.769902674194D0,  0.90411725637D0,    0.981560634247D0/ 
      DATA ABSCISSAS_14PT/
     &     -0.986283808697D0, -0.928434883664D0, -0.82720131507D0,
     &     -0.687292904812D0, -0.515248636358D0, -0.319112368928D0,	
     &     -0.108054948707D0,  0.108054948707D0,  0.319112368928D0,
     &      0.515248636358D0,  0.687292904812D0,  0.82720131507D0,
     &      0.928434883664D0,  0.986283808697D0/
      DATA ABSCISSAS_16PT/
     &     -0.989400934992D0, -0.944575023073D0, -0.865631202388D0,
     &     -0.755404408355D0, -0.617876244403D0, -0.458016777657D0,
     &     -0.281603550779D0, -0.095012509837D0,  0.095012509837D0,
     &      0.281603550779D0,  0.458016777657D0,  0.617876244403D0,
     &      0.755404408355D0,  0.865631202388D0,  0.944575023073D0,
     &      0.989400934992D0/                 	
      DATA WEIGHTS_6PT/
     &     0.171324492379D0,   0.360761573048D0,  0.467913934573D0,
     &     0.467913934573D0,   0.360761573048D0,  0.171324492379D0/
      DATA WEIGHTS_8PT/
     &     0.10122853629D0,    0.222381034453D0,  0.313706645878D0,
     &     0.362683783378D0,   0.362683783378D0,  0.313706645878D0,
     &     0.222381034453D0,   0.10122853629D0/
      DATA WEIGHTS_10PT/
     &     0.066671344308D0,   0.149451349151D0,  0.219086362516D0,
     &     0.26926671931D0,    0.295524224715D0,  0.295524224715D0,
     &     0.26926671931D0,    0.219086362516D0,  0.149451349151D0,
     &     0.066671344308D0/
      DATA WEIGHTS_12PT/
     &     0.047175336386D0,   0.106939325995D0, 0.160078328543D0,
     &     0.203167426723D0,   0.233492536538D0, 0.249147045813D0,
     &     0.249147045813D0,   0.233492536538D0, 0.203167426723D0,
     &     0.160078328543D0,   0.106939325995D0, 0.047175336386D0/
      DATA WEIGHTS_14PT/
     &     0.035119460332D0,   0.080158087159D0, 0.121518570688D0,
     &     0.157203167158D0,   0.185538397478D0, 0.205198463721D0,
     &     0.215263853463D0,   0.215263853463D0, 0.205198463721D0,
     &     0.185538397478D0,   0.157203167158D0, 0.121518570688D0,    
     &     0.080158087159D0,   0.035119460332D0/
      DATA WEIGHTS_16PT/
     &     0.027152459411D0,   0.062253523937D0, 0.095158511683D0,
     &     0.124628971256D0,   0.149595988817D0, 0.169156519395D0,
     &     0.182603415045D0,   0.189450610455D0, 0.189450610455D0,
     &     0.182603415045D0,   0.169156519395D0, 0.149595988817D0,
     &     0.124628971256D0,   0.095158511683D0, 0.062253523937D0,
     &     0.027152459411D0/
      END
      subroutine van_der_waals(nfreq,alpha,work,lfree)
C***********************************************************************
C
C     Written by R. Pilemalm and P. Norman 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C   
      real*8 C6
      dimension alpha(nfreq),work(lfree)
C
      kfree=1
      lwrk=lfree
C
      c6=0.0d0
      call gauss_legendre(nfreq,alpha,alpha,c6)
C
       CALL AROUND('Dispersion coefficient')
       WRITE(LUPRI,'(A)')
     &  'Please cite:',
     &  '   Sebastien Villaume, Trond Saue and Patrick Norman,',
     &  '   J. Chem. Phys. 133 (2010) 064105'
       WRITE(LUPRI,'(A)')
       WRITE(LUPRI,'(A)')'@'
       WRITE(LUPRI,'(A3,A10,F15.8,A5)') '@  ',
     &                         'C6       ', C6, 'a.u.'
       WRITE(LUPRI,'(A)')'@'
       WRITE(LUPRI,'(/A)')
     &  '----------------------------------------------'
       
C
      return 
      end
      subroutine gauss_legendre(N,A,B,C6)
C***********************************************************************
C
C     Written by R. Pilemalm and P. Norman 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "gauss_legendre.h"
C
C     Transformation of integration interval is made with w0=0.3 
C     in accordance with
C       R. D. Amos, N. C. Handy, P. J. Knowles, J. E. Rice and A. J. Stone
C       J. Phys. Chem. 89 (1985) 2186-2192
C     There is a slight typo in Eq. 9. The correct equation is Eq. 4 of
C       A. Jiemchooroj, P. Norman and B. E. Sernelius,
C       J. Chem. Phys. 125 (2006) 124306
C
      parameter ( w0 =0.3D0 )
C
      dimension A(N),B(N),abscissas(max_gl_pt),weights(max_gl_pt)
C
      if (N.eq.6) then
         do i=1,n
            abscissas(i)=abscissas_6pt(i)
            weights(i)=weights_6pt(i)
         end do
      else if (N.eq.8) then
         do i=1,n
            abscissas(i)=abscissas_8pt(i)
            weights(i)=weights_8pt(i)
         end do
      else if (N.eq.10) then
         do i=1,n
            abscissas(i)=abscissas_10pt(i)
            weights(i)=weights_10pt(i)
         end do
      else if (N.eq.12) then
         do i=1,n
            abscissas(i)=abscissas_12pt(i)
            weights(i)=weights_12pt(i)
         end do
      else if (N.eq.14) then
         do i=1,n
            abscissas(i)=abscissas_14pt(i)
            weights(i)=weights_14pt(i)
         end do
      else if (N.eq.16) then
         do i=1,n
            abscissas(i)=abscissas_16pt(i)
            weights(i)=weights_16pt(i)
         end do
      else
         write(lupri,*) 'Error in gauss-legendre integration'
         call quit('Error in gauss-legendre integration')
      endif
C     
      DO I=1,N
         C6 = C6 + weights(I)*A(I)*B(I)/(1 - abscissas(I))**2
      END DO
C
      C6 = 6.0D0/pi * w0 * C6
C     
      return
      end
