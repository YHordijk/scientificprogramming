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

C
C Purpose: General Purpose utility routines originally made for krmc,
C          but some are now also used other places (e.g. prp/).
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prdnz3 */
      SUBROUTINE PRDNZ3(AMAT,NDIM1,NDIM2,NZ,IP,LUNIT)
C***********************************************************************
C
C     Print a matrix in the Dirac (NZ,3) format.
C
C     Input:
C       AMAT  - the matrix
C       NDIM1, NDIM2, NZ - dimensions.
C       IP    - index vector from packed quarternion to quaternion.
C       LUNIT - print unit
C
C     Output:
C       None
C
C     Written by J. Thyssen - Nov 24 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NDIM1,NDIM1,NDIM2,NZ,3)
      CHARACTER M(4)*2
      DIMENSION IP(4)
      DATA M /'AR', 'AI', 'BR', 'BI' /
C
      IF (NDIM1 * (NDIM1+1) / 2 .NE. NDIM2) THEN
         WRITE(LUPRI,9000) NDIM1, NDIM2
         CALL QUIT('*** ERROR in PRDNZ3 ***')
      END IF
C
 9000 FORMAT(/,'*** ERROR in PRDNZ3 ***',/,
     &     'Inconsistent dimensions: NDIM1 = ',I5,
     &     ', NDIM2 = ',I5)
C
C     Print matrix.
C
 9001 FORMAT(/3X,'Indices',5X,'Matrix',4X,'Elements',
     &       /1X, 65('-'))
 9010 FORMAT(1P,1X,4I3,4X,A2,6X,3(D14.7,1X))
 9011 FORMAT(1P,17X,A2,6X,3(D14.7,1X))
C
      WRITE(LUNIT,9001)
C
      IJ = 0
      DO I = 1,NDIM1
         DO J = 1,I
            IJ = IJ + 1
            DO K = 1,NDIM1
               DO L = 1,NDIM1
                  WRITE(LUNIT,9010)
     &                 K,L,I,J,
     &                 M(IP(1)),
     &                 (AMAT(K,L,IJ,1,ICL),ICL = 1,3)
                  DO IZ = 2, NZ
                     WRITE(LUNIT,9011) M(IP(IZ)),
     &                    (AMAT(K,L,IJ,IZ,ICL),ICL = 1,3)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine mk34bindex_offset(indxb34_offset,nstr,trian)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
      logical trian
      integer nstr(2,0:2,2)
      integer indxb34_offset(2)
C
      ij = 0
      do irepij =1, nfsym
         indxb34_offset(irepij) = ij
         do irepj = 1, nfsym
            irepi = mod(irepj+irepij,2) + 1
            if(.not.trian.or.(trian.and.irepi.gt.irepj))then
              do j = 1, nstr(irepj,0,2)
                 do i = 1, nstr(irepi,0,1)
                    ij = ij + 1
                 end do
              end do
            else if(irepi.eq.irepj)then
!             same procedure for lower trangular gg or uu case..
              do j = 1, nstr(irepj,0,2)
                 do i = j, nstr(irepi,0,1)
                    ij = ij + 1
                 end do
              end do
            end if
         end do
      end do
!
!     write(lupri,'(2x,a,2i4)') '(mk34bindex_offset) offsets are',
!    &                          indxb34_offset(1),indxb34_offset(2)
!     write(lupri,'(2x,a,2i4)') '(mk34bindex_offset) input nstr(1) is',
!    &                           (nstr(i,0,1), i = 1, nfsym)
!     write(lupri,'(2x,a,2i4)') '(mk34bindex_offset) input nstr(2) is',
!    &                           (nstr(i,0,2), i = 1, nfsym)
!
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine int_irreps_b2x_2(irepij,indxb34_reod,orb_info_array)
C***********************************************************************
#include "implicit.h"
#include "priunit.h"

#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dcbopt.h"
#include "dcborb.h"
#include "dgroup.h"
      dimension indxb34_reod(2,*) 
      integer orb_info_array(2,norbt,norbt)

      iii = 0
      jjj = 0
      ikl = 0
      do irepl = 1, nfsym
        irepk = mod(irepl + irepij, 2) + 1
        do 10 icol = 1, nidx4(irepl)
          jj  = idxt2g(icol+iidx4(irepl),4)
          if(jj .le. 0) goto 10
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!         write(lupri,'(2x,a,i4,a,i4)')
!    &    ' (int_irreps_b2x_2) icol ==',icol, ': index jj =',jj
#endif
          do 20 irow = 1, nidx3(irepk)
            ii = idxt2g(irow+iidx3(irepk),3)
            if(ii .le. 0) goto 20
#if defined MOD_DEBUG || defined MCSCF_DEBUG_SPINFREE
!           write(lupri,'(2x,a,i4,a,i4)')
!    &      ' (int_irreps_b2x_2) irow ==',irow, ': index ii =',ii
!           write(lupri,'(2x,a,i6)')
!    &      '(int_irreps_b2x_2) reorder position - choice 1 ==',
!    &      (jj-1)*norbt+ii
!           write(lupri,'(2x,a,2i6)')
!    &      '(int_irreps_b2x_2) symmetry + quaternion phase factor ==',
!    &      orb_info_array(1,ii,jj), orb_info_array(2,ii,jj)
#endif
            indxb34_reod(1,(jj-1)*norbt+ii) = orb_info_array(1,ii,jj)
            indxb34_reod(2,(jj-1)*norbt+ii) = orb_info_array(2,ii,jj)
 20       continue
 10     continue
      end do

      end
!***********************************************************************
      subroutine setorb_attrib(orb_info_array,ibeig)
!
!     purpose: setup the orbital info array which will on output contain 
!              for the ( |kl) integral either 
!
!              - linear symmetry run : mj-values for each orbital k,l
!
!              - spinfree calculation: total boson symmetry of
!                                      kl-compound and the quaternion 
!                                      multiplication factor  
!
!     written by Stefan Knecht - Aug 2010
!     based on the routine MKINDXB originally written by Luuk Visscher
!
!***********************************************************************
      use memory_allocator
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
      integer orb_info_array(2,norbt,norbt)
      dimension ibeig(*)
      integer, allocatable  :: orb_sym_vec(:)
      character (len=8)     :: orbital_info_label
#include "ibtfun.h"
 
!     step 1: read mj-values/boson irreps from file KRMCSCF
      call alloc(orb_sym_vec,norbt,'orb_sym_vec-setorb_att')
      orb_sym_vec = 0

      call icopy(norbt,ibeig,1,orb_sym_vec,1)
C
C     debug print
!     write(lupri,'(2x,a)') 
!    & '(setorb_attrib) info: symmetry-info vector of orbitals'
!     call iwrtma(orb_sym_vec,1,norbt,1,norbt,lupri)

!     initialize
      call izero(orb_info_array,2*norbt**2)

!     step 2: setup the orbital info array
      do irepij =1, nfsym
         do irepl4 = 1, nfsym
            irepk3 = mod(irepl4+irepij,2) + 1
            do l4 = 1, norb(irepl4)
               index_orb4 = iorb(irepl4) + l4
               do k3 = 1, norb(irepk3)
                  index_orb3 = iorb(irepk3) + k3
                  irep3b = orb_sym_vec(index_orb3)
                  irep4b = orb_sym_vec(index_orb4)
!                 write(lupri,'(2x,a,2i5)') 
!    &            '(setorb_attrib) index_orb3 and index_orb4',
!    &                             index_orb3, index_orb4
                  if(linear)then
 
!                   store the mj values
                    orb_info_array(1,index_orb3,index_orb4) = irep3b
                    orb_info_array(2,index_orb3,index_orb4) = irep4b
                  else
 
!                   compound boson irrep can be obtained by an xor operation
                    irep34b                      = ibtxor(irep3b,irep4b)
!                   write(lupri,'(2x,a,3i5)') 
!    &              '(setorb_attrib) compound irrep for k,l',
!    &              irep34b,irep3b,irep4b
                    orb_info_array(1,index_orb3,
     &                               index_orb4) = irep34b
 
!                    We first determine which quaternion phase
!                    factors the individual functions receive and
!                    than obtain the compound phase factor by
!                    quaternion multiplication.
!
                     iq3 = ipqtoq(1,irep3b)
                     iq4 = ipqtoq(1,irep4b)
 
!                    The sign of the quaternion multiplication is given by iqphase,
!                    but since IQJ is conjugated (backtransform !) we also have to 
!                    multiply by IQSIGN (-1 for q-imaginary, 1 for real)
!                    write(lupri,'(2x,a,3i5)') 
!    &               '(setorb_attrib) phase factor for k,l',
!    &               iqphase(iq3,iq4,1)*iqsign(iq4,2,1),irep3b,irep4b
!
                     orb_info_array(2,index_orb3,index_orb4) =
     &                              iqphase(iq3,iq4,1)*iqsign(iq4,2,1)
                  end if
               end do
            end do
         end do
      end do

!     release scratch memory
      call dealloc(orb_sym_vec)

!     write(lupri,'(2x,a,2i4)') '(MKINDXB) input nstr(1) is',
!    &                           (nstr(i,0,1), i = 1, nfsym)
!     write(lupri,'(2x,a,2i4)') '(MKINDXB) input nstr(2) is',
!    &                           (nstr(i,0,2), i = 1, nfsym)

      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LINZERO */
      SUBROUTINE linzero_nz3(iclass,irepijb,irepklb,nkl,vnz3)
C
C     Zero integrals that should be zero due to symmetry
C     This means that the MJ values stored in IREPIJB and IREPKLB should match.
C
!     written by Stefan Knecht, August 2010
!     based on the routine linzero written by Luuk Visscher.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
      DIMENSION vnz3(nkl,nz,3)
      DIMENSION IREPIJB(2),IREPKLB(2,NKL)
C
C     to determine the MJ value we need to know the classes :
!     Class 1 : (i    j    | k    l   )
!     Class 2 : (i    jbar | k    lbar)
!     Class 3 : (ibar j    | kbar l   )
!
!     Luuks classes...
C     Class 1 : (i    j    | k    l   )
C     Class 2 : (i    j    | kbar lbar)
C     Class 3 : (i    jbar | k    lbar)
C     Class 4 : (ibar j    | k    lbar)
C
      if(nz.gt.2)then
        write(lupri,'(2x,a)') '*** error in (linzero_nz3):'//
     &  'code for nz > 2 is missing yet. ***'
        call quit('*** error in (linzero_nz3): code for nz > 2 is
     &  missing yet. ***')
      end if

        write(lupri,'(2x,a,1i6)') '(linzero_nz3): 
     &         iclass is:', iclass
      IF (ICLASS.eq.1) THEN
         MIJ = - IREPIJB(1) + IREPIJB(2)
      ELSEIF (ICLASS.EQ.2) THEN
         MIJ = - IREPIJB(1) - IREPIJB(2)
      ELSE
         MIJ =   IREPIJB(1) + IREPIJB(2)
      ENDIF
C
      DO I = 1, NKL
         IF (ICLASS.EQ.1) THEN
            MKL = - IREPKLB(1,I) + IREPKLB(2,I)
         ELSE IF(ICLASS.EQ.2)THEN
            MKL = - IREPKLB(1,I) - IREPKLB(2,I)
         ELSE
            MKL =   IREPKLB(1,I) + IREPKLB(2,I)
         ENDIF
         do j = 1, nz
           IF (MIJ+MKL.NE.0)then
           write(lupri,'(2x,a,i6,a,1p,d10.2,4i4)')
     &     '(linzero_nz3): removed integral(',i,')',vnz3(i,j,iclass),
     &     MIJ,MKL, IREPKLB(1,I), IREPKLB(2,I)
           else
           write(lupri,'(2x,a,i6,a,1p,d10.2,4i4)')
     &     '(linzero_nz3): kept integral(',i,')',vnz3(i,j,iclass),
     &     MIJ,MKL, IREPKLB(1,I), IREPKLB(2,I)
           end if
           IF (MIJ+MKL.NE.0) vnz3(i,j,iclass) = 0.D0
         end do
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rprkap */
      SUBROUTINE RPRKAP(XKAPEE,XKAPEP,JXOPE,NZXOPE,JXOPP,NZXOPP,
     &                  NZ,IP,LUPRI)
C***********************************************************************
C
C     Print Kappa vector
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
C
      DIMENSION XKAPEE(NZXOPE,NZ),XKAPEP(NZXOPP,NZ),
     &          JXOPE(2,NZXOPE),JXOPP(2,NZXOPP),IP(4)
      CHARACTER M(4)*1
      DATA M /'1','i','j','k'/
C
      IF (NZXOPE .EQ. 0) GOTO 100
      WRITE(LUPRI,9000) 'e-e'
      WRITE(LUPRI,9005) (M(IP(I)),I=1,NZ)
      WRITE(LUPRI,9006) ('---------------  ',I=1,NZ)
      DO IG = 1,NZXOPE
         WRITE(LUPRI,9010) IG,JXOPE(1,IG),JXOPE(2,IG),
     &                     (XKAPEE(IG,I),I=1,NZ)
      END DO
      DN = DNORM2(NZXOPE*NZ,XKAPEE,1)
      WRITE(LUPRI,9015) DN
  100 CONTINUE
      IF (NZXOPP .EQ. 0) GOTO 200
      WRITE(LUPRI,9000) 'e-p'
      WRITE(LUPRI,9005) (M(I),I=1,NZ)
      WRITE(LUPRI,9006) ('---------------  ',I=1,NZ)
      DO IG = 1,NZXOPP
         WRITE(LUPRI,9010) IG,JXOPP(1,IG),JXOPP(2,IG),
     &                     (XKAPEP(IG,I),I=1,NZ)
      END DO
      DN = DNORM2(NZXOPP*NZ,XKAPEP,1)
      WRITE(LUPRI,9015) DN
 9000 FORMAT(/,1X,A,' part of matrix:')
 9005 FORMAT(1P,/,1X,'Index(r,s)   r    s   ',4(A1,'-part           '))
 9006 FORMAT(1X, '----------  ---  ---  ',4A17)
 9010 FORMAT(1P,5X,I4,4X,I3,2X,I3,2X,4(D15.8,'  '))
 9015 FORMAT(/,1X,'Norm of vector: ',1P,D15.8)
C
  200 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrtkrmc */
      SUBROUTINE WRTKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Write VECTOR to LUKRMC under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8
C
C     Write label at end or - if already existing - after SODLABEL
C
      CALL LABKRMC(LUNIT,LABEL)
C
C     Write vector
C
      IF (ISIZE .GE. 4) THEN
         WRITE(LUNIT) VECTOR
      ELSE ! minimum record length of 32 bytes (for label searches)
         WRITE(LUNIT) VECTOR(1:ISIZE),(DUMMY,I=1,4)
      END IF
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck KRCI_PRPFILE */
      SUBROUTINE KRCI_PRPFILE(LUNIT,LABEL,VECTOR,ISIZE,IRWRT)
C***********************************************************************
C
C     Read/write VECTOR from/to LUNIT under label LABEL.
C
C     IRWRT controls reading/writing:
C                                     0: read vector
C                                     1: write vector
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C      
C
C     Output:
C
C        VECTOR  - vector of size ISIZE (if not input)   
C
C     Based on WRTKRMC.
C
C     Written by S. Knecht - Aug 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8
C
      NTEST = 0
C
      IF( IRWRT .eq. 0 )THEN
C
C       Write label at end or - if already existing - after SODLABEL
C
        CALL LABKRMC(LUNIT,LABEL)
C
C       debug print
        IF( NTEST .gt. 0)THEN
          WRITE(LUPRI,*) ' matrix put to file with label ',LABEL
          CALL WRTMATMN(VECTOR,1,ISIZE,1,ISIZE,LUPRI)
        END IF
C       ... write
        WRITE(LUNIT) VECTOR
        CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
        REWIND LUNIT
      ELSE
        IF( NTEST .gt. 0)THEN
          WRITE(LUPRI,*)' matrix to be read from file with label ',LABEL
        END IF
C       ... read from file
        REWIND LUNIT
        CALL MOLLAB(LABEL,LUNIT,LUPRI)
        CALL READT(LUNIT,ISIZE,VECTOR)
C
C       debug print
        IF( NTEST .gt. 0)THEN
          WRITE(LUPRI,*) ' matrix read from file with label ',LABEL
          CALL WRTMATMN(VECTOR,1,ISIZE,1,ISIZE,LUPRI)
        END IF
      END IF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck GET_KRCICVECS_E */
      SUBROUTINE GET_KRCICVECS_E(EVEC,NKRCI_CIROOTS,XSYMFLAB,
     &                           NKRCI_MAX_SYM,PRINT_EV)
C***********************************************************************
C
C     read eigenstate energies from file KRCI_CVECS.INFO
C
C     Input:
C        NKRCI_MAX_SYM  - highest active symmetry irrep
C        XSYMFLAB       - array containing file labels for each symmetry
C        NKRCI_CIROOTS  - # eigenstates per symmetry irrep
C
C     Output:
C        EVEC     - vector containing eigenvalues for all eigenstates
C        PRINT_EV - flag for printing eigenvalues in calling routine
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision: S. Knecht - Jan 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION EVEC(*), NKRCI_CIROOTS(*)
      CHARACTER*3 XSYMFLAB(*) 
      CHARACTER*8 FILELAB
      LOGICAL EX, PRINT_EV, FNDLAB
C
C     read final energies from file KRCI_CVECS.INFO
      INQUIRE(FILE = 'KRCI_CVECS.INFO',EXIST=EX)
      IF( .not. EX ) THEN
        WRITE(LUPRI,'(/A,2I4)')' *** ERROR in GET_KRCICVECS_E.'//
     &                         ' file KRCI_CVECS.INFO not found.'
        CALL QUIT(' *** ERROR in GET_KRCICVECS_E. file KRCI_CVECS.INFO
     &               not found.')
      END IF
      LU_INFO = 80
      CALL OPNFIL(LU_INFO,'KRCI_CVECS.INFO','OLD','XGETEN')
      NXSCR = 1
      DO I = 1, NKRCI_MAX_SYM
        IF(NKRCI_CIROOTS(I) .gt. 0 )THEN
          WRITE(FILELAB,'(A5,A3)') "eroot",XSYMFLAB(I)
          REWIND LU_INFO
          PRINT_EV = FNDLAB(FILELAB,LU_INFO)
          IF(PRINT_EV) THEN 
            CALL KRCI_PRPFILE(LU_INFO,FILELAB,EVEC(NXSCR),
     &                        NKRCI_CIROOTS(I),1)
          ELSE
            WRITE(LUPRI,'(/A,/A)')
     &      ' *** Warning: energies for requested eigenstates are not'//
     &      ' available. Applying a simple order per symmetry. ***',
     &      ' *** DO NOT RELY ON ANY ORDERING!'//
     &      ' DEGENERACIES MAY NOT BE DETECTED... ***'
            DO IVEC = 1,NKRCI_CIROOTS(I)
              EVEC(IVEC+NXSCR-1) = (IVEC+NXSCR-1) * 1.0D0
            END DO
          END IF
        END IF
        NXSCR = NXSCR + NKRCI_CIROOTS(I)
      END DO
C     ... close
      CLOSE(LU_INFO,STATUS='KEEP')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck iwrtkrmc */
      SUBROUTINE IWRTKRMC(LUNIT,LABEL,IVECTOR,ISIZE)
C***********************************************************************
C
C     Write VECTOR to LUKRMC under label LABEL.
C
C     Input:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Output:
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "iratdef.h"
C
      DIMENSION IVECTOR(ISIZE)
      CHARACTER LABEL*8
C
C     Write label at end or - if already existing - after SODLABEL
C
      CALL LABKRMC(LUNIT,LABEL)
C
C     Write vector
C
      LEN = MAX(4*IRAT,ISIZE)
      CALL WRITI(LUNIT,LEN,IVECTOR)
      CALL NEWLAB('EOFLABEL',LUNIT,LUPRI)
C
      REWIND LUNIT
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck labkrmc */
      SUBROUTINE LABKRMC(LUNIT,LABEL)
C***********************************************************************
C
C     Write label LABEL to LUNIT;
C       If label LABEL exists:
C         Seek to label SODLABEL and write from there.
C       Else:
C         Seek to label EOFLABEL, backspace and write from there
C       End If
C
C     Input:
C        LABEL   - the label
C
C     Output:
C
C     Written by J. Thyssen - Dec 28 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
      CHARACTER LABEL*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C     Label exists.
C     Seek to label SODLABEL and write from there.
C
         REWIND(LUNIT)
         IF (.NOT. FNDLAB('SODLABEL',LUNIT)) THEN
            CALL QUIT('LABKRMC: wrong structure in file - no SODLABEL')
         END IF
      ELSE
C
C     Label does not exist
C     Seek to label EOFLABEL and write from there
C
         REWIND(LUNIT)
         IF (FNDLAB('EOFLABEL',LUNIT))
     &      BACKSPACE LUNIT
C
      END IF
C
C     Write label
C
      CALL NEWLAB(LABEL,LUNIT,LUPRI)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reakrmc */
      SUBROUTINE REAKRMC(LUNIT,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUKRMC under label LABEL.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         CALL READT(LUNIT,ISIZE,VECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     &        '*** ERROR in REAKRMC ***',
     &        'Label <',LABEL,'> does not exist on unit ',LUNIT
         CALL QUIT('*** ERROR in REAKRMC ***')
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ireakrmc */
      SUBROUTINE IREAKRMC(LUNIT,LABEL,IVECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUKRMC under label LABEL.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IVECTOR(ISIZE)
      CHARACTER LABEL*8
      LOGICAL   FNDLAB
C
      REWIND(LUNIT)
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         CALL READI(LUNIT,ISIZE,IVECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     $        '*** ERROR in IREAKRMC ***',
     $        'Label <',LABEL,'> does not exists on unit ',LUNIT
         CALL QUIT('*** ERROR in IREAKRMC ***')
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reakrmc */
      SUBROUTINE REAKRMC2(LUNIT,FILENAME,LABEL,VECTOR,ISIZE)
C***********************************************************************
C
C     Read VECTOR from LUNIT under label LABEL.
C
C     Open & close files if necessary.
C
C     Input:
C
C     Output:
C        VECTOR  - vector of size ISIZE
C        LABEL   - the label
C        LUNIT   - unit number
C
C     Written by J. Thyssen - Jul 4 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
      DIMENSION VECTOR(ISIZE)
      CHARACTER LABEL*8
      CHARACTER FILENAME*(*)
      LOGICAL   FNDLAB, FEXIST, FOPEN
C
C
C     Check if file exists
C
      INQUIRE(FILE=FILENAME,EXIST=FEXIST)
      IF (FEXIST) THEN
         INQUIRE(FILE=FILENAME,OPENED=FOPEN)
         IF (.NOT. FOPEN) THEN
            OPEN(LUNIT,FILE=FILENAME,STATUS='UNKNOWN',
     &           FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         END IF
      ELSE
         FOPEN = .FALSE.
         OPEN(LUNIT,FILE=FILENAME,STATUS='NEW',
     &        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         REWIND(LUNIT)
         CALL NEWLAB('SODLABEL',LUNIT,LUPRI)
      END IF
C
      IF (FNDLAB(LABEL,LUNIT)) THEN
C
C        Label exists.
C
         LEN = MAX(4,ISIZE)
         CALL READT(LUNIT,ISIZE,VECTOR)
C
      ELSE
         WRITE(LUPRI,'(//A/,3A,I3)')
     $        '*** ERROR in REAKRMC ***',
     $        'Label <',LABEL,'> does not exists on unit ',LUNIT
         CALL QUIT('*** ERROR in REAKRMC ***')
      END IF
C
      IF (.NOT. FOPEN) CLOSE(LUNIT)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rupkwop */
      SUBROUTINE RUPKWOP(WOPEE,WOPEP,UNWOP,NORBT,
     &                   JXOPE,NZXOPE,JXOPP,NZXOPP,NZ)
C***********************************************************************
C
C     Unpack WOP into UNWOP
C     WOP is a rotation type vector and UNWOP is a square matrix.
C
C     Input:
C        WOPEE   - e-e orbital rotation type vector
C        WOPEP   - e-p orbital rotation type vector
C        JXOPE   - e-e rotations
C        JXOPP   - e-p rotations
C
C     Output:
C        UNWOP   - unpacked matrix
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION WOPEE(NZXOPE,NZ),WOPEP(NZXOPP,NZ),
     &          UNWOP(NORBT,NORBT,NZ),JXOPE(2,*),JXOPP(2,*)
C
C     e-e rotations
C
      DO IG = 1,NZXOPE
         K = JXOPE(1,IG)
         L = JXOPE(2,IG)
         UNWOP(K,L,1) =  WOPEE(IG,1)
         UNWOP(L,K,1) = -WOPEE(IG,1)
         DO IZ = 2,NZ
            UNWOP(K,L,IZ) =  WOPEE(IG,IZ)
            UNWOP(L,K,IZ) =  WOPEE(IG,IZ)
         END DO
      END DO
C
C     e-p rotations
C
      DO IG = 1,NZXOPP
         K = JXOPP(1,IG)
         L = JXOPP(2,IG)
         UNWOP(K,L,1) =  WOPEP(IG,1)
         UNWOP(L,K,1) = -WOPEP(IG,1)
         DO IZ = 2,NZ
            UNWOP(K,L,IZ) =  WOPEP(IG,IZ)
            UNWOP(L,K,IZ) =  WOPEP(IG,IZ)
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rotcmo */
      SUBROUTINE ROTCMO(CMO,NBASI,NORBI,NZ,XKAPPA,NORBT,IP,IPRINT,
     &                  WORK,LWORK)
C***********************************************************************
C
C     Rotate coefficients
C
C     CMO_{new} = CMO_{old} exp( -\kappa )
C     We use a trick a la Sostrup notes, chapter 3.1.5
C
C     Input:
C        CMO     - old coefficients
C        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
C        IPRINT  - print level
C
C     Output:
C        CMO     - new coefficients
C
C     Written by J. Thyssen - Nov 24 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
C
      PARAMETER (D1 = 1.00 D00, D0 = 0.00 D00, DM1 = -D1)
      PARAMETER (DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00)
      PARAMETER (D1P5 = 1.5D00, D2 = 2.0D00, THREQL = 1.0D-12)
C
      DIMENSION CMO(NBASI,NORBI,NZ), XKAPPA(NORBT,NORBT,NZ),
     &          IP(4), WORK(*)
      CHARACTER SECTID*12,CPUTID*12,WALLTID*12
C
#include "memint.h"
C
      CALL QENTER('ROTCMO')
      CALL GETTIM(CPU1,WALL1)
      IF (IPRINT .GE. 20) CALL HEADER('Output from ROTCMO',-1)
C
#ifndef OLDCODE
C
      NORBXQ = NORBI*NORBI*NZ
      CALL MEMGET('REAL',KEIG,NORBI,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KKSQ,NORBI*NBASI*NZ,WORK,KFREE,LFREE)
      KCOS = KKSQ
      CALL MEMGET('REAL',KSIN,NORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KV  ,NORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP,NORBXQ,WORK,KFREE,LFREE)
C
C     Calculate -\kappa^2 = \kappa^{\dagger} \kappa
C     ---------------------------------------------
C
      CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &           'H','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           'N','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           D0,IP,WORK(KKSQ),NORBI,NORBI,NZ)
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: -kappa^2',-1)
         CALL PRQMAT(WORK(KKSQ),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Diagonalize -\kappa^2
C     ---------------------------------------------
C
      CALL QDIAG(NZ,NORBI,WORK(KKSQ),NORBI,NORBI,
     &           WORK(KEIG),1,WORK(KV),NORBI,NORBI,
     &           WORK(KFREE),LFREE,IERR)
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: V matrix (eigenvector)',-1)
         CALL PRQMAT(WORK(KV),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO: eigenvalues',-1)
         CALL PRQMAT(WORK(keig),1,NORBI,1,NORBI,1,IP,LUPRI)
      END IF
C
C     Make cos(d) and sin(d)/d ; d = sqrt(eigenvalues)
C     ------------------------------------------------
C
      CALL GETTIM(CPU2,WALL2)
      CALL DZERO(WORK(KCOS),NORBXQ)
      CALL DZERO(WORK(KSIN),NORBXQ)
      DO I = 1,NORBI
         IADDR = (I-1) + (I-1)*NORBI
         EIG = WORK(KEIG+I-1)
         IF (EIG .LT. D0) THEN
            IF (IPRINT .GE. 3 .AND. ABS(EIG).GT. THRZER)
     &         WRITE(LUPRI,9010) EIG
            EIG = D0
            WORK(KEIG+I-1)   = D0
         END IF
         IF (EIG .LT. DTHRS) THEN
C
C           Use Taylor expansion to avoid division by zero and/or
C           to speed up evaluation (avoid call to sin/cos/sqrt).
C           For EIG < 1.0D-8 the 3. order Taylor expansion is correct
C             to machine precision.
C
            WORK(KCOS+IADDR) = D1-DP5*EIG
            WORK(KSIN+IADDR) = D1-DSIXTH*EIG
         ELSE
            SQEIG = SQRT(EIG)
            WORK(KCOS+IADDR) = COS(SQEIG)
            WORK(KSIN+IADDR) = SIN(SQEIG)/SQEIG
         END IF
      END DO
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: sin(d)/d',-1)
         CALL PRQMAT(WORK(KSIN),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO: cos(d)',-1)
         CALL PRQMAT(WORK(KCOS),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU3,WALL3)
 9010 FORMAT('WARNING from ROTCMO: Negative eigenvalue ',1P,D15.6)
C
C     Make V cos(d) V^{\dagger}
C     (use KTMP as temporatory storage)
C     ------------------------------------
C
      CALL QTRANS('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI,
     &            WORK(KTMP),NORBI,NORBI,NZ,IP,
     &            WORK(KCOS),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KFREE),LFREE,IPRINT)
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: V cos(d) V^{dagger}',-1)
         CALL PRQMAT(WORK(KTMP),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Make V sin(d)/d V^{\dagger}
C     (use KCOS as temporatory storage)
C     ------------------------------------
C
      CALL QTRANS('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI,
     &            WORK(KCOS),NORBI,NORBI,NZ,IP,
     &            WORK(KSIN),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KV),NORBI,NORBI,NZ,IP,
     &            WORK(KFREE),LFREE,IPRINT)
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(WORK(KCOS),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
C     ---------------------------------------------------------
C
      CALL QGEMM(NORBI,NORBI,NORBI,DM1,
     &           'N','N',IP,WORK(KCOS),NORBI,NORBI,NZ,
     &           'N','N',IP,XKAPPA,NORBT,NORBT,NZ,
     &           D1,IP,WORK(KTMP),NORBI,NORBI,NZ)
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: exp(-kappa)',-1)
         CALL PRQMAT(WORK(KTMP),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C     Check if unitary
C
c#ifdef KRMC_DEBUG
      CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &           'N','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           'H','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KCOS),NORBI,NORBI,NZ)
      ioff = 0
      i_max = 0
      j_max = 0
      d_max = D0
      do i = 1,norbi
      do j = 1,norbi
         if (i.eq.j) then
            if (abs(work(kcos+ioff)-d1).gt.d_max) then
               i_max = i
               j_max = j
               d_max = abs(work(kcos+ioff)-d1)
            end if
         else
            if (work(kcos+ioff).gt.d_max) then
               i_max = i
               j_max = j
               d_max = abs(work(kcos+ioff))
            end if
         end if
         ioff = ioff + 1
      end do
      end do
      if (d_max .gt. thrzer) then
         write(LUPRI,*) 'ROTCMO: maximum deviation from unitary matrix'
         write(LUPRI,*) 'Element ',j_max,i_max,' value :',d_max
      end if
C#endif
C
C
C     Calculate new orbitals: CMO_{old} exp(-kappa)
C     ---------------------------------------------
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
CSK   CALL DZERO(WORK(KCOS),NORBXQ)
      CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &           'N','N',IP,CMO,NBASI,NORBI,NZ,
     &           'N','N',IP,WORK(KTMP),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KCOS),NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,WORK(KCOS),1,CMO,1)
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('ROTCMO: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (IPRINT .GE. 20) THEN
         CPUTID = SECTID(CPU4-CPU1)
         WALLTID = SECTID(WALL4-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID
      END IF
 9020 FORMAT(/'ROTCMO: total CPU (WALL) time       : ',A12,' (',A12,')')
 9021 FORMAT( 'ROTCMO: sin/sqrt/cos CPU (WALL) time: ',A12,' (',A12,')')
C
#else
      call quit('this code does not work')
C     fixme: dette ser ikke ud til at virke for complex groups.
C
C     Memory alloc
C
      CALL MEMGET('REAL',KXKOLD,NORBI*NORBI*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXKNEW,NORBI*NORBI*NZ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP  ,NBASI*NBASI*NZ,WORK,KFREE,LFREE)
C
C
C     Calculate exp(-K) = 1 - K
C
      K = 0
      DO IZ = 1,NZ
         DO J = 1,NORBI
            DO I = 1,NORBI
               K = K + 1
               WORK(KXKNEW + K - 1) = -XKAPPA(I,J,IZ)
               IF ( I .EQ. J .AND. IZ .EQ. 1)
     $              WORK(KXKNEW + K - 1 ) = WORK(KXKNEW + K - 1 ) +  D1
            END DO
         END DO
      END DO
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('ROTCMO: 1 - Kappa',-1)
         CALL PRQMAT(WORK(KXKNEW),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      DMX = D0
      ITS = 0
C
C     Symmetric Orthogonalization with N-R iterative algorithm:
C
C     C(p+1) = 1/2( 3C(p) - C(p) S C^{\dagger}(p) C(p))
C
C            = 3/2 C(p) + ( -1/2 C(p) C^{\dagger}(p) ) C(p)
C
C     (here S = unit matrix)
C
 100  CONTINUE
         ITS = ITS + 1
         CALL DCOPY(NORBI*NORBI*NZ,WORK(KXKNEW),1,WORK(KXKOLD),1)
         CALL QGEMM(NORBI,NORBI,NORBI,D1,
     &        'N','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        'H','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        D0,     IP,WORK(KTMP),  NORBI,NORBI,NZ)
         CALL QGEMM(NORBI,NORBI,NORBI,-DP5,
     &        'N','N',IP,WORK(KTMP),  NORBI,NORBI,NZ,
     &        'N','N',IP,WORK(KXKOLD),NORBI,NORBI,NZ,
     &        D1P5,   IP,WORK(KXKNEW),NORBI,NORBI,NZ)
C
         DMX1 = DMX
         DMX = DABSMAXDIFF(WORK(KXKNEW),WORK(KXKOLD),NORBI*NORBI*NZ)
         IF (DMX .GT. D2 * DMX1 .AND. DMX1 .GT. D0) THEN
            WRITE(LUPRI,'(/1X,2A)')
     &           '*** ERROR in ROTCMO ***: ',
     &           'Divergence in symmetric ortho.'
            CALL QUIT('*** ERROR in ROTCMO: Divergence ***')
         END IF
         IF (DMX .GT. THREQL) GOTO 100
C
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: symm. ortho. rot. mat.',-1)
         CALL PRQMAT(WORK(KXKNEW),NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
C
C
C     Calculate new orbitals: CMO_{old} * (symm.ortho.rot.mat)
C     --------------------------------------------------------
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ROTCMO: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL QGEMM(NBASI,NORBI,NORBI,D1,
     &           'N','N',IP,CMO,NBASI,NORBI,NZ,
     &           'N','N',IP,WORK(KXKNEW),NORBI,NORBI,NZ,
     &           D0,IP,WORK(KTMP),NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,WORK(KTMP),1,CMO,1)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('ROTCMO: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
#endif
      CALL MEMREL('ROTCMO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('ROTCMO')
      RETURN
      END
      DOUBLE PRECISION FUNCTION DABSMAXDIFF(A,B,N)
#include "implicit.h"
      PARAMETER ( D0 = 0.0D00 )
      DOUBLE PRECISION A(N),B(N)
      DMX = D0
      DO I = 1,N
         DEL = ABS(A(I) - B(I))
         DMX = MAX(DMX,DEL)
      END DO
      DABSMAXDIFF = DMX
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setdbcidx */
      SUBROUTINE SETDCBIDX()
C***********************************************************************
C
C     Set common block DCBIDX (index vectors, orbital symmetries etc.)
C
C     See include/dcbidx.h for the definition of the index vectors.
C
C     Written by J. Thyssen - Nov 22 2000
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxash.h"
C
#include "maxorb.h"
#include "dcbidx.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbopt.h"
C
C
      IORBI = 0
C
      IESHI = 0
      IPSHI = 0
C
      IISHI = 0
      IASHI = 0
      ISSHI = 0
C
      DO IFSYM = 1, NFSYM
         DO I = 1, NPSH(IFSYM)
            IORBI = IORBI + 1
            IPSHI = IPSHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTPOSI
C
            IDXG2E(IORBI) = -1
            IDXG2P(IORBI) = IPSHI
            IDXG2I(IORBI) = -1
            IDXG2U(IORBI) = -1
            IDXG2S(IORBI) = -1
C
            IDXP2G(IPSHI) = IORBI
         END DO
         DO I = 1, NISH(IFSYM)
            IORBI = IORBI + 1
            IESHI = IESHI + 1
            IISHI = IISHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTINAC
C
            IDXG2E(IORBI) = IESHI
            IDXG2P(IORBI) = -1
            IDXG2I(IORBI) = IISHI
            IDXG2U(IORBI) = -1
            IDXG2S(IORBI) = -1
C
            IDXE2G(IESHI) = IORBI
            IDXI2G(IISHI) = IORBI
C
         END DO
         DO I = 1, NASH(IFSYM)
            IORBI = IORBI + 1
            IESHI = IESHI + 1
            IASHI = IASHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTACT
C
            IDXG2E(IORBI) = IESHI
            IDXG2P(IORBI) = -1
            IDXG2I(IORBI) = -1
            IDXG2U(IORBI) = IASHI
            IDXG2S(IORBI) = -1
C
            IDXU2G(IASHI) = IORBI
            IDXE2G(IESHI) = IORBI
C
         END DO
         DO I = 1, NSSH(IFSYM)
            IORBI = IORBI + 1
            IESHI = IESHI + 1
            ISSHI = ISSHI + 1
C
            IFSMO(IORBI) = IFSYM
            IOBTYP(IORBI) = JTSEC
C
            IDXG2E(IORBI) = IESHI
            IDXG2P(IORBI) = -1
            IDXG2I(IORBI) = -1
            IDXG2U(IORBI) = -1
            IDXG2S(IORBI) = ISSHI
C
            IDXE2G(IESHI) = IORBI
            IDXS2G(ISSHI) = IORBI
         END DO
      END DO
C
C
      IF (IPROPT .GE. 10) THEN
         WRITE(LUPRI,'(A)')
     &        ' (SETDCBIDX) I, IFSMO, IOBTYP, IDXG2U, IDXG2E'
         DO I = 1, NORBT
            WRITE(LUPRI,'(8X,I6,I7,I8,2I12)')
     &           I, IFSMO(I), IOBTYP(I), IDXG2U(I), IDXG2E(I)
         END DO
         WRITE(LUPRI,'(A)')
     &        ' (SETDCBIDX) E and IDXE2G'
         DO I = 1, NESHT
            WRITE(LUPRI,'(8X,I6,I11)') I, IDXE2G(I)
         END DO
         WRITE(LUPRI,'(A)')
     &        ' (SETDCBIDX) U and IDXU2G'
         DO I = 1, NASHT
            WRITE(LUPRI,'(8X,I6,I11)') I, IDXU2G(I)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck rprvec */
      SUBROUTINE RPRVEC(NDIM,VEC,INCVEC,NZ,IPQ,THRESH,MAXLIN,LUOUT)
C
C 19-Aug-1989 hjaaj
C   print VEC(1:NDIM*INCVEC:INCVEC)
C
C   NDIM      : Number of elements in vector VEC
C   VEC(:)    : Vector to be printed
C   INCVEC    : Increment between each element in vector VEC
C   NZ        : real, complex, or quaternion
C   IPQ       : index packed quaternion to quaternion
C   THRESH    : Print threshold for vector with unit norm
C               (if THRESH .lt. 0 then -THRESH is used without
C                renormalization).
C   MAXLIN    : max. lines of output with vector elements
C   LUOUT     : output unit
C
C
      use quaternion_algebra
#include "implicit.h"
      DIMENSION VEC(NDIM*INCVEC,NZ)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      PARAMETER (D1LOW = D1 - 1.0D-10, D1HGH = D1 + 1.0D-10)
      LOGICAL   VSCALE
      DIMENSION IPQ(4)
C
C     Test input
C
      NERR = 0
      IF (NDIM   .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &     'No print from PRVEC because NDIM .le. 0'
         NERR = NERR + 1
      END IF
      IF (INCVEC .LE. 0) THEN
         WRITE (LUOUT,'(/5X,A)')
     &      'No print from PRVEC because INCVEC .le. 0'
         NERR = NERR + 1
      END IF
      IF (NERR .GT. 0) RETURN
C
C
C
      C2NRM = DNORM2(NDIM*NZ,VEC,INCVEC)
      IF (THRESH .LE. D0) THEN
         THRES2 = -THRESH
         VSCALE = .FALSE.
      ELSE
         THRES2 =  THRESH * C2NRM
         VSCALE = (C2NRM .LT. D1LOW .OR. C2NRM .GT. D1HGH)
      END IF
      IF (VSCALE) THEN
         WRITE (LUOUT,'(/2A,1P,D12.4,/)')
     &      ' Print of vector elements (vector scaled to unit norm)',
     &      ' larger than',ABS(THRESH)
         SCALE = D1 / C2NRM
      ELSE
         WRITE (LUOUT,'(/A,1P,D12.4,/)')
     &      ' Print of vector elements larger than',ABS(THRESH)
         SCALE = D1
      END IF
C
      IPR   = 0
      IZER  = 0
      C2SUM = D0
C
      IF (VSCALE) THEN
         WRITE(LUOUT,9000)
      ELSE
         WRITE(LUOUT,9001)
      END IF
 9000 FORMAT(5X,'Element',3X,'Coefficient (scaled to unit norm)',
     &     /,1X,78('-'))
 9001 FORMAT(5X,'Element',3X,'Coefficient',/,1X,78('-'))
C
      DO 300 I = 1, NDIM
         NA = 1 + (I-1)*INCVEC
         DNNA = D0
         DO IZ = 1, NZ
            DNNA = DNNA + VEC(NA,IZ)**2
         END DO
         IF (SQRT(DNNA).LE.THRES2 .OR. IPR .GE. MAXLIN) THEN
            C2SUM = C2SUM + DNNA
            IF (DNNA .EQ. D0) IZER = IZER + 1
         ELSE
            IF (MOD(IPR,5) .EQ. 0 .AND. IPR .NE. 0) WRITE (LUOUT,'()')
            IPR = IPR + 1
            IF (VSCALE) THEN
               IF (NZ .EQ. 1) THEN
                  WRITE(LUOUT,50)
     &                 I,VEC(NA,1),
     &                 SCALE*VEC(NA,1)
               ELSE IF (NZ .EQ. 2) THEN
                  WRITE(LUOUT,51)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ),
     &                 SCALE*VEC(NA,1),
     &                 (SCALE*VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               ELSE
                  WRITE(LUOUT,52)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ),
     &                 SCALE*VEC(NA,1),
     &                 (SCALE*VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               END IF
            ELSE
               IF (NZ .EQ. 1) THEN
                  WRITE(LUOUT,60)I,VEC(NA,1)
               ELSE IF (NZ .EQ. 2) THEN
                  WRITE(LUOUT,61)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               ELSE
                  WRITE(LUOUT,62)
     &                 I,VEC(NA,1),
     &                 (VEC(NA,IZ),QUNIT(IPQ(IZ)),IZ=2,NZ)
               END IF
            END IF
 50         FORMAT(1X,I10,3X,1P,D10.2,0P,'(',F10.6,')')
 51         FORMAT(1X,I10,3X,1P,D10.2,2X,D10.2,1X,A1,2X,
     &           0P,'(',F10.6,2X,F10.6,1X,A1,')')
 52         FORMAT(1X,I10,3X,1P,D10.2,2X,3(D10.2,1X,A1,2X),
     &           0P,'(',F10.6,2X,3(F10.6,1X,A1,2X),')')
 60         FORMAT(1X,I10,3X,F20.10)
 61         FORMAT(1X,I10,3X,F20.10,2X,F20.10,1X,A1,2X)
 62         FORMAT(1X,I10,3X,F20.10,2X,3(F20.10,1X,A1,2X))
         END IF
  300 CONTINUE
      IF (IPR .GE. MAXLIN) THEN
C
C     *** We have reached the print limit
C
         WRITE (LUOUT,910) IPR
      END IF
  910 FORMAT(/' Print limit of',I6,' elements has been reached.')
      C2SUM = SQRT(C2SUM)
      IF (IZER .EQ. NDIM) THEN
         WRITE (LUOUT,920) NDIM
      ELSE
         WRITE (LUOUT,930) NDIM,NDIM-IPR,IZER,C2SUM,C2NRM
      END IF
  920 FORMAT(/' Length of vector                      :',I10,
     *       /' All elements are zero.',/)
  930 FORMAT(/' Length of vector                      :',I10,
     *       /' Number of elements not printed        :',I10,
     *       /' Number of zero elements               :',I10,
     *       /' Total norm of coefficients not printed:',F10.6,
     *       /' (the coefficients are normalized to    ',F10.6,')',/)
      RETURN
C
C End of PRVEC
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck qnorm */
      SUBROUTINE QNORM(S,VC,NBASI,NORBI,NZ,IPQTOQ,IRETUR,SCR1,SCR2)
C***********************************************************************
C
C     Gram-Schmidt orthogonalize coefficients.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Jan 23 2001
C
C     Last revision: HJAaJ - Aug 2008 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "thrzer.h"
#include "consts.h"
      PARAMETER (THRRND = 0.9D0, THNORM = 1.0D-6)
C
      DIMENSION S(NBASI,NBASI), VC(NBASI,NORBI,NZ)
      DIMENSION IPQTOQ(*), SCR1(NBASI,NZ), SCR2(NORBI)
C
      DO 20 I = 1, NORBI
         ITURN = 0
 1       ITURN = ITURN + 1
         IROUND = 0
C
C        Calculate S * VC(I)
C
         CALL QGEMM(NBASI,1,NBASI,D1,
     &        'N','N',IPQTOQ,S,NBASI,NBASI,1,
     &        'N','N',IPQTOQ,VC(1,I,1),NBASI,NORBI,NZ,
     &             D0,IPQTOQ,SCR1,NBASI,1,NZ)
C
         TNORM = D0
         DO IZ = 1, NZ
            TNORM = TNORM + DDOT(NBASI,VC(1,I,IZ),1,SCR1(1,IZ),1)
         END DO
C        write(lupri,*) 'tnorm = ',tnorm
         IF (TNORM .LT. THRZER) THEN
            IRETUR = -I
            RETURN
         END IF
C
C        Normalize vector
C
         TNORM = D1 / SQRT(TNORM)
         DO IZ = 1, NZ
            CALL DSCAL(NBASI,TNORM,VC(1,I,IZ),1)
            CALL DSCAL(NBASI,TNORM,SCR1(1,IZ),1)
         END DO
C
         TNORM = D0
         DO IZ = 1, NZ
            TNORM = TNORM + DDOT(NBASI,VC(1,I,IZ),1,SCR1(1,IZ),1)
         END DO
C
C        Begin normalization loop
C
         DO J = 1, I - 1
            T = D0
            DO IZ = 1, NZ
               T = T + DDOT(NBASI,VC(1,J,IZ),1,SCR1(1,IZ),1)
            END DO
            TNORM = TNORM - T * T
            SCR2(J) = -T
         END DO
C
         IF (TNORM .LT. THRZER) THEN
C           ... zero vector after orthogonalization
            IRETUR = -I * 1 000 000
            RETURN
         END IF
C
         IF (TNORM .LT. THNORM) IROUND = IROUND + 1
C
         TNORM = D1 / SQRT (TNORM)
         DO J = 1, I - 1
            SCR2(J) = SCR2(J) * TNORM
         END DO
         SCR2(I) = TNORM
C
C        Replace VC(*,I)
C
C        SCR1 = \sum_{k=1,i} SCR2(k) * VC(*,k)
C
         CALL QGEMM(NBASI,1,I,D1,
     &        'N','N',IPQTOQ,VC,NBASI,NORBI,NZ,
     &        'N','N',IPQTOQ,SCR2,NORBI,1,1,
     &             D0,IPQTOQ,SCR1,NBASI,1,NZ)
C
         DO IZ = 1, NZ
            CALL DCOPY(NBASI,SCR1(1,IZ),1,VC(1,I,IZ),1)
         END DO
C
         IF (ITURN .EQ. 1 .AND. IROUND .GT. 0) THEN
C           write(LUPRI,*) 'INFO: second round in QNORM'
            GOTO 1
         END IF
         IF (IROUND .GT. 0) CALL QUIT('QNORM internal error')
 20   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck matgat */
      SUBROUTINE MATGAT(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,NZ)
C***********************************************************************
C
C     Construct BMAT from AMAT elements using IDX.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ), BMAT(NRB,NCB,NZ), IDX(*)
C
      DO IZ = 1, NZ
         DO J = 1, NCB
            JJ = IDX(J)
            DO K = 1, NRB
               KK = IDX(K)
               BMAT(K,J,IZ) = AMAT(KK,JJ,IZ)
            END DO
         END DO
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATGAT_symm(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,
     &                       IBEIG,is_symm,NZ)
C***********************************************************************
C
C     Construct BMAT from AMAT elements using IDX and ibeig arrays.
C
C     Input:
C
C     Output:
C
C     Written by S. Knecht - August 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ), BMAT(NRB,NCB,NZ), IDX(*), IBEIG(*)
      real(8), parameter :: threshold = 1.0d-12
C
!     print *, 'MATGAT_symm start is_symm ==> ',is_symm
      DO IZ = 1, NZ
         do 5 J = 1, NCB
            JJ = IDX(J)
!           print *, ' JJ jz ==> ',JJ, ibeig(jj)
            if(IBEIG(JJ).ne.is_symm) goto 5
            do 10 K = 1, NRB
               KK = IDX(K)
!              print *, ' KK jz ==> ',JJ, ibeig(KK)
               if(IBEIG(kk).ne.is_symm) goto 10
               if(abs(AMAT(KK,JJ,IZ)).lt.threshold)then
                 BMAT(K,J,IZ) = 0.0d0
               else
                 BMAT(K,J,IZ) = AMAT(KK,JJ,IZ)
               end if
 10         continue
  5      continue
      END DO
!     print *, 'MATGAT_symm end  is_symm ==> ',is_symm
      RETURN
      END
!***********************************************************************

      subroutine create_jzbos_sym_vec(jzbos_vec,nvjz,jzbos_tot,nvjz_tot,
     &                                idx,joffset)
!***********************************************************************
!
!     Construct a jz/boson irrep sub-block vector out of the total 
!     jz/boson irrep vector (jzbos_tot) using the IDX array and starting 
!     at a given offset.
!
!     written by S. Knecht - june 2011
!
!***********************************************************************
      implicit none

      integer, intent(in)  :: nvjz, nvjz_tot
      integer, intent(out) :: jzbos_vec(nvjz) 
      integer, intent(in)  :: jzbos_tot(nvjz_tot)
      integer, intent(in)  :: idx(*)
      integer, intent(in)  :: joffset

      integer              :: j, jj
!----------------------------------------------------------------------

      call izero(jzbos_vec,nvjz)
#ifdef MCSCF_debug
      print *, 'create_jzbos_sym_vec start   ==> ',joffset
      print *, 'create_jzbos_sym_vec total # ==> ',nvjz_tot
      print *, 'create_jzbos_sym_vec sub   # ==> ',nvjz
#endif

      do j = 1, nvjz
        jj           = idx(j+joffset)
        jzbos_vec(j) = jzbos_tot(jj)
#ifdef MCSCF_debug
        print *, ' j JJ jz ==> ',j,jj, jzbos_tot(jj)
#endif
      end do

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE MATGATX(AMAT,NDIMA,BMAT,NDIMB,ISKIP,NZ)
C***********************************************************************
C
C     Construct BMAT from AMAT elements.
C
C     Input:
C
C     Output:
C
C     Written by S. Knecht - Nov 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NDIMA,NDIMA,NZ), BMAT(NDIMB,NDIMB,NZ)
C
      DO IZ = 1, NZ
        DO J = 1, NDIMB
          CALL DCOPY(NDIMB,AMAT(ISKIP+1,ISKIP+J,IZ),1,BMAT(1,J,IZ),1)
        END DO
      END DO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck matsct */
      SUBROUTINE MATSCT(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,NZ)
C***********************************************************************
C
C     Construct AMAT from BMAT elements using IDX.
C
C     Input:
C
C     Output:
C
C     Written by J. Thyssen - Mon Mar 19 16:41:48 "MET 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ), BMAT(NRB,NCB,NZ), IDX(*)
C
      DO IZ = 1, NZ
         DO J = 1, NCB
            JJ = IDX(J)
            DO K = 1, NRB
               KK = IDX(K)
               AMAT(KK,JJ,IZ) = BMAT(K,J,IZ)
            END DO
         END DO
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MATSCT_symm(AMAT,NRA,NCA,BMAT,NRB,NCB,IDX,
     &                       ibeig,is_symm,NZ)
C***********************************************************************
C
C     Construct AMAT from BMAT elements using IDX and ibeig.
C
C     Input:
C
C     Output:
C
C     Written by S. Knecht - August 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION AMAT(NRA,NCA,NZ), BMAT(NRB,NCB,NZ), IDX(*), ibeig(*)
C
      DO IZ = 1, NZ
         DO 5 J = 1, NCB
            JJ = IDX(J)
#ifdef MCSCF_debug
            print *, ' JJ jz ==> ',JJ, ibeig(jj)
#endif
            if(ibeig(jj).ne.is_symm) goto 5
            DO 10 K = 1, NRB
               KK = IDX(K)
#ifdef MCSCF_debug
               print *, ' KK jz ==> ',KK,ibeig(KK)
#endif
               if(ibeig(KK).ne.is_symm) goto 10
#ifdef MCSCF_debug
      print *,' put in b(k,j) value in a(KK,JJ) ==> ', 
     &                   k,j,BMAT(K,J,IZ),kk,jj
#endif
               AMAT(KK,JJ,IZ) = BMAT(K,J,IZ)
 10         continue
  5      continue
      END DO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      integer function iimax(n,iarr,inca)
C***********************************************************************
C
C     find largest absolute element in integer array iarr.
C
C     Input:
C
C     Output:
C
C     Written by S. Knecht - August 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION iarr(*)
C
      ismax = 0
      DO ielm = 1, n
        ismax = max(ismax,abs(iarr(ielm)))
      END DO
      iimax = ismax
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DETXFLAB */
      SUBROUTINE DETXFLAB(XSYMFLAB,IXSYM)
C***********************************************************************
C
C     determine KRCI_CVECS.x xSYMFLAB
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
C
      CHARACTER*3 XSYMFLAB(IXSYM) 
       
      do i = 1, IXSYM
        if(i.lt.10)then
          write(XSYMFLAB(I),'(i1,a2)') i,'  '
        else if(i.lt.100)then
          write(XSYMFLAB(I),'(i2,a1)') i,' '
        else
          write(XSYMFLAB(I),'(i3)')i
        end if
      end do
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck FIND_XPRPINDEX */
      INTEGER FUNCTION IXPRPINDEX(PRPNAMES,LPROP_KRCI,NPROP,XLABEL)
C***********************************************************************
C
C     find index of property XLABEL in KRCI format
C
C     Written by S. Knecht - Nov 2008
C
C     Last revision :
C
C***********************************************************************
#include "implicit.h"
C
      CHARACTER*16 PRPNAMES(*), XLABEL
      DIMENSION LPROP_KRCI(*)
C
      IXPRPINDEX = 0
      I_SCR      = 0
      DO I = 1, NPROP
        INDXPR = LPROP_KRCI(I)
        IF( PRPNAMES(INDXPR) .eq. XLABEL ) I_SCR = I
      END DO
      IXPRPINDEX = I_SCR
      END
