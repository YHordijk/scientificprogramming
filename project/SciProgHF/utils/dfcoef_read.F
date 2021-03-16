!      Copyright (c) 2018 by the authors of DIRAC.
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

c     program to dump data from a DFCOEF file.
c     Usage: dfcoef_read FILE
c     where FILE is the path to a DFCOEF file.
C     By Ulf Ekstrom 2009.

c     this is a cut down version of REACMO (dirac/dirgp.F)
c     -Ulf
      SUBROUTINE DUMPCMO(IUNIT)
      implicit none
      CHARACTER TEXT*74
      integer iunit,idim,i,j,k,nsym,pos1, totorbs,
     &     ftell, marker_size
      double precision toterg,eigs
      DIMENSION IDIM(3,2),eigs(100000)
C
      REWIND IUNIT
C
C     Read title line
C     ===============
C
      READ (IUNIT,END=10,ERR=20) TEXT,NSYM,
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
      print *,'TEXT ',TEXT
      print *,'NSYM',NSYM
      do i=1,NSYM
         print *,'In fermion symmetry ',i
         print *,'  Positronic orbitals', IDIM(1,I)
         print *,'  Electronic orbitals', IDIM(2,I)
         print *,'  AO basis function  ', IDIM(3,I)
      enddo
      print *,'TOTERG',TOTERG
      
      totorbs = 0
      do i=1,NSYM
         totorbs = totorbs + IDIM(1,I)+IDIM(2,I)
      enddo

      pos1 = ftell(iunit)
      if (pos1.ne.94+12*NSYM) then
         print *,'Sorry, unknown on-disk format, I cannot'//
     &        ' guess NZ value.'
         stop 'Unknown binary format'
      endif
c     Number of bytes per "marker" surrounding each file section
c     Typically 4 (or 8 in 64 bit formats). Will be used to guess NZ
c     for reading mo coefficients.
      marker_size = 4

c     for now we just skip down to the eigenvalues
      print *,'Skipping mo coefficients'
      READ (IUNIT,END=10,ERR=20)

      read (iunit) (eigs(i),i=1,totorbs)
      print *,'Eigenvalues:'
      k = 1
      do i=1,nsym       
         if (idim(1,i).gt.0) then
            print *,' Positronic (sym =',i,')'
            do j=1,idim(1,i)
               print *,j,eigs(k)
               k = k + 1
            enddo
         endif
         print *,' Electronic (sym =',i,')'
         do j=1,idim(2,i)
            print *,j,eigs(k)
            k = k + 1
         enddo
      enddo
      return
 10   continue
      print *,'End of file while looking for data.'
 20   continue
      print *,'Error reading DFCOEF file'
      return
      END

      PROGRAM DFCOEF_READ
      implicit none
      INTEGER   NARG,JARG,LUINP,LUPRI
      CHARACTER*80 ARGC
      PARAMETER (LUINP = 1, LUPRI = 6)
      NARG = command_argument_count()
      IF (NARG .NE. 1) THEN
         WRITE(6,'(/3A/)') 'Usage: dfcoef_read FILE'
         stop
      END IF
      CALL get_command_argument(1,ARGC)
      OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
      call dumpcmo(luinp)
      END
