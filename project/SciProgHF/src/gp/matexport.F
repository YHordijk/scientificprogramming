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

c
c     These subroutine write matrices to the file defined by
c     MATLAB_OUTFILE, appending new matrices to an existing file.  The
c     format is ascii with formatting suitable for automatic reading by
c     Matlab (or the free Octave) program. The output file is opened and
c     closed automatically. -ulfek
c


#define MATLAB_OUTFILE 'matexport.m'
#define MATLAB_LU 33
#define MATLAB_DFMT 'E24.16'
#define MATLAB_IFMT 'I16'

      subroutine matexport_text(text)
      implicit none
      character*(*) text
      open(unit=MATLAB_LU,file=MATLAB_OUTFILE,
     &     form='FORMATTED',position='append')
      write(MATLAB_LU,*) text
      close(MATLAB_LU)
      end

      subroutine matexport_integer(label,x,len)
      implicit none
      character*(*) label
      integer x,len,i
      dimension x(len)
      open(unit=MATLAB_LU,file=MATLAB_OUTFILE,
     &     form='FORMATTED',position='append')
      if (len.eq.1) then
         write (MATLAB_LU,'(1A,'//MATLAB_IFMT//',A)') label//'=',
     &        x(1),';'
      else
         write (MATLAB_LU,*) label//'= ['
         do i=1,len
            write (MATLAB_LU,'('//MATLAB_IFMT//')') x(i)
         enddo
         write (MATLAB_LU,*) '];'
      endif
      close(MATLAB_LU)
      end

      subroutine matexport_double(label,x,len)
      implicit none
      double precision x
      character*(*) label
      integer len,i
      dimension x(len)
      open(unit=MATLAB_LU,file=MATLAB_OUTFILE,
     &     form='FORMATTED',position='append')
      if (len.eq.1) then
         write (MATLAB_LU,'(1A,'//MATLAB_DFMT//',A)') label//'=',
     &        x(1),';'
      else
         write (MATLAB_LU,*) label//'= ['
         do i=1,len
            write (MATLAB_LU,'('//MATLAB_DFMT//')') x(i)
         enddo
         write (MATLAB_LU,*) '];'
      endif
      close(MATLAB_LU)
      end

      subroutine matexport_double2(label,x,lda,nrow,ncol)
      implicit none
      double precision x
      character*(*) label
      integer lda,nrow,ncol,i,j
      dimension x(lda,ncol)
      open(unit=MATLAB_LU,file=MATLAB_OUTFILE,
     &     form='FORMATTED',position='append')
C     Write as 1d vector and then reshape.
      if (nrow*ncol.eq.1) then
         write (MATLAB_LU,'(1A,'//MATLAB_DFMT//',A)') label//'=',
     &        x(1,1),';'
      else
         write (MATLAB_LU,*) label//'= ['
         do j=1,ncol
            do i=1,nrow
               write (MATLAB_LU,'('//MATLAB_DFMT//')') x(i,j)
            enddo
         enddo
         write (MATLAB_LU,*) '];'
      endif
      write (MATLAB_LU,*) label//'=reshape('//label//','
     &     ,nrow,',',ncol,');'
      close(MATLAB_LU)
      end

      subroutine matexport_double3(label,x,lda,ldb,nrow,ncol,nz)
      implicit none
      double precision x
      character*(*) label
      integer lda,ldb,nrow,ncol,nz,i,j,k
      dimension x(lda,ldb,nz)
      open(unit=MATLAB_LU,file=MATLAB_OUTFILE,
     &     form='FORMATTED',position='append')
C     Write as 1d vector and then reshape.
      if (nrow*ncol*nz.eq.1) then
         write (MATLAB_LU,'(1A,'//MATLAB_DFMT//',A)') label//'=',
     &        x(1,1,1),';'
      else
         write (MATLAB_LU,*) label//'= ['
         do k=1,nz
            do j=1,ncol
               do i=1,nrow
                  write (MATLAB_LU,'('//MATLAB_DFMT//')') x(i,j,k)
               enddo
            enddo
         enddo
         write (MATLAB_LU,*) '];'
      endif
      write (MATLAB_LU,*) label//'=reshape('//label//',',nrow,','
     &     ,ncol,',',nz,');'
      close(MATLAB_LU)
      end

#ifdef MATEXPORT_TEST
      program test
      implicit none
      double precision x,y,z
      integer iv
      dimension x(3),y(2,2),z(2,2,2),iv(2);
      iv(1) = 17
      iv(2) = 23

      x(1) = 7
      x(2) = 12
      x(3) = -1
      y(1,1) = 11;
      y(2,1) = 21;
      y(1,2) = 12;
      y(2,2) = 22;

      z(1,1,1) = 11;
      z(2,1,1) = 21;
      z(1,2,1) = 12;
      z(2,2,1) = 22;
      z(1,1,2) = 311;
      z(2,1,2) = 321;
      z(1,2,2) = 312;
      z(2,2,2) = 322;
      call matexport_integer('iv',iv,2)
      call matexport_double('x3',x,3)
      call matexport_double('x2',x,2)
      call matexport_double('x1',x,1)
      call matexport_double2('y',y,2,2,2)
      call matexport_double3('z',z,2,2,2,2,2)
      end
#endif
