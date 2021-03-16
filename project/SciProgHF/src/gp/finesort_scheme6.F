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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routines for "fine-sorting" of the 1HT data, prior to the 2HT step 
c in scheme 6
c
c initial implementation: Andre Gomes, Amsterdam/Lille, april/may 2009
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c sort_buffers: driver routine for the sorting. at this time it assumes 
c               that we are *not* using the "MPI-IO" code that is 
c               enabled with -DVAR_PFS
c
c               currently the functionality is as such that the data is
c               sorted one "buffer" at a time and completely in-memory;
c               this means that a whole buffer is a) read to memory, 
c               b) sorted, and c) written back to disk along with 
c               information on how its contents are now stored, so that 
c               during the 2HT only the necessary records will be read
c               again
c
c               in progress is the incorporation of functionality to 
c               perform the sorting out-of-memory, with the intention 
c               that the user may choose the strategy best suited for
c               a given machine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine sort_buffers(thisnode,i_file)
          use memory_allocator
#include "implicit.h"
#include "bucket_sort.h"
#include "infpar.h"

          integer :: i_file, thisnode, verbosity_factor
          integer :: data_chunksize, buffer, idim3

          if (thisnode.eq.0) call print_resort_banner()

          data_chunksize  = NGBFSZ 

c aspg
c with verbosity factor we want to control the printing of
c of the sorting progress, in intervals of 25% now
c 
          verbosity_factor = (N_BUFFER/(NUMNOD+1)/4)

c aspg, hardcoded! change this later
          i_metadata_file = 1

c we find the maximum buffer size in terms of written records, and
c use that as scratch space
          max_records_per_buffer = -1

          do buffer = 1, N_BUFFER 
             if (thisnode.eq.node_for_buffer(buffer)) then
                if (max_records_per_buffer.lt.
     &              records_in_buffer(buffer)) then
                   max_records_per_buffer = records_in_buffer(buffer) 
                endif 
             endif
          enddo

          if (max_records_per_buffer.eq.0) then
             call quit('All zero-sized buffers is not possible!')
          endif
       
c more scratch allocation
          imdsz = (KLpairsInBuffer+1)
          isiz0 = (max_records_per_buffer+1)*data_chunksize
          isiz1 = isiz0
          isiz2 = isiz1*2

          call init_metadata_file(i_metadata_file)  
       
          do buffer = 1, N_BUFFER 
            if (thisnode.eq.node_for_buffer(buffer)) then

               call read_and_resort_in_memory(buffer,
     &                              i_file, i_metadata_file,
     &                              data_chunksize,
     &                              KLpairsInBuffer,imdsz,isiz1)

c aspg: note to whoever wants to make this input "prettier":
c       please please *do not* make this write formatted unless 
c       you are *absolutely* sure that it is not broken for
c       compilers other than the one you are using!!! 

               if (mod(buffer,verbosity_factor).eq.0) then
                  write (6,*) '   > sort done up to buffer (@node ',
     &                        thisnode,')',buffer,
     &                        ' (out of',N_BUFFER,'buffers)'
                  call flush(6)
               endif
            endif
         enddo

c         if (thisnode.eq.0) then
            call print_resort_stats(thisnode)
            call timer('1HT_sort',TIMSTR,TIMEND)
c         endif

      end subroutine sort_buffers

c
c print_resort_banner: auxiliary routine, prints banner at start of fine-sorting
c
      subroutine print_resort_banner()
         write (*,*) ' '
         write (*,*) '<<< Starting sort of 1HT integral file(s)'
         write (*,*) ' '
      end subroutine print_resort_banner

c
c print_resort_stats
c
      subroutine print_resort_stats(thisnode)
      
         integer :: thisnode
         
         write (*,*) ' '
c         write (*,*) '+---> Fine-sort statistics not available yet'
         write (*,*) '>>> 1HT sorting finished for node ',thisnode
c         write (*,*) ' '
      end subroutine print_resort_stats

c
c read_and_resort_in_memory: as the name implies, reads the records related
c                           to a given buffer to memory, sorts them so that
c                           1HT data with same composite index KL (row/columns
c                           of the two-electron supermatrix) are grouped 
c                           together, and writes the sorted buffer back to file.
c                           in order to minimize the use of disk space, which
c                           can be an issue for very large calculations, the
c                           sorted buffer exactly replaces the unsorted one on
c                           disk - meaning that the same records written before
c                           are reused, and their content overwritten.
c      
      subroutine read_and_resort_in_memory(iactive_buffer,i_file,
     &                           i_metadata_file,
     &                           data_chunksize,
     &                           kl_in_buffer,
     &                           imdsz,isiz1)
      use memory_allocator
#include "implicit.h"
#include "bucket_sort.h"
#include "infpar.h"
       integer, intent(in) :: i_file, i_metadata_file
       integer :: iactive_buffer
       integer :: kl_in_buffer
       integer :: total_entries
       integer :: active_record
       integer :: mdata_idx, indexb
       REAL*8  :: RGBUF(NGBFSZ)
       DIMENSION IGBUF(2,NGBFSZ)
       integer, allocatable :: metadata(:,:), b(:,:), b2(:,:), count(:)
       real*8,  allocatable :: a(:), a2(:)
       integer, allocatable :: records_in_use(:,:)  ! keep track of which records are used in this buffer

c
c aspg, 20090422
c       for some reason, iactive_buffer was always getting overwritten, so
c       i save it and use iactive as working variable... seems like a compiler
c       problem... so, if you want to change this, please test carefully!!!
c
       iactive=iactive_buffer

       LGFIL    = I_BF_UNIT+I_FILE
       IREC     = LGREC(iactive_buffer,I_FILE)

       if (irec.eq.0) return ! if there is nothing in this buffer we leave

       LMFIL   = I_MD_UNIT+1
       imdrec  = irec_md_lastwritten(iactive_buffer)

       call alloc(records_in_use,3,max_records_per_buffer,              &
     &            id='records_in_use')
       call izero(records_in_use,3*max_records_per_buffer)

       call alloc(metadata,3,imdsz, id='metadata')
c 
c the tradeoff of the "suck all data into memory" sorting strategy
c is that while it avoids a larger number of reading operations, as
c the complete 1HT dataset is read only once, it may be that
c the a and b (and a2 and b2) temporary arrays below might grow
c too large - they should be of size:
c
c  size = total size of the dataset / (number of buffers)
c
       call alloc(a,isiz1, id="a")
       call izero(a,isiz1)

       call alloc(b,2,isiz1, id="b")
       call izero(b,2*isiz1)
c
c for the key-indexed counting sorting algorithm we need 
c the extra storage below (count, a2, b2 - a2 and b2 can perhaps be
c made redundant by an in-place exchange..?)
c
       call alloc(count,imdsz+1, id="count")
       count(1)=1
       forall (i=2:imdsz+1) count(i) = 0
c      
c we don't initialize a2, b2 to zero because there will be no accumulation done in them; 
c rather, they will be assigned the values in a,b
c
       call alloc(a2,isiz1,   id="a2")
       call alloc(b2,2,isiz1, id="b2")
c
c using the "suck all data into memory" sorting strategy: the whole buffer goes 
c in memory and is subsequently sorted, and broken down into the original record
c layout on file.
c
      call put_buffer_in_memory(lgfil,irec,active_record,
     &                          records_in_use,igbuf,b,rgbuf,a,
     &                          total_entries,NGBFSZ,count,kl_in_buffer)
c
c now for the actual sorting...
c
c key-indexed is linear in all operations, so it should show good 
c performance in practice. 
c
       call keyindexcount_sort_in_memory(a, b, total_entries,
     &                                   a2,b2,count,kl_in_buffer)
c
c and now we want to break down A and B in such a way that they fit the
c layout of the read records, and we then overwrite those in the file
c
      indexb = 0
      mdata_idx = 1 

      call initialize_metadata(metadata,imdsz)
      
      do j = 1, active_record
         n     = records_in_use(1,j)
         irec  = records_in_use(2,j)
         irecn = records_in_use(3,j)

         do k = 1, n            ! number of elements in the record
            indexb = indexb+1   ! indexb tracks total number of elements
            
            IGBUF(1,k) = b(1,indexb)
            IGBUF(2,k) = b(2,indexb)
            RGBUF(k)   = a(indexb)
c            
c we change value of icol, so we have to mark which record we are still in
c and if there is more then one column in this record, we mark its beginning as well
c
            call set_metadata(irec,b,indexb,metadata)
         enddo
c         
c overwrite the unsorted buffer with the sorted one
c
         write(LGFIL,REC=irec) irecn,n,rgbuf,igbuf
         
         if ( irecn .eq. 0 ) then
            call save_metadata(lmfil,iactive,metadata,
     &                         irec_md_global_lastwritten,
     &                         irec_md_lastwritten, kl_in_buffer)

            call dealloc(count)
            call dealloc(a2)
            call dealloc(b2)
            call dealloc(a)
            call dealloc(b)
            call dealloc(metadata)
            call dealloc(records_in_use)

            return
         endif
         
      enddo
c aspg, strange bug! see comments at the top of this routine...
      iactive_buffer=iactive
      
      call dealloc(count)
      call dealloc(a2)
      call dealloc(b2)
      call dealloc(a)
      call dealloc(b)
      call dealloc(metadata)
      call dealloc(records_in_use)

      end subroutine read_and_resort_in_memory
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routines dealing with the sorting "metadata" follow
c
c "metadata" here indicates the information of where a given
c row/column of the two-electron supermatrix is located in the
c file holding the corresponding buffer. 
c 
c such information (always related to the fully sorted buffer)
c is used to jump directly to the first record that holds such
c row/column data, and also to stop reading this file when 
c there is no data of this type in file anymore
c
c the metadata is written to a direct access file (one per mpi process)
c whose records have the following layout:
c
c [integer #1: pointer to next record]
c [integer #2: number of active buffer]
c [integer #3: index of a given row/column relative to this buffer]
c [integer #4: check on #3]
c [integer #5: start record of relative row/column] 
c [integer #6: end record of relative row/column]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c init_metadata_file: used to create the direct access file
c
      subroutine init_metadata_file(i_file)
      
         use memory_allocator
#include "implicit.h"
#include "bucket_sort.h"
#include "infpar.h"
         CHARACTER*10 NAME*5,FNODE,LFNAM

         IRECL  = IRECLEN(0,6,0)

         WRITE (NAME,'(A4,I1)') 'META',I_FILE
         FNODE = LFNAM(NAME)
         LMFIL = I_MD_UNIT+1

         OPEN (UNIT=LMFIL,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL,
     &         iostat=istat)

         if (istat.ne.0) call quit('Error opening metadata file!')

         irec_md_global_lastwritten = 0

         do i=1,N_BUFFER
            irec_md_lastwritten(i) = 0
         enddo
      
      end subroutine init_metadata_file

c
c set_metadata: keeps track of the start and end records fot a given
c               row/column of the 2-electron supermatrix
c
      subroutine set_metadata(irec,b,indexb,metadata)
         integer, intent(inout) :: metadata(3,*)
         integer, intent(in)  :: irec, indexb, b(2,*)
         integer :: iprev

         iprev = indexb - 1
         md_idx_current = b(2,indexb)

         if (indexb==1) then
            metadata(2,md_idx_current) = irec
         elseif ( (indexb > 1) .and. (b(2,iprev) .ne. b(2,indexb))) then
            md_idx_prev    = b(2,iprev) 
            metadata(3,md_idx_prev)    = irec
            metadata(2,md_idx_current) = irec
         endif
         metadata(3,md_idx_current) = irec

      end subroutine set_metadata

c
c print_metadata: utility to show the content of par of the 
c                 metadata record
c
      subroutine print_metadata(metadata,mdata_idx)
         integer, intent(in) :: mdata_idx, metadata(3,*)
         integer :: i

         write(6,*) 'metadata information stored so far'
         do i = 1, mdata_idx
            write(6,'(4I12)') i,(metadata(j,i),j=1,3)
         enddo

      end subroutine print_metadata

c
c initialize_metadata: utility to initialize the content of part of the 
c                      metadata record
c
      subroutine initialize_metadata(metadata_array,metadata_size)
         integer, intent(in) :: metadata_size
         integer :: i, metadata_array(3,*)
         
         do i = 1, metadata_size
            metadata_array(1,i) = i
            metadata_array(2,i) = 0
            metadata_array(3,i) = 0
         enddo
         
      end subroutine initialize_metadata

c
c save_metadata: after a buffer has been sorted and all the start,end
c                records for the corresponding rows/columns have been
c                properly identified, we write all the metadata records 
c                related to this particular buffer to the metada file 
c                associated to an mpi process
c
      subroutine save_metadata(lmfil,iactive_buffer,metadata_array,
     &                         irec_md_global_lastwritten,
     &                         irec_md_lastwritten, KLpairsInBuffer)
         integer, intent(in) :: lmfil,iactive_buffer, KLpairsInBuffer
         integer :: l, imdrec, imdrecn
         integer :: metadata_array(3,*)
         integer :: irec_md_global_lastwritten, irec_md_lastwritten(*)

         
         do l = 1, KLpairsInBuffer

            if (l.eq.1) then
               imdrec  = 0
            else
               imdrec  = irec_md_global_lastwritten
            endif
               imdrecn = irec_md_global_lastwritten+1

            if ( metadata_array(1,l).ne.l ) then
               write(6,*) 'inconsistent kl? ',l,
     &                     metadata_array(1,l),
     &                     metadata_array(2,l),
     &                     metadata_array(3,l)

            endif

            write(LMFIL,REC=imdrecn) imdrec,iactive_buffer,l,
     &                               metadata_array(1,l),
     &                               metadata_array(2,l),
     &                               metadata_array(3,l)

            irec_md_lastwritten(iactive_buffer) = imdrecn
            irec_md_global_lastwritten   = imdrecn

         enddo

      end subroutine save_metadata
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routines for the actual data sorting
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c put_buffer_in_memory: this routine will read a given buffer to
c                       temporary storage in memory, 
      subroutine put_buffer_in_memory(lgfil,start_record,active_record,
     &                                records_in_use,igbuf,b,rgbuf,a,
     &                                total_entries,max_size,count,
     &                                KLpairsInBuffer)
         integer, intent(in) :: max_size, lgfil
         integer :: records_in_use(3,*), count(KLpairsInBuffer)
         integer, intent(out) :: total_entries
         integer :: igbuf(2,max_size), b(2,*)
         real*8  :: rgbuf(max_size), a(*)
         integer :: irec, start_record, active_record

         total_entries = 0
         active_record = 0

         irec = start_record
C      first we want to have all the data from the buffer in memory...
    1    if (irec.eq.0) then
            goto 2
         endif

           READ (LGFIL,REC=IREC) IRECN,N,RGBUF,IGBUF
           active_record = active_record+1

C     we first keep track of how much data each read has, becase later we will write back
c     to the same locations..

           records_in_use(1,active_record) = N
           records_in_use(2,active_record) = IREC
           records_in_use(3,active_record) = IRECN

           IREC = IRECN
C
C     first we put everyone for this buffer in memory 
C
           DO I = 1, N
              a(i+total_entries)   = RGBUF(I)
              b(1,i+total_entries) = IGBUF(1,I)
              b(2,i+total_entries) = IGBUF(2,I)    ! this is the ICOL equivalent
              count(IGBUF(2,I) + 1) = count(IGBUF(2,I) + 1) + 1  ! for key-indexed counting sort
c              idxtable(IGBUF(2,I),count(IGBUF(2,I) + 1)) = i+total_entries
           ENDDO
           total_entries = total_entries + N
C
C     We read in reverse order, if we have the first block we're done
C
         goto 1

C     so here we read everything, and kept the structure of the buffer in terms
c     of the records written... now what we want to reorder the contents of temp
c     variables A and B

    2    continue

c for key-indexed counting sort, this will set the final offsets

         isum = count(1)
         do j = 2, KLpairsInBuffer
            isum = isum + count(j)
            count(j) = count(j) + count(j-1) 
         enddo

c aspg, debug
         if (.false.) then
            write (*,*) ' ---> records in use for this buffer '
            do iriu = 1, active_record
               write (*,*) '      ',iriu,
     &                        '  ',records_in_use(1,iriu),
     &                        '  ',records_in_use(2,iriu),
     &                        '  ',records_in_use(3,iriu) 
            enddo

            do i = 1, KLpairsInBuffer
               write (6,*) ' ---> offsets for IGBUF(2) =',
     &                      i,' is ',count(i)
            enddo
            write (6,*) ' ---> total number of elements is ',isum
c aspg, debug
            call flush(6)
         endif

      end subroutine put_buffer_in_memory 





C soubroutines to sort buffer in-memory

       subroutine keyindexcount_sort_in_memory(a,b,total_entries,
     &               a2,b2,count,KLpairsInBuffer)
          integer :: total_entries, i, j, itmp, min, KLpairsInBuffer
          integer :: b(2,total_entries)
          integer :: b2(2,total_entries), count(KLpairsInBuffer)
          real*8  :: a(total_entries),   a2(total_entries)

          do i = 1, total_entries
             a2(   count(b(2,i) ))   = a(i)
             b2(1, count(b(2,i) ))   = b(1,i)
             b2(2, count(b(2,i) ))   = b(2,i)
             count( b(2,i) ) = count( b(2,i) ) + 1
          enddo

          do i = 1, total_entries
             a(i)   = a2(i) 
             b(1,i) = b2(1,i)
             b(2,i) = b2(2,i)
          enddo


       end subroutine keyindexcount_sort_in_memory


       subroutine selection_sort_in_memory(a,b,total_entries)
          integer :: total_entries, i, j, itmp, min
          integer :: b(2,total_entries)
          real*8  :: rtmp
          real*8  :: a(total_entries)

c below is a selection sort implementation...
          do i = 1, total_entries-1
             min = i
             do j = i+1, total_entries
                if ( b(2,j) .lt. b(2,min) ) min = j
             enddo   

             itmp      = b(2,i)
             b(2,i)    = b(2,min)
             b(2,min) = itmp

             itmp      = b(1,i)
             b(1,i)    = b(1,min)
             b(1,min) = itmp

             rtmp    = a(i)
             a(i)    = a(min)
             a(min) = rtmp

          enddo
          
       end subroutine selection_sort_in_memory
      
      subroutine print_buffer(a,b,total_entries)
#include "priunit.h"
         integer :: total_entries, i
         integer :: b(2,total_entries)
         real*8  :: a(total_entries)

         do i = 1, total_entries
            write (LUPRI,'(2X,3I1,F12.8)') i,b(1,i),b(2,i),a(i)
         enddo
      
      end subroutine print_buffer



      subroutine init_finesort_commons(kl_in_buf,NumBuffers,BuffNodes)
#include "implicit.h"
#include "bucket_sort.h"

         integer :: kl_in_buf, NumBuffers, BuffNodes(NumBuffers)
         integer :: ibuf

         KLpairsInBuffer = KL_IN_BUF

         do ibuf=1,N_BUFFER
            records_in_buffer(ibuf) = 0
            node_for_buffer(ibuf)   = BuffNodes(ibuf)
         enddo

      end subroutine init_finesort_commons
c
#ifdef SCHEME6_FINESORT_INFILE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routines related to sorting the file to disk instead of in-memory
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c init_finesort_files is similar in spirit to buffer_file_init(), with
c the difference that we use a different seed for the file units, and
c we open all (=one for each kl in the buffer) files from this call
c
       subroutine create_finesort_files(buffer,
     &                                  kl_in_buffer,
     &                                  i_finesort_lastrec,
     &                                  i_finesort_rec)
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#if defined (VAR_MPI)
#include "mpif.h"
#endif
          CHARACTER*10 NAME*5,FNODE,LFNAM
          integer :: buffer, kl_in_buffer
          integer :: i_finesort_rec(kl_in_buffer)       ! saves current record of file holding a given kl
          integer :: i_finesort_lastrec(kl_in_buffer)   ! saves last record of file holding a given kl
C
C     The record size needs to be NGBFSZ real + (2*NGBFSZ + 2) integer
C     because we also need to store the labels
C
          IRECL = IRECLEN(NGBFSZ+1,2*NGBFSZ+2,0)
                  
          do i = 1, kl_in_buffer
             l_finesort_file = i_finesort_unit + i      

             WRITE (NAME,'(A2,I3)') 'fs',i_finesort_file
             FNODE = LFNAM(NAME)

             OPEN (UNIT=l_finesort_file,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
            
             i_finesort_lastrec(i) = 0
             i_finesort_name(i)    = FNODE
             i_finesort_rec(i)     = 0
            
          enddo

      end subroutine create_finesort_files

c
c write_finesort_file: analogous to buffer_file_write, as far as
c its task is to get a chunck of data and write to a given file.
c as with init_finesort_files, the major difference here is how
c one keeps track of the records written, the file unit and so on
c 
       subroutine write_finesort_file(buffer,
     &                                kl, kl_in_buffer,
     &                                kl_in_buffer,
     &                                i_finesort_lastrec,
     &                                i_finesort_rec,
     &                                rgbuf,
     &                                igbuf,
     &                                n)
      
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
          real(kind=8) :: rgbuf(NGBFSZ)
          integer      :: igbuf(2,NGBFSZ)
          integer      :: kl, buffer, kl_in_buffer, n
          integer      :: i_finesort_rec(kl_in_buffer)      ! next to last record written...
          integer      :: i_finesort_lastrec(kl_in_buffer)  ! last record written
C
C     Get file unit number and file pointer
C
          lgfil    = i_finesort_unit + kl
          irec     = i_finesort_lastrec(kl) + 1
          WRITE (LGFIL,REC=IREC) i_finesort_lastrec(kl),N,RGBUF,IGBUF
C
C     Update file pointers
C
          i_finesort_lastrec(kl) = irec
c         i_finesort_rec(kl)     = irec
         
       end subroutine write_finesort_file


c need a flush routine? i think not - aspg


       subroutine read_finesort_file(buffer,
     &                              kl, 
     &                              kl_in_buffer,
     &                              i_finesort_lastrec,
     &                              i_finesort_rec,
     &                              rgbuf,
     &                              igbuf,
     &                              n)
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
          real(kind=8) :: rgbuf(NGBFSZ)
          integer      :: igbuf(2,NGBFSZ)
          integer      :: kl, buffer, kl_in_buffer, n
          integer      :: i_finesort_rec(kl_in_buffer), i_finesort_lastrec(kl_in_buffer)
     
                      
          
       end subroutine read_finesort_file
      
      
       subroutine sort_buffer_on_disk(iactive_buffer,
     &                                i_file,
     &                                i_metadata_file,
     &                                data_chunksize,
     &                                kl_in_buffer,
     &                                imdsz,isiz1)
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
          integer      :: iactive_buffer
          integer,      allocatable :: irec(:), irecn(:)   

          
          call alloc(irec, kl_in_buffer, id="irec")
          call alloc(irecn,kl_in_buffer, id="irecn")
          
          call create_filesort_files(iactive_buffer,kl_in_buffer,irec,irecn)

          call process_unsorted_buffer(iactive_buffer,
     &                                 i_file,     
     &                                 kl_in_buffer,
     &                                 irec,
     &                                 irecn)
      
          call dealloc(finesort_buffer_real)
          call dealloc(finesort_buffer_int)
          
          call build_finesorted_buffer(iactive_buffer,kl_in_buffer)


c when we are done we get rid of these files...
          call destroy_filesort_files(iactive_buffer,kl_in_buffer)
          
          call dealloc(irec)
          call dealloc(irecn)
                             
       end subroutine sort_buffer_on_disk
       
       
       subroutine build_finesorted_buffer(iactive_buffer,kl_in_buffer)
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
          integer      :: iactive_buffer

          start from the top of the storage layout for this buffer
          set start,end metadata information 
          
          loop over all kl
             set start, end metadata for this kl
             
  22         read (data from splitfile for current)
  
             if (number of elements is the same as the size of storage)
                write value of temp storage to record belonging to final file
                update record count, that is, go to next record of final file
                update metadata for final file 
             else
                if (number of elements is larger than storage size) then
                   write part of the data to the record
                   move remainder to second temp storage
                   keep track of how many elements of this kl are in temp 
                   storage, and the start index 
                else
                   
                endif
             endif
   
  22         read (data from split file for kl=kl)
          
             if (last record of splot file is not zero)
            
                if (number of elements is the same as 
          
                goto 22
                
             endif   
             
             update metadata
          end loop over all kl       
          
       end subroutine build_finesorted_buffer
       
       
       subroutine process_unsorted_buffer(iactive_buffer,i_1ht_file,kl_in_buffer,irec_fs,irecn_fs)
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
          integer      :: iactive_buffer, kl_in_buffer, i_1th_file
          integer      :: irec_fs(kl_in_buffer),irecn_fs(kl_in_buffer)
          real(kind=8), allocatable :: sortbuff_real(:), rgbuf(:)
          integer     , allocatable :: sortbuff_int(:,:), igbuf(:,:)      

c we dynamically allocate the storage used here
          call alloc(rgbuf,ngbfsz, id="rgbuf")
          call alloc(igbuf,2,ngbfsz, id="igbuf")
          call alloc(sortbuff_real,kl_in_buffer,ngbfsz,id="sortbuff_real")
          call alloc(sortbuff_int ,kl_in_buffer,2,ngbfsz,id="sortbuff_int")
          call alloc(elements_in_sortbuff,kl_in_buffer,id="elements_in_sortbuff")
        
c initializations
        
          lgfil_1ht_buff = I_BF_UNIT+i_1ht_file
          irec           = LGREC(iactive_buffer,I_FILE)
          call izero(elements_in_sortbuff,kl_in_buffer)

c
c going over the 1HT buffer file now...
c         
   11     call dzero(rgbuf,ngbfsz)
          call izero(igbuf,2*ngbfsz)

c read buffer chunk, update next record pointer
          read (lgfil_1ht_buff,rec=irec) irecn,n,rgbuf,igbuf
          irec=irecn

c now we separate each kl into its "sub"-buffer, effectively sorting the 
c 1HT data as we move along...
          do i = 1, n
             kl = igbuf(2,1)
             elements_in_sortbuff(kl) = elements_in_sortbuff(kl) + 1
             sortbuff_real(kl,elements_in_sortbuff(kl))  = rgbuf(i)
             sortbuff_int(kl,1,elements_in_sortbuff(kl)) = igbuf(1,i)
             sortbuff_int(kl,2,elements_in_sortbuff(kl)) = igbuf(2,i)
          enddo
          
c if we reach the last record of the 1HT buffer, we dump things to file ...        
          if (irecn.eq.0) then
             do kl = 1,kl_in_buffer
                call write_finesort_file(kl,elements_in_sortbuff(kl),sortbuff_real(kl),sortbuff_int(kl))
             enddo
            
c ... and explicitly deallocate things, since we can't really trust the
c compilers to reclaim the allocatables on exit... blame the vendors...             
             call dealloc(rgbuf)
             call dealloc(igbuf)
             call dealloc(sortbuff_real)
             call dealloc(sortbuff_int)
             call dealloc(elements_in_sortbuff)
             
             return
          endif
          
c .. but otherwise, dump things to file as we reach the amount of data
c per record
          if (elements_in_sortbuff(kl).eq.ngbfsz) then
             call write_finesort_file(kl,ngbfsz,sortbuff_real(kl),sortbuff_int(kl))
          endif
          
          goto 11
       
          call quit("we should never reach this point!")
          
       end subroutine process_unsorted_buffer
#endif
