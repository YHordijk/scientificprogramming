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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FMTOFILE (PUT,FVO,FOO,FVV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Read/write of the Fock matrix: FOO, FVO and FVV parts
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 FVO(*),FOO(*),FVV(*)
      LOGICAL PUT
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ccpar.inc"
#include "waio.h"
C
C---------------Local variables--------------------------------------
C
      integer nsec, nsec0, nsec1, nsec2
C
C---------------Executable code--------------------------------------
C
C     Calculate start of sectors
C
      NSEC0 = 1
      NSEC = NFVO*IRW*RCW / NWORDS_ON_RECORD
      IF (MOD(NFVO*IRW*RCW,NWORDS_ON_RECORD).NE.0) NSEC = NSEC + 1
      NSEC1 = NSEC0 + NSEC
      NSEC = NFOO*IRW*RCW / NWORDS_ON_RECORD
      IF (MOD(NFOO*IRW*RCW,NWORDS_ON_RECORD).NE.0) NSEC = NSEC + 1
      NSEC2 = NSEC1 + NSEC
C
C     Read/write blocks
C
      IF (PUT) THEN
         CALL RWRIT(ITAPT+0,FVO,NFVO*RCW*IRW,NSEC0)
         CALL RWRIT(ITAPT+0,FOO,NFOO*RCW*IRW,NSEC1)
         CALL RWRIT(ITAPT+0,FVV,NFVV*RCW*IRW,NSEC2)
      ELSE
         CALL RREAD(ITAPT+0,FVO,NFVO*RCW*IRW,NSEC0)
         CALL RREAD(ITAPT+0,FOO,NFOO*RCW*IRW,NSEC1)
         CALL RREAD(ITAPT+0,FVV,NFVV*RCW*IRW,NSEC2)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETOOOO (VOOOO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get OOOO part of two-electron integrals
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VOOOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL RREAD(ITAPT+1,VOOOO,NV1*IRW*RCW,1)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETVOOO (VVOOO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOOO part of two-electron integrals
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVOOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL RREAD(ITAPT+2,VVOOO,NV2*IRW*RCW,1)
C
      RETURN
      END

    
!! Avijit Shee: 24th March, 2015.

      SUBROUTINE GETOOVO (VOOVO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOOO part of two-electron integrals and Sort it to OOVO
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8             :: VOOVO(*)
      REAL*8,ALLOCATABLE :: VVOOO(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
      allocate (VVOOO(NV2*IRW*RCW))

      CALL RREAD(ITAPT+2,VVOOO,NV2*IRW*RCW,1)
       
      CALL SRT1C1 (NREP,NVO,NOOT,VVOOO,VOOVO) 

      deallocate(VVOOO)

      RETURN
      END


!!---------------------------------------------------------------------------

      SUBROUTINE GETOVOO (VOVOO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOOO part of two-electron integrals and Sort it to OOVO
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8             :: VOVOO(*)
      REAL*8,ALLOCATABLE :: VVOOO(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
      allocate (VVOOO(NV2*IRW*RCW))

      CALL RREAD(ITAPT+2,VVOOO,NV2*IRW*RCW,1)

      CALL SRT1L1(NREP,MULTB,.FALSE.,NVO,NV,NO,NOOT,IOVOOT,IIOV,
     &              VVOOO,VOVOO)

      deallocate(VVOOO)

      RETURN
      END

!=======================================================================
      SUBROUTINE GETOOOV (VOOOV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOOO part of two-electron integrals and Sort it to OOVO
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8             :: VOOOV(*)
      REAL*8,ALLOCATABLE :: VOVOO(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
      allocate (VOVOO(IOVOOT(NREP+1)*IRW*RCW))

      CALL GETOVOO (VOVOO) 

      CALL SRT1C1 (NREP,NOV,NOOT,VOVOO,VOOOV) 

      deallocate(VOVOO)

      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETVVOO (VVVOO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VVOO part of two-electron integrals
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVVOO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL RREAD(ITAPT+3,VVVOO,NV3*IRW*RCW,1)
C
      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE GETOOVV (VOOVV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VVOO part of two-electron integrals
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 :: VOOVV(*)
      REAL*8,ALLOCATABLE :: VVVOO(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      ALLOCATE (VVVOO(NV3*IRW*RCW))

      CALL RREAD(ITAPT+3,VVVOO,NV3*IRW*RCW,1)

      CALL SRT1C1 (NREP,NVVT,NOOT,VVVOO,VOOVV) 

      DEALLOCATE (VVVOO)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETVOVO (VVOVO)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOVO part of two-electron integrals
C
C---------------Routines called----------------------------------------
C



C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVOVO(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      CALL RREAD(ITAPT+4,VVOVO,NV4*IRW*RCW,1)
C
      RETURN
      END


!=======================================================================

      SUBROUTINE GETVOOV (VVOOV)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOVO part of two-electron integrals
C
C---------------Routines called----------------------------------------


C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8             :: VVOOV(*)
      REAL*8,ALLOCATABLE :: VVOVO(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C

      ALLOCATE(VVOVO(NV4*IRW*RCW))

      CALL RREAD(ITAPT+4,VVOVO,NV4*IRW*RCW,1)

      CALL SRT1R1 (NREP,MULTB,.FALSE.,NVO,NV,NO,NVO,IVOOV,IIOV,
     &              VVOVO,VVOOV)

      DEALLOCATE(VVOVO)

      RETURN
      END


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETVOVV (IREP,ISTART,NINT,DONE,VVOVV,MINT)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOVV part of two-electron integrals in the
C     Representation                        IREP
C     with start address of integrals       ISTART (updated)
C     Number of integrals in batch          NINT (updated)
C     if there are read all in a batch      DONE=.TRUE.
C     Buffer:                               VVOVV(*)
C     Maximum number of VO batches that can
C     be buffered:                          MINT
C     
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVOVV(*)
      INTEGER IREP,ISTART,NINT,MINT
      LOGICAL DONE
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "waio.h"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
       INTEGER*8 JSTART
       integer jint,mtot,ntot
C
C---------------Executable code--------------------------------------
C
      DONE = .FALSE.
      MTOT = NVO(IREP)     ! one (VV)_ij batch has length NVO
      IF (IDIST(1,2,IREP).GT.ISTART) ISTART=IDIST(1,2,IREP)
      NTOT = IDIST(2,2,IREP)     ! local # of VO batches in this IRREP
      NINT = MIN0(MINT,NTOT-ISTART)  ! # of batches which will actually be read
      JSTART = IDIST(5,2,IREP)      ! compute relative local offset
      JSTART = JSTART + (ISTART-IDIST(1,2,IREP)) * MTOT
C
C     JINT : Total number of integrals to be read
C
      JINT = MTOT * NINT
C
C     The integrals may be complex, convert offset and length to real before reading
C
      JINT = JINT * RCW
      JSTART = JSTART * RCW
C
      CALL WAIO_REALIO (ITAPT+5,VVOVV,JINT,JSTART,.TRUE.)
C
C     Are we done with this IRREP ?
C
      IF (ISTART+NINT.EQ.NTOT) DONE = .TRUE.
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GETVVVO_incore (VVVVO)
C
      implicit none
C
C---------------Description--------------------------------------------
C     Get VOVO part of two-electron integrals

C---------------Calling variables--------------------------------------
C
      REAL*8             :: VVVVO(*)
      REAL*8, ALLOCATABLE:: VVOVV(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
      REAL*8, ALLOCATABLE :: VOVV_TEMP(:)
      INTEGER, ALLOCATABLE :: counter(:), displ(:)
      LOGICAL :: DONE
      INTEGER :: totsize
C---------------Executable code--------------------------------------
C
      ALLOCATE(VVOVV(NV5*IRW*RCW))

      CALL GETVOVV_incore (VVOVV)

      CALL SRT1C1 (NREP,NVO,NVVT,VVOVV,VVVVO) 

      DEALLOCATE(VVOVV)

      RETURN
      END
!=================================================================

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GETVOVV_incore (VVOVV)

      use interface_to_mpi
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VOVV part of two-electron integrals in core.
C     We first read the VOVV integrals on different nodes.
C     Then synchronize all of them to the master 
C
C---------------Routines called----------------------------------------


C---------------Last modified------------------------------------------
C
C     Author : Avijit Shee
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVOVV(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "param.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
      integer ISTAT(df_MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
      REAL*8, ALLOCATABLE :: VOVV_TEMP(:)
      INTEGER, ALLOCATABLE :: counter(:), displ(:)
      LOGICAL :: DONE
      INTEGER :: totsize, check, i, istart, mint, nint, off1, irp
C---------------Executable code--------------------------------------
C
       allocate(counter(0:nmproc-1))
       allocate(displ(0:nmproc-1))

         counter = 0

       OFF1 = 1
       DO  30 IRP = 1, NREP
         check = IRP
         totsize = 0
         IF ((NVO(IRP).EQ.0).OR.(NVVT(IRP).EQ.0)) GOTO 30
         ALLOCATE(VOVV_TEMP(NVO(IRP)*NVVT(IRP)*RCW))

         VOVV_TEMP = 0.0d0

         DONE = .FALSE.
         ISTART = 0
         MINT = NVVT(IRP)  ! max. # of NVO batches that fit in mem.
  10     CONTINUE
         CALL GETVOVV (IRP,ISTART,NINT,DONE,VOVV_TEMP,MINT)

#if defined (VAR_MPI)
       call interface_mpi_allgather(NINT,1,
     &   counter,1,global_communicator)
#endif

      do i = 0, nmproc-1
       counter(i) = counter(i)*nvo(irp)*rcw            
      enddo

         displ(0) = 0
         do i = 1, nmproc-1 
             displ(i) = displ(i-1) + counter(i-1)
         enddo

#if defined (VAR_MPI)
      IF (nmproc .GT. 1) THEN

        CALL XTIME(4,1,'-- combining VOVV via MPI_GATHER')

        call interface_mpi_allgatherv(vovv_temp(1),counter(myproc),
     &    VVOVV(totsize+off1),counter,displ,global_communicator)

        CALL XTIME(4,2,'-- combining VOVV via MPI_GATHER')

      ENDIF
#endif

      if (nmproc .eq. 1) then

       call xcopy(nvo(irp)*nvvt(irp),vovv_temp(1),1,VVOVV(off1),1)
    
      endif

         totsize = totsize + sum(counter)

         if (.not.done) then
            istart = istart + nint
            goto 10
         endif
         off1 = off1+nvo(irp)*nvvt(irp)*rcw

         deallocate(vovv_temp)

 30    continue

       deallocate(counter)
       deallocate(displ)

      RETURN
      END

      SUBROUTINE GETOVVV_incore (VOVVV)
C
      implicit none
C
C---------------Description--------------------------------------------

C     Get VOVO part of two-electron integrals

C---------------Calling variables--------------------------------------
C
      REAL*8             :: VOVVV(*)
      REAL*8, ALLOCATABLE:: VVOVV(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C

      ALLOCATE(VVOVV(NV5*IRW*RCW))

      CALL GETVOVV_incore (VVOVV)

      CALL SRT1L1(NREP,MULTB,.FALSE.,NVO,NV,NO,NVVT,IOVVVT,IIOV,
     &             VVOVV,VOVVV)

      DEALLOCATE(VVOVV)

      RETURN
      END


      SUBROUTINE GETVVOV_incore (VVVOV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------

C     Get VOVO part of two-electron integrals

C---------------Calling variables--------------------------------------
C
      REAL*8             :: VVVOV(*)
      REAL*8, ALLOCATABLE:: VOVVV(:)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "param.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C

      ALLOCATE(VOVVV(NV5*IRW*RCW))

      CALL GETOVVV_incore (VOVVV)

      CALL SRT1C1 (NREP,NOV,NVVT,VOVVV,VVVOV) 

      DEALLOCATE(VOVVV)

      RETURN
      END


      SUBROUTINE GETVVVV (IREP,ISTART,NINT,DONE,VVVVV,MINT)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Get VVVV part of two-electron integrals
C     Representation 			IREP
C     Start address of integrals 	ISTART
C     Maximum number of integrals 	MINT
C     Actual number of integrals  	NINT (updated on output)
C     Done with this representation ?	DONE
C     
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 VVVVV(*)
      INTEGER IREP,ISTART,NINT,MINT
      LOGICAL DONE
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "waio.h"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
       INTEGER*8 JSTART
       integer jint,mtot,ntot
C
C---------------Executable code--------------------------------------
C
C     Calculate start address. The batches are distributed over the
C     nodes for parallel runs. Take only the part that is on our node
C
      DONE = .FALSE.
      MTOT = NVVT(IREP)
      IF (IDIST(1,1,IREP).GT.ISTART) ISTART=IDIST(1,1,IREP)
      NTOT = IDIST(2,1,IREP)
      NINT = MIN0(MINT,NTOT-ISTART)
      JSTART = IDIST(5,1,IREP)
      JSTART = JSTART + (ISTART-IDIST(1,1,IREP)) * MTOT

C     JINT : Total number of integrals to be read

      JINT = MTOT * NINT

C     The integrals may be complex, convert offset and length to real before reading

      JINT = JINT * RCW
      JSTART = JSTART * RCW

      CALL WAIO_REALIO (ITAPT+6,VVVVV,JINT,JSTART,.TRUE.)

C     Are we done with this IRREP ?

      IF (ISTART+NINT.EQ.NTOT) DONE = .TRUE.

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PUTDZ (IREP,ISTART,NINT,DZ)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Puts D intermediate to file
C     Representation 			IREP
C     Start address of integrals 	ISTART 
C     Number of integrals in batch 	NINT
C     
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Avijit Shee
C
C---------------Calling variables--------------------------------------
C
      REAL*8 DZ(*)
      INTEGER IREP,ISTART,NINT
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "eqns.inc"
#include "ccpar.inc"
C---------------Local variables--------------------------------------
       INTEGER*8 JSTART
       integer jint,mtot
C---------------Executable code--------------------------------------
      MTOT = NVVT(IREP)
      JSTART = IDIST(5,1,IREP)
      JSTART = JSTART + (ISTART-IDIST(1,1,IREP)) * MTOT

C     JINT : Total number of integrals to write

      JINT = MTOT * NINT

C     The integrals may be complex, convert offset and length to real before reading

      JINT = JINT * RCW
      JSTART = JSTART * RCW

      CALL WAIO_REALIO (ITAPT+6,DZ,JINT,JSTART,.FALSE.)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETDZ (IREP,ISTART,NINT,DONE,DZ,MINT)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Gets D intermediate from file
C     Representation 			IREP
C     Start address of integrals 	ISTART 
C     Number of integrals in batch 	NINT
C     
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Avijit Shee
C
C---------------Calling variables--------------------------------------
C
      REAL*8 DZ(*)
      INTEGER IREP,ISTART,NINT,MINT
      LOGICAL DONE
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
#include "waio.h"
#include "ccpar.inc"
C
C---------------Local variables--------------------------------------
C
       INTEGER*8 JSTART
      integer jint,mtot,ntot
C
C---------------Executable code--------------------------------------
C
C     Calculate start address. The batches are distributed over the
C     nodes for parallel runs. Take only the part that is on our node
C
      DONE = .FALSE.
      MTOT = NVVT(IREP)
      IF (IDIST(1,1,IREP).GT.ISTART) ISTART=IDIST(1,1,IREP)
      NTOT = IDIST(2,1,IREP)
      NINT = MIN0(MINT,NTOT-ISTART)
      JSTART = IDIST(5,1,IREP)
      JSTART = JSTART + (ISTART-IDIST(1,1,IREP)) * MTOT
C
C     JINT : Total number of integrals to be read
C
      JINT = MTOT * NINT
C
C     The integrals may be complex, convert offset and length to real before reading
C
      JINT = JINT * RCW
      JSTART = JSTART * RCW
C
      CALL WAIO_REALIO (ITAPT+6,DZ,JINT,JSTART,.TRUE.)
C
C     Are we done with this IRREP ?
C
      IF (ISTART+NINT.EQ.NTOT) DONE = .TRUE.

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PUTAMPT (T1,T2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Write amplitudes to restart file
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 BUFF(RSRCLEN)
      DATA BUFF /RSRCLEN*0.D0/
      integer ioff,irec,irecl,n,nbuf,ireclen
C
C---------------Executable code--------------------------------------
C
      IRECL = IRECLEN (RSRCLEN,0,0)
      OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',STATUS='OLD',
     &      RECL=IRECL)
C
      NRECT1 = NDIMT1*RCW / RSRCLEN + 1
      IRECT2 = IRECT1 + NRECT1 ! T2 amplitudes
      NRECT2 = NDIMT2*RCW / RSRCLEN + 1
C
      IREC = IRECT1
      N = NDIMT1
      IOFF = 1
    1 NBUF = MIN0(N,(RSRCLEN/RCW))
      CALL XCOPY (NBUF,T1(IOFF),1,BUFF,1)
      WRITE (MCCRES,REC=IREC) BUFF
      IREC = IREC + 1
      N = N - NBUF
      IOFF = IOFF + NBUF * RCW
      IF (N.GT.0) GOTO 1
C
      IREC = IRECT2
      N = NDIMT2
      IOFF = 1
    2 NBUF = MIN0(N,(RSRCLEN/RCW))
      CALL XCOPY (NBUF,T2(IOFF),1,BUFF,1)
      WRITE (MCCRES,REC=IREC) BUFF
      IREC = IREC + 1
      N = N - NBUF
      IOFF = IOFF + NBUF * RCW
      IF (N.GT.0) GOTO 2
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GETAMPT (T1,T2)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Read amplitudes from restart file
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      REAL*8 T1(*),T2(*)
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "complex.inc"
#include "symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 BUFF(RSRCLEN)
      DATA BUFF /RSRCLEN*0.D0/
      integer ioff,irec,irecl,n,nbuf,ireclen
C
C---------------Executable code--------------------------------------
C
      IRECL = IRECLEN (RSRCLEN,0,0)
      OPEN (MCCRES,FILE=RSFILN(1:LRSFIL),ACCESS='DIRECT',STATUS='OLD',
     &      RECL=IRECL)
C
      NRECT1 = NDIMT1*RCW / RSRCLEN + 1
      IRECT2 = IRECT1 + NRECT1 ! T2 amplitudes
      NRECT2 = NDIMT2*RCW / RSRCLEN + 1
C
      IREC = IRECT1
      N = NDIMT1
      IOFF = 1
    1 NBUF = MIN0(N,(RSRCLEN/RCW))
      READ (MCCRES,REC=IREC) BUFF
      CALL XCOPY (NBUF,BUFF,1,T1(IOFF),1)
      IREC = IREC + 1
      N = N - NBUF
      IOFF = IOFF + NBUF * RCW
      IF (N.GT.0) GOTO 1
C
      IREC = IRECT2
      N = NDIMT2
      IOFF = 1
    2 NBUF = MIN0(N,(RSRCLEN/RCW))
      READ (MCCRES,REC=IREC) BUFF
      CALL XCOPY (NBUF,BUFF,1,T2(IOFF),1)
      IREC = IREC + 1
      N = N - NBUF
      IOFF = IOFF + NBUF * RCW
      IF (N.GT.0) GOTO 2
C
      CLOSE (MCCRES,STATUS='KEEP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
