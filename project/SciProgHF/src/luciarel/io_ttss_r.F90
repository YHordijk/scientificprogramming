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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE civcinf_mc(IBATS,LBATS,I1BATS,LEBATS,NBAT_MS2, &
                            NMS2VAL,NBLOCK)
!*********************************************************************
!
!     write TTSS-block information to file for KRMC - LUCIAREL
!     communication
!
!     written by S. Knecht - Jun 2009
!
!     last revision:
!
!**********************************************************************
      implicit none
!
      integer, intent(in) :: NMS2VAL, NBLOCK
      integer             :: LUX
      real(8), intent(in) :: IBATS(8,NBLOCK),LBATS(NBLOCK)
      real(8), intent(in) :: I1BATS(NBLOCK),LEBATS(NBLOCK)
      real(8), intent(in) :: NBAT_MS2(NMS2VAL)
!
!     write block information to file INFO.TTSS_BLOCK
      LUX = 95
      OPEN(UNIT=LUX,FILE='INFO.TTSS_BLOCK',STATUS='UNKNOWN', &
           FORM='UNFORMATTED')
      REWIND LUX
      WRITE(LUX) NBLOCK
      WRITE(LUX) NMS2VAL 
      WRITE(LUX) NBAT_MS2 
      WRITE(LUX) LEBATS 
      WRITE(LUX) LBATS 
      WRITE(LUX) I1BATS 
      WRITE(LUX) IBATS
      CLOSE(LUX)
!
      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SAVEVCDC(VECT,LU,IRILP,IPRNT)
!*********************************************************************
!
!     use vector partitioning information to save a vector
!     to LUCIAREL file.
!     driver routine 
!
!     written by S. Knecht - Jun 2009
!
!     last revision:
!
!**********************************************************************
use memory_allocator
      implicit none
!
      real(8), intent(inout) :: VECT(*)
      integer, intent(in)    :: LU, IRILP, IPRNT
      integer                :: NTESTL, NTEST
      integer                :: LUX, NMK2VAL, NBLOCK, LEN
      logical                :: EX
      integer, allocatable   :: NBAT_MS2_MC(:), LEBATS_MC(:)
      integer, allocatable   :: LBATS_MC(:), I1BATS_MC(:), IBATS_MC(:)
      character (LEN=5)      :: id
!
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
!
!     read block information from file INFO.TTSS_BLOCK
      LUX = 95
      INQUIRE (FILE="INFO.TTSS_BLOCK",EXIST=EX)
      IF( EX )THEN
        OPEN(UNIT=LUX,FILE='INFO.TTSS_BLOCK',STATUS='OLD', &
             FORM='UNFORMATTED')
        REWIND LUX
        READ(LUX) NBLOCK
        READ(LUX) NMK2VAL
!
!       allocate
        call alloc(NBAT_MS2_MC, NMK2VAL,  id="NBA  ")
        call alloc(LEBATS_MC  , NBLOCK,   id="LEB  ")
        call alloc(LBATS_MC   , NBLOCK,   id="LBA  ")
        call alloc(I1BATS_MC  , NBLOCK,   id="I1B  ")
        call alloc(IBATS_MC   , 8*NBLOCK, id="IBALL")
      ELSE
        write(6,'(/2A)') ' ** Error in SAVEVCDC: TTSS-block file ', &
                         ' INFO.TTSS_BLOCK not present **'
        CALL QUIT(' ** Error in SAVEVCDC: Block file not present **')
      END IF
      CALL read_ttss_blocks(NMK2VAL,NBAT_MS2_MC,LEBATS_MC, &
                            LBATS_MC,I1BATS_MC,IBATS_MC,   &
                            NBLOCK,LUX)
      CLOSE(LUX)
!
      CALL SAVEVCDC_RUN(VECT,LU,IRILP,NTEST,NBLOCK,NMK2VAL,& 
                        NBAT_MS2_MC,LEBATS_MC,LBATS_MC,I1BATS_MC,&
                        IBATS_MC)
!
!     deallocate
      call dealloc(NBAT_MS2_MC)
      call dealloc(LEBATS_MC)
      call dealloc(LBATS_MC)
      call dealloc(I1BATS_MC)
      call dealloc(IBATS_MC)
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SAVEVCDC_RUN(VECT,LU,IRILP,NTEST,NBLOCK,NMK2VAL,&
                              NBAT_MS2_MC,LEBATS_MC,LBATS_MC,I1BATS_MC,&
                              IBATS_MC)
!*********************************************************************
!
!     use vector partitioning information to save a vector
!     to LUCIAREL file.
!
!     written by S. Knecht - Jun 2009
!
!     last revision:
!
!**********************************************************************
      implicit none
!
      real(8), intent(inout) :: VECT(*)
      integer, intent(in)    :: LU, IRILP, NTEST
      integer                :: IMK2, JBATS, JBATABS
      integer                :: JOFF, LENGTH, IRI, LS, KBATSEND, IOFF
      integer                :: ISBLK, LEN_V
      integer, intent(in)    :: NMK2VAL, NBLOCK
      integer                :: NBAT_MS2_MC(NMK2VAL), LEBATS_MC(NBLOCK)
      integer                :: LBATS_MC(NBLOCK), I1BATS_MC(NBLOCK)
      integer                :: IBATS_MC(8,NBLOCK)

      if (NTEST.ge.100) then
        write(6,*) ' '
        write(6,*) ' +++++++++++++++++++++++++++++++ '
        write(6,*) '  SAVEVCDC called.               '
        write(6,*) '   (Save complex vec. to file.)  '
        write(6,*) ' +++++++++++++++++++++++++++++++ '
        write(6,*) '  File structure information:'
        write(6,*) '  Number of MK values :',NMK2VAL
        JBATABS = 0
        do IMK2 = 1,NMK2VAL,1
          write(6,*) '  MK value no.: ',IMK2
          do JBATS = 1,NBAT_MS2_MC(IMK2)
            JBATABS = JBATABS + 1
            write(6,*) ' Absolute batch no.: ',JBATABS
            write(6,*) ' No. of elements in batch :',LEBATS_MC(JBATABS)
          end do
        end do 
      end if
!
      JOFF = 1
      LENGTH = 0
!
      DO IRI = 1,IRILP,1
! Write end of vector mark between real and imaginary part.
        if (IRI.eq.2) call itods(-1,1,-1,LU)
!
        JBATABS = 0
        if (NTEST.ge.100) write(6,*) '  REAL/IMAG = ',IRI
        DO IMK2 = 1, NMK2VAL
!
! Loop over batches of sigma for this MK2 projection value
          do JBATS = 1,NBAT_MS2_MC(IMK2)
            JBATABS = JBATABS + 1
!
            if (NTEST.ge.100) then
               WRITE(6,*)
               WRITE(6,*) ' Start of writing of batch ', JBATABS
               write(6,'(2X,A,I3)') ' for MK2 value number ',IMK2
               WRITE(6,*) 
            end if
!
            LS = LEBATS_MC(JBATABS)
            LENGTH = LENGTH + LS
            if (NTEST.ge.100.and.LS.gt.0) then
              WRITE(6,*) ' **********************************'
              WRITE(6,*) ' * Array containing current batch *'
              WRITE(6,*) ' **********************************'
              CALL WRTMAT(VECT(JOFF),1,LS,1,LS)
            end if
            KBATSEND = I1BATS_MC(JBATABS) + LBATS_MC(JBATABS)
            do ISBLK = I1BATS_MC(JBATABS),KBATSEND-1,1
              IOFF = IBATS_MC(6,ISBLK)
              LEN_V = IBATS_MC(8,ISBLK)
              if (NTEST.ge.100) then
                write(6,*) '    Saving block no. ',ISBLK
                write(6,*) '     Block offset:   ',IOFF
                write(6,*) '     Block length:   ',LEN_V
              end if
              call itods(LEN_V,1,-1,LU)
              call todsc(VECT(JOFF+IOFF-1),LEN_V,-1,LU)
            end do
            JOFF = JOFF + LS
!
          end do
!         ^ End loop over S batches of given MK2 value
        end do
!       ^ End of loop over IMK2
!
      end do
!     ^ End of loop over IRI
!     write end of vector mark.
      CALL ITODS(-1,1,-1,LU)
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE READ_TTSS_BLOCKS(NMK2VAL,NBAT_MS2_MC,LEBATS_MC, &
                                  LBATS_MC,I1BATS_MC,IBATS_MC,   &
                                  NBLOCK,LUX)

!*********************************************************************
!
!     read TTSS-block information
!
!     written by S. Knecht - Jun 2009
!
!**********************************************************************
      implicit none
      integer, intent(in)  :: NMK2VAL, NBLOCK, LUX
      integer              :: NBLOCK_DUM, NMK2VAL_DUM
      integer              :: NBAT_MS2_MC(NMK2VAL)
      integer              :: LEBATS_MC(NBLOCK)
      integer              :: LBATS_MC(NBLOCK),I1BATS_MC(NBLOCK)
      integer              :: IBATS_MC(8,NBLOCK)
      REWIND LUX
      READ (LUX) NBLOCK_DUM
      READ (LUX) NMK2VAL_DUM
      READ (LUX) NBAT_MS2_MC 
      READ (LUX) LEBATS_MC 
      READ (LUX) LBATS_MC 
      READ (LUX) I1BATS_MC 
      READ (LUX) IBATS_MC
      END
