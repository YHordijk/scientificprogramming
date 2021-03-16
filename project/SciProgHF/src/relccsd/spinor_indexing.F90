module spinor_indexing 

 implicit none
 private

 logical, save              :: module_initialized = .false.
 integer, allocatable, save :: pair_index(:,:,:,:)

 public init_index, pair_index

 contains

      subroutine init_index(fspc,nsp,irepsp,irepspi)

      integer :: nsp, irepsp(:,:), irepspi(:,:,:)
      logical :: fspc
      
      if  (module_initialized) call delete_index

      if (fspc) then
         allocate(pair_index(nsp,nsp,4,4))
         call fkindex(nsp,irepsp,irepspi,pair_index)
      else
         allocate(pair_index(nsp,nsp,4,1))
         call mkindex(nsp,irepsp,irepspi,pair_index)
      end if

      module_initialized = .true.

      end subroutine init_index

      subroutine delete_index()

      if (module_initialized) then
         deallocate(pair_index)
         module_initialized = .false.
      end if

      end

      SUBROUTINE MKINDEX (NSP,IREPSP,IREPSPI,INDEX)
 
      implicit none
!
!---------------Description--------------------------------------------
!
!     Set up index array for initial integral sort
!     Index(i,j,1) : type (oo,vo,ov,vv)
!     Index(i,j,2) : position of fm(i,j) in appropriate part of fock matrix
!     Index(i,j,3) : symmetry of ij pair
!     Index(i,j,4) : position of ij pair in batch of this symmetry
!
!---------------Routines called----------------------------------------
!
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Common Blocks--------------------------------------
!
#include "param.inc"
#include "symm.inc"
!
!---------------Calling variables--------------------------------------
!
      INTEGER NSP,IREPSP(NSP,4),IREPSPI(NSP,MXREP,2),INDEX(NSP,NSP,4,1)
!
!---------------Local variables--------------------------------------
!
      integer i,ii,ij,ijrep,imin,irep,j,jj,jrep
!
!---------------Executable code--------------------------------------
!
  INDEX = 0
! Type of index pair (oo,vo,ov or vv)
      DO J = 1, NSP
         IF (IREPSP(J,3).GT.0) THEN
            JJ = 0 
         ELSE 
            JJ = 2
         ENDIF
         DO I = 1, NSP
            IF (IREPSP(I,3).GT.0) THEN
               II = 1 
            ELSE 
               II = 2
            ENDIF
            INDEX(I,J,1,1) = II + JJ
            INDEX(I,J,2,1) = 0
            INDEX(I,J,3,1) = 0
            INDEX(I,J,4,1) = 0
         ENDDO
      ENDDO
! Position in Fock matrix: oo part
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NO(IREP)
            J = IREPSPI(JJ,IREP,1)
            DO II = 1, NO(IREP)
               I = IREPSPI(II,IREP,1)
               IJ = IJ + 1
               INDEX(I,J,2,1) = IJ
            ENDDO
         ENDDO
      ENDDO
! Position in Fock matrix: vo part
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NO(IREP)
            J = IREPSPI(JJ,IREP,1)
            DO II = 1, NV(IREP)
               I = IREPSPI(II,IREP,2)
               IJ = IJ + 1
               INDEX(I,J,2,1) = IJ
            ENDDO
         ENDDO
      ENDDO
! Position in Fock matrix: ov part
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NV(IREP)
            J = IREPSPI(JJ,IREP,2)
            DO II = 1, NO(IREP)
               I = IREPSPI(II,IREP,1)
               IJ = IJ + 1
               INDEX(I,J,2,1) = IJ
            ENDDO
         ENDDO
      ENDDO
! Position in Fock matrix: vv part
      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NV(IREP)
            J = IREPSPI(JJ,IREP,2)
            DO II = 1, NV(IREP)
               I = IREPSPI(II,IREP,2)
               IJ = IJ + 1
               INDEX(I,J,2,1) = IJ
            ENDDO
         ENDDO
      ENDDO
! Position in |i,j> or <i,j| triangular lists: oo part
      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 10 IREP = JREP, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 10
               DO JJ = 1, NO(JREP)
                  J = IREPSPI(JJ,JREP,1)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = JJ + 1
                  DO II = IMIN, NO(IREP)
                     I = IREPSPI(II,IREP,1)
                     IJ = IJ + 1
                     INDEX(I,J,3,1) = IJREP
                     INDEX(I,J,4,1) = IJ
                  ENDDO
               ENDDO
 10         CONTINUE
         ENDDO
      ENDDO
! Position in |i,j> or <i,j| triangular lists: vo part
      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 20 IREP = 1, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 20
               DO JJ = 1, NO(JREP)
                  J = IREPSPI(JJ,JREP,1)
                  DO II = 1, NV(IREP)
                     I = IREPSPI(II,IREP,2)
                     IJ = IJ + 1
                     INDEX(I,J,3,1) = IJREP
                     INDEX(I,J,4,1) = IJ
                  ENDDO
               ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
! Position in |i,j> or <i,j| triangular lists: vv part
      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 40 IREP = JREP, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 40
               DO JJ = 1, NV(JREP)
                  J = IREPSPI(JJ,JREP,2)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = JJ + 1
                  DO II = IMIN, NV(IREP)
                     I = IREPSPI(II,IREP,2)
                     IJ = IJ + 1
                     INDEX(I,J,3,1) = IJREP
                     INDEX(I,J,4,1) = IJ
                  ENDDO
               ENDDO
 40         CONTINUE
         ENDDO
      ENDDO
 
      RETURN 
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FKINDEX (NSP,IREPSP,IREPSPI,INDEX)
 
 
!---------------Description--------------------------------------------
!
!     Set up index array for initial integral sort
!     Index(i,j,1) : type (oo,vo,ov,vv)
!     Index(i,j,2) : position of fm(i,j) in appropriate part of fock matrix
!     Index(i,j,3) : symmetry of ij pair
!     Index(i,j,4) : position of ij pair in batch of this symmetry
!
!---------------Routines called----------------------------------------
!
!
!---------------Last modified------------------------------------------
!
!     Author : Luuk Visscher
!
!---------------Common Blocks--------------------------------------
!
#include "param.inc"
#include "symm.inc"
!
!---------------Calling variables--------------------------------------
!
      INTEGER NSP,IREPSP(NSP,4),IREPSPI(NSP,MXREP,2),INDEX(NSP,NSP,4,4)
!
!---------------Local variables--------------------------------------
!
      integer i,icase,ii,ij,ijcase,ijrep,imin,irep,j,jcase,jj,jrep
!
!---------------Executable code--------------------------------------
!
      DO IJCASE = 1, 4
         JCASE = (IJCASE-1)/2
         ICASE = MOD(IJCASE-1,2)
         DO J = 1, NSP
            IF (IREPSP(J,3+JCASE).GT.0) THEN
               JJ = 0 
            ELSE 
               JJ = 2
            ENDIF
            DO I = 1, NSP
               IF (IREPSP(I,3+ICASE).GT.0) THEN
                  II = 1 
               ELSE 
                  II = 2
               ENDIF
               INDEX(I,J,1,IJCASE) = II + JJ
               INDEX(I,J,2,IJCASE) = 0
               INDEX(I,J,3,IJCASE) = 0
               INDEX(I,J,4,IJCASE) = 0
            ENDDO
         ENDDO
      ENDDO

      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NO(IREP)
            J = IREPSPI(JJ,IREP,1)
            DO II = 1, NO(IREP)
               I = IREPSPI(II,IREP,1)
               IJ = IJ + 1
               INDEX(I,J,2,1) = IJ
            ENDDO
         ENDDO
      ENDDO

      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NO(IREP)
            J = IREPSPI(JJ,IREP,1)
            DO II = 1, NV(IREP)
               I = IREPSPI(II,IREP,2)
               IJ = IJ + 1
               INDEX(I,J,2,2) = IJ
            ENDDO
         ENDDO
      ENDDO

      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NV(IREP)
            J = IREPSPI(JJ,IREP,2)
            DO II = 1, NO(IREP)
               I = IREPSPI(II,IREP,1)
               IJ = IJ + 1
               INDEX(I,J,2,3) = IJ
            ENDDO
         ENDDO
      ENDDO

      IJ = 0
      DO IREP = 1, NREP
         DO JJ = 1, NV(IREP)
            J = IREPSPI(JJ,IREP,2)
            DO II = 1, NV(IREP)
               I = IREPSPI(II,IREP,2)
               IJ = IJ + 1
               INDEX(I,J,2,4) = IJ
            ENDDO
         ENDDO
      ENDDO

      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 10 IREP = JREP, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 10
               DO JJ = 1, NO(JREP)
                  J = IREPSPI(JJ,JREP,1)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = JJ + 1
                  DO II = IMIN, NO(IREP)
                     I = IREPSPI(II,IREP,1)
                     IJ = IJ + 1
                     INDEX(I,J,3,1) = IJREP
                     INDEX(I,J,4,1) = IJ
                  ENDDO
               ENDDO
 10         CONTINUE
         ENDDO
      ENDDO

      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 20 IREP = 1, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 20
               DO JJ = 1, NO(JREP)
                  J = IREPSPI(JJ,JREP,1)
                  DO II = 1, NV(IREP)
                     I = IREPSPI(II,IREP,2)
                     IJ = IJ + 1
                     INDEX(I,J,3,2) = IJREP
                     INDEX(I,J,4,2) = IJ
                  ENDDO
               ENDDO
 20         CONTINUE
         ENDDO
      ENDDO

      DO IJREP = 1, NREP
         IJ = 0
         DO JREP = 1, NREP
            DO 40 IREP = JREP, NREP
               IF (MULTB(IREP,JREP,1).NE.IJREP) GOTO 40
               DO JJ = 1, NV(JREP)
                  J = IREPSPI(JJ,JREP,2)
                  IMIN = 1
                  IF (IREP.EQ.JREP) IMIN = JJ + 1
                  DO II = IMIN, NV(IREP)
                     I = IREPSPI(II,IREP,2)
                     IJ = IJ + 1
                     INDEX(I,J,3,4) = IJREP
                     INDEX(I,J,4,4) = IJ
                  ENDDO
               ENDDO
 40         CONTINUE
         ENDDO
      ENDDO

      RETURN 
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

end module 
