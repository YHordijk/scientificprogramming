!     LMOORD: Reorder orbitals?
!     MXREORD: maximum number of reorderings
!     IREORD, IMOORD: orbital reordering

      PARAMETER ( MXREORD = 250 )
      LOGICAL LMOORD
      COMMON /CBIREO/ IMOORD(MXREORD,2),IREORD(2)
      COMMON /CBLREO/ LMOORD
!
!     LMOORF: Post-DHF Reorder orbitals?
!     IREORF, IMOORF: orbital reordering

      LOGICAL LMOORF
      COMMON /CBI2RE/ IMOORF(MXREORD,2),IREORF(2)
      COMMON /CBL2RE/ LMOORF
