! File: energy.h
! Prefix -
!  ENER : energy contribution
!  GRAD : gradient contribution
! Suffix -
!  KE : kinetic energy
!  NA : nuclear attraction
!  EE : electron-electron repulsion
!  NN : nucleus-nucleus repulsion
!  FS : reorthonormalization ( trace(Fmat Smat) )
!
      REAL*8          ENERKE, ENERNA, ENEREE, ENERNN,                   &
     &                GRADKE(MXCOOR), GRADNA(MXCOOR), GRADEE(MXCOOR),   &
     &                GRADNN(MXCOOR), GRADFS(MXCOOR)
      COMMON /ENERGY/ ENERKE, ENERNA, ENEREE, ENERNN,                   &
     &                GRADKE, GRADNA, GRADEE, GRADNN, GRADFS
! --- end of energy.h ---
