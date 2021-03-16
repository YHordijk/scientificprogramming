 
 
!     Files to be used in RELCCSD
      INTEGER NCCFIL
      PARAMETER (NCCFIL=11) ! Number of word addressable files
!     Define the unit numbers

      INTEGER MSLVOUT,MRCONEE,MDINT,MDPROP,MCCRES,MCCDENS,ITAMPL,ITAMPD,                  &
     & IOAMPL,IOAMPS,ITAPT,ITAPL,ITAPADC,IUHEFF
      PARAMETER (MSLVOUT = 12) ! Output from slaves
      PARAMETER (MRCONEE = 13) ! Effective one-electron integrals
      PARAMETER (MDINT   = 14) ! Two-electron integrals
      PARAMETER (MDPROP  = 15) ! Property integrals
      PARAMETER (MCCRES  = 16) ! Restart information
      PARAMETER (MCCDENS = 17) ! Effective 1-electron density matrix
      PARAMETER (IUHEFF  = 18) ! Effective Hamiltonian (Fock space)
      PARAMETER (ITAMPL  = 26) ! T1 and T2 amplitudes (used in DIIS only)
      PARAMETER (ITAMPD  = 27) ! Error vector for T1 and T2 amplitudes
      PARAMETER (IOAMPL  = 28) ! Orbital gradients
      PARAMETER (IOAMPS  = 29) ! Sigma vector for orbital gradients
      PARAMETER (ITAPT   = 40) ! Sorted integral files
      PARAMETER (ITAPL   = 50) ! DZ intermediate for Lambda equations
      PARAMETER (ITAPADC = 60) ! Base for ADC-files
!     The standard input and output units may be redefined in the code
      COMMON /FILES/ IN,IW
      INTEGER IN,IW
!     Length of the records (integer words) on the restart file
      INTEGER RSRCLEN
      PARAMETER (RSRCLEN=1024)
      COMMON/RSFILE/ IRECIN,IRECST,IRECEN,IRECT1,NRECT1,IRECT2,NRECT2,                    &
     &               LRSFIL,LCCFIL
      INTEGER IRECIN,IRECST,IRECEN,IRECT1,NRECT1,IRECT2,NRECT2,LRSFIL,                    &
     & LCCFIL
      COMMON/PARNAM/ RSBASF,RSFILN,CCBASF,CCFILN
      CHARACTER*6  RSBASF      ! generic name for RESTART file
      CHARACTER*11 RSFILN      ! node-unique RESTART file name of length LRSFIL
      CHARACTER*11 CCBASF      ! generic name for RELCCSD output file
      CHARACTER*15 CCFILN      ! node-unique RELCCSD output file name
