!
!
!     Register the status of the independent parts of the code
!
!     Status given by integer :
!     0 - Never asked for
!     1 - On menu but not yet called
!     2 - Running
!     3 - Completed, restart information available
!     4 - Completed, no restart information
!
!     The indices of the different tasks are :
!      1 - Integral sort up to VOVO integrals
!      2 - Integral sort of VOVV integrals
!      3 - Integral sort of VVVV integrals
!      4 - Construction of Fock matrix
!      5 - MP2  energy and T1 + T2 amplitude calculation
!      6 - CCSD energy and T1 + T2 amplitude calculation
!      7 - CCSD(T) energy calculation
!      8 - fock-space sector 00 and T1 + T2 amplitudes
!      9 - fock-space sector 01 and T1 + T2 amplitudes
!     10 - fock-space sector 10 and T1 + T2 amplitudes
!     11 - fock-space sector 11 and T1 + T2 amplitudes
!     12 - fock-space sector 02 and T1 + T2 amplitudes
!     13 - fock-space sector 20 and T1 + T2 amplitudes
!     14 - 30 reserved for future use
!
      INTEGER :: NSTATS,IMSTAT
      PARAMETER (NSTATS=30)
      COMMON/STATUS/IMSTAT(NSTATS)
!
!     Gather results of the energy calculations
!     ETOT is the energy at the current highest level of theory
!     EFS is the energy of the lowest/selected eigenstate in the FSCC run (last sector)
!
      REAL*8 :: ECORE,ETOT,ESCF,EMP2,ECCSD,ECCSDIM,ET1,ET2,ET3,EFS(2)
      COMMON/ENERG/ECORE,ETOT,ESCF,EMP2,ECCSD,ECCSDIM,ET1,ET2,ET3,EFS
