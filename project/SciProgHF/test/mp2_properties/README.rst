MP2 first order property calculations of various molecular properties
=====================================================================

 The properties caclulated are:
 - PVC: Parity Violation (Chirality) energies
 - EFG: the electric field gradient
 - DIPOLE: the dipole moment
 - EFN: the electric field at nuclei 
 - NQCC: nuclear quadrupole couplings constants

 Note: - The geometrical parameters are the same as described in the 
         Chem. Phys paper about MP2 first order properties with
         dihedral angle 45 degrees


   Calculations are done with Dirac-Coulomb and two-component (BSS+MFSSO,BSS,DC2BSS)
 Hamiltonians. DC2BSS values are almost identical with DC as this is only expectation 
 value property with noncanonical MO's from picture change transformed Fock-Dirac Hamiltonian.


! Property-MP2/Hamilt.   DC           BSS+MFSSO         BSS            DC2BSS
! -------------------------------------------------------------------------------------
! dip.mom./D        3.00817730        3.00850091      3.00849150     3.00816434
!
! EPV*10E-19/au   -3.727081365648   -3.444869996424  -5.065054307153  -3.727775795683
!
! el.field H z/au   0.0523372481    0.0523444418      0.0523448074    0.0523383502
! el.field O z/au  -0.1003834332   -0.1003726433     -0.1003726009   -0.1003830833
!
! EFG H qzz/au    0.5614548739      0.5614888924      0.5614889372   0.5614550883
! EFG O qzz/au    0.2358820836      0.2356679609      0.2356709452   0.2358858419
!
! NQCC Eta 1      0.19973050          0.19972304      0.19972231     0.19973092
! NQCC Eta 2      0.27426600          0.27425784      0.27425723     0.27426631

The timings are 1min40.000s for the BSS+MFSSO and 5min46.000s for the DC Hamiltonian (PC Linux).


