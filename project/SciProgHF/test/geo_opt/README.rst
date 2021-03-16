Geometry optimization with DIRAC
================================


 This test performs geometry optimizations for different molecules.



 The test also illustrates calculation of properties at the final converged geometry:

 HeH+ molecule, comparison of BSS, BSS+MFSSO (.AMFICH=1)
   and DCrkb(lin.response with .SKIPEP) values

        Hamiltonian        DC_rkb      BSS+MFSSO      BSS
  =====================================================================================
               r(HeH+/A)  0.772857     0.772869     0.772869
    H isotr.shielding/au  11.8773      11.8779      11.8779
   He isotr.shielding/au  55.7332      55.7334      55.7334
   DSO isotr.contrib./au  3.9528       3.9526       3.9526
     Final spin-spin./au  -615.4195    -615.4744    -615.4744
         Total dip.mom/D  1.12682502   1.12682934   1.12682934
