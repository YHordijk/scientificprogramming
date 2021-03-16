=================== description of test =====================

The test calculates the ground state energy of SF_6 at
the HF level using the exact two-component (X2C) spinfree hamiltonian 
in the implementation interface to InteRest by M. Repisky.
the calculations are done using various X2C "flavors" listed in the following 
and described in more detail in the paper XXX:

--> molecular X2C approach (using the bare-nucleus h1 hamiltonian with "removal of spin-orbit terms at the end" 
    to derive the transformation matrix U) (large-molecular-x2c-spinfree.inp) [new X2C module implementation]

--> molecular X2C approach (using the bare-nucleus h1 hamiltonian with "removal of spin-orbit terms at the end" 
    to derive the transformation matrix U) (large-molecular-x2c-spinfree-end.inp) [old BSS/X2C implementation]

--> molecular X2C approach (using the bare-nucleus h1 hamiltonian with "removal of spin-orbit terms at the beginning" 
    to derive the transformation matrix U) (large-molecular-x2c-spinfree-beg.inp) [old BSS/X2C implementation]

--> atomic X2C approach using the transformation matrices U from atomic spinfree calculations 
                                           (large-atomic-x2c-spinfree.inp)

atomic X2C error wrt to full molecular X2C for SF6 (mH == milliHartree):
    -0.207 mH


    
    HF energy for molecular X2C - spinfree   (3-21G basis set) [new implementation]:
    -990.4068081192328 Hartree
    
    HF energy for molecular X2C - spinfree   (3-21G basis set) [old implementation] / end:
    -990.4068081183594 Hartree

    HF energy for molecular X2C - spinfree   (3-21G basis set) [old implementation] / beg:
    -990.4070480618358 Hartree
    
    HF energy for atomic    X2C + spinfree   (3-21G basis set):
    -990.4072548072289 Hartree

bug reports to dirac-users@googlegroups.com or stefan.knecht@phys.chem.ethz.ch
=================== =================== =====================
