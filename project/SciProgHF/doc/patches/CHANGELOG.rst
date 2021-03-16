:orphan:


20.0 (planned for 2021-01-31)
-----------------------------

New features:

- New parallel coupled cluster implementation ExaCorr. Contributors: J.V. Pototschnig,A. Papadopoulos,
   D. I. Lyakh, M. Repisky,L.Halbert,A.S.P. Gomes,H.J.Aa. Jensen, and L. Visscher.
   To be published.
- Improved start potential with far fewer functions and better accuracy. Contributor: Lucas Visscher.
   S. Lehtola,  L. Visscher and  E. Engel.
   Efficient implementation of the superposition of atomic potentials initial guess for electronic structure calculations in Gaussian basis sets.
   J. Chem. Phys. 152, 144105 (2020); http://dx.doi.org/10.1063/5.0004046 doi: 10.1063/5.0004046
- Calculation of Rotational Molecular g-tensors. Activated with ".ROTG". Contributor: I. Agustín Aucar.
   I. A. Aucar, S. S. Gómez, C. G. Giribet and M. C. Ruiz de Azúa.
   Theoretical study of the relativistic molecular rotational g-tensor.
   J. Chem. Phys. 141, 194103 (2014); http://dx.doi.org/10.1063/1.4901422 doi: 10.1063/1.4901422
- Improved reading of basis sets to be able to read all EMSL basis sets of type DALTON. Contributor: Hans Joergen Aa. Jensen.
- "make latex" for making a latex version of the documenation in "build/latex", used e.g. pdflatex to make a pdf version of the manual. Contributor: Hans Joergen Aa. Jensen

Revised features:

- Gauge origin (.GAUGEO and .GO ANG) now can only be set under \*\*HAMILTONIAN.
- Phase and dipole origins (.PHASEO and .DIPORG, respectively) now can only be set under \*\*HAMILTONIAN.
- EOM-CC excitation energies collected for all symmetries and summarized at the end  

New defaults:

- One-electron operator ANGMOM's origin was moved from the gauge-origin to the molecular center of mass.
- The CODATA2018 set of physical constants is now used as default (taken from http://physics.nist.gov/constants). Before this version, CODATA1998 was default.
- Non-relativistic and relativistic effective core-potential calculations (".NONREL" and ".ECP") now have nuclear point nucleus as default.
  For non-relativistic calculations this makes it easier to compare to other codes, and for RECP this is the intended behavior which by mistake was not activated.

-----------------------------------------------------------------------------------------------------

19.0 (2019-12-12)
-----------------

New features:

- Calculation of Nuclear Spin-Rotation tensors. Activated with ".SPINRO". Contributors: Trond Saue and I. Agustín Aucar.
   I. A. Aucar, S. S. Gómez, M. C. Ruiz de Azúa, and C. G. Giribet.
   Theoretical study of the nuclear spin-molecular rotation coupling for relativistic electrons and non-relativistic nuclei.
   J. Chem. Phys. 136, 204119 (2012); http://dx.doi.org/10.1063/1.4721627 doi: 10.1063/1.4721627
- Interface with the Quantum Computing package OpenFermion. https://github.com/bsenjean/Openfermion-Dirac. Contributors: Bruno Senjean and Luuk Visscher.
   T. E. O'Brien, B. Senjean, R. Sagastizabal, X. Bonet-Monroig, A. Dutkiewicz, F. Buda, L. DiCarlo and L. Visscher
   Calculating energy derivatives for quantum chemistry on a quantum computer.
   Accepted to npj: Quantum Information. https://arxiv.org/abs/1905.03742
- Expectation value of Nuclear Magnetic-Quadrupole-Moment interaction constant using ".OPERATOR" and ".NMQM" keyword in KRCI. Contributor: Malaya K. Nayak
   T. Fleig, M. K. Nayak and M. G. Kozlov.
   TaN, a molecular system for probing P,T-violating hadron physics.
   Phys. Rev. A 93, 012505, 2016; http://dx.doi.org/10.1103/PhysRevA.93.012505 doi: 10.1103/PhysRevA.93.012505

Revised features:

- Maximum number of basis functions increased from 30k to 40k.

Error corrections:

- Fixed use of point charges (do not try to obtain e.g. mass of a point charge).
- Fixed spelling error in keyword .GAUNTSCALE (was .GAUTSCALE).
- Fixed T1EQN1 function call in src/relccsd/ccgrad.F (thanks to Martin Siegert).
- Fixed wrong memory allocation for spin-free or non-relativistic KRMCSCF.
- Fixed runtime error for XPROPMAT in src/krmc/krmccan.F property modules.
- Fixed runtime error for calling XPRIND in src/prp/pamprp.F property modules.

Basis set corrections:

- Fixed Be basis set in cc-pV5Z.
- Fixed Al basis set in aug-cc-pV(D+d)Z, aug-cc-pV(T+d)Z and aug-cc-pV(Q+d)Z.
- Fixed Si basis set in aug-cc-pV(D+d)Z and aug-cc-pV(Q+d)Z.

Performance improvements:

- Reduced static memory requirement from common blocks.

Documentation:

- Restore local build of Sphinx documentation.

-----------------------------------------------------------------------------------------------------

18.0 (2018-12-13)
-----------------

New features:

- The visualization module has been extended to generate various radial distributions using the keyword .RADIAL. Contributor: Trond Saue.

Error corrections:

- Fix f correlating exponents for 3d in dyall.cv2z and dyall.ae2z basis sets. Contributor: Ken G. Dyall.
- Fix f polarizing exponents for 3d in dyall.cv3z and dyall.ae3z basis sets. Contributor: Ken G. Dyall.
- Bugfix in KRCI module (subroutine STRTYP_GAS_REL), reported by Shigeyoshi Yamamoto.
- Bugfix in polarizable embedding (the printed electronic part of the energy was not corrected with PE contributions).
- Fixed the formatting of the carbon potential in the file basis_ecp/ECPCE2SO.

Performance improvements:

- New much more efficient quaternion diagonalization routine using openMP for performance. Contributor: Hans Joergen Aa. Jensen.
- New much more efficient complex diagonalization routine using openMP for performance. Contributor: Hans Joergen Aa. Jensen.
- Use MPI_REDUCE for adding Fock matrix contributions from MPI worker nodes on MPI master node. Contributor: Hans Joergen Aa. Jensen.
- Save memory in direct Fock matrix constructions, now >13000 basis functions can be run in 32GB RAM. Contributor: Hans Joergen Aa. Jensen.
- (With these performance improvements, a 4c HF test calculation with 13000 basis functions, 3300 postive energy orbitals (6600 in total) has been run on Titan at ORNL with 800 nodes, each with 8 openMP threads, at 1h59m-2h06m per HF iteration.)

Framework changes:

- Allow to compile with OpenMP support.
- Update in CMake framework: pass compilers as cmake variables, not as environment variables.

-----------------------------------------------------------------------------------------------------

17.0 (2017-12-12)
-----------------

- DIRAC17 release.

New features:

- Kramers restricted Polarization Propagator in the ADC framework for electronic excitations, activated with ".POLPRP".
   M. Pernpointner. The relativistic polarization propagator for the calculation of electronic excitations in heavy systems.
   J. Chem. Phys. 140, 084108 (2014); http://dx.doi.org/10.1063/1.4865964

   M. Pernpointner, L. Visscher and A. B. Trofimov. 
   Four-component Polarization Propagator Calculations of Electron Excitations: Spectroscopic Implications of Spin-Orbit Coupling   Effects. J. Chem. Theory Comput. (2017) submitted.

- Polarizable embedding using PElib (https://gitlab.com/pe-software/pelib-public). Activated with ".PEQM", additional options under \*PEQM.

 Reference:
 E. D. Hedegård, R. Bast, J. Kongsted, J. M. H. Olsen, and H. J. Aa. Jensen:  Relativistic Polarizable Embedding., J. Chem. Theory Comput. 13, 2870-2880 (2017); http://doi.org/10.1021/acs.jctc.7b00162

- New and numerically stable procedure for elimination/freezing of orbitals at SCF level. Author: T. Saue.
- New ".MVOFAC" option in \*KRMC input section. Author: H. J. Aa. Jensen.
- New easier options for point charges in the .mol file: "LARGE POINTCHARGE" or "LARGE NOBASIS" (the two choices are equivalent).
- Added the RPF-4Z and aug-RPF-4Z basis sets for f-elements to the already existing files with sets for s, p and d elements. Deleted the aug-RPF-3Z set as that was not an official set.
- Write out effective Hamiltonian in Fock space coupled cluster to a file for post processing. Can be used with external code of Andrei Zaitsevskii (St. Petersburg).
- Provided memory counter for RelCC calculations, what is suitable for memory consuming large scale Coupled Cluster calculations. See the web-documentation.
- Updated basis_dalton/ with basis set updates in the Dalton distribution:
- fix of errors in Ahlrichs-pVDZ (several diffuse exponents were a factor 10 too big)
- fix of errors for 2. row atoms in aug-cc-pCV5Z
- added many atoms to aug-cc-PVTZ_J
- added many Frank Jensen "pc" type basis sets
- added Turbomole "def2" type basis sets

Error corrections:

- Fixed errors for quaternion symmetries in 2-electron MO integrals used in CI calculations with GASCIP.
  It is now possible to do CI calculations with GASCIP for C1 symmetry (i.e. no symmetry).
- Fixed the p exponents for Na in the dyall 4z basis sets to match the archive.
  The changes are small so should not significantly affect results.
- Fixed compilation of XCFun on Mac OS X High Sierra.
- Fixed error for parallel complex CI or MCSCF with GASCIP

New defaults:

- .SKIPEP is now default for KR-MCSCF, new keyword .WITHEP to include e-p rotations

Performance improvements:

- restored integral screening behavior from DIRAC15

-----------------------------------------------------------------------------------------------------

16.0 (2016-12-12)
-----------------

- DIRAC16 release.

New features:

- RELCCSD expectation values aka first-order unrelaxed analytic energy derivative.
  For more information, see J. Chem. Phys. 145 (2016) 184107 as well as test/cc_gradient for an example.
- Calculation of approximate dipole intensities together with the excitation energies in Fock space CC.
  See test/fscc_intensities for an example.
- Improved start potential for SCF: sum of atomic LDA potentials, generated by GRASP.

New defaults:

- Negative denominators (e.g. appearing in core ionized systems) accepted in RELCCSD
- AOFOCK is now default if at least 25 MPI nodes (parallelizes better than SOFOCK). .AOFOCK documented.

Error corrections and updates in isotope properties for the following atoms:
 -    Br isotope 2: quadrupole moment  .2620 ->  .2615
 -    Ag isotope 2: magnetic moment    .130563 -> -.130691 (note sign change)
 -    In isotope 2: quadrupole moment   .790 ->  .799
 -    Nd magnetic moments of isotopes 4 and 5 were interchanged: -0.065 -> -1.065 and -1.065 -> -0.065
 -    Gd: quadrupole meoments of isotopes 4 and 5 updated: 1.36 -> 1.35 and 1.30 -> 1.27
 -    Ho isotope 1: quadrupole moment updated 3.49 -> 3.58
 -    Lu isotope 2: quadrupole moment updtaed 4.92 -> 4.97
 -    Hf isotope 1: mass was real*4, not real*8, thus 7 digits instead of 179.9465457D0 (i.e. approx 179.9465)
 -    Ta isotope 1: quadrupole moment added 0.00 -> 3,17
 -    Tl isotope 1: nuclear moment 1.63831461D0 -> 1.63821461D0 (typo, error 1.d-4)
 -    Pb isotope 3: nuclear moment 0.582583D0 -> 0.592583D0 (typo, error 1.d-2)
 -    Po isotope 1: nuclear moment added: 0.000 -> 0.688

Framework improvements:

- Faster Dirac compilation due to implemented OBJECT libraries (requires CMake of version at least 2.8.10).
- Run benchmark tests and only benchmark tests with "ctest -C benchmarks -L benchmark".
- Change any tab character in coordinate specifications in .mol file to a blank.
  (tab characters have been reported to cause errors in reading coordinates).
- Added OPENBLAS in search list for blas and lapack libraries.
- Made keyword .X2Cmmf case insensitive.
- Cache compiler flags.

Bug fixes:

- Bugfix: MCSCF and GASCIP was not correct for odd number of electrons.
- Bugfix: now .CASSCF works again for KRMSCF input.
- Bugfix: KRMSCF was not converging is some cases (problem reported by Miro Ilias for Os atom).
- Bugfix: "all" keyword for ".GAS SHELLS" was not recognized (problem reported by Miro Ilias).
- Bugfix for numerical gradient when automatic symmetry detection: do no move or rotate molecule
  for finding highest possible symmetry because then the numerical gradient will not be correct
  (problem reported by Miro Ilias).
- Bugfix for spinfree X2Cmmf and linear symmetry.
- Bugfix: restored "make install" target.
- Bugfix: Test reladc_fano does not report failure if Stieltjes is deactivated.

-----------------------------------------------------------------------------------------------------

15.0 (2015-12-16)
-----------------

- FanoADC-Stieltjes: Calculation of decay widths of electronic decay processes.
- Better convergence for open-shell calculations, using .OPENFAC 0.5 as default.
  Final orbital energies are recalculated with .OPENFAC 1.0, for IP interpretation.
- Performance improvement for MP2 NO: only transform (G O|G O) integrals instead of (G G|G G).
- Geometry optimization with xyz input.
- Performance improvements for determinant generation in GASCIP.
- Use Kramers conjugation on CI vectors (cuts time for CI in half for ESR doublets).
- ANO-RCC basis: Fixed Carbon basis set (wrong contraction coefficients, see [MOLCAS ANO-RCC](http://www.molcas.org/ANO/).
- ANO-RCC basis: Modified the 3 Th h-functions by replacing them with the 3 Ac h-functions to Th.
  (A mistake was made in the original work when the 3 Th h-functions were made,
  and this has never been redone. They are too diffuse, exponents
  (0.3140887600, 0.1256355100, 0.0502542000) for Th, Z=90, compared to
  (0.7947153600, 0.3149038200, 0.1259615200) for Ac, Z=89, and
  (0.8411791300, 0.3310795400, 0.1324318200) for Pa, Z=91.
  We have selected to just replace the 3 Th h-functions with those from the Ac basis set,
  because the Ac g-functions are quite close tot he Th g-functions, closer than Ac g-functions,
  and therefore differences in results compared to optimized Th h-functions should be minimal.)
  Thanks to Kirk Peterson for pointing out the Th problem on http://daltonforum.org.
- Fixed reading of ANO-RCC and ANO-DK3 basis sets from the included basis set library.
- Update PCMSolver module to its [1.0.4 version](https://github.com/PCMSolver/pcmsolver/releases).
- Do not write Hartree-Fock in output if it is a DFT calculation, write DFT.
- Fix a bunch of problems where the number of arguments of the call does not match the subroutine definition (submitted by Martin Siegert).
- Fixed plotting of the electrostatic potential (electronic part was factor 2 too large).
- Configuration framework uses [Autocmake](http://autocmake.org).

-----------------------------------------------------------------------------------------------------

14.1 (2015-07-20)
-----------------

- Added OpenBLAS to default blas search list.
- Fixed misleading error message for .OPEN SHELL input.
- Workaround for CMake warning about nonexistent run_pamadm target.
- Intel compilers xHost flag automatic detection only upon request (-D ENABLE_XHOST_FLAG_DETECTION=ON).
- Fixed minor issues to allow compilation using IBM XL Fortran compiler (thanks to Bob Walkup, IBM).
- Update PCMSolver module to its [1.0.3 version](https://github.com/PCMSolver/pcmsolver/releases).


14.0 (2014-12-12)
-----------------

- DIRAC14 release, see doc/release/release-statement.txt
