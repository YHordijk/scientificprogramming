

.. raw:: html

  <a class="btn btn-info"
     href="http://diracprogram.org"
     data-mode="1"
     style="margin-bottom: 50px"
     target="_blank">Project website: http://diracprogram.org</a>


Setting up DIRAC
================

* **Installation:**
  :doc:`Requirements <installation/requirements>` |
  :doc:`Linux, Unix, Mac <installation/general>` |
  :doc:`Windows <installation/windows>` |
  :doc:`Expert options <installation/expert>` |
  :doc:`System administrators <installation/sysadmins>`

* **Math libraries:**
  :doc:`Detection and linking <installation/math>` |
  :doc:`MKL environment variables <installation/mkl>`

* **MPI:**
  :doc:`Forwarding environment variables <installation/mpi>`

* **64-bit integer support:**
  :doc:`Do I need it? <installation/int64/faq>` |
  :doc:`Math libraries <installation/int64/math>` |
  :doc:`64-bit OpenMPI <installation/int64/mpi64>` |
  :doc:`Troubleshooting <installation/int64/troubleshooting>`

* **Testing:**
  :doc:`Running the test set <installation/testing>`

* **Example installations and run scripts:**
  :doc:`HPC cluster (Strasbourg) <installation/examples/hpc-unistra>`

* **pam script:**
  :doc:`Setting the scratch directory <pam/scratch>` |
  :doc:`Alternative MPI launcher and passing arguments for it <pam/mpi_command>` |
  :doc:`String replacement <pam/replace>`


Updates and patches
===================

 :doc:`ChangeLog <patches/CHANGELOG>` |
 :doc:`Known problems <patches/KNOWN-PROBLEMS>`

.. _FirstSteps:

First steps in DIRAC territory
==============================

* **Getting started:**
  :doc:`First calculation <tutorials/getting_started>`

* **Molecule input formats:**
  :doc:`mol format <molecule_and_basis/molecule_using_mol>` |
  :doc:`xyz format <molecule_and_basis/molecule_using_xyz>` |
  :doc:`ecp input <molecule_and_basis/molecule_with_ecp>`   |
  :doc:`Troubleshooting <molecule_and_basis/troubleshooting>`

* **Basis sets:**
  :doc:`General information <molecule_and_basis/basis>` |
  :doc:`Howto uncontract basis sets <molecule_and_basis/howto_uncontract>`

* **Troubleshooting:**
  :doc:`Known problems <getting_help/known_problems>` |
  :doc:`Memory problems <getting_help/memory>`


Reference manual
================

.. toctree::
   :maxdepth: 1

   manual/index.rst

Tutorials and walkthrus
=======================

* **Basis sets:**
  :doc:`Basis sets for relativistic calculations <tutorials/relbas>` |
  :doc:`Augmenting basis sets <tutorials/augmentation>`

* **SCF start guess :**
  :doc:`Atomic start <tutorials/start_guess/atomic_start>` |
  :doc:`Extended Hückel start <tutorials/start_guess/atomic_huckel>` |

* **Restarting and multi-step jobs:**
  :doc:`SCF <tutorials/restart/scf>` |
  :doc:`X2C <tutorials/restart/x2c>` |
  :doc:`Coupled Cluster restart <tutorials/restart/coupled_cluster>` |
  :doc:`DFCOEF and DFPCMO <tutorials/restart/dfcoef_and_dfpcmo>` |
  :doc:`Troubleshooting <tutorials/restart/troubleshooting>`

* **2-component Hamiltonians:**
  :doc:`X2C and local X2C <tutorials/two_component_hamiltonians/x2c_mol_loc>` |
  :doc:`Molecular mean-field X2C <tutorials/two_component_hamiltonians/molecular_mean_field>` |
  :doc:`Selecting a 2-component Hamiltonian other than X2C <tutorials/two_component_hamiltonians/hamiltonians>` |
  :doc:`Case study <tutorials/two_component_hamiltonians/example>` 

* **Relativistic effective core potentials:**
  :doc:`Getting started <tutorials/ecp/first>`

* **Nonrelativistic limit:**
  :doc:`Reproducing nonrelativistic results <tutorials/reproducing_nr_results>`

* **DFT:**
  :doc:`TDDFT <tutorials/dft/tddft>` |
  :doc:`BSSE <tutorials/dft/bsse>` |
  :doc:`CAM functional <tutorials/dft/cam>` |
  :doc:`Troubleshooting <tutorials/dft/troubleshooting>`
  
* **Frozen orbitals:**
  :doc:`Frozen orbitals <tutorials/frozen_orbitals/tutorial>` 

* **Long-range WFT/short-range DFT:**
  :doc:`General <tutorials/srDFT/general>`

* **Property calculations:**
  :doc:`Calculation of NMR shieldings using simple magnetic balance <tutorials/simple_magnetic_balance/tutorial>` |
  :doc:`An introduction to complex reponse <tutorials/properties/complex_response>` |
  :doc:`Magnetizabilities with London Atomic Orbitals <tutorials/lao_magnetizabilities/lao-magnetizabilities>` |
  :doc:`Dipole moment and polarizability of open-shell molecule <tutorials/properties/dipmom-polariz>` |
  :doc:`Calculation of nuclear spin-rotation constants <tutorials/spinrot/tutorial>` |
  :doc:`Calculation of molecular rotational g-tensors <tutorials/rotg/tutorial>`

* **X-ray spectroscopy:**
  :doc:`Core ionization in the CO and N2 molecules <tutorials/x_ray/CO_N2_IP1s/tutorial>` |
  :doc:`Core electron excitations and ionization in water at the HF and DFT levels <tutorials/x_ray/water_Kedge/tutorial>`

* **Spectroscopy:**
  :doc:`Electronic excitations using the POLPRP module <tutorials/polprp/general>` |
  :doc:`Full scope application: Excitation spectrum of an osmium complex <tutorials/polprp/osmium>`

* **Analysis:**
  :doc:`Projection analysis <tutorials/analysis/projection_analysis>`

* **Visualization:**
  :doc:`General overview <tutorials/visual/general/tutorial>` |
  :doc:`Orbital densities <tutorials/visual/orbital_densities/tutorial>` |
  :doc:`Magnetizability density <tutorials/visual/magnetizability_density/tutorial>` |
  :doc:`Molecular electrostatic potential <tutorials/visual/mep/tutorial>` |
  :doc:`Radial distributions <tutorials/visual/radial_distributions/tutorial>`     

* **Open-shell SCF:**
  :doc:`Basics <tutorials/open_shell_scf/aoc>` |
  :doc:`Converging atoms <tutorials/open_shell_scf/converging_atoms>` 

* **Coupled-Cluster:**
  :doc:`The high spin oxygen molecule <tutorials/highspin_cc/O2>` |
  :doc:`Coupled Cluster memory count <tutorials/cc_memory_count/count_cc_memory>` |
  :doc:`Hybrid-parallel run <tutorials/hybrid_parallel_cc_run/hybrid_parallel_cc_run>` 

* **Case studies:**
  :doc:`Ir(16+) cation <tutorials/case_studies/Ir_16plus>` |
  :doc:`CmF molecule <tutorials/case_studies/CmF/CmF>` |
  :doc:`MnO6 system <tutorials/case_studies/MnO6>` |
  :doc:`UF6 molecule <tutorials/case_studies/UF6_molecule/UF6>` |
  :doc:`UF6(-) anion <tutorials/case_studies/UF6_anion/UF6_anion>` |
  :doc:`UO6(-6) anion <tutorials/case_studies/UO6_6minus/UO6_6minus>` |
  :doc:`LuF3 molecule <tutorials/case_studies/LuF3/LuF3>`

* **ADC:**
  :doc:`Triatomic molecule <tutorials/reladc/triatomic_molecule>` |
  :doc:`Atom <tutorials/reladc/atom>`

* **ECP:**
  :doc:`First calculation <tutorials/ecp/first>` |
  :doc:`Correlation calculations <tutorials/ecp/combine>`

* **Polarizable continuum model:**
  :doc:`Basics <tutorials/pcm/pcm_basics>` |
  :doc:`Hartree-Fock and DFT calculations in solution with the polarizable continuum model <tutorials/pcm/pcm_scf>` |
  :doc:`Calculation of polarizabilities in solution: response theory approach <tutorials/pcm/pcm_response>`

* **Polarizable embedding (PE) model:**
  :doc:`PE-HF calculations on micro-solvated H2O <tutorials/polarizable_embedding/pe_scf>` |
  :doc:`PE-TDDFT calculations on micro-solvated H2O <tutorials/polarizable_embedding/pe_response>`

* **Frozen density embedding (FDE):**
  :doc:`NMR shieldings in Frozen Density Embedding (FDE) scheme with London atomic orbitals (LAOs) <tutorials/fde_nmr/tutorial>`
 
* **Davidson corrections for relCI (LUCITA and KRCI):**
  :doc:`+Q corrections <tutorials/relci_q_corrections/q_corrections>`

* **Utility programs:**
  :doc:`TWOFIT <tutorials/utils/twofit>` |
  :doc:`VIBCAL <tutorials/utils/vibcal>` |
  :doc:`CFREAD <tutorials/utils/cfread>` |
  :doc:`LABREAD <tutorials/utils/labread>`

* **Python interface with OpenFermion:**
  :doc:`Openfermion-Dirac <tutorials/openfermion/openfermion-dirac>`

* **Outdated tutorials (need update):**
  :doc:`DIRRCI <tutorials/dusty_attic/dirrci>` |
  :doc:`GOSCIP <tutorials/dusty_attic/goscip>` |
  :doc:`LUCITA <tutorials/dusty_attic/lucita>` |
  :doc:`MOLTRA <tutorials/dusty_attic/moltra>` |
  :doc:`MP2 <tutorials/dusty_attic/mp2>` |
  :doc:`Open shells <tutorials/dusty_attic/open_shells>`


Exercises
=========

.. toctree::
   :maxdepth: 1

   exercises/index.rst


Developers
==========

* **Code review:**
  :doc:`Code review workflow <programmers/code_review>`

* **Releasing:**
  :doc:`Release preparation <programmers/release_preparation>` |
  :doc:`Beta testing <programmers/beta_testing>` |
  :doc:`Where to commit after the release is out <programmers/commit_policy>`

* **Programming:**
  :doc:`Rules <programmers/rules>` |
  :doc:`How to add new tests <programmers/testing>` |
  :doc:`How to add/move/remove sources <programmers/adding_moving_removing_sources>` |
  `runtest_dirac.py <http://runtest.readthedocs.org/en/latest/doc/adding_tests/dirac.html>`_ |
  :doc:`Input parsing <programmers/input_reading>` |
  :doc:`Git submodules <programmers/external_projects>` |
  :doc:`Nightly tests <programmers/cdash>` |
  :doc:`Dirac on Windows <programmers/windows>` |
  :doc:`Profiling <programmers/profiling>` |
  :doc:`Debugging <programmers/debugging>` |
  :doc:`History <programmers/museum>` |
  :doc:`Further development <programmers/further_development>` |
  :doc:`Good Fortran 90 practices <programmers/good_fortran_90_practices>` |
  :doc:`FAQ <programmers/faq>`

* **Basis sets:**
  :doc:`Testing basis sets in DIRAC <programmers/test_basis_sets>` 

* **Moving code between machines:**
  :doc:`Transfering uncommitted code <programmers/transfering_uncommitted_code>`

* **Notes:**
  :doc:`DFCOEF <programmers/dfcoef>` |
  :doc:`Screening <programmers/screening>` |
  :doc:`64bit integers <programmers/int8>` |
  :doc:`Numerical constants <programmers/num_constants>` |
  :doc:`XML <programmers/xml>` |
  :doc:`Static linking <programmers/static_linking>` |
  :doc:`Problems with lapack's DSYEVR <programmers/diag>` |
  :doc:`How this documentation works <programmers/documentation_howto>`


Bibliography
============

* :doc:`References <zreferences>`
