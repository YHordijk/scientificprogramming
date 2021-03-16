:orphan:

Keyword reference
=================

* :doc:`**DIRAC <dirac>` --- Job specification
   - :doc:`*OPTIMIZE <optimize>` --- Geometry optimization directives

* :doc:`**GENERAL <general>` --- General input
   - :doc:`*PARALLEL <parallel>` --- Parallel run directives

* :doc:`**MOLECULE<molecule>` --- Specification of basis set and other atom-specific information
   - :doc:`*CHARGE <molecule>` --- Define the molecular charge
   - :doc:`*BASIS <molecule>` --- Define default and atom-specific basis sets
   - :doc:`*CENTERS <molecule>` --- Define specific properties for atoms (e.g. ghost centers, point charges)
   - :doc:`*SYMMETRY <molecule>` --- Set symmetry manually (when reading from xyz file)
   - :doc:`*COORDINATES <molecule>` --- Set unit for the input coordinates

* :doc:`**HAMILTONIAN <hamiltonian>` --- Specify the Hamiltonian
   - :doc:`*DFT <dft>` --- DFT directives
   - :doc:`*AMFI <amfi>` --- AMFI directives
   - :doc:`*X2C <x2c>` --- X2C directives
   - :doc:`*FDE <fde>` --- Frozen Density Embedding directives
   - :doc:`*PCM <pcm>` --- Polarizable Continuum Model directives
   - :doc:`*PCMSOL <pcmsolver>` --- PCMSolver module directives
   - :doc:`*PEQM <peqm>` --- Polarizable embedding model directives
    
* :doc:`**WAVE FUNCTION <wave_function>` --- Method specification
   - :doc:`*SCF <wave_function/scf>` --- SCF module (Hartree-Fock/Kohn-Sham)
   - :doc:`*MP2CAL <wave_function/mp2cal>` --- second-order Møller-Plesset perturbation theory
   - :doc:`*RESOLVE <wave_function/resolve>` --- resolve open-shell states
   - :doc:`GOSCIP <wave_function/goscip>`  --- Complete Open-Shell module
   - :doc:`*COSCI  <wave_function/cosci>`   --- Complete Open-Shell CI module
   - :doc:`DIRRCI <wave_function/dirrci>` --- Direct Relativistic CI module
   - :doc:`*LUCITA <wave_function/lucita>` --- Spinfree large-scale CI module
   - :doc:`*KR CI <wave_function/krci>` --- Kramers-restricted relativistic large-scale CI module
   - :doc:`*KRMCSCF <wave_function/krmcscf>` --- Kramers-restricted relativistic large-scale MCSCF module
   - :doc:`**RELCCSD <wave_function/relccsd>` --- Coupled cluster module
   - :doc:`**RELADC <wave_function/reladc>` --- Propagator module (ADC) for single and double ionizations
   - :doc:`POLPRP <wave_function/polprp>` --- Polarization Propagator module (ADC) for excitations
   - :doc:`*MVOCAL <wave_function/mvocal>` --- Modified virtual orbitals
   - :doc:`*MP2 NO <wave_function/mp2no>` --- MP2 natural orbitals module
   - :doc:`*LAPLCE <wave_function/laplce>` --- Laplace transformation of orbital-energy denominators
   - :doc:`**EXACC <wave_function/exacorr>` --- Parallel implementation of coupled cluster methods based on ExaTensor library.

* :doc:`**ANALYZE <analyze>` --- Analyze the wave function
   - :doc:`*MULPOP <analyze/mulpop>` --- Mulliken population analysis
   - :doc:`*PRIVEC <analyze/privec>` --- Print coefficients
   - :doc:`*PROJECTION <analyze/prjana>` --- Projection analysis
   - :doc:`*LOCALIZATION <analyze/localization>` --- Localization 
   - :doc:`*DENSITY <analyze/density>` --- Density
   - :doc:`*RHO1 <analyze/density>` --- Rho1

* :doc:`**PROPERTIES <properties>` --- Property module
   - :doc:`*EXPECTATION VALUES <properties/expectation>` --- Expectation values
   - :doc:`*EXCITATION ENERGIES <properties/excitation>` --- Excitation energies
   - :doc:`*LINEAR RESPONSE <properties/linear_response>` --- Linear response
   - :doc:`*QUADRATIC RESPONSE <properties/quadratic_response>` --- Quadratic response
   - :doc:`*MOLGRD <properties/molgrd>` --- Molecular gradient
   - :doc:`*NMR <properties/nmr>` --- NMR directives
   - :doc:`*STEX <properties/stex>` --- Static exchange (STEX) directives


* :doc:`**VISUAL <visual>` --- Visualization module

* :doc:`**INTEGRALS <integrals>` --- Integral directives

* :doc:`**GRID <grid>` --- Numerical integration grid

* :doc:`**MOLTRA <moltra>` --- Integral transformation module

Notes:
------

  :doc:`One-electron operators <one_electron_operators>` |
  :doc:`Orbital strings <orbital_strings>` |
  :doc:`Symmetry-handling at the correlated level <groupchain>`
