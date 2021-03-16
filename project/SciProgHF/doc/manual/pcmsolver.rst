:orphan:
 

star(PCMSOL)

PCMSolver module configuration directives.

General information on the PCM can be found in :cite:`Tomasi2005`.
Details on the module can be found `here <http://pcmsolver.github.io/pcmsolver-doc/>`_
The module source code is freely available under the LGPLv3 license `on this URL <https://github.com/PCMSolver/pcmsolver>`_

There are two ways of providing input to the PCMSolver module:

1. by means of an additional input file, parsed by the `pcmsolver.py` script;
2. by means of a special section in the DIRAC input.

Method 1 is more flexible: all parameters that can be modified by the user are available.
Method 2 just gives access to the core parameters. We will here describe Method 2, for a description
of Method 1 refer to the `PCMSolver input documentation <http://pcmsolver.github.io/pcmsolver-doc/_input_description.html>`_
Examples of both methods are given in the PCM tutorial section.
All keywords expecting a numerical value are intended to be in Angstrom. The module will perform conversion
to atomic units internally.

.. warning:: Do **NOT** truncate keyword values to six characters. It will result in a crash with an obscure error message.

.. warning:: Some keywords can be set but won't have any effect on the module as their respective functionality is not in the released version of PCMSolver.

**Cavity directives**
=====================

The generating spheres for the cavity are obtained by centering spheres on the atoms making up the
molecule. Radii are automatically selected from a built-in set and are scaled by 1.2.
For a finer control over the cavity creation details use Method 1.

keyword(CAVTYP)

The type of the cavity. Completely specifies type of molecular surface and its discretization.
Allowed values are GEPOL and RESTART.
If RESTART is given, teh module will read the file specified by the RESTAR keyword 
and create a GePol cavity from that. Default::
 
  .CAVTYP
  GEPOL

keyword(PATCHL)

.. warning:: not available in the current release of the module.

keyword(COARSI)

.. warning:: not available in the current release of the module.

keyword(AREATS)

Average area of the surface partition for the GePol cavity. Default::

  .AREATS
  0.3

keyword(MINDIS)

.. warning:: not available in the current release of the module.

keyword(DERORD)

.. warning:: not available in the current release of the module.

keyword(NOSCAL)

If present the radii for the spheres will **not** be scaled by 1.2. Radii are scaled by default.  

keyword(RADIIS)

The set of built-in radii to be used. Allowed values are BONDI and UFF. 
Refer to the `PCMSolver input documentation <http://pcmsolver.github.io/pcmsolver-doc/_input_description.html>`_ 
for further details about the radii sets. Default::

  .RADIIS
  BONDI

keyword(RESTAR)

The name of the `.npz` file to be used for the GePol cavity restart.

keyword(MINRAD)

Minimal radius for additional spheres not centered on atoms. An arbitrary big value is equivalent
to switching off the use of added spheres. By default, no spheres not centered on atoms are added. Default::

  .MINRAD
  100.0

**Solver directives**
=====================

keyword(SOLVER)

Type of solver to be used. All solvers are based on the Integral Equation Formulation of
the Polarizable Continuum Model:

* IEFPCM. Collocation solver for a general dielectric medium;
* CPCM. Collocation solver for a conductor-like approximation to the dielectric medium.

keyword(SOLVNT)

Specification of the dielectric medium outside the cavity. This keyword **must always** be
given a value. Allowed values are:

* WATER               
* METHANOL                                    
* ETHANOL                                  
* CHLOROFORM                                  
* METHYLENECHLORIDE                          
* 1,2-DICHLOROETHANE                        
* CARBON TETRACHLORIDE                         
* BENZENE                                      
* TOLUENE                                   
* CHLOROBENZENE                              
* NITROMETHANE                               
* N-HEPTANE                                   
* CYCLOHEXANE                                 
* ANILINE                                   
* ACETONE                                    
* TETRAHYDROFURANE                              
* DIMETHYLSULFOXIDE                            
* ACETONITRILE                                
* EXPLICIT 

For further details on the dielectric properties of the available built-in solvents refer to the 
`PCMSolver input documentation <http://pcmsolver.github.io/pcmsolver-doc/_input_description.html>`_ 
If the solvent name given is different from Explicit any other settings in the Green's function
section will be overridden by the built-in values for the solvent specified.
Solvent = Explicit, triggers parsing of the Green's function directives.

keyword(EQNTYP)

.. warning:: not available in the current release of the module.

keyword(CORREC)

Correction, :math:`k` for the apparent surface charge scaling factor in the CPCM solver
:math:`f(\varepsilon) = \frac{\varepsilon - 1}{\varepsilon + k}`. Default::

  .CORREC
  0.0

keyword(PROBER)

Radius of the spherical probe approximating a solvent molecule. Used for generating the
solvent-excluded surface (SES) or an approximation of it. Overridden by the built-in value
for the chosen solvent. Default::

  .PROBER
  1.0

**Green's function directives**
===============================

If solvent is Explicit, **both** the Green's function inside and outside
must be specified.
The Green's function inside will always be the vacuum, while the 
Green's function outside might vary. Currently, only isotropic uniform dielectrics
are implemented.

keyword(GINSID)

Type of Green's function inside the cavity. Default::

  .GINSID
  VACUUM

keyword(GOUTSI)

Type of Green's function inside the cavity. Default::

  .GINSID
  UNIFORMDIELECTRIC

keyword(EPSILO)

Static dielectric permittivity for the medium outside the cavity. Default::

  .EPSILO
  1.0
