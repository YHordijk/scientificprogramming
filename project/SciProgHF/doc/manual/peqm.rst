:orphan:

star(PEQM)

Polarizable embedding model

This input section controls calculations using the polarizable embedding (PE) model.
In DIRAC it is possible to include the effects from a structured environment on a core
molecular system using the polarizable embedding (PE) model. The current implemen-
tation is a layered QM/MM-type embedding model capable of using advanced potentials
that include an electrostatic component as well as an induction (polarization) component.
The effects of the environment are included through effective operators that contain an em-
bedding potential, which is a representation of the environment, thereby directly affecting
the molecular properties of the core system. The wave function of the core system is opti-
mized while taking into account the explicit electrostatic interactions and induction effects
from the environment in a fully self-consistent manner. The electrostatic and induction
components are modeled using Cartesian multipole moments and anisotropic dipole-dipole
polarizabilities, respectively. The electrostatic part models the permanent charge distribu-
tion of the environment and will polarize the core system, while the induction part also
allows polarization of the environment. The environment response is included in the effec-
tive operator as induced dipoles, which arise due to the electric fields from the electrons and
nuclei in the core system as well as from the environment itself. It is therefore necessary
to recalculate the induced dipoles according to the changes in the electron density of the
core system as they occur in a wave function optimization. Furthermore, since the induced
dipoles are coupled through the electric fields, it is necessary to solve a set of coupled linear
equations. This can be done using either an iterative or a direct solver. This also means
that we include many-body effects of the total system.

The multipoles and polarizabilities can be obtained in many different ways. It is
possible to use the molecular properties, however, usually distributed/localized properties
are used because of the better convergence of the multipole expansion. These are typically
centered on all atomic sites in the environment (and sometimes also bond-midpoints), how-
ever, the implementation is general in this sense so they can be placed anywhere in space.
Currently, the PE library supports multipole moments up to fifth order and anisotropic
dipole-dipole polarizabilities are supported. For multipoles up to and including third order
(octopoles) the trace will be removed if present. Note, that the fourth and fifth order multi-
pole moments are expected to be traceless. In case polarizabilities are included it might be
necessary to use an exclusion list to ensure that only relevant sites can polarize each other.
The format of the POTENTIAL.INP file is demonstrated below.

Keywords
========

keyword(POTENTIAL)

This option can be used to specify a non-default name of the POTENTIAL.INP in-
put file that contains the embedding potential parameters, e.g.::

   .POTENT
   peqm.pot

The default name is POTENTIAL.INP.

keyword(DIRECT)

Use the direct solver to determine the induced dipole moments. It will explicitly build
a classical response matrix of size :math:`3S(3S + 1)/2`, where :math:`S` is the number of polarizable
sites and is therefore only recommendable for smaller molecular systems. Note that
this solver is not parallelized. The default is to use the iterative solver

keyword(ITERATIVE)

Use the iterative solver to determine the induced dipole moments. This is the
default. The convergence threshold defaults to :math:`1.0\cdot 10^{-8}>\sum^S_{s=1}|\mu^{(k)}_s-\mu^{(k-1)}_s|`, where :math:`k` is the 
current iteration, but can also be provided with this option::

  READ (LUCMD, *) THRITER (optional)

keyword(BORDER)

Controls the handling of the border between the core molecular system and its en-
vironment described by the embedding potential.::

   1) READ (LUCMD, *) BORDER_TYPE, RMIN, AUORAA
   2) READ (LUCMD, *) BORDER_TYPE, REDIST_ORDER, RMIN, AUORAA, NREDIST

There are two mutually exclusive schemes:
   1) BORDER_TYPE = REMOVE
   2) BORDER_TYPE = REDIST

The first option will remove all multipoles and polarizabilities that are within the given distance "RMIN" from any atom in the core molecular system. The "AUORAA" variable specifies
whether the distance threshold is given in ångström ("AA") or bohr ("AU"). The second
option will redistribute parameters that are within the given threshold "RMIN" from any
atom in the core system to nearest sites in the environment. The order of multipoles
up to which will be redistributed is determined by the "REDIST_ORDER" variable, e.g.
"REDIST_ORDER = 1" means that only charges will be redistributed and all other pa-
rameters removed, "REDIST_ORDER = 2" means charges and dipoles are redistributed
and so on. The sign of "REDIST_ORDER" specifies if the polarizabilities are redistributed.
Positive means that the polarizabilities are removed and negative means redistributed.
The number of sites that parameters on a given site are redistributed to is determined
by the "NREDIST" variable which can be between 1 and 3. The default is to redistribute
charges within 2.2 bohr to its nearest site and removing all other parameters.

keyword(DAMP INDUCED)

Damp the electric fields from induced dipole moments using Thole’s exponential damp-
ing scheme.::

  READ (LUCMD, *) IND_DAMP (optional)
  
The default damping coefficient is the standard 2.1304.

keyword(DAMP MULTIPOLE)

Damp the electric fields from permanent multipole moments using Thole’s expo-
nential damping scheme.::

  READ (LUCMD, *) MUL_DAMP (optional)

This option requires polarizabilities on all sites with permanent multipole moments. The default damping coefficient is the standard 2.1304.

keyword(DAMP CORE)

Damp the electric fields from the electrons and nuclei in the core region based on
Thole’s exponential damping scheme::

   READ (LUCMD, *) CORE_DAMP (optional)
   READ (LUCMD, *) NALPHAS
   DO I = 1, NALPHAS
      READ (LUCMD, *) ISO_ALPHA(I)
   END DO

The damping coefficient is optional unless isotropic polarizabilities are present. The default damping coefficient is the standard 2.1304. Standard polarizabilities from :cite:`vanduijnen1998`
have been implemented and will be used if none are given as input. However, only H, C, N, O, F, S, Cl, Br and I are available.

keyword(GSPOL)

Activate the ground-state polarization approximation, i.e. freeze the embedding po-
tential according to the ground-state density. This means that the polarizable envi-
ronment is self-consistently relaxed during the optimization of the ground-state den-
sity/wave function of the core molecular system and then kept frozen in any following
response calculations

keyword(NOMB)

Remove many-body effects in the environment. This is done by deactivating interac-
tions between inducible dipole moments.

keyword(RESTART)

Use any existing files to restart calculation.

keyword(CUBE)

Create cube file for the core molecular system containing the electrostatic potential due
to the final converged polarizable embedding potential.::

   1) READ (LUCMD, *) STD_GRID
   2) READ (LUCMD, *) OPTION
   2) IF (OPTION == ’GRID’) THEN
   2)    READ (LUCMD, *) XSIZE, XGRID, YSIZE, YGRID, ZSIZE, ZGRID
   2) END IF
   IF (OPTION == ’FIELD’) THEN
      FIELD = .TRUE.
   END IF

The grid density can be spec-
ified either using 1) standard grids (COARSE (3 points/bohr), MEDIUM (6 points/bohr)
or FINE (12 points/bohr) and in all cases 8.0 bohr are added in each direction) or
2) the GRID option which gives full control of the cube size and density, and requires 
an additional input line specifying the extent added (in bohr) in each direction and
density (in points/bohr). The default grid density is MEDIUM and default cube size
is the extent of the molecule plus 8.0 bohr in plus and minus each Cartesian coordi-
nate. If the FIELD option is given then also three cube files containing the Cartesian
components of the electric field from the embedding potential will be created.  

keyword(ISOPOL)

Converts all anisotropic polarizabilities into isotropic polarizabilities.

keyword(ZEROPOL)

Remove all polarizabilities

keyword(ZEROMUL)

Remove all multipoles of order ZEROMUL_ORDER and up.::

   READ (LUCMD, *) ZEROMUL_ORDER (optional)

Remove all multipoles of order ZEROMUL_ORDER and up. The default is to remove all
dipoles and higher-order multipoles (i.e. ZEROMUL_ORDER = 1).

keyword(VERBOSE)

Verbose output. Currently this will print the final converged induced dipole moments.

keyword(DEBUG)

Debug output. Prints the total electric field and the induced dipole moments in each
iteration. WARNING: for large systems this will produce very large output files.


   
The potential input format
==========================

The POTENTIAL.INP file is split into three sections: @COORDINATES, @MULTIPOLES and
@POLARIZABILITIES. The format is perhaps best illustrated using an example:

.. literalinclude:: ../tutorials/polarizable_embedding/3_h2o.pot

@COORDINATES
------------

The coordinates section follows the standard XYZ file format so that the environment can
be easily visualized using standard programs. The first line in gives the total number of
sites in the environment and the second line specifies whether the coordinates are given in
ångström (AA) or bohr (AU). The rest of the coordinates section is a list of the sites in the
environment where each line contains the element symbol and x-, y- and z-coordinates of
a site. If a site is not located on an atom, e.g. if it is a bond-midpoint, then the element
symbol should be specified as X. The listing also gives an implicit numbering of the sites,
so that the first line is site number one, the second line is site number two and so on. This
numbering is important and used in the following sections.

@MULTIPOLES
-----------
The multipoles section is subdivided into the orders of the multipoles, i.e. ORDER 0
for monopoles/charges, ORDER 1 for dipoles and so on. For each order there is a number
specifying the number of multipoles of that specific order. Note, that this number does
not have to be equal to the total number of sites. This is followed by a list of multipoles
where each line gives the multipole of a site. The lines begin with a number that specifies
which site the multipole is placed. Only the symmetry-independent Cartesian multipoles
(given in a.u.) should be provided using an ordering such that the components are stepped
from the right, e.g. xx xy xz yy yz zz or xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz.
Note, that the multipoles should in general be traceless, however, for multipoles up to and
including octopoles (i.e. ORDER 3) the trace is removed if present. Furthermore, the current
implementation is limited to fifth order multipoles.

@POLARIZABILITIES
-----------------

The polarizabilities section is also subdivided into orders, i.e. ORDER 1 1 for dipole-dipole
polarizabilities, which is the only type supported in the current release. The format is the
same as for multipoles, i.e. first line contains number of polarizabilities which is followed by
a list of the polarizabilities using the same ordering as the multipoles. The polarizabilities
should also be given in a.u. In addition, there is also the exclusion lists (EXCLISTS section).
Here the first line gives the number of lists (i.e. the number of lines) and the length of the
exclusion lists (i.e. the number of entries per line). The exclusion lists specify the polar-
ization rules. There is a list attached to each polarizable site that specifies which sites are
not allowed to polarize it, e.g. 1 2 3 4 5 means that site number 1 cannot be polarized by
sites 2, 3, 4 and 5.
