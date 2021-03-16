Examples of xyz-based symmetry detection
========================================

We provide  (in test/xyz_symmetry_recognition) collection of geometries aimed 
for the DIRAC automatic symmetry recognition
of molecular structures given in the xyz-coordinates format.

Sources of xyz-inputs
---------------------

Molecules were handled by a suitable software and 
their geometries (coordinates), reflecting imposed symmetry, were exported,
preferably in the xyz-format. 

If the GUI software does not save/export the geometry coordinates in xyz-format, try to save (export) them in another
common format, and read it in an other software, capable to export the desired xyz file.
Note that sometimes the resulting xyz-geometry file had to be manually controlled and coordinates 
corrected, if necessary, to keep the symmetry of the system.

Very good choice for getting symmetry imposed xyz coordinates is the GUI of the commercial ADF software (www.scm.com) .

List of symmetries
------------------

The DIRAC's *hersym.F/FIND_PGROUP* subroutine can handle the following list of point (full) groups:

* O(h) 
* I(h) 
* T(d)
* D(oo,h) 
* D(nh)
* D(nd)
* D(n)
* C(oo,v)
* C(nv)
* C(nh)
* C(n)
* C(s)
* C(i)
* C(1)
* S(n)

The detected full group is then turned into lower computational group,
represented as the point group.


Examples
--------

The following table contains examples of 
the symmetry detection from the xyz coordinates.

==========   ==========================     ==========    ==============
xyz input    molecule; imposed symmetry     full group    represented as 
==========   ==========================     ==========    ==============
molecule01   heterocycycle;     D(2h)       D(2h)         D2h
molecule02   Re2(CO)10;         D(4d)       C(4v)  (!)    C2v
molecule03   Os(CO)5;           D(3h)       C(2v)  (!)    C2v
molecule04   acetylene;         D(oo,h)     D(oo,h)       D2h
molecule05   SF6;               O(h)        O(h)          D2h
molecule06   allene;            S(4)        D(2d)  (!)    D2
molecule07   NH3;               C(3v)       C(3v)         Cs
molecule08   PF5;               D(3h)       D(3h)         C2v
molecule09   Hs(CO)4; {close to D(2d)}      D(2)   (!)    D2
molecule10   ethane;            D(3d)       C(2h)  (!)    C2h
molecule11   dodecahedrane;     I(h)        C(2h)  (!)    C2h
molecule12   dodecahedrane;     D(5d)       D(5d)         C2h
molecule13   benzene;           D(6h)       D(6h)         D2h 
molecule14   methane;           T(d)        T(d)          D2
molecule15   H2O2;              C(2h)       C(2h)         C2h
molecule16   ferrocene;         D(5h)       C(s)   (!)    C2v
molecule17   ferrocene-stagg;   D(5d)       C(2h)  (!)    C2h
molecule18   Hs(CO)4;           D(2d)       D(2d)         D2
==========   ==========================     ==========    ==============

The (!) symbol means that the (ADF's determined) symmetry does not correspond with the DIRAC's detected full group.

More examples have to be be added here to demostrate (and test) the full variety of symmetry detections in DIRAC.
