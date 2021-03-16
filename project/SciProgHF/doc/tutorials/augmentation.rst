:orphan:
 

Basis sets augmentation
=======================

DIRAC suite offers possibility to easily modify employed basis sets.

The purpose of this tutorial, which is reflected in the test *basis_automatic_augmentation*,  
is to show nonstandard automatic augmentation of provided standard basis sets, like
cc-pVXZ, Turbomole and Dyall's.

In the following table we demonstrate the influence of the (contracted) basis set extension on the total 
nonrelativistic energy and on the dipole moment, which is simply calculated through SCF wave-function expectation value.

Two molecules were chosen: HF and HBr. For the former molecule, the pre-cc-pVXZ and Turbomole basis sets are used for both F and H atoms.
In the latter molecule, the heavy Br atom is described with Dyall's relativistic basis set, while the H atom is provided
with the nonrelativistic basis set.

+---------------------+--------------------------------+--------------------+---------------+-------------+
|  Basis set name     | Atom    Basis size             | Atom   Basis size  | E(NR-SCF)/a.u.| dip. mom/D  |
+=====================+================================+====================+===============+=============+
|    aug-cc-pVDZ      | F      [10s5p2d|4s3p2d]        |  H    [5s2p|3s2p]  | -100.0337931  | -1.89747844 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
|  d-aug-cc-pVDZ      | F      [11s6p3d|5s4p3d]        |  H    [6s3p|4s3p]  | -100.0339787  | -1.88741792 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
|  t-aug-cc-pVDZ      | F      [12s7p4d|6s5p4d]        |  H    [7s4p|5s4p]  | -100.0341087  | -1.88839962 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
| Turbomole-DZP       | F      [8s4p1d|4s2p1d]         |  H    [4s1p|2s1p]  | -100.0028441  | -1.94418011 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
| s-a-Turbomole-DZP   | F      [9s5p2d|5s3p2d]         |  H    [5s2p|3s2p]  | -100.0180563  | -1.91810891 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
|  Turbomole-TZVPP    | F      [11s6p2d1f|5s3p2d1f]    |  H [5s2p1d|3s2p1d] | -100.0657904  | -1.91940698 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
| s-a-Turbomole-TZVPP | F      [12s7p3d2f|6s4p3d2f]    |  H [6s3p2d|4s3p2d] | -100.0668454  | -1.89121376 | 
+---------------------+--------------------------------+--------------------+---------------+-------------+
|       dyall.v2z     | Br     [15s11p7d1f|15s11p7d1f] | (H [4s1p|2s1p])    | -2572.6361195 | -0.97693933 |
+---------------------+--------------------------------+--------------------+---------------+-------------+
| d-aug-dyall.v2z     | Br     [17s13p9d3f|17s13p9d3f] | (H [4s1p|2s1p])    | -2572.6425365 | -0.76521484 |
+---------------------+--------------------------------+--------------------+---------------+-------------+

(H-atom in parenthesis has cc-pVDZ basis sets.)

Note that the nonrelativistic Hamiltonian with contracted basis sets was used in all these examples.
Contracted basis sets were set for light elements (F,H), while for Br atom uncontracted scheme with Dyall's basis
was preferred.

For relativistic calculation with either 2-component of 4-component Hamiltonians the nonrelativistic contraction
of light elements basis sets is no longer suitable, especially if the molecule contains heavy elements.
The user is adviced to resort to uncontracting his basis sets. 
For example, for the HBr molecule you can combine
decontracted cc-pVDZ basis set for hydrogen with the decontracted (as is) dyall.v2z basis for bromine.

Dyall's basis sets, which are constructed without contractions, can be used also with the X2C Hamiltonian.
