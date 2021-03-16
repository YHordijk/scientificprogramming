:orphan:
 

star(X2C)

Feature and debug options for the Exact 2-Component
relativistic Hamiltonian module, Refs. :cite:`Knecht2010` and :cite:`Knecht2014`.

keyword(fragx2c)

use the local X2C scheme (DLU approach) of Peng and Reiher :cite:`Peng2012` to construct an approximate molecular U matrix 
from a superposition of atomic U matrices. This requires to compute in advance atomic U matrices in atomic calculations 
which are stored in the file X2CMAT. For each atom type (C, H, F, U, ...) the module expects an X2CMAT file 
with the atomic number as appropriate ending, for example: X2CMAT.092 (for U), X2CMAT.006 (C), X2CMAT.001 (H) etc.. 
See also the tutorial section for a working example.    

keyword(skippc)

Skip the picture change transformation of four-component operators. This is for restarting X2c calculations after SCF, SCF+MOLTRA steps,
when such operators transformations were already carried out.

keyword(DEBUG)

raises the print level in the X2C module to print intermediate matrices and arrays. Use with care as the output might
become filled with lots of numbers.

keyword(freeh1)

use the free-particle Hamiltonian as defining Hamiltonian to solve the initial matrix problem in the X2C algorithm. for
debug purposes mostly.
