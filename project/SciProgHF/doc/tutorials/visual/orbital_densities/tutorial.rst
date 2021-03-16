:orphan:
 

Plotting orbital densities
==========================

In this tutorial we will plot the HOMO density of the water molecule.  We will
do this in two steps: first we will run the SCF and save the DFCOEF file, then
we will restart from this and calculate the orbital density on a grid of points
and write this into a cube file which can be opened by your favorite molecular
visualization program (the author of this tutorial likes Avogadro).


Getting the wave function
-------------------------

Let us use the following molecule input (h2o.mol):

.. literalinclude:: h2o.mol

Together with a simple job input (scf.inp):

.. literalinclude:: scf.inp

We can get the wave function with the following run script:

.. literalinclude:: scf.run

After running the script you should see the file DFCOEF in your submit
directory.


Plotting the orbital density
----------------------------

Water has 5 orbitals and we will plot the HOMO, orbital nr. 5 (density.inp):

.. literalinclude:: density.inp

The orbital occupation is controlled by the following keyword::

  .OCCUPATION
   1
   1 5 1.0

This means that we have one occupation set (first "1"), this occupation set is
for Fermion irrep 1 (second "1"), the range is orbital 5 ("5") and it is fully
occupied. To get the total density we could remove this keyword or alternatively::

  .OCCUPATION
   1
   1 1-5 1.0

We will use the following run script to generate the cube file:

.. literalinclude:: density.run

Finally you can open plot.3d.cube with your favorite program.
