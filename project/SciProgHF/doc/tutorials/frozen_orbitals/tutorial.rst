:orphan:
 

Frozen orbitals
===============

DIRAC allows the freezing (or elimination) of orbitals during Hartree-Fock or Kohn-Sham calculations.
We shall illustrate this functionality using bonding in the water molecule as an example and essentially
reproducing a calculation reported by Jarvie and co-workers in 1973 :cite:`Jarvie1973`.

The water molecule has a bond angle of 104.5 :math:`^{\circ}`. This is sometimes rationalized by invoking
:math:`sp^3` hybridized oxygen. Perfect hybridization would, however, lead to a bond angle of
109.5 :math:`^{\circ}` :math:`=\mathrm{arccos}(-1/3)`, and so the deviation is rationalized by invoking
valence shell electron pair repulsion (VSEPR) theory and in particular the stronger repulsion between
lone pairs.

Jarvie and co-workers therefore proposed to freeze the oxygen :math:`2s` orbital in water and see how
this affected the bond angle. Let us do this calculation the DIRAC way:

Atomic calculation
------------------

First we need an oxygen :math:`2s` orbital. We therefore set up an atomic calculation using the input files *O.inp*

.. literalinclude:: O.inp

and *O.mol*

.. literalinclude:: O.mol

Our plan now is to import this orbital into a subsequent molecular calculation. There is, however, one important restriction: The atomic calculation has to be carried out in the same symmetry as the molecular one, which in this case will be :math:`C_{2v}`, hence the explicit symmetry specification in the *O.mol* input file (for more details, see :ref:`molecule_using_xyz`)	    .
We run the calculation using the command::

     pam --inp=O --mol=O --get "DFCOEF=cf.O"

Geometry optimization with frozen orbitals: how it works
--------------------------------------------------------

For the molecular calculation we use the input files *H2O.inp*

.. literalinclude:: H2O.inp

and *H2O.mol*

.. literalinclude:: H2O.mol

The latter file specifies the experimental geometry of water with bond angle 104.5 :math:`^{\circ}` and bond distance 0.972 Å.
The former file includes the :ref:`OPTIMIZE_.NUMGRA` keyword to activate geometry optimization using a numerical gradient.
This will be needed here since we have included the constraint of a frozen orbital.

Now let us look at the specification of the frozen orbital: We have to tell DIRAC that the orbital was not calculated in
the full molecular basis, rather using only the basis of the oxygen atom. This information is indicated by the keyword
:ref:`SCF_.OWNBAS`. Next, we have to tell where to find the orbital. This is specified using the keyword :ref:`SCF_.PROJECT`
and its arguments. The oxygen :math:`2s` orbital is specified by the coefficients found on the file *DFCOEF* that we recovered
in the atomic calculation and renamed *cf.O*. The oxygen atom is considered a fragment of the water molecule, and the coefficient
file therefore constitute a fragment file. The first argument of the :ref:`SCF_.PROJECT` keyword is the specification of the
number of fragment files to read, which in this case is just one.

For each fragment we now have to specify the name of the fragment file. As seen in the FORTRAN snippet in the manual under the :ref:`SCF_.PROJECT` keyword, the fragment file should be given a six-character name, so we in this case choose *DFOXXX*. 
When the fragment has been calculated in its own basis, as specified by the :ref:`SCF_.OWNBAS` keyword, we specify it
by indicating the number of *symmetry-independent nuclei* it contains. Let us see how this works: Going back to the molecular input file *H2O.mol* file
we see that it gives the coordinate of one oxygen atom and *one* hydrogen atom. We do not need to specify the coordinates of the
second hydrogen atom because it is then known from symmetry. The two hydrogen atoms are symmetry-dependent and have to be treated
together as a fragment. We have, however, selected the oxygen atom as our fragment and specify that our fragment consists of a
singe symmetry-independent center. Oxygen will be chosen since it comes first in the molecular input file.

We now have to specify what orbital(s) to select on the fragment file. This is done through an orbital string (see the section
:ref:`orbital_strings`). The oxygen atom has five occupied orbitals and the :math:`2s` orbital is the second one, as we can for 
instance infer from looking at the Mulliken population analysis in the output of the atomic calculation.

Now we have to understand the :ref:`SCF_.FROZEN` keyword. This will take a little bit more theory:
At each iteration in the SCF-cycle we have to solve the generalized eigenvalue equation

.. math::

   F\mathbf{c} = S\mathbf{c}\varepsilon

where :math:`F` is the Fock (or Kohn-Sham) matrix in the AO-basis (symmetry-adapted or not) and :math:`\mathbf{c}` and :math:`\varepsilon` are the coefficients and eigenvalues we seek. The equation also features the overlap matrix :math:`S`
since the AO-basis is generally not orthogonal. A first step towards the solution of this equation is to introduce a
non-unitary transformation :math:`V` that will take us to an orthonormal basis, that is

.. math::

   V^{\dagger}SV = I

where :math:`I` is the identity matrix. We insert :math:`VV^{-1}=I` in front of the coefficients and pre-multiply with :math:`V^{\dagger}` to give
   
.. math::   
     
   V^{\dagger}FVV^{-1}\mathbf{c} = V^{\dagger}SVV^{-1}\mathbf{c}\varepsilon

Since the matrix :math:`V` transforms the overlap matrix to the identity matrix, the transformed equation simplifies to

.. math::

   F^{\prime}\mathbf{c}^{\prime}=\mathbf{c}^{\prime}\varepsilon;\quad F^{\prime}=V^{\dagger}FV;\quad \mathbf{c}^{\prime}=V^{-1}\mathbf{c} .

This is a normal eigenvalue problem and it suffices now to give the transformed Fock matrix :math:`F^{\prime}` to a standard
routine for matrix diagonalization in order to obtain eigenvalues :math:`\varepsilon` and the transformed eigenvectors
:math:`\mathbf{c}^{\prime}`. In a final step we get the MO-coefficients :math:`\mathbf{c}` by the backtransformation

.. math::

   \mathbf{c}=V\mathbf{c}^{\prime}

This diagonalization procedure is repeated in each iteration of the SCF-cycle and provides the optimization of the molecular orbitals.

The freezing of orbitals is done in the following manner: The imported fragment orbitals are projected out of the
transformation matrix :math:`V`, which means that they are not present in the orthonormal basis of the transformed
Fock matrix :math:`F^{\prime}`. They have in other words been removed from the variational space. If this was our
goal we could in fact have stopped here. This is for example what was done in :cite:`Fossgaard2003b` where the
effect of the lanthanide contraction on bonding in the CsAu molecule was investigated by projecting out the gold
:math:`4f` orbitals from the variational space (hence the name of the :ref:`SCF_.PROJECT` keyword).

In the present case we do not want to delete the oxygen :math:`2s` orbital, rather freeze it.
Therefore we have to put the oxygen :math:`2s` orbital back when backtransforming coefficients.
The :ref:`SCF_.FROZEN` keyword takes as argument an orbital string which indicates in what position one would
like to put the frozen orbitals in the total set of orbitals.

Geometry optimization of the water molecule with a frozen oxygen 2s orbital
---------------------------------------------------------------------------

We run the geometry optmization with the command::

  pam --inp=H2O --mol=H2O --put "cf.O=DFOXXX"

After optimization we find that the bond distance is 0.972 Å and the bond angle 95.9 :math:`^{\circ}`.
This is perhaps a surprising result since one would expect that breaking hybdridization would lead to
a bond angle of 90.0 :math:`^{\circ}`. Our result suggests that there must be some additional repulsion
between the hydrogens of water. The nature of this repulsion (electrostatic vs. Pauli repulsion) has been
a subject of considerable discussion in the literature. A full account is found in :cite:`Dubillard2006`.
