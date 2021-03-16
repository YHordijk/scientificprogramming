:orphan:
 

starstar(GRID)

The numerical integration scheme uses the Becke partitioning :cite:`Becke1988a`
of the molecular volume
into atomic ones, for which the quadrature is performed in spherical
coordinates. Radial integration is carried out using the scheme proposed
by Lindh *et al.* :cite:`Lindh2001`,
while the angular integration is handled by a set of highly accurate
Lebedev grids. Note that the Becke partitioning scheme generally uses
Slater-Bragg radii for atomic size adjustments. For the heavier
elements, with their more variable oxidation states, this can lead to
errors. If the user invokes the :ref:`GRID_.ATSIZE` keyword, DIRAC
tries to deduce relative atomic sizes from available densities, if it
can find coefficients.

The  **radial integration** employs an exponential grid

.. math::

   r_k=e^{kh};\quad w_k = r_k;\quad S=\sym_{k=-\infty}^{\infty}=w_kf\left(r_k\right)

and is specified in terms of step size :math:`h`, the inner point :math:`r_1` 
and the outermost point :math:`r_k_H`, chosen to provide the required relative 
precision :math:`R` (equal to the discretization error :math:`R_D`. 

keyword(RADINT)

Specify the maximum error in the radial integration.

*Default:*

::

    .RADINT
     1.0D-13

keyword(ANGINT)

Specify the precision of the Lebedev angular grid. The angular
integration of spherical harmonics will be exact to the L-value given by
the user. By default the grid will be pruned. The highest implemented
value is 64.

*Default:*

::

    .ANGINT
     41

keyword(NOPRUN)

Turn off the pruning of the angular grid.

keyword(ANGMIN)

Specify the minimum precision of the Lebedev angular grid after pruning.

::

    .ANGMIN
     LEBMIN

Close to the nucleus, the precision of the angular grid will be less
than L-value given by :ref:`GRID_.ANGINT`, but it will not be less
than the integer LEBMIN. Note that giving a LEBMIN value greater than or
equal to the L-value given by :ref:`GRID_.ANGINT` is equivalent to
turning the pruning off.

*Default:*

::

    .ANGMIN
     15

keyword(ATSIZE)

Generate new estimates for atomic size ratios for use in the Becke
partitioning scheme. Relative sizes for a pair of atoms A and B are
calculated from their contribution to the large component density along
the line connecting A and B.

keyword(IMPORT)

::

    .IMPORT
    numerical_grid

Import previously exported numerical grid.

For debugging you can also create your own grid file. The grid file is
formatted - in free format. A grid file for three points could look like
this:

::

    3
    0.1  0.1 0.1    1.0
    0.01 0.2 0.4    0.9
    9.9  9.9 9.0    0.8
    -1

The first integer is the number of points, then come the x-, y-, and
z-coordinate, and the weight for each point. The last line is a negative
integer. Works also in parallel.

You can export the DIRAC grid simply by copying back the file
"numerical\_grid" from the scratch directory.

keyword(NOZIP)

DIRAC will try to remove redundant grid points when symmetry is present.
Each reflection plane reduces the number of grid points by a factor of
two. With this keyword you can turn this default symmetry grid
compression off.

keyword(4CGRID)

Include the small component basis in the generation of the DFT grid even
if you run a 1- or 2-component calculation.

This can be useful if you want to compare to a 4-component calculation
using the identical integration grid.

By default the small component is ignored in the grid generation when
running 1- or 2-component DFT calculations.

keyword(DEBUG)

Very poor grid - corresponds to

::

    .RADINT
     1.0D-3
    .ANGINT
     10

keyword(COARSE)

Coarse grid - corresponds to

::

    .RADINT
     1.0D-11
    .ANGINT
     35

keyword(ULTRAFINE)

A better grid than default - corresponds to

::

    .RADINT
     2.0D-15
    .ANGINT
     64

keyword(INTCHK)

Test the performance of the grid by computing the overlap matrix
numerically and analyzing the errors. In addition, the error matrix can
be printed.

*Default (no test):*

::

    .INTCHK
     0

*Error analysis:*

::

    .INTCHK
     1

*Print error matrix:*

::

    .INTCHK
     2

