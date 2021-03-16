:orphan:
 

Troubleshooting
===============

WARNING: error in the number of electrons
-----------------------------------------

Sometimes you may see the following warning::

   number of grid points                          =     54836
   DFT exchange-correlation energy:               =      -230.0368431371773283
   number of electrons from numerical integration =       107.9973057254599667
   number of electrons from orbital occupations   =       108
   WARNING: error in the number of electrons      =        -0.0026942745400333
            is larger than 1.0d-3
            this can happen when starting from coefficients
            from a different geometry
            or it can mean that the quadrature grid is inappropriate

What does this mean?

When you run a Kohn--Sham (KS) density functional theory (DFT) calculation, the
exchange-correlation (XC) energy and XC potential contributions are integrated
numerically. This means that for KS DFT calculations the numerical quadrature
grid enters as an additional parameter which you control.

In every iteration, when integrating the XC contributions, the code will also
integrate the density at the same time and compare the number of electrons from
numerical integration with the number of electrons from orbital occupations as
a sanity check and issue the above warning if the difference is larger than a
threshold (here 0.001).

If you see this warning, then the number of electrons from numerical
integration deviates significantly from the target number of electrons.  This
is a sanity check and it does not necessarily mean that your results will be
useless, since the integrated number of electrons is not a sufficient measure
for the quality of the grid for the property that you calculate.  If you see
this warning in every iteration of your calculation, then it may be a good idea
to investigate the quality of the grid.  You can either select a finer (and
more expensive) grid (see :ref:`**GRID`), and/or tighten the screening
threshold (:ref:`DFT_.SCREENING` under :ref:`*DFT`).

When restarting from a ``DFCOEF`` file from a different geometry, the warning
can appear during the first few iterations and then this warning disappears
(because the density of the old geometry is then displaced with respect to the
grid at new geometry).

If the integrated number of electrons is significantly off, then there is a
serious problem, possibly a bug. In this case please contact the developers and
send your output.

If you do not see the above warning, then it does not guarantee that your
integration grid is useful for the property that you study. It is always a good
practice to verify that your results are converged with respect to the grid
parameters. When you calculate a table of functionals, it is a good practice
for one of the results to re-run with the best integration grid and to check
that your results do not change.

DFT nonconverging...
--------------------

What to do if your DFT run is not converging ? Instead of DFT try first to perform  the SCF method.
If the previous SCF step does converge, save the MO-coeffient file and use them as starting set for
your DFT method. 

This trick can help because the HOMO-LUMO gap is higher in the SCF than in the DFT. Restarting the DFT
from converged SCF coefficients makes the DFT method converging in many cases.
