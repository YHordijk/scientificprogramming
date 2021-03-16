:orphan:
 

star(PCM)

Polarizable continuum model (PCM) directives. 

General information on the PCM can be found in :cite:`Tomasi2005`, whereas details 
concerning the Dirac implementation are found in :cite:`DiRemigio2015`.

.. warning:: Calculations with the polarizable continuum model **cannot** exploit molecular point group symmetry and/or 2-component Hamiltonians.

**General**
===========

keyword(SKIPSS)

Performs a PCM-SCF calculation with the PCM-SCC approximation. The small-small block
of the electrostatic potential is neglected, in the spirit of the SCC approximation 
described in :cite:`Visscher1997a`. It is not employed by default.

.. warning:: Currently not working for response calculations.

keyword(PRINT )

Print level in the PCM subroutines. Default::

  .PRINT
   0

**Advanced/Debug**
==================

keyword(SEPARA)

Performs PCM-SCF separating the nuclear and electronic electrostatic potential and apparent surface
charge. Results are unaffected. *For debug purpose only*

keyword(DOSPF )

Remove spin-orbit dependent part from PCM potential. *For debug purpose only*

keyword(SKIPOI)

Skips the calculation of the one-index transformed apparent surface charge in the
linear response function. *For debug purpose only* 

