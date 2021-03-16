:orphan:

Combining other methods with ECP
================================

In DIRAC program, ECP can be incorporated with several ground and excited state calculation methods.
The calculation methods can be set in input file. Below is the few examples of hydrogen iodide calculations.
(See :ref:`ecp_input` for the molecular input) 

DFT calculation
---------------

See the quick Bi2 molecule test (:download:`DFT.inp <../../../../test/ecp/DFT.inp>` and :download:`Bi2.xyz <../../../../test/ecp/Bi2.xyz>`)

.. literalinclude:: ../../../../test/ecp/DFT.inp

COSCI calculation
-----------------

.. literalinclude:: COSCI.inp
   
MP2 and CC calculation
----------------------

.. literalinclude:: CC.inp

