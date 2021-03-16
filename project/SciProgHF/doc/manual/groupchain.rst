:orphan:

Symmetry-handling at the correlated level
=========================================

Introduction
------------

The DIRAC code can handle symmetries corresponding to :math:`D_{2h}` and subgroups (denoted binary groups) as well as linear supersymmetry.
At the SCF level DIRAC employs a unique quaternion symmetry scheme which combines time reversal and spatial symmetry (:cite:`Saue1999`).
A particularity of this scheme is that symmetry reductions due to spatial symmetry is translated into a reduction of algebra, from
quaternion down to complex and possibly real algebra. This leads to a classification of the binary groups as:

- Quaterion groups: :math:`C_1` , :math:`C_i`
- Complex groups:   :math:`C_2` , :math:`C_s` , :math:`C_{2h}`
- Real group:       :math:`D_2` , :math:`C_{2v}` , :math:`D_{2h}`

In general, at the correlated level, the highest Abelian subgroup of the point group under consideration is used.
All quaternion and complex subgroups are Abelian.
The real groups are non-Abelian and so an Abelian subgroup is chosen


.. math::  \begin{array}{lcl}D_2, C_{2v}&\rightarrow &C_2\\D_{2h}&\rightarrow &C_{2h}\end{array}

Character tables for Abelian double groups
------------------------------------------

Below we give character tables for the Abelian double groups handled by DIRAC. The irreps are given in the same order as internally in the correlation modules
of the program. This means that fermion irreps come before boson irreps. 
Otherwise, the tables are given according to the conventions of S. L. Altmann and P. Herzig, *Point-Group Theory Tables*, Clarendon Press, Oxford, 1994 
(A second corrected edition is now available free of charge `here <http://phaidra.univie.ac.at/o:104731>`_.). The final column of the tables give the irrep labels employed by DIRAC. 

- Real groups :

+----------------------+-----------+-----------+
| :math:`\mathbf{C_1}` | :math:`E` |           |
+----------------------+-----------+-----------+
| :math:`A_{1/2}`      | :math:`1` | :math:`A` |
+----------------------+-----------+-----------+
| :math:`A`            | :math:`1` | :math:`a` |
+----------------------+-----------+-----------+

+----------------------+-----------+----------------------+------------+
| :math:`\mathbf{C_i}` | :math:`E` | :math:`i`            |            |
+----------------------+-----------+----------------------+------------+
| :math:`A_{1/2,g}`    | :math:`1` | :math:`\phantom{-}1` | :math:`AG` |
+----------------------+-----------+----------------------+------------+
| :math:`A_{1/2,u}`    | :math:`1` |:math:`-1`            | :math:`AU` |
+----------------------+-----------+----------------------+------------+
| :math:`A_g`          | :math:`1` | :math:`\phantom{-}1` | :math:`ag` |
+----------------------+-----------+----------------------+------------+
| :math:`A_u`          | :math:`1` | :math:`-1`           | :math:`au` |
+----------------------+-----------+----------------------+------------+

- Complex groups:

+----------------------+-----------+----------------------+------------+
| :math:`\mathbf{C_2}` | :math:`E` | :math:`C_2`          |            |
+----------------------+-----------+----------------------+------------+
| :math:`\,^1E_{1/2}`  | :math:`1` | :math:`\phantom{-}i` | :math:`1E` |
+----------------------+-----------+----------------------+------------+
| :math:`\,^2E_{1/2}`  | :math:`1` |:math:`-i`            | :math:`2E` |
+----------------------+-----------+----------------------+------------+
| :math:`A`            | :math:`1` | :math:`\phantom{-}1` | :math:`a`  |
+----------------------+-----------+----------------------+------------+
| :math:`B`            | :math:`1` | :math:`-1`           | :math:`b`  |
+----------------------+-----------+----------------------+------------+

+----------------------------+-----------+----------------------+------------+
| :math:`\mathbf{C_s}`       | :math:`E` | :math:`\sigma_h`     |            |
+----------------------------+-----------+----------------------+------------+
| :math:`\,^1E_{1/2}`        | :math:`1` | :math:`\phantom{-}i` | :math:`1E` |
+----------------------------+-----------+----------------------+------------+
| :math:`\,^2E_{1/2}`        | :math:`1` |:math:`-i`            | :math:`2E` |
+----------------------------+-----------+----------------------+------------+
| :math:`A^{\prime}`         | :math:`1` | :math:`\phantom{-}1` | :math:`a`  |
+----------------------------+-----------+----------------------+------------+
| :math:`A^{\prime\prime}`   | :math:`1` | :math:`-1`           | :math:`b`  |
+----------------------------+-----------+----------------------+------------+


+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`\mathbf{C_{2h}}` | :math:`E` | :math:`C_2`          | :math:`i`            | :math:`\sigma_h`     |             |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`\,^1E_{1/2,g}`   | :math:`1` | :math:`\phantom{-}i` | :math:`\phantom{-}1` | :math:`\phantom{-}i` | :math:`1Eg` |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`\,^2E_{1/2,g}`   | :math:`1` | :math:`-i`           | :math:`\phantom{-}1` | :math:`-i`           | :math:`2Eg` |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`\,^1E_{1/2,u}`   | :math:`1` | :math:`\phantom{-}i` | :math:`-1`           | :math:`-i`           | :math:`1Eu` |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`\,^2E_{1/2,u}`   | :math:`1` | :math:`-i`           | :math:`-1`           | :math:`\phantom{-}i` | :math:`2Eu` |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`A_g`             | :math:`1` | :math:`\phantom{-}1` | :math:`\phantom{-}1` | :math:`\phantom{-}1` | :math:`ag`  |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`B_g`             | :math:`1` | :math:`-1`           | :math:`\phantom{-}1` | :math:`-1`           | :math:`bg`  |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`A_u`             | :math:`1` | :math:`\phantom{-}1` | :math:`-1`           | :math:`-1`           | :math:`au`  |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+
| :math:`B_u`             | :math:`1` | :math:`-1`           | :math:`-1`           | :math:`\phantom{-}1` | :math:`bu`  |
+-------------------------+-----------+----------------------+----------------------+----------------------+-------------+


