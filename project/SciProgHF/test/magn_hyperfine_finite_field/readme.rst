Cl(2P3/2) atom 
==============

Magnetic hyperfine operator as finite-field perturbation within Coupled-Cluster method.
Used ccpVDZ decontracted basis set; C2 symmetry due to the A_x operator

CCSD(T)
~~~~~~~
The X 2P3/2 state of Cl has with the occupation 1s(2) 2s(2) 2p(6) 3s(2) 3p12(2) 3p32(3).
One has to use orbital energies in RelCC to get CC convergence.

However, this CC treatment numerically unstable due the p32 degeneracy (the same for 2Paver state).
We have to resort to Fock-space CC method.

Fock-space CCSD
~~~~~~~~~~~~~~~

Numerically stable method, starting from closed-shell Cl-, removing one electron.

Energy of states - unperturbed:

 -461.212033602450    0.000000000000000 (   4 * ) ... 2P3/2
 -461.207835968605   -0.000000000000000 (   2 * ) ... 2P1/2

Perturbation is removing degeneracy; however, perturbed FSCC is also numerically unstable !

   1     0.0000000000      0.092613516014     -461.212063508532   -0.000000043449273 (   1 * )
   2     0.0000022743      0.092615790321     -461.212061234225    0.000000007898751 (   1 * )
   3     0.0000171691      0.092630685125     -461.212046339421   -0.000000027993950 (   1 * )
   4     0.0000227558      0.092636271842     -461.212040752704   -0.000000033090040 (   1 * )

   5     0.0042341723      0.096847688277     -461.207829336270   -0.000000050155310 (   1 * )
   6     0.0042345301      0.096848046142     -461.207828978405   -0.000000002909875 (   1 * )

Conclusion
----------

Perturbed CCSD energies are numerically unstable.
