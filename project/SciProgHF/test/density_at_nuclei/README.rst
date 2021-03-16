Electronic density at nuclei
============================

This test will do a LDA-DFT calculation of carbon monoxide at the experimental
bond length and calculates the electronic density at nuclei using the
expectation value module.

Two different schemes are employed: i) Fermi-contact (FC) operator (.RHONUC),
and ii) effective density using the PVC operator (.EFFDENS).

We assume that the basis sets are available in your standard path.
These sets can be obtained from the Dalton distribution and may later be included in the DIRAC distribution as well.

Note that a coarse DFT grid is used and the calculation is stopped after 5 iterations.
