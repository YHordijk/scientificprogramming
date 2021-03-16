:orphan:
 

star(NMR)

This section gives directives for the calculation of NMR parameters.

If common gauge origin (CGO) is used, i.e. :LONDON not specified,
then the user can define the gauge origin used for the external
magnetic field with :ref:`HAMILTONIAN_.GAUGEORIGIN` or :ref:`HAMILTONIAN_.GO ANG`
under :ref:`**HAMILTONIAN`.
If :ref:`NMR_.USECM` is specified then center-of-mass is used as gauge origin.
Default is (0, 0, 0).


**Advanced options**

keyword(LONDON)

Activate calculations of magnetic properties (NMR shielding constants
and magnetizabilities) with London atomic orbitals.

*Default:* Use conventional atomic orbitals.

keyword(USECM)

Use the center of mass as the gauge origin.

keyword(INTFLG)

Specify what two-electron integrals to include in the two-electron
London contributions to the magnetic field property gradient
(default: :ref:`HAMILTONIAN_.INTFLG` under :ref:`**HAMILTONIAN`).

keyword(NOTWO)

Do not calculate the two-electron London contributions for the magnetic
field property gradient when London atomic orbitals are used.

keyword(NOONEI)

Do not calculate the {H(0),T(B)} reorthonormalization terms for the
magnetic field property gradient when London atomic orbitals are used.

keyword(NOORTH)

Do not calculate the {T(B),h(mK)} reorthonormalization contributions for
the expectation value term when London atomic orbitals are used.

keyword(SYMCON)

Employ the symmetric connection for reorthonormalization terms when
using London atomic orbitals.

*Default:* Use the natural connection.

keyword(EXPPED)

Keyword used in the DFT calculations only. The contributions to the density perturbed by an external magnetic field
in LAO basis ("direct" LAO term and "reorthonormalization" term) are exported on files, 
`pertden_direct_lao.FINAL` and `pertden_reorth_lao.FINAL`.

