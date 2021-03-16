======================================================
The Nuclear Magnetic Quadrupole Moment (NMQM) operator 
======================================================

The NMQM operator is built of the EFG integrals that are already coded in DIRAC program package.

The NMQM operator is the cross product of the alpha matrix and the EFG integrals.
The EFG integrals by definition contain an extra factor of "3",
which to be taken care in the form of COMFACTOR or at the end with a
proper sign, if you use the .OPERATOR technique.

As you can see from the input, there are YZ and XZ component of EFG because 
one considers the Z-axis as the molecular axis of the MgF testing molecule. 
The ZAVECTOR (i.e. \alpha_z) will
couple Y and X component of \vec{r}. Hence we need to consider YZ and XZ
components and cross multiplied with (\alpha_z). 

See the paper PRA93, 012505_2016, and the 
Eq(s). 13, 17, 18, 19. You can see that the EFG integral gives an 
extra factor of "3" and with an opposite sign. 
