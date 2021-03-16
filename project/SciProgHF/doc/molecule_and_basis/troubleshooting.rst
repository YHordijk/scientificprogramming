:orphan:
 

GETPOT: Nuclei too close
------------------------

If you see this error then one of the following situations is likely:
 * You have two or more atoms on top of each other in your mol or xyz file.
 * You specify symmetry operations explicitly and have too many atoms (remove symmetry-generated centers).
 * You use .OPTIMIZE together with xyz format (use mol format instead).
