      parameter (Lmax=4)   ! max. angular momentum of basis functions  DO NOT INCREASE
cbs                          TO MORE THAN four  !!!!!!!!!!!!!!!!!!!
cbs                          if you do, you will have to edit the ixyzpow array by hand...............
cbs                          (in datapow.h)
cbs  you will also have to add the correct normalizations when using 
CBS  the old DALTON Normalisation..
      parameter (Lmax_occ=3)   ! highest L-value for occupied orbitals                     
      parameter (MxprimL=38) ! max. of primitives per angular momentum
      parameter (MxcontL=38) ! max. of contracted functions per angular momentum
      parameter (MxCart=520) ! max. number of contracted functions in the atom

      parameter (ndfmx=4*Lmax+4) ! dimension of precomputed double factorials     
