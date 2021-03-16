!     # msu    = max number of symmetry inequivalent types of atoms.
!     msup increased to 1000 to allow for large number of PCs (PSB 05/03/20)
!          also change recp_lnk.F line 145 for local definition of msup
!        number of centers after symmetry operations  
      integer    msup
      parameter (msup=203)
!      parameter (msup=1000) ! PSB 05/03/20

!     # mstu   = max number of irreps.
      integer    mstup
      parameter (mstup=8)

!     # mrcru  = max number of members in a contraction set.
      integer    mrcrup
      parameter (mrcrup=50) !v3 20->50

!     # mconu  = max number of primitives in a contraction.
      integer    mconup
      parameter (mconup=50) !v3 24->50

!     # mcu    = max number of symmetry equivalent atoms.
!     mcup increased to 200 to allow for large number of PCs (PSB 05/03/20)
!        number of centers before symmetry operations  
      integer    mcup
      parameter (mcup=50)   !v2 12->50
!      parameter (mcup=200)  ! PSB 05/03/20

!     # mconsu = max number of contraction sets.
      integer    mconsp
      parameter (mconsp=200) !v2 65->100

!     # msfu   = max number of function sets.
      integer    msfup
      parameter (msfup=200)  ! 52->v1(70),v2(100),1000

      integer    mnsfup
      parameter (mnsfup=(msfup*(msfup+1))/2)

!     # mgu    = max number of operators in the
!     #          nuclear interchange group.
      integer    mgup
      parameter (mgup=49)

!     # mccu   = max number of symmetry unique center
!     #          combinations from 4 nuclei.
      integer    mccup
      parameter (mccup=182)
 
! ### from (subroutine seg1mn) ###

!     # kaords = max number of ao reduction sets.
      integer    kaordp
      parameter (kaordp=240)

!     # mgcsu  = max number of symmetry orbital transformation sets.
      integer    mgcsup
      parameter (mgcsup=240)

!     # mru    = max number of irreps in an ao reduction set.
      integer    mrup
      parameter (mrup=100) !v2 24->50

!     # mcru   = max number of core potential expansion functions.
      integer    mcrup
      parameter (mcrup=140)

!     # msfru  = max number of so sets.
!     msfru increased to 4000 max for large component bfs (PSB 05/03/20)
      integer    msfrup
      parameter (msfrup=1000) !312,500
!      parameter (msfrup=4000) ! PSB 05/03/20  

!     # mnrup  = max number of charge dists. from a function set pair.
      integer    mnrup  ! used only to dimension three arrays (1 dimension)
!                         mau(mnrup) in recp1.F & maxics and maxjcs in recp_socfpd.F90
      parameter (mnrup=5000) !576->1000 - PSB(06/19/17) to 5000
!     parameter (mnrup=10000) !576->1000 - PSB(05/04/20) to 10000

