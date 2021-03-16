! FILE    : include/dcbidx.h
!
! Needs: maxorb.h, maxash.h
!
! PURPOSE: index arrays for the krmc module (and prp module when krmc)
!
! Originally designed by HJAaJ 2001 based on similar data structures
! in SIRIUS in DALTON.
!
!     IFSMO : Fermion symmetry of orbital (1 for gerade, 2 for ungerade)
!     IOBTYP: type of orbital:
!        - frozen orbital         : 0
!        - positronic orbital     : 1
!        - inactive orbital       : 2
!        - active orbital         : 3
!        - secondary orbital      : 4
!
!     IDXa2b index vector gives index of "a" to "b" renumbering
!     where "a" and "b" are one of:
!       G : general orbital
!       U : active orbital
!       E : electronic orbital ("positive energy" orbital)
!       P : positronic orbital ("negative energy" orbital)
!       I : inactive orbital (doubly occupied orbital)
!       F : frozen inactive orbital
!       S : secondary orbital (empty electronic orbital)
!       V : active secondary orbital (non-deleted secondary orbitals)
!       D : deleted empty orbitals
!     IDXT2G(:,1:4) :
!       general orbital index of each of the four "transformation" orbitals
!       for the four indices in the integral transformation (1 2 | 3 4)
!     NIDX*, IIDX* :
!       info related to the four indices in the integral transformation
!
!     
      INTEGER JTFRO, JTPOSI, JTINAC, JTACT, JTSEC
      PARAMETER (JTFRO=0, JTPOSI=1, JTINAC=2, JTACT=3, JTSEC=4)

      INTEGER IFSMO(MXCORB) , IOBTYP(MXCORB), IDXG2U(MXCORB),           &
     &        IDXU2G(MAXASH), IDXE2G(MXCORB), IDXG2E(MXCORB),           &
     &        IDXP2G(MXCORB), IDXG2P(MXCORB), IDXI2G(MXCORB),           &
     &        IDXG2I(MXCORB), IDXS2G(MXCORB), IDXG2S(MXCORB),           &
     &        IDXV2G(MXCORB), IDXG2V(MXCORB), IDXG2D(MXCORB),           &
     &        IDXD2G(MXCORB), IDXT2G(MXCORB, 4),                        &
     &        NIDXT(4), IIDX1(2), NIDX1(2), IIDX2(2), NIDX2(2),         &
     &                  IIDX3(2), NIDX3(2), IIDX4(2), NIDX4(2)
      COMMON /DCBIDX/ IFSMO,  IOBTYP, IDXG2U, IDXU2G, IDXE2G, IDXG2E,   &
     &                IDXP2G, IDXG2P, IDXI2G, IDXG2I, IDXS2G, IDXG2S,   &
     &                IDXV2G, IDXG2V, IDXG2D, IDXD2G, IDXT2G,           &
     &                NIDXT, IIDX1, NIDX1, IIDX2, NIDX2,                &
     &                       IIDX3, NIDX3, IIDX4, NIDX4
! end of include/dcbidx.h
