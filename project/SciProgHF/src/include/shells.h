!*****************************************************************************
!     Data on the AO-block associated with a given shell(ISHELL):
!       NUCO    - number of uncontracted functions
!       NRCO    - number of contracted functions
!     Data on a given shell in an AO-block:
!       NCENT   - index of symmetry independent center
!       NUMCF   - index of shell in AO-block
!       NBCH    - index of block in AO-vector
!       LAMN    - name of center
!       ISTBAO  - stabiliser: basic sym. op. that do not move center
!       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
!       KHKT    - number of spherical (Cartesian) components
!       KCKT    - number of Cartesian components
!       SEGM    - segmented contraction
!       LCLASS  - class: large component (1), small or Huckel(2), density fitting (0)
!       CENT    - coordinates of center
!       JSTRT is an offset for the shell
!       MBSID and MBIDBT have been added for multiple basis sets (WK/UniKA/31-10-2002).
!*****************************************************************************
      LOGICAL SHARE(MXSHEL), SEGM(MXSHEL), SPHR(MXSHEL)
      REAL*8 CENT(MXSHEL, 3, 2)
      INTEGER NHKT(MXSHEL), KHKT(MXSHEL), KCKT(MXSHEL), ISTBAO(MXSHEL), &
     &        NUCO(MXSHEL), JSTRT(MXSHEL), NSTRT(MXSHEL), NCENT(MXSHEL),&
     &        NRCO(MXSHEL), NUMCF(MXSHEL), NBCH(MXSHEL), KSTRT(MXSHEL), &
     &        LCLASS(MXSHEL), IPTSHL(MXSHEL), NUMCFT(MXSHEL), KMAX,     &
     &        NLRGSH, NSMLSH, NORBS,                                    &
     &        MST(MXSHEL), MBSID(MXSHEL), MBIDBT(MXSHEL)
      COMMON /SHELLS/ CENT, NHKT, KHKT, KCKT, ISTBAO, NUCO, JSTRT,      &
     &                NSTRT, NCENT, SHARE, NRCO, NUMCF, NBCH, KSTRT,    &
     &                SEGM, LCLASS, IPTSHL, NUMCFT, SPHR, KMAX, NLRGSH, &
     &                NSMLSH, NORBS,                                    &
     &                MST, MBSID, MBIDBT
