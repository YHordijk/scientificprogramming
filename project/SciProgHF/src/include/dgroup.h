! FILE: dgroup.h
!
!     NZ    - number of fermion matrices (1,2,4)
!     NZC1  - number of fermion matrices in C1 symmetry (1 or 4)
!     NZ_in_CI = min(2,NZ), if 2 then complex CI vectors, if 1 then real CI vectors
!     NFSYM - number of fermion irreps
!     NBSYM - number of boson irreps
!     JSPINR(IM,IC,IFRP) - boson irrep for a given spinor position
!     JBTOF(IBRP,IC) - points to fermion irrep for given boson 
!                      irrep/component
!     JFSYM(NZ,IFRP) - points to boson symmetries associated with 
!                      given fermion ircop
!     IQDEF(4)
!     JQBAS(IBRP,IC) - points to quaternionic position
!     IPQTOQ(4,IBRP) - pointer from packed quaternion to quaternion unit
!     IQTOPQ(4,IBRP) - pointer from quaternion to packed quaternion unit
!     IHQMAT(4,IH, ) - matrix symmetry(symmetric/antisymmetric)
!     IRQMAT(4,IBRP) - irrep of matrix
!     IBTOSF(0:7,2)  - pointer from list of boson irreps to list of symmetries
!                      sorted on fermion irreps
!
      CHARACTER*3 FREP(2)
      COMMON /DNAME/ FREP

      INTEGER NDGROUP
      PARAMETER (NDGROUP = 207)

      INTEGER I1_DGROUP, NZ, NFSYM, NBSYM, IFDIRP(4, 2), JBTOF(0:7, 2), &
     &        JMROI(4), JSPINR(4, 2, 2), JQBAS(0:7, 2), JFSYM(4, 2),    &
     &        IPQTOQ(4, 0:7), IQTOPQ(4, 0:7), IQPH(4, 2),IBTOSF(0:7,2), &
     &        IHQMAT(4,-1:1), IRQMAT(4, 0:7), IQDEF(4), ITOT,           &
     &        NZC1, NZ_in_CI, IQMAP(4), I2_DGROUP, ITERSEOUT
      COMMON /DGROUP/ I1_DGROUP, NZ, NFSYM, NBSYM, IFDIRP, JBTOF, JMROI,&
     &                JSPINR, JQBAS, JFSYM, IPQTOQ, IQTOPQ, IQPH,       &
     &                IHQMAT, IRQMAT, IQDEF, ITOT, IBTOSF,              &
     &                NZC1, NZ_in_CI, IQMAP, I2_DGROUP, ITERSEOUT


!     Do NOT add any variables after the end tag !
      LOGICAL LINEAR, ATOMIC
      COMMON /DGRPL/ LINEAR, ATOMIC

! =========================================================
! ... number of operator types in DIRAC is now in NOPTYP
! These are: 
! 1. P, 2. i[alpha_x]P,3. i[alpha_y]P,...19. beta, 20. iS.P,21 iBETAGAMA5
! 22 (iS.P)d
! =========================================================
      INTEGER NOPTYP,NMATTYP
      PARAMETER (NOPTYP = 22, NMATTYP=16)

! *** 4-component operators : Matrix part ***
!
!     M = I_4,iA_z,iA_y,iA_x
!
!     JM4REP(0:16) - bosonirrep of matrix operator
!     JM4TIM(0:16) - tmie reversal symmetry of matrix operator
!     JM4POS(0:16) - quaternionic position of matrix operator 
!                   (1,2(i),3(j),4(k))
!
!     OPERATORS (NOPTYP types):
!
!     MCMP  - number of components for a given operator
!     JM4   - pointer to matrix operator for a given component
!     JCOM  - coefficient of component
!
      INTEGER MCMP(NOPTYP), JM4(3,NOPTYP), JCOM(3,NOPTYP),              &
     &        JM4REP(0:NMATTYP),JM4POS(0:NMATTYP),JM4TIM(0:NMATTYP)
      CHARACTER*4 M4COMB
      COMMON /ONEOPS/ MCMP, JM4, JCOM, JM4REP, JM4POS, JM4TIM
      COMMON /ONECMB/ M4COMB(0:NMATTYP)
! -- end of dgroup.h --
