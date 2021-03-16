!C FILE: moltra/dcbtr3.h
!C
!C*** Symmetry packing information for 4-index transformation program **
!C
!C     (NSPCK and ISPCK in common block COMDIS)
!C
!C     1, 2, 3, 4 in variable names refer to index 1, 2, 3, 4, respectively,
!C        in the integral transformation.
!C     N* variables provide number of element for each index - in principle
!C        this makes it possible to have different basis sets for different indices.
!C     I* variables provide off-set to the corresponding N* elements in the relevant index
!C
!C     NDMOQC : used for density matrix to Fock core matrix
!C              (in fact, only NDMOQC(2,:,1) used in current code)
!C     NDMOQR : 1st dimension 1:2: row/column dimension of MO coefficient array
!C                                 (e.g. NFBAS(IFRP,0), NASH(IFRP) )
!C              2nd dimension 1:2: fermion irrep IFRP
!C              3rd dimension 1:4: 2-el transformation index
!C     NDMOQS and NDMOQT: as NDMOQR but for reduced coefficient spaces for
!C              index 3 and 4 when more than one pass is needed (see moltra/tradr4.F)
!C
      COMMON/DCBTR3/NSPCK12(0:7,0:3),  NSPCK34(0:7,0:3),                &
     &              ISPCK12(0:7,0:7,3),ISPCK34(0:7,0:7,3),              &
     &              NBBAS1(0:7,0:2),NBBAS2(0:7,0:2),                    &
     &              NBBAS3(0:7,0:2),NBBAS4(0:7,0:2),                    &
     &              IBBAS1(0:7,0:2),IBBAS2(0:7,0:2),                    &
     &              IBBAS3(0:7,0:2),IBBAS4(0:7,0:2),                    &
     &              IBAS1(2),IBAS2(2),IBAS3(2),IBAS4(2),                &
     &              NDMOQC(2,2,2),ICMOQC(2,2),                          &
     &              NDMOQR(2,2,4),ICMOQR(2,4),                          &
     &              NDMOQS(2,2,4),ICMOQS(2,4),                          &
     &              NDMOQT(2,2,4),ICMOQT(2,4)                 
!C -- end ofmoltra/dcbtr3.h
