!     Inclusion of this common block requires also
!        maxaqn.h
!        maxorb.h
!        mxcent.h

      INTEGER I1_SYMMTI, MAXREP, MAXOPR, MULT(0:7), ISYMAX(3, 2),       &
     &        ISYMAO(MXQN, MXAQN), NPARSU(8), NAOS(8), NPARNU(8, 8),    &
     &        IPTSYM(MXCORB, 0:7), IPTCNT(3*MXCENT, 0:7, 2),            &
     &        NCRREP(0:7, 2), IPTCOR(3*MXCENT, 2), NAXREP(0:7, 2),      &
     &        IPTAX(3, 2), IPTXYZ(3, 0:7, 2), IPTNUC(MXCENT, 0:7),      &
     &        ISOP(0:7), NROTS, NINVC, NREFL, IXVAL(0:7, 0:7),          &
     &        NCOS(8, 2), ICLASS(MXCORB), ICNTAO(MXCORB), IPAR(0:7),    &
     &        JSOP(0:7), ICOS(8, 2), I2COSX(8, 8, 2, 2), I2_SYMMTI

      COMMON /SYMMTI/ I1_SYMMTI, MAXREP, MAXOPR, MULT, ISYMAX, ISYMAO,  &
     &                NPARSU, NAOS, NPARNU, IPTSYM, IPTCNT, NCRREP,     &
     &                IPTCOR, NAXREP, IPTAX, IPTXYZ, IPTNUC, ISOP,      &
     &                NROTS, NINVC, NREFL, IXVAL, NCOS, ICLASS, ICNTAO, &
     &                IPAR, JSOP, ICOS, I2COSX, I2_SYMMTI

!     Do NOT add any variable after the end tag
      INTEGER NSYMMTR
      PARAMETER (NSYMMTR = 16)

      REAL*8 FMULT(0:7), PT(0:7)
      COMMON /SYMMTR/ FMULT, PT
