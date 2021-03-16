C-----For the ADC module

      COMMON/CCADCI/ADCLEVEL,ADCPRINT,SIPREPS,DIPREPS
      INTEGER ADCLEVEL,ADCPRINT(8),SIPREPS(32),DIPREPS(32)

      COMMON/CCADCR/VCONV,ADCTHR
      REAL*8 VCONV,ADCTHR

      COMMON/CCADCL/ DOSIPS, DODIPS, READQKL, 
     &               DOCONST, DOADCPOP,
     &               DOFULL, DOLANC

      LOGICAL        DOSIPS, DODIPS, READQKL, 
     &               DOCONST, DOADCPOP,
     &               DOFULL, DOLANC

C-----For the LANCZOS diagonalizer 

      COMMON/LANCZI/SIPITER,DIPITER,
     &              SIPEIGV,DIPEIGV,
     &              LANCMEM
      INTEGER SIPITER,DIPITER,
     &        SIPEIGV(32),DIPEIGV(32),
     &        LANCMEM
      COMMON/LANCZR/SIPPRNT,DIPPRNT
      REAL*8 SIPPRNT,DIPPRNT
      COMMON/LANCZL/DOINCORE
      LOGICAL DOINCORE

C-----For the DAVIDSON diagonalizer data are obtained from the module directly

