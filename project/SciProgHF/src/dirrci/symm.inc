C
C SOURCE  : $Source$
C
      PARAMETER (MXREP=32)
      COMMON/RAS/NORB,NORBR(3),NELEC,MINH1,MAXH1,MAXE3,IORBR(0:N3)
      COMMON/ABEL/NREP,IRA1,IREP,IRPAMO(N3),
     &            MULTB(MXREP*2,MXREP*2),MULTBI(MXREP*2,MXREP*2)
      COMMON/RECIND/IGENEX,NTVERT,NCOREX,LNORBC,NORBC,NORBD
C
      COMMON/SPSYMI/NSYMRP,IRPMO(N3)
      COMMON/SPSYMR/EPS(N3)
      COMMON/SPSYMC/REPN(16)
      COMMON/CHR/REPNA(MXREP*2)
      CHARACTER*4 REPNA,REPN*14