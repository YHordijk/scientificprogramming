!----------------------------------------------------------------------------------------
!
!
!            Common block for the relscf program. Note that it can not be used
!  in the DIRAC due to overlapping with its variables.
!
!
! Some parameters:
!
!          MXXTRP - maximum number of iterations
!
!
!----------------------------------------------------------------------------------------
      COMMON /DIMI/ NSIZE(10),NOSH(4),NCSH(4),NBCON(4),NBAS(4),           
     & NSAVF(4),N1(6)                                                   
      COMMON /DIMR/ FACTO(16),FATT(21),AJMN(24),OCCUCS(4),OCCUP(4),      
     & RADINC(10)

      COMMON /NONDIMI/  
     &  NBAS1,NVAR,NBAST,NBVAR1,NDIAG,NXTRP,MXXTRP,NBLOCK,NSYM,NUMVAR,  
     &  NRUN,NITSCF,NITSC1,NITSC2,NCONV,NALARM,N0,NSTEP,IIII,NSHT,NOSHT,
     &  NPRINT,NBASM,MXVAR,MVAR,NEXTRA,IEX, 
     &  NEWBAS,N1T,NSVD,IDATA,NDIMPQ,NORBIT,MXBF,MXCC,                  
     &  NPVEFA,NPOEFA,NPENFA,NPBAFI,NPCOFI,IPRINT

      COMMON /NONDIMR/ ZN,BIAS,DGATH,SCFAT,EXPMN,DMNXP,RMNED,RMXED, 
     &  ENERG,POT,CIN,VIR,THRSCF,THRSH,ENORM,RED,ZETDIF,CON1,CON2,ENCOR
                       
