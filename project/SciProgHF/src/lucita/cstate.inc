!comdeck cstate.inc $Revision$ 
!. CSTATE

      REAL(8) ETHRES,CTHRES,EROOT
      INTEGER NROOT,IROOT,IDC,PSSIGN,PLSIGN,INTSEL,IALIGN,              &
     &        IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS

      COMMON/CSTATE/IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,               &
     &              NROOT,IROOT(100),IDC,PSSIGN,PLSIGN,INTSEL,IALIGN,   &
     &              ETHRES,CTHRES,EROOT(100)
      
