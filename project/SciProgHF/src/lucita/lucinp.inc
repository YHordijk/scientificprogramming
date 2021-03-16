!comdeck lucinp.inc $Revision$ 
!./LUCINP/
      INTEGER PNTGRP,EXTSPC


      INTEGER NIRREP,NSMCMP,MAXML,MAXL,                                 &
     &              INTSPC,NRSSH,                                       &
     &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,                 &
     &              NSMOB,NRS0SH,NRS4SH,                                &
     &              MXR4TP,MXHR0,MXER4,                                 &
     &              NINASH,                                             &
     &              INTXCI,NDELSH,MNRS10,MXRS30,                        &
     &              MNRS1RE,MXRS1RE,MNRS3RE,MXRS3RE,                    &
     &              MNRS1ZE,MXRS1ZE,MNRS3ZE,MXRS3ZE,                    &
     &              MPORENP,MPORENR,MPORENZ

      COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,                    &
     &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),                      &
     &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,                 &
     &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPIRR,MXPR4T),       &
     &              MXR4TP, MXHR0,MXER4,                                &
     &              NINASH(MXPIRR),                                     &
     &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30,                &
     &              MNRS1RE,MXRS1RE,MNRS3RE,MXRS3RE,                    &
     &              MNRS1ZE,MXRS1ZE,MNRS3ZE,MXRS3ZE,                    &
     &              MPORENP,MPORENR,MPORENZ
