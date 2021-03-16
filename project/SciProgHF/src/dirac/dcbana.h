C     wave function analysis module
C     =============================
C
C     DOPOP  - perform Mulliken population analysis
C     DOPRJ  - perform projection onto another solution
C     DOVEC  - print MO-coefficients
C     DORHO1 - density plot in one dimension
C     DO3RHO - write density to unformatted file (Gaussian cube format)
C              for visualization
C     DOLOC  - perform orbital localization
C
C
      LOGICAL         DOPOP, DOPRJ, DOVEC, DO1RHO, DO3RHO, DO1WT,       &
     &                DOLOC
C
      COMMON /DCLPRP/ DOPOP, DOPRJ, DOVEC, DO1RHO, DO3RHO, DO1WT,       &
     &                DOLOC
C -- end of dcbana.h --
