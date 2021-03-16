C **************************************************************
C ***********    Gauss-Legendre integration       **************
C **************************************************************
C
      PARAMETER ( MAX_GL_PT = 16 )
C
      REAL*8 FREQ_6PT(6),FREQ_8PT(8),FREQ_10PT(10),
     &     FREQ_12PT(12),FREQ_14PT(14),FREQ_16PT(16)
      REAL*8 ABSCISSAS_6PT(6),ABSCISSAS_8PT(8),ABSCISSAS_10PT(10),
     &     ABSCISSAS_12PT(12),ABSCISSAS_14PT(14),ABSCISSAS_16PT(16)
      REAL*8 WEIGHTS_6PT(6),WEIGHTS_8PT(8),WEIGHTS_10PT(10),
     &     WEIGHTS_12PT(12),WEIGHTS_14PT(14),WEIGHTS_16PT(16)
C
      COMMON /GAUSS_LEGENDRE_PARAMETERS/ 
     &     FREQ_6PT,FREQ_8PT,FREQ_10PT,FREQ_12PT,FREQ_14PT,FREQ_16PT,
     &     ABSCISSAS_6PT,ABSCISSAS_8PT,ABSCISSAS_10PT,ABSCISSAS_12PT,
     &     ABSCISSAS_14PT,ABSCISSAS_16PT,
     &     WEIGHTS_6PT,WEIGHTS_8PT,WEIGHTS_10PT,WEIGHTS_12PT,
     &     WEIGHTS_14PT,WEIGHTS_16PT
C
C End of gauss_legendre.h
C

