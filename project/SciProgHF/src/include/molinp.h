!  FILE: molinp.h
!
!  input .mol file is saved in MLINE(1:NMLINE)
!
      INTEGER KMLINE, LEN_MLINE
      PARAMETER (KMLINE = 1000, LEN_MLINE = 80)

      CHARACTER*(LEN_MLINE) MLINE(KMLINE)
      COMMON /MOLINP_C/ MLINE

      INTEGER NMLINE, NCLINE(MXCENT), NMLAU
      COMMON /MOLINP/ NMLINE, NCLINE, NMLAU
! -- end of molinp.h --
