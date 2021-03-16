MODULE input_routines

CONTAINS

SUBROUTINE read_pseudo_spec(filename,e_init,length,energies,points)

  IMPLICIT NONE
! Data dictionary
  CHARACTER(len=100), INTENT(IN) :: filename
  INTEGER, INTENT(IN)            :: length
  REAL(8), INTENT(OUT)           :: e_init
  CHARACTER(len=200) :: buffer
  CHARACTER(len=1) :: dummy
  INTEGER :: pos
  INTEGER, PARAMETER :: fh = 17
  INTEGER :: ierror = 0
  INTEGER :: line = 0
  INTEGER :: ixx
  INTEGER :: input ! counter for input parameters
  REAL(8), DIMENSION(length), INTENT(OUT) :: energies
  REAL(8), DIMENSION(length), INTENT(OUT) :: points
  REAL(8), DIMENSION(length,2)            :: temp

  OPEN(fh, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
!  WRITE(*,*) 'file opening: ', ierror

  READ(fh,*, IOSTAT=ierror) e_init
  WRITE(*,*) 'e_init = ', e_init
  DO WHILE (ierror == 0)
    READ(fh, '(A)', IOSTAT=ierror) buffer

    IF (ierror == 0) THEN
! Find first white space
! If the first character is #, then don't consider this line
      pos = scan(buffer, '   ')
      dummy = TRIM(buffer(1:pos))
      buffer = buffer(pos+1:)

      IF (LGE(dummy,'#').AND.LLE(dummy,'#')) THEN
      ELSE
        line = line + 1
!        WRITE(*,*) 'line: ', line
        READ(buffer, *, IOSTAT=ierror) (temp(line,input), input=1,2)
!        WRITE(*,*) 'Channel ', (channels(line,input), input=1,2)
      END IF
    END IF
  END DO

  energies(1:length) = temp(1:length,1)
  points(1:length)   = temp(1:length,2)

  DO ixx = 1, length
    WRITE(*,*) energies(ixx), points(ixx)
  END DO

  CLOSE(fh)

END SUBROUTINE read_pseudo_spec




INTEGER FUNCTION len_file(filename)
! Purpose: determine the number of non commented lines in a given file
 IMPLICIT NONE
! Data dictionary
 CHARACTER(len=100), INTENT(IN) :: filename
 CHARACTER(len=100) :: buffer
 CHARACTER(len=1) :: dummy
 INTEGER :: pos
 INTEGER, PARAMETER :: fh = 16
 INTEGER :: ierror = 0
 INTEGER :: line

 OPEN(fh, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
! WRITE(*,*) 'ierror: ', ierror

 line = 0

 DO WHILE (ierror == 0)
   READ(fh, '(A)', IOSTAT=ierror) buffer
   IF (ierror == 0) THEN
! Find first white space
! If the first character is #, then don't consider this line
     pos = scan(buffer, '    ')
     dummy = TRIM(buffer(1:pos))
     
     IF (LGE(dummy,'#').AND.LLE(dummy,'#')) THEN
     ELSE
       line = line + 1
     END IF
   END IF
 END DO

 CLOSE(fh)

 len_file = line

END FUNCTION len_file

END MODULE input_routines
