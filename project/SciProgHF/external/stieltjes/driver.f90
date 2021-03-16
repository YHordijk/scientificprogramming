PROGRAM stieltjes_driver

  use StieltjesMod
  use input_routines

  IMPLICIT NONE

!
!  Input variables
!
  CHARACTER(len=100) :: input_file
  INTEGER            :: nr_points

!
!  Data dictionary: Pseudo-spectrum
!
  REAL(8)                           :: e_init
  REAL(8),ALLOCATABLE,DIMENSION(:)  :: energies
  REAL(8),ALLOCATABLE,DIMENSION(:)  :: points

!
!  File handling
!
  INTEGER                           :: fh_moments=2004
  CHARACTER(11)                     :: outfile="inp.dat"

!
!  Results
!
  REAL(8)                           :: gammae

!
!  Code
!
CALL get_command_argument(1, input_file)

nr_points = len_file(input_file) - 1

ALLOCATE(energies(nr_points))
ALLOCATE(points(nr_points))

CALL read_pseudo_spec(input_file,e_init,nr_points,energies,points)

CALL stieltjes(nr_points,energies,points,e_init,outfile,gammae,fh_moments,2)

DEALLOCATE(energies)
DEALLOCATE(points)



END PROGRAM stieltjes_driver
