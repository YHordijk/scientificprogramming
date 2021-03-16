!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

!     !!! MODULE FOR HANDLING FILES !!!
!     Contains
!       - OPEN_FILE(filename) or OPEN_FILE(filename, filenum,formatted,output)
!         Opens a file and returns the id. Which is a new number if only the
!         name is given, otherwise it is equal to filenum.
!         A negative or zero indicates an error.
!       - GET_EXTENSION(filename)
!         Returns the first extension in the array filetypes for which
!         filename.extension is an existing filename.
!    
!     JHS 05-07-07

      MODULE FILE
      IMPLICIT NONE

      PRIVATE
      PUBLIC OPEN_FILE, GET_EXTENSION, filter_comments

        CHARACTER(LEN=8), DIMENSION(2), PARAMETER :: filetypes=(/'MOL', &
     &                                                           'XYZ'/)
        INTEGER :: linewidth=0

        INTERFACE OPEN_FILE
          MODULE PROCEDURE OPEN_FILE_NAME, OPEN_FILE_NUM, OPEN_SCRATCH
        END INTERFACE
        
        
      CONTAINS

        INTEGER FUNCTION OPEN_FILE_NAME(filename)
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0  -   File does not exists.
!     -1  -   To many files opened.
!     -99  -   Error while executing the open statement.
          CHARACTER(LEN=*), INTENT(IN) :: filename
          INTEGER, parameter :: maxfile=99
          LOGICAL :: ex
          INTEGER :: n

          INQUIRE(file=filename, exist=ex)
          IF (ex) THEN
            DO n=1,maxfile
              INQUIRE(n, opened=ex)
              IF (.NOT.(ex)) THEN
                OPEN(n, file=filename, err=563)   
                EXIT
              END IF
            END DO
 563        open_file_name=-99
            IF (n>maxfile) THEN
              open_file_name=-1 ! To many files opened
            ELSE
               open_file_name=n ! Open file
            END IF                  
          ELSE
            open_file_name=0 ! Filename does not exist
          END IF

        END FUNCTION OPEN_FILE_NAME

        INTEGER FUNCTION OPEN_SCRATCH()
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0  -   File does not exists.
!     -1  -   To many files opened.
!     -99  -   Error while executing the open statement.
          INTEGER, parameter :: maxfile=99
          LOGICAL :: ex
          INTEGER :: n

          DO n=1,maxfile
            INQUIRE(n, opened=ex)
            IF (.NOT.(ex)) THEN
              OPEN(n, status='SCRATCH')
              EXIT
            END IF
          END DO
 563      open_scratch=-99
          IF (n>maxfile) THEN
            open_scratch=-1 ! To many files opened
          ELSE
             open_scratch=n ! Open file
          END IF                  

        END FUNCTION OPEN_SCRATCH

        INTEGER FUNCTION OPEN_FILE_NUM(filename, filenum, formatted,    &
     &                                 lupri)
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0   -   File does not exists.
!     -1   -   To many files opened.
!     -2   -   Number already assigned
!     -3   -   File already opened with a different number.
!     -99  -   Error while executing the open statement.
          CHARACTER(LEN=*), INTENT(IN)  :: filename
          INTEGER, INTENT(IN)           :: filenum
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: formatted
          INTEGER, OPTIONAL, INTENT(IN) :: lupri
          INTEGER, PARAMETER            :: maxfile=99
          LOGICAL                       :: ex
          INTEGER                       :: n

          INQUIRE(file=filename, exist=ex)
          IF (.NOT. ex) THEN
            open_file_num=0
            RETURN
          END IF
          INQUIRE(filenum, exist=ex)
          IF (.NOT. ex) THEN
            open_file_num=-2
            RETURN
          END IF
          OPEN(filenum, file=filename, FORM=formatted, err=463)   
          open_file_num=filenum

          RETURN
 463      open_file_num=-99

        END FUNCTION OPEN_FILE_NUM

        FUNCTION GET_EXTENSION(filename)
          CHARACTER(LEN=*), INTENT(IN) :: filename
          CHARACTER(LEN=LEN(filetypes(1))) :: get_extension
          INTEGER :: i, size_ft
          LOGICAL :: ex
  
          size_ft = SIZE(filetypes)
          DO i=1, size_ft
            INQUIRE(file=filename//'.'//TRIM(filetypes(i)), exist=ex)
            IF (ex) THEN
              get_extension=TRIM(filetypes(i))
              EXIT
            END IF
          END DO
          IF (i>size_ft)  get_extension=' '
        END FUNCTION
               
        subroutine filter_comments(filein,fileout)
          integer             :: filein,fileout
          character(len=192)  :: str
          character(len=1)    :: chr
          integer             :: i,ios

          rewind(filein)
          fileout=OPEN_SCRATCH()

          do 
            read(filein,'(192a)',iostat=ios) str
            if (ios.ne.0) exit
            do i=1,len(str)
              select case (str(i:i))
              case ('#','!')               
                exit
              end select
            end do      
            write(fileout,*) str(1:i-1)
          end do
          rewind(filein)
          rewind(fileout)
          

        end subroutine

      END MODULE FILE
