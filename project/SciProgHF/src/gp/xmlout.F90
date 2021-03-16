

! XML output module, Ulf Ekstrom 2011
! It's always safe to call the output functions even if no file has been opened,
! in that case there is simply no output.

module xmlout
  implicit none

  logical,public :: doxml
  data doxml/.false./
  integer, public :: luxml
  integer, private :: xml_id ! sometimes it's useful to have a unique id for some elements
  data xml_id/0/
  private

  public xml_open_file
  public xml_close_file
  public xml_begin
  public xml_end
  public xml_comment
  public xml_log
  public xml_newid
  public xml_matrix
  public xml_quantity
  public xml_embed_formatted

  integer, parameter :: max_depth = 32
  integer, parameter :: max_tag_len = 256
  integer :: depth
  character(max_tag_len) :: tag_stack(max_depth) = ' '

  interface xml_quantity
      module procedure xml_real
      module procedure xml_integer
  end interface

  contains
    function xml_newid()
      integer xml_newid
      xml_id = xml_id + 1
      xml_newid = xml_id
      return
    end function xml_newid

    subroutine push_tag(tag)
      character *(*) tag
      integer n
      if (depth.ge.max_depth) then
         print *,'ERROR in xmlout: max stack depth reached',depth
         if (doxml) write(luxml,*) 'ERROR in xmlout: max stack depth reached',depth
         return
      endif
      n = len_trim(tag)
      if (n.gt.max_tag_len) then
         print *,'ERROR in xmlout: tag too long',tag
         if (doxml) write(luxml,*) 'ERROR in xmlout: tag too long',tag
         return
      endif      
      depth = depth + 1
      tag_stack(depth) = tag(1:n)
    end subroutine push_tag

    subroutine pop_tag()
      if (depth.gt.0) then
         if (doxml) write (luxml,'(A,A,A)') '</',trim(tag_stack(depth)),'>'
         depth = depth - 1
      endif
    end subroutine pop_tag

    subroutine xml_begin(tag,attribs)
      character *(*) tag
      character *(*), optional :: attribs
      call push_tag(tag)
      if (doxml) then
         if (present(attribs)) then
            write (luxml,'(A,A,A,A,A)') '<',tag,' ',attribs,'>'
         else
            write (luxml,'(A,A,A)') '<',tag,'>'
         endif
      endif
    end subroutine xml_begin

    subroutine xml_end(tag)
      character *(*) tag
      if (depth.lt.1) then
         print *,'ERROR xmlout: stack empty when closing',tag
         if (doxml) write(luxml,*) 'ERROR xmlout: stack empty when closing',tag
         return
      endif
      if (trim(tag).ne.trim(tag_stack(depth))) then
         print *,'ERROR xmlout: attempting to close <',trim(tag_stack(depth)),'> with </',trim(tag),'>'
         if (doxml) write(luxml,*) &
              'ERROR xmlout: attempting to close <',trim(tag_stack(depth)),'> with </',trim(tag),'>'
         return
      endif
      call pop_tag()
    end subroutine xml_end

    integer function getUnit()
! find a free unit number (go fortran!)
      implicit none
      integer :: unit
      logical :: isOpen
      
      integer, parameter :: MIN_UNIT_NUMBER = 10
      integer, parameter :: MAX_UNIT_NUMBER = 399
      
      do unit = MAX_UNIT_NUMBER, MIN_UNIT_NUMBER,-1
         inquire(unit = unit, opened = isOpen)
         if (.not. isOpen) then
            getUnit = unit
            return
         end if
      end do
    end function getUnit

    subroutine xml_open_file(path)
      character *(*) path
      luxml = getUnit()
      open (unit=luxml,file=path,status='replace',form='formatted')
      write(luxml,'(A)') '<?xml version="1.0" encoding="UTF-8"?>'
      depth = 0
      doxml = .true.
    end subroutine xml_open_file

    subroutine xml_close_file()      
      integer i
      if (doxml) then
         if (depth.gt.0) then
	    call xml_comment('WARNING: open tag(s) found, closing..')
            do i=1,depth
	       call xml_comment('closing '//tag_stack(depth-i+1))
               call pop_tag()
            enddo
         endif
         close(luxml)
         doxml=.false.
      endif
    end subroutine xml_close_file

    subroutine print_datetime()
! prints a string to luxml with current date and time in xsd:dateTime format
! ulf: Did I get this right for time zones west of utc?
      integer, dimension(8) :: v
      integer :: utch,utcm
      character pm
      call date_and_time(VALUES=v)
      utch = v(4)/60
      utcm = v(4) - 60*utch
      if (utch.lt.0) then
         pm = '-'
         utch = -utch
      else
         pm = '+'
      endif
      write(luxml,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)',advance='no') &
           v(1),'-',v(2),'-', v(3),'T',v(5),':',v(6),':',v(7),pm,utch,':',utcm
    end subroutine print_datetime

    subroutine xml_log(note)
      character *(*) note
      if (doxml) then
         write(luxml,'(A)',advance='no') '<log time="'
         call print_datetime()
         write(luxml,'(A,A)',advance='no') '">',note
         write(luxml,'(A)') '</log>'
      endif
    end subroutine xml_log

    subroutine xml_comment(text)
      character *(*) text
      if (doxml) then
         write(luxml,'(A,A,A)') '<!--',text,'-->'
      endif
    end subroutine xml_comment

    subroutine xml_matrix(A,rows,lda,cols)
      double precision A
      integer rows,lda,cols,i
      dimension A(lda,cols)
      call push_tag('matrix')
      if (doxml) then
         write(luxml,'(A,I0.1,A,I0.1,A)') '<matrix type="real" rows="',rows,'" columns="', cols,'" packing="column">'
         do i=1,cols
            write(luxml,*) A(1:rows,i)
         enddo
      endif
      call xml_end('matrix')
    end subroutine xml_matrix

    subroutine xml_real(label,value,unit)
      character *(*) label,unit
      real(8),intent(in) :: value
      if (doxml) then
         write(luxml,'(A,A,A,A,A,F30.15,A)') '<quantity label="',trim(label),'" unit="',trim(unit),'">',value,'</quantity>'
      endif
    end subroutine xml_real

    subroutine xml_integer(label,value)
      character *(*) label
      integer,intent(in) :: value
      if (doxml) then
         write(luxml,'(A,A,A,I16,A)') '<quantity label="',trim(label),'">',value,'</quantity>'
      endif
    end subroutine xml_integer

! Read file from file_path and put it in the xml file.
! This is not very robust, one should use base64 encoding
    subroutine xml_embed_formatted(file_path,file_type)
      character file_path *(*), file_type *(*)
      character linebuf*(256)
      integer ios, lu
      write(LUXML,*) '<file name="'//trim(file_path)//'" type="'//trim(file_type)//'" encoding="UTF-8">'
      lu = getunit()
      OPEN(lu,FILE = file_path,IOSTAT=IOS)
      IF (IOS.NE.0) THEN
         CALL QUIT('Error in opening file for xml embed!')
      ENDIF
      do while (.true.)
         read(lu,'(A)',end=666) linebuf
         write(luxml,'(A)') trim(linebuf)
      enddo
666   continue
      write(LUXML,*) '</file>'
      close(lu)
    end subroutine xml_embed_formatted
end module xmlout
