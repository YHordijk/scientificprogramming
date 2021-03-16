module include_dcbtra_h
!> module encapsulating the src/moltra/dcbtra.h include file
!> written by Miro Ilias, GSI.de, July 2016
implicit none
private
#include "dcbtra.h"
!> public method: set ASCII
public :: set_TRA_ASCII
contains
subroutine set_tra_ascii(tra_ascii_value)
!> sets the logical TRA_ASCII variable
logical, intent(in) :: tra_ascii_value
TRA_ASCII = tra_ascii_value
end subroutine set_tra_ascii

end module include_dcbtra_h
