module file_units

   implicit none

   character(*), parameter, public :: interface_file_name_ao      = 'interface_ao'
   character(*), parameter, public :: interface_file_name_mo      = 'interface_mo'
   character(*), parameter, public :: interface_file_name_mo_coef = 'interface_mo_coef'
   integer,      parameter, public :: interface_file_unit = 137
   logical                , public :: interface_file_exists = .false.
   logical                , public :: interface_file_open   = .false.

!export perturbed densities (for fde)
   character(*), parameter, public :: pertden_direct_lao_file_name   = 'pertden_direct_lao.FINAL'
   integer,      parameter, public :: pertden_direct_lao_file_unit   = 111
   logical                , public :: pertden_direct_lao_file_exists = .false.
   logical                , public :: pertden_direct_lao_file_open   = .false.

   character(*), parameter, public :: pertden_direct_lao_par_file_name   = 'pertden_direct_lao'
   integer,      parameter, public :: pertden_direct_lao_par_file_unit   = 333
   logical                , public :: pertden_direct_lao_par_file_exists = .false.
   logical                , public :: pertden_direct_lao_par_file_open   = .false.

   character(*), parameter, public :: pertden_reorth_lao_file_name   = 'pertden_reorth_lao.FINAL'
   integer,      parameter, public :: pertden_reorth_lao_file_unit   = 222
   logical                , public :: pertden_reorth_lao_file_exists = .false.
   logical                , public :: pertden_reorth_lao_file_open   = .false.

   character(*), parameter, public :: pertden_reorth_lao_par_file_name   = 'pertden_reorth_lao'
   integer,      parameter, public :: pertden_reorth_lao_par_file_unit   = 444
   logical                , public :: pertden_reorth_lao_par_file_exists = .false.
   logical                , public :: pertden_reorth_lao_par_file_open   = .false.

!  nr of columns for exported density
   integer                , public :: pertden_nr_columns      = 52
   character(*),parameter , public :: pertden_file_format     = '(52d20.12)'

   private

end module
