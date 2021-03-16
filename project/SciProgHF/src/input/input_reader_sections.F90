   subroutine read_input_sections(word, kw_section)

      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      select case (kw_section)

!       case ('**METHO', '**WAVE ')
!          call read_input_method(word, kw_section)

        case ('*DFT   ')
           call read_input_dft(word, kw_section)

        case ('**GRID ')
           call read_input_grid(word, kw_section)

!       *VISUAL is for backward compatibility
!       *VISUAL used to be under **ANALYZE
        case ('**VISUA', '*VISUAL')
           call read_input_visual(word, kw_section)

        case ('*OPENRS')
           call read_input_openrsp(word, kw_section)

        case ('*X2C')
           call read_input_x2c(word, kw_section)

        case ('*aooSOC')
           call read_input_aoosoc(word, kw_section)

        case ('**EXACC')
           call read_input_exacc(word, kw_section)

        case ('**RELCC')
           call read_input_relcc(word, kw_section)

        case ('*CCREST')
           call read_input_ccrestart(word, kw_section)

        case ('*CCENER')
           call read_input_ccener(word, kw_section)

        case ('*CCFOPR')
           call read_input_ccfopr(word, kw_section)

        case ('*CCSORT')
           call read_input_ccsort(word, kw_section)

        case ('*CCFSPC')
           call read_input_ccfspc(word, kw_section)

        case ('*CCIH  ')
           call read_input_ccih(word, kw_section)

        case ('*EOMCC')
           call read_input_cceom(word, kw_section)

        case ('*EOMPROP')
           call read_input_eomprop(word, kw_section)

        case ('*CCPROJ')
           call read_input_ccproj(word, kw_section)

        case ('*CCDIAG')
           call read_input_ccdiag(word, kw_section)

#ifdef HAS_PCMSOLVER
        case ('*PCM   ')
           call read_input_pcm(word, kw_section)

        case ('*PCMSOL')
           call read_input_pcmsolver(word, kw_section)
#endif

        case ('**RELAD')
           call read_input_reladc(word, kw_section)

        case ('**POLPR')
           call read_input_polprp(word, kw_section)

        case ('**LANCZ')
           call read_input_lanczos(word, kw_section)

        case ('**DAVID')
           call read_input_davidson(word, kw_section)

        case ('*QCORR ')
           call read_input_qcorr(word, kw_section)

        case ('**DIAG ')  ! Miro's standalone program
           call read_input_diagonalization_tests(word, kw_section)

        case default
!          activate as soon as everything is merged to here
!          call quit('section '//kw_section//' not recognized')

      end select

   end subroutine

   subroutine read_input_grid(word, kw_section)

      use input_reader
      use num_grid_cfg

      implicit none

!     ---------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ---------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.RADINT')) then
         call kw_read(word, num_grid_cfg_radint)
      end if

      if (kw_matches(word, '.ANGINT')) then
         call kw_read(word, num_grid_cfg_angint)
      end if

      if (kw_matches(word, '.NOPRUN')) then
!        turn off pruning of angular grid
         num_grid_cfg_no_pruning = .true.
      end if

      if (kw_matches(word, '.ANGMIN')) then
         call kw_read(word, num_grid_cfg_angmin)
      end if

      if (kw_matches(word, '.ATSIZE')) then
!        estimate relative atomic sizes for use in the Becke
!        partitioning scheme from atomic contributions
         num_grid_cfg_estimate_radii = .true.
      end if

      if (kw_matches(word, '.IMPORT')) then
         num_grid_cfg_import_grid = .true.
         call kw_read(word, num_grid_cfg_gridfile)
      end if

!radovan: not documented on wiki
      if (kw_matches(word, '.ZIPGRI')) then
!        is true by default, keep for backward compatibimity
         num_grid_cfg_zipgrid = .true.
      end if

      if (kw_matches(word, '.NOZIP ')) then
         num_grid_cfg_zipgrid = .false.
      end if

      if (kw_matches(word, '.4CGRID')) then
!        include small component in the grid generation
!        also for 1- and 2-c dft calculations
         num_grid_cfg_force_4c_grid = .true.
      end if

      if (kw_matches(word, '.DEBUG ')) then
         num_grid_cfg_radint = 1.0d-3
         num_grid_cfg_angint = 10
      end if

      if (kw_matches(word, '.COARSE')) then
         num_grid_cfg_radint = 1.0d-11
         num_grid_cfg_angint = 35
      end if

      if (kw_matches(word, '.ULTRAF')) then
         num_grid_cfg_radint = 2.0d-15
         num_grid_cfg_angint = 64
      end if

      if (kw_matches(word, '.INTCHK')) then
         call kw_read(word, num_grid_cfg_integration_check_level)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_visual(word, kw_section)

      use dirac_cfg
      use visual_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i, j
      integer                          :: isym, istart, iend
      real(8)                          :: f, s, v(3), v_12(3), v_13(3)
      integer                          :: nr_lines
      character(100)                   :: line
      character(80), allocatable       :: word_array(:)
      integer                          :: nr_words
!     --------------------------------------------------------------------------

#include "dgroup.h"
#include "maxorb.h"
! orgcom.h : GAGORG
#include "orgcom.h"

      call reset_available_kw_list()

      visual_cfg_gauge_origin(1:3) = GAGORG(1:3) ! use magnetic gauge origin defined globally in Dirac
      dirac_cfg_visual = .true.

      if (kw_matches(word, '.LIST  ')) then
         visual_cfg_list = .true.
         call kw_read(word, visual_cfg_nr_points_in_list)
         allocate(visual_cfg_xyz_list(visual_cfg_nr_points_in_list, 3))
         do i = 1, visual_cfg_nr_points_in_list
            call kw_read(word, visual_cfg_xyz_list(i, 1), &
                               visual_cfg_xyz_list(i, 2), &
                               visual_cfg_xyz_list(i, 3))
         end do
      end if

      if (kw_matches(word, '.LINE  ')) then
         visual_cfg_line = .true.
         call kw_read(word, visual_cfg_line_from(1), &
                            visual_cfg_line_from(2), &
                            visual_cfg_line_from(3))
         call kw_read(word, visual_cfg_line_to(1), &
                            visual_cfg_line_to(2), &
                            visual_cfg_line_to(3))
         call kw_read(word, visual_cfg_line_nr_steps)
      end if

      if (kw_matches(word, '.2D    ')) then
         visual_cfg_2d = .true.
         read(get_file_unit(), *) visual_cfg_2d_p_origin(1:3)
         read(get_file_unit(), *) visual_cfg_2d_p_right(1:3)
         read(get_file_unit(), *) visual_cfg_2d_nr_right
         read(get_file_unit(), *) visual_cfg_2d_p_top(1:3)
         read(get_file_unit(), *) visual_cfg_2d_nr_top
!        check that angle between 1-2 and 1-3 is 90 degs
         v_12 = visual_cfg_2d_p_right - visual_cfg_2d_p_origin
         v_13 = visual_cfg_2d_p_top   - visual_cfg_2d_p_origin
         s = v_12(1)*v_13(1) + v_12(2)*v_13(2) + v_12(3)*v_13(3)
         if (s > tiny(0.0d0)) then
            call quit('.2D: lines connecting points 1-2 and 1-3 do not enclose 90 degrees')
         end if
      end if

      if (kw_matches(word, '.2D_INT')) then
         visual_cfg_2d_integration = .true.
         read(get_file_unit(), *) visual_cfg_2d_integration_p_origin(1:3)
         read(get_file_unit(), *) visual_cfg_2d_integration_p_right(1:3)
         read(get_file_unit(), *) visual_cfg_2d_integration_nr_right
         read(get_file_unit(), *) visual_cfg_2d_integration_p_top(1:3)
         read(get_file_unit(), *) visual_cfg_2d_integration_nr_top
         read(get_file_unit(), *) visual_cfg_2d_integration_order
!        check that angle between 1-2 and 1-3 is 90 degs
         v_12 = visual_cfg_2d_integration_p_right - visual_cfg_2d_integration_p_origin
         v_13 = visual_cfg_2d_integration_p_top   - visual_cfg_2d_integration_p_origin
         s = v_12(1)*v_13(1) + v_12(2)*v_13(2) + v_12(3)*v_13(3)
         if (s > tiny(0.0d0)) then
            call quit('.2D_INT: lines connecting points 1-2 and 1-3 do not enclose 90 degrees')
         end if
      end if

      if (kw_matches(word, '.3D    ')) then
         visual_cfg_3d = .true.
         call kw_read(word, visual_cfg_ncube(1), &
                            visual_cfg_ncube(2), &
                            visual_cfg_ncube(3))
      end if

      if (kw_matches(word, '.3DFAST')) then
         visual_cfg_3d_fast = .true.
         call kw_read(word, visual_cfg_ncube(1), &
                            visual_cfg_ncube(2), &
                            visual_cfg_ncube(3))
      end if

      if (kw_matches(word, '.3D_ADD')) then
         call kw_read(word, visual_cfg_add_3d)
      end if

      if (kw_matches(word, '.3D_INT')) then
         visual_cfg_3d_integration = .true.
      end if

      if (kw_matches(word, '.3D_IMP')) then
         visual_cfg_3d_gridimp = .true.
         read(get_file_unit(), *) visual_cfg_3d_gridfil
      end if

      if (kw_matches(word, '.SCALE ')) then
         call kw_read(word, visual_cfg_scale)
      end if

      if (kw_matches(word, '.CARPOW')) then
         call kw_read(word, visual_cfg_cartesian_power(1), &
                            visual_cfg_cartesian_power(2), &
                            visual_cfg_cartesian_power(3))
      end if

      if (kw_matches(word, '.DSCALE ')) then
         call kw_read(word, visual_cfg_scale)
         visual_cfg_scale=1.0D0/visual_cfg_scale
      end if

      if (kw_matches(word, '.DENSIT')) then
         call read_input_visual_helper(iq_density)
      end if

      if (kw_matches(word, '.ELF   ')) then
         call read_input_visual_helper(iq_elf)
      end if

      if (kw_matches(word, '.GAMMA5')) then
         call read_input_visual_helper(iq_gamma5)
      end if

      if (kw_matches(word, '.J     ')) then
         call read_input_visual_helper(iq_j)
      end if

      if (kw_matches(word, '.JDIA  ')) then
         call read_input_visual_helper(iq_jdia)
      end if

      if (kw_matches(word, '.JX    ')) then
         call read_input_visual_helper(iq_jx)
      end if

      if (kw_matches(word, '.JY    ')) then
         call read_input_visual_helper(iq_jy)
      end if

      if (kw_matches(word, '.JZ    ')) then
         call read_input_visual_helper(iq_jz)
      end if

      if (kw_matches(word, '.KIN   ')) then
         call read_input_visual_helper(iq_kin)
      end if

      if (kw_matches(word, '.KINLS ')) then
         call read_input_visual_helper(iq_kin_ls)
      end if

      if (kw_matches(word, '.KINSL ')) then
         call read_input_visual_helper(iq_kin_sl)
      end if

      if (kw_matches(word, '.KINTAU')) then
         call read_input_visual_helper(iq_kin_tau)
      end if

      if (kw_matches(word, '.KINLAP')) then
         call read_input_visual_helper(iq_kin_lap)
      end if

      if (kw_matches(word, '.KINNR ')) then
         call read_input_visual_helper(iq_kin_nr)
      end if

      if (kw_matches(word, '.DIVJ  ')) then
         call read_input_visual_helper(iq_divj)
      end if

      if (kw_matches(word, '.ROTJ  ')) then
         call read_input_visual_helper(iq_rotj)
      end if

      if (kw_matches(word, '.S     ')) then
         call read_input_visual_helper(iq_s)
      end if

      if (kw_matches(word, '.DIVS  ')) then
         call read_input_visual_helper(iq_divs)
      end if

      if (kw_matches(word, '.ROTS  ')) then
         call read_input_visual_helper(iq_rots)
      end if

      if (kw_matches(word, '.EDIPX ')) then
         call read_input_visual_helper(iq_edipx)
      end if
      if (kw_matches(word, '.EDIPY ')) then
         call read_input_visual_helper(iq_edipy)
      end if
      if (kw_matches(word, '.EDIPZ ')) then
         call read_input_visual_helper(iq_edipz)
      end if

      if (kw_matches(word, '.BDIPX ')) then
         call read_input_visual_helper(iq_bdipx)
      end if
      if (kw_matches(word, '.BDIPY ')) then
         call read_input_visual_helper(iq_bdipy)
      end if
      if (kw_matches(word, '.BDIPZ ')) then
         call read_input_visual_helper(iq_bdipz)
      end if

      if (kw_matches(word, '.BEDCOS')) then
         call read_input_visual_helper(iq_tcos)
      end if
      if (kw_matches(word, '.BEDSIN')) then
         call read_input_visual_helper(iq_tsin)
      end if
      if (kw_matches(word, '.BEDFIX')) then
         call kw_read(word, v(1),v(2),v(3),visual_cfg_freq)
         CALL ANGLES(v(1),v(2),v(3),visual_cfg_wave_vector,visual_cfg_pol_vector)
      end if

      if (kw_matches(word, '.NDIPX ')) then
         call read_input_visual_helper2(iq_ndipx)
      end if
      if (kw_matches(word, '.NDIPY ')) then
         call read_input_visual_helper2(iq_ndipy)
      end if
      if (kw_matches(word, '.NDIPZ ')) then
         call read_input_visual_helper2(iq_ndipz)
      end if

      if (kw_matches(word, '.BDIPXD')) then
         call read_input_visual_helper(iq_bdipxdia)
      end if
      if (kw_matches(word, '.BDIPYD')) then
         call read_input_visual_helper(iq_bdipydia)
      end if
      if (kw_matches(word, '.BDIPZD')) then
         call read_input_visual_helper(iq_bdipzdia)
      end if

      if (kw_matches(word, '.NDIPXD')) then
         call read_input_visual_helper2(iq_ndipxdia)
      end if
      if (kw_matches(word, '.NDIPYD')) then
         call read_input_visual_helper2(iq_ndipydia)
      end if
      if (kw_matches(word, '.NDIPZD')) then
         call read_input_visual_helper2(iq_ndipzdia)
      end if

      if (kw_matches(word, '.ESP   ')) then
         call read_input_visual_helper(iq_esp)
      end if
      if (kw_matches(word, '.ESPE  ')) then
         call read_input_visual_helper(iq_espe)
      end if
      if (kw_matches(word, '.ESPN  ')) then
         call read_input_visual_helper(iq_espn)
      end if
      if (kw_matches(word, '.ESPRHO')) then
         call read_input_visual_helper(iq_esprho)
      end if
      if (kw_matches(word, '.ESPERH')) then
         call read_input_visual_helper(iq_esperho)
      end if
      if (kw_matches(word, '.ESPNRH')) then
         call read_input_visual_helper(iq_espnrho)
      end if

      if (kw_matches(word, '.NICS  ')) then
         visual_cfg_nics = .true.
         call kw_read(word, visual_cfg_nics_origin(1), &
                            visual_cfg_nics_origin(2), &
                            visual_cfg_nics_origin(3))
      end if

      if (kw_matches(word, '.READJB')) then
         visual_cfg_3d_readjb = .true.
         read(get_file_unit(), *) visual_cfg_3d_jbfile
      end if

      if (kw_matches(word, '.RADIAL')) then
         visual_cfg_radial = .true.
         call kw_read(word, visual_cfg_radial_from(1), &
                            visual_cfg_radial_from(2), &
                            visual_cfg_radial_from(3))
         call kw_read(word, visual_cfg_radial_length)         
         call kw_read(word, visual_cfg_radial_nr_steps)
      end if

      if (kw_matches(word, '.SMALLA')) then
         visual_cfg_force_small_ao = .true.
      end if

      if (kw_matches(word, '.OCCUPA')) then
         visual_cfg_use_orbital_string = .true.
         allocate(visual_cfg_occupation(mxprim, 2))
         visual_cfg_occupation = 0.0d0
         ! read nr of lines
         call kw_read(word, nr_lines)
         print *, 'reading .OCCUPATION, nr lines:', nr_lines
         do i = 1, nr_lines
            ! read each line
            read(get_file_unit(), '(a)') line
            nr_words = word_count(line)
            allocate(word_array(nr_words))
            ! cut into words
            read(line, *) (word_array(j), j = 1, nr_words)
            ! read ifsym
            read(word_array(1), *) isym
            ! find '-'
            j = index(word_array(2), '-')
            if (j > 0) then
               ! read start orbital
               read(word_array(2)(1:j-1), *) istart
               ! read end orbital
               read(word_array(2)(j+1:len(word_array(2))), *) iend
            else
               ! start and end orbital are the same
               read(word_array(2), *) istart
               iend = istart
            end if
            ! sanity check
            if (istart < 1) then
               print *, 'ERROR: start orbital should be larger than 0'
               stop 1
            end if
            if (iend > mxprim) then
               print *, 'ERROR: end orbital higher than max nr of orbitals'
               stop 1
            end if
            ! read occupation
            read(word_array(3), *) f
            ! update occupation vector
            do j = istart, iend
               visual_cfg_occupation(j, isym) = f
            end do
            ! print confirmation to screen
            print *, 'occupy symmetry      = ', isym
            print *, '       start orbital = ', istart
            print *, '       end orbital   = ', iend
            print *, '       occupation    = ', f
            deallocate(word_array)
         end do
      end if

      if (kw_matches(word, '.LONDON')) then
         visual_cfg_london = .true.
         read(get_file_unit(), *) visual_cfg_london_component
      end if

      if (kw_matches(word, '.NONE  ')) then
         visual_cfg_london_none      = .true.
      end if

      if (kw_matches(word, '.NODIRE')) then
         visual_cfg_london_skip_direct = .true.
      end if

      if (kw_matches(word, '.NOREOR')) then
         visual_cfg_london_skip_ro = .true.
      end if

      if (kw_matches(word, '.NOKAPP')) then
         visual_cfg_london_skip_kappa = .true.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_visual_helper(j)

      use visual_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      integer, intent(in) :: j
!     --------------------------------------------------------------------------
      integer             :: i
      character(80)       :: line
!     --------------------------------------------------------------------------

      visual_cfg_nr_dmat = visual_cfg_nr_dmat + 1
      i = visual_cfg_nr_dmat
      visual_cfg_property(i) = j
      read(get_file_unit(), *) line
      backspace(get_file_unit())
      if (word_contains(line, 'DFCOEF')) then
         read(get_file_unit(), *) visual_cfg_dmat_file(i)
      else
         read(get_file_unit(), *) visual_cfg_dmat_file(i), &
                          visual_cfg_dmat_file_record(i)
      end if

   end subroutine

   subroutine read_input_visual_helper2(j)

      use visual_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      integer, intent(in) :: j
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      visual_cfg_nr_dmat = visual_cfg_nr_dmat + 1
      i = visual_cfg_nr_dmat
      visual_cfg_property(i) = j
      read(get_file_unit(), *) visual_cfg_ref_nucleus(i), &
                       visual_cfg_dmat_file(i),   &
                       visual_cfg_dmat_file_record(i)

   end subroutine

   subroutine read_input_openrsp(word, kw_section)

      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i
!     --------------------------------------------------------------------------

      call quit(kw_section//' not available in this version')

   end subroutine

   subroutine read_input_exacc(word, kw_section)

#ifdef MOD_EXACORR
      use exacc_cfg
#endif
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
      integer :: i,j
      real(8) :: val_real, val_imag
!     --------------------------------------------------------------------------

#ifndef MOD_EXACORR
      call quit(kw_section//' not available in this version')
#else
      call reset_available_kw_list()

      if (kw_matches(word, '.OCCUPIED')) then
         call kw_read(word, string_occupied)
      end if

      if (kw_matches(word, '.VIRTUAL')) then
         call kw_read(word, string_virtual)
      end if

      if (kw_matches(word, '.EXATENSOR')) then
         exa_input%talsh = .false.
      end if

      if (kw_matches(word, '.CCDOUBLES')) then
         exa_input%ccd = .true.
      end if

      if (kw_matches(word, '.CC2         ')) then 		!L
         exa_input%cc2 = .true.
      end if

      if (kw_matches(word, '.TCONVERG')) then
         call kw_read(word, exa_input%t_econv)
      end if

      if (kw_matches(word, '.LSHIFT')) then
         call kw_read(word, exa_input%level_shift)
      end if

      if (kw_matches(word, '.CHOLESKY')) then
         call kw_read(word, exa_input%t_cholesky)
      end if

      if (kw_matches(word, '.LAMBDA')) then
         exa_input%lambda = .true.
      end if

      if (kw_matches(word, '.NOTRIPLES')) then
         exa_input%do_triples = .false.
      end if

      if (kw_matches(word, '.TRIPL_BLOCK')) then
         call kw_read(word, exa_input%tripl_block)
      end if

      if (kw_matches(word, '.EXA_BLOCKSIZE')) then
         call kw_read(word, exa_input%exa_blocksize)
      end if

      if (kw_matches(word, '.MOINT_SCHEME')) then
         call kw_read(word, exa_input%moint_scheme)
      end if

      if (kw_matches(word, '.NCYCLES')) then
         call kw_read(word, exa_input%ncycles)
      end if

      if (kw_matches(word, '.OCC_BETA')) then
         exa_input%beta_occ = .true.
         call kw_read(word, string_occ_beta)
      end if

      if (kw_matches(word, '.VIR_BETA')) then
         exa_input%beta_vir = .true.
         call kw_read(word, string_vir_beta)
      end if

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, exa_input%print_level)
      end if

      if (kw_matches(word, '.TALSH_BUFF')) then
         call kw_read(word, exa_input%talsh_buff)
      end if

      if (kw_matches(word, '.FF_PROP')) then
         call kw_read(word, exa_input%nff)
         allocate(exa_input%ff_names(exa_input%nff(1)))
         allocate(exa_input%ff(exa_input%nff(1),exa_input%nff(2)))
         do i=1,exa_input%nff(1)
           call kw_read(word, exa_input%ff_names(i))
           do j=1,exa_input%nff(2)
             call kw_read(word, val_real,val_imag)
             exa_input%ff(i,j)=dcmplx(val_real,val_imag)
           enddo
         enddo
      end if

      call check_whether_kw_found(word, kw_section)

#endif

   end subroutine

   subroutine read_input_relcc(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
      logical :: explicitly_set_energy = .false.
!     --------------------------------------------------------------------------

#include "dgroup.h"

      call reset_available_kw_list()

      if (kw_matches(word, '.NELEC ')) then
         nelec_input = .true.
         call kw_read(word, relcc_nelec, 2*nfsym)
      end if

      if (kw_matches(word, '.NELEC_OPEN ')) then
         nelec_open_input = .true.
         call kw_read(word, relcc_nelec_open, 2*nfsym)
      end if

      if (kw_matches(word, '.NEL_F1')) then
         call kw_read(word, relcc_nelec_f1)
      end if

      if (kw_matches(word, '.NEL_F2')) then
         call kw_read(word, relcc_nelec_f2)
      end if

      if (kw_matches(word, '.NFROZ ')) then
         call kw_read(word, relcc_nfroz, 2*nfsym)
      end if

      if (kw_matches(word, '.DEBUG ')) then
         relcc_debug = .true.
      end if

      if (kw_matches(word, '.CARITH')) then
         relcc_carith = .true.
      end if

      if (kw_matches(word, '.TIMING')) then
         relcc_timing = .true.
      end if

      if (kw_matches(word, '.RESTART')) then
         relcc_do_restart = .true.
      end if

      if (kw_matches(word, '.NORESTART')) then
         relcc_do_restart = .false.
      end if

      if (kw_matches(word, '.NOSORT')) then
         relcc_do_sort = .false.
      end if

      if (kw_matches(word, '.COUNTMEM')) then
         relcc_do_count_memory = .true.
      end if

      if (kw_matches(word, '.ENERGY')) then
         relcc_do_energy = .true.
         explicitly_set_energy = .true.
      end if

      if (kw_matches(word, '.GRADIE')) then
         relcc_do_gradient = .true.
!         relcc_do_mp2gradient = .true.
      end if

      if (kw_matches(word, '.HESSIA')) then
         relcc_do_hessian  = .true.
      end if

      if (kw_matches(word, '.EOMCC ')) then
         relcc_do_eomcc  = .true.
      end if

      if (kw_matches(word, '.EOMPROP')) then
         relcc_do_eomprop  = .true.
      end if

      if (kw_matches(word, '.FOCKSP')) then
         relcc_do_fspc     = .true.
         if (.not. explicitly_set_energy) relcc_do_energy = .false.
      end if

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, relcc_print)
      end if

      if (kw_matches(word, '.INTERFACE')) then
         call kw_read(word, relcc_integral_interface)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccener(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.NOMP2 ')) then
         relcc_do_mp2 = .false.
      end if

      if (kw_matches(word, '.NOSING')) then
         relcc_no_singles = .true.
      end if

      if (kw_matches(word, '.NODOUB')) then
         relcc_no_doubles = .true.
      end if

      if (kw_matches(word, '.NOSD  ')) then
         relcc_do_ccsd = .false.
      end if

      if (kw_matches(word, '.NOSDT ')) then
         relcc_do_ccsd_t = .false.
      end if

      if (kw_matches(word, '.MAXIT ')) then
         read(get_file_unit(),'(i5)') relcc_ccener_max_iterations
      end if

      if (kw_matches(word, '.MAXDIM')) then
         read(get_file_unit(),'(i5)') relcc_ccener_max_dimension_diis
      end if

      if (kw_matches(word, '.NTOL  ')) then
         read(get_file_unit(),'(i5)') relcc_ccener_desired_convergence
      end if

! miro (luuk) - expert's keyword
      if (kw_matches(word, '.DHOLU ')) then
         call kw_read(word, relcc_ccener_dholu_limit)
         relcc_ccener_dholu_limit_set = .true.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccfopr(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.MP2G   ')) then
         relcc_do_mp2gradient = .true.
      end if

      if (kw_matches(word, '.MP2G_O')) then
         relcc_do_oldmp2gradient = .true.
      end if

      if (kw_matches(word, '.CCSDG ')) then
         relcc_do_ccsdgradient = .true.
      end if

      if (kw_matches(word, '.CCSDTG')) then
         relcc_do_ccsdtgradient = .true.
      end if

      if (kw_matches(word, '.NATORB')) then
         relcc_do_naturalorbitals = .true.
      end if

      if (kw_matches(word, '.RELAXED')) then
         relcc_do_relaxed = .true.
      end if

      if (kw_matches(word, '.NEOPER')) then
         read(get_file_unit(),'(i5)') relcc_ne_oper
      end if

      if (kw_matches(word, '.MAXIT ')) then
         read(get_file_unit(),'(i5)') relcc_fopr_max_iterations
      end if

      if (kw_matches(word, '.MAXDIM')) then
         read(get_file_unit(),'(i5)') relcc_fopr_max_dimension_diis
      end if

      if (kw_matches(word, '.NTOL  ')) then
         read(get_file_unit(),'(i5)') relcc_fopr_desired_convergence
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine


   subroutine read_input_cceom(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
       integer                          :: nr_roots(2) ! number of roots in each symmetry 
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.EE    ')) then
         relcc_do_eomee=.true.
         call kw_read(word,  nr_roots)
         relcc_eom_nroots(nr_roots(1)) = nr_roots(2)           
      end if

      if (kw_matches(word, '.IP    ')) then
         relcc_do_eomip=.true.
         call kw_read(word,  nr_roots)
         relcc_eom_nroots(nr_roots(1)) = nr_roots(2)           
      end if

      if (kw_matches(word, '.EA    ')) then
         relcc_do_eomea=.true.
         call kw_read(word,  nr_roots)
         relcc_eom_nroots(nr_roots(1)) = nr_roots(2)           
      end if

      call check_whether_kw_found(word, kw_section)

      if (    (relcc_do_eomip .and. relcc_do_eomea) &
          .or.(relcc_do_eomip .and. relcc_do_eomee) &
          .or.(relcc_do_eomea .and. relcc_do_eomee)) &
          call quit('Only one EOM kind of calculation at a time (IP/EE/EA) allowed !')


   end subroutine


   subroutine read_input_eomprop(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
       integer   :: nr_roots_prop(2),nr_roots_left(2),nr_roots_right(2) ! number of roots in each symmetry 
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.EXCPRP')) then
         relcc_do_excprp=.true.
         call kw_read(word,  nr_roots_prop)
         relcc_eom_nroots_prop(nr_roots_prop(1)) = nr_roots_prop(2)           
      end if

      if (kw_matches(word, '.LEFT_STATE')) then
         relcc_do_excprp=.true.
         call kw_read(word,  nr_roots_left)
         relcc_eom_nroots_left(nr_roots_left(1)) = nr_roots_left(2)           
      end if

      if (kw_matches(word, '.RIGHT_STATE')) then
         relcc_do_excprp=.true.
         call kw_read(word,  nr_roots_right)
         relcc_eom_nroots_right(nr_roots_right(1)) = nr_roots_right(2)           
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine


   subroutine read_input_ccproj(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.REW_STRICT')) then
         relcc_projectors_rew_strict = .true.
      end if

      if (kw_matches(word, '.REW_NOSTRICT')) then
         relcc_projectors_rew_strict = .false.
      end if

      if (kw_matches(word, '.REW_OCC')) then
         relcc_projectors_do_restricted_excitation_window = .true.
         call kw_read(word,  relcc_projectors_rew_occ_min_energy)
         call kw_read(word,  relcc_projectors_rew_occ_max_energy)
      end if

      if (kw_matches(word, '.REW_VIRT')) then
         relcc_projectors_do_restricted_excitation_window = .true.
         call kw_read(word,  relcc_projectors_rew_virt_min_energy)
         call kw_read(word,  relcc_projectors_rew_virt_max_energy)
      end if

      if (kw_matches(word, '.CVS_CORE')) then
         relcc_projectors_do_core_valence_separation = .true.
         call kw_read(word,  relcc_projectors_rew_occ_max_energy)
         relcc_projectors_rew_occ_min_energy = -huge(relcc_projectors_rew_occ_min_energy)
         relcc_projectors_rew_virt_min_energy = -huge(relcc_projectors_rew_virt_min_energy)
         relcc_projectors_rew_virt_max_energy = huge(relcc_projectors_rew_virt_min_energy)
      end if

      if (kw_matches(word, '.FCORE ')) then
         relcc_projectors_frozen_core = .true.
         call kw_read(word,  relcc_projectors_frozen_core_max_energy)
      end if

      if (kw_matches(word, '.NODOCC')) then
         relcc_projectors_rew_remove_double_occupied = .true.
      end if

      call check_whether_kw_found(word, kw_section)

      if (relcc_projectors_do_restricted_excitation_window.and.relcc_projectors_do_core_valence_separation) & 
          call quit('Only one projection scheme (CVS/REW) allowed at a time !')


   end subroutine


   subroutine read_input_ccdiag(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.DEBUG ')) then
         relcc_mfd_verbose = .true.
      end if

      if (kw_matches(word, '.MAXSIZE')) then
         call kw_read(word,  relcc_mfd_max_subspace_size)
      end if

      if (kw_matches(word, '.REFRESH')) then
         call kw_read(word, relcc_mfd_refresh_rate)
      end if

      if (kw_matches(word, '.MAXITER')) then
         call kw_read(word, relcc_mfd_max_iterations)
      end if

      if (kw_matches(word, '.CONVERG')) then
         call kw_read(word, relcc_mfd_convergence_threshold)
      end if

      if (kw_matches(word, '.TRV_CCS')) then
         relcc_mfd_trial_ccs      = .true.
         relcc_mfd_trial_diagonal = .not.relcc_mfd_trial_ccs
         relcc_mfd_trial_full_matrix = .not.relcc_mfd_trial_ccs
         relcc_mfd_trial_restart = .not.relcc_mfd_trial_ccs
      end if

      if (kw_matches(word, '.TRV_I  ')) then
         relcc_mfd_trial_diagonal = .true.
         relcc_mfd_trial_ccs      = .not.relcc_mfd_trial_diagonal 
         relcc_mfd_trial_full_matrix = .not.relcc_mfd_trial_diagonal
         relcc_mfd_trial_restart  = .not.relcc_mfd_trial_diagonal
      end if

      if (kw_matches(word, '.TRV_FULLMATRIX')) then
         relcc_mfd_trial_full_matrix = .true.
         relcc_mfd_trial_diagonal = .not.relcc_mfd_trial_full_matrix
         relcc_mfd_trial_ccs      = .not.relcc_mfd_trial_full_matrix
         relcc_mfd_trial_restart  = .not.relcc_mfd_trial_full_matrix
      end if

      if (kw_matches(word, '.TRV_RESTART')) then
         relcc_mfd_trial_restart = .true.
         relcc_mfd_trial_diagonal = .not.relcc_mfd_trial_restart
         relcc_mfd_trial_ccs      = .not.relcc_mfd_trial_restart
         relcc_mfd_trial_restart  = .not.relcc_mfd_trial_restart
      end if

      if (kw_matches(word, '.TRV_R2L ')) then
         relcc_mfd_trial_lhs_use_rhs = .true.
      end if

      if (kw_matches(word, '.TRV_NOR2L')) then
         relcc_mfd_trial_lhs_use_rhs = .false.
      end if

      if (kw_matches(word, '.OVERLAP')) then
         relcc_mfd_overlap_sorting = .true.
      end if

      if (kw_matches(word, '.NOOVERLAP')) then
         relcc_mfd_overlap_sorting = .false.
      end if

      if (kw_matches(word, '.EV_Re_SHIFT')) then
         call kw_read(word, relcc_mfd_eigenvalues_energy_shift(1))
      end if

      if (kw_matches(word, '.EV_Im_SHIFT')) then
         call kw_read(word, relcc_mfd_eigenvalues_energy_shift(2))
      end if

      call check_whether_kw_found(word, kw_section)


   end subroutine



   subroutine read_input_ccsort(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
      integer :: i
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

!     Use HF orbital energies
      if (kw_matches(word, '.USEOE ')) then
         relcc_use_orbital_energies = .true.
      end if

!     do not re-compute orbital energies
      if (kw_matches(word, '.NORECMP')) then
         relcc_no_recompute   = .true.
      end if

! miro&luuk: expert option - adding finite field perturbation operators
      if (kw_matches(word, '.NFFOPER')) then
#ifdef MOD_UNRELEASED
         call kw_read(word, relcc_nffoper)
         relcc_ADD_FINITE_FIELD = .TRUE.
         ! after reading the number of perturbative operators
         ! read each one after one  - property name and its complex field strength
         do i=1,relcc_nffoper
           call kw_read(word, relcc_FF_PROP_NAMES(i))
           call kw_read(word, relcc_FF_PROP_STRENGTHS(1,i),relcc_FF_PROP_STRENGTHS(2,i))
         enddo
#else
         call quit('.NFFOPER keyword not in the Dirac release version !')
#endif
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccrestart(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i
      integer                          :: j
      integer                          :: nsect
      character(2)                     :: ifss(6)
      character(2)                     :: read_sect(6)
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

!     forces results from unconverged iterative procedures to be taken as converged
      if (kw_matches(word, '.UNCONV')) then
         relcc_restart_unconverged = .true.
      end if

!     in the case of fock-space calculations, indicate which sectors are to be
!     skipped (e.g. because converged, or otherwise)
      if (kw_matches(word, '.SKIPSECT')) then

         do i = 1, 6
            read_sect(i) = '  '
         enddo

         read(get_file_unit(), *) nsect
         read(get_file_unit(), *) (read_sect(i), i = 1, nsect)

         ifss(1) = '00'
         ifss(2) = '01'
         ifss(3) = '10'
         ifss(4) = '11'
         ifss(5) = '02'
         ifss(6) = '20'

         do i = 1, nsect
            do j = 1, 6
               if (read_sect(i).eq.IFSS(j)) then
                  relcc_restart_skipsect(j) = .true.
               endif
            enddo
         enddo
      end if

!     forces restart even if setup is potentially different (number of electrons, active/inactive etc)
      if (kw_matches(word, '.FORCER')) then
         relcc_restart_ignore_check = .true.
      end if

!     forces sorting to be done even if it has been completed before.
      if (kw_matches(word, '.REDOCCSD')) then
         relcc_restart_redo_ccsd = .true.
      end if

!     forces sorting to be done even if it has been completed before.
      if (kw_matches(word, '.REDOSORT')) then
         relcc_restart_redo_sorting = .true.
      end if

!     in the case of fock-space calculations, indicate which sectors are to be
!     recalculated even if converged.
      if (kw_matches(word, '.REDOSECT')) then

         do i = 1, 6
            read_sect(i) = '  '
         enddo

         read(get_file_unit(), *) nsect
         read(get_file_unit(), *) (read_sect(i), i = 1, nsect)

         ifss(1) = '00'
         ifss(2) = '01'
         ifss(3) = '10'
         ifss(4) = '11'
         ifss(5) = '02'
         ifss(6) = '20'

         do i = 1, nsect
            do j = 1, 6
               if (read_sect(i).eq.IFSS(j)) then
                  relcc_restart_redosect(j) = .true.
               endif
            enddo
         enddo
      end if



      call check_whether_kw_found(word, kw_section)

   end subroutine


   subroutine read_input_ccfspc(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i
!     --------------------------------------------------------------------------

#include "dgroup.h"

      call reset_available_kw_list()

      if (kw_matches(word, '.DOEA  ')) then
         relcc_fs_do_ea = .true.
      end if

      if (kw_matches(word, '.DOIE  ')) then
         relcc_fs_do_ie = .true.
      end if

      if (kw_matches(word, '.DOEA2 ')) then
         relcc_fs_do_ea2 = .true.
      end if

      if (kw_matches(word, '.DOIE2 ')) then
         relcc_fs_do_ie2 = .true.
      end if

      if (kw_matches(word, '.DOEXC ')) then
         relcc_fs_do_exc = .true.
      end if

      if (kw_matches(word, '.FSSECT ')) then
         read(get_file_unit(),*) (relcc_fs_fssect(i), i=1,6)
      end if

      if (kw_matches(word, '.DOIH  ')) then
         relcc_fs_do_ih = .true.
      end if

      if (kw_matches(word, '.NACTH ')) then
         call kw_read(word, relcc_fs_nacth, 2*NFSYM)
      end if

      if (kw_matches(word, '.NACTP ')) then
         call kw_read(word, relcc_fs_nactp, 2*NFSYM)
      end if

      if (kw_matches(word, '.MAXIT ')) then
         read(get_file_unit(),'(i5)') relcc_fs_max_iterations
      end if

      if (kw_matches(word, '.MXIT02')) then
         read(get_file_unit(),'(i5)') relcc_fs_max02_iterations
      end if

      if (kw_matches(word, '.MXIT01')) then
         read(get_file_unit(),'(i5)') relcc_fs_max01_iterations
      end if

      if (kw_matches(word, '.MXIT20')) then
         read(get_file_unit(),'(i5)') relcc_fs_max20_iterations
      end if

      if (kw_matches(word, '.MXIT10')) then
         read(get_file_unit(),'(i5)') relcc_fs_max10_iterations
      end if

      if (kw_matches(word, '.MXIT11')) then
         read(get_file_unit(),'(i5)') relcc_fs_max11_iterations
      end if

      if (kw_matches(word, '.MXIT00')) then
         read(get_file_unit(),'(i5)') relcc_fs_max00_iterations
      end if

      if (kw_matches(word, '.MAXDIM')) then
         read(get_file_unit(),'(i5)') relcc_fs_max_dimension_diis
      end if

      if (kw_matches(word, '.NTOL  ')) then
         read(get_file_unit(),'(i5)') relcc_fs_desired_convergence
      end if

      if (kw_matches(word, '.TSHOLD')) then
         read(get_file_unit(),'(i5)') relcc_fs_tshold
      end if

      if (kw_matches(word, '.GESTAT')) then
         read(get_file_unit(),'(i5)') relcc_fs_select_state_for_numgrad_energy
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccih(word, kw_section)

      use relcc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

#include "dgroup.h"

      call reset_available_kw_list()

      if (kw_matches(word, '.IHSCHE')) then
         read(get_file_unit(),'(i5)') relcc_ih_scheme
      end if

      if (kw_matches(word, '.SH_H11')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_h11
      end if

      if (kw_matches(word, '.SH_H12')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_h12
      end if

      if (kw_matches(word, '.SH_P11')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_p11
      end if

      if (kw_matches(word, '.SH_P12')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_p12
      end if

      if (kw_matches(word, '.SH_P2 ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_p2
      end if

      if (kw_matches(word, '.SH_H2 ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_h2
      end if

      if (kw_matches(word, '.SH_HP ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_shift_hp
      end if

      if (kw_matches(word, '.AIH   ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_aih
      end if

      if (kw_matches(word, '.EHMIN ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_eh_min
      end if

      if (kw_matches(word, '.EHMAX ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_eh_max
      end if

      if (kw_matches(word, '.EPMIN ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_ep_min
      end if

      if (kw_matches(word, '.EPMAX ')) then
         read(get_file_unit(),'(f15.8)') relcc_ih_ep_max
      end if

      if (kw_matches(word, '.NIH   ')) then
         read(get_file_unit(),'(i5)') relcc_ih_nih
      end if

      if (kw_matches(word, '.NACTHI')) then
         call kw_read(word, relcc_ih_nacthi, 2*NFSYM)
      end if

      if (kw_matches(word, '.NACTPI')) then
         call kw_read(word, relcc_ih_nactpi, 2*NFSYM)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine
!
!  reladc input subroutine
!
   subroutine read_input_reladc(word, kw_section)

      use adc_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i, j, k !counters
      character(100)                   :: line
      integer                          :: temp, already
      integer                          :: user_req_symmetries
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.ADCLEV')) then
         read(get_file_unit(),'(i5)') reladc_adclevel
      end if

      if (kw_matches(word, '.DOSIPS')) then
         reladc_dosips = .true.
      end if

      if (kw_matches(word, '.DODIPS')) then
         reladc_dodips = .true.
      end if

      if (kw_matches(word, '.FANO  ')) then
         reladc_dofano = .true.
      end if

      if (kw_matches(word, '.SIPREP')) then
         call kw_read(word, user_req_symmetries)
         !call kw_read(word, reladc_sipreps, user_req_symmetries)
         call kw_read(word, reladc_sipreps, reladc_no_sipreps, user_req_symmetries)
      end if

      if (kw_matches(word, '.DIPREP')) then
         call kw_read(word, user_req_symmetries)
         call kw_read(word, reladc_dipreps, user_req_symmetries)
      end if

      if (kw_matches(word, '.READQK')) then
         reladc_readqkl = .true.
      end if

      if (kw_matches(word, '.NOCONS')) then
         reladc_doconst = .false.
      end if

      if (kw_matches(word, '.DOADCP')) then
         reladc_doadcpop = .true.
      end if

      if (kw_matches(word, '.DOFULL')) then
         reladc_dofull = .true.
      end if

      if (kw_matches(word, '.DOLANC')) then
         reladc_dolanc = .true.
      end if

      if (kw_matches(word, '.VCONV ')) then
         read(get_file_unit(),'(f15.8)') reladc_vconv
      end if

      if (kw_matches(word, '.ADCTHR')) then
         read(get_file_unit(),'(f18.6)') reladc_adcthr
      end if

      if (kw_matches(word, '.FANOIN')) then
         call kw_read(word, reladc_fano_inrep)
         reladc_no_sipreps = reladc_no_sipreps + 1
         reladc_sipreps(reladc_no_sipreps) = reladc_fano_inrep
         reladc_sipeigv(reladc_no_sipreps) = 0
         call kw_read(word, reladc_fano_inrelsp) ! use relative spinor numbers
      end if

      if (kw_matches(word, '.FANOCHNL')) then
         ! read number of channels
         call kw_read(word, reladc_fano_nrgroups)
         allocate(reladc_fano_groups(reladc_fano_nrgroups))
         allocate(reladc_fano_labels(reladc_fano_nrgroups))
         reladc_fano_groups = 0
         call kw_read(word, reladc_fano_groups, reladc_fano_nrgroups)
         reladc_fano_nrchannels = SUM(reladc_fano_groups)
         allocate(reladc_fano_channels(reladc_fano_nrchannels,2))
         reladc_fano_channels = 0
         already = 0
         print *, ''
         print *, 'Absolute spinor numbers of Fano final states: '
         do k = 1, reladc_fano_nrgroups
           call kw_read(word,reladc_fano_labels(k))
           do i = 1, reladc_fano_groups(k)
              ! read every line
              read(get_file_unit(), '(a)') line
              ! read the absolute spinor numbers of the final states
              read(line,*) (reladc_fano_channels(i+already,j), j = 1, 2)
              if (reladc_fano_channels(i+already,1) < reladc_fano_channels(i+already,2)) then
                temp = reladc_fano_channels(i+already,1)
                reladc_fano_channels(i+already,1) = reladc_fano_channels(i+already,2)
                reladc_fano_channels(i+already,2) = temp
              end if 
              print *, (reladc_fano_channels(i+already,j), j=1,2)
              if (reladc_fano_channels(i+already,1) == reladc_fano_channels(i+already,2)) then
                call quit ( 'Equal spinors in one channel -> Unphysical choice')
              end if
           end do
           already = already + reladc_fano_groups(k)
         end do
      end if

      if (kw_matches(word, '.FANOONLY')) then
         reladc_dosips     = .true.
         reladc_sipreps    = 0
         reladc_sipreps(1) = reladc_fano_inrep
         reladc_dodiag     = .false.
         reladc_fanoonly   = .true.
      end if

      if (kw_matches(word, '.FANOFOVL')) then
         call kw_read(word, reladc_fano_fovl)
      end if

      if (kw_matches(word, '.FANOFCHK')) then
         reladc_fano_checkfin = .true.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine
!
!  polprp input subroutine
!
   subroutine read_input_polprp(word, kw_section)

      use polprp_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: user_req_symmetries
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.STATES')) then
         call kw_read(word, user_req_symmetries)
         call kw_read(word, polprp_statesym, user_req_symmetries)
      end if

      if (kw_matches(word, '.DOEXTE')) then
         polprp_doextended = .true.
      end if

      if (kw_matches(word, '.WRITET')) then
         read(get_file_unit(),'(f15.8)') polprp_writethr
      end if

      if (kw_matches(word, '.NODIAG')) then
         polprp_dodiag = .false.
      end if

      if (kw_matches(word, '.DOTRMO')) then
         polprp_dotrmo = .true.
      end if

      if (kw_matches(word, '.SKIPCC')) then
         polprp_skipccseti = .true.
      end if

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, polprp_printlev)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine
!
!  lanczos input subroutine
!
   subroutine read_input_lanczos(word, kw_section)

      use adc_cfg   ! lanczos only used in combination with reladc. same module
      use adc_mat
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: user_req_symmetries
      character(100)                   :: line
      integer                          :: i, j
!     --------------------------------------------------------------------------


      call reset_available_kw_list()

      if (kw_matches(word, '.SIPITER')) then
         read(get_file_unit(),'(i8)') reladc_sipiter
      end if

      if (kw_matches(word, '.DIPITER')) then
         read(get_file_unit(),'(i8)') reladc_dipiter
      end if

      if (kw_matches(word, '.SIPEIGV')) then
         call kw_read(word, user_req_symmetries)
         print *, 'lower and upper energy boundaries for ADC single ion. eigenvectors'
         do i = 1, user_req_symmetries
            ! read every line
            read(get_file_unit(), '(a)') line
            read(line,*) (reladc_md_sip_eeigv(i,j), j = 1, 2)
            print *, (reladc_md_sip_eeigv(i,j), j=1,2)
         end do
      end if

      if (kw_matches(word, '.DIPEIGV')) then
         call kw_read(word, user_req_symmetries)
         print *, 'lower and upper energy boundaries for ADC double ion. eigenvectors'
         do i = 1, user_req_symmetries
            ! read every line
            read(get_file_unit(), '(a)') line
            read(line,*) (reladc_md_dip_eeigv(i,j), j = 1, 2)
            print *, (reladc_md_dip_eeigv(i,j), j=1,2)
         end do
         !read(get_file_unit(),'(i8)') user_req_symmetries
         !call kw_read(word, reladc_dipeigv, user_req_symmetries)
      end if

      if (kw_matches(word, '.SIPPRNT')) then
         read(get_file_unit(),'(f18.6)') reladc_sipprnt
      end if

      if (kw_matches(word, '.DIPPRNT')) then
         read(get_file_unit(),'(f18.6)') reladc_dipprnt
      end if

      if (kw_matches(word, '.DOINCORE')) then
         reladc_doincore = .true.
      end if

      if (kw_matches(word, '.LANCMEM')) then
         read(get_file_unit(),'(i12)') reladc_lancmem
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine
!
!  davidson input subroutine (now separate from Lanczos for consistency)
!
   subroutine read_input_davidson(word, kw_section)

      use adc_cfg   ! lanczos only used in combination with reladc. same module
      use adc_mat
      use polprp_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.DVROOTS')) then
         read(get_file_unit(),'(i12)') reladc_davroots
         polprp_davroots = reladc_davroots
      end if

      if (kw_matches(word, '.DVMAXSP')) then
         read(get_file_unit(),'(i12)') reladc_davmaxsp
         polprp_davmaxsp = reladc_davmaxsp
      end if

      if (kw_matches(word, '.DVMAXIT')) then
         read(get_file_unit(),'(i12)') reladc_davmaxit
         polprp_davmaxit = reladc_davmaxit
      end if

      if (kw_matches(word, '.DVCONV')) then
         read(get_file_unit(),'(f18.6)') reladc_davconv
         polprp_davconv = reladc_davconv
      end if

      if (kw_matches(word, '.DVREORT')) then
         polprp_davreort = .true.
      end if

      if (kw_matches(word, '.DVINCORE')) then
         reladc_davooc = .false.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine
!
!  x2c 
!
   subroutine read_input_x2c(word, kw_section)

      use x2cmod_cfg
      use fragment_x2c_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.fragx2')) then
         x2cmod_fragment_x2c                            = .true.
         fragment_x2c_info%fragment_approach_enabled    = .true.
      end if

      if (kw_matches(word, '.fragmo')) then
         x2cmod_fragment_x2c                            = .true.
         fragment_x2c_info%fragment_approach_enabled    = .true.
         fragment_x2c_info%fragment_approach_ismolecule = .true.
      end if

      if (kw_matches(word, '.DEBUG ')) then
         x2cmod_debug = .true.
      end if

      if (kw_matches(word, '.hdirty')) then
         x2cmod_h1dirty = .true.
      end if

!     free-particle matrix as defining h1 hamiltonian
      if (kw_matches(word, '.freeh1')) then
         x2c_free_part_mtx_defh1 = .true.
      end if

      if (kw_matches(word, '.fockao')) then
         x2c_fock_saao_basis = .true.
      end if

!     completely skip any PCT of property operators in restart X2C energy calcs
      if (kw_matches(word, '.skippc')) then
         x2cmod_skip_op_pct = .true.
      end if

!     restart mmf calculation
      if (kw_matches(word, '.mmf-re')) then
         x2cmod_mmf_restart = .true.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_aoosoc(word, kw_section)

      use x2cmod_cfg
      use input_reader
#ifdef MOD_AOOSOC
      use atomic_oo_order_so_correction_cfg
#endif

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

#ifdef MOD_AOOSOC
      call reset_available_kw_list()


      if (kw_matches(word, '.DEBUG ')) then
         aoomod_debug = .true.
         aoomod       = .true.
         x2c_add_amfi = 2
      end if
      if (kw_matches(word, '.atomic')) then
         aoomod       = .true.
         x2c_add_amfi = 2
      end if

      call check_whether_kw_found(word, kw_section)

#else
      call quit('Atomic oo-order spin-orbit correction not included in this version')
#endif

   end subroutine

   subroutine read_input_qcorr(word, kw_section)

      use input_reader
      use qcorr_cfg

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      character(len=400)               :: input_line
      integer                          :: i, j, la, lb
      integer                          :: ios, islash, ibar
!     --------------------------------------------------------------------------

#include "dgroup.h"

      luci_cfg_qcorr = .true.

      call reset_available_kw_list()

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, print_qcorr)
      end if

      if (kw_matches(word, '.CVORB ')) then
         call kw_read(word, cvorb_qcorr, 2)
      end if

      if (kw_matches(word, '.REFFIL')) then

        allocate(nash_qcorr(nfsym))

        open(1234,file='refvec.luci',status='old',form='unformatted',    &
        access='sequential',action="read",position='rewind')

        read(1234) nref_qcorr, nash_qcorr(1:nfsym)
        read(1234) nref_e_qcorr
        allocate(reference_energy_qcorr(nref_e_qcorr))
        reference_energy_qcorr = 0
        read(1234) reference_energy_qcorr(1:nref_e_qcorr)
!#define DEBUG_QCORR
#ifdef DEBUG_QCORR
        print *, 'reference energy/ies: #',nref_e_qcorr, '--> energies: ',reference_energy_qcorr(1:nref_e_qcorr)
#undef DEBUG_QCORR
#endif
!       read # reference configurations
        allocate(ref_wavefunction_qcorr(nref_qcorr,2))
        allocate(ref_wavefunction_coeff_qcorr(nref_qcorr))
        do i = 1, nref_qcorr
          read(1234) la,lb,ref_wavefunction_qcorr(i,1)(1:la), ref_wavefunction_qcorr(i,2)(1:lb),  &
                           ref_wavefunction_coeff_qcorr(i)

!#define DEBUG_QCORR
#ifdef DEBUG_QCORR
          print *, ref_wavefunction_qcorr(i,1)(1:la)
          print *, ref_wavefunction_qcorr(i,2)(1:lb)
          print *, ref_wavefunction_coeff_qcorr(i)
#undef DEBUG_QCORR
#endif

        end do

        close(1234,status="keep")
        deallocate(ref_wavefunction_qcorr)
        deallocate(ref_wavefunction_coeff_qcorr)
        deallocate(reference_energy_qcorr)
        deallocate(nash_qcorr)

      end if

      if(kw_matches(word, '.REFDAT')) then
        set_nash_q_qcorr = .false.
        allocate(nash_qcorr(nfsym))
        read(get_file_unit(), *) nash_qcorr(1:nfsym)
!       read # reference energies (if AOC-HF and more than 1 root)
        call kw_read(word, nref_e_qcorr)
        allocate(reference_energy_qcorr(nref_e_qcorr))
        reference_energy_qcorr = 0
        do i = 1, nref_e_qcorr
          call kw_read(word, reference_energy_qcorr(i))
        end do
!       read # reference configurations
        call kw_read(word, nref_qcorr)
        allocate(ref_wavefunction_qcorr(nref_qcorr,2))
        allocate(ref_wavefunction_coeff_qcorr(nref_qcorr))
!       read reference configurations + coefficients

!       process the reference configuration / coefficient
        do i = 1, nref_qcorr

          read(get_file_unit(),'(a)') input_line
          call upcase(input_line)

!         print *, 'ref input line is',input_line

          islash = index(input_line,'/')
          ibar   = index(input_line,'|')
          if(islash <= 1)then
            write(*,'(/a,i2,a,i2/a/2a)') 'ERROR for *QCORR .REFDAT no.',                                        &
                  i,'/',nref_qcorr,'- the defining input line does not contain a "/":',                         &
                  '- the bad line : ',input_line
            write(*,'(a)') '- Or the specification of the number of reference configurations might be wrong!'// &
                           ' Check the input.'
            call quit('Input error for .REFDAT under *QCORR')
          end if
          if(ibar <= 1)then
            write(*,'(/a,i2,a,i2/a/2a)') 'ERROR for *QCORR .REFDAT no.',                                        &
                  i,'/',nref_qcorr,'- the defining input line does not contain a "|":',                         &
                  '- the bad line : ',input_line
            write(*,'(a)') '- Or the specification of the number of reference configurations might be wrong!'// &
                           ' Check the input.'
            call quit('Input error for .REFDAT under *QCORR')
          end if

!         alpha string
          read(input_line(1:ibar-1),'(a)',iostat=ios) ref_wavefunction_qcorr(i,1)

          if(ios /= 0)then
            write(*,'(/a,i2,a,i2/a/2a)') 'ERROR for *QCORR .REFDAT no.',                                        &
                  i,'/',nref_qcorr,'- the defining input line does not contain an alpha string',                &
                  '- the bad line : ',input_line
            write(*,'(a)') '- Or the specification of the number of reference configurations might be wrong!'// &
                           ' Check the input.'
            call quit('Input error for .REFDAT under *QCORR')
          end if

!         beta string
          read(input_line(ibar+1:islash-1),'(a)',iostat=ios) ref_wavefunction_qcorr(i,2)

          if(ios /= 0)then
            write(*,'(/a,i2,a,i2/a/2a)') 'ERROR for *QCORR .REFDAT no.',                                        &
                  i,'/',nref_qcorr,'- the defining input line does not contain a beta string',                  &
                  '- the bad line : ',input_line
            write(*,'(a)') '- Or the specification of the number of reference configurations might be wrong!'// &
                           ' Check the input.'
            call quit('Input error for .REFDAT under *QCORR')
          end if

!         coefficient of reference configuration
          read(input_line(islash+1:),*,iostat=ios) ref_wavefunction_coeff_qcorr(i)

          if(ios /= 0)then
            write(*,'(/a,i2,a,i2/a/2a)') 'ERROR for *QCORR .REFDAT no.',                                        &
                  i,'/',nref_qcorr,'- the defining input line does not contain a coefficient',                  &
                  '- the bad line : ',input_line
            write(*,'(a)') '- Or the specification of the number of reference configurations might be wrong!'// &
                           ' Check the input.'
            call quit('Input error for .REFDAT under *QCORR')
          end if

        end do

        open(1234,file='refvec.luci',status='replace',form='unformatted',    &
        access='sequential',action="readwrite",position='rewind')

        write(1234) nref_qcorr, nash_qcorr(1:nfsym)
        write(1234) nref_e_qcorr
        write(1234) reference_energy_qcorr(1:nref_e_qcorr)
!       write # reference configurations
        do i = 1, nref_qcorr

          la = len(trim(ref_wavefunction_qcorr(i,1)))
          lb = len(trim(ref_wavefunction_qcorr(i,2)))

          write(1234) la,lb,ref_wavefunction_qcorr(i,1)(1:la), ref_wavefunction_qcorr(i,2)(1:lb),               &
                      ref_wavefunction_coeff_qcorr(i)

!#define DEBUG_QCORR
#ifdef DEBUG_QCORR
          print *, ref_wavefunction_qcorr(i,1)(1:la)
          print *, ref_wavefunction_qcorr(i,2)(1:lb)
          print *, ref_wavefunction_coeff_qcorr(i)
#undef DEBUG_QCORR
#endif

        end do

        close(1234,status="keep")
        deallocate(ref_wavefunction_qcorr)
        deallocate(ref_wavefunction_coeff_qcorr)
        deallocate(reference_energy_qcorr)
        deallocate(nash_qcorr)

      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

#ifdef HAS_PCMSOLVER
   subroutine read_input_pcm(word, kw_section)

      use dirac_cfg
      use pcmmod_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      dirac_cfg_pcm = .true.

      call reset_available_kw_list()

      if (kw_matches(word, '.SKIPSS')) then
!         Use precomputed small charges
         pcmmod_skipss = .true.
      end if

      if (kw_matches(word, '.SEPARA')) then
!         .SEPARATE option
!         Split potentials and polarization charges in nuclear and electronic
         pcmmod_separate = .true.
      end if

      if (kw_matches(word, '.DOSPF ')) then
!        Perform elimination of spin-orbit dependent terms in the PCM operator
         pcmmod_dospf = .true.
      end if

      if (kw_matches(word, '.SKIPOI')) then
!         .SKIPOIT option
!         Do not compute One-Index Transformed ASC (only for debugging)
         pcmmod_skipoit = .true.
      end if

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, pcmmod_print)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_pcmsolver(word, kw_section)

      use pcmmod_cfg
      use input_reader

      implicit none

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      pcmmod_host_provides_input = .true.

      call reset_available_kw_list()

      if (kw_matches(word, '.CAVTYP')) then
         call kw_read(word, pcmmod_cavity_type)
      end if

      if (kw_matches(word, '.PATCHL')) then
         call kw_read(word, pcmmod_patch_level)
      end if

      if (kw_matches(word, '.COARSI')) then
         call kw_read(word, pcmmod_coarsity)
      end if

      if (kw_matches(word, '.AREATS')) then
         call kw_read(word, pcmmod_cavity_area)
      end if

      if (kw_matches(word, '.MINDIS')) then
         call kw_read(word, pcmmod_min_distance)
      end if

      if (kw_matches(word, '.DERORD')) then
         call kw_read(word, pcmmod_der_order)
      end if

      if (kw_matches(word, '.NOSCAL')) then
              pcmmod_scaling = .false.
      end if

      if (kw_matches(word, '.RADIIS')) then
         call kw_read(word, pcmmod_radii_set)
      end if

      if (kw_matches(word, '.RESTAR')) then
         call kw_read(word, pcmmod_restart_name)
      end if

      if (kw_matches(word, '.MINRAD')) then
         call kw_read(word, pcmmod_min_radius)
      end if

      if (kw_matches(word, '.SOLVER')) then
         call kw_read(word, pcmmod_solver_type)
      end if

      if (kw_matches(word, '.SOLVNT')) then
         call kw_read(word, pcmmod_solvent)
      end if

      if (kw_matches(word, '.EQNTYP')) then
         call kw_read(word, pcmmod_equation_type)
      end if

      if (kw_matches(word, '.CORREC')) then
         call kw_read(word, pcmmod_correction)
      end if

      if (kw_matches(word, '.PROBER')) then
         call kw_read(word, pcmmod_probe_radius)
      end if

      if (kw_matches(word, '.GINSID')) then
         call kw_read(word, pcmmod_inside_type)
      end if

      if (kw_matches(word, '.GOUTSI')) then
         call kw_read(word, pcmmod_outside_type)
      end if

      if (kw_matches(word, '.EPSILO')) then
         call kw_read(word, pcmmod_outside_epsilon)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine
#endif

   subroutine read_input_dft(word, kw_section)

      use dft_cfg
      use sigma_prefactor_setting
      use input_reader

      implicit none

!     ---------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ---------------------------------------------------------------------------
      character(80)                    :: line
      character(4)                     :: first_4
      integer                          :: i, n
!     ---------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.TINYDE')) then
         call kw_read(word, dft_cfg_tinydens)
      end if

      if (kw_matches(word, '.SCREEN')) then
         call kw_read(word, dft_cfg_screening)
      end if

      if (kw_matches(word, '.GRAC  ')) then
!        gradient regulated asymptotic correction
         dft_cfg_grac = .true.
         read(get_file_unit(), '(a4)') first_4
         if      (lowercase(first_4) == 'lb94') then
            dft_cfg_asymptote_is_lb94 = .true.
         else if (lowercase(first_4) == 'lbal') then
            dft_cfg_asymptote_is_lbalpha = .true.
         else
            print *, 'keyword following .GRAC not recognized'
            stop
         end if
         call kw_read(word,               &
                      dft_cfg_grac_alpha, &
                      dft_cfg_grac_beta,  &
                      dft_cfg_ac_ip,      &
                      dft_cfg_ac_threshold)
      end if

      if (kw_matches(word, '.SAOP  ')) then
!        statistical averaging of (model) orbital potentials
         dft_cfg_saop = .true.
         read(get_file_unit(), '(a4)') first_4
         if      (lowercase(first_4) == 'lb94') then
            dft_cfg_asymptote_is_lb94 = .true.
         else if (lowercase(first_4) == 'lbal') then
            dft_cfg_asymptote_is_lbalpha = .true.
         else
            print *, 'keyword following .SAOP not recognized'
            stop
         end if
      end if

#ifdef PRG_DIRAC
      if (kw_matches(word, '.SAOP! ')) then
!        "original" SAOP as in JCP 112, 1344 (2000).
         dft_cfg_saop                    = .true.
         dft_cfg_saop_with_response_part = .true.
         dft_cfg_asymptote_is_lbalpha    = .true.
         dft_cfg_alda_hs                 = .true.
         dft_cfg_alda_ha                 = .true.
         call setaldahs(.true.) !activate alda for the h+ part
         call setaldaha(.true.) !activate alda for the h- part
      end if
#endif

      if (kw_matches(word, '.NOSDFT')) then
         dft_cfg_no_sdft = .true.
      end if

      if (kw_matches(word, '.COLLIN')) then
         dft_cfg_sdft_collinear = .true.
      end if

      if (kw_matches(word, '.BETASI')) then
         sigma_prefactor = -1.0d0
      end if

#ifdef PRG_DIRAC
      if (kw_matches(word, '.ALDA  ')) then
         dft_cfg_alda_hs = .true.
         dft_cfg_alda_ha = .true.
         call setaldahs(.true.) !activate alda for the h+ part
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.XALDA ')) then
!        use        (1-f)*S + VWN + f*HFx
!        instead of       S + VWN
         dft_cfg_xalda   = .true.
         dft_cfg_alda_hs = .true.
         dft_cfg_alda_ha = .true.
         call setxalda(.true.)
         call setaldahs(.true.) !activate alda for the h+ part
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.ALDA+ ')) then
!        use alda only for the      hermitian part
         dft_cfg_alda_hs = .true.
         call setaldahs(.true.) !activate alda for the h+ part
      end if

      if (kw_matches(word, '.ALDA- ')) then
!        use alda only for the anti-hermitian part
         dft_cfg_alda_ha = .true.
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.XALDA+')) then
!        use xalda only for the      hermitian part
         dft_cfg_xalda   = .true.
         dft_cfg_alda_hs = .true.
         call setxalda(.true.)
         call setaldahs(.true.) !activate alda for the h+ part
      end if

      if (kw_matches(word, '.XALDA-')) then
!        use xalda only for the anti-hermitian part
         dft_cfg_xalda   = .true.
         dft_cfg_alda_ha = .true.
         call setxalda(.true.)
         call setaldaha(.true.) !activate alda for the h- part
      end if

      if (kw_matches(word, '.GAUNTS')) then
         call set_lscale_dft_gaunt()
      end if
#endif

!radovan: not documented on wiki, code not perfectly tested
      if (kw_matches(word, '.BLOCKE')) then
         dft_cfg_blocked = .true.
      end if

      if (kw_matches(word, '.POINTW')) then
         dft_cfg_pointwise = .true.
      end if

#ifdef MOD_UNRELEASED
      if (kw_matches(word, '.OVLDIA')) then
         dft_cfg_overlap_diagnostic = .true.
      end if
#endif

      call check_whether_kw_found(word, kw_section)

   end subroutine


! reading subroutine for the standalone testing program of diagonalization
! routines (miro ilias)
     subroutine read_input_diagonalization_tests(word, kw_section)
#ifdef MOD_UNRELEASED
      use input_reader
      use diagmod
      implicit none
!     ---------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ---------------------------------------------------------------------------
      call reset_available_kw_list()

     if (kw_matches(word, '.TITLE ')) then
         call kw_read(word, title_text)
     end if

     if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, print_level)
     end if

     if (kw_matches(word, '.MTXFIL')) then
         call kw_read(word, matrix_file_name)
     end if

     if (kw_matches(word, '.DSYEVR')) then
        do_dsyevr = .true.
     end if

     ! Ulf's simple routine
     if (kw_matches(word, '.ULF   ')) then
        do_dsyevr_ulf = .true.
     end if

     if (kw_matches(word, '.RS    ')) then
        do_rs = .true.
     end if

     if (kw_matches(word, '.RSJACO')) then
        do_rsjaco = .true.
     end if

     if (kw_matches(word, '.QJACO ')) then
        do_qjaco = .true.
     end if

     if (kw_matches(word, '.PAUL  ')) then
        do_householder_psb = .true.
     end if

     if (kw_matches(word, '.EIGVTS')) then
        do_eigv_check = .true.
     end if

     call check_whether_kw_found(word, kw_section)
#else
     call quit('read_input_diagonalization_tests not in release version')
#endif  /* MOD_UNRELEASED */
   end subroutine
