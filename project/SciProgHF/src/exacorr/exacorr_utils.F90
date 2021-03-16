module exacorr_utils

!  Collection of simple utility functions used in exacorr

  implicit none
  private

  public relative_error
  public print_date
  public get_free_fileunit
  public print_exacorr_logo
  public print_iteration
#ifdef VAR_MPI
  public write_tensor
#endif
  public quicksort
  public print_exainput
  public print_exaoutput
  public shift_orbital_energy
  public print_orbital_energy

  contains

    real(8) function relative_error(quantity_1,quantity_2)
    ! Gives percentage error for values > 1, otherwise absolute error

      real(8),intent(in):: quantity_1,quantity_2
      real(8)           :: error

      error = quantity_1 - quantity_2

      if (abs(quantity_1).gt.1.0) then
        relative_error = abs (error / quantity_1)
      else
        relative_error = abs (error)
      end if

    end function relative_error
       
    subroutine print_date(text)
    ! Prints time and date with a string
      character(len=10) :: datex
      character(len=8)  :: timex
      character(len=*)  :: text
      integer           :: lentext
      lentext = len(text)
      call daytime(datex,timex)
      write (*,'(1x,a10,1x,a8,1x,a)') datex,timex,text(1:lentext)
    end subroutine print_date

    subroutine get_free_fileunit(unit)
    ! idea from http://fortranwiki.org, will become obsolete once fortran2008 is supported everywhere
    ! gets free file unit
      integer, intent(out) :: unit
      integer, parameter   :: LUN_MIN=100, LUN_MAX=1000
      logical :: opened
      integer :: lun
      unit=-1
      do lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        if (.not. opened) then
          unit=lun
          exit
        end if
      end do
    end subroutine get_free_fileunit

    subroutine print_exacorr_logo()
    ! Prints logo and authors of exacorr

      print*, ""
      print*, ""
      print*, "                   *****************************************************"
      print*, "                   ***     Entering the ExaCorr module in DIRAC      ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Authors:         - Lucas Visscher            ***"
      print*, "                   ***                   - Anastasios Papadopoulos   ***"
      print*, "                   ***  Contributors:    - Johann Pototschnig        ***"
      print*, "                   ***                   - Loic Halbert              ***"
      print*, "                   ***                   - Andre S. P. Gomes         ***"
      print*, "                   ***                   - Michal Repisky            ***"
      print*, "                   ***  Features:        - CCD/CC2/CCSD              ***"
      print*, "                   ***                   - 1DM                       ***"
      print*, "                   ***                   - ReSpect interface         ***"
      print*, "                   ***                                               ***"
      print*, "                   ***  Relativistic Coupled Cluster code using the  ***"
      print*, "                   ***  TALSH and ExaTensor libraries developed by   ***"
      print*, "                   ***  Dmitry Lyakh                                 ***"
      print*, "                   *****************************************************"
      print*, ""
      print*, ""

    end subroutine print_exacorr_logo

    subroutine print_exainput(exa_input,nao)

      use exacorr_datatypes,   only : exacc_input

      type(exacc_input), intent(in) :: exa_input
      integer, intent(in) :: nao
      integer             :: i

      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,'(A)') "      Characteristics of this calculation      "
      write(*,'(A)') "-----------------------------------------------"
      write(*,'(A,L1)') " - Solve CCD equations :                ", exa_input%ccd
      write(*,'(A,L1)') " - Solve Lambda equations :             ", exa_input%lambda
      write(*,'(A,L1)') " - Compute perturbative Triples :       ", exa_input%do_triples
      write(*,'(A,I5)') " - Integral transformation scheme : ", exa_input%moint_scheme
      write(*,'(A,I5)') " - Maximum number of cycles :       ", exa_input%ncycles
      write(*,'(A,ES10.1)') " - Target precision of amplitudes :  ", exa_input%t_econv
      if (exa_input%moint_scheme.eq.42) then
        write(*,'(A,ES10.1)') " - Cholesky decomposition threshold :",exa_input%t_cholesky
      end if
      if (exa_input%level_shift.gt.1.0D-14) then
        write(*,'(A,ES10.1)') " - Level Shift :                      ",exa_input%level_shift
        write(*,'(A)')        " WARNING: MP2 energies are wrong due to level shift"
      end if
      if (exa_input%nff(1).gt.0) then
        write(*,'(A,L1)') " - Adding finite fields                 "
        do i=1,exa_input%nff(1)
          write(*,'(A,A8)') " - - Property :                         ",exa_input%ff_names(i)
          write(*,*) "- - Field strengths :                  ",exa_input%ff(i,1:exa_input%nff(2))
        end do
      end if
      write(*,*) "" 
      write(*,'(A,I5)') " - Number of atomic orbitals :      ", nao
      write(*,'(A,I5)') " - Number of occupied spinors :     ", exa_input%nocc
      write(*,'(A,I5)') " - Number of virtual spinors :      ", exa_input%nvir
      write(*,*) ""
      if (exa_input%talsh) then 
        write(*,'(A)') " Using the single-node TALSH library "
        write(*,'(A,I5,A)') " - TALSH buffer size:               ",exa_input%talsh_buff," GB"
      else 
        write(*,'(A)') " Using the multi-node ExaTensor library (MPI) "
        write(*,'(A,I5)') " - ExaTensor block size :           ", exa_input%exa_blocksize
        write(*,'(A)') "!! ATTENTION !! Properties require TALSH, density matrix must be smaller than: "
        write(*,'(A,I5,A)') " - TALSH buffer size:              ",exa_input%talsh_buff," GB"
      end if
      write(*,*) ""
      write(*,'(A,I5)') " - Print Level :                    ", exa_input%print_level
      write(*,*) ""
      write(*,'(A)') " - Estimating memory requirements"
      write(*,'(A,I9,A)') " -- nao  ( scheme 1-3  ) : ",3*16*((nao*nao)/1000)**2/1000," GB" 
      write(*,'(A,I9,A)') " -- nvir ( scheme 4/42 ) : ",3*16*((exa_input%nvir*exa_input%nvir)/1000)**2/1000," GB"
      write(*,'(A)') "-----------------------------------------------"
      write(*,*) ""

    end subroutine print_exainput

    subroutine print_exaoutput(talsh, ccd, cc2, scf_energy, mp2_energy, cc_energy, t1diag, t_energy)

      logical, intent(in)           :: talsh, ccd, cc2
      real(8), intent(in)           :: scf_energy, mp2_energy, cc_energy
      real(8), intent(in), optional :: t1diag
      real(8), intent(in), optional :: t_energy(3)
      
      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,'(A)') "     -     Final results from EXACORR     -    "
      write(*,'(A)') "-----------------------------------------------"
      if (talsh) then 
        write(*,'(A)') " Using the single-node TALSH library "
      else 
        write(*,'(A)') " Using the multi-node ExaTensor library (MPI) "
      end if
      write(*,*) ""
      write (*,'(A,F20.15)') "  MP2 energy :                      ",mp2_energy
      if (ccd) then 
        write (*,'(A,F20.15)') "  Final  CCD energy :               ",cc_energy
      else 
        if (cc2) then
          write (*,'(A,F20.15)') "  Final  CC2 energy :               ",cc_energy
        else
          write (*,'(A,F20.15)') "  Final CCSD energy :               ",cc_energy
          if (present(t_energy)) then 
            write (*,'(A,F20.15)') "  4th order triples correction :    ",t_energy(1)
            write (*,'(A,F20.15)') "  5th order triples (T) correction :",t_energy(2)
            write (*,'(A,F20.15)') "  5th order triples -T  correction :",t_energy(3)
            write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",cc_energy+t_energy(1)
            write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",cc_energy+t_energy(1)+t_energy(2)
            write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",cc_energy+t_energy(1)+t_energy(3)
          end if
        end if
        if (present(t1diag)) then
          write (*,'(A,F20.15)') "  T1 diagnostic      :              ",t1diag
        end if
      end if
      write(*,*) ""
      write (*,'(A,F25.15)') "  SCF energy :                 ",scf_energy
      write (*,'(A,F25.15)') "  Total MP2 energy :           ",scf_energy+mp2_energy
      if (ccd) then 
        write (*,'(A,F25.15)') "  Total  CCD energy :          ",scf_energy+cc_energy
      else 
        if (cc2) then
          write (*,'(A,F25.15)') "  Total  CC2 energy :          ",scf_energy+cc_energy
        else
          write (*,'(A,F25.15)') "  Total CCSD energy :          ",scf_energy+cc_energy
          if (present(t_energy)) then 
            write (*,'(A,F25.15)') "  Total CCSD+T  energy :       ",scf_energy+cc_energy+t_energy(1)
            write (*,'(A,F25.15)') "  Total CCSD(T) energy :       ",scf_energy+cc_energy+t_energy(1)+t_energy(2)
            write (*,'(A,F25.15)') "  Total CCSD-T  energy :       ",scf_energy+cc_energy+t_energy(1)+t_energy(3)
          end if
        end if
      end if
      write(*,*) ""
      write(*,'(A)') "-----------------------------------------------"
      write(*,*) ""

    end subroutine print_exaoutput

    subroutine print_iteration(iteration,convergence,energy,print_level)
    ! Formatted print for iterations

      ! input variables
      integer, intent(in)           :: iteration
      real(8), intent(in)           :: convergence
      real(8), intent(in), optional :: energy
      integer, intent(in), optional :: print_level

      if (iteration.eq.1 .and. present(energy)) then
        write(*,*) ""
        write(*,*) "-----------------------------------------------"
        write(*,*) "  It.     Energy                  Convergence  "
        write(*,*) "-----------------------------------------------"
        write(*,'(A3,I2,A4,F20.15,A3,E10.1)') "   ", iteration, "    ", energy, "   ", convergence
      else if (iteration.eq.1) then
        write(*,*) ""
        write(*,*) "------------------------"
        write(*,*) "  It.      Convergence  "
        write(*,*) "------------------------"
        write(*,'(A3,I2,A4,E10.1)') "   ", iteration, "    ", convergence
      else if (iteration.ne.1 .and. present(energy)) then
        write(*,'(A3,I2,A4,F20.15,A3,E10.1)') "   ", iteration, "    ", energy, "   ", convergence
      else if (iteration.ne.1) then
        write(*,'(A3,I2,A4,E10.1)') "   ", iteration, "    ", convergence
      end if
      
      if (present(print_level)) then
        if (print_level.gt.4) then
          call print_date(' CC iteration done')
        end if
      end if

    end subroutine print_iteration

#ifdef VAR_MPI
    subroutine write_tensor(exa_tensor, i_dims)

      use exatensor
      use talsh
      use tensor_algebra

      type(tens_rcrsv_t)              :: exa_tensor
      integer(INTL), dimension(2)     :: i_dims
      integer(INTD), dimension(2)     :: h_dims
      type(talsh_tens_t)              :: talsh_tensor
      complex(8), pointer, contiguous :: talsh_tens(:,:)
      type(C_PTR)                     :: body_p
      integer                         :: i,j
      integer(INTD)                   :: ierr
      integer(C_SIZE_T)               :: buf_size=1_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
      integer(C_INT)                  :: host_arg_max

      ierr=talsh_init(buf_size,host_arg_max)
      call print_date('Initialized talsh library')
      write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

      h_dims=i_dims
      write(*,*) '+++ write to file'
      ierr=talsh_tensor_construct(talsh_tensor,C8,h_dims,init_val=(0.D0,0.D0))
      write(*,*) 'construct err:', ierr
      ierr=exatns_tensor_get_slice(exa_tensor,talsh_tensor)
      write(*,*) 'get_slice err:', ierr
      ierr=talsh_tensor_get_body_access(talsh_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p,talsh_tens,h_dims) 
      do i=1,h_dims(2)
          write(*,*) (talsh_tens(i,j), j=1,h_dims(1))
      enddo
      ierr=talsh_tensor_destruct(talsh_tensor)

      ierr = talsh_shutdown()
      
    end subroutine write_tensor
#endif

    recursive subroutine quicksort(mo_list, i_beg, i_end)

        implicit none

        integer,intent(inout) :: mo_list(:)
        integer, intent(in) :: i_beg, i_end
        integer :: i,j
        integer :: ref, tem

        i=floor((float(i_end)+float(i_beg))/2)
        ref=mo_list(i)
        i=i_beg
        j=i_end

        do
          do while (mo_list(i)<ref)
            i=i+1
          end do
          do while (mo_list(j)>ref)
            j=j-1
          end do
          if (i<j) then
            tem=mo_list(i)
            mo_list(i)=mo_list(j)
            mo_list(j)=tem
            i=i+1
            j=j-1
          else if(i.eq.j) then
            i=i+1
            EXIT
          else
            EXIT
          end if
        end do

        if (i_beg<j) call quicksort(mo_list, i_beg, j)
        if (i<i_end) call quicksort(mo_list, i, i_end)

    end subroutine quicksort

    subroutine shift_orbital_energy(eps_vir,eps_occ,level_shift)
          implicit none

          real(8), intent(inout)   :: eps_vir(:),eps_occ(:)
          real(8), intent(in)      :: level_shift
          
          integer            :: print_level=1
          integer            :: lumo,i
          real(8)            :: e_homo, e_lumo
          real(8), parameter :: THRESHOLD=1.D-8
          
          e_homo = maxval(eps_occ)
          lumo   = minloc(eps_vir,1)


          if (level_shift.gt.1.0D-14) then
            e_lumo = minval(eps_vir)

            do i=1,size(eps_vir)
              eps_vir(i)=eps_vir(i)+level_shift
              
              if (print_level.gt.-1) then
                  write(*,*) "LUMO shifted from  ", e_lumo," to ",eps_vir(lumo)
              end if

            end do
          end if
          
          if (eps_vir(lumo)-e_homo<THRESHOLD) then
            write(*,*) 'WARNING: Negative HOMO-LUMO gap (.LSHIFT possible)'
          end if

    end subroutine shift_orbital_energy

    subroutine print_orbital_energy(eps_occ,nocc,eps_vir,nvir,info)
      implicit none

      real(8), intent(in)   :: eps_vir(:),eps_occ(:)
      integer, intent(in)   :: nocc,nvir,info

      integer :: i
      
      if (info==1) then
        write(*,*) ' --- shifted active occupied spinors --- '
      else
        write(*,*) ' --- recomputed active occupied spinors --- '
      end if 
      write(*,*) '   #             E   '
      do i=1,nocc
        write(*,'(I8,E18.9)') i, eps_occ(i)
      end do

      if (info==1) then
        write(*,*) ' --- shifted active virtual spinors --- '
      else
        write(*,*) ' --- recomputed active virtual spinors --- '
      end if
      write(*,*) '   #             E   '
      do i=1,nvir
        write(*,'(I8,E18.9)') i, eps_vir(i)
      end do

    end subroutine print_orbital_energy

end module exacorr_utils
