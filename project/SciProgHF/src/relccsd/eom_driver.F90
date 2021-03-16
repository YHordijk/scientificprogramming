  module eom_driver
  use projectors
  use sigma_eom
  use trial_vectors_eom
  use relcc_cfg
  implicit none
!this module calculates the sigma vectors for IP(2h-1p), EA(1h-2p) and EE(2h-2p)  

!    public eom_ip  
!    public eom_ea
    public eom_ee
    public number_of_states

    private

    character(8) :: eomtype
    integer :: nroot,state_sym
    integer :: irpoff = 0
 

   contains

  subroutine set_continuous_orbital(active_irrep)
#include "symm.inc"
#include "dcbham.h"
#include "dgroup.h"
    integer, intent(in) :: active_irrep

!Assigning Dimension to the continuous orbital
!===============================================
     ncont(1:nrep) = 0
     ncont(active_irrep) = 1

  end subroutine

  subroutine verify_input()
#include "symm.inc"
#include "dcbham.h"
#include "dgroup.h"
#include "files.inc"
    logical :: reset_roots_in_symmetry = .false.
    
    if (    (relcc_do_eomip .and. relcc_do_eomea) &
        .or.(relcc_do_eomip .and. relcc_do_eomee) &
        .or.(relcc_do_eomea .and. relcc_do_eomee)) &
          call quit('Specifying more than one EOM kind of calculation at a time (IP/EE/EA)  not allowed !') 
    if (relcc_do_eomip) eomtype = 'IP'
    if (relcc_do_eomea) eomtype = 'EA'
    if (relcc_do_eomee) then
       eomtype = 'EE'
       irpoff = nrep
     endif

     do state_sym = 1, nrep 
        if (relcc_eom_nroots(state_sym).ne.0) then
! checking whether or not, for ip/ea, we can actually ask for this many roots 
!
! we choose to restrict the number of roots one cas ask for in the eom-ee/ip/ea calculationsnos to
! avoid setting up guess trial vectors which would have to sample the 2h1p (ip, ea) or 2h2p subblocks:
!
! 1. IP: maximum number of roots for irrepN is reset back to the number of occupied spinors (nocc) in irrepN
! 2. EA: maximum number of roots for irrepN is reset back to the number of virtual  spinors (nvirt), per symmetry 
! 3. EE: maximum number of roots for irrepN is reset back to  nocc(irrepM)*nvirt(irrepO) for irrepN = irrepM * irrepO

           write (iw,'(2x,A,I6,2x,A,2x,A4,A1,I2,A1)') "Requested",relcc_eom_nroots(state_sym),"root(s) for symmetry: ", &
                & REPNA(irpoff+state_sym),"(",state_sym,")'"
          select case (eomtype)
             case ('IP')
                if (relcc_eom_nroots(state_sym).gt.no(state_sym)) then
                   relcc_eom_nroots(state_sym) = no(state_sym)
                   reset_roots_in_symmetry = .true.
                else if (relcc_eom_nroots(state_sym).eq.-1) then
                   reset_roots_in_symmetry = .true.
                else
                   reset_roots_in_symmetry = .false.
                end if
             case ('EA')
                if (relcc_eom_nroots(state_sym).gt.nv(state_sym)) then
                    relcc_eom_nroots(state_sym) = nv(state_sym)
                   reset_roots_in_symmetry = .true.
                else if (relcc_eom_nroots(state_sym).eq.-1) then
                   reset_roots_in_symmetry = .true.
                else
                   reset_roots_in_symmetry = .false.
                end if
             case ('EE')
                if (relcc_eom_nroots(state_sym).gt.nvo(state_sym)) then
                    relcc_eom_nroots(state_sym) = nvo(state_sym)
                   reset_roots_in_symmetry = .true.
                else if (relcc_eom_nroots(state_sym).eq.-1) then
                   reset_roots_in_symmetry = .true.
                else
                   reset_roots_in_symmetry = .false.
                end if
          end select

          if (reset_roots_in_symmetry) then
             if(relcc_do_eomip) then
                write (iw,'(4x,A)') "Warning: more EOM-IP roots than there are occupied orbitals."
             endif
             if(relcc_do_eomea) then
                write (iw,'(4x,A)') "Warning: more EOM-EA roots than there are virtual orbitals."
             endif
             if(relcc_do_eomee) then
                write (iw,'(4x,A)') "Warning: more EOM-EE roots than there are occupied/virtual pairs."
             endif
             write (iw,'(13x,A,I3,2x,A,2x,I3)') "Number or roots for symmetry",state_sym,"set to",relcc_eom_nroots(state_sym)
          end if
       end if
    end do

  end subroutine

   subroutine input_processor_right(Fbar_HH,Fbar_PP,Vbar_SS,R1,R2,H_II)
      use davidson 

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"

!description
!===========
!1)  this subroutine processes the guess vectors for eomcc.
!2)  it also defines the Preconditioner for the error vectors, which is relevant at
!    the Davidson step.

!calling variables
!================

   real*8, intent(in)  :: Fbar_HH(:),Fbar_PP(:),Vbar_SS(:)
   real*8, intent(inout) :: R1(:,:), R2(:,:) ! (ndimt1,nroot) and (ndimt2,nroot)
   real*8, intent(inout),target :: H_II(:)

!Local Variables
!==============

   integer :: i,ii,irep,j,iaai,aoff,a,aa,ai,arep,iroot,k(1),jroot,newarep
   real*8  :: fac1
   real*8  :: eps((iv(nrep+1)+io(nrep+1))) 
   real*8, pointer :: diag_H(:)
   integer :: icount, acount

   integer :: abij, ijrp, abrp, jrp, jj, irp, ioff, imin, bb, amin, arp, brp, b, iarep, airep
   real(kind=8) :: fac, fac3, fac2
   integer :: ndimr1, ndimr2, search_dimension, xroot, N1, N2, N
   real*8,allocatable :: R1_dum(:,:) 
   integer :: fudge_factor
   real*8 :: in_window_HbarSS_threshold
   integer :: number_diagonal_elements_in_window
   logical, allocatable :: index_in_window(:)
   integer :: ai_norcw

!initialize the R-vector
!=======================

    N1 = size(R1,1)
    N2 = size(R2,1)
    ndimr1 = N1 / rcw
    ndimr2 = N2 / rcw

    R1(1:N1,1:nroot) = 0.0d0 
    R2(1:N2,1:nroot) = 0.0d0 

! aspg, WIP

    allocate(index_in_window(ndimr1))
    index_in_window = .false.
    number_diagonal_elements_in_window = 0
    in_window_HbarSS_threshold = huge(in_window_HbarSS_threshold)
! aspg, WIP

! aspg, WIP

! 1> Determine the minimum energy determinants for IP, EA & EE. The number of them is equal to the number of roots. 

     i = 1
     ii = 1
     do irep = 1, nrep
        if (no(irep) .gt. 0) call dcopy (no(irep),Fbar_HH(ii),((no(irep)+1)*rcw),eps(i),1)
        i  =  i + no(irep)
        ii =  ii + NO(irep)* no(irep) * rcw
     enddo

     ii = 1
     do irep = 1, nrep
        if (nv(irep) .gt. 0) call dcopy (nv(irep),Fbar_PP(ii),((nv(irep)+1)*rcw),eps(i),1)
        i  =  i + nv(irep)
        ii =  ii + NV(irep)* nv(irep) * rcw
     enddo

   select case (eomtype)

   case ("EE")
! aspg, WIP
! whereas for IP/EA one can use the spinor energies as a guide for which
! diagonal element of Hbar relates to the core-core part of the diagonal, 
! for EE we have to look at the orbital energy differences and therefore
! use a helper variable to map between the core-core part of the diagonal
! and the occupied spinor energy threshold given as input 
   in_window_HbarSS_threshold = huge(in_window_HbarSS_threshold)
   ai_norcw = 1
! aspg, WIP

    ai = 1
    ii = 0
    do irep = 1, nrep
        arep = multb (state_sym+nrep,irep,2)

        iarep = multb(irep,arep,2)
        airep = multb(arep,irep,2)

        do i = 1, no(irep)
         ii = ii + 1
         fac1 = eps(ii)
         aoff = io(nrep+1) + iv(arep)
          do a = 1, nv(arep)
             aa = aoff + a
              iaai = (iovvo(state_sym) + (iivo(arep,irep)+(i-1)*nv(arep)+a-1)*nov(iarep) + &
       &      iiov(irep,arep) + (a-1)*no(irep)+ i - 1)*rcw+1   
              H_II(ai) = -fac1 + eps(aa)! - vbar_ss (iaai)
! aspg, WIP
! whereas for IP/EA one can use the spinor energies as a guide for which
! diagonal element of Hbar relates to the core-core part of the diagonal, 
! for EE we have to look at the orbital energy differences and therefore
! use a helper variable to map between the core-core part of the diagonal Hbar(ai,ai)
! and the occupied spinor energy threshold given as input. We will use this
! variable to determine the trial vectors below 
              if (((eps(ii).le.relcc_projectors_rew_occ_max_energy).and.&
                 & (eps(ii).ge.relcc_projectors_rew_occ_min_energy).and.&
                 & ((eps(aa).ge.relcc_projectors_rew_virt_min_energy).and.&
                 &  (eps(aa).le.relcc_projectors_rew_virt_max_energy))) &
                 ) then
                 index_in_window(ai_norcw) = .true. 
                 if (H_II(ai).lt.in_window_HbarSS_threshold) then
                    in_window_HbarSS_threshold = H_II(ai)
                 end if
                 number_diagonal_elements_in_window = number_diagonal_elements_in_window + 1
              end if
! aspg, WIP
              ai_norcw = ai_norcw + 1
              ai = ai + rcw
          enddo
        enddo
    enddo

!aspg+as
      ABIJ = ai 
      DO IJRP = 1, NREP
       abrp = multb (state_sym+nrep,ijrp+nrep,2)
      DO JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) CYCLE
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            IF (ARP.LT.BRP) CYCLE
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + rcw
               ENDDO
            ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
!aspg+as

      if (relcc_mfd_trial_full_matrix) then
! New trial vectors for R1 
          do i = 1, min(ndimr1,nroot)
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1(iroot,jroot) = 1.0d0
          end do
          if (nroot.gt.ndimr1) then
              do i = 1, min(ndimr2,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
                  R2(iroot,ndimr1+jroot) = 1.0d0
             end do
          end if

      else if (relcc_mfd_trial_ccs) then
          write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
          allocate (R1_dum(N1,ndimr1))
          R1_dum = 0.0d0

          do i = 1, ndimr1
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1_dum(iroot,jroot) = 1.0d0
          end do
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector         ")

          call trial_vectors_ee(R1_dum,Fbar_HH,Fbar_PP,Vbar_SS)
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector  CD1    ",complex_dimension=1)

          R1 = R1_dum(:,1:nroot)
          if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")

          deallocate(R1_dum)


      else if (relcc_mfd_trial_diagonal) then
          write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a
! particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

          allocate(diag_H(ndimr1))
          diag_H = H_II(1:N1:rcw)

! aspg, WIP june/july 2018
! in principle we can eigenvalue shift with any of the schemes below 
          if ((relcc_mfd_eigenvalues_energy_shift(1).gt.(0.0d0))) then
              do i = 1, size(diag_H,1)
                 diag_H(i) = abs(diag_H(i) - relcc_mfd_eigenvalues_energy_shift(1))
              end do
          end if
! aspg, WIP june/july 2018

          if (relcc_projectors_do_restricted_excitation_window &
             & .or.relcc_projectors_do_core_valence_separation) then ! find the lowest elements of the diagonal within the given excitation window

              if (iprnt.ge.1) then
                  write (iw,*) ''
                  if (relcc_projectors_do_restricted_excitation_window) then 
                      write (iw,*) '   Restricted excitation window (REW) via projectors will be invoked'
                      write (iw,*) '      - 1h1p determinants with all (quasi)particles within' 
                      write (iw,*) '        the occupied/virtual windows will be retained'
                  else if (relcc_projectors_do_core_valence_separation) then
                      write (iw,*) '   Core valence separation (CVS) via projectors will be invoked'
                      write (iw,*) '      - 1h1p determinants with all particles within' 
                      write (iw,*) '        the occupied window will be retained'
                  end if
                  if (relcc_projectors_rew_strict) then
                      write (iw,*) '      - 2h2p/... determinants with all (quasi)particles within' 
                      write (iw,*) '        the occupied/virtual windows will be retained'
                  else if (relcc_projectors_rew_remove_double_occupied) then
                      write (iw,*) '      - 2h2p/... determinants with only one (quasi)particle within'            
                      write (iw,*) '        the occupied window will be retained'
                  else
                      write (iw,*) '      - 2h2p/... determinants with at least one (quasi)particle within' 
                      write (iw,*) '        the occupied/virtual windows will be retained'
                  end if
                  write (iw,*) '   Spinors/orbitals in the following window(s) will be active'
                  write (iw,*) '      occupied min     ',relcc_projectors_rew_occ_min_energy,'a.u.'
                  write (iw,*) '               max     ',relcc_projectors_rew_occ_max_energy,'a.u.'
                  write (iw,*) '      virtual  min     ',relcc_projectors_rew_virt_min_energy,'a.u.'
                  write (iw,*) '               max     ',relcc_projectors_rew_virt_max_energy,'a.u.'
                  write (iw,*) '   There are',number_diagonal_elements_in_window,'possible 1h1p determinants within the window'
                  write (iw,*) ''
              end if
              if (number_diagonal_elements_in_window.lt.nroot) then
                  call quit ('Requested more roots than there are 1h1p determinants within the window. Stopping.')
              end if

              if (iprnt.ge.2) write (iw,*) 'diag(Hbar_SS): ',diag_H
              if (iprnt.ge.2) write (iw,*) 'in_window_HbarSS_threshold',in_window_HbarSS_threshold
! aspg, WIP june/july 2018
              iroot = 1
              do while (iroot.le.nroot)
                 k = minloc(diag_H, MASK=(diag_H .ge. in_window_HbarSS_threshold))
                 i = k(1)
                 if (iprnt.ge.2) write(iw,*)'roots',iroot,diag_H(i),i
                 if (iprnt.ge.4) write(iw,*) 'making',i,'huge'
                 diag_H(i) = huge(diag_H)
! we may have diagonal elements which do not correspond to an element within the windows, so we check that
! this is the case
                 if (index_in_window(i).eqv..true.) then
                     if (iprnt.ge.2) write(iw,*)'determinant in window',i
                     jroot = i
                     if (rcw.eq.2) jroot = 2*i - 1
                     R1(jroot,iroot) = 1.0d0 
                     iroot = iroot + 1
                 end if
              end do

              if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")
! aspg, WIP june/july 2018

          else ! this is the standard way of doing it, find the lowest elements of the diagonal

              do iroot = 1, nroot
                  k = minloc(diag_H)
                  if (iprnt.ge.3) write(iw,*)'roots',iroot,diag_H(k(1)),k(1)
                  i = k(1)
                  jroot = i
                  if (rcw.eq.2) jroot = 2*i - 1
                  R1(jroot,iroot) = 1.0d0
                  diag_H(i) = huge(diag_H)
              enddo
          end if
          deallocate(diag_H)

      else if (relcc_mfd_trial_restart) then
          call quit ('restart from prior RHS vectors not implemented')
      else
          call quit ('unknown RHS trial vectors generation scheme')
      end if ! if chain for trial vectors 


   case ('IP')

!As of now, I am targetting the lowest energy determinants across the symmetries.

    ai_norcw = 1
    ai = 1
    ii = 0

    do irep = 1, nrep
        arep = multb (1+nrep,irep,2)
        do i = 1, no(irep)
         ii = ii + 1
         fac1 = eps(ii)
          do a = 1, ncont(arep)
             H_II(ai) = -fac1 
! aspg, WIP
              if ( ((eps(ii).le.relcc_projectors_rew_occ_max_energy).and.&
                 &  (eps(ii).ge.relcc_projectors_rew_occ_min_energy)) &
                 ) then
                 index_in_window(ai_norcw) = .true.
                 if (H_II(ai).lt.in_window_HbarSS_threshold) then
                    in_window_HbarSS_threshold = H_II(ai)
                 end if
                 number_diagonal_elements_in_window = number_diagonal_elements_in_window + 1
              end if
! aspg, WIP
             ai_norcw = ai_norcw + 1
             ai = ai + rcw
          enddo

        enddo
    enddo

!      H_II = -eps(1:(ndimr1*rcw)) !+ 100.0d0

      ABIJ = ai 
      DO IJRP = 1, NREP
       abrp = multb (1+nrep,ijrp+nrep,2)
      DO JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) CYCLE
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = EPS(II) + FAC1
            DO BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               DO A = 1, NCONT(ARP)
                  FAC = FAC3! - 100.0d0
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + rcw
               ENDDO
            ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      if (relcc_mfd_trial_full_matrix) then
! New trial vectors for R1 
          do i = 1, min(ndimr1,nroot)
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1 
              jroot = i 
              R1(iroot,jroot) = 1.0d0
          end do
          if (nroot.gt.ndimr1) then
              do i = 1, min(ndimr2,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1 
                  jroot = i 
                  R2(iroot,ndimr1+jroot) = 1.0d0
             end do
          end if

      else if (relcc_mfd_trial_ccs) then
          write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
          allocate (R1_dum(N1,ndimr1))
          R1_dum = 0.0d0

          do i = 1, ndimr1
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1_dum(iroot,jroot) = 1.0d0
          end do
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector         ")

          call trial_vectors_ip(R1_dum,Fbar_HH)
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector  CD1    ",complex_dimension=1)
      
! aspg 19/04/2017 : this if structure should't be here, but it appears that the
!                   unit trial vectors constructed for irrep 1 in effect span other irreps (at least for water with DC)
!                   is this due to issues setting the continuous orbitals for other irreps ?
!         if (nrep.eq.1) then
             R1 = R1_dum(:,1:nroot) 
!         else
!            jroot = 1
!            do i = 1, nroot
!               R1(:,i) = R1_dum(:,jroot)
!               jroot = jroot + 2
!            end do
!         end if
          if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")

          deallocate(R1_dum)

      else if (relcc_mfd_trial_diagonal) then
          write (iw,*) '  Using (reordered) unit vectors as trial vectors'

          allocate(diag_H(ndimr1))
          diag_H = H_II(1:N1:rcw)

! aspg, WIP june/july 2018
! in principle we can eigenvalue shift with any of the schemes below 
          if ((relcc_mfd_eigenvalues_energy_shift(1).gt.(0.0d0))) then
              do i = 1, size(diag_H,1)
                 diag_H(i) = abs(diag_H(i) - relcc_mfd_eigenvalues_energy_shift(1))
              end do
          end if
! aspg, WIP june/july 2018


          if (relcc_projectors_do_restricted_excitation_window &
             & .or.relcc_projectors_do_core_valence_separation) then ! find the lowest elements of the diagonal within the given excitation window

              if (iprnt.ge.1) then
                  write (iw,*) ''
                  if (relcc_projectors_do_restricted_excitation_window) then
                      write (iw,*) '   Restricted ionization window (RIW) via projectors will be invoked'
                      write (iw,*) '      - 1h determinants with all (quasi)particles within'
                      write (iw,*) '        the occupied/virtual windows will be retained'
                  else if (relcc_projectors_do_core_valence_separation) then
                      write (iw,*) '   Core valence separation (CVS) via projectors will be invoked'
                      write (iw,*) '      - 1h  determinants with all particles within'
                      write (iw,*) '        the occupied window will be retained'
                  end if
                  if (relcc_projectors_rew_strict) then
                      write (iw,*) '      - 2h1p/... determinants with all (quasi)particles within'
                      write (iw,*) '        the occupied/virtual windows will be retained'
                  else if (relcc_projectors_rew_remove_double_occupied) then
                      write (iw,*) '      - 2h1p/... determinants with only one (quasi)particle within'
                      write (iw,*) '        the occupied window will be retained'
                  else
                      write (iw,*) '      - 2h1p/... determinants with at least one (quasi)particle within'
                      write (iw,*) '        the occupied/virtual windows will be retained'
                  end if
                  write (iw,*) '   Spinors/orbitals in the following window(s) will be active'
                  write (iw,*) '      occupied min     ',relcc_projectors_rew_occ_min_energy,'a.u.'
                  write (iw,*) '               max     ',relcc_projectors_rew_occ_max_energy,'a.u.'
                  write (iw,*) '      virtual  min     ',relcc_projectors_rew_virt_min_energy,'a.u.'
                  write (iw,*) '               max     ',relcc_projectors_rew_virt_max_energy,'a.u.'
                  write (iw,*) '   There are',number_diagonal_elements_in_window,'possible 1h determinans within the window'
                  write (iw,*) ''
              end if

              if (number_diagonal_elements_in_window.lt.nroot) then
                  call quit ('Requested more roots than there are 1h determinants within the window. Stopping.')
              end if

              if (iprnt.ge.2) write (iw,*) 'diag(Hbar_SS): ',diag_H
              if (iprnt.ge.2) write (iw,*) 'in_window_HbarSS_threshold',in_window_HbarSS_threshold
! aspg, WIP june/july 2018
              iroot = 1
              do while (iroot.le.nroot)
                 k = minloc(diag_H, MASK=(diag_H .ge. in_window_HbarSS_threshold))
                 if (iprnt.ge.2) write(iw,*)'roots',iroot,diag_H(k(1)),k(1)
                 if (iprnt.ge.3) write(iw,*) diag_H
                 i = k(1)
                 if (iprnt.ge.3) write(iw,*) 'making',i,'huge'
                 diag_H(i) = huge(diag_H)
! we may have diagonal elements which do not correspond to an element within the windows, so we check that
! this is the case
                 if (index_in_window(i).eqv..true.) then
                     jroot = i
                     if (rcw.eq.2) jroot = 2*i - 1
                     R1(jroot,iroot) = 1.0d0
                     iroot = iroot + 1
                 end if
              end do

              if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")
! aspg WIP june/july 2018

          else ! this is the standard way of doing it, find the lowest elements of the diagonal

              do iroot = 1, nroot
                  k = minloc(diag_H) 
                  if (iprnt.ge.1) write(iw,*)'roots',iroot,diag_H(k(1)),k(1) 
                  i = k(1)
                  diag_H(i) = huge(diag_H)
                  jroot = i
                  if (rcw.eq.2) jroot = 2*i - 1
                  R1(jroot,iroot) = 1.0d0
              enddo
          end if
          deallocate(diag_H)

      else if (relcc_mfd_trial_restart) then
          call quit ('restart from prior RHS vectors not implemented')
      else
          call quit ('unknown RHS trial vectors generation scheme')
      end if ! if chain for trial vectors 


   case ("EA")

    ai = 1
    do irep = 1, nrep
       arep = multb (irep,1+nrep,1)
        do i = 1, ncont(irep)
         aoff = io(nrep+1) + iv(arep)
          do a = 1, nv(arep)
             aa = aoff + a
             H_II(ai) =  eps(aa) 
             ai = ai + rcw
          enddo
        enddo
    enddo

      ABIJ = ai 
      DO IJRP = 1, NREP
       abrp = multb (ijrp+nrep,1+nrep,1)
      DO JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = EPS(JJ)
         DO I = 1, NCONT(IRP)
            FAC2 =  FAC1
            DO BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,ABRP+NREP,2)
            IF (ARP.LT.BRP) CYCLE
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = FAC2 - EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = FAC3 - EPS(AA)
                  H_II(ABIJ) = -FAC
                  ABIJ = ABIJ + rcw
               ENDDO
            ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      if (relcc_mfd_trial_full_matrix) then
! New trial vectors for R1 
          do i = 1, min(ndimr1,nroot)
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1(iroot,jroot) = 1.0d0
          end do
          if (nroot.gt.ndimr1) then
              do i = 1, min(ndimr2,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
                  R2(iroot,ndimr1+jroot) = 1.0d0
             end do
          end if

      else if (relcc_mfd_trial_ccs) then
          write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
          allocate (R1_dum(N1,ndimr1))
          R1_dum = 0.0d0

          do i = 1, ndimr1
              iroot = i
              if (rcw.eq.2) iroot = 2*i - 1
              jroot = i
              R1_dum(iroot,jroot) = 1.0d0
          end do
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector         ")

          call trial_vectors_EA(R1_dum,Fbar_PP)
          if (iprnt.ge.3) call davidson_print_vectors(R1_dum, label="debug:  R1_dum vector  CD1    ",complex_dimension=1)

! aspg 19/04/2017 : this if structure should't be here, but it appears that the
!                   unit trial vectors constructed for irrep 1 in effect span other irreps (at least for water with DC)
!                   is this due to issues setting the continuous orbitals for other irreps ?
!         if (nrep.eq.1) then
             R1 = R1_dum(:,1:nroot)
!         else
!            jroot = 1
!            do i = 1, nroot
!               R1(:,i) = R1_dum(:,jroot)
!               jroot = jroot + 2
!            end do
!         end if
          if (iprnt.ge.3) call davidson_print_vectors(R1, label="debug:  R1 vector             ")

          deallocate(R1_dum)

! calculate the orbital dimension within which our desired root should belong..

      else if (relcc_mfd_trial_diagonal) then
          write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

          search_dimension = ncont(state_sym)*nv(state_sym)

          allocate(diag_H(ndimr1))
          diag_H = H_II(1:N1:rcw)

! aspg, WIP june/july 2018
          if ((relcc_mfd_eigenvalues_energy_shift(1).gt.(0.0d0))) then
              do i = 1, size(diag_H,1)
                 diag_H(i) = abs(diag_H(i) - relcc_mfd_eigenvalues_energy_shift(1))
              end do
          end if
! aspg, WIP june/july 2018

          fudge_factor = 2
! here we set the factor to 1 because for C1, with only one fermion irrep,
! we can't get around to working with the degenerate roots
          if (nrep.eq.1) fudge_factor = 1

          do iroot = 1, nroot
              k = minloc(diag_H) 
              if (iprnt.ge.2) write(iw,*)'roots',iroot,diag_H(k(1)),k(1) 
              i = k(1)
              diag_H(i) = huge(diag_H)
              jroot = i
              if (rcw.eq.2) jroot = 2*i - 1
              R1(jroot,iroot) = 1.0d0
          enddo
 
          deallocate(diag_H)

      else if (relcc_mfd_trial_restart) then
          call quit ('restart from prior RHS vectors not implemented')
      else
          call quit ('unknown RHS trial vectors generation scheme')
      end if ! if chain for trial vectors 

  end select 
 
   deallocate(index_in_window)

   end subroutine

   subroutine input_processor_left(Fbar_HH,Fbar_PP,Vbar_SS,eVectorsR,L1,L2,H_II)

  use symmetry_offset
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"


!description
!===========
!1)  this subroutine processes the left hand guess vectors for eomcc.
!2)  it also defines the Preconditioner for the error vectors, which is relevant at
!    the Davidson step.

!calling variables
!================

   real*8, intent(in)  :: Fbar_HH(:),Fbar_PP(:),Vbar_SS(:),eVectorsR(:,:) 
   real*8, intent(inout) :: L1(:,:), L2(:,:) ! (ndimt1,nroot) and (ndimt2,nroot)
   real*8, intent(inout),target :: H_II(:)

!Local Variables
!==============

   integer :: i,ii,irep,j,iaai,aoff,a,aa,ia,arep,iroot,k(1)
   real*8  :: fac1
   real*8  :: eps((iv(nrep+1)+io(nrep+1))) 
   real*8, pointer :: diag_H(:)
   integer :: icount, acount

   integer :: ijab, ijrp, abrp, jrp, jj, irp, ioff, imin, bb, amin, arp, brp, b
   real(kind=8) :: fac, fac3, fac2
   integer :: ndimL1, ndimL2, search_dimension, xroot, N, N1, N2, fudge_factor, jroot
   type(Offset)         :: e
   real*8,allocatable :: L1_dum(:,:) 


!initialize the R-vector
!=======================

    N1 = size(L1,1)
    N2 = size(L2,1)
    N  = N1 + N2

    ndimL1 = N1 / rcw
    ndimL2 = N2 / rcw

!   L1(1:N1,1:nroot) = 0.0d0 
!   L2(1:N2,1:nroot) = 0.0d0 
    L1 = 0.0d0 
    L2 = 0.0d0 

! 1> store the diagonal elements of the transformed Hamiltonian..  

     i = 1
     ii = 1
     do irep = 1, nrep
        if (no(irep) .gt. 0) call dcopy (no(irep),Fbar_HH(ii),((no(irep)+1)*rcw),eps(i),1)
        i  =  i + no(irep)
        ii =  ii + NO(irep)* no(irep) * rcw
     enddo

     ii = 1
     do irep = 1, nrep
        if (nv(irep) .gt. 0) call dcopy (nv(irep),Fbar_PP(ii),((nv(irep)+1)*rcw),eps(i),1)
        i  =  i + nv(irep)
        ii =  ii + NV(irep)* nv(irep) * rcw
     enddo

   select case (eomtype)

   case ("EE")

    ia = 1
    ii = 0

      do arp = 1, nrep
         irp = multb (state_sym+nrep,arp,2)
         aa = io(nrep+1) + iv(arp)
      do a = 1, nv(arp)
         aa = aa + 1
         fac1 = eps(aa)
         ioff = io(irp)
         do i = 1, no(irp)
             ii = ioff + i
             h_ii(ia) = fac1 - eps(ii) ! I will add the missing two-electron term later 
             ia = ia + 1
         enddo
      enddo
      enddo


      ijab = ia
      do abrp = 1, nrep
         ijrp = multb (state_sym+nrep,abrp,2)
      do brp = 1, nrep
      bb = iv(brp) + io(nrep+1)
      arp = multb(brp,abrp+nrep,2)
      if (arp.lt.brp) cycle
      do b = 1, nv(brp)
         bb = bb + 1
         fac1 = eps(bb)
            aoff = iv(arp) + io(nrep+1)
           amin = 1
           if (arp.eq.brp) amin = b + 1
         do a = amin, nv(arp)
            aa = aoff + a
            fac2 = eps(aa) + fac1
            do jrp = 1, nrep
            jj = io(jrp)
            irp = multb(jrp,ijrp+nrep,2)
            if (irp.lt.jrp) cycle
            ioff = io(irp)
            do j = 1, no(jrp)
               jj = jj + 1
               fac3 = fac2 - eps(jj)
               imin = 1
               if (irp.eq.jrp) imin = j + 1
               do i = imin, no(irp)
                  ii = ioff + i
                  fac = fac3 - eps(ii)
                  h_ii(ijab) = fac
                  ijab = ijab + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo


    do iroot = 1, nroot

      call srt1c1 (nrep,nv,no,eVectorsR(1:N1,iroot), &
     &             L1(1:N1,iroot)) 
      call srt1c1 (nrep,nvvt,noot,eVectorsR(N1+1:N,iroot), &
     &             L2(1:N2,iroot)) 

    enddo

!         allocate(diag_H(N1))
!         diag_H = H_II(1:N1:rcw)

!         do iroot = 1, nroot
!             k = minloc(diag_H)
!             if (iprnt.ge.3) write(iw,*)'roots',iroot,diag_H(k(1)),k(1)
!             i = k(1)
!             jroot = i
!             if (rcw.eq.2) jroot = 2*i - 1
!             L1(jroot,iroot) = 1.0d0
!             diag_H(i) = 8000.0d0
!         enddo
!         deallocate(diag_H)


   case ('IP')

    ia = 1
    ii = 0

      do arp = 1, nrep
         irep = multb (arp,1+nrep,1)
         aa = io(nrep+1) + iv(arp)
      do a = 1, ncont(arp)
         ioff = io(irep)
         do i = 1, no(irep)
             ii = ioff + i
             h_ii(ia) = -eps(ii)  
             ia = ia + 1
         enddo
      enddo
      enddo


      ijab = ia
      do abrp = 1, nrep
      do brp = 1, nrep
      bb = iv(brp) + io(nrep+1)
      arp = multb(brp,abrp+nrep,2)
      do b = 1, nv(brp)
         bb = bb + 1
         fac1 = eps(bb)
         do a = 1, ncont(arp)
            fac2 =  fac1
            do jrp = 1, nrep
            jj = io(jrp)
            irp = multb(jrp,abrp+nrep,2)
            if (irp.lt.jrp) cycle
            ioff = io(irp)
            do j = 1, no(jrp)
               jj = jj + 1
               fac3 = fac2 - eps(jj)
               imin = 1
               if (irp.eq.jrp) imin = j + 1
               do i = imin, no(irp)
                  ii = ioff + i
                  fac = fac3 - eps(ii)
                  h_ii(ijab) = fac
                  ijab = ijab + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo

    call alloc_array(e,nrep)
    call auto_symmetry_offset(e,nv,ncont,.false.,.false.)

     allocate(diag_H(ndiml1))
     diag_H = H_II(1:ndiml1:rcw)

        xroot = 0 

       do iroot = 1, nroot
        k = minloc(diag_H) 
        i = k(1)
        diag_H(i) = 8000.0d0

         xroot = xroot + 1
         jroot = i
         if (rcw.eq.2) jroot = 2*i - 1

         L1(jroot,xroot) = 1.0d0
       enddo
     deallocate(diag_H)

      if (relcc_mfd_trial_lhs_use_rhs) then
         write (iw,*) '  Using RHS eigenvectors as LHS trial vectors'
! if we are reusing the right-hand side eigenvectors, we resort them

!         do iroot = 1, nroot
!             call srt1c1 (nrep,ncont,no,eVectorsR(1:N1,iroot:iroot), &
!    &                     L1(1:N1,iroot:iroot))
!             call srt1c1 (nrep,e%oneNonDirac,noot,eVectorsR(N1+1:N,iroot:iroot), &
!    &                     L2(1:N2,iroot:iroot))
!         enddo

      else
! otherwise, we are going to generate trial vectors for the lhs
          if (relcc_mfd_trial_full_matrix) then
! New trial vectors for L1 
              do i = 1, min(ndiml1,nroot)
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
!                 print *, 'iroot',iroot,' jroot',jroot
                  L1(iroot,jroot) = 1.0d0
              end do
              if (nroot.gt.ndiml1) then
                  do i = 1, min(ndiml2,nroot)
                      iroot = i
                      if (rcw.eq.2) iroot = 2*i - 1
                      jroot = i
                      L2(iroot,ndiml1+jroot) = 1.0d0
                 end do
              end if

          else if (relcc_mfd_trial_ccs) then
              write (iw,*) '  Using eigenvalues of singles-singles block as trial vectors'
              allocate (L1_dum(N1,ndiml1))
              L1_dum = 0.0d0

              do i = 1, ndiml1
                  iroot = i
                  if (rcw.eq.2) iroot = 2*i - 1
                  jroot = i
                  L1_dum(iroot,jroot) = 1.0d0
              end do
              if (iprnt.ge.3) call davidson_print_vectors(L1_dum, label="debug: L1_dum vector          ")

! this is a temporary hack (as of 31/03/2017), since the routine below only returns the right-hand vectors
! to generalize
              call trial_vectors_ip(L1_dum,Fbar_HH) 
              if (iprnt.ge.3) call davidson_print_vectors(L1_dum, label="debug: L1_dum vector  CD1     ",complex_dimension=1)

              L1 = L1_dum(:,1:nroot)

              if (iprnt.ge.3) call davidson_print_vectors(L1, label="debug:  R1 vector             ")

              deallocate(L1_dum)

          else if (relcc_mfd_trial_diagonal) then
              write (iw,*) '  Using (reordered) unit vectors as trial vectors'

! In the following section we generate the guess vectors belonging to a
! particular state symmetry.
! Thereby we avoid the sigma vector calculations of the degenerate roots.  

              search_dimension = ncont(state_sym)*no(state_sym)

              allocate(diag_H(N1))
              diag_H = H_II(1:N1)

          xroot = 0 

          do iroot = 1, nroot
              k = minloc(diag_H) 
              if (iprnt.ge.2) write(iw,*)'roots',iroot,diag_H(k(1)),k(1) 
              i = k(1)
              diag_H(i) = 8000.0d0
                   xroot = xroot + 1
                   jroot = i
                   if (rcw.eq.2) jroot = 2*i - 1
                   L1(jroot,xroot) = 1.0d0
          enddo
          deallocate(diag_H)

          else if (relcc_mfd_trial_restart) then
              call quit ('restart from prior LHS vectors not implemented')
          else
              call quit ('unknown LHS trial vectors generation scheme')
          end if ! if chain for trial vectors 

      end if ! relcc_mfd_trial_lhs_use_rhs

  call dealloc_array(e)

   case ("EA")

    ia = 1
    aa = 0

      do arp = 1, nrep
         irep = multb (arep,state_sym+nrep,1)
         aa = io(nrep+1) + iv(arp)
      do a = 1, nv(arp)
         aa = aa + 1
         fac1 = eps(aa)
         do i = 1, ncont(arp)
             h_ii(ia) = fac1  
             ia = ia + 1
         enddo
      enddo
      enddo

      ijab = ia
      do abrp = 1, nrep
      do brp = 1, nrep
      bb = iv(brp) + io(nrep+1)
      arp = multb(brp,abrp+nrep,2)
      if (arp.lt.brp) cycle
      do b = 1, nv(brp)
         bb = bb + 1
         fac1 = eps(bb)
            aoff = iv(arp) + io(nrep+1)
           amin = 1
           if (arp.eq.brp) amin = b + 1
         do a = amin, nv(arp)
            aa = aoff + a
            fac2 = eps(aa) + fac1
            do jrp = 1, nrep
            jj = io(jrp)
            irp = multb(jrp,abrp+nrep,2)
            do j = 1, no(jrp)
               jj = jj + 1
               fac3 = fac2 - eps(jj)
               do i = 1, ncont(irp)
                  fac = fac3  
                  h_ii(ijab) = -fac
                  ijab = ijab + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo

! calculate the orbital dimension within which our desired root should belong..

  search_dimension = ncont(state_sym)*nv(state_sym) 

     allocate(diag_H(N1))
     diag_H = H_II(1:N1)
        xroot = 1
       do iroot = 1, 2*nroot
        k = minloc(diag_H) 
        i = k(1)
        diag_H(i) = 8000.0d0
        if (i <= search_dimension) then
         L1(i,xroot) = 1.0d0
         xroot = xroot + 1
        endif
       enddo
     deallocate(diag_H)
!! end if 

  end select 

   end subroutine


   subroutine eom_ee(icalc,mxcorr,eps)
   use davidson
   use memory_allocator
   use allocator_parameters, only : klongint
   use eom_density
   use symmetry_offset 

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dgroup.h"
!--------------------------------------------------
 integer, intent(in)                   :: icalc
 integer(kind=klongint), intent(inout) :: mxcorr
 real*8                                :: eps(*)
!--------------------------------------------------
! what is ndt here should be the dimension of t1+t2 for the irrep in question
   type(intermediates) :: B
   real*8, allocatable, target :: buf1(:),buf2(:),buf3(:) ! may be not all of them are necessary. 
   type(davidson_results) :: davResultsR(nrep), davResultsL(nrep)
   integer :: irp

!memory allocation
!------------------

  integer(kind=8)        :: current_mem_use = 0
  integer(kind=8)        :: start_mem_use   = 0
  integer                :: freemem = 0
  integer                :: work_memsize
  integer                :: nbuf1,nbuf2
  integer                :: nbuf3
  integer                :: state
  integer*8              :: nv5irp,nv6irp,nvirp

    call allocator_get_words_inuse(start_mem_use)
    call legacy_lwork_get(work_memsize)

    nbuf1 = max(joooo(nrep+1),jvooo(nrep+1),jvvoo(nrep+1),jvovo(nrep+1))
    nbuf2 = nbuf1

    call alloc( buf1  , nbuf1*rcw             , id="buf1" )
    call alloc( buf2  , nbuf2*rcw             , id="buf2" )

!   For the moment I am using as much in-core memory as possible. However we
!   may think of writing them in files, replacing the integral ones...  
! aspg: we should collect these actions into a routine, like
!    call initialize_intermediates() 
    call alloc(B%Fbar_mi,   (nfoo*rcw), id="Fbar_mi")  
    call alloc(B%Fbar_me,   (nfvo*rcw), id="Fbar_me")  
    call alloc(B%Fbar_ae,   (nfvv*rcw), id="Fbar_ae")  

    call alloc(B%W_ijmn,    (nv1*rcw),            id="W_ijmn") 
    call alloc(B%W_ejmb,    (iovvo(nrep+1)*rcw),  id="W_ejmb")
    call alloc(B%W_mbej,    (iovvo(nrep+1)*rcw),  id="W_mbej")
    call alloc(B%W_iemn,    (iovoot(nrep+1)*rcw), id="W_iemn")

    call alloc(B%W_efam,    (nv5*rcw),            id="W_efam")
    call alloc(B%W_amef,    (nv5*rcw),            id="W_amef")

    call alloc(B%W_mnie,    (iovoot(nrep+1)*rcw), id="W_mnie")
    call alloc(B%Wbar_mbej, (iovvo(nrep+1)*rcw),  id="Wbar_mbej")

    call alloc(B%t1, (ndimt1*rcw), id="t1")
    call alloc(B%t2, (ndimt2*rcw), id="t2")
! end of intermediates, now passing on to the last local buffer
       
    call allocator_get_maxbuff(freemem,kind(BUF3))

    NBUF3 = freemem / RCW

    NV5IRP = MAX(JVOVV(NREP+1),JOVVV(NREP+1))
    NV6IRP = 1
    DO IRP = 1, NREP
       NV6IRP = MAX(NV6IRP,INT(NVVT(IRP),8)*INT(NVVT(IRP),8))
    ENDDO

!   Check the minimum size for the general use of this array

    NVIRP = MAX(NV5IRP,NV6IRP)

    NBUF3 = MIN(NVIRP,INT(NBUF3,8))
    NBUF3 = MAX(NBUF1,NBUF3)

!     Allocate now also this remaining work array

    call alloc( BUF3  , NBUF3*RCW             , id="buf3" )


    call allocator_get_words_inuse(current_mem_use)
    mxcorr = current_mem_use - start_mem_use

    if (icalc.eq.1) then
       write (iw,1001) "EOMCC calculation",mxcorr
       goto 451
    endif

    write (iw,*) " "
    write (iw,*) "Equation of Motion CC module setup"
    write (iw,*) " "

    call verify_input()
    call initialize_buffers(buf1, buf2, buf3)
    call ccdump
    call getampt(B%t1,B%t2)
    call IM_container(eomtype,B)
    call free_buffers(buf1,buf2,buf3)

    do state_sym = 1, nrep ! 16 
       call eom_solve_for_left_or_right(state_sym, .false., .true., davResultsL, davResultsR, B, eps)

       if (relcc_do_eomprop) then
          call eom_solve_for_left_or_right(state_sym, .true., .false., davResultsL, davResultsR, B, eps)
          call eom_check_left_right_roots(state_sym, davResultsL,davResultsR)
          call eom_renormalize_birthogonal(state_sym, davResultsL, davResultsR)
        endif ! relcc_do_eomprop
   end do

! if the analysis is only being performed for the eigenvalues, and from above we ensured the RHS and LHS 
! have the same energy and are biorthogonal, we only need to get the energies from one, so we use the RHS
   write (iw,*) ""
   write (iw,'(2x,A)') "Final results for EOM-"//eomtype(1:2)
   write (iw,*) ""
   call eom_analysis_allirreps(davResultsR, eps)

! now we calculate all density and transition density matrices with the states at our disposal ?
   if (relcc_do_eomprop) call eom_calculate_store_density_matrices(B, davResultsL, davResultsR) 

   do state_sym = 1, nrep 
      if (relcc_eom_nroots(state_sym) <= 0) cycle
      call davidson_cleanup(davResultsR(state_sym))
      if (relcc_do_eomprop)  call davidson_cleanup(davResultsL(state_sym))
   end do

 451 continue

!
! final cleanup
!
     call free_intermediates(B)
     if (icalc.eq.1) then
        call free_buffers(buf1,buf2,buf3)
     end if

1001  FORMAT (" Core used for ",A," :",T50,I15," 8-byte words")
   end subroutine

   subroutine eom_check_left_right_roots(state_sym,davResultsL,davResultsR)
! calculate density matrices of sym target symmetry. which combinations of L and R ?
      use davidson
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dgroup.h"
!--------------------------------------------------
      integer :: iroot, i
      integer, intent(in) :: state_sym
      type(davidson_results) :: davResultsR(nrep), davResultsL(nrep)
      real*8 :: difference_in_evalue
      real*8 :: rhs_lhs_evalue(2)
      real*8 :: threshold
      double precision dznrm2


      if (relcc_eom_nroots(state_sym) <= 0) return

      write (iw,*) " "
      do iroot = 1, relcc_eom_nroots(state_sym)
         rhs_lhs_evalue(1) = davResultsR(state_sym)%eValues(1,iroot) - davResultsL(state_sym)%eValues(1,iroot)
         rhs_lhs_evalue(2) = davResultsR(state_sym)%eValues(2,iroot) - davResultsL(state_sym)%eValues(2,iroot)
         difference_in_evalue = dznrm2(1,rhs_lhs_evalue,1)
         write (iw,'(2X,A,I4,A,I3,2X,A,2X,E16.8)') &
                "INFO: difference between RHS and LHS eigenvalues for root ",iroot,&
     &          " irrep ",state_sym,":",difference_in_evalue
         write (iw,'(2X,A)') "INFO:              Re                  Im"
         write (iw,'(2X,A,2X,E16.8,4X,E16.8)') "INFO: E(RHS)   ",davResultsR(state_sym)%eValues(:,iroot)
         write (iw,'(2X,A,2X,E16.8,4X,E16.8)') "INFO: E(LHS)   ",davResultsL(state_sym)%eValues(:,iroot)
         write (iw,*) " "
      end do
   end subroutine

subroutine eom_renormalize_birthogonal(state_sym, davResultsL, davResultsR)
! calculate density matrices of sym target symmetry. which combinations of L and R ?
   use davidson
   use memory_allocator
   use allocator_parameters, only : klongint
   use eom_density
   use symmetry_offset

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dgroup.h"
!--------------------------------------------------
   type(davidson_results), intent(inout) :: davResultsR(nrep), davResultsL(nrep)
   integer, intent(in) :: state_sym

   integer :: iroot, i, j, irep
   integer :: ndimr1,ndimr2

   integer :: N, N1, N2 ! size of determinantal space 
   type(Offset)         :: e
   real*8               :: overlap

!memory allocation
!------------------
  real*8,allocatable     :: l1_temp(:),l2_temp(:)

           if (relcc_eom_nroots(state_sym) <= 0) return 

           do iroot = 1, relcc_eom_nroots(state_sym)
               select case(eomtype)
                  case ('IP')
!                    call set_continuous_orbital(state_sym)
                     ndimr1 = 0

                     do irep = 1, nrep
                         ndimr1 = ndimr1 + ncont(irep)*no(irep)
                     enddo

                     call alloc_array(e,nrep)
                     call auto_symmetry_offset(e,ncont,nv,.false.,.false.)

                     ndimr2 = 0

                     do irep = 1, nrep
                         ndimr2 = ndimr2 + e%oneNonDirac(irep)*noot(irep)
                     enddo

                     N1 = ndimr1*rcw
                     N2 = ndimr2*rcw
                     N  = N1 + N2

                     allocate(l1_temp(N1))
                     allocate(l2_temp(N2))

                     call srt1c1(nrep,no,ncont,davResultsL(state_sym)%eVectorsL(1:N1,iroot),l1_temp)
                     call srt1c1(nrep,noot,e%oneNonDirac,davResultsL(state_sym)%eVectorsL(1+N1:N,iroot),l2_temp)

                     call dealloc_array(e)

                     overlap = dot_product(l1_temp,davResultsR(state_sym)%eVectorsR(1:N1,iroot)) + &
    &                          dot_product(l2_temp,davResultsR(state_sym)%eVectorsR(1+N1:N,iroot))

                     overlap = 1.0d0/overlap

                     call dscal(N,overlap,davResultsL(state_sym)%eVectorsL(:,iroot),1)
                     call dscal(N1,overlap,l1_temp,1)
                     call dscal(N2,overlap,l2_temp,1)

                     deallocate(l1_temp,l2_temp)

                  case ('EE')

                     ndimr1 = ndimt1
                     ndimr2 = ndimt2

                     N1 = ndimr1*rcw
                     N2 = ndimr2*rcw
                     N  = N1 + N2

                     allocate(l1_temp(N1))
                     allocate(l2_temp(N2))

                     call srt1c1(nrep,no,nv,davResultsL(state_sym)%eVectorsL(1:N1,iroot),l1_temp)
                     call srt1c1(nrep,noot,nvvt,davResultsL(state_sym)%eVectorsL(1+N1:N,iroot),l2_temp)

                     overlap = dot_product(l1_temp,davResultsR(state_sym)%eVectorsR(1:N1,iroot)) + &
    &                          dot_product(l2_temp,davResultsR(state_sym)%eVectorsR(1+N1:N,iroot))

                     overlap = 1.0d0/overlap

                     call dscal(N,overlap,davResultsL(state_sym)%eVectorsL(1:N,iroot:iroot),1)
                     call dscal(N1,overlap,l1_temp,1)
                     call dscal(N2,overlap,l2_temp,1)

                     deallocate(l1_temp,l2_temp)

                 case ('EA')
                    call quit ('LHS EOM-EA Not yet implemented')

               end select

           enddo ! loop over roots 
   end subroutine

subroutine eom_calculate_store_density_matrices(B, davResultsL, davResultsR)
! calculate density matrices of sym target symmetry. which combinations of L and R ?
   use davidson
   use memory_allocator
   use allocator_parameters, only : klongint
   use eom_density
   use symmetry_offset

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "dgroup.h"
!--------------------------------------------------
   type(davidson_results), intent(inout) :: davResultsR(nrep), davResultsL(nrep)
   type(intermediates), intent(inout) :: B

   integer :: iroot, i, j, irep
   integer :: ndimr1,ndimr2,sym
   integer :: N, N1, N2 ! size of determinantal space 
   type(Offset)         :: e

!memory allocation
!------------------
  integer                :: state
  real*8,allocatable     :: dmo(:)

!!!!!
       write (iw,*) ""
       write (iw,*) "Determining EOM density matrices"
       write (iw,*) ""

       allocate(dmo(norbt*norbt*nz))

       dmo =0.0d0
       state =1
       do sym = 1, nrep

           if (relcc_eom_nroots(sym) <= 0) cycle

           do iroot = 1, relcc_eom_nroots(sym)
               select case(eomtype)
                  case ('IP')
!                    call set_continuous_orbital(state_sym)
                     ndimr1 = 0

                     do irep = 1, nrep
                         ndimr1 = ndimr1 + ncont(irep)*no(irep)
                     enddo

                     call alloc_array(e,nrep)
                     call auto_symmetry_offset(e,ncont,nv,.false.,.false.)

                     ndimr2 = 0

                     do irep = 1, nrep
                         ndimr2 = ndimr2 + e%oneNonDirac(irep)*noot(irep)
                     enddo
                     call dealloc_array(e)

                     N1 = ndimr1*rcw
                     N2 = ndimr2*rcw
                     N  = N1 + N2

                     call density_IP(davResultsL(sym)%eVectorsL(1:N1,iroot),&
  &                                  davResultsL(sym)%eVectorsL(1+N1:N,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1:N1,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1+N1:N,iroot),B%t1,B%t2,dmo)

                  case ('EE')

                     ndimr1 = ndimt1
                     ndimr2 = ndimt2

                     N1 = ndimr1*rcw
                     N2 = ndimr2*rcw
                     N  = N1 + N2

                     call density_EE(davResultsL(sym)%eVectorsL(1:N1,iroot),&
  &                                  davResultsL(sym)%eVectorsL(1+N1:N,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1:N1,iroot), &
  &                                  davResultsR(sym)%eVectorsR(1+N1:N,iroot),B%t1,B%t2,dmo)

                 case ('EA')
                    call quit ('LHS EOM-EA Not yet implemented')

               end select

               call store_eom_density ('EOM ',state,dmo)

               state = state + 1
           enddo ! loop over roots 
       enddo ! loop over symmetries
       deallocate(dmo)
   end subroutine

   subroutine number_of_states(state)
#include "symm.inc"

   integer, intent(out) :: state 
   integer              :: sym

   state = 0
   do sym = 1, nrep  

     state = state + relcc_eom_nroots(sym) 

   enddo

   end subroutine  

   function merge_arrays(dimension_to_merge, a,b)
       real(kind=8), pointer :: merge_arrays(:,:)
       real(kind=8), intent(in) :: a(:,:), b(:,:)
       integer, intent(in) :: dimension_to_merge

       integer :: Na, Nb, La, Lb, Ns, Ls

       Na = size(a,1)
       La = size(a,2)

       Nb = size(b,1)
       Lb = size(b,2)

       if ((dimension_to_merge .eq. 1) .and. (La .eq. Lb)) then
           Ns = Na + Nb
           Ls = La
           allocate(merge_arrays(Ns,Ls))
           merge_arrays(1:Na, :) = a(1:Na, :) 
           merge_arrays(Na+1:Ns, :) = b(1:Nb, :)
       else if ((dimension_to_merge .eq. 2) .and. (Na .eq. Nb)) then
           Ls = La + Lb
           Ns = Na
           allocate(merge_arrays(Ns,Ls))
           merge_arrays(:, 1:La) = a(:, 1:La) 
           merge_arrays(:, La+1:Ls) = b(:, 1:Lb)
       else
           call quit ('mismatch in array dimensions')
       end if

   end function

   subroutine dump_to_file(fileno, filename, matrix)
       integer, intent(in) :: fileno
       character(len=14), intent(in) :: filename
       real(kind=8) :: matrix(:,:)
       integer :: irec
  
       inquire (iolength = irec) matrix
       open (unit=fileno, file=filename, form="unformatted", access="direct", recl=irec)
       write(unit=fileno, rec=1) matrix
       close(unit=fileno)
   end subroutine


   subroutine eom_analysis(eigenvalues,eigenvectors,eps) 

! Analyse the R-vectors of EOM-IP,EE,EA. Based on FS_ANALYSIS1 routine by Luuk

#include "param.inc"
#include "inpt.inc"
#include "complex.inc"
#include "symm.inc"
#include "results.inc"
#include "ccpar.inc"
#include "files.inc"

!---------------Calling variables--------------------------------------

     real*8,intent(inout) :: eigenvalues(rcw,1:*),eigenvectors(*)
     real*8,intent(in)    :: eps(*)

!---------------Local variables -----------------------------------
      INTEGER :: IERR, ISIZE, select_eom,iroot
      REAL*8 :: E_OFF
      INTEGER :: nstate_sym(nrep,4)
!---------------Executable code--------------------------------------

   nstate_sym(1:nrep,1:4) = 0
   CALL PRSYMB(iw,'=',80,0)

   select case(eomtype)

    case('IP')

      nstate_sym(state_sym,3) = nroot
      WRITE(iw,'(/A,A4,A,I2,A/)') ' *** EXCITATIONS OF FERMION IRREP ',REPNA(irpoff+state_sym),"(",state_sym,") ***"
      select_eom = 3

    case('EA')

      nstate_sym(state_sym,2) = nroot
      WRITE(iw,'(/A,A4,A,I2,A/)') ' *** EXCITATIONS OF FERMION IRREP ',REPNA(irpoff+state_sym),"(",state_sym,") ***"
      select_eom = 2

    case('EE')

      nstate_sym(state_sym,4) = nroot
      WRITE(iw,'(/A,A4,A,I2,A/)') ' *** EXCITATIONS OF BOSON IRREP ',REPNA(irpoff+state_sym),"(",state_sym,") ***"
      select_eom = 4

  end select  

!     Initialize the reference energy

      E_OFF = ESCF + ECCSD

!     Print an ordered list of eigenvalues
      CALL PRSYMB(iw,'-',80,0)
       IF      (RCW.EQ.1) THEN
          WRITE (iw,1005)
          DO iroot = 1,nroot
           WRITE(iw,1010)   &
     &       iroot,eigenvalues(1,iroot),E_OFF+eigenvalues(1,iroot)
         ENDDO

       ELSE IF (RCW.EQ.2) THEN
         WRITE (iw,1007)
         DO iroot = 1,nroot
           WRITE(iw,1015)   &
     &       iroot,eigenvalues(1,iroot),E_OFF+eigenvalues(1,iroot),ECCSDIM+eigenvalues(2,iroot)
         ENDDO
       ELSE
         call quit('eom_anaysis: Unknown algebra (RCW) !')
       ENDIF
      CALL PRSYMB(iw,'-',80,0)
      
      call evector_analysis (state_sym,select_eom,nstate_sym,eps,iprnt,e_off,eigenvalues(1:rcw,1:nroot),eigenvectors)

      CALL PRSYMB(iw,'=',80,0)
 1005 FORMAT(//' Level   Abs eigenvalue      Total Energy   (atomic units)'/)
 1010 FORMAT(I5,F17.10,F20.12,'  ',F20.12)
 1007 FORMAT(//' Level   Abs eigenvalue      Total Energy       Total imag.ener.   (atomic units)'/)
 1015 FORMAT(I5,F17.10,F20.12,'  ',F20.12,' ',F20.15)
   end subroutine


   subroutine evector_analysis (state_sym,select_eom,nstate_sym,eps,ipr,e_off,ev_all,eigenvectors)

!---------------Description--------------------------------------------
!
!     Analyze the excited states and calculate transition moments
!           SELECT_EOM : selection of eomtype
!                  EPS : orbital energies
!                  IPR : print level
!                E_OFF : CCSD energy to be added to eigenvalue
!               EV_ALL : eigenvalues for all irreps
!         eigenvectors : R-Vectors for all irreps (packed)
!
!---------------Last modified------------------------------------------
!
!     Author : Avijit Shee from FS_ANALYSIS routine by Luuk.
!
!---------------Calling variables--------------------------------------

#include "complex.inc"
#include "symm.inc"

   real*8, intent(in) :: eps(*),e_off
   real*8, intent(in) :: ev_all(rcw,nroot),eigenvectors(rcw,*) 
   integer            :: state_sym  
   integer            :: select_eom,ipr,nstate_sym(nrep,4)

!---------------Common Blocks--------------------------------------

#include "files.inc"
#include "param.inc"
#include "ccpar.inc"

!---------------Local variables--------------------------------------

   integer       :: i_eigenvalue, j_eigenvalue,ndt
   character*100 :: label
   real*8        :: value,tres
   integer       :: i,ii,irep,j,aoff,a,aa,ai,arep,ssym
   integer       :: ijrp, abrp, jrp, jj, kk, irp, ioff, imin, bb, amin, arp, brp, b

!---------------Executable code--------------------------------------
 
!     Now print an analysis of the vectors

      if (ipr.le.1) then
         tres = 0.1d0
      elseif (ipr.le.2) then
         tres = 0.01d0
      else
         tres = 0.01d0
      end if

      i_eigenvalue = 0

      write (iw,1000) tres
      do ssym = 1,nrep
 
!     get the dimension of the space
 
         ndt = nstate_sym(ssym,select_eom)
         if (ndt.eq.0) cycle

!   Loop over the eigenvalues

     ii = 1
     do j_eigenvalue = 1, nstate_sym(ssym,select_eom) 
          i_eigenvalue = i_eigenvalue + 1 
!    I have seen crazy ouputs with large fock space calcs, put a limit
          if (j_eigenvalue.lt.100000) then
          write (iw,1010) repna(irpoff+ssym),j_eigenvalue,  &
    &             e_off+ev_all(1,i_eigenvalue),ev_all(1,i_eigenvalue)
!    Loop over the R-vectors and print the interesting ones

    select case(eomtype)

      case('EA')

      do irep=1,nrep
        arep = multb (irep,1+nrep,1)
        do i = 1, ncont(irep)
          aoff = io(nrep+1) + iv(arep)
         do a = 1, nv(arep)
          aa = aoff + a
          write(label,1001) repna(arep),(a+no(arep)),eps(aa)
            value = abs(eigenvectors(1,ii))
            if (carith) value = value + abs(eigenvectors(2,ii))
            if (carith .and. value.gt.tres) then
                write (iw,'(2f15.5,3x,a50)') & 
  &             eigenvectors(1,ii),eigenvectors(2,ii),label
            else if (value.gt.tres) then              
            write (iw,'(f15.5,3x,a50)') eigenvectors(1,ii),label
            end if
          ii = ii + 1
         enddo
        enddo
      enddo

      do irep=1,nrep
      ijrp = irep
      abrp = multb (ijrp+nrep,1+nrep,1)
      do jrp = 1, nrep
      jj = io(jrp)
      irp = multb(jrp,ijrp+nrep,2)
      do j = 1, no(jrp)
         jj = jj + 1
         do i = 1, ncont(irp)
            do brp = 1, nrep
            bb = iv(brp) + io(nrep+1)
            arp = multb(brp,abrp+nrep,2)
            if (arp.lt.brp) cycle
            aoff = iv(arp) + io(nrep+1)
            do b = 1, nv(brp)
               bb = bb + 1
               amin = 1
               if (arp.eq.brp) amin = b + 1
               do a = amin, nv(arp)
                  aa = aoff + a
                  write(label,1002) repna(jrp),j,eps(jj), &
                &          repna(arp),(a+no(arp)),eps(aa),&
                &          repna(brp),(b+no(brp)),eps(bb)

                  value = abs(eigenvectors(1,ii))
                  if (carith) value = value + abs(eigenvectors(2,ii))
                  if (carith .and. value.gt.tres) then
                      write (iw,'(2f15.5,3x,a70)') & 
            &             eigenvectors(1,ii),eigenvectors(2,ii),label
                  else if (value.gt.tres) then              
                  write (iw,'(f15.5,3x,a70)') eigenvectors(1,ii),label
                  end if
                  ii = ii + 1
               enddo
            enddo
            enddo
         enddo
      enddo
      enddo
      enddo
      case('IP')

    jj = 0
  do irep = 1, nrep
     arep = multb (irep,1+nrep,1)
      do i = 1, no(irep)
       jj = jj + 1
        do a = 1, ncont(arep)
          write(label,1001) repna(irep),i,eps(jj)
            value = abs(eigenvectors(1,ii))
            if (carith) value = value + abs(eigenvectors(2,ii))
            if (carith .and. value.gt.tres) then
                write (iw,'(2f15.5,3x,a50)') & 
  &             eigenvectors(1,ii),eigenvectors(2,ii),label
            else if (value.gt.tres) then              
            write (iw,'(f15.5,3x,a50)') eigenvectors(1,ii),label
            end if
          ii = ii + 1
        enddo
      enddo
  enddo


  do ijrp = 1, nrep
!       abrp = multb (ijrp+nrep,1+nrep,1)
       abrp = multb (1+nrep,ijrp+nrep,2)
    do jrp = 1, nrep
      jj = io(jrp)
      irp = multb(jrp,ijrp+nrep,2)
      if (irp.lt.jrp) cycle
      ioff = io(irp)
      do j = 1, no(jrp)
         jj = jj + 1
         imin = 1
         if (irp.eq.jrp) imin = j + 1
         do i = imin, no(irp)
            kk = ioff + i
            do brp = 1, nrep
            bb = iv(brp) + io(nrep+1)
            arp = multb(brp,abrp+nrep,2)
            aoff = iv(arp) + io(nrep+1)
            do b = 1, nv(brp)
               bb = bb + 1
               do a = 1, ncont(arp)

                  write(label,1003) repna(irp),i,eps(kk), &
                &          repna(jrp),j,eps(jj),&
                &          repna(brp),(b+no(brp)),eps(bb)

                  value = abs(eigenvectors(1,ii))
                  if (carith) value = value + abs(eigenvectors(2,ii))
                  if (carith .and. value.gt.tres) then
                      write (iw,'(2f15.5,3x,a70)') & 
            &             eigenvectors(1,ii),eigenvectors(2,ii),label
                  else if (value.gt.tres) then              
                  write (iw,'(f15.5,3x,a70)') eigenvectors(1,ii),label
                  end if
                  ii = ii + 1
               enddo
            enddo
            enddo
         enddo
      enddo
    enddo
  enddo

      case('EE')

    jj = 0
  do irep = 1, nrep
     arep = multb (state_sym+nrep,irep,2)
      do i = 1, no(irep)
       jj = jj + 1
          aoff = io(nrep+1) + iv(arep)
        do a = 1, nv(arep)
          aa = aoff + a
          write(label,1004)repna(irep),i,eps(jj), &
            repna(arep),(a+no(arep)),eps(aa)
            value = abs(eigenvectors(1,ii))
            if (carith) value = value + abs(eigenvectors(2,ii))
            if (carith .and. value.gt.tres) then
                write (iw,'(2f15.5,3x,a50)') & 
  &             eigenvectors(1,ii),eigenvectors(2,ii),label
            else if (value.gt.tres) then              
            write (iw,'(f15.5,3x,a50)') eigenvectors(1,ii),label
            end if
          ii = ii + 1
        enddo
      enddo
  enddo


  do ijrp = 1, nrep
       abrp = multb (state_sym+nrep,ijrp+nrep,2)
    do jrp = 1, nrep
      jj = io(jrp)
      irp = multb(jrp,ijrp+nrep,2)
      if (irp.lt.jrp) cycle
      ioff = io(irp)
      do j = 1, no(jrp)
         jj = jj + 1
         imin = 1
         if (irp.eq.jrp) imin = j + 1
         do i = imin, no(irp)
            kk = ioff + i
            do brp = 1, nrep
            bb = iv(brp) + io(nrep+1)
            arp = multb(brp,abrp+nrep,2)
            if (arp.lt.brp) cycle
            aoff = iv(arp) + io(nrep+1)
            do b = 1, nv(brp)
               bb = bb + 1
               amin = 1
               if (arp.eq.brp) amin = b + 1
               do a = amin, nv(arp)
                  aa = aoff + a
                  write(label,1005) repna(irp),i,eps(kk), &
                &          repna(jrp),j,eps(jj),&
                &          repna(arp),(a+no(arp)),eps(aa),&
                &          repna(brp),(b+no(brp)),eps(bb)

                  value = abs(eigenvectors(1,ii))
                  if (carith) value = value + abs(eigenvectors(2,ii))
                  if (carith .and. value.gt.tres) then
                      write (iw,'(2f15.5,3x,a100)') & 
            &             eigenvectors(1,ii),eigenvectors(2,ii),label
                  else if (value.gt.tres) then              
                  write (iw,'(f15.5,3x,a100)') eigenvectors(1,ii),label
                  end if
                  ii = ii + 1
               enddo
            enddo
            enddo
         enddo
      enddo
    enddo
  enddo

      end select
      endif
         end do
      end do
 
 1000 FORMAT (//' Analysis of EOM eigenvectors' &
     &/' First line  :  State Energy, Eigenvalue' &
     &/' other lines : Largest R1 and R2 Vectors',&
     & ' (above a threshold of ',E6.1,')')
 1010 FORMAT (/' Irrep ',A4,' State ',I5,2F16.8)

!for R1 of IP and EA
 1001 format ('| ',a4,' #',i4,' (',f8.3,') |',25x)

!for R2 of EA
 1002 format (1x,a4,' #',i4,' (',f8.3,') -> ', &
     &           a4,' #',i4,' (',f8.3,'), ',   &
     &           a4,' #',i4,' (',f8.3,')  ')


!for R2 of IP
 1003 format (1x,a4,' #',i4,' (',f8.3,'), ',   &
                 a4,' #',i4,' (',f8.3,') -> ', &
     &           a4,' #',i4,' (',f8.3,')  ')

!for R1 of EE
 1004 format (1x,a4,' #',i4,' (',f8.3,') -> ', &
     &             a4,' #',i4,' (',f8.3,')  ')

!for R2 of EE
 1005 format (1x,a4,' #',i4,' (',f8.3,'), ',   &
                 a4,' #',i4,' (',f8.3,') -> ', &
     &           a4,' #',i4,' (',f8.3,'), ',   &
     &           a4,' #',i4,' (',f8.3,')  ')

  end subroutine


subroutine eom_set_davidson_mode_energies(left_vectors, right_vectors)
   use davidson
   use allocator_parameters, only : klongint
   use symmetry_offset
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"
#include "ccpar.inc"
   logical, intent(in) :: left_vectors
   logical, intent(in) :: right_vectors
   logical, save :: davidson_splash = .true.
! 
!  we set the state for the diagonalizer, as this should not change within a given type of calculation,
! e.g. excitation energies
!
    call davidson_setup(unit_output= iw,   &
                        solve_right=right_vectors, &
                        solve_left=left_vectors,   &
                        symmetric=.false., &
                        verbose=relcc_mfd_verbose,   &
                        refresh_trial_rate=relcc_mfd_refresh_rate,  &
                        complex_mode=(rcw.eq.2), &
                        convergence_threshold = relcc_mfd_convergence_threshold, &
                        max_subspace_size = relcc_mfd_max_subspace_size, &
                        max_iterations = relcc_mfd_max_iterations, &
                        overlap_sorting = relcc_mfd_overlap_sorting, &
                        energy_shift=relcc_mfd_eigenvalues_energy_shift, &
                        print_config=davidson_splash)
! aspg: this is to silence davidson setup after the first call. we will know from the
! rest of the output whether we are doing left or right, and all other variables will not
! change for a given run, so we need to only print them once 
    if (davidson_splash) davidson_splash = .false.
    if (left_vectors .and. right_vectors) call quit('EOM does left OR right vectors for now')

end subroutine

subroutine eom_solve_for_left_or_right(state_sym, left_vectors, right_vectors, davResultsL, davResultsR, B, eps)
   use davidson
   use allocator_parameters, only : klongint
   use symmetry_offset

#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "inpt.inc"
#include "ccpar.inc"
!--------------------------------------------------
   integer, intent(in) :: state_sym
   type(davidson_results), intent(inout) :: davResultsL(nrep)
   type(davidson_results), intent(inout) :: davResultsR(nrep)
   type(intermediates), intent (inout)   :: B
   real(kind=8)                          :: eps(*)
   logical, intent(in)                   :: left_vectors
   logical, intent(in)                   :: right_vectors
!--------------------------------------------------
   real*8, pointer :: trial_vector(:,:)
   integer :: irp,iroot,jroot,i,j,irep,jrep,ijrep
   integer :: ndimr1,ndimr2
   integer :: AllocateStatus
   integer :: N, N1, N2 ! size of determinantal space 
   character(len=3) :: left_or_right
   type(Offset) :: e,f,g,h
   logical, pointer :: projector(:) => NULL()

   call eom_set_davidson_mode_energies(left_vectors, right_vectors)

    if (right_vectors) then
        left_or_right = "RHS"
    else if (left_vectors) then
        left_or_right = "LHS"
    end if
!
! be very careful with state_sym and nroot, now these are global variables in this module !!
!
    call create_eom_ndet(2)
    nroot = relcc_eom_nroots(state_sym)

    if (nroot.ne.0) then
       write (iw,*) ""
       write (iw,'(2x,A,2x,A4,A1,I2,A5)') "<<< SOLVING "//left_or_right(1:3)//" EOM-"//eomtype(1:2)// &
       &   " EQUATIONS FOR SYMMETRY",REPNA(irpoff+state_sym),"(",state_sym,") >>>"
       write (iw,*) ""

! setting up information for continuous orbital. this is a trick to perform
! EOM-IP and EOM-EA calculations within the same framework as EOM-EE 
     call set_continuous_orbital(state_sym)

     select case (eomtype)

        case ('EE')

    call alloc_array(e,nrep,f,g,h)
    call auto_symmetry_offset(e,nv,no,.false.,.false.)
    call auto_symmetry_offset_triangular(f,nv,nv)
    call auto_symmetry_offset_triangular(g,no,no)
    call auto_symmetry_offset(h,f%oneNonDirac,g%oneNonDirac,.true.,.true.)

       ndimr1 = e%oneNonDirac(state_sym)
       ndimr2 = h%oneNonDirac(state_sym)
 
    call dealloc_array(e,f,g,h)

           if (nroot.eq.-1) nroot = ndimr1 + ndimr2

           N1 = ndimr1*rcw
           N2 = ndimr2*rcw
           N =  N1 + N2 
           call set_eom_ndet(1, N1)
           call set_eom_ndet(2, N2)

           call set_eom_ndet(1, N1)
           call set_eom_ndet(2, N2)

           allocate (trial_vector(N,nroot), stat = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
           allocate (B%H_II(N),  stat = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

           trial_vector = 0.0d0
           call eom_initialize_HII(B%H_II)
!          B%H_II = 1.0d0

           if (right_vectors.and..not.left_vectors) then
               call input_processor_right(B%Fbar_mi, &
    &                               B%Fbar_ae, &
    &                               B%W_mbej,  &
    &                               trial_vector(1:N1,:), &
    &                               trial_vector(N1+1:N,:),  &
    &                               B%H_II)
           else if (left_vectors.and..not.right_vectors) then    
               call input_processor_left(B%Fbar_mi, &
    &                               B%Fbar_ae, &
    &                               B%W_mbej,davResultsR(state_sym)%eVectorsR(1:N,:),  &
    &                               trial_vector(1:N1,:), &
    &                               trial_vector(N1+1:N,:),  &
    &                               B%H_II)
 
           end if
           if (.not.(relcc_projectors_do_core_valence_separation.or.relcc_projectors_do_restricted_excitation_window)) then
              projector => NULL()
           else
              projector => create_projector(4, (ndimr1 + ndimr2), state_sym, eps)
           end if

        case ('IP')

           ndimr2 = 0
           ndimr1 = 0
     
           call alloc_array(e,nrep)
           call auto_symmetry_offset(e,nv,ncont,.false.,.false.)

           do irep = 1, nrep
              if (iprnt.ge.3) write(*,*)'continuous orbital ',irep,ncont(irep),nroot
              ndimr1 = ndimr1 + ncont(irep)*no(irep)
           enddo

           do irep = 1, nrep
               ndimr2 = ndimr2 + e%oneNonDirac(irep)*noot(irep)
           enddo
          call dealloc_array(e)

           N1 = ndimr1*rcw
           N2 = ndimr2*rcw
           N = N1 + N2 

           call set_eom_ndet(1, N1)
           call set_eom_ndet(2, N2)

           if (nroot.eq.-1) nroot = ndimr1 + ndimr2

           allocate (trial_vector(N,nroot), stat = AllocateStatus)
           IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
           allocate (B%H_II(N)) 
     
           trial_vector = 0.0d0
           call eom_initialize_HII(B%H_II)
 !         B%H_II = 1.0d0

           if (right_vectors.and..not.left_vectors) then
               call input_processor_right(B%Fbar_mi, &
    &                               B%Fbar_ae, &
    &                               B%W_mbej,  &
    &                               trial_vector(1:N1,:), &
    &                               trial_vector(N1+1:N,:),  &
    &                               B%H_II)
           else if (left_vectors.and..not.right_vectors) then    
               call input_processor_left(B%Fbar_mi, &
    &                               B%Fbar_ae, &
    &                               B%W_mbej,davResultsR(state_sym)%eVectorsR(1:N,:),  &
    &                               trial_vector(1:N1,:), &
    &                               trial_vector(N1+1:N,:),  &
    &                               B%H_II)

           end if

           if (.not.(relcc_projectors_do_core_valence_separation.or.relcc_projectors_do_restricted_excitation_window)) then
              projector => NULL()
           else
              projector => create_projector(3, (ndimr1 + ndimr2), state_sym, eps)
           end if

        case('EA')

           ndimr2 = 0
           ndimr1 = 0

           call alloc_array(f,nrep)
           call auto_symmetry_offset(f,no,ncont,.false.,.false.)


            do irep = 1, nrep
             ndimr1 = ndimr1 + ncont(irep)*nv(irep)
            enddo

          do irep = 1, nrep
              ndimr2 = ndimr2 + f%oneNonDirac(irep)*nvvt(irep)
          enddo

          call dealloc_array(f)

           N1 = ndimr1*rcw
           N2 = ndimr2*rcw
           N = N1 + N2

           call set_eom_ndet(1, N1)
           call set_eom_ndet(2, N2)

           if (nroot.eq.-1) nroot = ndimr1 + ndimr2

           allocate (trial_vector(N,nroot))
           allocate (B%H_II(N)) 
     
           trial_vector = 0.0d0

           call eom_initialize_HII(B%H_II)
!          B%H_II = 1.0d0

           if (right_vectors.and..not.left_vectors) then 
               call input_processor_right(B%Fbar_mi, &
    &                               B%Fbar_ae, &
    &                               B%W_mbej,  &
    &                               trial_vector(1:N1,:), &
    &                               trial_vector(N1+1:N,:),  &
    &                               B%H_II)
           else if (left_vectors.and..not.right_vectors) then 
               call input_processor_left(B%Fbar_mi, &
    &                               B%Fbar_ae, &
    &                               B%W_mbej,davResultsR(state_sym)%eVectorsR,  &
    &                               trial_vector(1:N1,:), &
    &                               trial_vector(N1+1:N,:),  &
    &                               B%H_II)
 
           end if
           
           if (.not.(relcc_projectors_do_core_valence_separation.or.relcc_projectors_do_restricted_excitation_window)) then
              projector => NULL()
           else
              projector => create_projector(2, (ndimr1 + ndimr2), state_sym, eps)
           end if

       end select 

       if (right_vectors) then
           call davidson_driver(nroot, B, davResultsR(state_sym), trial_vector,state_sym,projector)
           call eom_analysis(davResultsR(state_sym)%eValues(1:rcw,:),davResultsR(state_sym)%eVectorsR,eps) ! Fock-Space like output for EA, EE and IP

       else if (left_vectors) then

           call davidson_driver(nroot, B, davResultsL(state_sym), trial_vector,state_sym,projector)

       end if

       call free_eom_ndet()

       if (allocated(B%H_II)) deallocate(B%H_II)
       if (associated(projector)) deallocate(projector)
! aspg: trial_vector should have been nullified inside davidson_driver, check
! why this crashes
!        if (associated(trial_vector)) deallocate(trial_vector)

      end if ! this ends the if checking whether there are any roots required for this symmetry

end subroutine

subroutine eom_initialize_HII(H_II)
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"

    real(kind=8), intent(inout) :: H_II(:)
    integer :: i

    H_II = 0.0d0
    do i = 1, size(H_II,1), rcw
        H_II(i) = 1.0d0
    end do

end subroutine




subroutine eom_analysis_allirreps(davResults, eps)
    use davidson
#include "symm.inc"
#include "complex.inc"
#include "files.inc"
#include "ccpar.inc"
#include "results.inc"
!--------------------------------------------------
   type(davidson_results), intent(inout) :: davResults(nrep)
   real(kind=8)                          :: eps(*)
!---------------Local variables -----------------------------------
   integer :: iroot,ioff,isize,dum2
   real*8, allocatable  :: eigenvalues(:,:)
   real*8               :: E_OFF,dum1(rcw)
   integer, allocatable :: nstate_sym(:),ind(:),irps(:),indx(:)
!---------------Executable code--------------------------------------
   allocate(nstate_sym(nrep))
   allocate(ind(nrep))   
   nstate_sym(1:nrep)=relcc_eom_nroots(1:nrep)
   if(eomtype.eq.'EE') then
      isize=sum(relcc_eom_nroots)+1
      allocate(indx(isize))
      allocate(irps(isize))      
      allocate(eigenvalues(rcw,isize))
      nstate_sym(1)=nstate_sym(1)+1
      ioff=2
   else
     isize=sum(relcc_eom_nroots)
     allocate(indx(isize))
     allocate(irps(isize))      
     allocate(eigenvalues(rcw,isize))
     ioff=1
   endif
   eigenvalues=0.0D0
   do state_sym = 1, nrep
   if(relcc_eom_nroots(state_sym).gt.0) then
     eigenvalues(1:RCW,ioff:ioff+relcc_eom_nroots(state_sym)-1)=davResults(state_sym)%eValues(1:RCW,1:relcc_eom_nroots(state_sym))
     ioff=ioff+relcc_eom_nroots(state_sym)
   endif
   enddo
   dum1(1:rcw) = 0.0d0
   dum2 = 0
   E_OFF = ESCF + ECCSD

   CALL PRT_EV (IW,NREP,nstate_sym,REPNA(1+irpoff), &
     &                E_OFF,ECCSDIM,RCW, &
     &                eigenvalues(1,1),isize, &
     &                indx,irps,ind,dum1,dum2)
   deallocate(nstate_sym)
   deallocate(ind)   
   deallocate(indx)
   deallocate(irps)
   deallocate(eigenvalues)
end subroutine

#if 0
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMF (EPS,T1,T2,CT1,CT2,OPTION)
C
      implicit none
C
C---------------Description--------------------------------------------
C
C     Option 1) :
C     Divide parts of T1 and T2 that correspond to WF by denominators
C     Parts that correspond to the effective Hamiltonian are untouched.
C     In the intermediate Hamiltonian IH-2 formalism we zero out Pi->Q
C     excitations (but keep Pm->Q) and corresponding Hamiltonian parts.
C
C     Option 2) :
C     Delete parts of T1 and T2 that correspond to the Hamiltonian.
C
C     Option 3) :
C     In the intermediate Hamiltonian (IH-2) formalism we zero out Pi->Q
C     excitations (but keep Pm->Q) and corresponding Hamiltonian parts.
C     NO division of the T1 and T2 by dets!
C
C---------------Routines called----------------------------------------
C
C     XCOPY
C
C---------------Last modified------------------------------------------
C
C     Authors : Lucas Visscher (LV), Ephraim Eliav (EE) & Andre Zaitsevskii (AZ)
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T1(*),T2(*)
      COMPLEX*16 CT1(*),CT2(*)
      INTEGER OPTION
C
C---------------Common Blocks----------------------------------------
C
#include "symm.inc"
#include "files.inc"
#include "complex.inc"
#include "param.inc"
#include "ihm.inc"
#include "eqns.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3,VALUE,FACIH
      complex*16 CFACIH
      REAL*8 AFAC
      logical yimshift
      INTEGER IOI(0:3), IVA(0:3)
      INTEGER IOJ(0:3), IVB(0:3)
      CHARACTER*2 CLASS(3,2)
      LOGICAL DENO,ZERO,Z1,Z2,Z3,Z4,IHZERO
      LOGICAL DOIHDN
      INTEGER A,AA,ABIJ,acv,acvMIN,AI,AMIN,AOFF,ARP,B,BB,bcv,BRP,i,ico
      INTEGER icomin,ijrp,ii,imin,ioff,irp,isecth,isecth1,isecth2
      INTEGER isectp,isectp1,isectp2,j,jj,jco,jrp
C
C---------------Executable code--------------------------------------
C
      CLASS(1,1) = 'Oi'
      CLASS(2,1) = 'Oa'
      CLASS(3,1) = 'Va'
      CLASS(1,2) = 'Va'
      CLASS(2,2) = 'Vi'
      CLASS(3,2) = 'Oa'

      yimshift= CARITH .and. ( AIH .lt. -1.d-7 )
    
c
c
c The Fock Space blocks construction:
c
c
c      #       O         V
c     -----------------------
c      1      IO        AV
c      2      AO        IV
c      3      AV        AO
C
C T1 part:
C
c
      II = 0
      AI = 0
C
      DO IRP = 1, NREP
c
      IOI(0)=0
      IOI(1)=NIO(IRP)+IOI(0)
      IOI(2)=NAO(IRP)+IOI(1)
      IOI(3)=NAV(IRP)+IOI(2)
c
      do ico=1,3
      DO I = IOI(ico-1)+1, IOI(ico)
         Z1 = .FALSE.
         II = II + 1
         FAC1 = EPS(II)
         ISECTP=0
         IF (DOIH.AND.ICO.EQ.3) THEN
             ISECTP=1
             IF (IPIORB(II).EQ.1) Z1=.TRUE.
         ENDIF
         AOFF = IO(NREP+1) + IV(IRP)
c
         IVA(0)=0
         IVA(1)=NAV(IRP)+IVA(0)
         IVA(2)=NIV(IRP)+IVA(1)
         IVA(3)=NAO(IRP)+IVA(2)
c
         do acv=1,3
C        Find out whether we need to divide by the denominators
C        leave Hamiltonian matrix elements untouched.

C        Find out whether we need to zero non-wf parts
C        This will zero Hamiltonian elements

         IF (OPTION.EQ.2) THEN
            IF ((ICO.EQ.3.AND.ACV.EQ.3).OR.
     &          (ICO.EQ.3.AND.ACV.EQ.1).OR.
     &          (ICO.EQ.2.AND.ACV.EQ.3)) THEN
                ZERO = .TRUE.
            ELSE
                ZERO = .FALSE.
            ENDIF
         ELSE
            ZERO = .FALSE. 
         ENDIF
C
         DO A = IVA(acv-1)+1, IVA(acv)
             Z2 = .FALSE.
             AA = AOFF + A
             FAC = FAC1 - EPS(AA)
             FACIH = FAC
             CFACIH=cmplx(FAC,0.d00)
             ISECTH=0
             IF (DOIH.AND.ACV.EQ.3) THEN
                ISECTH=1
                IF (IPIORB(AA).EQ.1) Z2=.TRUE.
             ENDIF
             AI = AI + 1
             DOIHDN=.FALSE.
             IHZERO=.FALSE.
             DOIHDN = (Z1.OR.Z2).AND.IHSCHEME.EQ.1
             IHZERO = (Z1.OR.Z2).AND.IHSCHEME.EQ.2
             IF (DENO) THEN
                IF (DOIHDN) THEN
c                  Compute modified denominators for extrapolated intermediate Hamiltonian scheme
caz                EEs version of intermediate denominators is assumed 
caz                for positive aih parameter. otherwise, AZs version works 
                   if ( AIH  .gt. 1.d-7) then
c                     positive aih
                      IF (ISECTP.EQ.1.AND.ISECTH.EQ.0) 
     &                   FACIH=FAC+SHIFT_IH(1,2)
                      IF (ISECTP.EQ.0.AND.ISECTH.EQ.1) 
     &                   FACIH=FAC+SHIFT_IH(1,3)
                   else
c                     zero or negative aih
                      if (yimshift) then
c                        imaginary shift                  
                         if (ISECTP.EQ.1 .AND. ISECTH.EQ.0) then
c                           0h1p sector        
                            AFAC = dabs( SHIFT_IH(1,2) ) / 
     &                            ( dsqrt( FAC**2+SHIFT_IH(1,2)**2 ) ) 
                            CFACIH=cmplx(FAC,SHIFT_IH(1,2)*AFAC**NIH)
                        
                         endif    
                         if (ISECTP.EQ.0.AND.ISECTH.EQ.1) then 
c                           1h0p sector        
                             AFAC = dabs( SHIFT_IH(1,3) ) / 
     &                            ( dsqrt( FAC**2+SHIFT_IH(1,3)**2 ) ) 
                             CFACIH=cmplx(FAC,SHIFT_IH(1,3)*AFAC**NIH)
                         endif
                      else
c                        real shift
                         if (ISECTP.EQ.1 .AND. ISECTH.EQ.0) then
c                           0h1p sector        
                            AFAC=1.d0 - ( FAC / (FAC+SHIFT_IH(1,2)) )
                            FACIH=FAC+SHIFT_IH(1,2)*AFAC**NIH
                         endif    
                         if (ISECTP.EQ.0.AND.ISECTH.EQ.1) then 
c                           1h0p sector        
                            AFAC=1.d0 - ( FAC / (FAC+SHIFT_IH(1,3)) )
                            FACIH=FAC+SHIFT_IH(1,3)*AFAC**NIH
                         endif
                      endif ! end imaginary / real shift branching
                   endif ! end positve / negative shift branching
                ENDIF ! end special code for extrapolated IH 

!               Start division of amplitudes by (modified) denominators
                IF (CARITH) THEN
!                  complex arithmetic variant (CARITH True)
                   IF (IHZERO) CT1(AI) = A0
                   if ( yimshift ) then
!                     imaginary shift is used
                      IF (DOIHDN.AND.ABS(CFACIH-cmplx(FAC,0.d+0))
     &                    .GT.1.D-7)  THEN
                          CT1(AI)=CT1(AI)/CFACIH
                      ELSE
                          CT1(AI) = CT1(AI)/FAC
                      ENDIF
                   else                    
!                     real shift is used
                      IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                         CT1(AI)=CT1(AI)/FACIH
caz                      aih is strictly positive because the ee denoms are not active
                         IF (AIH.GT.1.D-7) THEN
                            IF (NIH.GT.0.AND.NIH.LT.100) THEN
                                CT1(AI)=CT1(AI)*
     &                          (1.d0-(AIH*(FACIH-FAC)/FACIH)**NIH)/
     &                          (1.d0-AIH*(FACIH-FAC)/FACIH)
                            ELSE IF (NIH.GE.100) THEN
                                CT1(AI)=CT1(AI)/
     &                          (1.d0 - AIH*(FACIH-FAC)/FACIH)
                            ENDIF
                         ENDIF
                      ELSE
                         CT1(AI) = CT1(AI)/FAC
                      ENDIF
                   endif ! end of real/imaginary shift branching
                ELSE
!                   real arithmetic variant (CARITH False)
                    IF (IHZERO) T1(AI) = A0
                    IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                       T1(AI)=T1(AI)/FACIH
                       IF (AIH.GT.1.D-7) THEN
                          IF (NIH.GT.0.AND.NIH.LT.100) THEN
                             T1(AI)=T1(AI)*
     &                         (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     &                         (1.d0 - AIH*(FACIH-FAC)/FACIH)
                           ELSE IF (NIH.GE.100) THEN
                             T1(AI)=T1(AI)/
     &                         (1.d0 - AIH*(FACIH-FAC)/FACIH)
                           ENDIF
                       ENDIF
                    ELSE
                       T1(AI) = T1(AI)/FAC
                    ENDIF
                ENDIF ! end of real/complex arithmetic branching

             ENDIF ! end of code for division of T1-amplitudes by denominators
C
             IF (ZERO.OR.(OPTION.EQ.3.AND.IHZERO)) THEN
!               ampltitudes will be zeroed (sector is not yet active, or
!               non-dressed, CI-like, IH variant is used)
                IF (CARITH) THEN
                   CT1(AI) = A0
                ELSE
                   T1(AI) = A0
                ENDIF
             ENDIF

             ISECTP=0
             ISECTH=0
         ENDDO ! A (spinor index within specific ACV)
         enddo ! ACV (types of virtual spinors)
      ENDDO ! I (spinor index within specific ICO)
      enddo ! ICO (types of occupied spinors)
      ENDDO ! IRP (symmetry irreps)
C

c T2 part:
      ABIJ = 0
c
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
         JJ = IO(JRP)
         IRP = MULTB(JRP,IJRP+NREP,2)
         IF (IRP.LT.JRP) GOTO 10
c
         IOJ(0)=0
         IOJ(1)=NIO(JRP)+IOJ(0)
         IOJ(2)=NAO(JRP)+IOJ(1)
         IOJ(3)=NAV(JRP)+IOJ(2)
c
         IOFF = IO(IRP) 
c
         IOI(0)=0
         IOI(1)=NIO(IRP)+IOI(0)
         IOI(2)=NAO(IRP)+IOI(1)
         IOI(3)=NAV(IRP)+IOI(2)
c

         do jco=1,3
         DO J = IOJ(jco-1)+1, IOJ(jco)
           Z1 = .FALSE.
           JJ = JJ + 1
           FAC1 = EPS(JJ)
           ISECTP1=0
           IF (DOIH.AND.JCO.EQ.3) THEN
                ISECTP1=1
                IF (IPIORB(JJ).EQ.1) Z1 = .TRUE.
           ENDIF
           icomin=1
           IF (IRP.EQ.JRP) icoMIN = jco
           do ico=icoMIN,3
           IMIN = IOI(ico-1)+1
           IF (IRP.EQ.JRP.and.ico.eq.jco) IMIN = J + 1
           DO I = IMIN, IOI(ico)
             Z2 = .FALSE.
             II = IOFF + I
             FAC2 = EPS(II) + FAC1
             ISECTP2=0
             IF (DOIH.AND.ICO.EQ.3) THEN
                ISECTP2=1
                IF (IPIORB(II).EQ.1) Z2 = .TRUE.
             ENDIF
c
             DO 20 BRP = 1, NREP
                BB = IV(BRP) + IO(NREP+1)
                ARP = MULTB(BRP,IJRP+NREP,2)
                IF (ARP.LT.BRP) GOTO 20
c
                IVB(0)=0
                IVB(1)=NAV(BRP)+IVB(0)
                IVB(2)=NIV(BRP)+IVB(1)
                IVB(3)=NAO(BRP)+IVB(2)
c
                AOFF = IV(ARP) + IO(NREP+1)
c
                IVA(0)=0
                IVA(1)=NAV(ARP)+IVA(0)
                IVA(2)=NIV(ARP)+IVA(1)
                IVA(3)=NAO(ARP)+IVA(2)
c
                do bcv=1,3
                DO B = IVB(bcv-1)+1, IVB(bcv)
                   Z3 = .FALSE.
                   BB = BB + 1
                   FAC3 = FAC2 - EPS(BB)
                   ISECTH1=0
                   IF (DOIH.AND.BCV.EQ.3) THEN
                      ISECTH1=1
                      IF (IPIORB(BB).EQ.1) Z3 = .TRUE.
                   ENDIF
                   acvMIN=1
                   IF (ARP.EQ.BRP) acvMIN = Bcv
                   do acv=acvMIN,3
C                    Find out whether we need to divide by the denoms
C                    leave Hamiltonian matrix elements untouched.
                     IF (OPTION.EQ.1) THEN
                 IF ( (ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.1)
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.1)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.1.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.1)
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.3)
     &            .OR.(             JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.             AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.             AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3             ))
     &                  THEN
                            DENO = .FALSE.
                        ELSE
                            DENO = .TRUE.
                        ENDIF
                     ELSE
                        DENO = .FALSE. 
                     ENDIF
C                    Find out whether we need to zero non-wf parts
C                    This will zero out Hamiltonian elements
                     IF (OPTION.EQ.2) THEN
!lv: Indentation is broken below due to 72 character F77 format limitation
                 IF ( (ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.1)      ! H2 (EF,GH) (0,2) sector
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.3)      ! H2 (MN,OP) (2,0) sector
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.3.AND.BCV.EQ.1)      ! H2 (MF,GP) (1,1) sector
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.2.AND.ACV.EQ.1.AND.BCV.EQ.3)      ! H2 (EN,GP) (1,1) sector
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.1)      ! H2 (MF,OH) (1,1) sector
     &            .OR.(ICO.EQ.2.AND.JCO.EQ.3.AND.ACV.EQ.1.AND.BCV.EQ.3)      ! H2 (EN,OH) (1,1) sector
     &            .OR.(             JCO.EQ.3.AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.             AND.ACV.EQ.3.AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.             AND.BCV.EQ.3)
     &            .OR.(ICO.EQ.3.AND.JCO.EQ.3.AND.ACV.EQ.3             ))
     &                  THEN
!lv: Resume proper indentation
                            ZERO = .TRUE.
                        ELSE
                            ZERO = .FALSE.
                        ENDIF
                     ELSE
                        ZERO = .FALSE. 
                     ENDIF
                     AMIN = IVA(acv-1)+1
                     IF (ARP.EQ.BRP.and.acv.eq.bcv) AMIN = B + 1
                        DO A = AMIN, IVA(acv)
                           Z4 = .FALSE.
                           AA = AOFF + A
                           FAC = FAC3 - EPS(AA)
                           FACIH = FAC
                           CFACIH = cmplx( FAC, 0.d+0)
                           ISECTH2=0
                           IF (DOIH.AND.ACV.EQ.3) THEN
                            ISECTH2=1
                            IF (IPIORB(AA).EQ.1) Z4 = .TRUE.
                           ENDIF
                           ABIJ = ABIJ + 1
                           DOIHDN = .FALSE.
                           IHZERO = .FALSE.
                           DOIHDN = (Z1.OR.Z2.OR.Z3.OR.Z4)
     &                              .AND.IHSCHEME.EQ.1
                           IHZERO = (Z1.OR.Z2.OR.Z3.OR.Z4)
     &                              .AND.IHSCHEME.EQ.2
                           ISECTP=ISECTP1+ISECTP2
                           ISECTH=ISECTH1+ISECTH2
                           IF (DENO) THEN
!lv: Indentation is resetbelow due to 72 character F77 format limitation
! this whole DOIHDN block has to do with the extrapolated IH scheme of AZ & EE
             IF (DOIHDN) THEN
                if (AIH .gt. 1.d-7) then
czai               EEs version .......
                   IF (ISECTP.EQ.1.AND.ISECTH.EQ.0) THEN
                       FACIH=FAC+SHIFT_IH(2,2)
                   ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.1) THEN
                       FACIH=FAC+SHIFT_IH(2,3)
                   ELSEIF (ISECTP.EQ.2.AND.ISECTH.EQ.0) THEN
                      FACIH=FAC+SHIFT_IH(2,5)
                      IF (Z1.AND.Z2) FACIH=FACIH+SHIFT_IH(2,5)
                   ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.2) THEN
                      FACIH=FAC+SHIFT_IH(2,6)
                      IF (Z3.AND.Z4) FACIH=FACIH+SHIFT_IH(2,6)
                   ELSEIF (ISECTP.EQ.1.AND.ISECTH.EQ.1) THEN
                      FACIH=FAC+SHIFT_IH(2,4)
                      IF ((Z1.AND.Z4).OR.(Z1.AND.Z3).OR.
     &                    (Z2.AND.Z4).OR.(Z2.AND.Z3)) 
     &                FACIH=FACIH+SHIFT_IH(2,4)
                   ENDIF
                else
c                  ...... AZs version 
                  if ( yimshift ) then
c                     imaginary shift
                      IF (ISECTP.EQ.1.AND.ISECTH.EQ.0) THEN
                         AFAC=dabs(SHIFT_IH(2,2))/
     &                        dsqrt(FAC**2+SHIFT_IH(2,2)**2)
                         CFACIH=cmplx(FAC,SHIFT_IH(2,2)*AFAC**NIH)
                      ELSE IF (ISECTP.EQ.0.AND.ISECTH.EQ.1) THEN
                         AFAC=dabs(SHIFT_IH(2,3))/
     &                        dsqrt(FAC**2+SHIFT_IH(2,3)**2)
                         CFACIH=cmplx(FAC,SHIFT_IH(2,3)*AFAC**NIH)
                      ELSEIF (ISECTP.EQ.2.AND.ISECTH.EQ.0) THEN
                         AFAC=SHIFT_IH(2,5)
                         IF (Z1.AND.Z2) AFAC=AFAC+SHIFT_IH(2,5) 
                         AFAC=dabs(AFAC)/sqrt(FAC**2+AFAC**2)
                         CFACIH=cmplx(FAC,SHIFT_IH(2,5)*AFAC**NIH) 
                      ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.2) THEN
                         AFAC=SHIFT_IH(2,6)
                         IF (Z3.AND.Z4) AFAC=AFAC+SHIFT_IH(2,6)
                         AFAC=dabs(AFAC)/sqrt(FAC**2+AFAC**2)
                         CFACIH=cmplx(FAC,SHIFT_IH(2,6)*AFAC**NIH)  
                      ELSEIF (ISECTP.EQ.1.AND.ISECTH.EQ.1) THEN
                          AFAC=SHIFT_IH(2,4)
                          IF ((Z1.AND.Z4).OR.(Z1.AND.Z3).OR.
     &                    (Z2.AND.Z4).OR.(Z2.AND.Z3))
     &                    AFAC=AFAC+SHIFT_IH(2,4)
                          AFAC=dabs(AFAC)/sqrt(FAC**2+AFAC**2)
                          CFACIH=cmplx(FAC,SHIFT_IH(2,4)*AFAC**NIH)
                       ENDIF
                  else
c                   real shift
                    IF (ISECTP.EQ.1.AND.ISECTH.EQ.0) THEN
                        AFAC=1.d0 - ( FAC / (FAC+SHIFT_IH(2,2)) )
                        FACIH=FAC+SHIFT_IH(2,2)*AFAC**NIH
                    ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.1) THEN
                        AFAC=1.d0 - ( FAC / (FAC+SHIFT_IH(2,3)) )
                        FACIH=FAC+SHIFT_IH(2,3)*AFAC**NIH
                    ELSEIF (ISECTP.EQ.2.AND.ISECTH.EQ.0) THEN
                         AFAC=SHIFT_IH(2,5)
                         if (Z1.and.Z2) AFAC=AFAC+SHIFT_IH(2,5)
                         FACIH=FACIH+AFAC*(AFAC/(FAC+AFAC))**NIH
                    ELSEIF (ISECTP.EQ.0.AND.ISECTH.EQ.2) THEN
                        AFAC=SHIFT_IH(2,6)
                        IF (Z3.AND.Z4) AFAC=AFAC+SHIFT_IH(2,6)
                        FACIH=FACIH+AFAC*(AFAC/(FAC+AFAC))**NIH
                    ELSEIF (ISECTP.EQ.1.AND.ISECTH.EQ.1) THEN
                        AFAC=SHIFT_IH(2,4)
                        IF ((Z1.AND.Z4).OR.(Z1.AND.Z3).OR.
     &                     (Z2.AND.Z4).OR.(Z2.AND.Z3)) 
     &                  AFAC=AFAC+SHIFT_IH(2,4)
                        FACIH=FACIH+AFAC*(AFAC/(FAC+AFAC))**NIH
                    ENDIF
                  endif ! end of real/imaginary shift branching
                endif ! end positive / negative shift branching
c
             ENDIF ! end of DOIHDN block but we are still inside the T2 DENOM block

!            Check for zero denominators and abort if necessary
             IF( ((.not.yimshift).and.DABS(FACIH).LT.1.D-7).or. 
     &           (DOIHDN.and.yimshift.and.ABS(CFACIH).LT.1.D-7) ) THEN
                        WRITE (IW,1000) IRP,JRP,ARP,BRP,
     &                  I,J,A,B,
     &                  CLASS(ICO,1),CLASS(JCO,1),
     &                  CLASS(ACV,2),CLASS(BCV,2),
     &                  EPS(II),EPS(JJ),EPS(AA),EPS(BB),FAC,FACIH
                        write (IW,*) ' CFACIH' , CFACIH
                        CALL QUIT('Zero denominator in DENOMF')
             ENDIF

             IF (CARITH) THEN
!               complex arithmetics branch
                if( yimshift ) then
!                  complex shift
                   IF (DOIHDN.AND.ABS(CFACIH-FAC).GT.1.D-7) THEN
                     CT2(ABIJ)=CT2(ABIJ)/CFACIH
                   ELSE
                     CT2(ABIJ) = CT2(ABIJ)/FAC
                   ENDIF
                else 
!                  real shift
                   IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                      CT2(ABIJ)=CT2(ABIJ)/FACIH
                      IF ( AIH .GT. 1.D-7) THEN
                         IF (NIH.GT.0.AND.NIH.LT.100) THEN
                              CT2(ABIJ)=CT2(ABIJ)*
     &                        (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     &                        (1.d0 - AIH*(FACIH-FAC)/FACIH)
                         ELSEIF (NIH.GE.100) THEN
                              CT2(ABIJ)=CT2(ABIJ)/
     &                        (1.d0 - AIH*(FACIH-FAC)/FACIH)
                         ENDIF
                       ENDIF
                   ELSE
                     CT2(ABIJ) = CT2(ABIJ)/FAC
                   ENDIF
                endif ! end of real/imaginary shift branching
c
             ELSE
!               real arithmetics branch
                IF (IHZERO) T2(ABIJ) = A0
                IF (DOIHDN.AND.DABS(FACIH-FAC).GT.1.D-7) THEN
                   T2(ABIJ)=T2(ABIJ)/FACIH
                    IF (AIH .GT. 1.D-7) THEN
                       IF (NIH.GT.0.AND.NIH.LT.100) THEN
                          T2(ABIJ)=T2(ABIJ)*
     &                    (1.d0 - (AIH*(FACIH-FAC)/FACIH)**NIH)/
     &                    (1.d0 - AIH*(FACIH-FAC)/FACIH)
                       ELSEIF (NIH.GE.100) THEN
                          T2(ABIJ)=T2(ABIJ)/
     &                    (1.d0 - AIH*(FACIH-FAC)/FACIH)
                       ENDIF
                    ENDIF
                ELSE
                   T2(ABIJ) = T2(ABIJ)/FAC
                ENDIF
             ENDIF ! end of real/complex arithmetic branching
    
                           ENDIF ! End of T2-denominator (DENO True) code block, resume original indentation

                           IF (ZERO.OR.(OPTION.EQ.3.AND.IHZERO)) THEN
                             IF (CARITH) THEN
                                CT2(ABIJ) = A0
                             ELSE
                                T2(ABIJ) = A0
                             ENDIF
                           ENDIF
C-----------------------------------------------------------------------
CThis may cause large output : activate only when debugging toy systems
C                            Check for large amplitudes
C                            IF (CARITH) THEN
C                               VALUE = CDABS(CT2(ABIJ))
C                            ELSE
C                               VALUE = ABS(T2(ABIJ))
C                            ENDIF
C                            IF (VALUE.GT.0.1D0.AND..OPTION.EQ.2) THEN
C                               WRITE (IW,1001) IRP,JRP,ARP,BRP,
C    &                          I,J,A,B,
C    &                          CLASS(ICO,1),CLASS(JCO,1),
C    &                          CLASS(ACV,2),CLASS(BCV,2),
C    &                          EPS(II),EPS(JJ),EPS(AA),EPS(BB),VALUE
C                            ENDIF
C-----------------------------------------------------------------------
                        ENDDO ! A (spinor index within specific ACV)
                        enddo ! ACV (types of virtual spinors)
                  ENDDO ! B (spinor index within specific BCV)
                  enddo ! BCV (types of virtual spinors)
 20          CONTINUE ! BRP (irreps of second virtual spinor)
             ENDDO ! I (spinor index within specific ICO)
             enddo ! ICO (types of occupied spinors)
          ENDDO ! J (spinor index within specific JCO)
         enddo ! JCO ((types of occupied spinors)
 10   CONTINUE ! JRP (irreps of second occupied spinor)
      ENDDO ! IJRP (compound irreps: Gamma(IRP)*Gamma(JRP) = Gamma (ARP) * Gamma(BRP) = Gamma(IJRP)
C
      IF (ABIJ.NE.NDIMT2) CALL QUIT ('ERROR in DENOMF') ! Sanity check: we should have looped over all array elements
C
      RETURN
 1000 FORMAT (//' Zero denominator found in DENOMF',
     &/T30,'I',14X,'J',14x,'A',14X,'B',
     &/' Irreps      : ',4I15,
     &/' Indices     : ',4I15,
     &/' Classes     : ',4(13X,A2),
     &/' Energies    : ',4F15.5,
     &/' Denominators: ',2E15.5,
     &//' Please choose a better model space ! We will stop now.')
 1001 FORMAT (//' Large amplitude in DENOMF',
     &/T30,'I',14X,'J',14x,'A',14X,'B',
     &/' Irreps      : ',4I15,
     &/' Indices     : ',4I15,
     &/' Classes     : ',4(13X,A2),
     &/' Energies    : ',4F15.5,
     &/' Absolute value ',E15.5)
      END
#endif
   end module 
