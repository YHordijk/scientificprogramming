module xc_ac

   use dft_cfg
   use interface_functional_read
   use extra
   use xc_derv
   use xc_max_block_length
   use interface_ao
   use interface_mo

   implicit none

   public dftac
   public dftac_saop_get_dmats
   public dftac_eigenvalues

   private

   real(8) :: homo_eigenvalue
   real(8) :: homo_eigenvalue_previous


contains

    subroutine get_safe_minimum(safe_minimum)
        real(8) :: safe_minimum
        double precision dlamch
        safe_minimum = dlamch('S')
    end subroutine

   subroutine dftac(              &
                    block_length, &
                    derv,         &
                    weight,       &
                    rho,          &
                    rho_resp,     &
                    rho_outer,    &
                    z_0           &
                   )

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(inout) :: derv(max_block_length, 0:99)
      real(8), intent(in)    :: weight(*)
      real(8), intent(in)    :: rho(*)
      real(8), intent(in)    :: rho_resp(*)
      real(8), intent(in)    :: rho_outer(*)
      real(8), intent(in)    :: z_0(*)
!     --------------------------------------------------------------------------
      integer                :: k
      real(8)                :: rho43, f
      real(8)                :: vx_bulk(3), vx_outer(3), vx_slater(3), vx_resp(3)
      real(8)                :: x, shiftbulk
      real(8)                :: safe_minimum
!     --------------------------------------------------------------------------

#ifdef PRG_DIRAC
#include "../cfun/functionals.h"
#else
      print *, 'adapt standalone AC to xcfun'
      stop
#endif
      call get_safe_minimum(safe_minimum)

!     shift bulk?
      if (dft_cfg_grac) then
         shiftbulk = dft_cfg_ac_ip + homo_eigenvalue
      else
         shiftbulk = 0.0d0
      end if

      do k = 1, block_length

!        bulk potential
         vx_bulk = 0.0d0
         if (dft_cfg_saop_with_response_part) then
            vx_slater = 0.0d0
            vx_resp = 0.0d0
#ifdef PRG_DIRAC
            if (rho(k).gt.safe_minimum) &
            call gllbholepot(vx_slater, weight(k), rho(k), dsqrt(z_0(k)))
#endif
            if (rho(k).gt.safe_minimum) &
              vx_resp(1) = 0.42d0*weight(k)*rho_resp(k)/rho(k)

            vx_bulk(1) = vx_slater(1) + vx_resp(1)
         else
            vx_bulk(1) = derv(k, d1000000)
         end if
         vx_bulk(2) = derv(k, d0010000)
        
!        asymptotic potential
         vx_outer = 0.0d0
#ifdef PRG_DIRAC
         if (dft_cfg_asymptote_is_lb94) then
            call lb94pot(vx_outer, weight(k), rho(k), dsqrt(z_0(k)))
         end if
         if (dft_cfg_asymptote_is_lbalpha) then
            call lbalphapot(vx_outer, weight(k), rho(k), dsqrt(z_0(k)))
         end if
#endif

!        interpolation
         if (dft_cfg_saop) then
            if (rho(k).gt.safe_minimum) &
               f = rho_outer(k)/rho(k)
         end if
         if (dft_cfg_grac) then
            rho43 = rho(k)**(4.0d0/3.0d0)
            if (rho43.gt.safe_minimum) then
               x = dsqrt(z_0(k))/rho43
            else
               x = 0.0d0
               return
            end if
            f = 1.0d0/(1.0d0 + exp(-dft_cfg_grac_alpha*(x - dft_cfg_grac_beta)))
         end if
       
!        connection
         derv(k, d0010000) = (1.0d0-f)* vx_bulk(2)
         derv(k, d1000000) = (1.0d0-f)*(vx_bulk(1) - weight(k)*shiftbulk) + f*vx_outer(1)*(1.0d0 - hf_exchange_factor)
      end do

   end subroutine


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* deck dftac_eigenvalues */
      SUBROUTINE DFTAC_EIGENVALUES()
!=======================================================================
!     check orbital energy convergence
!     for asymptotic correction
!     modelled after Andreas Hesselmann's routine HOMO
!     which has been incorporated into this one, so thanks Andreas!
!-----------------------------------------------------------------------
!     radovan bast                         last revision: september 2006
!=======================================================================

      if (.not. interface_ao_set) then
         return
      end if
      if (.not. mo_eigenvalues_available) then
         return
      end if

      if (dft_cfg_grac) then
         if (.not. dft_cfg_grac_is_active) then
            call get_homo_eigenvalue()
            if (dabs(homo_eigenvalue_previous - homo_eigenvalue) <= dft_cfg_ac_threshold) then
               dft_cfg_grac_is_active = .true.
            end if
            homo_eigenvalue_previous = homo_eigenvalue
         end if
      end if

      if (dft_cfg_saop) then
         dft_cfg_saop_is_active = .true.
      end if

      if (dft_cfg_grac_is_active) then
         print *, 'GRAC asymptotic correction active'
      end if
      if (dft_cfg_saop_is_active) then
         print *, 'SAOP asymptotic correction active'
      end if

   end subroutine


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     /* deck dftac_saop_get_dmats */
      SUBROUTINE DFTAC_SAOP_GET_DMATS(DMAT_OUTER_AO,                    &
     &                                DMAT_RESP_AO)
      LOGICAL   FILE_IS_THERE,                                          &
     &          FILE_IS_OPENED
      real(8) ::  DMAT_OUTER_AO(nr_ao*nr_ao*nr_quaternion_blocks),                               &
     &          DMAT_RESP_AO(nr_ao*nr_ao*nr_quaternion_blocks)
      real(8), allocatable :: dmat_outer_mo(:)
      real(8), allocatable :: dmat_resp_mo(:)
      real(8), allocatable :: saop_energies(:)
      integer :: i, j, joff, jhelp, joffnext, np, ni, no
      integer :: k, kcheck, i2basx(2, 2), i2orbx(2, 2), icmoq(2)
      real(8) :: somediff
      real(8) :: toterg, e
      integer :: lusaop, itemp(2)
      integer :: nr_ao_in_fermion_corep(2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      if (.not. mo_eigenvalues_available) then
         print *, 'error: dft ac code cannot find mo_eigenvalues'
         stop
      end if
      if (.not. mo_coef_available) then
         print *, 'error: dft ac code cannot find mo_coef'
         stop
      end if

      allocate(dmat_outer_mo(nr_mo*nr_mo))
      if (dft_cfg_saop_with_response_part) then
         allocate(dmat_resp_mo(nr_mo*nr_mo))
      end if
      allocate(saop_energies(nr_mo))


      call get_homo_eigenvalue()
!
!
!     calculate energy differences
!     ============================
!
      K = 0
      saop_energies = 0.0d0
      DO I = 1,nr_fermion_coreps
        if (i == 1) then
           np = nr_mo_gerade_negative_secondary
           ni = nr_mo_gerade_positive_inactive
        else
           np = nr_mo_ungerade_negative_secondary
           ni = nr_mo_ungerade_positive_inactive
        end if
        IF(NI.GT.0) THEN
          DO J = 1,NI
            K = K+1
            saop_energies(k) =                                          &
     &            homo_eigenvalue                               &
     &           -mo_eigenvalues(J+(I-1)*nr_mo_gerade+NP)
          ENDDO
        ENDIF
      ENDDO
!
!
!     get orbital density matrices
!     ============================
!
      dmat_outer_ao = 0.0d0
      dmat_outer_mo = 0.0d0
      if (dft_cfg_saop_with_response_part) then
        dmat_resp_ao = 0.0d0
        dmat_resp_mo = 0.0d0
      end if
      K = 0
      DO I = 1,nr_fermion_coreps
        if (i == 1) then
           np = nr_mo_gerade_negative_secondary
           ni = nr_mo_gerade_positive_inactive
        else
           np = nr_mo_ungerade_negative_secondary
           ni = nr_mo_ungerade_positive_inactive
        end if
        IF(NI.GT.0) THEN
          DO J = 1,NI
            K = K+1
            JOFF = (nr_mo_gerade*(I-1)+NP)*nr_mo
            JOFF = JOFF+nr_mo_gerade*(I-1)+NP
            JOFF = JOFF+(J-1)*nr_mo+J-1
            IF(dft_cfg_saop_with_response_part) THEN
              DMAT_RESP_MO(1+JOFF)                                      &
     &                                 = SQRT(saop_energies(K))
            ENDIF
            DMAT_OUTER_MO(1+JOFF)                                       &
     &                        = EXP(-2.0d0*(saop_energies(K)**2.0d0))
          ENDDO
        ENDIF
      ENDDO
!
!
!     to avoid numerical problems for degenerate orbitals
!     ===================================================
!
      DO I = 1,nr_fermion_coreps
        if (i == 1) then
           np = nr_mo_gerade_negative_secondary
           ni = nr_mo_gerade_positive_inactive
        else
           np = nr_mo_ungerade_negative_secondary
           ni = nr_mo_ungerade_positive_inactive
        end if
        IF(NI.GT.0) THEN
          DO JHELP = 1,(NI-1)
            J = NI-JHELP
            KCHECK = nr_mo_gerade_positive_inactive*(I-1)+J
!
            JOFF = (nr_mo_gerade*(I-1)+NP)*nr_mo
            JOFF = JOFF+nr_mo_gerade*(I-1)+NP
            JOFF = JOFF+(J-1)*nr_mo+J-1
!
            JOFFNEXT = JOFF+nr_mo+1
!
            SOMEDIFF = ABS(saop_energies(1+KCHECK)                      &
     &                                  -saop_energies(KCHECK))
            IF(SOMEDIFF.LT.1.0D-5) THEN
              IF(dft_cfg_saop_with_response_part) THEN
                dmat_resp_mo(1+joff) = dmat_resp_mo(1+joffnext)
              ENDIF
           dmat_outer_mo(1+joff) = dmat_outer_mo(1+joffnext)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      i2basx = 0
      if (nr_fermion_coreps > 1) then
         i2basx(2, 1) = nr_ao_gerade
         i2basx(1, 2) = nr_ao_gerade*nr_ao_gerade &
                      + nr_ao_gerade*nr_ao_ungerade
         i2basx(2, 2) = i2basx(2, 1) + i2basx(1, 2)
      end if

      i2orbx = 0
      if (nr_fermion_coreps > 1) then
         i2orbx(2, 1) = nr_mo_gerade
         i2orbx(1, 2) = nr_mo_gerade*nr_mo_gerade &
                      + nr_mo_gerade*nr_mo_ungerade
         i2orbx(2, 2) = i2orbx(2, 1) + i2orbx(1, 2)
      end if

      icmoq = 0
      if (nr_fermion_coreps > 1) then
         icmoq(2) = nr_mo_gerade*nr_ao_gerade
      end if

      nr_ao_in_fermion_corep(1) = nr_ao_gerade
      nr_ao_in_fermion_corep(2) = nr_ao_ungerade
!
!
!     transform density matrices
!     ==========================
!
      DO I = 1,nr_fermion_coreps
        if (i == 1) then
           np = nr_mo_gerade_negative_secondary
           ni = nr_mo_gerade_positive_inactive
           no = nr_mo_gerade
        else
           np = nr_mo_ungerade_negative_secondary
           ni = nr_mo_ungerade_positive_inactive
           no = nr_mo_ungerade
        end if
        IF(NI.GT.0) THEN
          IF(dft_cfg_saop_with_response_part) THEN
            CALL QTRANS90('MOAO','S',0.0d0,                             &
     &                  nr_ao_in_fermion_corep(i),nr_ao_in_fermion_corep(i),NO, no, &
     &                  DMAT_RESP_AO(1+I2BASX(I,I)),                    &
     &                  nr_ao, nr_ao,1,pq_to_uq,                            &
     &                  DMAT_RESP_MO(1+I2ORBX(I,I)),                    &
     &                  nr_mo,nr_mo,1,pq_to_uq,                         &
     &                  mo_coef(1+ICMOQ(I)),                           &
     &                  nr_ao_in_fermion_corep(i),no,nr_quaternion_blocks,pq_to_uq,                 &
     &                  mo_coef(1+ICMOQ(I)),                           &
     &                  nr_ao_in_fermion_corep(i),no,nr_quaternion_blocks,pq_to_uq,                 &
     &                  0)
          ENDIF
          CALL QTRANS90('MOAO','S',0.0d0,                               &
     &                  nr_ao_in_fermion_corep(i),nr_ao_in_fermion_corep(i),NO, no, &
     &                DMAT_OUTER_AO(1+I2BASX(I,I)),                     &
     &                  nr_ao, nr_ao,1,pq_to_uq,                            &
     &                DMAT_OUTER_MO(1+I2ORBX(I,I)),                     &
     &                nr_mo,nr_mo,1,pq_to_uq,                           &
     &                mo_coef(1+ICMOQ(I)),                             &
     &                nr_ao_in_fermion_corep(i),no,nr_quaternion_blocks,pq_to_uq,                   &
     &                mo_coef(1+ICMOQ(I)),                             &
     &                nr_ao_in_fermion_corep(i),no,nr_quaternion_blocks,pq_to_uq,                   &
     &                0)
        ENDIF
      ENDDO





      dmat_outer_ao = 2.0d0*dmat_outer_ao
      if (dft_cfg_saop_with_response_part) then
         dmat_resp_ao = 2.0d0*dmat_resp_ao
      end if

      deallocate(dmat_outer_mo)
      if (dft_cfg_saop_with_response_part) then
         deallocate(dmat_resp_mo)
      end if
      deallocate(saop_energies)

      end subroutine

   subroutine get_homo_eigenvalue()

!     --------------------------------------------------------------------------
      integer :: i, j
!     --------------------------------------------------------------------------

      homo_eigenvalue = -1.0d20
      j = 0
      do i = 1, nr_mo_gerade_negative_secondary
         j = j + 1
      end do
      do i = 1, nr_mo_gerade_positive_inactive
         j = j + 1
         if (mo_eigenvalues(j) > homo_eigenvalue) then
            homo_eigenvalue = mo_eigenvalues(j)
         end if
      end do
      do i = 1, nr_mo_gerade_positive_secondary
         j = j + 1
      end do
      do i = 1, nr_mo_ungerade_negative_secondary
         j = j + 1
      end do
      do i = 1, nr_mo_ungerade_positive_inactive
         j = j + 1
         if (mo_eigenvalues(j) > homo_eigenvalue) then
            homo_eigenvalue = mo_eigenvalues(j)
         end if
      end do
      do i = 1, nr_mo_ungerade_positive_secondary
         j = j + 1
      end do

   end subroutine

end module
