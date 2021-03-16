module pcm_integrals

    use, intrinsic :: iso_c_binding, only: c_size_t

    implicit none

    public get_mep
    public get_nuclear_mep
    public get_electronic_mep

    private

contains

    subroutine get_mep(nr_points, nr_points_irr, points, mep, dmat, work, lwork, iprint, skipss)

#include "dcbbas.h"
#include "dgroup.h"

        ! Passed variables
        integer(c_size_t), intent(in)  :: nr_points
        integer(c_size_t), intent(in)  :: nr_points_irr
        real(8), intent(in)  :: points(3, nr_points)
        real(8), intent(in)  :: dmat(ntbas(0), ntbas(0), nz)
        real(8), intent(out) :: mep(nr_points)
        real(8)              :: work(*)
        integer              :: iprint, lwork
        logical, optional    :: skipss
        ! Local variables
        logical              :: skipss_block

        ! Skip SS block?
        skipss_block = .false.
        if (present(skipss)) then
            if (skipss) then
                skipss_block = .true.
            end if
        end if

        ! Get nuclear and electronic contribution.
        call get_nuclear_mep(nr_points, nr_points_irr, points, mep, iprint)
        call get_electronic_mep(nr_points, nr_points_irr, dmat, points,       &
            mep, work, lwork, iprint, .false., skipss_block)

    end subroutine

    subroutine get_nuclear_mep(nr_points, nr_points_irr, points, nuc_mep, iprint)
        !
        ! Routine for the calculation of the nuclear part of the molecular
        ! electrostatic potential on a certain grid of points {r_i}:
        !     V_nuc(r_i) = sum_K Z_K/|R_K - r_i|
        !
        !  array nuc_mep(r_i) contains the nuclear MEP at point r_i.
        !
        ! We also take care of symmetry, since it is straightforward
        ! for the nuclei.
        !
        ! RDR 0512.
        !

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"

        !     Passed variables
        integer(c_size_t), intent(in)  :: nr_points
        integer(c_size_t), intent(in)  :: nr_points_irr
        real(8), intent(in)  :: points(3, nr_points)
        real(8), intent(out) :: nuc_mep(nr_points)
        integer :: iprint

        !     Local variables
        real(8) :: coora(3, nucdep), charges(nucdep), dist
        integer :: ipoint, i, j, k
        real(8) :: renorm, temp

        !     Get coordinates and charges of all centers, not only the ones
        !     that are symmetry independent
        call getacord(coora)
        i = 0
        do j=1,nucind
            do k = 1,nucdeg(j)
                i = i + 1
                charges(i) = charge(j)
            enddo
        enddo

        nuc_mep = 0.0d0
        do i = 1, nucdep
            do ipoint = 1, nr_points_irr
                dist = (coora(1, i) - points(1, ipoint))**2 +               &
                    (coora(2, i) - points(2, ipoint))**2 +               &
                    (coora(3, i) - points(3, ipoint))**2
                nuc_mep(ipoint) = nuc_mep(ipoint) + charges(i) / (sqrt(dist))
            end do
        end do

        ! Renormalization
        renorm = dble(maxrep + 1)
        temp = 0.0d0
        do ipoint = 1, nr_points_irr
            temp = nuc_mep(ipoint) * renorm
            nuc_mep(ipoint) = temp
        end do

        if (iprint > 5) then
            do ipoint = 1, nr_points
                print *, "Nuclear ESP @point", ipoint, nuc_mep(ipoint)
            end do
        end if

    end subroutine

    subroutine get_scc_electronic_mep(nr_points, nr_points_irr, points, ele_mep, iprint)
        !
        ! Routine for the calculation of the electronic part of the molecular
        ! electrostatic potential on a certain grid of points {r_i} for the
        ! Small-Small block, using the Simple Coulombic Correction (SCC):
        !
        !     V_nuc(r_i) = sum_K Z_K/|R_K - r_i|
        !
        !  Z_K is now the tabulated small charge for nucleus K.
        !  array nuc_mep(r_i) contains the nuclear MEP at point r_i.
        !
        ! We also take care of symmetry, since it is straightforward
        ! for the nuclei.
        ! The ele_mep vector is not zeroed out in this subroutine!
        !
        ! RDR 0612.
        !

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"

        !     Passed variables
        integer(c_size_t), intent(in)  :: nr_points
        integer(c_size_t), intent(in)  :: nr_points_irr
        real(8), intent(in)  :: points(3, nr_points)
        real(8), intent(out) :: ele_mep(nr_points)
        integer :: iprint

        !     Local variables
        real(8) :: coora(3, nucdep), charges(nucdep), dist
        integer :: ipoint, i, j, k
        real(8) :: renorm, temp
        real(8), allocatable :: scratch(:)

        !     Get coordinates and charges of all centers, not only the ones
        !     that are symmetry independent
        call getacord(coora)
        i = 1
        do j=1, nucind
            !     Get small component charges
            call scdens_2(charges(i), j, j, 0)
            do k = 1, nucdeg(j)-1
              charges(i+k) = charges(i)
            enddo
            i = i + nucdeg(j) 
        enddo

        allocate(scratch(nr_points_irr))
        scratch = 0.0d0
        do i = 1, nucdep
            do ipoint = 1, nr_points_irr
                dist = (coora(1, i) - points(1, ipoint))**2 +               &
                    (coora(2, i) - points(2, ipoint))**2 +               &
                    (coora(3, i) - points(3, ipoint))**2
                scratch(ipoint) = scratch(ipoint) + charges(i) / (sqrt(dist))
            end do
        end do
        ! Renormalize the scratch vector
        renorm = dble(maxrep + 1)
        temp = 0.0d0
        do ipoint = 1, nr_points_irr
            temp = scratch(ipoint) * renorm
            scratch(ipoint) = temp
        end do
        ! Accumulate on top of ele_mep, in the totally symmetric block
        do ipoint = 1, nr_points_irr
            ele_mep(ipoint) = ele_mep(ipoint) + scratch(ipoint)
        end do
        deallocate(scratch)

        if (iprint > 5) then
            do ipoint = 1, nr_points
                print *, "Electronic SS ESP @point", ipoint, ele_mep(ipoint)
            end do
        end if

    end subroutine

    subroutine get_electronic_mep(nr_points, nr_points_irr, matrix, points, vector,    &
            work, lwork, iprint, get_matrix, skipss_block)
        !
        ! Driver routine for the calculation of the electronic part of the molecular
        ! electrostatic potential on a certain grid of points {r_i}:
        !     V_el(r_i) = tr(DV_i)
        ! tr is the trace operator, D is the density matrix, V^i is the matrix of the
        ! "nuclear attraction" integrals calculated at point r_i of the grid:
        !     V_mu,nu,i =  - <mu|1/|r-r_i||nu>
        !
        !  array ADER(mu, nu, r_i) contains these integrals.
        !  array vc(r_i) contains the electronic MEP at point r_i.
        !
        ! RDR 0512.
        !
        ! RDR 010512 CANNOT yet handle 2-component Hamiltonians.
        ! RDR 050312 CANNOT yet handle symmetry.
        ! RDR 220312 This routine will be used to form both potentials and Fock
        !            matrix contribution for PCM.
        !            matrix is Fock or density matrix, vector is potentials
        !            or charges vector.
        !            get_matrix logical is present and TRUE:
        !            charges vector as input, Fock matrix contribution as output.
        !            get_matrix logical is absent or is present and FALSE:
        !            density matrix as input, potentials vector as output.
        !

#include "pi.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#endif
#include "ccom.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
#include "cbisol.h"
#include "shells.h"
#include "onecom.h"
#include "primit.h"
#include "symind.h"
#include "symmet.h"
#include "lmns.h"
#include "nuclei.h"

        ! Passed variables
        integer(c_size_t), intent(in) :: nr_points
        integer(c_size_t), intent(in) :: nr_points_irr
        real(8)             :: matrix(ntbas(0), ntbas(0), nz)
        ! Triggers whether we want PCM Fock matrix contribution or just the
        ! potential.
        logical, optional   :: get_matrix
        real(8)             :: vector(nr_points)
        real(8), intent(in) :: points(3, nr_points)
        real(8)             :: work(*)
        integer             :: iprint, lwork
        logical             :: skipss_block

        ! Local variables
        real(8), allocatable :: ader(:, :, :)
        real(8) :: tolog, tols, factor = 1.0d0
        integer :: ishela, ishelb, ica, icb, ia, ib
        integer :: multa, multb, nhktab, kab
        integer :: ipoint, idxmax
        logical :: do_matrix

        if (nbsym > 1) then
            call quit('ESP calculation cannot handle symmetry!')
        end if
        if (BSS.or.x2c) then
            call quit('ESP calculation not working with X2C Hamiltonian!')
        end if

        ! Fock matrix contribution or potential calculation?
        do_matrix = .false.
        if (present(get_matrix)) then
            if (get_matrix) then
                do_matrix = .true.
            end if
        end if

        tols = thrs**2
        tolog = 2 * log(thrs)
        ! Decide whether we loop on LL block only or not
        if (nosmlv) then
            ! Looping only over LL block only
            idxmax = nlrgsh
        else if  (skipss_block) then
            ! Looping only over LL block only
            idxmax = nlrgsh
        else
            idxmax = kmax
        end if
        ! Loop over bras <mu|
        idena = 0
        do ishela = 1, idxmax
            ica   = lclass(ishela)
            nhkta = nhkt(ishela)
            khkta = khkt(ishela)
            kckta = kckt(ishela)
            call lmnval(nhkta, kckta, lvalua, mvalua, nvalua)
            ncenta = ncent(ishela)
            icenta = nucnum(ncenta, 1)
            mula   = istbao(ishela)
            multa  = mult(mula)
            nuca   = nuco(ishela)
            numcfa = numcf(ishela)
            jsta   = jstrt(ishela)
            corax  = cent(ishela, 1, 1)
            coray  = cent(ishela, 2, 1)
            coraz  = cent(ishela, 3, 1)
            ! Loop over kets |nu>
            idenb = 0
            do ishelb = 1, ishela
                icb   = lclass(ishelb)
                ldiag = ishela .eq. ishelb
                nhktb = nhkt(ishelb)
                khktb = khkt(ishelb)
                kcktb = kckt(ishelb)
                call lmnval(nhktb, kcktb, lvalub, mvalub, nvalub)
                ncentb = ncent(ishelb)
                nhktab = nhkta + nhktb
                mulb   = istbao(ishelb)
                multb  = mult(mulb)
                nucb   = nuco(ishelb)
                numcfb = numcf(ishelb)
                jstb   = jstrt(ishelb)
                corbx  = cent(ishelb, 1, 1)
                corby  = cent(ishelb, 2, 1)
                corbz  = cent(ishelb, 3, 1)
                khktab = khkta * khktb
                kcktab = kckta * kcktb
                mab    = ior(mula, mulb)
                kab    = iand(mula ,mulb)
                hkab   = fmult(kab)
                ! Calculate -<mu|1/|r-r_i||nu> integrals for each shell pair mu, nu and point r_i,
                ! then contract with the density matrix. The minus sign accounts for the
                ! charge of the electron.
                if (ica == icb) then
                    allocate(ader(nr_points, kckta, kcktb))
                    ader = 0.0d0
                    call vc_shell(nr_points, ader, tolog, tols,               &
                        points, iprint, work, lwork)
                    if (ishela /= ishelb) then
                        ader = 2.0d0 * ader
                    end if
                    ! We multiply by +2.0d0. The D matrix refers to the alpha part only.
                    !     do ipoint = 1, nr_points
                    !        print *, "ELECTROSTATIC_POTENTIAL_MATRIX @point", ipoint
                    !        do ib = 1, kcktb
                    !          do ia = 1, kckta
                    !             print *, "ELEMENT", ia, ib
                    !             print *, ader(ipoint, ia, ib)
                    !        call output(ader(ipoint,:,:),1,kckta,1,kcktb,kckta,kcktb,2,6)
                    !          enddo
                    !        enddo
                    !     enddo
                    do ipoint = 1, nr_points
                        if (do_matrix) then
                            do ib = 1, kcktb
                                do ia = 1, kckta
                                    matrix(idena+ia,idenb+ib, 1) = matrix(idena+ia,idenb+ib, 1) &
                                        - ader(ipoint, ia, ib) * vector(ipoint)
                                end do
                            end do
                        else
                            vector(ipoint) = vector(ipoint) + 2.0d0                    * &
                                sum(matrix(idena+1:idena+kckta, idenb+1:idenb+kcktb, 1) * &
                                ader(ipoint, 1:kckta, 1:kcktb))
                        end if
                    end do
                    deallocate(ader)
                end if
                idenb = idenb + khktb * multb
            end do
            idena = idena + khkta * multa
        end do

        if (skipss_block .and. (.not. do_matrix) .and. (.not.nosmlv)) then
            ! Simple Coulombic Correction for the MEP when skipss is selected.
            ! This must NOT be done when:
            ! 1. we want a Fock matrix contribution;
            ! 2. NOSMLV, i.e. the Hamiltonian specified does not make use of the SS-block
            call get_scc_electronic_mep(nr_points, nr_points_irr, points, vector, iprint)
        end if

        !          do ipoint = 1, nr_points
        !            print *, "Electronic ESP @point", ipoint, vector(ipoint)
        !          end do
        if (iprint > 5) then
            do ipoint = 1, nr_points
                print *, "Electronic ESP @point", ipoint, vector(ipoint)
            end do
        end if

    end subroutine


    subroutine vc_shell(nr_points, ader, tolog, tols,                 &
            points, iprint, work, lwork)
        !
        ! Calculates the contribution for one primitive orbital set.
        !
        ! RDR 090312 Clean-up.
        !
#include "pi.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#endif
#include "cbisol.h"
#include "onecom.h"
#include "ader.h"
#include "primit.h"

        ! Passed variables
        integer(c_size_t), intent(in) :: nr_points
        integer, intent(in) :: iprint, lwork
        real(8) :: points(3, nr_points), ader(nr_points, kckta, kcktb)
        real(8) :: work(*)
        real(8) :: tolog, tols, factor=1.0

        ! Local variables
        real(8), allocatable :: ahgtf(:, :, :, :), odc(:, :, :, :, :, :)
        real(8), allocatable :: r(:, :, :, :, :)
        real(8) :: difab(3), corp(3), origin(3), cora(3), corb(3)
        real(8) :: difcp(3, nr_points)
        real(8) :: distab, conta, expa, contb, expb, expp, exppi,         &
            &           expabq, saab, asaab, saab13, expapi, expbpi
        integer :: jmaxa, jmaxb, jmaxd, jmaxt, jmaxm,               &
            &           ipoint, iprima, iprimb, jprima, jprimb
        integer :: idummy
        real(8) :: pval

        ! Allocation
        jmaxd = 2
        jmaxa = nhkta - 1
        jmaxb = nhktb - 1
        jmaxt = jmaxa + jmaxb + jmaxd
        jmaxm = 0
        jmax = jmaxa + jmaxb
        ! Initialization
        allocate(odc(0:jmaxa,0:jmaxb,0:jmaxt,0:jmaxd,0:jmaxm,3))
        allocate(ahgtf(nr_points, 0:jmax, 0:jmax, 0:jmax))
        allocate(r(nr_points, 0:jmax, 0:jmax, 0:jmax, 0:jmax))

        cora(1) = corax
        cora(2) = coray
        cora(3) = coraz
        corb(1) = corbx
        corb(2) = corby
        corb(3) = corbz
        difab(:) = cora(:) - corb(:)
        distab = difab(1) * difab(1) + difab(2) * difab(2) + difab(3) * difab(3)

        ! Loop over primitive orbitals
        ! Shell a
        do iprima = 1, nuca
            jprima = jsta + iprima
            conta = priccf(jprima, numcfa)
            expa = priexp(jprima)
            ! Shell b
            do iprimb = 1, nucb
                jprimb = jstb + iprimb
                contb = priccf(jprimb, numcfb)
                expb = priexp(jprimb)
                expp = expa + expb
                exppi = 1.0d0 / expp
                expabq = expa * expb * exppi
                saab = conta * contb * exp(-expabq * distab)
                asaab = abs(saab)
                if (expabq * distab < tolog) then
                    cycle
                end if
                if (asaab < tols) then
                    cycle
                end if
                saab13 = sign(asaab**(1.0d0/3.0d0), saab)
                expapi  = expa * exppi
                expbpi  = expb * exppi
                corp(:) = expapi * cora(:) + expbpi * corb(:)
                do ipoint = 1, nr_points
                    difcp(:, ipoint) =  points(:, ipoint) - corp(:)
                end do
                ! Calculate the Overlap Distribution Coefficients
                idummy = 0
                call getodc(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm, .false.,  &
                    &           .false., onecen, expa, expb, iprint, saab13, exppi,    &
                    &           work, lwork, corp(1), corp(2), corp(3),                &
                    &           .true.,.false.,origin,idummy)
                call vnuc_vec(ahgtf, r, nr_points, jmax, expp, difcp)
                call cart_vc_vec(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm,        &
                    &            ader, ahgtf, nr_points)
            end do ! Close loop over second shell
        end do ! Close loop over first shell
        deallocate(odc)
        deallocate(ahgtf)
        deallocate(r)

    end subroutine


    subroutine vnuc_vec(ahgtf, r, nr_points, jmax, pval, cp)
        !
        ! This subroutine calculates the R integrals as defined by
        ! McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
        ! The recursion formulas (4.6) - (4.8) are used.
        !
        ! JHS 260308 Only slightly slower than the implementation hernai
        !            in abacus/her1car.f of TUH
        ! RDR 090312 Clean-up.
        !
#include "maxaqn.h"
#include "gamcom.h"

        ! Parameters
        real(8), parameter :: pi = acos(-1.0d0)

        ! Passed variables
        integer(c_size_t) :: nr_points
        integer :: jmax
        real(8), intent(out) :: ahgtf(nr_points, 0:jmax, 0:jmax, 0:jmax)
        real(8), intent(in)  :: cp(3, nr_points)
        real(8) :: r(nr_points, 0:jmax, 0:jmax, 0:jmax, 0:jmax)
        real(8) :: pval

        ! Local variables
        real(8) :: factor, prod
        integer :: jval, t, u, v, ipoint

        do ipoint = 1, nr_points
            ! Incomplete gamma function
            wval = pval * (cp(1,ipoint)**2 + cp(2,ipoint)**2 + cp(3,ipoint)**2)
            jmax0 = jmax
            call gamma_function
            ! Calculate r(ipoint, 0, 0, 0, jval)
            factor = (2.0d0 * pi) / pval
            do jval = 0, jmax
                fjw(jval)         =   factor * fjw(jval)
                factor            = - 2.0d0 * pval * factor
                r(ipoint, 0, 0, 0, jval)  =   fjw(jval)
            end do

            ! Calculate r(t, u, v, jval)
            do jval = jmax, 1, -1
                do v = 0, jmax - jval
                    prod = -cp(3, ipoint) * r(ipoint, 0, 0, v, jval)
                    if (v > 0) then
                        prod = prod + v * r(ipoint, 0, 0, v - 1, jval)
                    end if
                    r(ipoint, 0, 0, v + 1, jval - 1) = prod
                    do u = 0, jmax - jval - v
                        prod = -cp(2, ipoint) * r(ipoint, 0, u, v, jval)
                        if (u > 0) then
                            prod = prod + u * r(ipoint, 0, u - 1, v, jval)
                        end if
                        r(ipoint, 0, u + 1, v, jval - 1) = prod
                        do t = 0, jmax - jval - u - v
                            prod = -cp(1, ipoint) * r(ipoint, t, u, v, jval)
                            if (t > 0) then
                                prod = prod + t * r(ipoint, t - 1, u, v, jval)
                            end if
                            r(ipoint, t + 1, u, v, jval - 1) = prod
                        end do
                    end do
                end do
            end do
        end do
        ! The nuclear attraction integrals are given as r(ipoint, t, u, v, 0)
        do v = 0, jmax
            do u = 0, jmax
                do t = 0, jmax
                    do ipoint = 1, nr_points
                        ahgtf(ipoint, t, u, v) = - r(ipoint, t, u, v, 0)
                    end do
                end do
            end do
        end do

    end subroutine

    subroutine cart_vc_vec(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm,        &
            &                 ader, ahgtf, nr_points)
        !
        ! Originally written by JHS.
        !
        ! RDR 060312 Clean-up.
        !

#include "maxaqn.h"
#include "onecom.h"
#include "lmns.h"

        ! Passed variables
        integer :: jmaxa, jmaxb, jmaxt, jmaxd, jmaxm
        integer(c_size_t) :: nr_points
        real(8) :: ahgtf(nr_points, 0:jmax, 0:jmax, 0:jmax)
        real(8) :: ader(nr_points, kckta, kcktb),                                 &
            &           odc(0:jmaxa, 0:jmaxb, 0:jmaxt, 0:jmaxd, 0:jmaxm, 3)

        ! Local variables
        real(8) :: ev, ee, eee
        integer :: icompa, lvala, mvala, nvala,                           &
            &           icompb, lvalb, mvalb, nvalb,                           &
            &           t, u, v, ipoint

        do icompa = 1,kckta
            lvala = lvalua(icompa)
            mvala = mvalua(icompa)
            nvala = nvalua(icompa)
            do icompb = 1,kcktb
                lvalb = lvalub(icompb)
                mvalb = mvalub(icompb)
                nvalb = nvalub(icompb)
                do v = 0, nvala + nvalb
                    ev = odc(nvala,nvalb,v,0,0,3)
                    do u = 0, mvala + mvalb
                        ee = odc(mvala,mvalb,u,0,0,2)*ev
                        do t = 0, lvala + lvalb
                            eee = odc(lvala,lvalb,t,0,0,1)*ee
                            do ipoint = 1, nr_points
                                ader(ipoint, icompa, icompb) = ader(ipoint, icompa, icompb) + eee * ahgtf(ipoint, t, u, v)
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine

    subroutine gamma_function
        !
        ! Trygve Ulf Helgaker fall 1984
        !
        ! This subroutine calculates the incomplete gamma function as
        ! described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
        !
        ! Roberto Di Remigio May 2012
        ! Purified from the evil implicit.h and all the other common blocks.
        !
#include "maxaqn.h"

        ! A long and boring parameter list
        real(8), parameter :: d1 = 1.0d0,  d2 = 2.0d0, d10 = 10.0d0
        real(8), parameter :: half = 0.5d0, tenth = 0.1d0, ten6 = 1.0d6
        real(8), parameter :: pi = acos(-1.0d0)
        real(8), parameter :: sqrtpi = sqrt(pi)
        real(8), parameter :: pi2 = pi * pi
        real(8), parameter :: sqrpih = sqrtpi/d2
        real(8), parameter :: coef2 = half,  coef3 = - d1/6.0d0, coef4 = d1/24.0d0
        real(8), parameter :: coef5 = - d1/120.0d0, coef6 = d1/720.0d0
        real(8), parameter :: gfac30 = 0.4999489092d0, gfac31 = -0.2473631686d0,       &
            &   gfac32 = 0.321180909d0, gfac33 = -0.3811559346d0, gfac20 = 0.4998436875d0,  &
            &   gfac21 = -0.24249438d0, gfac22 = 0.24642845d0, gfac10 = 0.499093162d0,      &
            &   gfac11 = -0.2152832d0, gfac00 = 0.490d0

        ! Passed variables

        ! Local variables
        real(8)            :: tabjfw, wdif, d2wal, rexpw, denom, rwval, summ, term
        real(8)            :: r2max1, d2max1, gval, factor
        integer            :: jmax, jmx, j, iadr, jadr, maxj0, istart, ipoint, iorder
        !
#include "gamcom.h"
        !
        save maxj0
        data maxj0 /-1/
        !
        ipoint = d10 * min(wval, ten6) + half
        !     have seen problems with NINT intrinsic function here (rarely)
        !     therefore the "+ HALF" before integer truncation
        if (ipoint < 0) then
            call quit('Fatal error in gammafun')
        else if (ipoint < 120) then
            istart = 1 + 121 * jmax0 + ipoint
            wdif = wval - tenth * ipoint
            fjw(jmax0) = (((((coef6 * tabfjw(istart + 726) * wdif    &  ! 726 = 6*121
                &                   + coef5 * tabfjw(istart + 605)) * wdif   &
                &                    + coef4 * tabfjw(istart + 484)) * wdif  &
                &                     + coef3 * tabfjw(istart + 363)) * wdif   &
                &                      + coef2 * tabfjw(istart + 242)) * wdif  &
                &                       - tabfjw(istart + 121)) * wdif       &
                &                        + tabfjw(istart)
            d2wal = d2 * wval
            rexpw = exp(-wval)
            denom = 2.0d0 * jmax0 + 1.0d0
            do j = jmax0, 1, -1
                denom = denom - d2
                fjw(j - 1) = (d2wal * fjw(j) + rexpw) / denom
            end do
        else if (ipoint <= (20 * jmax0 + 360)) then
            rwval = d1 / wval
            rexpw = exp(-wval)
            gval = gfac30 + rwval * (gfac31 + rwval * (gfac32 + rwval * gfac33))
            fjw(0) = sqrpih * sqrt(rwval) - rexpw * gval * rwval
            factor = half * rwval
            term = factor * rexpw
            do j = 1, jmax0
                fjw(j) = factor * fjw(j - 1) - term
                factor = rwval + factor
            end do
        else
            rwval  = d1 / wval
            fjw(0) = sqrpih * sqrt(rwval)
            factor = half * rwval
            do j = 1, jmax0
                fjw(j) = factor * fjw(j - 1)
                factor = rwval + factor
            end do
        end if
        return
        !
        !     ***** Tabulation of incomplete gamma function *****
        !
        entry gamtab(jmx)
        !
        !     For j = jmx a power series expansion is used, see for
        !     example Eq.(39) given by V. Saunders in "Computational
        !     Techniques in Quantum Chemistry and Molecular Physics",
        !     Reidel 1975.  For j < jmx the values are calculated
        !     using downward recursion in j.
        !
        !
        if (jmx > maxj) then
            call quit('Gamtab error: jmx greater than limit.')
        end if
        jmax = jmx + 6
        maxj0 = jmax
        !
        !     WVAL = 0.0
        !
        iadr = 1
        denom = d1
        do j = 0, jmax
            tabfjw(iadr) = d1 / denom
            iadr = iadr + 121
            denom = denom + d2
        end do
        !
        !     WVAL = 0.1, 0.2, 0.3,... 12.0
        !
        iadr = iadr - 121
        d2max1 = 2.0d0 * jmax + 1.0d0
        r2max1 = d1 / d2max1
        do ipoint = 1, 120
            wval = tenth * ipoint
            d2wal = wval + wval
            iadr = iadr + 1
            term = r2max1
            summ = term
            denom = d2max1
            do iorder = 2, 200
                denom = denom + d2
                term = term * d2wal / denom
                summ = summ + term
                if (term .le. 1.0d-15) exit
            end do
            rexpw = exp(-wval)
            tabfjw(iadr) = rexpw * summ
            denom = d2max1
            jadr = iadr
            do j = 1, jmax
                denom = denom - d2
                tabfjw(jadr - 121) = (tabfjw(jadr) * d2wal + rexpw) / denom
                jadr = jadr - 121
            end do
        end do

    end subroutine

end module pcm_integrals
