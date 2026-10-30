! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Calculates a global lightning flash rate parametrization
!          based on inputs from the convection scheme. The flash rates
!          generated can then be made available to other sections (e.g.
!          UKCA or INFERNO (Fire modelling)).
!
! Method: This routine uses the method of Allen and Pickering (2002,
!         JGR), with modifications to match the UKCA lightning code
!         to ensure the flash rates remain the same (or similar). A
!         full description of the methodology with links to references
!         can be found inline with the calculations below.
!
!----------------------------------------------------------------------
! References:
!
! Allen, D. J. and Pickering, K. E. (2002)
! Evaluation of lightning flash rate parameterizations for use in a
! global chemical transport model.
! Journal of Geophysical Research Vol 107 (D23), 4711,
! https://dx.doi.org/10.1029/2002JD002066
! [Hereafter referred to as AP02]
!
! Price, C. and Rind, D. (1992)
! A Simple Lightning Parametrization for Calculating Global Lightning
! Distributions.
! Journal of Geophysical Research Vol 97 (D9), pp. 9919-9933,
! https://dx.doi.org/10.1029/92JD00719
! [Hereafter referred to as PR92]
!
! Price, C. and Rind, D. (1993)
! What Determines the Cloud-to-Ground Lightning Fraction in
! Thunderstorms?
! Geophysical Research Letters, Vol 20 pp. 463-466.
! https://doi.org/10.1029/93GL00226
! [Hereafter referred to as PR93]
!
! Price, C. and Rind, D. (1994)
! Modeling Global Lightning Distributions in a General Circulation
! Model
! Monthly Weather Review Vol 122, pp 1930-1939.
! https://doi.org/10.1175/1520-0493(1994)122<1930:MGLDIA>2.0.CO;2
! [Hereafter referred to as PR94]
!
! Luhar, A. K., Galbally, I. E., Woodhouse, M. T., and Abraham, N. L.
! Assessing and improving global lightning flash-rate parameterisations
! based on convective cloud-top height, and their impact on NOx
! production and tropospheric composition
! To be submitted to Atmos. Chem. Phys.
!
!-----------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric
!
! Code written in Fortran to conform with UMDP3 programming standards
!-----------------------------------------------------------------------

module fr_conv_price_rind_mod

use um_types,               only: real_umphys
use atm_fields_bounds_mod,  only: tdims

implicit none

character(len=*), parameter, private :: ModuleName = 'FR_CONV_PRICE_RIND_MOD'

contains

subroutine fr_conv_price_rind(lcbase, lctop, land_sea_mask,                    &
                              r_rho_levels, r_theta_levels,                    &
                              latitude, cos_theta_latitude,                    &
                              tot_flash_rate, cg_flash_rate, ic_flash_rate,    &
                              frac_cg, storm_field)

use cderived_mod,           only: delta_lambda, delta_phi
use conversions_mod,        only: recip_pi_over_180
use electric_constants_mod, only: min2sec, zero, one
use electric_inputs_mod,    only: l_pr_update

! Dr Hook modules
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: lcbase(tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end)
! Integer level of the convective cloud base, from the convection scheme.

integer, intent(in) :: lctop( tdims%i_start:tdims%i_end,                       &
                              tdims%j_start:tdims%j_end)
! Integer level of the convective cloud top, from the convection scheme.

logical, intent(in) :: land_sea_mask( tdims%i_start:tdims%i_end,               &
                                      tdims%j_start:tdims%j_end )
! Land/Sea mask (Land points are .true.).

real(kind=real_umphys), intent(in) :: r_rho_levels(                            &
                                          tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end,           &
                                                      1:tdims%k_end)
! Distance from the centre of the planet to the rho levels [m]

real(kind=real_umphys), intent(in) :: r_theta_levels(                          &
                                          tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end,           &
                                                      0:tdims%k_end)
! Distance from the centre of the planet to the theta levels [m]

real(kind=real_umphys), intent(in) :: latitude(tdims%i_start:tdims%i_end,      &
                                               tdims%j_start:tdims%j_end)
! True latitude [radians]

real(kind=real_umphys), intent(in) :: cos_theta_latitude(                      &
                                          tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end)
! Cosine of the latitude, used in the calculation of the grid box area

real(kind=real_umphys), intent(out) ::                                         &
     tot_flash_rate(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end)
! Total lightning flash rate [s-1]

real(kind=real_umphys), intent(out) ::                                         &
     cg_flash_rate(tdims%i_start:tdims%i_end,                                  &
                   tdims%j_start:tdims%j_end)
! Cloud-to-ground flash rate [s-1]

real(kind=real_umphys), intent(out) ::                                         &
     ic_flash_rate(tdims%i_start:tdims%i_end,                                  &
                   tdims%j_start:tdims%j_end)
! Intracloud flash rate [s-1]

real(kind=real_umphys), intent(out) :: frac_cg(tdims%i_start:tdims%i_end,      &
                                               tdims%j_start:tdims%j_end)
! Fraction of lightning which is cloud-to-ground []
! (i.e. cg_flash_rate / tot_flash_rate), with point where no lightning
! occurs taking a value of zero.

logical, intent(out) :: storm_field(tdims%i_start:tdims%i_end,                 &
                                    tdims%j_start:tdims%j_end)
! Logical storm field (points where lightning is calculated are .true.)

! Local variables
real(kind=real_umphys) :: top_alt     ! Cloud top altitude [km]
real(kind=real_umphys) :: base_alt    ! Cloud base altitude [km]
real(kind=real_umphys) :: x_res_deg   ! Resolution of x-grid [degrees long.]
real(kind=real_umphys) :: y_res_deg   ! Resolution of y-grid [degrees lat.]

real(kind=real_umphys) :: gb_norm     ! Grid box normalisation factor [m-1]
real(kind=real_umphys) :: cfact       ! PR grid resolution conversion factor []
                                      ! (See PR94 Figure 1)

real(kind=real_umphys) :: surf_area    ! Grid box surface area [m2]
real(kind=real_umphys) :: norm_fact    ! Flash rate grid box normalisation
                                       ! factor [ ]
real(kind=real_umphys) :: tot_fr1      ! Temporary flash rate [s-1]
real(kind=real_umphys) :: var_cc_thick ! Variation in the cold cloud thickness
                                       ! as a function of latitude [km]
real(kind=real_umphys) :: f_cg_tmp     ! Temporary variable used in the
                                       ! calculation of the cloud-to-ground
                                       ! lightning fraction.
real(kind=real_umphys) :: w_lat        ! Working ABS(latitude) [degrees lat.]

! Variables used in flash rate calculation
! Values chosen via namelist logical l_pr_update
! Flash_rate = land_flash_pref * ctop_height ** land_flash_exp
real(kind=real_umphys) :: land_flash_pref
real(kind=real_umphys) :: land_flash_exp
! Flash_rate = ocean_flash_pref * ctop_height ** ocean_flash_exp
real(kind=real_umphys) :: ocean_flash_pref
real(kind=real_umphys) :: ocean_flash_exp

! PR92 equation 6 (land) constants:
! Flash_rate = land_flash_pref * ctop_height ** land_flash_exp
real(kind=real_umphys), parameter :: land_flash_pref_pr  = 3.44e-05
real(kind=real_umphys), parameter :: land_flash_exp_pr   = 4.90
! Luhar et al. updated constants
real(kind=real_umphys), parameter :: land_flash_pref_upd = 2.40e-05
real(kind=real_umphys), parameter :: land_flash_exp_upd  = 5.09

! PR92 equation 8 (ocean) constants:
! Flash_rate = ocean_flash_pref * ctop_height ** ocean_flash_exp
real(kind=real_umphys), parameter :: ocean_flash_pref_pr  = 6.40e-04
real(kind=real_umphys), parameter :: ocean_flash_exp_pr   = 1.73
! Luhar et al. updated constants
real(kind=real_umphys), parameter :: ocean_flash_pref_upd = 2.00e-05
real(kind=real_umphys), parameter :: ocean_flash_exp_upd  = 4.38

real(kind=real_umphys), parameter :: m2km = 0.001
! Conversion factor from metres to kilometres

real(kind=real_umphys), parameter :: min_conv_depth = 5.0
! Minimum convective cloud depth [km]

real(kind=real_umphys), parameter :: lat2km           = 111.11
! Latitude to km conversion factor

real(kind=real_umphys), parameter :: sqkm2sqm         = 1.0e6
! Square km to square m conversion factor

! Parameters from PR94 Figure 1 to link the grid box size
! to a lightning calibration factor (cfact)
real(kind=real_umphys), parameter :: pr_conv_pref     = 0.97241
real(kind=real_umphys), parameter :: pr_conv_exp      = 0.048203

! Parameters to link the cold cloud thickness variation
! with the cloud-to-ground lightning fraction following
! PR93 Eq (1):
real(kind=real_umphys), parameter :: fcg_a = 0.021
real(kind=real_umphys), parameter :: fcg_b = -0.648
real(kind=real_umphys), parameter :: fcg_c = 7.493
real(kind=real_umphys), parameter :: fcg_d = -36.54
real(kind=real_umphys), parameter :: fcg_e = 63.09

! Parameters to link the cold-cloud thickness with latitude
! PR93 Eq (3):
real(kind=real_umphys), parameter :: cc_a = -6.64e-5
real(kind=real_umphys), parameter :: cc_b = -4.73e-3
real(kind=real_umphys), parameter :: cc_c = 7.34

integer :: i, j ! Loop counters

! Dr Hook variables
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
character(len=*), parameter   :: RoutineName='FR_CONV_PRICE_RIND'

!===============================================================================
! End of declarations and start of the subroutine code proper
!===============================================================================

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------
! Calculate initial constants for scheme
!-----------------------------------------------------------------

! Calculate the resolution of a grid box in degrees
! (This is constant everywhere)
x_res_deg = delta_lambda * recip_pi_over_180
y_res_deg = delta_phi    * recip_pi_over_180

! Calculate the grid box normalisation factor and the Price and Rind
! spatial calibration factor.
gb_norm = one / (sqkm2sqm * ((x_res_deg * lat2km * 0.87_real_umphys)           &
                          * (y_res_deg * lat2km)))

! PR94 Figure 1, repeated as AP02 Eq (8)
cfact = pr_conv_pref * exp(pr_conv_exp * (x_res_deg * y_res_deg))

! Set the values of the constants used in the flash rate calculation
if (l_pr_update) then
  ! use updated values from Luhar et al.
  land_flash_pref  = land_flash_pref_upd
  land_flash_exp   = land_flash_exp_upd
  ocean_flash_pref = ocean_flash_pref_upd
  ocean_flash_exp  = ocean_flash_exp_upd
else
  ! use original PR92 values
  land_flash_pref  = land_flash_pref_pr
  land_flash_exp   = land_flash_exp_pr
  ocean_flash_pref = ocean_flash_pref_pr
  ocean_flash_exp  = ocean_flash_exp_pr
end if

!===============================================================================
! Main Lightning Calculation Loop
!===============================================================================

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP SHARED( cg_flash_rate, ic_flash_rate, tot_flash_rate, frac_cg,           &
!$OMP         storm_field, lcbase, lctop, r_rho_levels, r_theta_levels,        &
!$OMP         delta_lambda, delta_phi, cos_theta_latitude, gb_norm, cfact,     &
!$OMP         land_sea_mask, latitude, tdims, land_flash_pref,                 &
!$OMP         land_flash_exp, ocean_flash_pref, ocean_flash_exp )              &
!$OMP PRIVATE( i, j, base_alt, top_alt, surf_area, tot_fr1, norm_fact,         &
!$OMP          w_lat, var_cc_thick, f_cg_tmp )
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    ! Initialize flash rate arrays to zero
    cg_flash_rate(i,j)  = zero
    ic_flash_rate(i,j)  = zero
    tot_flash_rate(i,j) = zero
    tot_fr1             = zero
    frac_cg(i,j)        = zero
    storm_field(i,j)    = .false.

    if (lcbase(i,j) > 0 .and. lctop(i,j) > 0) then

      !---------------------------------------------------------------------
      ! Calculate convective cloud base and top altitude in km, based on
      ! the cloud base and top model levels provided by the parametrized
      ! (subgrid) convection scheme.
      !---------------------------------------------------------------------
      base_alt = m2km * (r_rho_levels(i,j,lcbase(i,j)) - r_theta_levels(i,j,0))
      top_alt  = m2km * (r_rho_levels(i,j,lctop(i,j))  - r_theta_levels(i,j,0))

      if ( (top_alt - base_alt) > min_conv_depth) then

        ! Set storm field to .true. as calculations for lightning are taking
        ! place.
        storm_field(i,j) = .true.

        !---------------------------------------------------------------------
        ! Calculate the surface area of this grid, along with the
        ! normalisation factor for the lightning flash rates based on
        ! the default grid box area.
        !---------------------------------------------------------------------

        surf_area = delta_lambda * delta_phi * (r_theta_levels(i,j,0)**2) *    &
                    cos_theta_latitude(i,j)

        norm_fact = surf_area * gb_norm * cfact

        !-------------------------------------------------------------------
        ! Cloud is deeper than the minimum convective cloud depth, so
        ! calculate the lightning flash rate appropriately.
        ! Lightning flash rate is related to the cloud top height
        ! in PR92 equation 6 (land) and PR92 equation 8 (ocean).
        ! Note that these equations are repeated in PR94 as equations 1
        ! (land) and 2 (ocean) and again in AP02 as equations 6 (land)
        ! and 8 (ocean).
        !-------------------------------------------------------------------

        if (land_sea_mask(i,j)) then
          ! Land point; use the land equations:
          ! PR92 Eq (6), repeated as AP02 Eq (6)
          tot_fr1 = min2sec * land_flash_pref * ( top_alt ** land_flash_exp )
        else
          ! Ocean/sea point; use the ocean equations:
          ! PR92 Eq (8), repeated as AP02 Eq (7)
          tot_fr1 = min2sec * ocean_flash_pref * ( top_alt ** ocean_flash_exp )
        end if ! land_sea_mask

        !-------------------------------------------------------------------
        ! Normalise the calculation to generate the total flash rate
        !-------------------------------------------------------------------
        tot_flash_rate(i,j) = tot_fr1 * norm_fact

        !-------------------------------------------------------------------
        ! Calculate the cloud-to-ground lightning fraction as a function
        ! of latitude. This is based on the cold-cloud thickness from PR93.
        !-------------------------------------------------------------------

        w_lat = abs(recip_pi_over_180 * latitude(i,j))

        ! Calculate the variation in the cold cloud thickness (in km) as a
        ! function of latitude.
        ! PR93 Eq (3)
        var_cc_thick = ( cc_a * (w_lat**2 )) + ( cc_b * w_lat) + cc_c

        ! PR93 Eq (1), repeated as AP02 Eq (5 (a))
        f_cg_tmp = (fcg_a * (var_cc_thick**4)) + (fcg_b * (var_cc_thick**3))   &
                 + (fcg_c * (var_cc_thick**2)) + (fcg_d *  var_cc_thick)       &
                 + fcg_e

        !-------------------------------------------------------------------
        ! In AP02 equation 5, the value of var_cc_thick is set at a
        ! minimum value of 5.5 km. Note that PR93 state that the
        ! polynomial used to calculate f_cg_tmp above should only be
        ! valid between 5.5 and 14 km. However, for PR93 Eq (3), the cold cloud
        ! thickness can only vary between 6.38 km (poles) and 7.34 km (equator)
        ! so no limits need to be applied here. PR94 explicitly state that for
        ! clouds with cold sectors less than 5.5 km, no CG lightning should
        ! occur.
        !-------------------------------------------------------------------

        ! PR93 Eq (2), repeated as AP02 Eq (5 (b))
        frac_cg(i,j) = one / (f_cg_tmp + one)

        !-------------------------------------------------------------------
        ! Using the cloud-to-ground fraction, we can then simply calculate
        ! the intracloud and cloud-to-ground flash rates
        !-------------------------------------------------------------------
        cg_flash_rate(i,j) = frac_cg(i,j) * tot_flash_rate(i,j)
        ic_flash_rate(i,j) = tot_flash_rate(i,j) - cg_flash_rate(i,j)

      end if ! Convective cloud depth above minimum
    end if ! lcbase and lctop > 0
  end do ! i
end do  ! j
!$OMP END PARALLEL DO

!===============================================================================
! End of the subroutine; no more scientific code below here.
!===============================================================================
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine fr_conv_price_rind
end module fr_conv_price_rind_mod
