!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the convection-based electric (lightning) scheme.

module conv_electric_main_kernel_mod

use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_SCALAR, GH_REAL,                     &
                             GH_READ, GH_WRITE,                                &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             CELL_COLUMN
use fs_continuity_mod, only: W3, WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: conv_electric_main_kernel_type
  private
  type(arg_type) :: meta_args(15) = (/                                         &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                             & ! height_w3
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  WTHETA),                         & ! height_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),      & ! latitude
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),      & ! lowest_cv_base
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),      & ! lowest_cv_top
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),      & ! tile_fraction
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                                  & ! dt
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! tot_flash_rate
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! cg_flash_rate
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! ic_flash_rate
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! frac_cg
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! num_cg_flashes
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! num_ic_flashes
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),      & ! tot_flashes
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)       & ! storm_field
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: conv_electric_main_code
end type

public :: conv_electric_main_code

contains

!> @brief Interface to the UM electric (lightning) scheme driven by convection
!> @param[in]     nlayers          Number of layers
!> @param[in]     height_w3        Height of density space levels above surface
!> @param[in]     height_wth       Height of theta space levels above surface
!> @param[in]     latitude         Latitude (radians)
!> @param[in]     lowest_cv_base   Level number of the base of the lowest
!!                                 convection in the column
!> @param[in]     lowest_cv_top    Level number of the top of the lowest
!!                                 convection in the column
!> @param[in]     tile_fraction    Surface tile fractions
!> @param[in]     dt               The model timestep length (s)
!> @param[in,out] tot_flash_rate   Total flash rate (s-1)
!> @param[in,out] cg_flash_rate    Cloud-to-ground flash rate (s-1)
!> @param[in,out] ic_flash_rate    Intracloud flash rate (s-1)
!> @param[in,out] frac_cg          Fraction of flashes which are
!!                                 cloud-to-ground
!> @param[in,out] num_cg_flashes   Number of cloud-to-ground flashes this
!!                                 timestep
!> @param[in,out] num_ic_flashes   Number of intracloud flashes this timestep
!> @param[in,out] tot_flashes      Total number of flashes this timestep
!> @param[in,out] storm_field      Flag (1/0) to indicate whether a grid column
!!                                 contains a thunderstorm
!> @param[in]     ndf_w3           Number of degrees of freedom per cell for
!!                                 density space
!> @param[in]     undf_w3          Number unique of degrees of freedom in
!!                                 segment for density space
!> @param[in]     map_w3           Dofmap for a segment of any density space
!!                                 field
!> @param[in]     ndf_wth          Number of degrees of freedom per cell for
!!                                 potential temperature space
!> @param[in]     undf_wth         Number unique of degrees of freedom in
!!                                 segment for potential temperature space
!> @param[in]     map_wth          Dofmap for a segment of any potential
!!                                 temperature space field
!> @param[in]     ndf_2d           Number of degrees of freedom per cell for
!!                                 2D fields
!> @param[in]     undf_2d          Number unique of degrees of freedom in
!!                                 segment for 2D fields
!> @param[in]     map_2d           Dofmap for a segment of any 2D field
!> @param[in]     ndf_tile         Number of degrees of freedom per cell for
!!                                 tiled fields
!> @param[in]     undf_tile        Number unique of degrees of freedom in
!!                                 segment for tiled fields
!> @param[in]     map_tile         Dofmap for a segment of any tiled field

subroutine conv_electric_main_code(nlayers, height_w3, height_wth,             &
                                   latitude, lowest_cv_base, lowest_cv_top,    &
                                   tile_fraction, dt,                          &
                                   tot_flash_rate, cg_flash_rate,              &
                                   ic_flash_rate, frac_cg, num_cg_flashes,     &
                                   num_ic_flashes, tot_flashes, storm_field,   &
                                   ndf_w3, undf_w3, map_w3,                    &
                                   ndf_wth, undf_wth, map_wth,                 &
                                   ndf_2d, undf_2d, map_2d,                    &
                                   ndf_tile, undf_tile, map_tile)

    use constants_mod,          only: r_def, r_um, i_def, i_um
    use jules_control_init_mod, only: n_land_tile
    use planet_constants_mod,   only: planet_radius

    !------------------
    ! UM modules
    !------------------
    use fr_conv_price_rind_mod, only: fr_conv_price_rind

    implicit none

    !-----------------------
    ! Subroutine arguments
    !-----------------------
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3, undf_wth, undf_2d, undf_tile
    integer(kind=i_def), intent(in) :: ndf_w3, ndf_wth, ndf_2d, ndf_tile
    integer(kind=i_def), intent(in), dimension(ndf_w3)   :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_2d)   :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_tile) :: map_tile

    real(kind=r_def), intent(in), dimension(undf_w3)   :: height_w3
    real(kind=r_def), intent(in), dimension(undf_wth)  :: height_wth
    real(kind=r_def), intent(in), dimension(undf_2d)   :: latitude
    real(kind=r_def), intent(in), dimension(undf_2d)   :: lowest_cv_base
    real(kind=r_def), intent(in), dimension(undf_2d)   :: lowest_cv_top
    real(kind=r_def), intent(in), dimension(undf_tile) :: tile_fraction
    real(kind=r_def), intent(in) :: dt

    real(kind=r_def), intent(inout), dimension(undf_2d) :: tot_flash_rate
    real(kind=r_def), intent(inout), dimension(undf_2d) :: cg_flash_rate
    real(kind=r_def), intent(inout), dimension(undf_2d) :: ic_flash_rate
    real(kind=r_def), intent(inout), dimension(undf_2d) :: frac_cg
    real(kind=r_def), intent(inout), dimension(undf_2d) :: num_cg_flashes
    real(kind=r_def), intent(inout), dimension(undf_2d) :: num_ic_flashes
    real(kind=r_def), intent(inout), dimension(undf_2d) :: tot_flashes
    real(kind=r_def), intent(inout), dimension(undf_2d) :: storm_field

    !-------------------------------------
    ! Local variables for the kernel
    !-------------------------------------
    integer(i_um), dimension(1,1) :: lcbase_work, lctop_work

    logical, dimension(1,1) :: land_sea_mask, storm_field_work

    real(r_um), dimension(1,1,nlayers)   :: r_rho_levels
    real(r_um), dimension(1,1,0:nlayers) :: r_theta_levels

    real(r_um), dimension(1,1) :: latitude_work, cos_theta_latitude
    real(r_um), dimension(1,1) :: tot_fr_work, cg_fr_work, ic_fr_work
    real(r_um), dimension(1,1) :: frac_cg_work

    real(r_um) :: flandg
    real(r_um) :: num_cg, num_ic

    integer(i_um) :: k, n

    !-------------------------------------------------------------
    ! Convert inputs from the algorithm layer into UM-dimensions
    !-------------------------------------------------------------
    do k = 1, nlayers
      r_rho_levels(1,1,k) = real(height_w3(map_w3(1) + k-1) + planet_radius,   &
                                 r_um)
    end do

    do k = 0, nlayers
      r_theta_levels(1,1,k) = real(height_wth(map_wth(1) + k) + planet_radius, &
                                   r_um)
    end do

    latitude_work(1,1) = real(latitude(map_2d(1)), r_um)
    ! Set to 1 to account for the quasi-uniform grid
    cos_theta_latitude(1,1) = 1.0_r_um
    lcbase_work(1,1) = nint(lowest_cv_base(map_2d(1)), i_um)
    lctop_work(1,1)  = nint(lowest_cv_top(map_2d(1)), i_um)

    flandg = 0.0_r_um
    do n = 1, n_land_tile
      flandg = flandg + real(tile_fraction(map_tile(1) + n-1), r_um)
    end do
    land_sea_mask(1,1) = flandg > 0.0_r_um

    !-------------------------------------------------------------
    ! Calculate the flash rates from the convective cloud
    !-------------------------------------------------------------
    call fr_conv_price_rind( lcbase_work, lctop_work, land_sea_mask,           &
                             r_rho_levels, r_theta_levels,                     &
                             latitude_work, cos_theta_latitude,                &
                             tot_fr_work, cg_fr_work, ic_fr_work,              &
                             frac_cg_work, storm_field_work )

    ! Number of flashes in this timestep
    num_cg = cg_fr_work(1,1) * real(dt, r_um)
    num_ic = ic_fr_work(1,1) * real(dt, r_um)

    !--------------------------------------------------
    ! Pass the results back to the algorithm layer
    !--------------------------------------------------
    tot_flash_rate(map_2d(1)) = real(tot_fr_work(1,1), r_def)
    cg_flash_rate(map_2d(1))  = real(cg_fr_work(1,1), r_def)
    ic_flash_rate(map_2d(1))  = real(ic_fr_work(1,1), r_def)
    frac_cg(map_2d(1))        = real(frac_cg_work(1,1), r_def)
    num_cg_flashes(map_2d(1)) = real(num_cg, r_def)
    num_ic_flashes(map_2d(1)) = real(num_ic, r_def)
    tot_flashes(map_2d(1))    = real(num_cg + num_ic, r_def)

    ! Convert storm field from logical to real type
    if (storm_field_work(1,1)) then
      storm_field(map_2d(1)) = 1.0_r_def
    else
      storm_field(map_2d(1)) = 0.0_r_def
    end if

end subroutine conv_electric_main_code
end module conv_electric_main_kernel_mod
