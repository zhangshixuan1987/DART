! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use        types_mod, only : r8, i8, MISSING_R8, digits12

use time_manager_mod, only : time_type, set_time, GREGORIAN

use     location_mod, only : location_type, get_close_type, &
                             set_location, set_location_missing, &
                             set_vertical_localization_coord, &
                             VERTISHEIGHT, VERTISLEVEL, &
                             loc_get_close => get_close

use    utilities_mod, only : register_module, error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             to_upper

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 NF90_MAX_NAME

use state_structure_mod, only : add_domain, get_domain_size, get_model_variable_indices

use obs_kind_mod, only : get_index_for_quantity

use ensemble_manager_mod, only : ensemble_type

use default_model_mod, only : read_model_time, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time

! module variables
character(len=256), parameter :: source   = "wrf/model_mod.f90"
logical :: module_initialized = .false.
type(time_type) :: assimilation_time_step

integer, parameter :: MAX_STATE_VARIABLES = 100
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 4
integer, parameter :: NUM_BOUNDS_TABLE_COLUMNS = 4

integer, allocatable :: dom_id(:)


!-- Namelist with default values --
logical :: default_state_variables = .true.
character(len=NF90_MAX_NAME) :: wrf_state_variables(MAX_STATE_VARIABLES*NUM_STATE_TABLE_COLUMNS) = 'NULL'
character(len=NF90_MAX_NAME) :: wrf_state_bounds(num_bounds_table_columns,max_state_variables) = 'NULL'
integer :: num_domains = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
! Max height a surface obs can be away from the actual model surface
! and still be accepted (in meters)
real (kind=r8) :: sfc_elev_max_diff  = -1.0_r8   ! could be something like 200.0_r8

real (kind=r8) :: center_search_half_length = 500000.0_r8
real(r8) :: circulation_pres_level = 80000.0_r8
real(r8) :: circulation_radius     = 108000.0_r8
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord = VERTISHEIGHT

! Allow observations above the surface but below the lowest sigma level.
logical :: allow_obs_below_vol = .false.

! Do the interpolation of pressure values only after taking the log (.true.)
! vs doing a linear interpolation directly in pressure units (.false.)
logical :: lallow_obs_below_vol = .true.
logical :: log_horz_interpM = .false.
logical :: log_horz_interpQ = .false.


logical :: allow_perturbed_ics = .false.
!-------------------------------


namelist /model_nml/ &
default_state_variables, &
wrf_state_variables, &
wrf_state_bounds, &
num_domains, &
calendar_type, &
assimilation_period_seconds, &
sfc_elev_max_diff, &
center_search_half_length, &
circulation_pres_level, &
circulation_radius, &
center_spline_grid_scale, &
vert_localization_coord, &
allow_perturbed_ics, &
allow_obs_below_vol, &
allow_obs_below_vol, &
log_horz_interpM, &
log_horz_interpQ

type grid
  real(r8), dimension(:,:),   allocatable :: latitude, latitude_u, latitude_v
  real(r8), dimension(:,:),   allocatable :: longitude, longitude_u, longitude_v
end type grid

integer(i8) :: model_size



contains

!------------------------------------------------------------------

subroutine static_init_model()

integer  :: iunit, io

character(len=NF90_MAX_NAME) :: varname(MAX_STATE_VARIABLES)
integer :: state_qty(MAX_STATE_VARIABLES)
logical :: update_var(MAX_STATE_VARIABLES)
character(len=9) :: in_domain(MAX_STATE_VARIABLES) ! assumes <=9 or 999

integer :: nfields
logical, allocatable :: domain_mask(:)
integer :: i, field ! loop indices
integer :: model_dt, assim_dt
character (len=1)     :: idom ! assumes <=9

module_initialized = .true.

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.
model_dt = 1
print*, 'FAKE model_dt', model_dt
assim_dt = (assimilation_period_seconds / model_dt) * model_dt
assimilation_time_step = set_time(assim_dt)

allocate(dom_id(num_domains))

call verify_state_variables(nfields, varname, state_qty, update_var, in_domain)

model_size = 0
allocate(domain_mask(nfields))

do i = 1, num_domains

  do field = 1, nfields
     domain_mask(field) = variable_is_on_domain(in_domain(field), i)
  end do

  write( idom , '(I1)') i
  dom_id(i) = add_domain('wrfinput_d0'//idom, &
                          num_vars=count(domain_mask), &
                          var_names = pack(varname(1:nfields), domain_mask), &
                          kind_list = pack(state_qty(1:nfields), domain_mask), &
                          !clamp_vals  = &
                          update_list = pack(update_var(1:nfields), domain_mask) )
   
  model_size = model_size + get_domain_size(dom_id(i))
enddo

call set_vertical_localization_coord(vert_localization_coord)
 
deallocate(domain_mask)

end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

if ( .not. module_initialized ) call static_init_model

expected_obs(:) = MISSING_R8

istatus(:) = 1

end subroutine model_interpolate



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty_out)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: qty_out

integer :: i, j, k, id, var_id, state_id, qty

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, i, j, k, var_id=var_id, dom_id=state_id, kind_index=qty)

location = convert_indices_to_lon_lat_lev(i, j, k, qty)

! return DART variable qty if requested
if(present(qty_out)) qty_out = qty

end subroutine get_state_meta_data

!------------------------------------------------------------------
function convert_indices_to_lon_lat_lev(i, j, k, qty)

integer, intent(in) :: i, j, k
integer, intent(in) :: qty
type(location_type) :: convert_indices_to_lon_lat_lev

print*, 'not done'
convert_indices_to_lon_lat_lev = set_location(1.0_r8,2.0_r8,3.0_r8, VERTISLEVEL)

end function convert_indices_to_lon_lat_lev


!------------------------------------------------------------------
! obs have a type and qty
!  observation type not taken in to account for wrf get close calulations
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:), loc_types(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, state_handle)

end subroutine get_close_obs

!------------------------------------------------------------------
! state only has qty
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer(i8),                   intent(in)     :: loc_indx(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

call get_close(gc, base_loc, base_type, locs, loc_qtys, &
               num_close, close_ind, dist, state_handle)

end subroutine get_close_state

!------------------------------------------------------------------
subroutine get_close(gc, base_loc, base_type, locs, loc_qtys, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs/state locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs/state
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close'
logical :: fail

call convert_base_vertical(base_loc, fail)
if (fail) then
   num_close = 0
   return
endif

call loc_get_close(gc, base_loc, base_type, locs, loc_qtys, &
                          num_close, close_ind)

if (.not. present(dist)) return

!call convert_vertical()

!call calculate_distances()

end subroutine get_close

!------------------------------------------------------------------
subroutine convert_base_vertical(base_loc, fail)

type(location_type), intent(inout) :: base_loc
logical, intent(out) :: fail

fail = .true.

end subroutine convert_base_vertical
!------------------------------------------------------------------

subroutine end_model()


end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to netcdf files
subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!------------------------------------------------------------------
subroutine pert_model_copies(ens_handle, ens_size,  dummy_pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: dummy_pert_amp ! not used
logical,             intent(out)   :: interf_provided

if (.not. allow_perturbed_ics) then
call error_handler(E_ERR,'pert_model_copies', &
                     'starting WRF model from a single vector requires additional steps', &
                      text2='see comments in wrf/model_mod.f90::pert_model_copies()')
endif

print*, 'not done'

end subroutine pert_model_copies

!------------------------------------------------------------------
!------------------------------------------------------------------
! Verify that the namelist was filled in correctly, and check
! that there are valid entries for the dart_kind.

subroutine verify_state_variables(nvar, varname, qty, update, in_domain)

integer,           intent(out) :: nvar
character(len=NF90_MAX_NAME), intent(out) :: varname(MAX_STATE_VARIABLES)
integer,           intent(out) :: qty(MAX_STATE_VARIABLES)
logical,           intent(out) :: update(MAX_STATE_VARIABLES)
character(len=9),  intent(out) :: in_domain(MAX_STATE_VARIABLES) ! assumes <=9 or 999

integer :: i
character(len=NF90_MAX_NAME) :: qty_str, update_str
character(len=256) :: string1, string2

if ( .not. module_initialized ) call static_init_model

nvar = 0
varloop: do i = 1, MAX_STATE_VARIABLES

   if ( wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -3) == 'NULL ' ) exit varloop ! Found end of list. !HK do you need to test all columns for ' '?

   varname(i)   = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -3))
   qty_str      = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -2))
   update_str   = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -1))
   in_domain(i) = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i   ))

   call to_upper(update_str)
   

   ! Make sure DART qty is valid
   qty(i) = get_index_for_quantity(qty_str)
   if( qty(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(qty_str)
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif
   
   ! Make sure the update variable has a valid name
   select case (update_str)
      case ('UPDATE')
         update(i) = .true.
      case ('NO_COPY_BACK')
         update(i) = .false.
      case default
         write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(string2,'(6A)') 'you provided : ', trim(varname(i)), ', ', trim(qty_str), ', ', trim(update_str)
         call error_handler(E_ERR,'verify_state_variables',string1, text2=string2)
   end select
   
   nvar = nvar + 1

enddo varloop

end subroutine verify_state_variables

!------------------------------------------------------------------
! This is assuming that the number of domains <=9
! do while loop could be replaced with intrinsic scan
logical function variable_is_on_domain(domain_id_string, id)

integer,           intent(in) :: id
character(len=*),  intent(in) :: domain_id_string

integer                       :: domain_int, i

variable_is_on_domain = .false.

! if '999' then counts all domains
if ( trim(domain_id_string) == '999' ) then
   variable_is_on_domain = .true.
else
i = 1
   do while ( domain_id_string(i:i) /= ' ' )
   read(domain_id_string(i:i),'(i1)') domain_int
   if ( domain_int == id ) variable_is_on_domain = .true.
   i = i+1
   enddo
endif

end function variable_is_on_domain

!------------------------------------------------------------------
function compute_geometric_height(geopot, lat)

real(r8), intent(in)  :: geopot
real(r8), intent(in)  :: lat
real(r8)              :: compute_geometric_height


real(digits12) :: pi2, latr
real(digits12) :: semi_major_axis, semi_minor_axis, grav_polar, grav_equator
real(digits12) :: earth_omega, grav_constant, flattening, somigliana
real(digits12) :: grav_ratio, sin2, termg, termr, grav, eccentricity

!  Parameters below from WGS-84 model software inside GPS receivers.
parameter(semi_major_axis = 6378.1370d3)    ! (m)
parameter(semi_minor_axis = 6356.7523142d3) ! (m)
parameter(grav_polar = 9.8321849378)        ! (m/s2)
parameter(grav_equator = 9.7803253359)      ! (m/s2)
parameter(earth_omega = 7.292115d-5)        ! (rad/s)
parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
parameter(eccentricity = 0.081819d0)        ! unitless
parameter(pi2 = 3.14159265358979d0/180.d0)

!  Derived geophysical constants
parameter(flattening = (semi_major_axis-semi_minor_axis) / semi_major_axis)

parameter(somigliana = (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

parameter(grav_ratio = (earth_omega*earth_omega * &
                semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)

! To use geopotential height uncomment the following two lines:
!compute_geometric_height = geopot
!return

latr = lat * (pi2)        ! in radians
sin2  = sin(latr) * sin(latr)
termg = grav_equator * ( (1.d0+somigliana*sin2) / &
        sqrt(1.d0-eccentricity*eccentricity*sin2) )
termr = semi_major_axis / (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

compute_geometric_height = (termr*geopot) / ( (termg/grav) * termr - geopot )


end function compute_geometric_height






!===================================================================
! End of model_mod
!===================================================================
end module model_mod

