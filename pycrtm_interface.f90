module pycrtm_interface

use iso_c_binding, only: c_char, c_int, c_double, c_null_char, &
 c_long, c_ptr, c_loc, c_f_pointer
use crtm_module, only: crtm_init, crtm_destroy, crtm_channelinfo_type, &
 success, strlen, crtm_channelinfo_inspect, crtm_geometry_inspect, &
 crtm_geometry_setvalue, crtm_geometry_type, crtm_geometry_destroy, &
 crtm_options_type,crtm_options_create,crtm_options_inspect
implicit none 

contains

! get CRTM string length
subroutine get_strlen(lenstr) bind(c)
  integer(c_int), intent(out) :: lenstr
  lenstr = strlen
end subroutine get_strlen

! initialize crtm_channelinfo_type, read in coefficient files for a given instrument/sensor
subroutine init_crtm(nchanl,isis,iload_cloudcoeff,iload_aerosolcoeff,&
                     crtm_coeffs_path,channel_infop) bind(c)
!   input argument list:
!     nchanl - (int) number of channels 
!     isis   - (char*strlen) instrument/sensor character string 
!     iload_cloudcoeff - (int) 1 to load cloud coeffs
!     iload_aerosolcoeff - (int) 1 to load aerosol coeffs
!     crtm_coeffs_path - (char*256) path to CRTM coeffs files
!   output:
!     channel_infop - (c_ptr) opaque pointer to channel_info type
! input variables.
  integer(c_int),intent(in) :: nchanl,iload_cloudcoeff,iload_aerosolcoeff
  character(c_char), intent(in) :: isis(strlen)
  character(c_char), intent(in) :: crtm_coeffs_path(257)
! output variables
  type(c_ptr),intent(out) :: channel_infop
! local variables.
  character(len=strlen) :: isis_f
  integer :: error_status
  logical :: ice,Load_AerosolCoeff,Load_CloudCoeff
  type(crtm_channelinfo_type), pointer :: channel_info
  character(len=256) :: crtm_coeffs_path_f
! local parameters
  character(len=*), parameter :: myname_='pycrtm_interface*init_crtm'

  ! convert integer to logical
  Load_CloudCoeff=iload_cloudcoeff /= 0
  Load_AerosolCoeff=iload_aerosolcoeff /= 0
  call copy_string_ctof(isis,isis_f)
  call copy_string_ctof(crtm_coeffs_path,crtm_coeffs_path_f)
  allocate(channel_info)

! Initialize radiative transfer
  write(6,*)myname_,': crtm_init() on path "'//trim(crtm_coeffs_path_f)//'"'
  error_status = crtm_init_wrap(isis_f,channel_info,&
     Load_CloudCoeff=Load_CloudCoeff,Load_AerosolCoeff=Load_AerosolCoeff, &
     File_Path = crtm_coeffs_path_f )
  if (error_status /= success) then
     print *,myname_,':  ***ERROR*** crtm_init error_status=',error_status,&
        '   TERMINATE PROGRAM EXECUTION'
     stop
  endif
  channel_infop = c_loc(channel_info)
end subroutine init_crtm

! wrapper that accepts scalar derived type instead of array of derived types
FUNCTION crtm_init_wrap( &
  Sensor_ID         , &  ! Input
  ChannelInfo       , &  ! Output
  File_Path         , &  ! Optional input
  Load_CloudCoeff   , &  ! Optional input
  Load_AerosolCoeff ) &  ! Optional input
RESULT( error_status )
  ! Arguments
  CHARACTER(*)               , INTENT(IN)  :: Sensor_ID
  TYPE(CRTM_ChannelInfo_type), INTENT(OUT) :: ChannelInfo
  CHARACTER(*),      OPTIONAL, INTENT(IN)  :: File_Path
  LOGICAL     ,      OPTIONAL, INTENT(IN)  :: Load_CloudCoeff
  LOGICAL     ,      OPTIONAL, INTENT(IN)  :: Load_AerosolCoeff
  ! Function result
  INTEGER :: error_status
  ! local variables.
  type(crtm_channelinfo_type), dimension(1) :: channelinfo_array
  character(len=len(Sensor_ID)), dimension(1) :: sensor_id_array
  sensor_id_array(1) = sensor_id
  error_status = crtm_init(sensor_id_array,channelinfo_array,&
     Load_CloudCoeff=Load_CloudCoeff,Load_AerosolCoeff=Load_AerosolCoeff, &
     File_Path = file_path )
  channelinfo = channelinfo_array(1)
END FUNCTION crtm_init_wrap

FUNCTION crtm_destroy_wrap( &
  ChannelInfo       ) &  ! Input/Output
RESULT( error_status )
  ! Arguments
  TYPE(CRTM_ChannelInfo_type), INTENT(INOUT) :: ChannelInfo
  ! Function result
  INTEGER :: error_status
  ! local variables.
  type(crtm_channelinfo_type), dimension(1) :: channelinfo_array
  channelinfo_array(1) = ChannelInfo
  error_status = crtm_destroy(channelinfo_array)
END FUNCTION crtm_destroy_wrap

! print info in Channel_Info type
subroutine print_channelinfo(channel_infop) bind(c)
  type(c_ptr), intent(in) :: channel_infop
  type (crtm_channelinfo_type), pointer :: channel_info
  call c_f_pointer(channel_infop, channel_info)
  call crtm_channelinfo_inspect( channel_info )
end subroutine print_channelinfo

! set crtm_channel_info derived type member n_Channels
subroutine channelinfo_set_n_channels(channel_infop, n_Channels) bind(c)
   type(c_ptr), intent(out) :: channel_infop
   integer(c_int), intent(in) :: n_Channels
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % n_Channels = n_Channels
end subroutine channelinfo_set_n_channels

! get crtm_channel_info derived type member n_Channels
subroutine channelinfo_get_n_channels(channel_infop, n_Channels) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out) :: n_Channels
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   n_Channels = channel_info % n_Channels 
end subroutine channelinfo_get_n_channels

! set crtm_channel_info derived type member Sensor_Type
subroutine channelinfo_set_sensor_type(channel_infop, Sensor_Type) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in) :: Sensor_Type
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % Sensor_Type = Sensor_Type
end subroutine channelinfo_set_sensor_type

! get crtm_channel_info derived type member Sensor_Type
subroutine channelinfo_get_sensor_type(channel_infop, Sensor_Type) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out) :: Sensor_Type
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   Sensor_Type = channel_info % Sensor_Type 
end subroutine channelinfo_get_sensor_type

! set crtm_channel_info derived type member Sensor_Index
subroutine channelinfo_set_sensor_index(channel_infop, Sensor_Index) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in) :: Sensor_Index
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % Sensor_Index = Sensor_Index
end subroutine channelinfo_set_sensor_index

! get crtm_channel_info derived type member Sensor_Index
subroutine channelinfo_get_sensor_index(channel_infop, Sensor_Index) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out) :: Sensor_Index
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   Sensor_Index = channel_info % Sensor_Index 
end subroutine channelinfo_get_sensor_index

! set crtm_channel_info derived type member WMO_Satellite_ID
subroutine channelinfo_set_wmo_satellite_id(channel_infop, WMO_Satellite_ID) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in) :: WMO_Satellite_ID
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % WMO_Satellite_ID = WMO_Satellite_ID
end subroutine channelinfo_set_wmo_satellite_id

! get crtm_channel_info derived type member WMO_Satellite_ID
subroutine channelinfo_get_wmo_satellite_id(channel_infop, WMO_Satellite_ID) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out) :: WMO_Satellite_ID
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   WMO_Satellite_ID = channel_info % WMO_Satellite_ID 
end subroutine channelinfo_get_wmo_satellite_id

! set crtm_channel_info crtm_channel_info derived type member WMO_Sensor_ID
subroutine channelinfo_set_wmo_sensor_id(channel_infop, WMO_Sensor_ID) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in) :: WMO_Sensor_ID
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % WMO_Sensor_ID = WMO_Sensor_ID
end subroutine channelinfo_set_wmo_sensor_id

! get crtm_channel_info derived type member WMO_Sensor_ID
subroutine channelinfo_get_wmo_sensor_id(channel_infop, WMO_Sensor_ID) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out) :: WMO_Sensor_ID
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   WMO_Sensor_ID = channel_info % WMO_Sensor_ID 
end subroutine channelinfo_get_wmo_sensor_id

! get crtm_channel_info  derived type member Sensor_ID
subroutine channelinfo_get_sensor_id(channel_infop,name) bind (c)
   type(c_ptr), intent(in) :: channel_infop
   character(c_char), intent(out) :: name(strlen+1)
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   call copy_string_ftoc(channel_info % Sensor_ID,name)
end subroutine channelinfo_get_sensor_id

! set crtm_channel_info derived type member Sensor_ID
subroutine channelinfo_set_sensor_id(channel_infop, name) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   character(c_char), intent(in) :: name(strlen+1)
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   call copy_string_ctof(name,channel_info % Sensor_ID)
end subroutine channelinfo_set_sensor_id

! set crtm_channel_info derived type member Sensor_Channel
subroutine channelinfo_set_sensor_channel(channel_infop, sensor_channel, n) bind(c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in), dimension(n) :: sensor_channel
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % Sensor_Channel = sensor_channel
end subroutine channelinfo_set_sensor_channel

! get crtm_channel_info derived type member Sensor_Channel
subroutine channelinfo_get_sensor_channel(channel_infop,sensor_channel, n) bind (c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out), dimension(n) :: sensor_channel
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   sensor_channel = channel_info % Sensor_Channel
end subroutine channelinfo_get_sensor_channel

! set crtm_channel_info derived type member Channel_Index
subroutine channelinfo_set_channel_index(channel_infop, channel_index, n) bind(c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in), dimension(n) :: channel_index
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % Channel_Index = channel_index
end subroutine channelinfo_set_channel_index

! get crtm_channel_info derived type member Channel_Index
subroutine channelinfo_get_channel_index(channel_infop,channel_index, n) bind (c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out), dimension(n) :: channel_index
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_index = channel_info % Channel_Index
end subroutine channelinfo_get_channel_index

! set crtm_channel_info derived type member Process_Channel
subroutine channelinfo_set_process_channel(channel_infop, process_channel, n) bind(c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(in), dimension(n) :: process_channel
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   channel_info % Process_Channel = process_channel /= 0
end subroutine channelinfo_set_process_channel

! get crtm_channel_info derived type member Process_Channel
subroutine channelinfo_get_process_channel(channel_infop,process_channel, n) bind (c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out), dimension(n) :: process_channel
   type (crtm_channelinfo_type), pointer :: channel_info
   integer i
   call c_f_pointer(channel_infop, channel_info)
   do i=1,n
      if (channel_info % Process_Channel(i)) then
          process_channel(i) = 1
      else     
          process_channel(i) = 0
      endif
   enddo
end subroutine channelinfo_get_process_channel

! deallocate crtm_channelinfo_type
subroutine destroy_channelinfo(channel_infop) bind(c)
   type(c_ptr), intent(in) :: channel_infop
   integer error_status
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   error_status = crtm_destroy_wrap(channel_info)
   if (error_status /= success) then
      write(6,*) ' ***ERROR*** crtm_destroy,error_status=',error_status
     stop
   endif
   deallocate(channel_info % process_channel)
   deallocate(channel_info % sensor_channel)
   deallocate(channel_info % channel_index)
   deallocate(channel_info)
end subroutine destroy_channelinfo

! initialize crtm_geometry_type
subroutine init_geometry(ifov,longitude,latitude,&
                surface_altitude,sensor_scan_angle,&
                sensor_zenith_angle,sensor_azimuth_angle,&
                source_zenith_angle,source_aziumth_angle,flux_zenith_angle,&
                year,month,day,geometryp) bind(c)
  integer(c_int), intent(in) :: ifov,year,month,day              
  type(c_ptr), intent(out) :: geometryp
  real(c_double), intent(in) :: longitude,latitude,surface_altitude,&
  sensor_scan_angle,sensor_zenith_angle,sensor_azimuth_angle,&
  source_zenith_angle,source_aziumth_angle,flux_zenith_angle
  type (crtm_geometry_type), pointer :: geometry
  allocate(geometry)
  call crtm_geometry_setvalue(geometry,ifov,longitude,latitude,&
               surface_altitude,sensor_scan_angle,&
               sensor_zenith_angle,sensor_azimuth_angle,&
               source_zenith_angle,source_aziumth_angle,flux_zenith_angle,&
               year,month,day)
  geometryp = c_loc(geometry)
end subroutine init_geometry

! print info in Geometry type
subroutine print_geometry(geometryp) bind(c)
  type(c_ptr), intent(in) :: geometryp
  type (crtm_geometry_type), pointer :: geometry
  call c_f_pointer(geometryp, geometry)
  call crtm_geometry_inspect( geometry )
end subroutine print_geometry

! get crtm_geometry_type derived type member ifov
subroutine geometry_get_ifov(geometryp, ifov) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(out) :: ifov
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   ifov = geometry % ifov
end subroutine geometry_get_ifov

! set crtm_geometry_type derived type member ifov
subroutine geometry_set_ifov(geometryp, ifov) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(in) :: ifov
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % ifov = ifov
end subroutine geometry_set_ifov

! get crtm_geometry_type derived type member year
subroutine geometry_get_year(geometryp, year) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(out) :: year
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   year = geometry % year
end subroutine geometry_get_year

! set crtm_geometry_type derived type member year
subroutine geometry_set_year(geometryp, year) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(in) :: year
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % year = year
end subroutine geometry_set_year

! get crtm_geometry_type derived type member month
subroutine geometry_get_month(geometryp, month) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(out) :: month
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   month = geometry % month
end subroutine geometry_get_month

! set crtm_geometry_type derived type member month
subroutine geometry_set_month(geometryp, month) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(in) :: month
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % month = month
end subroutine geometry_set_month

! get crtm_geometry_type derived type member day
subroutine geometry_get_day(geometryp, day) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(out) :: day
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   day = geometry % day
end subroutine geometry_get_day

! set crtm_geometry_type derived type member day
subroutine geometry_set_day(geometryp, day) bind(c)
   type(c_ptr), intent(in) :: geometryp
   integer(c_int), intent(in) :: day
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % day = day
end subroutine geometry_set_day

! get crtm_geometry_type derived type member latitude
subroutine geometry_get_latitude(geometryp, latitude) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: latitude
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   latitude = geometry % latitude
end subroutine geometry_get_latitude

! set crtm_geometry_type derived type member latitude
subroutine geometry_set_latitude(geometryp, latitude) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: latitude
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % latitude = latitude
end subroutine geometry_set_latitude

! deallocate crtm_geometry_type
subroutine destroy_geometry(geometryp) bind(c)
   type(c_ptr), intent(in) :: geometryp
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   call crtm_geometry_destroy(geometry)
   deallocate(geometry)
end subroutine destroy_geometry

! get crtm_geometry_type derived type member longitude
subroutine geometry_get_longitude(geometryp, longitude) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: longitude
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   longitude = geometry % longitude
end subroutine geometry_get_longitude

! set crtm_geometry_type derived type member longitude
subroutine geometry_set_longitude(geometryp, longitude) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: longitude
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % longitude = longitude
end subroutine geometry_set_longitude

! get crtm_geometry_type derived type member surface_altitude
subroutine geometry_get_surface_altitude(geometryp, surface_altitude) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: surface_altitude
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   surface_altitude = geometry % surface_altitude
end subroutine geometry_get_surface_altitude

! set crtm_geometry_type derived type member surface_altitude
subroutine geometry_set_surface_altitude(geometryp, surface_altitude) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: surface_altitude
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % surface_altitude = surface_altitude
end subroutine geometry_set_surface_altitude

! get crtm_geometry_type derived type member sensor_scan_angle
subroutine geometry_get_sensor_scan_angle(geometryp, sensor_scan_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: sensor_scan_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   sensor_scan_angle = geometry % sensor_scan_angle
end subroutine geometry_get_sensor_scan_angle

! set crtm_geometry_type derived type member sensor_scan_angle
subroutine geometry_set_sensor_scan_angle(geometryp, sensor_scan_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: sensor_scan_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % sensor_scan_angle = sensor_scan_angle
end subroutine geometry_set_sensor_scan_angle

! get crtm_geometry_type derived type member sensor_zenith_angle
subroutine geometry_get_sensor_zenith_angle(geometryp, sensor_zenith_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: sensor_zenith_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   sensor_zenith_angle = geometry % sensor_zenith_angle
end subroutine geometry_get_sensor_zenith_angle

! set crtm_geometry_type derived type member sensor_zenith_angle
subroutine geometry_set_sensor_zenith_angle(geometryp, sensor_zenith_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: sensor_zenith_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % sensor_zenith_angle = sensor_zenith_angle
end subroutine geometry_set_sensor_zenith_angle

! get crtm_geometry_type derived type member sensor_azimuth_angle
subroutine geometry_get_sensor_azimuth_angle(geometryp, sensor_azimuth_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: sensor_azimuth_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   sensor_azimuth_angle = geometry % sensor_azimuth_angle
end subroutine geometry_get_sensor_azimuth_angle

! set crtm_geometry_type derived type member sensor_azimuth_angle
subroutine geometry_set_sensor_azimuth_angle(geometryp, sensor_azimuth_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: sensor_azimuth_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % sensor_azimuth_angle = sensor_azimuth_angle
end subroutine geometry_set_sensor_azimuth_angle

! get crtm_geometry_type derived type member source_zenith_angle
subroutine geometry_get_source_zenith_angle(geometryp, source_zenith_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: source_zenith_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   source_zenith_angle = geometry % source_zenith_angle
end subroutine geometry_get_source_zenith_angle

! set crtm_geometry_type derived type member source_zenith_angle
subroutine geometry_set_source_zenith_angle(geometryp, source_zenith_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: source_zenith_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % source_zenith_angle = source_zenith_angle
end subroutine geometry_set_source_zenith_angle

! get crtm_geometry_type derived type member source_azimuth_angle
subroutine geometry_get_source_azimuth_angle(geometryp, source_azimuth_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: source_azimuth_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   source_azimuth_angle = geometry % source_azimuth_angle
end subroutine geometry_get_source_azimuth_angle

! set crtm_geometry_type derived type member source_azimuth_angle
subroutine geometry_set_source_azimuth_angle(geometryp, source_azimuth_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: source_azimuth_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % source_azimuth_angle = source_azimuth_angle
end subroutine geometry_set_source_azimuth_angle

! get crtm_geometry_type derived type member flux_zenith_angle
subroutine geometry_get_flux_zenith_angle(geometryp, flux_zenith_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(out) :: flux_zenith_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   flux_zenith_angle = geometry % flux_zenith_angle
end subroutine geometry_get_flux_zenith_angle

! set crtm_geometry_type derived type member flux_zenith_angle
subroutine geometry_set_flux_zenith_angle(geometryp, flux_zenith_angle) bind(c)
   type(c_ptr), intent(in) :: geometryp
   real(c_double), intent(in) :: flux_zenith_angle
   type (crtm_geometry_type), pointer :: geometry
   call c_f_pointer(geometryp, geometry)
   geometry % flux_zenith_angle = flux_zenith_angle
end subroutine geometry_set_flux_zenith_angle

! initialize crtm_options_type
subroutine init_options(n_Channels,icheck_input, &
                iUse_Old_MWSSEM,iUse_Antenna_Correction, &
                iApply_NLTE_Correction,RT_Algorithm_Id, &
                Aircraft_Pressure,iUse_n_Streams, &
                n_Streams,iInclude_Scattering, &
                Channel,iUse_Emissivity,optionsp) bind(c)
   integer(c_long), intent(in) :: n_Channels,RT_Algorithm_Id,n_Streams,Channel
   integer(c_int), intent(in) :: icheck_input,iUse_Old_MWSSEM, &
      iUse_Antenna_Correction,iApply_NLTE_Correction,iInclude_Scattering, &
      iUse_n_Streams,iUse_Emissivity
   real(c_double), intent(in) :: Aircraft_Pressure
   type(c_ptr), intent(out) :: optionsp
   logical :: check_input,Use_Old_MWSSEM, &
      Use_Antenna_Correction,Apply_NLTE_Correction,Include_Scattering, &
      Use_Emissivity
   type (crtm_options_type), pointer :: options
   allocate(options)
   ! cast ints to fortran logicals
   check_input = icheck_input /= 0
   Use_Old_MWSSEM = iUse_Old_MWSSEM /= 0
   Use_Antenna_Correction = iUse_Antenna_Correction /= 0
   Apply_NLTE_Correction = iApply_NLTE_Correction /= 0
   Include_Scattering = iInclude_Scattering /= 0
   Use_Emissivity = iUse_Emissivity /= 0
   call crtm_options_create(options,int(n_Channels))
   options%check_input = check_input
   options%Use_Old_MWSSEM = Use_Old_MWSSEM
   options%Use_Antenna_Correction = Use_Antenna_Correction
   options%Apply_NLTE_Correction = Apply_NLTE_Correction
   options%Include_Scattering = Include_Scattering
   options%Use_Emissivity = Use_Emissivity
   options%RT_Algorithm_Id = RT_Algorithm_Id
   options%n_Streams = n_Streams
   options%Channel = Channel
   options%Aircraft_Pressure = Aircraft_Pressure
   optionsp = c_loc(options)
end subroutine init_options

! print info in Options type
subroutine print_options(optionsp) bind(c)
  type(c_ptr), intent(in) :: optionsp
  type (crtm_options_type), pointer :: options
  call c_f_pointer(optionsp, options)
  call crtm_options_inspect( options )
end subroutine print_options

! deallocate crtm_options_type
subroutine destroy_options(optionsp) bind(c)
   type(c_ptr), intent(in) :: optionsp
   type (crtm_options_type), pointer :: options
   call c_f_pointer(optionsp, options)
   !call crtm_options_destroy(options)
   deallocate(options)
end subroutine destroy_options

! utility functions
subroutine copy_string_ctof(stringc,stringf)
  ! utility function to convert c string to fortran string
  character(len=*), intent(out) :: stringf
  character(c_char), intent(in) :: stringc(:)
  integer j
  stringf = ''
  char_loop: do j=1,min(size(stringc),len(stringf))
     if (stringc(j)==c_null_char) exit char_loop
     stringf(j:j) = stringc(j)
  end do char_loop
end subroutine copy_string_ctof

subroutine copy_string_ftoc(stringf,stringc)
  ! utility function to convert c string to fortran string
  character(len=*), intent(in) :: stringf
  character(c_char), intent(out) :: stringc(:)
  integer j,n
  n = len_trim(stringf)   
  do j=1,n    
    stringc(j) = stringf(j:j)   
  end do
  stringc(n+1) = c_null_char
end subroutine copy_string_ftoc

end module pycrtm_interface
