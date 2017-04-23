module pycrtm_interface

use iso_c_binding, only: c_char, c_int, c_double, c_null_char, &
 c_ptr, c_loc, c_f_pointer
use crtm_module, only: crtm_init, crtm_destroy, crtm_channelinfo_type, &
 success, strlen, crtm_channelinfo_inspect, crtm_geometry_inspect, &
 crtm_geometry_setvalue, crtm_geometry_type
implicit none 

contains

! get CRTM string length
subroutine get_strlen(lenstr) bind(c)
  integer(c_int), intent(out) :: lenstr
  lenstr = strlen
end subroutine get_strlen

! initialize Channel_Info type, read in coefficient files for a given instrument/sensor
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

  Load_CloudCoeff=iload_cloudcoeff
  Load_AerosolCoeff=iload_aerosolcoeff
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
   channel_info % Process_Channel = process_channel
end subroutine channelinfo_set_process_channel

! get crtm_channel_info derived type member Process_Channel
subroutine channelinfo_get_process_channel(channel_infop,process_channel, n) bind (c)
   integer(c_int), intent(in) :: n
   type(c_ptr), intent(in) :: channel_infop
   integer(c_int), intent(out), dimension(n) :: process_channel
   type (crtm_channelinfo_type), pointer :: channel_info
   call c_f_pointer(channel_infop, channel_info)
   process_channel = channel_info % Process_Channel
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

! initialize Geometry_Info type
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
  character(c_char), intent(out) :: stringc(strlen+1)
  integer j,n
  n = len_trim(stringf)   
  do j=1,n    
    stringc(j) = stringf(j:j)   
  end do
  stringc(n+1) = c_null_char
end subroutine copy_string_ftoc

end module pycrtm_interface
