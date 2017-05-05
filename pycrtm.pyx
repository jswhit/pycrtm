import numpy as np
from numpy cimport ndarray 

# fortran functions with iso_c_binding interfaces.
cdef extern int get_strlen(int *strlen);
cdef extern int init_crtm(int *nchanl, char *isis, int *iload_cloudcoeffs, int *iload_aerosolcoeffs, char *crtm_coeffs_path, void *channelinfop);
cdef extern int destroy_channelinfo(void *channelinfop);
cdef extern int print_channelinfo(void *channelinfop);
cdef extern int channelinfo_set_n_channels(void *channelinfop, int *n);
cdef extern int channelinfo_get_n_channels(void *channelinfop, int *n);
cdef extern int channelinfo_set_sensor_index(void *channelinfop, int *n);
cdef extern int channelinfo_get_sensor_index(void *channelinfop, int *n);
cdef extern int channelinfo_set_sensor_type(void *channelinfop, int *n);
cdef extern int channelinfo_get_sensor_type(void *channelinfop, int *n);
cdef extern int channelinfo_set_wmo_satellite_id(void *channelinfop, int *n);
cdef extern int channelinfo_get_wmo_satellite_id(void *channelinfop, int *n);
cdef extern int channelinfo_set_wmo_sensor_id(void *channelinfop, int *n);
cdef extern int channelinfo_get_wmo_sensor_id(void *channelinfop, int *n);
cdef extern int channelinfo_get_sensor_id(void *channelinfop, char *name);
cdef extern int channelinfo_set_sensor_id(void *channelinfop, char *name);
cdef extern int channelinfo_get_sensor_channel(void *channelinfop, int *value, int *nchanl);
cdef extern int channelinfo_set_sensor_channel(void *channelinfop, int *value, int *nchanl);
cdef extern int channelinfo_get_channel_index(void *channelinfop, int *value, int *nchanl);
cdef extern int channelinfo_set_channel_index(void *channelinfop, int *value, int *nchanl);
cdef extern int channelinfo_get_process_channel(void *channelinfop, int *value, int *nchanl);
cdef extern int channelinfo_set_process_channel(void *channelinfop, int *value, int *nchanl);
cdef extern int print_geometry(void *geometryp);
cdef extern int init_geometry(int *ifov,double *longitude,double *latitude,
                double* surface_altitude,double *sensor_scan_angle,
                double *sensor_zenith_angle,double *sensor_azimuth_angle,
                double *source_zenith_angle,double *source_aziumth_angle,
                double *flux_zenith_angle,
                int *year,int *month,int *day, void *geometryp);
cdef extern int geometry_set_ifov(void *geometryp, int *n);
cdef extern int geometry_get_ifov(void *geometryp, int *n);
cdef extern int geometry_set_year(void *geometryp, int *n);
cdef extern int geometry_get_year(void *geometryp, int *n);
cdef extern int geometry_set_month(void *geometryp, int *n);
cdef extern int geometry_get_month(void *geometryp, int *n);
cdef extern int geometry_set_day(void *geometryp, int *n);
cdef extern int geometry_get_day(void *geometryp, int *n);
cdef extern int geometry_set_latitude(void *geometryp, double *x);
cdef extern int geometry_get_latitude(void *geometryp, double *x);
cdef extern int geometry_set_longitude(void *geometryp, double *x);
cdef extern int geometry_get_longitude(void *geometryp, double *x);
cdef extern int geometry_set_surface_altitude(void *geometryp, double *x);
cdef extern int geometry_get_surface_altitude(void *geometryp, double *x);
cdef extern int geometry_set_sensor_scan_angle(void *geometryp, double *x);
cdef extern int geometry_get_sensor_scan_angle(void *geometryp, double *x);
cdef extern int geometry_set_sensor_zenith_angle(void *geometryp, double *x);
cdef extern int geometry_get_sensor_zenith_angle(void *geometryp, double *x);
cdef extern int geometry_set_sensor_azimuth_angle(void *geometryp, double *x);
cdef extern int geometry_get_sensor_azimuth_angle(void *geometryp, double *x);
cdef extern int geometry_set_source_zenith_angle(void *geometryp, double *x);
cdef extern int geometry_get_source_zenith_angle(void *geometryp, double *x);
cdef extern int geometry_set_source_azimuth_angle(void *geometryp, double *x);
cdef extern int geometry_get_source_azimuth_angle(void *geometryp, double *x);
cdef extern int geometry_set_flux_zenith_angle(void *geometryp, double *x);
cdef extern int geometry_get_flux_zenith_angle(void *geometryp, double *x);
cdef extern int destroy_geometry(void *geometryp);
cdef extern int init_options(long *n_Channels,int *check_input,
                int *Use_Old_MWSSEM,int *Use_Antenna_Correction,
                int *Apply_NLTE_Correction,long *RT_Algorithm_Id,
                double *Aircraft_Pressure,int *Use_n_Streams,
                long *n_Streams,int *Include_Scattering,
                long *Channel,int *Use_Emissivity, void *optionsp);
cdef extern int print_options(void *optionsp);
cdef extern int destroy_options(void *optionsp);
cdef extern int set_zeeman_input(void *optionsp, double *field_strength,
        double *cos_theta8, double *cos_phi8, double *doppler_shift);
cdef extern int set_ssu_input(void *optionsp, double *time, double *cell_pressure,
        int *nchannel);

# header file with constants
cdef extern from 'pycrtm_interface.h':
   enum: CRTM_STRLEN

# from CRTM_Parameters.f90 (used for default initialization of CRTM types)
cdef double DIFFUSIVITY_ANGLE  = 53.130102354156 
cdef int RT_ADA = 56
cdef int RT_SOI = 168

#  When interfacing between Fortran and C, you will have to pass pointers to all
#  the variables you send to the Fortran function as arguments. Passing a variable
#  directly will probably crash Python.

def crtm_strlen():
    cdef int strlen
    get_strlen(&strlen)
    return strlen

# make sure constant in header file consistent with constant in library.
_crtm_strlen = crtm_strlen()
if _crtm_strlen != CRTM_STRLEN:
    raise ValueError('inconsistent value of CRTM_STRLEN')

# python version of crtm_channelinfo_type fortran derived type
cdef class ChannelInfo:
    cdef void *ptr
    cdef int nchanl
    def __init__(self, int nchanl, char *isis, int iload_cloudcoeff, int iload_aerosolcoeff, char *crtm_coeffs_path):
        self.nchanl = nchanl
        init_crtm(&nchanl, isis, &iload_cloudcoeff, &iload_aerosolcoeff,
                crtm_coeffs_path, &self.ptr)
    def show(self):
        print_channelinfo(&self.ptr)
    property n_Channels:
        """get and set n_Channels member of derived type"""
        def __get__(self):
            cdef int i
            channelinfo_get_n_channels(&self.ptr, &i)
            return i
        def __set__(self,int value):
            channelinfo_set_n_channels(&self.ptr, &value)
    property Sensor_Type:
        """get and set Sensor_Type member of derived type"""
        def __get__(self):
            cdef int i
            channelinfo_get_sensor_type(&self.ptr, &i)
            return i
        def __set__(self,int value):
            channelinfo_set_sensor_type(&self.ptr, &value)
    property WMO_Satellite_ID:
        """get and set WMO_Satellite_ID member of derived type"""
        def __get__(self):
            cdef int i
            channelinfo_get_wmo_satellite_id(&self.ptr, &i)
            return i
        def __set__(self,int value):
            channelinfo_set_wmo_satellite_id(&self.ptr, &value)
    property WMO_Sensor_ID:
        """get and set WMO_Sensor_ID member of derived type"""
        def __get__(self):
            cdef int i
            channelinfo_get_wmo_sensor_id(&self.ptr, &i)
            return i
        def __set__(self,int value):
            channelinfo_set_wmo_sensor_id(&self.ptr, &value)
    property Sensor_ID:
        """get and set Sensor_ID member of derived type"""
        def __get__(self):
            cdef char name[CRTM_STRLEN+1] # null char will be added
            channelinfo_get_sensor_id(&self.ptr, name)
            return name
        def __set__(self,char *value):
            channelinfo_set_sensor_id(&self.ptr, value)
    property Sensor_Index:
        """get and set Sensor_Index member of derived type"""
        def __get__(self):
            cdef int i
            channelinfo_get_sensor_index(&self.ptr, &i)
            return i
        def __set__(self,int value):
            channelinfo_set_sensor_index(&self.ptr, &value)
    property Sensor_Channel:
        """get and set Sensor_Channel member of derived type"""
        def __get__(self):
            cdef ndarray iarr = np.empty(self.nchanl,np.intc)
            channelinfo_get_sensor_channel(&self.ptr, <int *>iarr.data, &self.nchanl)
            return iarr
        def __set__(self,ndarray value):
            value = value.astype(np.intc)
            if value.size != self.nchanl:
                raise ValueError('cannot change the size of Sensor_Channel member')
            channelinfo_set_sensor_channel(&self.ptr, <int *>value.data, &self.nchanl)
    property Process_Channel:
        """get and set Process_Channel member of derived type"""
        def __get__(self):
            cdef ndarray iarr = np.empty(self.nchanl,np.intc)
            channelinfo_get_process_channel(&self.ptr, <int *>iarr.data, &self.nchanl)
            return iarr.astype(np.bool)
        def __set__(self,ndarray value):
            value = value.astype(np.intc)
            if value.size != self.nchanl:
                raise ValueError('cannot change the size of Process_Channel member')
            channelinfo_set_process_channel(&self.ptr, <int *>value.data, &self.nchanl)
    property Channel_Index:
        """get and set Channel_Index member of derived type"""
        def __get__(self):
            cdef ndarray iarr = np.empty(self.nchanl,np.intc)
            channelinfo_get_channel_index(&self.ptr, <int *>iarr.data, &self.nchanl)
            return iarr
        def __set__(self,ndarray value):
            value = value.astype(np.intc)
            if value.size != self.nchanl:
                raise ValueError('cannot change the size of Channel_Index member')
            channelinfo_set_channel_index(&self.ptr, <int *>value.data, &self.nchanl)
    def __dealloc__(self):
        destroy_channelinfo(&self.ptr)
    def __repr__(self):
        printlist = [' ChannelInfo OBJECT:\n']
        printlist.append('   n_Channels       : %s\n' % self.n_Channels)
        printlist.append('   Sensor_ID        : %s\n' % self.Sensor_ID)
        printlist.append('   Sensor_Type      : %s\n' % self.Sensor_Type)
        printlist.append('   WMO_Satellite_ID : %s\n' % self.WMO_Satellite_ID)
        printlist.append('   WMO_Sensor_ID    : %s\n' % self.WMO_Sensor_ID)
        printlist.append('   Sensor_Index     : %s\n' % self.Sensor_Index)
        printlist.append('   Channel#     Index     Process?\n')
        for n in range(self.n_Channels):
            printlist.append('        %s           %s        %s\n' % (self.Sensor_Channel[n],\
            self.Channel_Index[n],self.Process_Channel[n]))
        return ''.join(printlist)

# python version of crtm_geometry_type fortran derived type
cdef class Geometry:
    cdef void *ptr
    def __init__(self, int ifov=0, double longitude=0, double latitude=0,
                 double surface_altitude=0, double sensor_scan_angle=0,
                 double sensor_zenith_angle=0,
                 double sensor_azimuth_angle=999.9,
                 double source_zenith_angle=100, double source_aziumth_angle=0,
                 double flux_zenith_angle=DIFFUSIVITY_ANGLE, int year=2001,
                 int month=1, int day=1):
        init_geometry(&ifov,&longitude,&latitude,
                &surface_altitude,&sensor_scan_angle,
                &sensor_zenith_angle,&sensor_azimuth_angle,
                &source_zenith_angle,&source_aziumth_angle,&flux_zenith_angle,
                &year,&month,&day,&self.ptr)
    def show(self):
        print_geometry(&self.ptr)
    property ifov:
        """get and set ifov member of derived type"""
        def __get__(self):
            cdef int i
            geometry_get_ifov(&self.ptr, &i)
            return i
        def __set__(self,int value):
            geometry_set_ifov(&self.ptr, &value)
    property year:
        """get and set year member of derived type"""
        def __get__(self):
            cdef int i
            geometry_get_year(&self.ptr, &i)
            return i
        def __set__(self,int value):
            geometry_set_year(&self.ptr, &value)
    property month:
        """get and set month member of derived type"""
        def __get__(self):
            cdef int i
            geometry_get_month(&self.ptr, &i)
            return i
        def __set__(self,int value):
            geometry_set_month(&self.ptr, &value)
    property day:
        """get and set day member of derived type"""
        def __get__(self):
            cdef int i
            geometry_get_day(&self.ptr, &i)
            return i
        def __set__(self,int value):
            geometry_set_day(&self.ptr, &value)
    property latitude:
        """get and set latitude member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_latitude(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_latitude(&self.ptr, &value)
    property longitude:
        """get and set longitude member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_longitude(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_longitude(&self.ptr, &value)
    property surface_altitude:
        """get and set surface_altitude member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_surface_altitude(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_surface_altitude(&self.ptr, &value)
    property sensor_scan_angle:
        """get and set sensor_scan_angle member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_sensor_scan_angle(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_sensor_scan_angle(&self.ptr, &value)
    property sensor_zenith_angle:
        """get and set sensor_zenith_angle member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_sensor_zenith_angle(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_sensor_zenith_angle(&self.ptr, &value)
    property sensor_azimuth_angle:
        """get and set sensor_azimuth_angle member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_sensor_azimuth_angle(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_sensor_azimuth_angle(&self.ptr, &value)
    property source_zenith_angle:
        """get and set source_zenith_angle member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_source_zenith_angle(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_source_zenith_angle(&self.ptr, &value)
    property source_azimuth_angle:
        """get and set source_azimuth_angle member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_source_azimuth_angle(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_source_azimuth_angle(&self.ptr, &value)
    property flux_zenith_angle:
        """get and set flux_zenith_angle member of derived type"""
        def __get__(self):
            cdef double x
            geometry_get_flux_zenith_angle(&self.ptr, &x)
            return x
        def __set__(self,double value):
            geometry_set_flux_zenith_angle(&self.ptr, &value)
    def __repr__(self):
        printlist = [' Geometry OBJECT:\n']
        printlist.append('   FOV index           : %s\n' % self.ifov)
        printlist.append('   Longitude           : %s\n' % self.longitude)
        printlist.append('   Latitude            : %s\n' % self.latitude)
        printlist.append('   Surface altitude    : %s\n' %\
                self.surface_altitude)
        printlist.append('   Sensor scan angle   : %s\n' %\
                self.sensor_scan_angle)
        printlist.append('   Sensor zenith angle : %s\n' %\
                self.sensor_zenith_angle)
        printlist.append('   Sensor azimuth angle: %s\n' %\
                self.sensor_azimuth_angle)
        printlist.append('   Source zenith angle : %s\n' %\
                self.source_zenith_angle)
        printlist.append('   Source azimuth angle: %s\n' %\
                self.source_azimuth_angle)
        printlist.append('   Flux zenith angle   : %s\n' %\
                self.flux_zenith_angle)
        printlist.append('   Year                : %s\n' % self.year)
        printlist.append('   Month               : %s\n' % self.month)
        printlist.append('   Day                 : %s\n' % self.day)
        return ''.join(printlist)
    def __dealloc__(self):
        destroy_geometry(&self.ptr)

# python version of crtm_options_type fortran derived type
cdef class Options:
    cdef void *ptr
    def __init__(self,long n_Channels,check_input=True,Use_Old_MWSSEM=False,
             Use_Antenna_Correction=False,Apply_NLTE_Correction=True,
             long RT_Algorithm_Id=RT_ADA,double Aircraft_Pressure=-1,Use_n_Streams=False,
             long n_Streams=0,Include_Scattering=True,long Channel=0,Use_Emissivity=False):
        cdef int icheck_input,iUse_Old_MWSSEM,iUse_Antenna_Correction,\
                 iUse_n_Streams,iInclude_Scattering,iUse_Emissivity,iApply_NLTE_Correction
        # cast python bools to ints
        icheck_input = check_input
        iUse_Old_MWSSEM = Use_Old_MWSSEM
        iUse_Antenna_Correction = Use_Antenna_Correction
        iApply_NLTE_Correction = Apply_NLTE_Correction
        iUse_n_Streams = Use_n_Streams
        iInclude_Scattering = Include_Scattering
        iUse_Emissivity = Use_Emissivity
        init_options(&n_Channels,&icheck_input,
                &iUse_Old_MWSSEM,&iUse_Antenna_Correction,
                &iApply_NLTE_Correction,&RT_Algorithm_Id,
                &Aircraft_Pressure,&iUse_n_Streams,
                &n_Streams,&iInclude_Scattering,
                &Channel,&iUse_Emissivity, &self.ptr)
    def set_zeeman(self, field_strength=None, cos_theta8=None, cos_phi8=None,
                   doppler_shift=None):
        cdef double dfield_strength
        cdef double dcos_theta8
        cdef double dcos_phi8
        cdef double ddoppler_shift
        if field_strength is not None:
            dfield_strength = field_strength
            set_zeeman_input(&self.ptr,&dfield_strength,NULL,NULL,NULL)
        if cos_theta8 is not None:
            dcos_theta8 = cos_theta8
            set_zeeman_input(&self.ptr,NULL,&dfield_strength,NULL,NULL)
        if cos_phi8 is not None:
            dcos_phi8 = cos_phi8
            set_zeeman_input(&self.ptr,NULL,NULL,&dcos_phi8,NULL)
        if doppler_shift is not None:
            ddoppler_shift = doppler_shift
            set_zeeman_input(&self.ptr,NULL,NULL,NULL,&ddoppler_shift)
    def set_ssu(self, time=None, cell_pressure=None, nchannel=None):
        cdef double dcell_pressure, dtime
        cdef int ichannel
        if (cell_pressure is not None and nchannel is None) or\
           (cell_pressure is None and nchannel is not None):
            msg = 'must specify both cell_pressure and nchannel, or neither'
            raise ValueError(msg)
        if time is not None:
            dtime = time
            set_ssu_input(&self.ptr,&dtime,NULL,NULL)
        if cell_pressure is not None and nchannel is not None:
            dcell_pressure = cell_pressure; ichannel = nchannel
            set_ssu_input(&self.ptr,NULL,&dcell_pressure,&ichannel)
    def show(self):
        print_options(&self.ptr)
    def show(self):
        print_options(&self.ptr)
    def __dealloc__(self):
        destroy_options(&self.ptr)

# TODO:
# crtm_surface_type
# crtm_atmosphere_type
# crtm_rtsolution_type
# crtm_cloud_type ? (part of crtm_atmosphere_type)
# crtm_aerosol_type ? (part of crtm_atmosphere_type)
# crtm_SensorData_type ? (part of crtm_surface_type)
