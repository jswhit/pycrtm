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

# header file with constants
cdef extern from 'pycrtm_interface.h':
   enum: CRTM_STRLEN

cdef double DIFFUSIVITY_ANGLE  = 53.130102354156 # from CRTM_Parameters.f90

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

# TODO:
# crtm_surface_type
# crtm_atmosphere_type
# crtm_rtsolution_type
# crtm_options_type
# crtm_cloud_type ? (part of crtm_atmosphere_type)
# crtm_aerosol_type ? (part of crtm_atmosphere_type)
# crtm_SensorData_type ? (part of crtm_surface_type)
# ssu_input_type? (part of crtm_options_type)
# zeeman_input_type? (part of crtm_options_type)
