import pycrtm
import numpy as np
isis='amsua_n15'
crtm_coeffs_path = '/Users/jsw/python/fixcrtm-2.2.3/'
nchanl = 15
iload_cloudcoeffs=1
iload_aerosolcoeffs=1
channel_info = pycrtm.Channel_Info(nchanl,isis,iload_cloudcoeffs,iload_aerosolcoeffs,crtm_coeffs_path)
print 
channel_info.show() # fortran side
print
print channel_info # python side
geometry = pycrtm.Geometry()
geometry.ifov = 10
geometry.year = 1993
geometry.month = 11
geometry.day = 21
geometry.latitude = -3.1415
geometry.show()
print geometry
