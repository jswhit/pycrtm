import pycrtm
import numpy as np
isis='amsua_n15'
crtm_coeffs_path = '/Users/jsw/python/fixcrtm-2.2.3/'
crtm_coeffs_path = '/Volumes/Drobo/fixcrtm-2.2.3/'
nchanl = 15
iload_cloudcoeffs=1
iload_aerosolcoeffs=1
channel_info = pycrtm.ChannelInfo(nchanl,isis,iload_cloudcoeffs,iload_aerosolcoeffs,crtm_coeffs_path)
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
options = pycrtm.Options(15)
options.set_zeeman(field_strength=-999,cos_phi8=3)
options.set_ssu(time=-999)
options.show()
surface = pycrtm.Surface(15)
surface.Land_Coverage=1
surface.Land_Type=9
surface.show()
