import crtm_interface
isis='amsua_n15'
crtm_coeffs_path = '/scratch3/BMC/gsienkf/whitaker/gsi/EXP-enkflinhx/fix/crtm-2.2.3/'
n_channels = crtm_interface.crtm_get_nchannels(isis,crtm_coeffs_path)
print 'n_channels = ',n_channels
iload_cloudcoeffs=1
iload_aerosolcoeffs=1
sensor_id,sensor_type,wmo_satid,wmo_sensor_id = crtm_interface.crtm_initialize(isis,iload_cloudcoeffs,iload_aerosolcoeffs,crtm_coeffs_path)
print 'back in python:'
print 'sensor_id,sensor_type,wmo_satid,wmo_sensor_id',sensor_id,sensor_type,wmo_satid,wmo_sensor_id
print 'crtm strlen:',crtm_interface.crtm_strlen()
