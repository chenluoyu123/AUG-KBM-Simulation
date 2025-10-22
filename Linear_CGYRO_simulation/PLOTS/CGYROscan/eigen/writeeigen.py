# write a script to deploy the outputs
# the root['PLOTS']['CGYROscan']['eigencal.py'] needs to be run before runing of this script
# the data format is col
# casename,
# ky, kys, n_theta, n_r
# theta, phi, apar,  btheta, jpar
# r, phi_r, apar_r, btheta_r, apar_r
rootdir=OMFIT['CGYRO_SCAN']['OUTPUTS']['20_0.2']
plt=OMFIT['CGYRO_SCAN']['SETTINGS']['PLOTS']
#eigenout='/home/jianx/176125@2600/jx_forinputpro/deploy/ky0.2eigen.txt'
eigenout=root['SETTINGS']['DEPENDENCIES']['eigenout']
fid=open(eigenout,'w')
fid.write(str(plt['paraval_eigen'])+'_'+str(plt['ky_eigen']))
fid.write('\n')
fid.write(str(plt['ky_eigen'])+'    '+str(rootdir['kys'])+'    '+str(len(rootdir['theta']))+'    '+str(len(rootdir['r_arr'])))
fid.write('\n')
for k in range(len(rootdir['theta'])):
    line=str(rootdir['theta'][k])+'    '+\
         str(real(rootdir['phi'][k]))+'    '+str(imag(rootdir['phi'][k]))+'    '+\
         str(real(rootdir['apar'][k]))+'    '+str(imag(rootdir['apar'][k]))+'    '+\
         str(real(rootdir['btheta'][k]))+'    '+str(imag(rootdir['btheta'][k]))+'    '+\
         str(real(rootdir['jpar'][k]))+'    '+str(imag(rootdir['jpar'][k]))
    fid.write(line)
    fid.write('\n')
for k in range(len(rootdir['r_arr'])):
    line=str(rootdir['r_arr'][k])+'    '+\
         str(real(rootdir['phi_r_arr'][k]))+'    '+str(imag(rootdir['phi_r_arr'][k]))+'    '+\
         str(real(rootdir['apar_r_arr'][k]))+'    '+str(imag(rootdir['apar_r_arr'][k]))+'    '+\
         str(real(rootdir['btheta_r_arr'][k]))+'    '+str(imag(rootdir['btheta_r_arr'][k]))+'    '+\
         str(real(rootdir['jpar_r_arr'][k]))+'    '+str(imag(rootdir['jpar_r_arr'][k]))
    fid.write(line)
    fid.write('\n')
fid.close()
