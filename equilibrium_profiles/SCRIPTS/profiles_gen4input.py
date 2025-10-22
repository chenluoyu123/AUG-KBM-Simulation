##--------------------------
# What are the input files
##--------------------------
inputs=[(root['INPUTS']['input.gacode'],'input.gacode')]
##----------------------
### output
##----------------------
outputs=['input.tglf.locpargen','input.cgyro.locpargen']
#-----------------------
# Execute Profile_gen
#-----------------------
# executable ='profiles_gen -i pfile -g gfile'
workdir=root['SETTINGS']['SETUP']['workDir']
rmtworkdir=root['SETTINGS']['REMOTE_SETUP']['workDir']
rmtserver=root['SETTINGS']['REMOTE_SETUP']['server']
tunnel=root['SETTINGS']['REMOTE_SETUP']['tunnel']
tgyroout=root['OUTPUTS']['TGYRO']
rho_arr=tgyroout['rho'][0][1:]
num_rho=len(rho_arr)
for k in range(num_rho):
    executable = 'profiles_gen -i input.gacode -loc_rho '+ str(rho_arr[k])
    ret_code=OMFITx.executable(root, inputs=inputs, outputs=outputs, workdir=workdir,remotedir=rmtworkdir,executable=executable,server=rmtserver,tunnel=tunnel,clean=True)
    root['OUTPUTS']['Profiles_gen']['input.tglf_'+str(k+1)]=OMFITgacode('input.tglf.locpargen')
    root['OUTPUTS']['Profiles_gen']['input.cgyro_' + str(k + 1)] = OMFITgacode('input.cgyro.locpargen')
