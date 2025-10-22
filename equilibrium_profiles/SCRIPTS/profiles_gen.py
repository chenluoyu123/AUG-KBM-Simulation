##--------------------------
# What are the input files
##--------------------------
inputs=[]
if 'gfile' in root['INPUTS'].keys():
    inputs.append((root['INPUTS']['gfile'],'gfile'))
if 'statefile.nc' in root['INPUTS'].keys():
    inputs.append((root['INPUTS']['statefile.nc'],'statefile.nc'))
    executable ='profiles_gen -i statefile.nc -g gfile'
elif 'pfile' in root['INPUTS'].keys():
    inputs.append((root['INPUTS']['pfile'],'pfile'))
    executable ='profiles_gen -i pfile -g gfile'
elif 'input.profiles' in root['INPUTS'].keys():
    inputs.append((root['INPUTS']['input.profiles'],'input.profiles'))
    executable ='profiles_gen -i input.profiles'
else:
    raise('You must have statefile.nc or pfile or input.profiles for generation of input.gacode!!')
##----------------------
### output
##----------------------
# outputs=['input.profiles','input.profiles.geo']
outputs=['input.gacode']
#executable ='profiles_gen -i statefile.nc -g gfile'
#-----------------------
# Execute Profile_gen
#-----------------------
rmtdir=root['SETTINGS']['REMOTE_SETUP']['workDir']
ret_code=OMFITx.executable(root, inputs=inputs, outputs=outputs, executable=executable，remotedir=rmtdir)
#-----------------------
# load the results
#-----------------------
root['OUTPUTS']['Profiles_gen']['input.gacode']=OMFITgacode('input.gacode')
# root['OUTPUTS']['Profiles_gen']['input.profiles']=OMFITgacode('input.profiles')
# root['OUTPUTS']['Profiles_gen']['input.profiles.geo']=OMFITgacode('input.profiles.geo')
# # determine whether to generate the Er profile in the input.profiles
# if root['SETTINGS']['PHYSICS']['iGenEr']==1:
#     inputs=[root['OUTPUTS']['Profiles_gen']['input.profiles'],'input.profiles']
#     executable ='pbsMonitor -jq short -wt 00:30:00 -cn 8 -exe profiles_gen -vgen -i input.profiles -er 2 -vel 1 -in DC -ix 2'
#     ret_code=OMFITx.executable(root, inputs=inputs, outputs=outputs, executable=executable)
#     root['OUTPUTS']['Profiles_gen']['input.profiles']=OMFITgacode('vgen/input.profiles')
