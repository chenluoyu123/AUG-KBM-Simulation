# this script is used to load the nonlin data
# to the output tree
# setup for cori
#sever='cori'
#server='xiang@cori.nersc.gov'
#tunnel='jianx@cybele.gat.com:22'
server='mingyue' #can do nersc and shenma as well
# setup for mingyue-hf
if server=='mingyue':
    servername='mingyue-hf'
    server='jianx@hfeshell.nscc-hf.cn:65062'
    tunnel=''
elif server=='nersc':
    servername='perlmutter'
    server='xiang@perlmutter.nersc.gov'
    tunnel=''
else:
    servername='shenma'
pltnl=root['SETTINGS']['PLOTS']['NL']
npath=len(pltnl['case_load'].keys())
outputs=root['OUTPUTS']
icgyro=root['SETTINGS']['SETUP']['icgyro']
# set the default server
for k in range(npath):
    keyname=pltnl['case_load'].keys()[k]
#    outputs[keyname]=OMFITcgyro([pltnl['path'][keyname],server,tunnel],extra_files=['bin.cgyro.kxky_apar','bin.cgyro.kxky_bpar'])
    if servername!='shenma':
        if icgyro==1:
            outputs[keyname]=OMFITcgyro([pltnl['case_load'][keyname],server,tunnel],extra_files=['bin.cgyro.kxky_apar','bin.cgyro.kxky_bpar','bin.cgyro.kxky_n'])
        else:
            outputs[keyname]=OMFITgyro([pltnl['case_load'][keyname],server,tunnel],extra_files=['bin.gyro.kxkyspec','bin.gyro.moment_n','bin.gyro.moment_u','bin.gyro.moment_e','bin.gyro.moment_0'])
    
    else:
        if icgyro==1:
            outputs[keyname]=OMFITcgyro(pltnl['case_load'][keyname])
        else:
            outputs[keyname]=OMFITgyro([pltnl['case_load'][keyname]],extra_files=['bin.gyro.kxkyspec','bin.gyro.moment_n','bin.gyro.moment_u','bin.gyro.moment_e','bin.gyro.moment_0'])
