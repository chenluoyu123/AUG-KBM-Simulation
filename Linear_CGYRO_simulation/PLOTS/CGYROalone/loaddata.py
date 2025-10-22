# this script is used to load the nonlin data
# to the output tree
# setup for cori
#sever='cori'
#server='xiang@cori.nersc.gov'
#tunnel='jianx@cybele.gat.com:22'
# setup for mingyue-hf
servername='mingyue-hf'
server='jianx@hfeshell.nscc-hf.cn:65062'
tunnel=''
pltnl=root['SETTINGS']['PLOTS']['NL']
npath=len(pltnl['case_load'].keys())
outputs=root['OUTPUTS']
# set the default server

for k in range(npath):
    keyname=pltnl['case_load'].keys()[k]
#    outputs[keyname]=OMFITcgyro([pltnl['path'][keyname],server,tunnel],extra_files=['bin.cgyro.kxky_apar','bin.cgyro.kxky_bpar'])
    if servername!='shenma':
        outputs[keyname]=OMFITcgyro([pltnl['case_load'][keyname],server,tunnel],extra_files=['bin.cgyro.kxky_apar','bin.cgyro.kxky_bpar','bin.cgyro.kxky_n'])
    else:
        outputs[keyname]=OMFITcgyro(pltnl['case_load'][keyname])
