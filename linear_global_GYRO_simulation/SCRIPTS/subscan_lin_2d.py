import numpy as np
import sys
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROscan/assist')
from cgyro_read_xj import *
print ("============== scan_cgyro ===========")
# specify the caseroot and casetag, note usually caseroot needs to have the key caseName
caseRoot=root['Cases'] # the root of the tag, defines where to store the scanning data
physics=root['SETTINGS']['PHYSICS']
phy2d=physics['2d']
caseName=physics['case_tag'] # the name of the tag, which indicates the name of the scanning
if not caseName in caseRoot.keys():
    caseRoot[caseName]=OMFITtree()
#else:
#    for item in caseRoot[caseName].keys():     #cover the previous calculation
#        del caseRoot[caseName][item]
kyarr=physics['kyarr']
Para_x=phy2d['Para_x']
Para_y=phy2d['Para_y']
Range_x=phy2d['Range_x']
Range_y=phy2d['Range_y']
setup=root['SETTINGS']['SETUP']
effnum=setup['effnum']
wkdir=setup['workDir'] # the local directory(not the remote one)
restart_mode=physics['restart_mode']
dir_list=[]
inputs_node_names={'input.cgyro':OMFITgacode}
inputs_node=root['INPUTS']
inputcgyro=inputs_node['input.cgyro']
if inputcgyro['PROFILE_MODEL']==2:
    inputs_node_names={'input.cgyro':OMFITgacode,'input.profiles':OMFITgacode,'input.profiles.geo':OMFITgacode}
inputs=[]
# OK , before doing the scan, we need to keep the original input.cgyro
inputs_node['input.cgyro_orig']=inputcgyro.duplicate()
# Gamma_E need to be set to 0 in the linear run
inputcgyro['GAMMA_E']=0
inputcgyro['NONLINEAR_FLAG']=0
if inputcgyro['PROFILE_MODEL']==2:
    inputcgyro['GAMMA_E_SCALE']=0.
# the outputs needs to be updated due to the update of CGYRO
base_loadOutputs={'bin.cgyro.aparb':OMFITgacode,'bin.cgyro.phib':OMFITgacode,'bin.cgyro.bparb':OMFITgacode,\
                  'bin.cgyro.geo':OMFITgacode,'bin.cgyro.kxky_phi':OMFITgacode,'bin.cgyro.ky_flux':OMFITgacode,\
                  'bin.cgyro.restart':OMFITgacode,'bin.cgyro.restart.old':OMFITgacode,\
                  'input.cgyro':OMFITgacode,'input.cgyro.gen':OMFITgacode,\
                  'out.cgyro.freq':OMFITgacode,'out.cgyro.info':OMFITgacode,'out.cgyro.time':OMFITgacode,'out.cgyro.timing':OMFITgacode,\
                  'out.cgyro.egrid':OMFITgacode,'out.cgyro.equilibrium':OMFITgacode,'out.cgyro.grids':OMFITgacode,\
                  'out.cgyro.hosts':OMFITgacode,'out.cgyro.memory':OMFITgacode,'out.cgyro.mpi':OMFITgacode,\
                  'out.cgyro.prec':OMFITgacode,'out.cgyro.version':OMFITgacode,\
                  'out.cgyro.tag':OMFITgacode,'run_log':OMFITgacode
                 }
loadOutputs={}
ncount_x=-1
# get the original delta_t and max_time
delta_t_orig=inputcgyro['DELTA_T']
max_time_orig=inputcgyro['MAX_TIME']
inputcgyro['NONLINEAR_FLAG']=0
if restart_mode==1:
    root['PLOTS']['CGYROscan']['assist']['changeOMFITcgyrogacode.py'].run()
for para_x_val in Range_x:
    ncount_x=ncount_x+1
    inputcgyro[Para_x]=para_x_val
    nPara_x=2
    while 'Para_x'+str(nPara_x) in phy2d.keys() and 'Range_x'+str(nPara_x) in phy2d.keys():
        inputcgyro[phy2d['Para_x'+str(nPara_x)]]=phy2d['Range_x'+str(nPara_x)][ncount_x]
        nPara_x=nPara_x+1
    ncount_y=-1
    for para_y_val in Range_y:
        printi('running the calculation on '+Para_x+ '='+num2str_xj(para_x_val,effnum)+' and '+Para_y+' = '+num2str_xj(para_y_val,effnum))
        ncount_y=ncount_y+1
        inputcgyro[Para_y]=para_y_val
        nPara_y=2
        while 'Para_y'+str(nPara_y) in phy2d.keys() and 'Range_y'+str(nPara_y) in phy2d.keys():
            inputcgyro[phy2d['Para_y'+str(nPara_y)]]=phy2d['Range_y'+str(nPara_y)][ncount_y]
            nPara_y=nPara_y+1
        for ky in kyarr:
            new_dir=Para_x+'~'+num2str_xj(para_x_val,effnum)+'~'+Para_y+'~'+num2str_xj(para_y_val,effnum)+'~ky~'+num2str_xj(ky,effnum)
            dir_list.append(new_dir)
            inputcgyro['KY']=ky
            time_scheme=physics['time_scheme']
            if time_scheme[0]==1:  # determine whether to use the time scheme
                if ky>1:
                    inputcgyro['DELTA_T']=time_scheme[1]/ky
                    inputcgyro['MAX_TIME']=time_scheme[2]/ky
                else:
                    inputcgyro['DELTA_T']=delta_t_orig
                    inputcgyro['MAX_TIME']=max_time_orig
            ## ----------------- ##
            if not new_dir in caseRoot[caseName].keys():
                caseRoot[caseName][new_dir]=OMFITtree()
            if os.path.isfile(new_dir+r'.zip'):
                os.remove(new_dir+r'.zip')
            tmp_zip=zipfile.ZipFile(new_dir+r'.zip', mode='w') # construct a new zip file, and waiting for files to be written in 
            if restart_mode==0:
                inputs_node_names_key=inputs_node_names.keys()
            else:
#                inputs_node_names_key=inputs_node_names.keys()+base_loadOutputs.keys()
                inputs_node_names_key=caseRoot[caseName][new_dir].keys()
                # we need to consider the condition that we may add new paramters regime to scan or there is no data in previous scan
                for item in inputs_node_names.keys():
                    if not item in inputs_node_names_key:
                        inputs_node_names_key=append(inputs_node_names_key,item)
                        # inputs_node[item].save()
                        # caseRoot[caseName][new_dir][item]=copy.deepcopy(inputs_node[item])
            caseRoot[caseName][new_dir]['input.cgyro'] = inputcgyro.duplicate()  # update the input.cgyro in CaseRoot by root['INPUTS']['input.cgyro']
            # for item  in inputs_node_names_key:
            #     if restart_mode==0:
            #         inputs_node[item].save()
            #         caseRoot[caseName][new_dir][item]=copy.deepcopy(inputs_node[item])
    #        item='input.cgyro'
    #        inputcgyro.save()
    #        caseRoot[caseName][new_dir][item]=copy.deepcopy(inputcgyro)
            for item in caseRoot[caseName][new_dir].keys():
                tmp_zip.write(caseRoot[caseName][new_dir][item].filename, caseRoot[caseName][new_dir][item].filename.split(os.sep)[-1])  # os.sep is / for linux and \ for windows
            tmp_zip.close()
            caseRoot[caseName][new_dir]['zip']=OMFITpath(tmp_zip.filename) # the inputs are compressed as an zip file and loaded in the output directory
            inputs.append(caseRoot[caseName][new_dir]['zip'])              # then the zipfile are import as input
            # set OUTPUTS( a lot of dir ) for loadOutputs, all the parameters ranges are covered
            for item in base_loadOutputs:
                loadOutputs[new_dir+os.sep+item]=base_loadOutputs[item]
outputs=loadOutputs.keys()
#
executable='mv input.cgyro* input.cgyro;\n '+setup['executable']

#set job scipts
pbs_file=setup['workDir']+os.sep+'scan.pbs'
username=root['SETTINGS']['REMOTE_SETUP']['server'].split('@')[0]
ps_name='cgyro'

rmt_setup=root['SETTINGS']['REMOTE_SETUP']
## the bash_head below is used for runing on SHNEMA
if rmt_setup['serverPicker'][0:5]=='login' or rmt_setup['serverPicker']=='shenma':
    bash_head= \
    r'#!/bin/sh '+'\n'+ \
    r'#PBS -N '+ps_name +'\n'+ \
    r'#PBS -l nodes='+str(setup['num_nodes'])+':ppn='+str(setup['num_cores']) +'\n'+ \
    r'#PBS -j oe' +'\n'+ \
    r'#PBS -l walltime='+str(setup['wall_time']) +'\n'+ \
    r'#PBS -q '+setup['pbs_queue'] +'\n'+ \
    r'cd ${PBS_O_WORKDIR}' +'\n'+ \
    r'pwd'  +'\n'+ \
    r'NP=`cat ${PBS_NODEFILE}|wc -l`' +'\n'+ \
    '\n'+ \
    r'JOBID_FILE="JOBID_${PBS_JOBID}"' +'\n'+ \
    r'touch ${JOBID_FILE}' +'\n'
#########
# the bash_head below is suitable for kuafu
#    r'#SBATCH --workdir=./'+'\n' +\
else:
    num_nodes=setup['num_nodes']
    num_cores=setup['num_cores']
    bash_head= \
    r'#!/bin/sh '+'\n' + \
    r'#SBATCH -p '+setup['pbs_queue'] +'\n' +\
    r'#SBATCH -J '+ps_name +'\n' +\
    r'#SBATCH -t '+str(setup['wall_time']) +'\n' +\
    r'#SBATCH -o -o.out'+'\n' +\
    r'#SBATCH -e -e.out'+'\n' +\
    r'#SBATCH --ntasks '+str(num_nodes*num_cores) +'\n'
    if setup['server']=='cori':
        bash_head=bash_head + r'#SBATCH -C knl,quad,cache' +'\n' 
###########
bash_content= \
'dir_list=('+' '.join(dir_list)+')' +'\n'+ \
r'for i in ${dir_list[@]}; ' +'\n'+  \
'do '  +'\n'+  \
r'  while [[ $( ps -u '+username+r' |grep '+ps_name+r' | wc -l ) -gt '+str(setup['num_cores']-1)+' ]]; do ' +'\n' + \
r'    sleep 2s' +'\n' +\
r'  done ' +'\n'+ \
r'  if [ -d ${i} ]; then  ' +'\n'+\
r'    echo "${i} exist. Do not unzip this case"   ' +'\n'+\
r'  cd $i' +'\n'+\
r'  else'   +'\n'+\
r'    unzip ${i}.zip -d $i'+'\n'+\
r'    cd $i' +'\n'+\
r'    mv input.cgyro* input.cgyro;'+'\n'+\
r'  fi'+'\n'+\
executable+r' > run_log & ' +'\n'+ \
r"  cd .." +'\n'+ \
r'done' +'\n'

bash_tail= \
r'while [[ $( ps -u '+username+r' |grep '+ps_name+r' | wc -l ) -gt 0 ]]; do ' +'\n' + \
r'  sleep 5s' +'\n'+ \
r'done ' +'\n'+ \
r'rm ${JOBID_FILE}' +'\n'



##############################################
if not os.path.exists(str(setup['workDir'])):
    os.makedirs(str(setup['workDir']))
##############################################
with open(pbs_file,'w') as f1:
    f1.write(bash_head+'\n'+bash_content+'\n'+bash_tail)
caseRoot[caseName][pbs_file.split(os.sep)[-1]]=OMFITascii(pbs_file)
inputs.append(caseRoot[caseName][pbs_file.split(os.sep)[-1]])
if setup['irun']==1:
    # sub jobs
    executable='pbsMonitor '+'scan.pbs'
    rmtserver=root['SETTINGS']['REMOTE_SETUP']['server']
    rmttunnel=root['SETTINGS']['REMOTE_SETUP']['tunnel']
    workdir=root['SETTINGS']['SETUP']['workDir']
    rmtworkdir=root['SETTINGS']['REMOTE_SETUP']['workDir']
    ret_code=OMFITx.executable(root, inputs=inputs, outputs=[],  \
                               server=rmtserver, \
                               tunnel=rmttunnel, \
                               workdir=workdir,\
                               remotedir=rmtworkdir,\
                               executable=executable,clean=True)