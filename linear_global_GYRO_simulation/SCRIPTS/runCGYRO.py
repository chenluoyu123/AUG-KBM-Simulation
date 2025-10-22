# this script is used to do the overall call for the 1D and 2D scan
rmt_setup=root['SETTINGS']['REMOTE_SETUP']
setup=root['SETTINGS']['SETUP']
server_list=OMFIT['MainSettings']['SERVER']
icgyro=root['SETTINGS']['SETUP']['icgyro']
user_name={'iris':server_list['GA_USERNAME'],\
           'shenma':server_list['ASIPP_username'],\
           'kuafu':server_list['ASIPP_username']}
if rmt_setup['serverPicker']=='iris':
    rmt_setup['server']=user_name['iris']+'@iris'
    rmt_setup['tunnel']=''
#    setup['pbs_queue']='medium'
    if len(rmt_setup['workDir'])<6:
        rmt_setup['workDir']='/cluster-scratch/'+user_name['iris']+'/OMFIT/runs/CGYRO_Scan/'+rmt_setup['workDir']
elif rmt_setup['serverPicker']=='omega':
    rmt_setup['server']='jianx@omega'
    rmt_setup['tunnel']='jianx@cybele.gat.com:2039'
    if len(rmt_setup['workDir'])<6:
        rmt_setup['workDir']='/cscratch/jianx/OMFIT/runs/CGYRO_Scan/'+rmt_setup['workDir']
elif rmt_setup['serverPicker'][0:5]=='login':
#    rmt_setup['server']=user_name['shenma']+'@shenma.ipp.ac.cn:22'
    rmt_setup['server']=user_name['shenma']+'@'+rmt_setup['serverPicker']
    rmt_setup['tunnel']=''
    if rmt_setup['serverPicker']=='login110':
        setup['pbs_queue']='parallel21'
#    else:
#        setup['pbs_queue']='parallel01'
    if len(rmt_setup['workDir']) < 6:
        if icgyro==1:
            rmt_setup['workDir']='/scratch/'+user_name['shenma']+'/OMFIT/runs/CGYRO_Scan/'+rmt_setup['workDir']
        else:
            rmt_setup['workDir']='/scratch/'+user_name['shenma']+'/OMFIT/runs/GYRO_Scan/'+rmt_setup['workDir']
else:
    rmt_setup['server']=user_name['kuafu']+'@202.141.160.120:15522'
    rmt_setup['tunnel']=user_name['iris']+'@cybele.gat.com:2039'
    setup['pbs_queue']='batch'
    if len(rmt_setup['workDir']) < 6:
        rmt_setup['workDir']='/scratch/'+user_name['kuafu']+'/OMFIT/runs/CGYRO_Scan/'+rmt_setup['workDir']
idimrun=setup['idimrun']
if idimrun==1:
    root['SCRIPTS']['CGYROScan.py'].run()
else:
    root['SCRIPTS']['CGYROScan_2d.py'].run()
