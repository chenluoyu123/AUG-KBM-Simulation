# this script is used to scan the Para of CGYRO using initial value method
setup=root['SETTINGS']['SETUP']
icgyro=setup['icgyro']
ieigensolver=root['SETTINGS']['PHYSICS']['1d']['gyroeigen']['ieigensolver']
# run
if setup['irun']==1:
    if icgyro==1:
        root['SCRIPTS']['subscan_lin.py'].run()
    else:
        root['SCRIPTS']['set_resolution.py'].run()
        if ieigensolver == 0:
            root['SCRIPTS']['subscan_lin_gyro.py'].run()
        else:
            print('GYRO eigensolver Comes!')
            setup['idownsync']=0
            root['SCRIPTS']['LM3Scan.py'].run()
if setup['idownsync']==1:
    root['SCRIPTS']['downsync.py'].run()
# rearange those results
case_tag=root['SETTINGS']['PHYSICS']['case_tag']
if 'scan.pbs' in root['Cases'][case_tag].keys():
    del root['Cases'][case_tag]['scan.pbs']
for item in root['Cases'][case_tag].keys():
    item_temp=item.split('~')
    if not  item_temp[0] in root['OUTPUTScan'].keys():
        root['OUTPUTScan'][item_temp[0]]=OMFITtree()
    if not  item_temp[1] in root['OUTPUTScan'][item_temp[0]].keys():
        root['OUTPUTScan'][item_temp[0]][item_temp[1]]=OMFITtree()
    if not 'lin' in root['OUTPUTScan'][item_temp[0]][item_temp[1]].keys():
        root['OUTPUTScan'][item_temp[0]][item_temp[1]]['lin']=OMFITtree()
    if isinstance(root['Cases'][case_tag][item],OMFITtree):
        root['OUTPUTScan'][item_temp[0]][item_temp[1]]['lin'][item_temp[3]]=OMFITtree()
        for files in root['Cases'][case_tag][item]:
            root['OUTPUTScan'][item_temp[0]][item_temp[1]]['lin'][item_temp[3]][files]=root['Cases'][case_tag][item][files]
    else:
        root['OUTPUTScan'][item_temp[0]][item_temp[1]]['lin'][item_temp[3]]=copy.deepcopy(root['Cases'][case_tag][item])
## delete those results which has already been aranged
# this is not required since we may want the restart mode
#        del root['Cases'][case_tag][item]
