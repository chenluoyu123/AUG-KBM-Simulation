# this script is used to collect the CGYRO result of 2D scan in a good format
setup=root['SETTINGS']['SETUP']
# start to scan
# run
if setup['irun']==1:
    if setup['icgyro']==1:
        root['SCRIPTS']['subscan_lin_2d.py'].run()
    else:
        root['SCRIPTS']['set_resolution.py'].run()
        root['SCRIPTS']['subscan_lin_2d_gyro.py'].run()
# load the results from the remote sever
if setup['idownsync']==1:
    root['SCRIPTS']['downsync.py'].run()
case_tag=root['SETTINGS']['PHYSICS']['case_tag']
if 'scan.pbs' in root['Cases'][case_tag].keys():
    del root['Cases'][case_tag]['scan.pbs']
for item in root['Cases'][case_tag].keys():
#        print(item)
    item_temp=item.split('~')
    if not  item_temp[0] in root['OUTPUTScan'].keys():
        root['OUTPUTScan'][item_temp[0]]=OMFITtree()
    if not  item_temp[1] in root['OUTPUTScan'][item_temp[0]].keys():
        root['OUTPUTScan'][item_temp[0]][item_temp[1]]=OMFITtree()
    if not item_temp[2] in root['OUTPUTScan'][item_temp[0]][item_temp[1]].keys():
        root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]]=OMFITtree()
    if not  item_temp[3] in root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]].keys():
        root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]][item_temp[3]]=OMFITtree()
    if not 'lin' in root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]][item_temp[3]].keys():
        root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]][item_temp[3]]['lin']=OMFITtree()
    if isinstance(root['Cases'][case_tag][item],OMFITtree):
        root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]][item_temp[3]]['lin'][item_temp[5]]=OMFITtree()
        for files in root['Cases'][case_tag][item]:
            root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]][item_temp[3]]['lin'][item_temp[5]][files]=root['Cases'][case_tag][item][files].duplicate()
    else:
        root['OUTPUTScan'][item_temp[0]][item_temp[1]][item_temp[2]][item_temp[3]]['lin'][item_temp[5]]=copy.deepcopy(root['Cases'][case_tag][item])
