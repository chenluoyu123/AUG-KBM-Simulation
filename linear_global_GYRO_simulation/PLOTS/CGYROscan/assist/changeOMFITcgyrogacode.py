# this script is used to change the format in the caseRoot in between of OMFITcgyro and OMFITgacode
import os
CaseRoot=root['Cases']
case_tag=root['SETTINGS']['PHYSICS']['case_tag']
setup=root['SETTINGS']['SETUP']
deploypath='/home/jianx/temp/'
for item in CaseRoot[case_tag]:
    if isinstance(CaseRoot[case_tag][item],OMFITcgyro):
        try:
            CaseRoot[case_tag][item].deploy(deploypath+item)
            CaseRoot[case_tag][item]=OMFITtree()
            for (roott, dirs, filenames) in os.walk(deploypath+item):
                for filename in filenames:
                    CaseRoot[case_tag][item][filename]=OMFITgacode(roott+'/'+filename)
        except:
            CaseRoot[case_tag][item]=OMFITtree()
    else:
         CaseRoot[case_tag][item].deploy(deploypath+item)
         CaseRoot[case_tag][item]=OMFITcgyro(deploypath+item)

# turn from Caseroot to OUTPUT
irun=setup['irun']
idownsync=setup['idownsync']
setup['irun']=0
setup['idownsync']=0
if len(item.split('~'))<5:
    print('Redistribute the 1D scan data!')
    root['SCRIPTS']['CGYROScan.py'].run()
else:
    print('Redistribute the 2D scan data!')
    root['SCRIPTS']['CGYROScan_2d.py'].run()
# recover
setup['irun']=irun
setup['idownsync']=idownsync
