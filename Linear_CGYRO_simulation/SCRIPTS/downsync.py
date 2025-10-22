# this script is used to downsync those results after the calculation running is finished.
# we will have two kinds of loading method, determined by iloadmthd
setup=root['SETTINGS']['SETUP']
# iloadmthd=setup['iloadmthd'] # 0(default): OMFITcgyro, 1: OMFITgacode
# outputs={'bin.cgyro.aparb':OMFITgacode,'bin.cgyro.phib':OMFITgacode,'bin.cgyro.bparb':OMFITgacode,\
#                   'bin.cgyro.geo':OMFITgacode,'bin.cgyro.kxky_phi':OMFITgacode,'bin.cgyro.ky_flux':OMFITgacode,\
#                   'bin.cgyro.restart':OMFITgacode,'bin.cgyro.restart.old':OMFITgacode,\
#                   'input.cgyro':OMFITgacode,'input.cgyro.gen':OMFITgacode,\
#                   'out.cgyro.freq':OMFITgacode,'out.cgyro.info':OMFITgacode,'out.cgyro.time':OMFITgacode,'out.cgyro.timing':OMFITgacode,\
#                   'out.cgyro.egrid':OMFITgacode,'out.cgyro.equilibrium':OMFITgacode,'out.cgyro.grids':OMFITgacode,\
#                   'out.cgyro.hosts':OMFITgacode,'out.cgyro.memory':OMFITgacode,'out.cgyro.mpi':OMFITgacode,\
#                   'out.cgyro.prec':OMFITgacode,'out.cgyro.version':OMFITgacode,\
#                   'out.cgyro.tag':OMFITgacode,'run_log':OMFITgacode
#                  }
icgyro=setup['icgyro']
rmtsetup=root['SETTINGS']['REMOTE_SETUP']
rmtserver=rmtsetup['server']
rmtworkdir=rmtsetup['workDir']
rmttunnel=rmtsetup['tunnel']
physics=root['SETTINGS']['PHYSICS']
caseTag=physics['case_tag']
caseRoot=root['Cases']
# before going on, we may want to judge whether it's a local or remote downsync
local=False  # by defaults, the calculation is not expected to run on iris due to its limited core number
if 'scan.pbs' in caseRoot[caseTag].keys():
    del caseRoot[caseTag]['scan.pbs']
parentdir=caseRoot[caseTag].keys()
if rmtserver.split('@')[1]=='iris':
    local=True
for itemdir in parentdir:
    print("Loading "+itemdir)
    try:
        if icgyro==0:
            if local:
                caseRoot[caseTag][itemdir]=OMFITgyro(rmtworkdir+'/'+itemdir,extra_files=['RESTART_0','RESTART_tag_0','RESTART_1','RESTART_tag_1','restart.dat'])
            else:
                caseRoot[caseTag][itemdir]=OMFITgyro([rmtworkdir+'/'+itemdir,rmtserver,rmttunnel],extra_files=['RESTART_0','RESTART_tag_0','RESTART_1','RESTART_tag_1','restart.dat'])
        else:
            if local:
                caseRoot[caseTag][itemdir]=OMFITcgyro(rmtworkdir+'/'+itemdir,extra_files=['bin.cgyro.restart','bin.cgyro.restart.old'])
            else:
                caseRoot[caseTag][itemdir]=OMFITcgyro([rmtworkdir+'/'+itemdir,rmtserver,rmttunnel],extra_files=['bin.cgyro.restart','bin.cgyro.restart.old'])
    except:
        continue
#if iloadmthd==1:
for item in caseRoot[caseTag].keys():
    if 'zip' in caseRoot[caseTag].keys():
        del caseRoot[caseTag][item]['zip']
# recover the original input.cgyro
inputs_node=root['INPUTS']
if 'input.cgyro_orig' in inputs_node.keys():
    inputs_node['input.cgyro']=inputs_node['input.cgyro_orig'].duplicate()
    del inputs_node['input.cgyro_orig']
if 'input.gyro_orig' in inputs_node.keys():
    inputs_node['input.gyro']=inputs_node['input.gyro_orig'].duplicate()
    del inputs_node['input.gyro_orig']
