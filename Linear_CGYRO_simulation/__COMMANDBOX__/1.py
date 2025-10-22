# AUG, rho=0.35, fast ion on
# nfield=3
# refined scan of ky 
# scan of fast ion fraction
for k in array([9]):
    OMFIT.load('/home/users/xiangjian/mymodule/CGYRO_SCAN/OMFITsave.txt')
    dirpath='/home/users/xiangjian/augwork/rho0.35/'
    inputs=OMFIT['INPUTS']
    inputs['input.cgyro']=OMFITgacode(dirpath+'input.cgyro_2ion')
    inputcgyro=inputs['input.cgyro']
    inputs['input.cgyro_bak']=inputcgyro.duplicate()
    inputcgyro_bak=inputs['input.cgyro_bak']
    inputcgyro['PROFILE_MODEL']=1
    inputcgyro['EQUILIBRIUM_MODEL']=2
    # node and directory setup
    settings=OMFIT['SETTINGS']
    setup=settings['SETUP']
    setup['num_nodes']=1
    setup['num_cores']=24
    setup['wall_time']='96:00:00'
    setup['pbs_queue']='parallel01'
    setup['executable']='cgyro -e . -n 4'
    setup['effnum']=5
    setup['idimrun']=1
    rmtsetup=settings['REMOTE_SETUP']
    rmtsetup['serverPicker']='login112'
    rmtsetup['workDir']='10'
    # prepare for linear calculation
    inputcgyro['NONLINEAR_FLAG']=0
    inputcgyro['GAMMA_E_SCALE']=0.
    inputcgyro['N_TOROIDAL']=1
    # resolution
    inputcgyro['DELTA_T']=0.05
    inputcgyro['PRINT_STEP']=20
    inputcgyro['DELTA_T_METHOD']=1
    inputcgyro['ERROR_TOL']=1.e-5
    inputcgyro['MAX_TIME']=500.
    inputcgyro['N_RADIAL']=6
    inputcgyro['N_THETA']=32
    inputcgyro['THETA_PLOT']=1
    inputcgyro['FREQ_TOL']=1.e-3
    # some other important setting
    inputcgyro['N_FIELD']=3
    # start to scan
    physics=OMFIT['SETTINGS']['PHYSICS']
    scale_fac=linspace(0,1,6)
    phy1d=physics['1d']
    phy1d['Para']='DENS_2'
    phy1d['Range']=inputcgyro['DENS_2']*array([0.001,0.3,0.7,1])
    physics['kyarr']=concatenate((linspace(0.03,0.3,10),linspace(0.4,0.8,4)))
    physics['mode']='lin'
    OMFIT['SCRIPTS']['runCGYRO.py'].run()
    OMFIT.saveas('/scratch/xiangjian/OMFITdata/project/cgyro/AUG_cgyro_rho0.35_nfscan_2ion_refine.zip')