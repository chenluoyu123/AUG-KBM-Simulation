# try a global gyro run
for k in array([8]):
    OMFIT.load('/home/users/xiangjian/mymodule/CGYRO_SCAN/OMFITsave.txt')
    inputs=OMFIT['INPUTS']
    dirpath='/home/users/xiangjian/augwork/rho0.35/'
    inputs['input.gyro']=OMFITgacode(dirpath+'input.gyro_2ion')
    inputs['input.gacode']=OMFITgacode(dirpath+'input.gacode')
    inputgyro=inputs['input.gyro']
    inputs['input.gyro_bak']=inputgyro.duplicate()
    inputgyro_bak=inputs['input.gyro_bak']
    inputgyro['RADIAL_PROFILE_METHOD']=5 # 5 for model equilibirium with input.gyro as input
    # node and directory setup
    settings=OMFIT['SETTINGS']
    setup=settings['SETUP']
    setup['icgyro']=0
    setup['num_nodes']=1
    setup['num_cores']=24
    setup['wall_time']='24:00:00'
    setup['pbs_queue']='parallel01'
    setup['executable']='gyro -e . -n 4'
    setup['effnum']=4
    setup['idimrun']=1
    rmtsetup=settings['REMOTE_SETUP']
    rmtsetup['serverPicker']='login112'
    rmtsetup['workDir']='11'
    # prepare for linear calculation
    inputgyro['NONLINEAR_FLAG']=0
    # remove the fast ions
#    para_arr=array(['Z_3','MU_3','NI_OVER_NE_3','TI_OVER_TE_3','DLNNDR_3','DLNTDR_3'])
#    for para in para_arr:
#        del inputgyro[para]
    # some other important setting
    inputgyro['N_FIELD']=3
    inputgyro['ELECTRON_METHOD']=2
    inputgyro['BOX_MULTIPLIER']=2
    inputgyro['RADIAL_PROFILE_METHOD']=3
#    inputgyro['ORD_RBF']=5
    # resolution
    physics=OMFIT['SETTINGS']['PHYSICS']
    physics['mnp']=array([0,3,4,0])
    inputgyro['TIME_STEP']=0.02
    inputgyro['TIME_SKIP']=50
    inputgyro['TIME_MAX']=-800.001
    inputgyro['THETA_PLOT']=24
    # change the input.gyro for radial_profile_method =3
    inputgyro['NI_OVER_NE']=1
    inputgyro['Z_2']=49
    inputgyro['MU_2']=0.1
    inputgyro['NI_OVER_NE_2']=1
    inputgyro['Z_3']=1
    inputgyro['MU_3']=1
    inputgyro['NI_OVER_NE_3']=1
    inputgyro['BOUNDARY_METHOD']=2
    inputgyro['TOROIDAL_REFERENCE']=2
    inputgyro['BOX_MULTIPLIER']=2
    inputgyro['RADIAL_GRID']=128
    inputgyro['EXPLICIT_DAMP_GRID']=32
    inputgyro['RADIAL_GYRO_BAND']=16
    inputgyro['ORBIT_GRID']=8    
    inputgyro['BLEND_GRID']=8
    inputgyro['PASS_GRID']=4
    inputgyro['TRAP_GRID']=4
    inputgyro['ENERGY_GRID']=16
    inputgyro['ENERGY_MAX']=8
    # start to scan
    phy1d=physics['1d']
    gyroeigen=phy1d['gyroeigen']
    phy1d['Para']='TOROIDAL_MIN'
#    phy1d['Range']=array([3,4,5,6,8,10])     # toroidal_mean =2 corresponds to kyrhos=0.04
#    phy1d['Range']=array([3,4,5,6,8,10,12,14,16,18,20,22,24,26,28,30])
    phy1d['Range']=array([3,4,5,6,8,10,12])
#    physics['kyarr']=concatenate((linspace(0.03,0.3,10),linspace(0.4,0.8,5)))
    physics['kyarr']=array([0.12]) # this setting does not make any sense if box_multiplier is not -1`
    gyroeigen['ieigensolver']=0
    physics['mode']='lin'
    OMFIT['SCRIPTS']['runCGYRO.py'].run()
    OMFIT.saveas('/scratch/xiangjian/OMFITdata/project/gyro/initialvalue/AUG_cgyro_rho0.35_toroidalminscan_2ion_extend_highresolution.zip')