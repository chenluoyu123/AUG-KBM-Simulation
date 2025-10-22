plots=root['SETTINGS']['PLOTS']
plt1d=plots['1d']  # type: object
plt2d=plots['2d']
physics=root['SETTINGS']['PHYSICS']
phy1d=physics['1d']
phy2d=physics['2d']
root['PLOTS']['CGYROscan']['assist']['starter.py'].run()
setup = root['SETTINGS']['SETUP']
Para = plt1d['Para']
Para_x = plt2d['Para_x']
Para_y = plt2d['Para_y']
Range=plt1d['Range']*plots['Range_scale']
Range_x=plt2d['Range_x']
Range_y=plt2d['Range_y']
Range_x = Range_x * plots['Range_xscale']
Range_y = Range_y * plots['Range_yscale']
nRange=len(Range)
nRange_x=len(Range_x)
nRange_y=len(Range_y)
kyarr = plots['kyarr']
num_ky=len(kyarr)
para_eigen=plt1d['para_eigen']
para_x_eigen=plt2d['para_x_eigen']
para_y_eigen=plt2d['para_y_eigen']
ky_eigen=plots['ky_eigen']
#setup
effnum = setup['effnum']
iloadmthd= setup['iloadmthd']
inputcgyro = root['INPUTS']['input.cgyro']
ns=inputcgyro['N_SPECIES']
# for plotting
ilogx=plots['ilogx']
ilogy=plots['ilogy']
scale_dic={0:'linear',1:'log'}
idimplt=plots['idimplt']
iflwphy=plots['iflwphy']
ibelow0=plots['ibelow0']
ilim_err=plots['ilim_err']
Paraname_dic={'DLNTDR_1':'$a/L_{Ti}$', 'DLNTDR_'+str(ns):'$a/L_{Te}$','BETA_STAR_SCALE':'$\\alpha_{scale}$','BETAE_UNIT':'$\\beta_{E}$','GAMMA_E':'$\gamma_{E}$','S':'s','Q':'q'}
mounttree=root['OUTPUTS']
#deploypath='/home/task1/xiangjian/temp/'
deploypath='/home/jianx/temp/'
# specifically for 1d
powky=plt1d['powky']
ind_kyfit=plt1d['ind_kyfit']
ifit_mthd=plt1d['ifit_mthd']
ipltExB=plt1d['ipltExB']
# specifically for 2d
Channel_name=['Pe','Pi','Ge','Gi','Qe','Qi']
Para_x_dic={0:Para_x, 1: Para_y}
Para_y_dic={0:Para_y, 1: Para_x}
##
ms=12
lw=2
fs1=24
fs2=20
fs3=16
lab=['-kd','-md','-bd','-rd','-cd','-gd','-yd','-k' ,'-m', '-b', '-r', '-c', '-g', '-y','-ko','-mo','-bo','-ro','-co','-go','-yo','-k*','-m*','-b*','-r*','-c*','-g*','-y*','-ks','-ms','-bs','-rs','-cs','-gs','-ys','--k','--m','--b','--r','--c','--g','--y' ];

