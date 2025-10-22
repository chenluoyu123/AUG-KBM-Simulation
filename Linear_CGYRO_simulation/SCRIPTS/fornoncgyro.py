#dirpath='/global/cfs/cdirs/atom/users/xiang/highgradient/ESRun/highntheta/box4/dky0.2/box8/box12'
#OMFIT['scratch']['input.cgyro']=OMFITgacode([dirpath+'/input.cgyro','xiang@cori.nersc.gov','jianx@cybele.gat.com:2039'])
#inputcgyro=OMFIT['scratch']['input.cgyro']
inputcgyro=root['INPUTS']['input.cgyro']
#root['INPUTS']['input.cgyro_bak']=inputcgyro.duplicate()
# addjust line
#inputcgyro['BOX_SIZE']=1
# calculate
q=inputcgyro['Q']
s=inputcgyro['S']
ky=inputcgyro['KY']
box_size=inputcgyro['BOX_SIZE']
rmin=inputcgyro['RMIN']
# get the parameters for kx, donated by _x
n_x=inputcgyro['N_RADIAL']
delta_x=2*pi*ky*s/box_size
max_x=(n_x/2.-1)*delta_x
Lx=box_size/ky/s
# get the parameters for ky, donated by _y
n_y=inputcgyro['N_TOROIDAL']
delta_y=ky
max_y=(n_y-1)*ky
Ly=2*pi/ky
# print
print ('          n     Delta    Max    L/rho')
print ('kx : %7.0f %7.3f %7.3f %7.1f'%(n_x,delta_x,max_x,Lx))
print ('ky : %7.0f %7.3f %7.3f %7.1f'%(n_y,delta_y,max_y,Ly))
print('dx/rho_s=%7.2f'%(Lx/n_x))
# we will also have some suggestions for the number of cores
nv=inputcgyro['N_ENERGY']*inputcgyro['N_XI']*inputcgyro['N_SPECIES']
nc=inputcgyro['N_RADIAL']*inputcgyro['N_THETA']
#n_proc1=n_proc/inputcgyro['N_TOROIDAL']
# common divisor
print('you can choose the node number,total cores to be:')
ncore_per_node=64
for k in arange(1,min([nv,nc])):
    if ((nv%k==0) and (nc%k==0)):
#        print('common diviror=')
#        print(k)
        if (k*inputcgyro['N_TOROIDAL'])%ncore_per_node==0:
            print(k*inputcgyro['N_TOROIDAL']/ncore_per_node,k*inputcgyro['N_TOROIDAL'])
# add some nonlinear run cases
inputcgyro['NONLINEAR_FLAG']=1
inputcgyro['THETA_PLOT']=inputcgyro['N_THETA']/2
inputcgyro['FIELD_PRINT_FLAG']=1
inputcgyro['MOMENT_PRINT_FLAG']=1
inputcgyro['KXKYFLUX_PRINT_FLAG']=1
