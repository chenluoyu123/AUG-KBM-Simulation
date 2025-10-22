import matplotlib.pyplot as plt
n_arr_gyro=array([4,5,6,8,10,12,14,16,18,20,22,24,26,28,30])
omega_gyro=array([  -8.26826875e-01, -9.17233333e-01, -1.19134286e+00,\
  -1.06360000e+00, -1.22252500e+00, -1.69310000e+00, -1.41940000e+00,\
  -1.88583333e+00, -1.56195000e+00,  4.48075000e+01,  4.98460000e+01,\
   5.49160000e+01,  5.99450000e+01,  6.48730000e+01,  6.96660000e+01])
gamma_gyro=array([ 0.09426975, 0.092292  , 0.06826571 ,0.06398544 ,0.05971256,\
  0.06142975, 0.07930933 ,0.07603133, 0.084169  , 0.201915 ,  1.3311,\
  2.7601     ,4.5182  ,   6.5653   ,  8.9278    ])
n_arr_nlt=array([4,5,6,7,8,10])
omega_nlt=array([-848757,-933051,-761556,-880731,-984234,-1047459])
gamma_nlt=array([95784,107854,37801,84118,117845,55390])
kyarr_cgyro=concatenate((linspace(0.03,0.3,10),linspace(0.4,0.8,4)))
omega_cgyro=array( [-0.256915 ,  -0.46844   , -0.61442 ,   -0.73648 ,   -0.849717 ,  -0.939657,\
  -1.0607 ,    -1.14795  ,  -0.166069  , -0.187176  , -0.24213 ,   -0.30098,\
  -0.358736 ,  -0.420711  ])
gamma_cgyro=array( [0.0484965,  0.157911 ,  0.19624  ,  0.20193   , 0.183617  , 0.145541,\
  0.10172   , 0.0553525 , 0.0456209 , 0.0437758,  0.064151,   0.0666611,\
  0.04696666, 0.0200055 ])
csovera=144000*2*pi
ind_end_cgyro=12
n_cgyro=linspace(3,12,10)
f_linear=interp1d(kyarr_cgyro[0:ind_end_cgyro]*50,omega_cgyro[0:ind_end_cgyro],'cubic')
omega_cgyro_intp=f_linear(n_cgyro)
f_linear=interp1d(kyarr_cgyro[0:ind_end_cgyro]*50,gamma_cgyro[0:ind_end_cgyro],'cubic')
gamma_cgyro_intp=f_linear(n_cgyro)
figure(figsize=[6,8])
lw=3
fs2=24
ind_end_gyro=6
ind_end_cgyro=8
subplot(211)
plot(n_arr_gyro[0:ind_end_gyro],omega_gyro[0:ind_end_gyro],'-bo',linewidth=lw,label='GYRO-Global')
plot(n_arr_nlt,omega_nlt/csovera,'-ko',linewidth=lw,label='NLT-Global')
#plot(kyarr_cgyro[0:ind_end_cgyro]*50,omega_cgyro[0:ind_end_cgyro],'-ro',linewidth=lw,label='CGYRO-Local($\\rho$=0.35)')
plot(n_cgyro,omega_cgyro_intp,'-ro',linewidth=lw,label='CGYRO-Local($\\rho$=0.35)')
xticks([],fontsize=fs2)
yticks(linspace(-2,0,3),fontsize=fs2)
xlim([2,12])
ylim([-2.2,0])
text(3,-0.2,'(e)',fontsize=fs2-4)
text(6,-0.2,'$\\omega(c_s/a)$',fontsize=fs2)
#title('$\\omega(c_s/a)$',fontsize=fs2)
legend(loc=0,fontsize=fs2-4).draggable(True)
subplot(212)
plot(n_arr_gyro[0:ind_end_gyro],gamma_gyro[0:ind_end_gyro],'-bo',linewidth=lw)
plot(n_arr_nlt,gamma_nlt/csovera,'-ko',linewidth=lw)
#plot(kyarr_cgyro[0:ind_end_cgyro]*50,gamma_cgyro[0:ind_end_cgyro],'-ro',linewidth=lw)
plot(n_cgyro,gamma_cgyro_intp,'-ro',linewidth=lw)
xlabel('n',fontsize=fs2)
xlim([2,12])
xticks(linspace(3,12,4),fontsize=fs2)
yticks(linspace(0,0.2,3),fontsize=fs2)
ylim([0,0.25])
#title('$\\gamma(c_s/a)$',fontsize=fs2)
text(3,0.22,'(f)',fontsize=fs2-4)
text(6,0.22,'$\\gamma(c_s/a)$',fontsize=fs2)
plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.1, wspace=0.4, hspace=0.05)