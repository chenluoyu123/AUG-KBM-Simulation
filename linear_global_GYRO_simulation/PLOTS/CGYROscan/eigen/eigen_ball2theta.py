# this is used to turn the eigenfunction from ballooning space to RZ plane
# no eikon function so that we can have an intuitive feeling about the eigenfunction distribution along theta
import numpy as np
import sys
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROscan/assist')
from cgyro_read_xj import *
from cgyro_ball import *
from cgyro_ball_class import OMFITcgyro_eigen
f = open(root['PLOTS']['CGYROscan']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
# definition of ball2theta

if setup['icgyro']==1:
    root['PLOTS']['CGYROscan']['assist']['collect.py'].run()
else:
    root['PLOTS']['CGYROscan']['assist']['collect_gyro.py'].run()

# for plot
iabs=1   # 0: real part; 1: absolute value
ncount=0
inorm=2  #normalized to be better visulized, 1: normalized to average value; 2: normalized to phi(0)
# RZ contour plot
figure(figsize=[10,8])
for dirname in root['SETTINGS']['PLOTS']['dirname']:
#        dirname='para_'+num2str_xj(paraval_item,effnum)+'~ky_'+num2str_xj(ky_item,effnum)
    case=mounttree[dirname]
    if setup['icgyro']==1:
        rmin=case['input.cgyro.gen']['RMIN']
    else:
        rmin = case['input.gyro.gen']['RADIUS']
    case.turn_ball2theta()
    phi_p0=interp(0,case.theta_p,case.phi_p)
    if case.iapar==1:
        apar_p0=interp(0,case.theta_p,case.apar_p)
    if inorm==1:
        phi_p=case.phi_p/mean(abs(case.phi_p))
        if case.iapar==1:
            apar_p = case.apar_p / mean(abs(case.apar_p))
    elif inorm==2:
        phi_p=case.phi_p/phi_p0
        if case.iapar==1:
            apar_p = case.apar_p / apar_p0
    ax=subplot(1,2,1)
    if iabs==0:
        ax.plot(case.theta_p/pi,real(case.phi_p),lab[ncount],linewidth=lw,label=dirname)
    else:
        ax.plot(case.theta_p / pi, abs(case.phi_p), lab[ncount], linewidth=lw, label=dirname)
#    print(abs(phi_theta))
    title('$\phi$',fontsize=fs2,family='serif')
    xlabel('$\\theta\pi$',fontsize=fs2,family='serif')
    if case.iapar==1:
        ax=subplot(1,2,2)
        if iabs==0:
            ax.plot(case.theta_p/pi,real(case.apar_p),lab[ncount],linewidth=lw,label=dirname)
        else:
            ax.plot(case.theta_p / pi, abs(case.apar_p), lab[ncount], linewidth=lw, label=dirname)
        title('$A_{||}$',fontsize=fs2,family='serif')
        xlabel('$\\theta\pi$',fontsize=fs2,family='serif')
    ncount=ncount+1
xticks(fontsize=fs2,family='serif')
yticks(fontsize=fs2,family='serif')
legend(loc=0).draggable(True)

for item in mounttree.keys():
    del mounttree[item]
