# this script is used to plot the zonal flow shearing rate versus theta
import sys
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROalone/assist')
f = open(root['PLOTS']['CGYROalone']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
icgyro=root['SETTINGS']['SETUP']['icgyro']
if icgyro==1:
    root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
else:
    root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()  # currently the gyro part is NOT supported
B_phi0=1 # used for normlization of all the other magnetic related quantities
rhos_over_a=5.e-3 # subject to change depending on the concrete experimental conditions
box_size=40       # box_size * rhos_over_a determines the radial range of interests
num_grid=128
kx_pvt=0.06
ishape=1
isign=-1 # might be 1 as well
# plot the zonal field
figure(figsize=[18,12])
for k_case in range(n_case):
    casename=case_plot[k_case]
    casek = outputs[casename]
# construct the equilibrium profile
    B_theta_0=casek.rmin*B_phi0/casek.Rmaj/casek.q
    rmin_0=casek.rmin
    rmin=linspace(rmin_0-box_size*rhos_over_a/2, rmin_0+box_size*rhos_over_a/2, num_grid)
    rmin_over_Btheta=casek.shear/B_theta_0*(rmin-rmin_0)+casek.rmin/B_theta_0
    B_theta=rmin/rmin_over_Btheta  # equilibrium B_theta profile
#  plot the profiles
    subplot(3,2,1)
    plot(rmin,B_theta, linewidth=lw,label='Eq'+casename)
    legend(loc=0,fontsize=fs2).draggable(True)
    ylabel('$B_{\\theta}/B_{\phi}$',fontsize=fs1,family='serif')
    # print(casek.freq_n[i_n].data)
    xticks([],fontsize=fs1)
    xlim(array([min(rmin),max(rmin)]))
    yticks(fontsize=fs1)
    subplot(3,2,2)
    casek.get_zonal_field(theta=0,field=1,imthd=1,kx_pvt=kx_pvt,ishape=ishape)
    print('energy_diff=',casek.energy_diff)
    f_cubic=interp1d(casek.x*rhos_over_a+rmin_0,real(casek.B_theta_x),kind='cubic')
    B_theta_Z=f_cubic(rmin)
    plot(rmin, B_theta_Z,lab[k_case],linewidth=lw,label=casename)
    xticks([],fontsize=fs1)
    yticks(fontsize=fs1)
    subplot(3,2,3)
    q=rmin*B_phi0/casek.Rmaj/B_theta # equilibrium q
    q_tot=rmin*B_phi0/casek.Rmaj/(B_theta+B_theta_Z)
    plot(rmin,q,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,q_tot,labo[k_case],linewidth=lw,label='Eq+Zonal_'+casename)
    ylabel('$q$', fontsize=fs1, family='serif')
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    legend(loc=0,fontsize=fs2).draggable(True)
    subplot(3,2,4)
    plot(rmin,rmin/q*np.gradient(q)/np.gradient(rmin),lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,rmin/q_tot*np.gradient(q_tot)/np.gradient(rmin),labo[k_case],linewidth=lw,label='Eq+Zonal_'+casename)
    ylabel('$s$',fontsize=fs1)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    subplot(3,2,5)
    plot(rmin,1./rmin*gradient(rmin*B_theta)/gradient(rmin),linewidth=lw,label='Eq_'+casename)
    plot(rmin,1./rmin*gradient(rmin*(B_theta+B_theta_Z))/gradient(rmin),linewidth=lw,label='Eq+Zonal_'+casename)
    legend(loc=0,fontsize=fs2).draggable(True)
    xlim(array([min(rmin),max(rmin)]))
    xlabel('$r_{min}/a$',fontsize=fs2)
    ylabel('$J_{\phi}(.a.u)$',fontsize=fs2)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    subplot(3,2,6)
    casek.get_zonal_field(theta=0,field=0,imthd=1,kx_pvt=kx_pvt,ishape=ishape)
    print(casek.energy_diff)
    f_cubic=interp1d(casek.x*rhos_over_a+rmin_0,casek.field_x,kind='cubic')
    phi_Z_x=f_cubic(rmin)
    plot(rmin, abs(phi_Z_x)**2,linewidth=lw,label=casename)
    legend(loc=0,fontsize=fs2).draggable(True)
    xlim(array([min(rmin),max(rmin)]))
    xlabel('$r_{min}/a$',fontsize=fs2)
    ylabel('$\phi_{ZF}^2$',fontsize=fs2)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)

# plot the zonal profiles which can corrugate the equilibrium profiles
ion_species=0
figure(figsize=[24,12])
for k_case in range(n_case):
    casename=case_plot[k_case]
    casek = outputs[casename]
# construct the equilibrim profiles
    rmin_0=casek.rmin
    rmin=linspace(rmin_0-box_size*rhos_over_a/2, rmin_0+box_size*rhos_over_a/2, num_grid)
    inputcgyro=casek['input.cgyro.gen']
    n_species=inputcgyro['N_SPECIES']
    ni_over_ne=inputcgyro['DENS_'+str(ion_species+1)]
    Ti_over_Te=inputcgyro['TEMP_'+str(ion_species+1)]
    aLni=casek['input.cgyro.gen']['DLNNDR_'+str(ion_species+1)]
    aLTi=casek['input.cgyro.gen']['DLNTDR_'+str(ion_species+1)]
    aLne=casek['input.cgyro.gen']['DLNNDR_'+str(n_species)]
    aLTe=casek['input.cgyro.gen']['DLNTDR_'+str(n_species)]
# equilibrium n and P
    ni=(rmin-rmin_0)*ni_over_ne*aLni*(-1)+ni_over_ne
    ne=(rmin-rmin_0)*1*aLne*(-1)+1
    Pi=(rmin-rmin_0)*ni_over_ne*Ti_over_Te*(aLni+aLTi)*(-1)+ni_over_ne*Ti_over_Te
    Pe=(rmin-rmin_0)*1*1*(aLne+aLTe)*(-1)+1
    subplot(4,3,1)
    casek.get_zonal_pro(i_species=ion_species, i_moment=0,imthd=1,kx_pvt=kx_pvt,ishape=ishape) # 0 for n and 1 for e
    print(casek.energy_diff)
    f_cubic=interp1d(casek.x*rhos_over_a+rmin_0,real(casek.pro_x),kind='cubic')
    pro_x_ni=isign*f_cubic(rmin)
    plot(rmin,pro_x_ni,lab[k_case],linewidth=lw,label='Eq_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$n_{i,Z}$',fontsize=fs1)
    subplot(4,3,2)
    plot(rmin,ni,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,ni+pro_x_ni,labo[k_case],linewidth=lw,label='Eq+zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$n_i$',fontsize=fs1)
    subplot(4,3,3)
    plot(rmin,-1*gradient(ni)/gradient(rmin)/ni,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,-1*gradient(ni+pro_x_ni)/gradient(rmin)/(ni+pro_x_ni),labo[k_case],linewidth=lw,label='Eq+zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$a/L_{ni}$',fontsize=fs1)
    subplot(4,3,4)
    casek.get_zonal_pro(i_species=-1, i_moment=0,imthd=1,kx_pvt=kx_pvt,ishape=ishape) # 0 for n and 1 for e
    print(casek.energy_diff)
    f_cubic=interp1d(casek.x*rhos_over_a+rmin_0,real(casek.pro_x),kind='cubic')
    pro_x_ne=isign*f_cubic(rmin)
    plot(rmin,pro_x_ne,lab[k_case],linewidth=lw,label='zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$n_{e,Z}$',fontsize=fs1)
    subplot(4,3,5)
    plot(rmin,ne,lab[k_case],linewidth=lw,label='zonal_'+casename)
    plot(rmin,ne+pro_x_ne,labo[k_case],linewidth=lw,label='zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$n_e$',fontsize=fs1)
    subplot(4,3,6)
    plot(rmin,-1*gradient(ne)/gradient(rmin)/ne,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,-1*gradient(ni+pro_x_ne)/gradient(rmin)/(ni+pro_x_ne),labo[k_case],linewidth=lw,label='Eq+zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$a/L_{ne}$',fontsize=fs1)
    subplot(4,3,7)
    casek.get_zonal_pro(i_species=ion_species, i_moment=1,imthd=1,kx_pvt=kx_pvt,ishape=ishape) # 0 for n and 1 for e
    print(casek.energy_diff)
    f_cubic=interp1d(casek.x*rhos_over_a+rmin_0,real(casek.pro_x),kind='cubic')
    pro_x_pi=isign*f_cubic(rmin)
    plot(rmin,pro_x_pi,lab[k_case],linewidth=lw,label='zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$P_{i,Z}$',fontsize=fs1)
    subplot(4,3,8)
    plot(rmin,Pi,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,Pi+pro_x_pi,labo[k_case],linewidth=lw,label='Eq+zonal_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$Pi$',fontsize=fs1)
    subplot(4,3,9)
    plot(rmin,-1*gradient(Pi)/gradient(rmin)/Pi,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,-1*gradient(Pi+pro_x_pi)/gradient(rmin)/(Pi+pro_x_pi),labo[k_case],linewidth=lw,label='Eq_'+casename)
    xlim(array([min(rmin),max(rmin)]))
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$a/L_{Pi}$',fontsize=fs1)
    subplot(4,3,10)
    casek.get_zonal_pro(i_species=-1, i_moment=1,imthd=1,kx_pvt=kx_pvt,ishape=ishape) # 0 for n and 1 for e
    print(casek.energy_diff)
    f_cubic=interp1d(casek.x*rhos_over_a+rmin_0,real(casek.pro_x),kind='cubic')
    pro_x_pe=isign*f_cubic(rmin)
    plot(rmin,pro_x_pe,lab[k_case],linewidth=lw,label='zonal_'+casename)
    legend(loc=0,fontsize=fs2).draggable(True)
    xlim(array([min(rmin),max(rmin)]))
    xlabel('$r_{min}/a$',fontsize=fs2)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$P_{e,Z}$',fontsize=fs1)
    subplot(4,3,11)
    plot(rmin,Pe,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,Pe+pro_x_pe,labo[k_case],linewidth=lw,label='Eq+zonal_'+casename)
    legend(loc=0,fontsize=fs2).draggable(True)
    xlim(array([min(rmin),max(rmin)]))
    xlabel('$r_{min}/a$',fontsize=fs2)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$P_{e}$',fontsize=fs1)
    subplot(4,3,12)
    plot(rmin,-1*gradient(Pe)/gradient(rmin)/Pe,lab[k_case],linewidth=lw,label='Eq_'+casename)
    plot(rmin,-1*gradient(Pe+pro_x_pe)/gradient(rmin)/(Pe+pro_x_pe),labo[k_case],linewidth=lw,label='Eq+zonal_'+casename)
    legend(loc=0,fontsize=fs2).draggable(True)
    xlim(array([min(rmin),max(rmin)]))
    xlabel('$r_{min}/a$',fontsize=fs2)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    ylabel('$a/L_{Pe}$',fontsize=fs1)