# for plotting the energy transfer to a selected (kx_select, ky_select)
import sys
from matplotlib import ticker, cm
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROalone/assist')
f = open(root['PLOTS']['CGYROalone']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
# check to see how many fields that all the cases have in common
n_field=1
for k_case in range(n_case):
    casek = outputs[case_plot[k_case]]
    n_field=min(len(casek['field_tags']),n_field)
nlevel=11
kx_select=0.0
ky_select=20
# note that it could be better if we split the time window in to several sub time windows
plot_nl=root['SETTINGS']['PLOTS']['nl']
t_ave_orig=plot_nl['t_ave']
t_end_orig=plot_nl['t_end']
n_split=1
# # plot Energy transfer of phi (kinetic energy)
figure(figsize=[15,8])
for k_case in range(n_case):
    casek = outputs[case_plot[k_case]]
    t_ave = t_ave_orig / n_split
    print(case_plot[k_case])
    i_theta_plot=casek.theta_plot//2
    S_k_kp_split = zeros([n_split, casek.n_r, casek.n_n])
    S_k_kp_norm_split = zeros([n_split, casek.n_r, casek.n_n])
    T_phi_split = zeros([n_split, casek.n_r, casek.n_n])
    for i_split in arange(n_split):
        t_end=t_end_orig-t_ave*i_split
        plot_nl['t_ave']=t_ave
        plot_nl['t_end'] = t_end
        print(t_ave,t_end)
        # casek = outputs[case_plot[k_case]]  # datanode
        root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
        casek.Energy_transfer_phi(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
        S_k_kp_split[i_split]=casek.S_k_kp_phi
        S_k_kp_norm_split[i_split] = casek.S_k_kp_norm_phi
        T_phi_split[i_split] = casek.T_phi
    S_k_kp_ave=mean(S_k_kp_split,axis=0)
    S_k_kp_norm_ave = mean(S_k_kp_norm_split, axis=0)
    T_phi_ave = mean(T_phi_split, axis=0)
    kx_grid, ky_grid = meshgrid(casek.kx, casek.ky)
    axk=subplot(3,n_case,k_case+1)
    Lamda_abs_max=amax(abs(casek.Lamda_phi))
    contourf(kx_grid, ky_grid, casek.Lamda_phi.T,levels=linspace(-1*Lamda_abs_max,Lamda_abs_max,nlevel),cmap='seismic')
    colorbar()
    xticks([])
    if k_case == 0:
        ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    title('Coupling Coefficient # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.2f $k_y\\rho_s$=%.2f' % (casek.ftheta_plot[i_theta_plot],kx_select,ky_select) , fontsize=fs2, family='serif')

    axk=subplot(3,n_case,n_case+k_case+1)
    # S_k_kp_absmax=amax(abs(casek.S_k_kp))
    # contourf(kx_grid,ky_grid,casek.S_k_kp.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    S_k_kp_absmax=amax(abs(S_k_kp_ave))
    contourf(kx_grid,ky_grid,S_k_kp_ave.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    colorbar()
    xticks([])
    if k_case==0:
        ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    title('Bicoherence', fontsize=fs2, family='serif')

    # axk=subplot(4,n_case,2*n_case+k_case+1)
    # S_k_kp_absmax=amax(abs(casek.S_k_kp_norm))
    ## contourf(kx_grid,ky_grid,casek.S_k_kp_norm.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    # contourf(kx_grid,ky_grid,S_k_kp_norm_ave.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    # colorbar()
    # if k_case==0:
    #     ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    # title('Bicoherence(norm)', fontsize=fs2, family='serif')

    axk=subplot(3,n_case,2*n_case+k_case+1)
    # T_phi_absmax=amax(abs(casek.T_phi))
    # contourf(kx_grid,ky_grid,casek.T_phi.T,levels=linspace(-1*T_phi_absmax,T_phi_absmax,nlevel),cmap='seismic')
    T_phi_absmax=amax(abs(T_phi_ave))
    contourf(kx_grid,ky_grid,T_phi_ave.T,levels=linspace(-1*T_phi_absmax,T_phi_absmax,nlevel),cmap='seismic')
    colorbar()
    if k_case==0:
        ylabel('$k_y\\rho_s$',fontsize=fs2,family='serif')
    title('$T_{\phi}$',fontsize=fs2,family='serif')
    xlabel('$k_x\\rho_s$',fontsize=fs2,family='serif')

# plot Energy transfer of p (internal energy)
figure(figsize=[15,8])
for k_case in range(n_case):
    t_ave = t_ave_orig / n_split
    print(case_plot[k_case])
    i_theta_plot=casek.theta_plot//2
    S_k_kp_split = zeros([n_split, casek.n_r, casek.n_n])
    S_k_kp_norm_split = zeros([n_split, casek.n_r, casek.n_n])
    T_p_split = zeros([n_split, casek.n_r, casek.n_n])
    for i_split in arange(n_split):
        t_end=t_end_orig-t_ave*i_split
        plot_nl['t_ave']=t_ave
        plot_nl['t_end'] = t_end
        casek = outputs[case_plot[k_case]]
        root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
        casek.Energy_transfer_p(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
        S_k_kp_split[i_split]=casek.S_k_kp_p
        S_k_kp_norm_split[i_split] = casek.S_k_kp_norm_p
        T_p_split[i_split] = casek.T_p
    S_k_kp_ave=mean(S_k_kp_split,axis=0)
    S_k_kp_norm_ave = mean(S_k_kp_norm_split, axis=0)
    T_p_ave = mean(T_p_split, axis=0)
    kx_grid, ky_grid = meshgrid(casek.kx, casek.ky)

    axk=subplot(3,n_case,k_case+1)
    Lamda_abs_max=amax(abs(casek.Lamda_p))
    contourf(kx_grid, ky_grid, casek.Lamda_p.T,levels=linspace(-1*Lamda_abs_max,Lamda_abs_max,nlevel),cmap='seismic')
    colorbar()
    xticks([])
    if k_case == 0:
        ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    title('Coupling Coefficient # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.2f $k_y\\rho_s$=%.2f' % (casek.ftheta_plot[i_theta_plot],kx_select,ky_select) , fontsize=fs2, family='serif')
    #
    axk=subplot(3,n_case,n_case+k_case+1)
    S_k_kp_absmax=amax(abs(S_k_kp_ave))
    contourf(kx_grid,ky_grid,S_k_kp_ave.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    colorbar()
    xticks([])
    if k_case==0:
        ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    title('Bicoherence', fontsize=fs2, family='serif')

    axk=subplot(3,n_case,2*n_case+k_case+1)
    T_p_absmax=amax(abs(T_p_ave))
    contourf(kx_grid,ky_grid,T_p_ave.T,levels=linspace(-1*T_p_absmax,T_p_absmax,nlevel),cmap='seismic')
    colorbar()
    if k_case==0:
        ylabel('$k_y\\rho_s$',fontsize=fs2,family='serif')
    title('$T_{p}$',fontsize=fs2,family='serif')
    xlabel('$k_x\\rho_s$',fontsize=fs2,family='serif')
#
# return the original parameters
plot_nl['t_ave']=t_ave_orig
plot_nl['t_end']=t_end_orig
