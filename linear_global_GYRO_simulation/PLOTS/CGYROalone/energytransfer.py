# for plotting the energy transfer to a selected (kx_select, ky_select)
import sys
from matplotlib import ticker, cm
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROalone/assist')
f = open(root['PLOTS']['CGYROalone']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
icgyro=root['SETTINGS']['SETUP']['icgyro']
if icgyro==1:
    root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
else:
    root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
# check to see how many fields that all the cases have in common
n_field=3
for k_case in range(n_case):
    casek = outputs[case_plot[k_case]]
    if icgyro==1:
        n_field=min(len(casek['field_tags']),n_field)
    else:
        n_field=min(len(casek['tagfieldtext']),n_field)
nlevel=128
kx_select=0.0
ky_select=0.12
# note that it could be better if we split the time window in to several sub time windows
plot_nl=root['SETTINGS']['PLOTS']['nl']
t_ave_orig=plot_nl['t_ave']
t_end_orig=plot_nl['t_end']
n_split=2
ifilter=0
iplotp=0 # whether plot the internal energy transfer
i_improve=0 # if 1, use the improve version which can accelate the calculation but the exact value is not correct; 
i_cal_ratio=1 # calculate T_phi(ky'=0)/(ky**2*phi**2)
icontour=1 # choose to use contour or pcolor
ky_bdry=0.2
ms=2
# # plot Energy transfer of phi (kinetic energy)
figure(figsize=[15,16])
for k_case in range(n_case):
    casek = outputs[case_plot[k_case]]
    t_ave = t_ave_orig / n_split
    print(case_plot[k_case])
    if icgyro==1:
        i_theta_plot=casek.theta_plot//2
    else:
        i_theta_plot=casek.n_theta_plot//2
    S_k_kp_split = zeros([n_split, casek.n_r-1, casek.n_n])
#    S_k_kp_norm_split = zeros([n_split, casek.n_r-1, casek.n_n])
    T_phi_split = zeros([n_split, casek.n_r-1, casek.n_n])
    for i_split in arange(n_split):
        t_end=t_end_orig-t_ave*i_split
        plot_nl['t_ave']=t_ave
        plot_nl['t_end'] = t_end
        print(t_ave,t_end)
        # casek = outputs[case_plot[k_case]]  # datanode
        if icgyro==1:
            root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
        else:
            root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
        if i_improve==1:
            casek.Energy_transfer_phi_improve(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
        else:
            casek.Energy_transfer_phi(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
        S_k_kp_split[i_split]=casek.S_k_kp_phi
#        S_k_kp_norm_split[i_split] = casek.S_k_kp_norm_phi
        T_phi_split[i_split] = casek.T_phi
    S_k_kp_ave=mean(S_k_kp_split,axis=0)
#    S_k_kp_norm_ave = mean(S_k_kp_norm_split, axis=0)
    T_phi_ave = mean(T_phi_split, axis=0)
    kx_grid, ky_grid = meshgrid(casek.kx, casek.ky)
    ax1=subplot(3,n_case,k_case+1)
    Lamda_abs_max=amax(abs(casek.Lamda_phi))
    if icontour==1:
        contourf(kx_grid, ky_grid, casek.Lamda_phi.T,levels=linspace(-1*Lamda_abs_max,Lamda_abs_max,nlevel),cmap='seismic')
    else:
        pcolor(kx_grid, ky_grid, casek.Lamda_phi.T,cmap='seismic',vmin=-1*Lamda_abs_max,vmax=Lamda_abs_max)
    plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
    plot(casek.kx_select_incode,casek.ky_select_incode,'o',markersize=ms*6,color='k')
    colorbar()
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    xlim([-1*ky_bdry,ky_bdry])
    if k_case == 0:
        ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    if icgyro==1:
        title('Coupling Coefficient # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')
    else:
        title('Coupling Coefficient # $k_x\\rho_s$=%.2f $k_y\\rho_s$=%.2f' % (kx_select,ky_select) , fontsize=fs2, family='serif')

    ax2=subplot(3,n_case,n_case+k_case+1)
    S_k_kp_absmax=amax(abs(S_k_kp_ave))
    if icontour==1:
        contourf(kx_grid,ky_grid,S_k_kp_ave.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    else:
        pcolor(kx_grid,ky_grid,S_k_kp_ave.T,cmap='seismic',vmin=-1*S_k_kp_absmax,vmax=S_k_kp_absmax)
    colorbar()
    plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
    plot(casek.kx_select_incode,casek.ky_select_incode,'o',markersize=ms*6,color='k')
    xlim([-1*ky_bdry,ky_bdry])
    xticks([],fontsize=fs2)
    yticks(fontsize=fs2)
    if k_case==0:
        ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    title('Bicoherence # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')

    # axk=subplot(4,n_case,2*n_case+k_case+1)
    # S_k_kp_absmax=amax(abs(casek.S_k_kp_norm))
    ## contourf(kx_grid,ky_grid,casek.S_k_kp_norm.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    # contourf(kx_grid,ky_grid,S_k_kp_norm_ave.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
    # colorbar()
    # if k_case==0:
    #     ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
    # title('Bicoherence(norm)', fontsize=fs2, family='serif')

    ax3=subplot(3,n_case,2*n_case+k_case+1)
# in order to filter out the noise induced by high kx,ky, T_phi can be multiplied by the coherence factor
    if ifilter==1:
        T_phi_ave=T_phi_ave*abs(S_k_kp_ave)**2
    T_phi_absmax=amax(abs(T_phi_ave))
    if icontour==1:
        contourf(kx_grid,ky_grid,T_phi_ave.T,levels=linspace(-1*T_phi_absmax,T_phi_absmax,nlevel),cmap='seismic')    
    else:
        pcolor(kx_grid,ky_grid,T_phi_ave.T,cmap='seismic',vmin=-1*T_phi_absmax,vmax=T_phi_absmax)    
    colorbar()
    plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
    plot(casek.kx_select_incode,casek.ky_select_incode,'o',markersize=ms*6,color='k')
    xlim([-1*ky_bdry,ky_bdry])
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    if ifilter==1:
        title('T_{\phi}-Filter # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')
    else:
        title('$T_{\phi}$ # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')
    xlabel('$k_x\\rho_s$',fontsize=fs2,family='serif')
    ax1.get_shared_x_axes().join(ax2,ax3)
    ax1.get_shared_y_axes().join(ax2,ax3)
    ax2.get_shared_x_axes().join(ax1,ax3)
    ax2.get_shared_y_axes().join(ax1,ax3)
    ax3.get_shared_x_axes().join(ax1,ax2)
    ax3.get_shared_y_axes().join(ax1,ax2)
    # get the value of T_phi/(ky**2*phi_kx_0**2)
    if i_cal_ratio==1:
        print('please compare this value with the linear growth rate: ky=', ky_select)
        fk=casek.phi_cmplx_t[1:,i_theta_plot,:,:]
        phi_ave_t=mean(abs(fk[:,:,casek.ind_t_ave]),axis=-1)/casek.rho
        ky_select_incode=casek.ky.flat[np.abs(casek.ky-ky_select).argmin()]  # the real selected kx value in code
        ind_ky=where(casek.ky==ky_select_incode)[0][0]
        phi_n_kx0=phi_ave_t[len(casek.kx)//2,ind_ky]
        print(sum(T_phi_ave,axis=0)[0]/(ky_select**2*phi_n_kx0**2))
# plot Energy transfer of p (internal energy)
if iplotp==1:
    figure(figsize=[15,8])
    for k_case in range(n_case):
        casek = outputs[case_plot[k_case]]
        t_ave = t_ave_orig / n_split
        print(case_plot[k_case])
        if icgyro==1:
            i_theta_plot=casek.theta_plot//2
        else:
            i_theta_plot=casek.n_theta_plot//2
        S_k_kp_split = zeros([n_split, casek.n_r-1, casek.n_n])
        S_k_kp_norm_split = zeros([n_split, casek.n_r-1, casek.n_n])
        T_p_split = zeros([n_split, casek.n_r-1, casek.n_n])
        for i_split in arange(n_split):
            t_end=t_end_orig-t_ave*i_split
            plot_nl['t_ave']=t_ave
            plot_nl['t_end'] = t_end
    #        casek = outputs[case_plot[k_case]]
            if icgyro==1:
                root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
            else:
                root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
            casek.Energy_transfer_p_improve(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
            S_k_kp_split[i_split]=casek.S_k_kp_p
            S_k_kp_norm_split[i_split] = casek.S_k_kp_norm_p
            T_p_split[i_split] = casek.T_p
        S_k_kp_ave=mean(S_k_kp_split,axis=0)
        S_k_kp_norm_ave = mean(S_k_kp_norm_split, axis=0)
        T_p_ave = mean(T_p_split, axis=0)
        kx_grid, ky_grid = meshgrid(casek.kx, casek.ky)

        ax1=subplot(3,n_case,k_case+1)
        Lamda_abs_max=amax(abs(casek.Lamda_p))
        if icontour==1:
            contourf(kx_grid, ky_grid, casek.Lamda_p.T,levels=linspace(-1*Lamda_abs_max,Lamda_abs_max,nlevel),cmap='seismic')
        else:
            pcolor(kx_grid, ky_grid, casek.Lamda_p.T,cmap='seismic',vmin=-1*Lamda_abs_max,vmax=Lamda_abs_max)
        colorbar()
        plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
        plot(casek.kx_select_incode,casek.ky_select_incode,'o',markersize=ms*6,color='k')
        xlim([-1*ky_bdry,ky_bdry])
        xticks([],fontsize=fs2)
        yticks(fontsize=fs2)
        if k_case == 0:
            ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
        if icgyro==1:
            title('Coupling Coefficient # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')
        else:
            title('Coupling Coefficient # $$k_x\\rho_s$=%.2f $k_y\\rho_s$=%.2f' % (kx_select,ky_select) , fontsize=fs2, family='serif')
        #
        ax2=subplot(3,n_case,n_case+k_case+1)
        S_k_kp_absmax=amax(abs(S_k_kp_ave))
        if icontour==1:
            contourf(kx_grid,ky_grid,S_k_kp_ave.T,levels=linspace(-1*S_k_kp_absmax,S_k_kp_absmax,nlevel),cmap='seismic')
        else:
            pcolor(kx_grid,ky_grid,S_k_kp_ave.T,cmap='seismic',vmin=-1*S_k_kp_absmax,vmax=S_k_kp_absmax)
        colorbar()
        plot(casek.kx_select_incode,casek.ky_select_incode,'o',markersize=ms*6,color='k')
        xlim([-1*ky_bdry,ky_bdry])
        plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
        xticks([],fontsize=fs2)
        yticks(fontsize=fs2)
        if k_case==0:
            ylabel('$k_y\\rho_s$', fontsize=fs2, family='serif')
        title('Bicoherence # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')

        ax3=subplot(3,n_case,2*n_case+k_case+1)
        T_p_absmax=amax(abs(T_p_ave))
        if icontour==1:
            contourf(kx_grid,ky_grid,T_p_ave.T,levels=linspace(-1*T_p_absmax,T_p_absmax,nlevel),cmap='seismic')
        else:
            pcolor(kx_grid,ky_grid,T_p_ave.T,cmap='seismic',vmin=-1*T_p_absmax,vmax=T_p_absmax)
        colorbar()
        xlim([-1*ky_bdry,ky_bdry])
        xticks(fontsize=fs2)
        yticks(fontsize=fs2)
        plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
        plot(casek.kx_select_incode,casek.ky_select_incode,'o',markersize=ms*6,color='k')
        if k_case==0:
            ylabel('$k_y\\rho_s$',fontsize=fs2,family='serif')
        title('T_p # $\\theta_{plot}$=%.2f $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (casek.ftheta_plot[i_theta_plot],casek.kx_select_incode,casek.ky_select_incode) , fontsize=fs2, family='serif')
        xlabel('$k_x\\rho_s$',fontsize=fs2,family='serif')
        ax1.get_shared_x_axes().join(ax2,ax3)
        ax1.get_shared_y_axes().join(ax2,ax3)
        ax2.get_shared_x_axes().join(ax1,ax3)
        ax2.get_shared_y_axes().join(ax1,ax3)
        ax3.get_shared_x_axes().join(ax1,ax2)
        ax3.get_shared_y_axes().join(ax1,ax2)
#
# return the original parameters
plot_nl['t_ave']=t_ave_orig
plot_nl['t_end']=t_end_orig
