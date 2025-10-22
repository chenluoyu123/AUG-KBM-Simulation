# this script is used to plot the fluctuation intensity over different poloidal theta
# note the dky should be the same for all input cases
import sys
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
# plot
figure(figsize=[15,12])
theta_p=linspace(-pi,pi,37)
i_thetap0=list(theta_p).index(0)  # geometry theta discretation
field_name=['phi_abs','fluct_n_abs','fluct_e_abs','fluct_v_abs']
title_name=['\phi^2','n^2','p^2','v^2']
n_field_plot=3
for k_case in range(n_case):
    print(case_plot[k_case])
    casek = outputs[case_plot[k_case]]
    casek.getbigfield()
    # casek.miller_wd_s(theta_p=theta_p)
#   index specification
    i_s=1  # species
    i_r=list(casek.kx).index(0)  # find the index of kx=0
    i_n=2
    # i_theta0=list(casek.ftheta_plot).index(0)  # the i_theta for theta_plot=0
###
    n_t_ft = len(casek.ind_t_ave)
## caution, don't write to be : phi_cmplx = casek.kxky_phi[0, :, :, :, case.t_ind_ave] + 1j * casek.kxky_phi[1, :, :, :, case.t_ind_ave]
    phi_cmplx = casek.kxky_phi[0, :, :, :, :] + 1j * casek.kxky_phi[1, :, :, :,:] #phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
    fluct_n_cmplx=casek.kxky_n[0,:,:,:,:,:]+1j*casek.kxky_n[1,:,:,:,:,:]   #kxky_n(2, self.n_radial, self.theta_plot, self.n_species, self.n_n, nt), 'F')
    fluct_e_cmplx = casek.kxky_e[0, :, :, :, :, :] +1j* casek.kxky_e[1, :, :, :, :, :]
    # fluct_v_cmplx = casek.kxky_v[0, :, :, :, :, :] + casek.kxky_v[1, :, :, :, :, :]
    phi_cmplx=phi_cmplx[:,:,:,casek.ind_t_ave]/casek.dky/casek.dkx
    fluct_n_cmplx=fluct_n_cmplx[:,:,:,:,casek.ind_t_ave]/casek.dky/casek.dkx
    fluct_e_cmplx = fluct_e_cmplx[:,: , :, :, casek.ind_t_ave] / casek.dky/casek.dkx
    # fluct_v_cmplx = fluct_v_cmplx[:,: , :, :, casek.ind_t_ave] / casek.dky
    phi_abs=mean(abs(phi_cmplx),axis=-1) # phi2_abs[i_r,i_theta_plot,i_n]
    fluct_n_abs = mean(abs(fluct_n_cmplx), axis=-1)  # fluct_n_abs[i_r,i_theta_plot,i_s,i_n]
    fluct_e_abs = mean(abs(fluct_e_cmplx), axis=-1)
    # fluct_v_abs = mean(abs(fluct_v_cmplx), axis=-1)
    fluct_n_abs = fluct_n_abs[:,:,i_s,:]  # fluct_n_abs[i_r,i_theta_plot,i_n]
    fluct_e_abs = fluct_e_abs[:,:,i_s,:]
    # fluct_v_abs = fluct_v_abs[:, :, i_s, :]
    for k_field in arange(n_field_plot):
        subplot(n_field_plot,3,1+n_field_plot*k_field)
        # phi over i_theta for a given i_n and i_r
        cmd='plot(casek.ftheta_plot[0:-1]/pi,'+field_name[k_field]+'[i_r,:,i_n]**2,labo[k_case],linewidth=lw,label=case_plot[k_case])'
        exec(cmd)
        if k_field==0:
            legend(loc=0,fontsize=fs2).draggable(True)
        # plt.gca().set_ylim(bottom=0)
        xticks(fontsize=fs2)
        yticks(fontsize=fs2)
        # xlabel('$\\theta(\pi)$', fontsize=fs1, family='serif')
        # title('$\phi^2(k_y\\rho_s=%.2f,k_x\\rho_s=%.2f)$' % (casek.ky[i_n], casek.kx[i_r]),fontsize=fs1)
        if k_field==n_field_plot-1:
            xlabel('$\\theta(\pi)$',fontsize=fs1)
        if k_field==0:
            title('$k_y\\rho_s=%.2f,k_x\\rho_s=%.2f$' % (casek.ky[i_n], casek.kx[i_r]), fontsize=fs1)
        ylabel('$'+title_name[k_field]+'$',fontsize=fs2)
        # plt.gca().set_ylim(bottom=0)
        subplot(n_field_plot,3,2+n_field_plot*k_field)
        # phi over i_theta for a given i_n summed over kx
        # plot(casek.ftheta_plot[0:-1]/pi,sum(phi_abs[:,:,i_n]*casek.dkx,axis=0)**2,labo[k_case],linewidth=lw,label=case_plot[k_case])
        cmd='plot(casek.ftheta_plot[0:-1] / pi, sum('+field_name[k_field]+'[:, :, i_n] * casek.dkx, axis=0) ** 2, labo[k_case], linewidth=lw,label=case_plot[k_case])'
        exec (cmd)
        xticks(fontsize=fs2)
        yticks(fontsize=fs2)
        # xlabel('$\\theta(\pi)$', fontsize=fs1, family='serif')
        if k_field==n_field_plot-1:
            xlabel('$\\theta(\pi)$',fontsize=fs1)
        if k_field==0:
            title('$k_y\\rho_s=%.2f$' % casek.ky[i_n],fontsize=fs1)
        subplot(n_field_plot,3,3+n_field_plot*k_field)
        # phi over i_theta summed over i_n&i_r
        # can choose to remove the zonal part or not
        # plot(casek.ftheta_plot[0:-1]/pi,sum(phi_abs[:,:,:]*casek.dkx*casek.dky,axis=(0,2))**2,labo[k_case],linewidth=lw,label=case_plot[k_case])
        i_rmzf=0  # remove the zonal component
        if i_rmzf==0:
            cmd='plot(casek.ftheta_plot[0:-1]/pi,sum('+field_name[k_field]+'[:,:,:]*casek.dkx*casek.dky,axis=(0,2))**2,labo[k_case],linewidth=lw,label=case_plot[k_case])'
        else:
            cmd = 'plot(casek.ftheta_plot[0:-1]/pi,sum(' + field_name[\
                k_field] + '[:,:,1:]*casek.dkx*casek.dky,axis=(0,2))**2,labo[k_case],linewidth=lw,label=case_plot[k_case])'
        exec(cmd)
        # plt.gca().set_ylim(bottom=0)
        xticks(fontsize=fs2)
        yticks(fontsize=fs2)        # legend(loc=0, fontsize=fs2).draggable(True)
        if k_field==n_field_plot-1:
            xlabel('$\\theta(\pi)$',fontsize=fs1)
        if k_field==0:
            if i_rmzf==0:
                title('sum over all n',fontsize=fs1,family='serif')
            else:
                title('sum over finite n',fontsize=fs1,family='serif')