# this script is used to plot the comparison of flux spectrum versus the field amplitude spectrum
pltnl=root['SETTINGS']['PLOTS']['nl']
t_ave=pltnl['t_ave']  # average over [1-t_ave,1]*time_tot
path=pltnl['path']
casename=root['SETTINGS']['PLOTS']['nl']['case_t_trace']
outputs=root['OUTPUTS']
t_ind=int((1-t_ave)*outputs[casename]['n_time'])
fig=figure()
#ax=fig.add_subplot(1,1,1)
fs1=24
fs2=20
fs3=16
lw=2
lab=['-bo','-ro','-ko','-go','-mo','-co']
clr=['b','r','k','g','m','c']
lab2=['--b','--r','--k','--g','--m','--c']
chan=pltnl['chan']  # which channel to be compared, including  Q_*, G_*,P_*, multiple value is available
nchan=len(chan)
field=int(pltnl['field'])  # 0,1,2, only single value is permitted
field_name={0:'phi',1:'apar',2:'bpar'}
dic1={'Q':'energy','G':'particle','P':'momentum'}
casek=outputs[casename]
kyrhos=sign(casek['kyrhos'][1])*casek['kyrhos']
#multfct=10.
for ic in range(nchan):
#    ax1=subplot(2,nchan,ic+1)
    ax1=fig.add_subplot(2,nchan,ic+1)
#    casek=outputs[casename]
    chan_name_in=dic1[chan[ic].split('_')[0]]
    speci=chan[ic].split('_')[1]
    speci=int(speci)
    dky=sign(casek['kyrhos'][1])*casek['kyrhos'][1]
# flux
    flux_ave=mean(casek['flux_ky'][chan_name_in][speci][field].T[t_ind:-1],0)
    flux_std=std(casek['flux_ky'][chan_name_in][speci][field].T[t_ind:-1],0)
    ax1.errorbar(kyrhos,flux_ave/dky,flux_std/dky,color=clr[0],linewidth=lw,label=chan[ic]+'_'+field_name[field])
    if root['SETTINGS']['PLOTS']['ilogx']==1:
        ax1.set_xscale('log')
        xlim([kyrhos[1],kyrhos[-1]])
# phi
    casek.getbigfield()
    fk,ftk = casek.kxky_select(0,field,'phi',0)
    pk=sum(fk[:,:,:],axis=0)/casek.rho
    phi_ave=mean(pk.T[t_ind:-1],0)
    phi_std=std(pk.T[t_ind:-1],0)
    multfct=mean(array(flux_ave))/dky/mean(array(phi_ave))
    errorbar(kyrhos[1:],multfct*phi_ave[1:],multfct*phi_std[1:],color=clr[1],linewidth=lw,label=str(multfct)+'*'+field_name[field])
#    errorbar(kyrhos[1:],phi_ave[1:],phi_std[1:],color=clr[1],linewidth=lw,label=str(multfct)+'*'+field_name[field])
    legend(loc=0,fontsize=fs2).draggable(True)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    title(casename,fontsize=fs1)
    xlabel('$k_y\\rho_s$',fontsize=fs1,family='serif')
#    ax2=subplot(2,nchan,nchan+ic+1)
    ax2=fig.add_subplot(2,nchan,nchan+ic+1)
# divide the flux by phi amplitude for each ky
#    plot(kyrhos,mean(casek['flux_ky'][chan_name_in][speci][field].T[t_ind:-1],0)/dky/mean(pk.T[t_ind:-1],0)**2,'-ro',linewidth=lw,label=chan[ic]+'_'+field_name[field]+'/'+field_name[field])
#    plot(kyrhos,mean(casek['flux_ky'][chan_name_in][speci][field].T[t_ind:-1],0)/dky/mean(pk.T[t_ind:-1],0)**2/kyrhos**2,'-bo',linewidth=lw,label=chan[ic]+'_'+field_name[field]+'/'+'ky/'+field_name[field])
    ax2.plot(kyrhos,mean(casek['flux_ky'][chan_name_in][speci][field].T[t_ind:-1],0)/dky/mean(pk.T[t_ind:-1],0),'-ro',linewidth=lw,label=chan[ic]+'_'+field_name[field]+'/'+field_name[field])
    ax2.plot(kyrhos,mean(casek['flux_ky'][chan_name_in][speci][field].T[t_ind:-1],0)/dky/mean(pk.T[t_ind:-1],0)/kyrhos,'-bo',linewidth=lw,label=chan[ic]+'_'+field_name[field]+'/'+'ky/'+field_name[field])
    if root['SETTINGS']['PLOTS']['ilogx']==1:
        ax2.set_xscale('log')
        xlim([kyrhos[1],kyrhos[-1]])
    legend(loc=0,fontsize=fs2).draggable(True)
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)
    xlabel('$k_y\\rho_s$',fontsize=fs1,family='serif')
