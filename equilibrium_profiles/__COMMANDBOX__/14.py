# plot the diamagnetic frequency
iptg=OMFIT['WholeInfo']['INPUTS']['input.gacode']
# diamagentic frequency= (a/Lne+a/LTe)*cs/a*rho_s
figure(figsize=[8,8])
lw=3
fs1=32
fs2=28
v_exp=60; #km/s ## get fro project: /fusion/projects/omfit-results/jianx/bispectrum/185959_bes_crossphase.zip
a=iptg['rmin'][-1]
BunitOverBt=abs(iptg['bunit']/iptg['bt0']) # bt0 is the toroidal field at theta=0
v_diag=(iptg['dlnnedr']+iptg['dlntedr'])*iptg['cs']*iptg['rhos']/1.e3
v_Er=iptg['rmin']*iptg['omega0']/1.e3
plot(iptg['rho'],BunitOverBt*v_diag,'-b',linewidth=lw,label='Diamagnetic')
plot(iptg['rho'],BunitOverBt*v_Er,'-k',linewidth=lw,label='Er')
plot(iptg['rho'],BunitOverBt*(v_diag+v_Er),'-r',linewidth=lw,label='Tot')
plot(array([0,1]),array([v_exp,v_exp]),'--m',linewidth=lw)
fill(array([0,1,1,0]),array([54,54,66,66]),color='m',alpha=0.2)
text(0.98,15,'$v_{*pe}$',fontsize=fs2,color='b')
text(0.97,30,'$v_{Er}$',fontsize=fs2,color='k')
text(0.955,42,'$v_{tot}$',fontsize=fs2,color='r')
text(0.98,62,'$v_{BES}$',fontsize=fs2,color='m')
text(0.945,70,'(d)',fontsize=fs2,color='k')
xlim([0.94,1.0])
ylim([0,80])
xticks(linspace(0.95,1.0,6),fontsize=fs2)
yticks(linspace(0,80,5),fontsize=fs2)
ylabel('km/s',fontsize=fs1)
xlabel('$\\rho$',fontsize=fs1)