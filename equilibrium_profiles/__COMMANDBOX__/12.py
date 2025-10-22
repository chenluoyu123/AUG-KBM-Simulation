# plot the s
figure(figsize=[8,8])
rct=[0.2,0.2,0.6,0.6]
ax=plt.axes(rct)
fs1=24
fs2=20
fs3=16
lw=2
scriptrun=OMFIT['scriptsRun']
ax.plot(scriptrun['rho_sparse'],scriptrun['Alpha'],'-bo',linewidth=lw,label='$\\alpha_{Exp}$')
ax.plot(scriptrun['rho_sparse'],0.6*scriptrun['s'],'-ro',linewidth=lw,label='$\\alpha_{criti,KBM}(0.6s)$')
ax.plot(scriptrun['rho_sparse'],scriptrun['s'],'-ko',linewidth=lw,label='s')
xlabel('$\\rho$',fontsize=fs1,family='serif')
#ylabel('$\\rho$',fontsize=24,family='serif')
xticks(fontsize=fs2,family='serif')
yticks(fontsize=fs2,family='serif')
legend(loc=0,fontsize=fs2).draggable(True)
ylim([0,0.5])