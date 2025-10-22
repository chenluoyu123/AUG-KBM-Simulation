# plot the gradient value
figure(figsize=[8,8])
inputgacode=OMFIT['INPUTS']['input.gacode']
a=inputgacode['rmin'][-1]
plot(inputgacode['rho'],a*inputgacode['dlnnedr'],'-k',linewidth=lw)
plot(inputgacode['rho'],a*inputgacode['dlntedr'],'--b',linewidth=lw)
plot(inputgacode['rho'],a*inputgacode['dlntidr_1'],'-b',linewidth=lw)
text(0.30,3.5,'$a/L_{Ti}$',color='b',fontsize=fs2)
text(0.30,1.5,'$a/L_{Te}$',color='b',fontsize=fs2)
text(0.20,0.5,'$a/L_{ne}$',color='k',fontsize=fs2)
xlabel('$\\rho$',fontsize=fs2)
xlim([0,0.9])
ylim([0,4])
xticks(array([0,0.35,0.6,0.9]),fontsize=fs2)
yticks(linspace(0,4,3),fontsize=fs2)