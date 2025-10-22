# plot vf/vA
# execulate normpara.py first
figure(figsize=[6,6])
fs2=24
case=OMFIT['scriptsRun']
plot(case['rho_sparse'],case['Vf']/case['VA'],'-ro',linewidth=lw)
xlim([0.1,0.9])
ylim([0,0.3])
xlabel('$\\rho$',fontsize=fs2)
title('$V_f/V_A$',fontsize=fs2)
xticks(fontsize=fs2)
yticks(linspace(0,0.3,4),fontsize=fs2)