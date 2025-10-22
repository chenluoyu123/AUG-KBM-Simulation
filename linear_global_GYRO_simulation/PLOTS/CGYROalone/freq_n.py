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
for k_case in range(n_case):
    print(case_plot[k_case])
    casek = outputs[case_plot[k_case]]
    casek.get_freq_n()
    # plot(abs(casek['kyrhos']),casek['freq']['omega'].isel(t=arange(t_ind[k_case],casek['n_time'])).mean(axis=1),lab[k_case],linewidth=lw,label='nonlin-'+case_plot[k_case])
    plot(casek.ky,casek.freq_n,labo[k_case],linewidth=lw,label='nonlin-'+case_plot[k_case])
    if case_plot[k_case] in freq_lin.keys():
        casek_lin=case_plot[k_case]
        plot(freq_lin[casek_lin]['ky'],freq_lin[casek_lin]['omega'],lab2[k_case],linewidth=lw,label='lin-'+casek_lin)
xticks(fontsize=fs2)
yticks(fontsize=fs2)
ylabel('$\omega$',fontsize=fs1)
legend(loc=0,fontsize=fs2).draggable(True)
xlabel('$k_y$',fontsize=fs1,family='serif')
ylabel('$\\omega(c_s/a)$',fontsize=fs1,family='serif')