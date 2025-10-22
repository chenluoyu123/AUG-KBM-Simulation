# this script is used to plot the quasilinear flux weight
import numpy as np
import sys
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROscan/assist')
from cgyro_read_xj import *
f = open(root['PLOTS']['CGYROscan']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
# get the data
for paraval_item in Range:
    for ky_item in kyarr:
        datanode=root['OUTPUTScan'][Para][num2str_xj(paraval_item,effnum)]['lin'][num2str_xj(ky_item,effnum)]
        dirname='para_'+num2str_xj(paraval_item,effnum)+'~ky_'+num2str_xj(ky_item,effnum)
        mounttree[dirname]=copy.deepcopy(datanode)

# get the value
nfield=3
Ge=zeros([3,num_ky,nRange]) # 3 is for the number of fields
Gi=zeros([3,num_ky,nRange]) # can be the sum of over species or just one species, depending on i_ion_species, shown below
Qe=zeros([3,num_ky,nRange]) #
Qi=zeros([3,num_ky,nRange]) #
Pe=zeros([3,num_ky,nRange]) #
Pi=zeros([3,num_ky,nRange]) #
i_ion_species=1;                   # -1: sum over the species; positiva value: the index of the ion species, 1 for the fist species
for p_ky in range(num_ky):
    for p_range in range(nRange):
        dirname='para_'+num2str_xj(Range[p_range],effnum)+'~ky_'+num2str_xj(kyarr[p_ky],effnum)
        case = mounttree[dirname]
        qlflux = case['qlflux_ky']
        if setup['icgyro']==1:
            ns = inputcgyro['N_SPECIES']
            inputgen = case['input.cgyro.gen']
        else:
            ns=len(case['tagspec'])+1
            inputgen = case['input.gyro.gen']
        for p_field in range(inputgen['N_FIELD']):
            try:
                if 'AE_FLAG' in inputgen.keys():
                    if inputgen['AE_FLAG'] == 1:
                        Ge[p_field][p_ky][p_range] = 0
                        Qe[p_field][p_ky][p_range] = 0
                        Pe[p_field][p_ky][p_range] = 0
                        if i_ion_species==-1:
                            if setup['icgyro']==0:
                                Gi[p_field][p_ky][p_range] = qlflux['particle'].T[-1][p_field][0]*inputgen['Z']
                                if ns>2:
                                    Gi[p_field][p_ky][p_range]=Gi[p_field][p_ky][p_range]+sum([qlflux['particle'].T[-1][p_field][ind_i]*inputgen['Z_'+str(ind_i+1)] for ind_i in range(1,ns)])
                            else:
                                Gi[p_field][p_ky][p_range] = sum([qlflux['particle'].T[-1][p_field][ind_i] * inputgen['Z_' + str(ind_i + 1)] for \
                                     ind_i in range(0, ns)])
                            Qi[p_field][p_ky][p_range] = sum(qlflux['energy'].T[-1][p_field])
                            Pi[p_field][p_ky][p_range] = sum(qlflux['momentum'].T[-1][p_field])
                        else:
                            Gi[p_field][p_ky][p_range] = qlflux['particle'].T[-1][p_field][i_ion_species - 1]
                            Qi[p_field][p_ky][p_range] = qlflux['energy'].T[-1][p_field][i_ion_species-1]
                            Pi[p_field][p_ky][p_range] = qlflux['momentum'].T[-1][p_field][i_ion_species-1]
                    else:
                        Ge[p_field][p_ky][p_range] = qlflux['particle'][-1][p_field][-1]
                        Qe[p_field][p_ky][p_range] = qlflux['energy'][-1][p_field][-1]
                        Pe[p_field][p_ky][p_range] = qlflux['momentum'][-1][p_field][-1]
                        if i_ion_species == -1:
                            Gi[p_field][p_ky][p_range] = sum(
                                [qlflux['particle'].T[-1][p_field][ind_i] * inputgen['Z_' + str(ind_i + 1)] for ind_i
                                 in range(ns - 1)])
                            Qi[p_field][p_ky][p_range] = sum(qlflux['energy'].T[-1][p_field][0:-1])
                            Pi[p_field][p_ky][p_range] = sum(qlflux['momentum'].T[-1][p_field][0:-1])
                        else:
                            Gi[p_field][p_ky][p_range] = qlflux['particle'].T[-1][p_field][i_ion_species - 1]
                            Qi[p_field][p_ky][p_range] = qlflux['energy'].T[-1][p_field][i_ion_species - 1]
                            Pi[p_field][p_ky][p_range] = qlflux['momentum'].T[-1][p_field][i_ion_species - 1]
                else:
                    Ge[p_field][p_ky][p_range] = qlflux['particle'][-1][p_field][-1]
                    Qe[p_field][p_ky][p_range]=qlflux['energy'][-1][p_field][-1]
                    Pe[p_field][p_ky][p_range] = qlflux['momentum'][-1][p_field][-1]
                    if i_ion_species==-1:
                        Gi[p_field][p_ky][p_range] = sum(
                            [qlflux['particle'].T[-1][p_field][ind_i] * inputgen['Z_' + str(ind_i+1)] for ind_i in
                             range(ns-1)])              # the Z_* weight is verified to be required to put on to calculate the corresponding electron flux
                        Qi[p_field][p_ky][p_range]=sum(qlflux['energy'].T[-1][p_field][0:-1])
                        Pi[p_field][p_ky][p_range]=sum(qlflux['momentum'].T[-1][p_field][0:-1])
                    else:
                        Gi[p_field][p_ky][p_range] = qlflux['particle'].T[-1][p_field][i_ion_species - 1]
                        Qi[p_field][p_ky][p_range]=qlflux['energy'].T[-1][p_field][i_ion_species-1]
                        Pi[p_field][p_ky][p_range]=qlflux['momentum'].T[-1][p_field][i_ion_species-1]
            except:
                Ge[p_field][p_ky][p_range] = 0.
                Qe[p_field][p_ky][p_range]=0.
                Gi[p_field][p_ky][p_range] = 0.
                Qi[p_field][p_ky][p_range]=0.
                Pe[p_field][p_ky][p_range] = 0.
                Pi[p_field][p_ky][p_range]=0.
Channel_arr=['Pe','Pi','Ge','Gi','Qe','Qi']
# plots
for p_field in range(nfield):
    figure('field' + str(p_field))
    for k in range(1, 7):
        axk = subplot(3, 2, k)
        if idimplt==0:
            for p_range in range(nRange):
                print(p_range)
                cmd='axk.plot(kyarr,'+Channel_arr[k-1]+'[p_field].T[p_range],lab[p_range],label=str(Range[p_range]))'
                exec (cmd)
                if k==5 or k==6:
                    xlabel('$k_y$',fontsize=fs1,family='serif')
        else:
            for p_ky in range(num_ky):
                cmd = 'axk.plot(Range,' + Channel_arr[k - 1] + '[p_field][p_ky],lab[p_ky],label=str(kyarr[p_ky]))'
                exec (cmd)
                if k==5 or k==6:
                    xlabel(Para,fontsize=fs1,family='serif')
        axk.set_xscale(scale_dic[ilogx])
        axk.set_yscale(scale_dic[ilogy])
        xticks(fontsize=fs2,family='serif')
        yticks(fontsize=fs2,family='serif')
        if k==1 or k==3 or k==5:
            title(Channel_arr[k-1],fontsize=fs1,family='serif')
        else:
            if i_ion_species==-1:
                tag='-sum'
            else:
                tag='-ion_'+str(i_ion_species)
            title(Channel_arr[k-1]+tag,fontsize=fs1,family='serif')
    legend(loc=0,fontsize=fs2).draggable(True)
# plot the sum over all fields
figure('field-sum')
for k in range(1, 7):
    axk = subplot(3, 2, k)
    if idimplt==0:
        for p_range in range(nRange):
            print(p_range)
            # cmd='axk.plot(kyarr,'+Channel_arr[k-1]+'[p_field].T[p_range],lab[p_range],label=str(Range[p_range]))'
            cmd='axk.plot(kyarr,sum('+Channel_arr[k-1]+', axis=0).T[p_range],lab[p_range],linewidth=lw,label=str(Range[p_range]))'
            exec (cmd)
            if k==5 or k==6:
                xlabel('$k_y$',fontsize=fs1,family='serif')
    else:
        for p_ky in range(num_ky):
            cmd = 'axk.plot(Range,sum(' + Channel_arr[k - 1] + ',axis=0)[p_ky],lab[p_ky],linewidth=lw,label=str(kyarr[p_ky]))'
            exec (cmd)
            if k==5 or k==6:
                xlabel(Para,fontsize=fs1,family='serif')
    axk.set_xscale(scale_dic[ilogx])
    axk.set_yscale(scale_dic[ilogy])
    xticks(fontsize=fs2,family='serif')
    yticks(fontsize=fs2,family='serif')
    if k==1 or k==3 or k==5:
        title(Channel_arr[k-1],fontsize=fs1,family='serif')
    else:
        if i_ion_species==-1:
            tag='-sum'
        else:
            tag='-ion_'+str(i_ion_species)
        title(Channel_arr[k-1]+tag,fontsize=fs1,family='serif')
legend(loc=0,fontsize=fs2).draggable(True)
#
figure('field-sum-ratio to Qe')
for k in range(1, 7):
    axk = subplot(3, 2, k)
    if idimplt==0:
        for p_range in range(nRange):
            print(p_range)
            # cmd='axk.plot(kyarr,'+Channel_arr[k-1]+'[p_field].T[p_range],lab[p_range],label=str(Range[p_range]))'
            cmd='axk.plot(kyarr,sum('+Channel_arr[k-1]+', axis=0).T[p_range]/sum(Qe,axis=0).T[p_range],lab[p_range],linewidth=lw,label=str(Range[p_range]))'
            exec (cmd)
            if k==5 or k==6:
                xlabel('$k_y$',fontsize=fs1,family='serif')
    else:
        for p_ky  in range(num_ky):
            cmd = 'axk.plot(Range,sum(' + Channel_arr[k - 1] + ',axis=0)[p_ky]/sum(Qe,axis=0)[p_ky],lab[p_ky],linewidth=lw,label=str(kyarr[p_ky]))'
            exec (cmd)
            if k==5 or k==6:
                xlabel(Para,fontsize=fs1,family='serif')
    axk.set_xscale(scale_dic[ilogx])
    axk.set_yscale(scale_dic[ilogy])
    xticks(fontsize=fs2,family='serif')
    yticks(fontsize=fs2,family='serif')
    if k==1 or k==3 or k==5:
        title(Channel_arr[k-1],fontsize=fs1,family='serif')
    else:
        if i_ion_species==-1:
            tag='-sum'
        else:
            tag='-ion_'+str(i_ion_species)
        title(Channel_arr[k-1]+tag,fontsize=fs1,family='serif')
legend(loc=0,fontsize=fs2).draggable(True)

# write the flux out
iwriteflux=plots['iwriteflux']
if iwriteflux==1:
    fluxout=root['SETTINGS']['DEPENDENCIES']['fluxout']
    fid=open(fluxout,'w')
# firstly write the scanned parameter name into the file
    fid.write(Para)
    fid.write('\n')
# write the nRange and the para_val into the file
    fid.write(str(nRange))
    fid.write('\n')
    line=''
    for m in range(nRange):
        line=line+str(Range[m])+'    '
    fid.write(line)
    fid.write('\n')
# write the number of poloidal modes into the file for flux summation
    fid.write(str(nfield)+'    '+str(num_ky))
    fid.write('\n')
# write the flux spectrum into the file
# fromation: row number, num_ky
# colum: ky, ((pflux,Qe, Qi, Pi)*nmodes)*nRange)
    for k in range(num_ky):
        line=str(kyarr[k])
        for p in range(nRange):
            for n in range(nfield):
                line=line+'    '+str(Ge[n][k][p])+'    '+\
                                 str(Qe[n][k][p])+'    '+\
                                 str(Qi[n][k][p])+'    '+\
                                 str(Pi[n][k][p])
        fid.write(line)
        fid.write('\n')
    fid.close()
    print('write data done!')
