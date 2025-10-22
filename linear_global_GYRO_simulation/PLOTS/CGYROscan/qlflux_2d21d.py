# this script is used to plot the quasilinear flux.
import sys
import xarray as xr
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROscan/assist')
from cgyro_read_xj import *
f = open(root['PLOTS']['CGYROscan']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()

for paraxval_item in Range_x:
    for parayval_item in Range_y:
        for ky_item in kyarr:
            datanode=root['OUTPUTScan'][Para_x][num2str_xj(paraxval_item,effnum)][Para_y][num2str_xj(parayval_item,effnum)]['lin'][num2str_xj(ky_item,effnum)]
            dirname='parax_'+num2str_xj(paraxval_item,effnum)+'~paray_'+num2str_xj(parayval_item,effnum)+'~ky_'+num2str_xj(ky_item,effnum)
            mounttree[dirname]=copy.deepcopy(datanode)

# get the value
nfield=3
Qe_data=zeros([nfield,num_ky,nRange_x,nRange_y]) # 3 is for the number of fields
Qi_data=zeros([nfield,num_ky,nRange_x,nRange_y]) # the results of a given ion species
Ge_data=zeros([nfield,num_ky,nRange_x,nRange_y])
Gi_data=zeros([nfield,num_ky,nRange_x,nRange_y])
Pe_data=zeros([nfield,num_ky,nRange_x,nRange_y])
Pi_data=zeros([nfield,num_ky,nRange_x,nRange_y]) # for a given ion species
i_ion_species=3
#
for p_ky in range(num_ky):
    for p_rangex in range(nRange_x):
        for p_rangey in range(nRange_y):
            dirname='parax_'+num2str_xj(Range_x[p_rangex],effnum)+'~paray_'+num2str_xj(Range_y[p_rangey],effnum)+'~ky_'+num2str_xj(kyarr[p_ky],effnum)
            case = mounttree[dirname]
            qlflux = case['qlflux_ky']
            if setup['icgyro'] == 1:
                inputgen = mounttree[dirname]['input.cgyro.gen']
            else:
                inputgen = mounttree[dirname]['input.gyro.gen']
            for p_field in range(nfield):
                try:
                    if 'AE_FLAG' in inputgen.keys():
                        if inputgen['AE_FLAG'] == 1:
                            Qe_data[p_field][p_ky][p_rangex][p_rangey] = 0
                            Qi_data[p_field][p_ky][p_rangex][p_rangey] = sum(qlflux['energy'].T[-1][p_field])
                            Ge_data[p_field][p_ky][p_rangex][p_rangey] = 0
                            Gi_data[p_field][p_ky][p_rangex][p_rangey] = 0
                            Pe_data[p_field][p_ky][p_rangex][p_rangey] = 0
                            Pi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['momentum'].T[-1][p_field][i_ion_species]
                        else:
                            Ge_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['particle'][-1][p_field][-1]
                            Qe_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['energy'][-1][p_field][-1]
                            # Qi_data[p_field][p_ky][p_rangex][p_rangey] = sum(qlflux['energy'].T[-1][p_field][0:-1])
                            # Ge_data[p_field][p_ky][p_rangex][p_rangey] = sum(qlflux['particle'][-1][p_field][-1])
                            # Pi_data[p_field][p_ky][p_rangex][p_rangey] = sum(qlflux['momentum'].T[-1][p_field][0:-1])
                            Gi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['particle'].T[-1][p_field][i_ion_species-1]
                            Qi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['energy'].T[-1][p_field][i_ion_species-1]
                            Pe_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['momentum'].T[-1][p_field][- 1]
                            Pi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['momentum'].T[-1][p_field][i_ion_species-1]
                    else:
                        Ge_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['particle'][-1][p_field][-1]
                        Gi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['particle'].T[-1][p_field][i_ion_species-1]
                        Qe_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['energy'][-1][p_field][-1]
                        # Qi_data[p_field][p_ky][p_rangex][p_rangey] = sum(qlflux['energy'].T[-1][p_field][0:-1])
                        Qi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['energy'][-1][p_field][i_ion_species-1]
                        Pe_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['momentum'].T[-1][p_field][-1]
                        # Pi_data[p_field][p_ky][p_rangex][p_rangey] = sum(qlflux['momentum'].T[-1][p_field][0:-1])
                        Pi_data[p_field][p_ky][p_rangex][p_rangey] = qlflux['momentum'].T[-1][p_field][i_ion_species-1]
                except:
                    Ge_data[p_field][p_ky][p_rangex][p_rangey] = 0.
                    Qe_data[p_field][p_ky][p_rangex][p_rangey] = 0.
                    Gi_data[p_field][p_ky][p_rangex][p_rangey] = 0.
                    Qi_data[p_field][p_ky][p_rangex][p_rangey] = 0.
                    Pe_data[p_field][p_ky][p_rangex][p_rangey] = 0.
                    Pi_data[p_field][p_ky][p_rangex][p_rangey] = 0.
# construct to xarray
Qe=xr.DataArray(data=Qe_data,
                coords={
                    "field":range(nfield),
                    "ky": kyarr,
                    "rangex": Range_x,
                    "rangey": Range_y,
                },
                dims=['field','ky','rangex','rangey'],
                )
Qi=xr.DataArray(data=Qi_data,
                coords={
                    "field":range(nfield),
                    "ky": kyarr,
                    "rangex": Range_x,
                    "rangey": Range_y,
                },
                dims=['field','ky','rangex','rangey'],
                )
Ge=xr.DataArray(data=Ge_data,
                coords={
                    "field":range(nfield),
                    "ky": kyarr,
                    "rangex": Range_x,
                    "rangey": Range_y,
                },
                dims=['field','ky','rangex','rangey'],
                )
Gi=xr.DataArray(data=Gi_data,
                coords={
                    "field":range(nfield),
                    "ky": kyarr,
                    "rangex": Range_x,
                    "rangey": Range_y,
                },
                dims=['field','ky','rangex','rangey'],
                )

Pe=xr.DataArray(data=Pe_data,
                coords={
                    "field":range(nfield),
                    "ky": kyarr,
                    "rangex": Range_x,
                    "rangey": Range_y,
                },
                dims=['field','ky','rangex','rangey'],
                )
Pi=xr.DataArray(data=Pi_data,
                coords={
                    "field":range(nfield),
                    "ky": kyarr,
                    "rangex": Range_x,
                    "rangey": Range_y,
                },
                dims=['field','ky','rangex','rangey'],
                )
Channel_arr=['Pe','Pi','Ge','Gi','Qe','Qi']
if Para_x in Paraname_dic.keys():
    plt2d['Para_x_display']=Paraname_dic[Para_x]
else:
    plt2d['Para_x_display']=Para_x
if Para_y in Paraname_dic.keys():
    plt2d['Para_y_display']=Paraname_dic[Para_y]
else:
    plt2d['Para_y_display']=Para_y
idisplaysum_only=1
# based on the profiles we get, then all the linear information can be plotted
Para_x_display=plt2d['Para_x_display']
Para_y_display=plt2d['Para_y_display']
Para_x_dic={0:Para_x_display, 1: Para_y_display}
Para_y_dic={0:Para_y_display, 1: Para_x_display}
if idisplaysum_only==0:
    for p_field in range(nfield):
        for p_ky in range(num_ky):
            figure('field'+str(p_field)+'~ky'+str(kyarr[p_ky]),figsize=[16,10])
            for k in arange(1,7):
                axk=subplot(3,2,k)
                if idimplt==0:
                    cmd='plot(Range_x,'+Channel_name[k-1]+'.isel(field=p_field,ky=p_ky,rangey=p_rangey),lab[p_rangey],label=str(Range_y[p_rangey]))'
                    for p_rangey in range(nRange_y):
                        exec(cmd)
                else:
                    cmd = 'plot(Range_y,' + Channel_name[k - 1] + '.isel(field=p_field,ky=p_ky,rangex=p_rangex),lab[p_rangex],label=str(Range_x[p_rangex]))'
                    for p_rangex in range(nRange_x):
                        exec(cmd)
                axk.set_xscale(scale_dic[ilogx])
                axk.set_yscale(scale_dic[ilogy])
                legend(loc=0, fontsize=fs2).draggable(True)
                if k==5 or k==6:
                    xlabel(Para_x_dic[idimplt], fontsize=fs1, family='serif')
                xticks(fontsize=fs2,family='serif')
                yticks(fontsize=fs2,family='serif')
                if k == 1 or k == 3 or k == 5:
                    title(Channel_arr[k - 1], fontsize=fs1, family='serif')
                else:
                    if i_ion_species == -1:
                        tag = '-sum'
                    else:
                        tag = '-ion_' + str(i_ion_species)
                    title(Channel_arr[k - 1] + tag, fontsize=fs1, family='serif')

# fieldsum
for p_ky in range(num_ky):
    figure('fieldsum~ky'+str(kyarr[p_ky]),figsize=[16,10])
    for k in arange(1,7):
        axk=subplot(3,2,k)
        if idimplt==0:
            cmd='plot(Range_x,sum('+Channel_name[k-1]+'.isel(ky=p_ky,rangey=p_rangey),axis=0),lab[p_rangey],label=str(Range_y[p_rangey]),linewidth=lw*2)'
            for p_rangey in range(nRange_y):
                exec(cmd)
        else:
            cmd = 'plot(Range_y,sum(' + Channel_name[k - 1] + '.isel(ky=p_ky,rangex=p_rangex),axis=0),lab[p_rangex],label=str(Range_x[p_rangex]),linewidth=lw*2)'
            for p_rangex in range(nRange_x):
                exec(cmd)
        axk.set_xscale(scale_dic[ilogx])
        axk.set_yscale(scale_dic[ilogy])
        legend(loc=0, fontsize=fs2).draggable(True)
        if k==5 or k==6:
            xlabel(Para_x_dic[idimplt], fontsize=fs1, family='serif')
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


# fieldsum to Qe ratio
for p_ky in range(num_ky):
    figure('fieldsum2Qe~ky'+str(kyarr[p_ky]),figsize=[16,10])
    for k in arange(1,7):
        axk=subplot(3,2,k)
        if idimplt==0:
            cmd='plot(Range_x,sum('+Channel_name[k-1]+'.isel(ky=p_ky,rangey=p_rangey),axis=0)/sum(Qe.isel(ky=p_ky,rangey=p_rangey),axis=0),lab[p_rangey],label=str(Range_y[p_rangey]),linewidth=lw*2)'
            for p_rangey in range(nRange_y):
                exec(cmd)
        else:
            cmd = 'plot(Range_y,sum(' + Channel_name[k - 1] + '.isel(ky=p_ky,rangex=p_rangex),axis=0)/sum(Qe.isel(ky=p_ky,rangex=p_rangex),axis=0),lab[p_rangex],label=str(Range_x[p_rangex]),linewidth=lw*2)'
            for p_rangex in range(nRange_x):
                exec(cmd)
        axk.set_xscale(scale_dic[ilogx])
        axk.set_yscale(scale_dic[ilogy])
        legend(loc=0, fontsize=fs2).draggable(True)
        if k==5 or k==6:
            xlabel(Para_x_dic[idimplt], fontsize=fs1, family='serif')
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


# write the flux out
iwriteflux=plots['iwriteflux']
if iwriteflux==1:
    fluxout=root['SETTINGS']['DEPENDENCIES']['fluxout']
    fid=open(fluxout,'w')
# firstly write the scanned parameter name into the file
    fid.write(para)
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
