# this script is used to plot the eigenfrequency and growth rate for all the scanning cases,
import sys
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROscan/assist')
from cgyro_read_xj import *
f = open(root['PLOTS']['CGYROscan']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()

# then read all the data
w_arr=zeros([nRange_x,nRange_y,num_ky])     # dominate mode frequency
err_w_arr=zeros([nRange_x,nRange_y,num_ky])     # error in dominate mode frequency
gamma_arr=zeros([nRange_x,nRange_y,num_ky]) # dominate mode growth rate
err_gamma_arr=zeros([nRange_x,nRange_y,num_ky]) # error in dominate mode growth rate
for k in range(nRange_x):
    for p in range(nRange_y):
        for q in range(num_ky):
            try:
                datanode=root['OUTPUTScan'][Para_x][num2str_xj(Range_x[k],effnum)][Para_y][num2str_xj(Range_y[p],effnum)]['lin'][num2str_xj(kyarr[q],effnum)]
                flag_read=1
                w,err_w=readfreq(datanode,flag_read)
                if ibelow0[0]==1:
                    if w[1]<ibelow0[1]:
                        w[1]=ibelow0[1]
            except:
#            print('not work.')
                w=array([0,0])
                err_w=array([1,1])
            w_arr[k][p][q]=w[0]
            gamma_arr[k][p][q]=w[1]
            err_w_arr[k][p][q]=err_w[0]
            err_gamma_arr[k][p][q]=err_w[1]
########################################
# based on the profiles we get, then all the linear information can be plotted
Para_x_dic={0:Para_x, 1: Para_y}
Para_y_dic={0:Para_y, 1: Para_x}
bdry_dic={1:0, 2:0, 3:1.e-4, 4:1.e-4}
freqarr=['w_arr','err_w_arr','gamma_arr','err_gamma_arr']
if Para_x in Paraname_dic.keys():
    plt2d['Para_x_display']=Paraname_dic[Para_x]
else:
    plt2d['Para_x_display']=Para_x
if Para_y in Paraname_dic.keys():
    plt2d['Para_y_display']=Paraname_dic[Para_y]
else:
    plt2d['Para_y_display']=Para_y
    
# allows for scale of range_x and range_y
x_scale=plt2d['x_scale']
y_scale=plt2d['y_scale']
Para_x_display=plt2d['Para_x_display']
Para_y_display=plt2d['Para_y_display']
# based on the profiles we get, then all the linear information can be plotted
Para_x_dic={0:Para_x_display, 1: Para_y_display}
Para_y_dic={0:Para_y_display, 1: Para_x_display}
bdry_dic=[0.3, 1.e-4, 0.10, 1.e-4]
freqarr=['w_arr','err_w_arr','gamma_arr','err_gamma_arr']

for k in range(num_ky):
    fig=figure(str(kyarr[k]),figsize=[10,10])
    for p in arange(1,5):
        axp=subplot(2,2,p)
        if idimplt==0:
            cmd='plot(Range_x,'+freqarr[p-1]+'.T[k][m],lab[m],label=num2str_xj(Range_y[m],effnum),linewidth=lw)'
            for m in range(0,nRange_y):
                exec(cmd)
        else:
            cmd='plot(Range_y,'+freqarr[p-1]+'.T[k].T[m],lab[m],label=num2str_xj(Range_x[m],effnum),linewidth=lw)'
            for m in range(0,nRange_x):
                exec(cmd)
        legend(loc=0,fontsize=fs2).draggable(True)
        if p==3 or p==4:
            xlabel(Para_x_dic[idimplt],fontsize=fs1,family='serif')
        axp.set_xscale(scale_dic[ilogx])
        axp.set_yscale(scale_dic[ilogy])
        if p==2 or p==4:
            ylim([1.e-6,1.0])
            axp.set_yscale('log')
        xticks(fontsize=fs2,family='serif')
        yticks(fontsize=fs2,family='serif')
        title(freqarr[p-1],fontsize=fs1,family='serif')
# # ===================================
iwritelin=root['SETTINGS']['PLOTS']['iwritelin']  # determine whether to write out to a file
numky=len(kyarr)
nmodes=1
if iwritelin==1 and 'linout' in root['SETTINGS']['DEPENDENCIES'].keys():
    eigenout=root['SETTINGS']['DEPENDENCIES']['linout']
    fid=open(eigenout,'w')
# first write the scanned parameter name into the file
    fid.write(Para_x+'    '+Para_y)
    fid.write('\n')
# write the nRange and the para_val into the file
    fid.write(str(len(Range_x))+'    '+str(len(Range_y)))
    fid.write('\n')
    line=''
    for m in range(len(Range_x)):
        line=line+str(Range_x[m])+'    '
    fid.write(line)
    fid.write('\n')
    line=''
    for m in range(len(Range_y)):
        line=line+str(Range_y[m])+'    '
    fid.write(line)
    fid.write('\n')
# write nmodes into the file
    fid.write(str(nmodes)+'    '+str(numky))
    fid.write('\n')
    for k in range(numky):
        line=str(kyarr[k])
        for p in range(len(Range_y)):
            for q in range(len(Range_x)):
# w_arr=zeros([nRange_x,nRange_y,num_ky])
#            line=line+'    '+str(w_arr[p][k])+'    '+str(gamma_arr[p][k])+'    '+str(err_w_arr[p][k])+'    '+str(err_gamma_arr[p][k])
                line=line+'    '+str(w_arr[q][p][k])+'    '+str(gamma_arr[q][p][k])
        fid.write(line)
        fid.write('\n')
    fid.close()
