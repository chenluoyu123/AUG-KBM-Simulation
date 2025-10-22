# plot the eigenfunction average parameters for 1d scan
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROscan/assist')
from cgyro_read_xj import *
from cgyro_ball import *
from cgyro_ball_class import OMFITcgyro_eigen
f = open(root['PLOTS']['CGYROscan']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
if setup['icgyro']==1:
    root['PLOTS']['CGYROscan']['assist']['collect.py'].run()
else:
    root['PLOTS']['CGYROscan']['assist']['collect_gyro.py'].run()

# then read all the data
nRange_x_eigen=len(para_x_eigen)
nRange_y_eigen=len(para_y_eigen)
num_ky_eigen=len(ky_eigen)
q_loc_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])
s_loc_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])
wd1_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])
k_par_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])
k_perp_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])
epar_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])
epar_esratio_ave=zeros([nRange_x_eigen,nRange_y_eigen,num_ky_eigen])

for k in range(nRange_x_eigen):
    for p in range(nRange_y_eigen):
        for q in range(num_ky_eigen):
            datanode=root['OUTPUTScan'][Para_x][num2str_xj(Range_x[k],effnum)][Para_y][num2str_xj(Range_y[p],effnum)]['lin'][num2str_xj(kyarr[q],effnum)]
            if setup['icgyro']==1:
                datanode=OMFITcgyro_eigen(datanode.filename)
            else:
                datanode = OMFITcgyro_eigen(datanode.filename)
            datanode.eigen_ave()
            q_loc_ave[k][p][q]=datanode.q_loc_ave
            s_loc_ave[k][p][q] = datanode.s_loc_ave
            wd1_ave[k][p][q] = datanode.wd1_ave
            k_par_ave[k][p][q] = datanode.k_par_ave
            k_perp_ave[k][p][q] = datanode.k_perp_ave
            epar_ave[k][p][q] = datanode.epar_ave
#            epar_esratio_ave[k][p][q] = datanode.epar_esratio_ave

########################################
# based on the profiles we get, then all the linear information can be plotted
Para_x_dic={0:Para_x, 1: Para_y}
Para_y_dic={0:Para_y, 1: Para_x}
bdry_dic=[0, 0, 1.e-4,1.e-4,1,1]
#avearr=['q_loc_ave','s_loc_ave','epar_ave','epar_esratio_ave','wd1_ave','k_par_ave','k_perp_ave']
avearr=['q_loc_ave','s_loc_ave','epar_ave','wd1_ave','k_par_ave','k_perp_ave']
#unitarr=['','','','','(c_s/a)','(a^{-1})','*rhos']
unitarr=['','','','(c_s/a)','(a^{-1})','*rhos']
if idimplt==0:
    Range_x_grid,Range_y_grid=meshgrid(para_x_eigen,para_y_eigen)
else:
    Range_x_grid,Range_y_grid=meshgrid(para_y_eigen,para_x_eigen)
for k in range(num_ky_eigen):
    for p in arange(1,len(avearr)+1):
        print(p)
        fig=figure(ky_eigen[k],figsize=[10,10])
        ax = fig.gca()
        axp=subplot(2,3,p)
        if idimplt==0:
            cmd="contourf(para_x_eigen,para_y_eigen,"+avearr[p-1]+".T[k],cmap='seismic')"
        else:
            cmd="contourf(para_y_eigen,para_x_eigen,"+avearr[p-1]+".T[k].T,cmap='seismic')"
        exec(cmd)
        colorbar()
        # try:
        if idimplt==0:
            cmd="cs=contour(Range_x_grid,Range_y_grid,"+avearr[p-1]+".T[k],["+str(bdry_dic[p-1])+"])"
        else:
            cmd="cs=contour(Range_y_grid,Range_x_grid,"+avearr[p-1]+".T[k].T, ["+str(bdry_dic[p-1])+"])"
        exec(cmd)

        try:
            pp = cs.collections[0].get_paths()[0]
            v = pp.vertices
            plot(v[:,0],v[:,1],'--r',linewidth=lw*2)
        except:
            print('not work at subplot '+str(p))
            continue
        axp.set_xscale(scale_dic[ilogx])
        axp.set_yscale(scale_dic[ilogy])
        if p>3 :
            xlabel(Para_x_dic[idimplt],fontsize=fs1,family='serif')
        if p==1 or p==4:
            ylabel(Para_y_dic[idimplt],fontsize=fs1,family='serif')
#        if p==2 or p==4:
#	    axp.set_zscale('log')
        xticks(fontsize=fs2,family='serif')
        yticks(fontsize=fs2,family='serif')
        title(avearr[p - 1] + unitarr[p - 1], fontsize=fs1, family='serif')
