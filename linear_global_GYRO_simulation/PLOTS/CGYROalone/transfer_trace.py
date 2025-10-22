# trace back the energy from a given (kx_end,ky_end), which is linearly unstable but has finite amplitute nonlinearly
# the goal is to understand the (kx_start,ky_start), which should have kx_start=0
#  (kx_start, ky_start) should be linearly unstable and of course nonlinearly can contribute flux as well
import sys
from matplotlib import ticker, cm
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROalone/assist')
f = open(root['PLOTS']['CGYROalone']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()
icgyro=root['SETTINGS']['SETUP']['icgyro']
if icgyro==1:
    root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
else:
    root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
# check to see how many fields that all the cases have in common
n_field=3
for k_case in range(n_case):
    casek = outputs[case_plot[k_case]]
    if icgyro==1:
        n_field=min(len(casek['field_tags']),n_field)
    else:
        n_field=min(len(casek['tagfieldtext']),n_field)

def find_min_index(arr):
    """
    此函数用于寻找二维数组中最小值的索引
    :param arr: 输入的二维数组
    :return: 最小值的索引，格式为 (行索引, 列索引)
    """
    # 转换为 numpy 数组
    arr = np.array(arr)
    # 找到最小值的一维索引
    flat_index = np.argmin(arr)
    # 将一维索引转换为二维索引
    rows, columns = arr.shape
    row_index = flat_index // columns
    col_index = flat_index % columns
    return row_index, col_index

def pltcase(T_phi_ave,kx_grid,ky_grid,kx_select_incode,ky_select_incode,ky_bdry):
    fs2=24
    figure(figsize=[8,8])
    T_phi_absmax=amax(abs(T_phi_ave))
    contourf(kx_grid,ky_grid,T_phi_ave.T,levels=linspace(-1*T_phi_absmax,T_phi_absmax,nlevel),cmap='seismic')    
    colorbar()
    plot(kx_grid,ky_grid,'o',color='b',markersize=ms)
    plot(kx_select_incode,ky_select_incode,'o',markersize=ms*6,color='k')
    title('Coupling Coefficient # $k_x\\rho_s$=%.3f $k_y\\rho_s$=%.3f' % (kx_select_incode,ky_select_incode) , fontsize=fs2, family='serif')
    xlim([-1*ky_bdry,ky_bdry])
    xticks(fontsize=fs2)
    yticks(fontsize=fs2)

plot_nl=root['SETTINGS']['PLOTS']['nl']
t_ave_orig=plot_nl['t_ave']
t_end_orig=plot_nl['t_end']
n_split=1
iplt=0 # choose to whether to plot in the iteration or not
ms=2
ky_bdry=0.2
nlevel=128
n_count_max=20
kx_end_arr=array([0])
ky_end_arr=linspace(0.3,0.84,10)
#ky_end_arr=array([0.48])
kx_select_arr_overall=zeros([len(kx_end_arr),len(ky_end_arr),n_count_max+2])
ky_select_arr_overall=zeros([len(kx_end_arr),len(ky_end_arr),n_count_max+2])
icount_kx=-1
for kx_end in kx_end_arr:
    icount_kx=icount_kx+1
    icount_ky=-1
    for ky_end in ky_end_arr:
        icount_ky=icount_ky+1
        iexit_status=0
        # the (kx_end, ky_end) is the k-pair that we want to trace its energy origin
        # initialize the kx_select and ky_select
        kx_select=kx_end
        ky_select=ky_end
        print('Starting, kx_select=%.3f,ky_select=%.3f' %(kx_select,ky_select))
        # note that it could be better if we split the time window in to several sub time windows
        kx_select_arr=zeros(n_count_max+2)
        ky_select_arr=zeros(n_count_max+2)
        kx_select_arr[0]=kx_end
        ky_select_arr[0]=ky_end
        # # plot Energy transfer of phi (kinetic energy)
        for k_case in range(n_case):
            casek = outputs[case_plot[k_case]]
            t_ave = t_ave_orig / n_split
            print(case_plot[k_case])
            if icgyro==1:
                i_theta_plot=casek.theta_plot//2
            else:
                i_theta_plot=casek.n_theta_plot//2
            T_phi_split = zeros([n_split, casek.n_r-1, casek.n_n])
            for i_split in arange(n_split):
                t_end=t_end_orig-t_ave*i_split
                plot_nl['t_ave']=t_ave
                plot_nl['t_end'] = t_end
                if icgyro==1:
                    root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
                else:
                    root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
                casek.Energy_transfer_phi_improve(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
                T_phi_split[i_split] = casek.T_phi
            T_phi_ave = mean(T_phi_split, axis=0)
            kx_grid, ky_grid = meshgrid(casek.kx, casek.ky)
            if iplt==1:
                pltcase(T_phi_ave,kx_grid,ky_grid,casek.kx_select_incode,casek.ky_select_incode,ky_bdry)
            row_index,col_index=find_min_index(T_phi_ave[:,1:])  # do note search the zonal one
            icount=1
            while abs(casek.kx[row_index])>1.e-3:
                kx_select=casek.kx[row_index]
                ky_select=casek.ky[col_index+1]
                kx_select_arr[icount]=kx_select
                ky_select_arr[icount]=ky_select
                print('Going On, kx_select=%.3f,ky_select=%.3f' %(kx_select,ky_select))
                T_phi_split = zeros([n_split, casek.n_r-1, casek.n_n])
                for i_split in arange(n_split):
                    t_end=t_end_orig-t_ave*i_split
                    plot_nl['t_ave']=t_ave
                    plot_nl['t_end'] = t_end
                    if icgyro==1:
                        root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
                    else:
                        root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
                    casek.Energy_transfer_phi_improve(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
                    T_phi_split[i_split] = casek.T_phi
                T_phi_ave = mean(T_phi_split, axis=0)
                row_index,col_index=find_min_index(T_phi_ave[:,1:])  # do note search the zonal one
                if iplt==1:
                    pltcase(T_phi_ave,kx_grid,ky_grid,casek.kx_select_incode,casek.ky_select_incode,ky_bdry)
                if icount>n_count_max:  # if after 20 iterations, still can not converge, then exist
                    iexit_status=1
                    break
                icount=icount+1
        # plot the last case
            kx_select=casek.kx[row_index]
            ky_select=casek.ky[col_index+1]
            kx_select_arr[icount]=kx_select
            ky_select_arr[icount]=ky_select
            print('Last case, kx_select=%.3f,ky_select=%.3f' %(kx_select,ky_select))
            T_phi_split = zeros([n_split, casek.n_r-1, casek.n_n])
            for i_split in arange(n_split):
                t_end=t_end_orig-t_ave*i_split
                plot_nl['t_ave']=t_ave
                plot_nl['t_end'] = t_end
                if icgyro==1:
                    root['PLOTS']['CGYROalone']['assist']['collect.py'].run()
                else:
                    root['PLOTS']['CGYROalone']['assist']['collect_gyro.py'].run()
                casek.Energy_transfer_phi_improve(i_theta_plot=i_theta_plot, kx_select=kx_select, ky_select=ky_select)
                T_phi_split[i_split] = casek.T_phi
            T_phi_ave = mean(T_phi_split, axis=0)
            row_index,col_index=find_min_index(T_phi_ave[:,1:])  # do note search the zonal one
            if iplt==1:
                pltcase(T_phi_ave,kx_grid,ky_grid,casek.kx_select_incode,casek.ky_select_incode,ky_bdry)
        if iexit_status==1:
            print('reached maximum interations '+str(n_count_max))
        else:
            print('Congratulations, convergence is reached!')
        kx_select_arr_overall[icount_kx,icount_ky]=kx_select_arr
        ky_select_arr_overall[icount_kx,icount_ky]=ky_select_arr
#
# return the original parameters
plot_nl['t_ave']=t_ave_orig
plot_nl['t_end']=t_end_orig