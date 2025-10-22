# used to plot the time trace of several important quantities
import scipy
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    leny=len(y)
    lenw=len(w)
    yy=y[int(lenw/2):int(leny-lenw/2)+1]
    return yy


server='DIII-D'
treename=None
shot=174783
fig=figure(figsize=[10,11])
tagname=array(['wmhd','h_thh98y2','taue','cerqrott20','tinj','pinj','prmtan_neped','prmtan_peped','MPI1A322D'])
scalefac=array([1/1.e3,1,1,-1,1,1/1.e3,1/1.e19,1,1])
ylabelarr=array(['$W_{MHD}(kJ)$','$H_{98,y2}$','$\\tau_E(s)$','$V_{\phi,\\rho=0.3}(kms^{-1})$','$T_{inj}(N.m)$','$P_{inj}(MW)$','$n_{e,ped}(10^{19}m^{-3})$','$p_{e,ped}(kPa)$','$\\tilde B_{\\theta}(.a.u)$'])
fs1=24
fs2=20
fs3=16
lw=2
nsubplot=len(tagname)
#nsubplot=6
char_arr=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(H)','(I)']
for k in range(nsubplot):
    print(ylabelarr[k])
    ax=subplot(nsubplot,1,k+1)
    mdsstring=tagname[k]
    paraval=OMFITmdsValue(server=server,treename=treename,shot=shot,TDI=mdsstring).data()
    time = OMFITmdsValue(server=server,treename=treename,shot=shot,TDI=mdsstring).dim_of(0)
    plot(time/1.e3,smooth(paraval*scalefac[k],16),'-k',linewidth=lw)
    plot([2.1,2.1],[ax.get_ylim()[0],ax.get_ylim()[1]],'--b',linewidth=lw+1)
    plot([2.8,2.8],[ax.get_ylim()[0],ax.get_ylim()[1]],'--r',linewidth=lw+1)
#    ylabel(ylabelarr[k],fontsize=fs2)
    if k==0:
        text(2.1,3200,'SH-HI',ha='center', fontsize=fs2,family='serif',color='b')
        text(2.8,3200,'SH',ha='center', fontsize=fs2,family='serif',color='r')
        text(4.0,3200,'#174783',ha='center', fontsize=fs2,family='serif',color='k')
    text(3.0,1./4*ax.get_ylim()[-1],ylabelarr[k],fontsize=fs1)
    text(4.0,1./4*ax.get_ylim()[-1],char_arr[k],fontsize=fs1)
#    ylabel('$'+ylabelarr[k]+'$',fontsize=fs2)
    xticks(linspace(1,5,5),fontsize=fs1,family='serif')
    yticks(linspace(ax.get_ylim()[0],ax.get_ylim()[1],3),fontsize=fs1,family='serif')
    xlim([1.5,4.5])
#    if k==1:
#        ylim([1,2.5])
    if k==nsubplot-1:
        xlabel('Time(s)',fontsize=fs1,family='serif')
    if k!=nsubplot-1:
        xticks([])
    xticks(fontsize=fs2,family='serif')
    yticks(fontsize=fs2,family='serif')
# next step, add a line to the whole plot