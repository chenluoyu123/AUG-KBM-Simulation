import sys
sys.path.append('/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/CGYROalone/assist')
f = open(root['PLOTS']['CGYROalone']['assist']['getglobal.py'].filename, 'Ur')
for line in f:
    exec(line)
f.close()

# this is a class inherient from OMFITcgyro_base
# this script will mostly focus on handling the nonlinear CGYRO object
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
from classes.omfit_gacode import OMFITcgyro

class OMFITcgyro_nonlin(OMFITcgyro):
    """
    Class used for handling the nonlinear behaviors, functions included:
    get_flux_n, get_phi_n, get_quasiweight_n,
    """
    def __init__(self,filename=None,extra_files=[],test_mode=False,window=0.3,t_end=1):
        OMFITcgyro.__init__(self,filename,extra_files,test_mode)
        print('You have entered OMFITcgyro_nonlin!')
    	# get the background parameters
        inputcgyro=self['input.cgyro.gen']
        # geometric parameters
        self.Rmaj=inputcgyro['RMAJ']
        self.rmin=inputcgyro['RMIN']
        self.shift=inputcgyro['SHIFT']
        self.kappa=inputcgyro['KAPPA']
        self.skappa=inputcgyro['S_KAPPA']
        self.delta=inputcgyro['DELTA']
        self.sdelta=inputcgyro['S_DELTA']
        self.q=inputcgyro['Q']
        self.shear=inputcgyro['S']
        alpha=np.sum([inputcgyro['DENS_'+str(m)]*inputcgyro['TEMP_'+str(m)]*\
              (inputcgyro['DLNNDR_'+str(m)]+inputcgyro['DLNTDR_'+str(m)]) for m in np.arange(1,inputcgyro['N_SPECIES']+1) ])
        alpha=alpha*inputcgyro['BETAE_UNIT']*inputcgyro['Q']**2*inputcgyro['RMAJ']
        self.alpha=alpha*inputcgyro['BETA_STAR_SCALE']  # this parameter will be used
        self.betastar=self.alpha/self.q**2/self.Rmaj
        # resolution parameters
        self.n_theta=inputcgyro['N_THETA']
        self.ky=abs(self['kyrhos'])    #note that we will always use positive kyrhos, regardless of IPCCW&BTCCT
        self.n_n=self['n_n']
        if self.n_n>1:
            self.dky = abs(self.ky[1]-self.ky[0])
        else:
            self.dky=abs(self.ky[0])
        self.kx = self['kxrhos']
        self.n_r=self['n_r']
        self.dkx = abs(self.kx[1]-self.kx[0])
        self.t=self['t']
        self.n_time=self['n_time']
        # others
        self.field_tags=self['field_tags']
        self.n_field=self['n_field']
        self.species_tags=self['species_tags']
        self.n_species=self['n_species']
        # time averaging related
        self.n_theta_plot = self['theta_plot']
        if self.n_theta_plot==1:
            self.ftheta_plot=array([0.,0.])
        else:
            self.ftheta_plot=np.linspace(-np.pi,np.pi,self.n_theta_plot+1)
        self.window=window
        self.ind_t_ave=np.arange(int((t_end-window)*self.n_time),int(self.n_time*t_end))
        n_ind_temp=len(self.ind_t_ave)
        n_precies=10
        self.ind_t_ave = np.arange(int(self.n_time * t_end)-int(n_precies*ceil(n_ind_temp/n_precies)), int(self.n_time * t_end))  # sparse outline
        # also add a function on fluctuating issues, which is very time consuming
        self.getbigfield()
        self.phi_cmplx_t = self.kxky_phi[0, :, :, :, :] + 1j * self.kxky_phi[1, :, :, :,:]  # phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave]), includes the time dynamics
        if self.n_field>1:
            self.apar_cmplx_t = self.kxky_apar[0, :, :, :, :] + 1j * self.kxky_apar[1, :, :, :, :]
        elif self.n_field>2:
            self.bpar_cmplx_t = self.kxky_bpar[0, :, :, :, :] + 1j * self.kxky_bpar[1, :, :, :, :]



    def get_flux_t(self,i_field=0,i_species=0):
        """
        Usage: get_flux_n(self)
        functionality: get the time averaged flux for all channels versus time
        output: self.Gamma_t, Q_t, Pi_t, Gamma_t_ave, Q_t_ave, Pi_t_ave
        :return:
        """
        flux_t=self['flux_t']
        # use i_field<0 to do the sum
        if i_field<0:
            Gamma_t = flux_t['particle'].isel(species=i_species).sum(axis=0)
            Q_t = flux_t['energy'].isel(species=i_species).sum(axis=0)
            Pi_t = flux_t['momentum'].isel(species=i_species).sum(axis=0)
        else:
            Gamma_t = flux_t['particle'].isel(species=i_species,field=i_field)
            Q_t     = flux_t['energy'].isel(species=i_species, field=i_field)
            Pi_t    = flux_t['momentum'].isel(species=i_species, field=i_field)
        self.Gamma_t=Gamma_t
        self.Q_t=Q_t
        self.Pi_t=Pi_t
        self.Gamma_t_ave=self.Gamma_t[self.ind_t_ave].mean().data
        self.Q_t_ave = self.Q_t[self.ind_t_ave].mean().data
        self.Pi_t_ave = self.Pi_t[self.ind_t_ave].mean().data
        self.Gamma_t_std=self.Gamma_t[self.ind_t_ave].std().data
        self.Q_t_std = self.Q_t[self.ind_t_ave].std().data
        self.Pi_t_std = self.Pi_t[self.ind_t_ave].std().data

    def get_flux_n(self,i_field=0,i_species=0):
        """
        Usage: get_flux_n(self)
        functionality: get the time averaged flux for all kinetic channels versus n
        output: self.Gamma_n_ave,self.Q_n_ave, self_Pi_n_ave
        :return:
        """
        flux_n=self['flux_ky']
        if i_field<0:
            Gamma_n_ave=flux_n['particle'].isel(species=i_species,t=self.ind_t_ave).sum(axis=0).mean(axis=1).data
            Q_n_ave = flux_n['energy'].isel(species=i_species,  t=self.ind_t_ave).sum(axis=0).mean(axis=1).data
            Pi_n_ave = flux_n['momentum'].isel(species=i_species,  t=self.ind_t_ave).sum(axis=0).mean(axis=1).data
        else:
            Gamma_n_ave=flux_n['particle'].isel(species=i_species,field=i_field,t=self.ind_t_ave).mean(axis=1).data
            Q_n_ave = flux_n['energy'].isel(species=i_species, field=i_field, t=self.ind_t_ave).mean(axis=1).data
            Pi_n_ave = flux_n['momentum'].isel(species=i_species, field=i_field, t=self.ind_t_ave).mean(axis=1).data
        self.Gamma_n_ave=Gamma_n_ave/self.dky
        self.Q_n_ave=Q_n_ave/self.dky
        self.Pi_n_ave=Pi_n_ave/self.dky

    def get_quasiweight_n(self,i_field=0,i_species=0):
        """
        Usage: get_quasiweight_n(self,i_field.i_species)
        functionality: get the quasi-linear weight versu n, the quasilinear weight is flux/dky/phi^2
        output: self.quasiweight_Gamma,self.quasiweight_Q,self.quasiweight_Pi
        :return:
        """
        self.get_flux_n(i_field,i_species)
        self.get_phi_n(i_field, i_species)
        self.quasiweight_Gamma=self.Gamma_n_ave/self.dky/self.phi_n_ave
        self.quasiweight_Q = self.Q_n_ave /self.dky/ self.phi_n_ave
        self.quasiweight_Pi = self.Pi_n_ave /self.dky/ self.phi_n_ave

    def get_freq_n(self):
        """
        Usage: get_freq_n(self)
        functionality: get frequency  versus n
        output: self.freq_n
        :return:
        """
        self.freq_n=self['freq']['omega'].isel(t=self.ind_t_ave).mean(axis=1)


    # def get_phi_freq(self,i_theta_plot=self.n_theta_plot//2+1,i_r=self.n_r//2,i_n=3):
    def get_phi_freq(self,i_theta_plot=1,i_r=1,i_n=3,i_field=0):
        """
        Usage: get_phi(i_theta_plot=self.n_theta_plot//2+1,i_r=self.n_r//2,i_n=2)
        functionality: get phi plot on different dimensions
        phi_omega: for a given (i_theta_plot,i_r,i_n)   # 1D: frequency spectrum
        phi_kx_omega: for a given (i_theta_plot, i_n)   # 2D:over kx and freq
        phi_ky_omega: for a given (i_theta_plot)        # 2D: averaged over kx
        phi_ithetaplot_omega: for a given(i_n)          # 2D: kx averaged
        :return:
        """
        # self.getbigfield()
        t_ft=self.t[self.ind_t_ave]   # the time series for doing the fourier transform, ft is short for fourier transform
        n_t_ft=len(self.ind_t_ave)
        ## caution, don't write to be : phi_cmplx = casek.kxky_phi[0, :, :, :, case.t_ind_ave] + 1j * casek.kxky_phi[1, :, :, :, case.t_ind_ave]
        if i_field==0:
            phi_cmplx=self.kxky_phi[0,:,:,:,:]+1j*self.kxky_phi[1,:,:,:,:]
        elif i_field==1:
            phi_cmplx=self.kxky_apar[0,:,:,:,:]+1j*self.kxky_apar[1,:,:,:,:]
        else:
            phi_cmplx = self.kxky_bpar[0, :, :, :, :] + 1j * self.kxky_bpar[1, :, :, :, :]
        # get the frequency phi_omega
        omega_ft, phi_omega= OMFITcgyro_nonlin.fft_jian(self,t_ft,phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
        self.omega_ft=omega_ft  # the omega_ft is determined by the time series that we choose
        self.phi_omega=phi_omega
        # get the frequency phi_kx_omega
        phi_kx_omega = np.zeros([self.n_r, n_t_ft],dtype='complex')
        for i_rr in range(self.n_r):
            data = phi_cmplx[i_rr, i_theta_plot, i_n, self.ind_t_ave]
            omega, phi_kx_omega[i_rr, :] =  OMFITcgyro_nonlin.fft_jian(self,t_ft, data)
        kx_omega_grid, omega_kx_grid = np.meshgrid(self.kx,omega)
        self.kx_omega_grid = kx_omega_grid
        self.omega_kx_grid=omega_kx_grid
        self.phi_kx_omega=phi_kx_omega
#        contourf(kx_omega_grid, omega_kx_grid, abs(phi_kx_omega).T)
        # get the phi_ky_omega & phi_ky_omega_kx0
        i_r_kx0 = list(self.kx).index(0)  # find the index of kx=0
        phi_ky_omega = np.zeros([self.n_n, n_t_ft], dtype='complex')
        phi_ky_omega_kx0 = np.zeros([self.n_n, n_t_ft], dtype='complex')
        phi_kx_omegaa = np.zeros([self.n_r, n_t_ft], dtype='complex')
        for i_nn in np.arange(self.n_n):
            for i_rr in np.arange(self.n_r):
                data = phi_cmplx[i_rr, i_theta_plot, i_nn, self.ind_t_ave]
                omega, phi_kx_omegaa[i_rr, :] =  OMFITcgyro_nonlin.fft_jian(self,t_ft, data)
                phi_f = np.sum(phi_kx_omegaa, axis=0)
                phi_f_kx0 = phi_kx_omegaa[i_r_kx0,:]
                phi_ky_omega[i_nn] = phi_f  / max(abs(phi_f))  # will be normalized
                phi_ky_omega_kx0[i_nn] = phi_f_kx0 / max(abs(phi_f_kx0))  # will be normalized
        ky_omega_grid, omega_ky_grid = np.meshgrid(self.ky, omega)
        # contourf(ky_omega_grid, omega_ky_grid, abs(phi_ky_omega).T)  # can be compared to the linear/nonlinear frequency
        self.ky_omega_grid=ky_omega_grid
        self.omega_ky_grid=omega_ky_grid
        self.phi_ky_omega=phi_ky_omega
        self.phi_ky_omega_kx0 = phi_ky_omega_kx0
        # get the phi_ithetaplot_omega & phi_ithetaplot_omega_kx0 for a given i_n
        phi_ithetaplot_omega = np.zeros([self.n_theta_plot, n_t_ft],dtype='complex')
        phi_ithetaplot_omega_kx0 = np.zeros([self.n_theta_plot, n_t_ft], dtype='complex')
        for i_theta_plott in np.arange(self.n_theta_plot):
            for i_rr in np.arange(self.n_r):
                data = phi_cmplx[i_rr, i_theta_plott, i_n, self.ind_t_ave]
                omega, phi_kx_omegaa[i_rr, :] = OMFITcgyro_nonlin.fft_jian(self,t_ft, data)
                phi_f = sum(phi_kx_omegaa, axis=0)
                phi_f_kx0 = phi_kx_omegaa[i_r_kx0,:]
                phi_ithetaplot_omega[i_theta_plott] = phi_f #/max(abs(phi_f))  # will not be normalized to show the strength over poloidal plane
                phi_ithetaplot_omega_kx0[i_theta_plott] = phi_f_kx0  # /max(abs(phi_f))
        ftheta_plot = self.ftheta_plot[0:-1]
        ftheta_omega_grid, omega_ftheta_grid = np.meshgrid(ftheta_plot, omega)
        self.ftheta_omega_grid=ftheta_omega_grid
        self.omega_ftheta_grid=omega_ftheta_grid
        self.phi_ithetaplot_omega=phi_ithetaplot_omega
        self.phi_ithetaplot_omega_kx0 = phi_ithetaplot_omega_kx0
        # contourf(theta_omega_grid, omega_theta_grid, abs(phi_ithetaplot_omega).T)
        # get the omega_ky by averaged over phi from phi_ky_omega
        omega_n = np.zeros(self.n_n)
        for i_nn in np.arange(self.n_n):
            omega_n[i_nn] = np.sum(abs(self.phi_ky_omega[i_nn]) * self.omega_ft) / np.sum(abs(self.phi_ky_omega[i_nn]))
        self.omega_n=omega_n

    def get_kxky_phi(self, i_theta_plot=1):
        """
        Usage: get_kxky_phi(i_theta_plot=1)
        functionality: get the kxky_phi_avt and kx_phi_avetky for a given i_theta_plot
        :return:
        """
        # self.getbigfield()
        # phi_cmplx = self.kxky_phi[0, :, :, :, :] + 1j * self.kxky_phi[1, :, :, :,:]  # phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
        phi_cmplx=self.phi_cmplx_t[:,i_theta_plot,:,self.ind_t_ave]
        kxky_phi_avet=mean(abs(phi_cmplx),axis=-1)   #kxky_phi[i_r,i_theta_plot,i_n]
        kx_phi_avetky=mean(kxky_phi_avetky,axis=-1)
        self.kx_phi_avetky=kxky_phi_avetky
        self.kxky_phi_avet=kxky_phi_avet

    def get_kx_ky(self,i_field=0,i_theta_plot=1):
        """
        Usage: get_kx_ky(i_field=0,i_theta_plot=1)
        functionality: get the potential amplitude averaged kx
        :return:
        """
        # self.getbigfield()
        # phi_cmplx = self.kxky_phi[0, :, :, :, :] + 1j * self.kxky_phi[1, :, :, :,:]  # phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
        if i_field==0:
            phi_cmplx=self.phi_cmplx_t[:,i_theta_plot,:,self.ind_t_ave]
        elif i_field==1:
            phi_cmplx = self.apar_cmplx_t[:, i_theta_plot, :, self.ind_t_ave]
        else:
            phi_cmplx = self.bpar_cmplx_t[:, i_theta_plot, :, self.ind_t_ave]
        print(phi_cmplx.shape)
        phi_abs = mean(abs(phi_cmplx), axis=0)
        print(phi_abs.shape)
        kx_over_ky=np.zeros(self.n_n)
        kx0 = np.zeros(self.n_n)
        for i_n in range(self.n_n):
            kx0[i_n]=sum(self.kx*abs(phi_abs[:,i_n])**2)/sum(abs(phi_abs[:,i_n])**2)
            kx_rms=sum((self.kx-kx0[i_n])**2*abs(phi_abs[:,i_n])**2)/sum(abs(phi_abs[:,i_n])**2)
#            kx_rms=sum(abs(self.kx-kx0[i_n])**0.5*abs(phi_abs[:,i_n])**2)/sum(abs(phi_abs[:,i_n])**2)
            kx_rms=kx_rms**0.5
#            kx_rms=kx_rms**2
        self.kx0=kx0
        self.kx_over_ky=kx_rms/self.ky


    def get_phi_n(self,i_field=0,theta=0,i_n=1):
        """
        Usage: get_phi_n(self,i_field,theta=0)
        functionality: get the time averaged fluctuation amplitude versus n
        output: self.phi_n_ave
        This is abs(phi), not phi^2
        theta is in unit of phi ,should be in the range of [-1,1],
        theta=0 (-1/1) :the outboard/inboard midplane
        :return:
        """
        # self.getbigfield()
        moment='phi'
        fk,ftk = self.kxky_select(theta,i_field,moment,0) # ft[i_r,i_n,i_t]
        phi_ave_t=np.mean(abs(fk[:,:,self.ind_t_ave]),axis=-1)/self.rho  # for phi over ky
        phi_m_ave = np.sum(phi_ave_t, axis=-1)
        phi_n_ave = np.sum(phi_ave_t, axis=0)
        phi_n_overkx = phi_ave_t.T[i_n]
        self.phi_n_ave=phi_n_ave
        self.phi_m_ave = phi_m_ave
        self.phi_n_overkx=phi_n_overkx

    def get_npv_n(self,theta=0,i_n=1, i_species=-1, i_moment=0):
        """
        Usage: get_npv_n(self,i_field,theta=0)
        functionality: get the time averaged intensity of fluctuating density/pressure/velocity versus toroidal mode number
        i_moment: 0: den; 1: pressure; 2: v
        output: self.den_n_ave, self_p_ave, self_v_ave
        theta is in unit of phi ,should be in the range of [-1,1],
        theta=0 (-1/1) :the outboard/inboard midplane
        :return:
        """
        # self.getbigfield()
        moment_arr={0:'n', 1: 'e', 2:'v'}
        moment=moment_arr[i_moment]
        fk,ftk = self.kxky_select(theta,0,moment, i_species) # ft[i_r,i_n,i_t]
        npv_ave_t=np.mean(abs(fk[:,:,self.ind_t_ave]),axis=-1)/self.rho  # for phi over ky
        npv_m_ave = np.sum(npv_ave_t, axis=-1)
        npv_n_ave = np.sum(npv_ave_t, axis=0)
        npv_n_overkx = npv_ave_t.T[i_n]
        self.npv_n_ave=npv_n_ave
        self.npv_m_ave = npv_m_ave
        self.npv_n_overkx=npv_n_overkx



    def get_zfshear(self,i_theta_plot=1):
        """
        Usage: get_zfshear(i_theta_plot=1)
        functionality: get the zonal flow shearing rate at differnet poloidal places &
         zonal flow kx spectrum for a given poloidal space
        :return:
        """
        phi_cmplx=self.phi_cmplx_t#/self.rho
        zf_shear=np.zeros(self.n_theta_plot)
        for i_theta_plott in np.arange(self.n_theta_plot):
            zf_cmplx=phi_cmplx[:,i_theta_plott,0,self.ind_t_ave]
            zf_kx=np.mean(abs(zf_cmplx),axis=1)
            if i_theta_plott == i_theta_plot:
                self.zf_kx=zf_kx
                self.gamma_zf=self.kx**2*zf_kx/self.rho/self.dkx  # in unit of cs/a
            # zf_shear[i_theta_plott]=sum(self.kx**2*zf_kx*self.dkx)**0.5
            zf_shear[i_theta_plott]=sum(self.kx**2*zf_kx)/self.rho
        self.zf_shear=zf_shear  # sum gamma_zf over kx

    def get_shear_stress(self,theta=0):
        """
        :param theta:
        :return:
        the shear stress(S_stress) across ky_arr, shear address is defined in Candy-PPCF-2007, secion 3.2
        """
        moment='phi'
        i_field=0
        fk,ftk = self.kxky_select(theta,i_field,moment,0) # ft[i_r,i_n,i_t]
        phi_ave_t=np.mean(abs(fk[:,:,self.ind_t_ave]),axis=-1)
        phi_norm=phi_ave_t/self.dkx
        S_stress=np.zeros(self.n_n) # shear stress
        sign=-1;
        for k_n in np.arange(self.n_n):
            S_stress[k_n]=abs(sum((self.ky[k_n]**2+sign*self.kx**2)*abs(phi_norm[:,k_n])))*self.dkx
        self.S_stree=S_stress/self.rho  # has the unit of c_s/a

    def Energy_transfer_phi(self,i_theta_plot=1, kx_select=0.0, ky_select=0.2):
    #     do the kinetic energy transfer analysis for a given k_select=[kx_select,ky_select] from the whole [kx,ky] space
    #     self.getbigfield()
    #     phi_cmplxx = self.kxky_phi[0, :, :, :, :] + 1j * self.kxky_phi[1, :, :, :,:]  # phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
        phi_cmplx = self.phi_cmplx_t[:,i_theta_plot,:,:]
        # print(np.size(phi_cmplxx,axis=-2))
        # print(np.size(phi_cmplx, axis=-2))
        # find the the index of the kx_select and ky_select
        delta_k_err=0.01  # used to avoid unneccesary error
        if kx_select>np.max(self.kx)+delta_k_err or kx_select<np.min(self.kx)-delta_k_err:
            raise('Please choose the kx value in between'+str(self.kx[0])+' and '+str(self.kx[-1]))
        else:
            kx_select_incode=self.kx.flat[np.abs(self.kx-kx_select).argmin()]  # the real selected kx value in code
            ind_kx=np. where(self.kx==kx_select_incode)[0][0]
            print('The kxrhos that you really select is '+str(kx_select_incode))    # find the index of the kx_select_incode
        if ky_select>np.max(self.ky)+delta_k_err or ky_select < np.min(self.ky)-delta_k_err:
            raise('Please choose the ky value in between' + str(self.ky[0]) + ' and ' + str(self.ky[-1]))
        else:
            ky_select_incode = self.ky.flat[np.abs(self.ky - ky_select).argmin()]  # the real selected ky value in code
            ind_ky = np.where(self.ky == ky_select_incode)[0][0]
            print('The kyrhos that you really select is ' + str(ky_select_incode))  # find the index of the kx_select_incode
        # do the calculation
        len_ind_t_ave=len(self.ind_t_ave)
        S_k_kp=np.zeros([self.n_r,self.n_n],dtype=complex)
        S_k_kp_norm = np.zeros([self.n_r, self.n_n], dtype=complex)
        Lamda=np.zeros([self.n_r,self.n_n]) # the coupling coefficient
        ind_kx_0=np.where(self.kx == 0)
        for i_m in np.arange(self.n_r):
            for i_n in np.arange(self.n_n):
                kxp=self.kx[i_m]   # kx_prime
                kyp=self.ky[i_n]   # ky_prime
                kx_kxp=kx_select_incode-kxp  # kx-kx_prime
                ky_kyp=ky_select_incode-kyp  # ky-ky_prime
                Lamda[i_m,i_n]=1./2*kx_select_incode*kyp-kxp*ky_select_incode         # coupling coefficient
                Lamda[i_m,i_n]=Lamda[i_m,i_n]*((kxp**2+kyp**2)-(kx_kxp**2+ky_kyp**2))
                S_k_kp_temp=np.zeros([len_ind_t_ave],dtype=complex)
                S_k_kp_temp_norm = np.zeros([len_ind_t_ave],dtype=complex)
                for i_t in np.arange(len_ind_t_ave):
                    # print(ind_kx,ind_ky,i_t)
                    phi_k_select=phi_cmplx[ind_kx,ind_ky,self.ind_t_ave[i_t]]  # phi(k)
                    phi_k_p=phi_cmplx[i_m,i_n,self.ind_t_ave[i_t]]             # phi(k_p)
                    # for phi(k-kp)
                    if abs(ind_kx-i_m)>self.n_r//2-1 or abs(ind_ky-i_n)>self.n_n:
                        phi_k_kp=0                                              # phi(k-k_p)
                    else:
                        ind_kx_m= ind_kx-i_m                                    #ind_kx-i_m
                        ind_ky_kyp= 0 + (ind_ky-i_n)
                    # f(-n,-m)=conj(f(n,m)) and f(-n,m)=conj(f(n,-m))
                        if ind_ky_kyp<0:                                        # n < 0, conjuction is required
                            ind_kx_kxp = ind_kx_0 - ind_kx_m                    # ind of kx-kxp
                            phi_k_kp = np.conj(phi_cmplx[ind_kx_kxp, -1*ind_ky_kyp, self.ind_t_ave[i_t]])
                        else:
                            ind_kx_kxp=ind_kx_0 + ind_kx_m                      # ind of kx-kxp
                            phi_k_kp=phi_cmplx[ind_kx_kxp,ind_ky_kyp,self.ind_t_ave[i_t]]
                    S_k_kp_temp[i_t]=np.conj(phi_k_select)*phi_k_p*phi_k_kp
                    if abs(S_k_kp_temp[i_t])!=0:
                        S_k_kp_temp_norm[i_t] = S_k_kp_temp[i_t]/(abs(phi_k_select)*abs(phi_k_p)*abs(phi_k_kp))
                    else:
                        S_k_kp_temp_norm[i_t] = 0
                S_k_kp[i_m,i_n]=np.real(np.mean(S_k_kp_temp))                  # the real part is required for the bicoherence
                S_k_kp_norm[i_m, i_n] = np.real(np.mean(S_k_kp_temp_norm))  # the real part is required for the bicoherence
        S_k_kp=-1*self['input.cgyro.gen']['IPCCW']*S_k_kp                       #I believe its sign should be affected by IPCCW
        S_k_kp_norm = -1 * self['input.cgyro.gen']['IPCCW'] * S_k_kp_norm
        T_phi=Lamda*S_k_kp
        self.Lamda_phi=Lamda                                                        # coupling coefficient
        self.S_k_kp_phi=S_k_kp                                                      # Bicoherence,
        self.S_k_kp_norm_phi =S_k_kp_norm                                           # normalized Bicoherence
        self.T_phi=T_phi                                                        # energy transfer function

    def Energy_transfer_p(self,i_theta_plot=1, i_s=0, kx_select=0.0, ky_select=0.2):
    #     do the internal energy transfer analysis for a given k_select=[kx_select,ky_select] from the whole [kx,ky] space
    #     self.getbigfield()
    #     phi_cmplxx = self.kxky_phi[0, :, :, :, :] + 1j * self.kxky_phi[1, :, :, :,:]  # phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
        phi_cmplx = self.phi_cmplx_t[:,i_theta_plot,:,:]
        p_cmplxx = self.kxky_e[0, :, :, :, :, :] + 1j * self.kxky_e[1, :, :, :,:, :]  # phi_cmplx[i_r,i_theta_plot,i_n,self.ind_t_ave])
        p_cmplx = p_cmplxx[:,i_theta_plot, i_s, : ,:]  # i_s for i_species
        # find the the index of the kx_select and ky_select
        if kx_select>np.max(self.kx) or kx_select<np.min(self.kx):
            exit('Please choose the kx value in between'+str(self.kx[0])+' and '+str(self.kx[-1]))
        else:
            kx_select_incode=self.kx.flat[np.abs(self.kx-kx_select).argmin()]  # the real selected kx value in code
            ind_kx=np. where(self.kx==kx_select_incode)[0][0]
            print('The kxrhos that you really select is '+str(kx_select_incode))    # find the index of the kx_select_incode
        if ky_select>np.max(self.ky) or ky_select < np.min(self.ky):
            exit('Please choose the ky value in between' + str(self.ky[0]) + ' and ' + str(self.ky[-1]))
        else:
            ky_select_incode = self.ky.flat[np.abs(self.ky - ky_select).argmin()]  # the real selected ky value in code
            ind_ky = np.where(self.ky == ky_select_incode)[0][0]
            print('The kyrhos that you really select is ' + str(ky_select_incode))  # find the index of the kx_select_incode
        # do the calculation
        len_ind_t_ave=len(self.ind_t_ave)
        S_k_kp=np.zeros([self.n_r,self.n_n],dtype=complex)
        S_k_kp_norm = np.zeros([self.n_r, self.n_n], dtype=complex)
        Lamda=np.zeros([self.n_r,self.n_n]) # the coupling coefficient
        ind_kx_0=np.where(self.kx == 0)
        for i_m in np.arange(self.n_r):
            for i_n in np.arange(self.n_n):
                kxp=self.kx[i_m]   # kx_prime
                kyp=self.ky[i_n]   # ky_prime
                # kx_kxp=kx_select_incode-kxp  # kx-kx_prime
                # ky_kyp=ky_select_incode-kyp  # ky-ky_prime
                Lamda[i_m,i_n]=kx_select_incode*kyp-kxp*ky_select_incode         # coupling coefficient
                # Lamda[i_m,i_n]=Lamda[i_m,i_n]*((kxp**2+kyp**2)-(kx_kxp**2+ky_kyp**2))
                S_k_kp_temp=np.zeros([len_ind_t_ave],dtype=complex)
                S_k_kp_temp_norm = np.zeros([len_ind_t_ave],dtype=complex)
                for i_t in np.arange(len_ind_t_ave):
                    # print(ind_kx,ind_ky,i_t)
                    p_k_select=p_cmplx[ind_kx,ind_ky,self.ind_t_ave[i_t]]  # p(k)
                    p_k_p=p_cmplx[i_m,i_n,self.ind_t_ave[i_t]]             # p(k_p)
                    # for phi(k-kp)
                    if abs(ind_kx-i_m)>self.n_r//2-1 or abs(ind_ky-i_n)>self.n_n:
                        phi_k_kp=0                                              # phi(k-k_p)
                    else:
                        ind_kx_m= ind_kx-i_m                                    #ind_kx-i_m
                        ind_ky_kyp= 0 + (ind_ky-i_n)
                    # f(-n,-m)=conj(f(n,m)) and f(-n,m)=conj(f(n,-m))
                        if ind_ky_kyp<0:                                        # n < 0, conjuction is required
                            ind_kx_kxp = ind_kx_0 - ind_kx_m                    # ind of kx-kxp
                            phi_k_kp = np.conj(phi_cmplx[ind_kx_kxp, -1*ind_ky_kyp, self.ind_t_ave[i_t]])
                        else:
                            ind_kx_kxp=ind_kx_0 + ind_kx_m                      # ind of kx-kxp
                            phi_k_kp=phi_cmplx[ind_kx_kxp,ind_ky_kyp,self.ind_t_ave[i_t]]
                    S_k_kp_temp[i_t]=np.conj(p_k_select)*p_k_p*phi_k_kp
                    if abs(S_k_kp_temp[i_t])!=0:
                        S_k_kp_temp_norm[i_t] = S_k_kp_temp[i_t]/(abs(p_k_select)*abs(p_k_p)*abs(phi_k_kp))
                    else:
                        S_k_kp_temp_norm[i_t] = 0
                S_k_kp[i_m,i_n]=np.real(np.mean(S_k_kp_temp))                  # the real part is required for the bicoherence
                S_k_kp_norm[i_m, i_n] = np.real(np.mean(S_k_kp_temp_norm))  # the real part is required for the bicoherence
        S_k_kp=-1*self['input.cgyro.gen']['IPCCW']*S_k_kp                       #I believe its sign should be affected by IPCCW
        S_k_kp_norm = -1 * self['input.cgyro.gen']['IPCCW'] * S_k_kp_norm
        T_p=Lamda*S_k_kp
        self.Lamda_p=Lamda                                                        # coupling coefficient
        self.S_k_kp_p=S_k_kp                                                      # Bicoherence,
        self.S_k_kp_norm_p =S_k_kp_norm                                           # normalized Bicoherence
        self.T_p=T_p                                                        # energy transfer function


# geometry
    def miller_wd_s(self,theta_p=np.linspace(-np.pi,np.pi,37)):
        """
        Usage: miller_wd_s(self)
        Functionality: get the self.wd1, wd2,wd3, Rs,Zs,Gq,rc_theta,l, BoverBunit,Gtheta, gcos1, gcos2, gsin, gradr, k_perp
        output: self.q_loc, self.s_loc as a function of theta_p
        :return:
        """
        n_theta_p=len(theta_p)
        xd = np.arcsin(self.delta)
        arg = theta_p + xd * np.sin(theta_p)
        Rs = self.Rmaj + self.rmin * np.cos(arg)  # R position
        Zs = self.kappa * self.rmin * np.sin(theta_p)  # Z Position
        #  calculate the Jocobian
        dRdtheta = -1 * self.rmin * np.sin(arg) * (1 + np.cos(theta_p) * xd)
        dZdtheta = self.kappa * self.rmin * np.cos(theta_p)
        dldtheta = (dRdtheta ** 2 + dZdtheta ** 2) ** 0.5
        dRdr = self.shift + np.cos(arg) - np.sin(theta_p) * np.sin(arg) * self.sdelta
        dZdr = self.kappa * np.sin(theta_p) * (1 + self.skappa)
        det = Rs * (dRdr * dZdtheta - dRdtheta * dZdr)
        # look at grad r
        gradr = dldtheta * Rs / det
        l = integrate.cumtrapz(dldtheta, theta_p, initial=0)
        d2Zdtheta2 = np.gradient(dZdtheta) / np.gradient(theta_p)
        d2Rdtheta2 = np.gradient(dRdtheta) / np.gradient(theta_p)
        rc_theta = dldtheta ** 3 / (dRdtheta * d2Zdtheta2 - dZdtheta * d2Rdtheta2)
        #        IoverBunit=2*np.pi*self.rmin/integrate.trapz(1/Rs/gradr,l) # scale
        IoverBunit = 2 * np.pi * self.rmin / np.trapz(1 / Rs / gradr, l)  # scale
        BtoverBunit = IoverBunit / Rs
        BpoverBunit = self.rmin / Rs * gradr / self.q
        BoverBunit = (BtoverBunit ** 2 + BpoverBunit ** 2) ** 0.5
        #    geometric components
        cosu = np.gradient(Zs) / np.gradient(l)  # dZ/dl
        sinu = -np.gradient(Rs) / np.gradient(l)  # - dR/dl
        gsin = BtoverBunit / BoverBunit * self.Rmaj / BoverBunit * np.gradient(BoverBunit) / np.gradient(l)
        gcos1 = (BtoverBunit / BoverBunit) ** 2 * self.Rmaj / Rs * cosu + (
                    BpoverBunit / BoverBunit) ** 2 * self.Rmaj / rc_theta
        gcos2 = -1. / 2. / BoverBunit ** 2 * self.Rmaj * gradr * self.betastar
        # E series for mu
        E1kernel = 2. / Rs / gradr * BtoverBunit / BpoverBunit * (self.rmin / rc_theta - self.rmin / Rs * cosu)
        E2kernel = 1. / Rs / gradr * (BoverBunit / BpoverBunit) ** 2
        E3kernel = 1. / 2. / Rs * BtoverBunit / BpoverBunit / BpoverBunit ** 2
        # chang the order[0~2*pi], denoted new
        E1kernel_new = OMFITcgyro_nonlin.changeorder(self, E1kernel)
        E2kernel_new = OMFITcgyro_nonlin.changeorder(self, E2kernel)
        E3kernel_new = OMFITcgyro_nonlin.changeorder(self, E3kernel)
        ntheta_half = np.int(np.round((n_theta_p + 1) / 2))
        l_new = np.zeros(n_theta_p)
        l_new[0:ntheta_half - 1] = l[ntheta_half - 1:n_theta_p - 1] - l[ntheta_half - 1]
        l_new[ntheta_half - 1:n_theta_p - 1] = l[0:ntheta_half - 1] + l[ntheta_half - 1]
        E1_new = integrate.cumtrapz(E1kernel_new, l_new, initial=0)
        E2_new = integrate.cumtrapz(E2kernel_new, l_new, initial=0)
        E3_new = integrate.cumtrapz(E3kernel_new, l_new, initial=0)  # in the order of 0~2*pi
        # change back to [-pi,pi]
        E1 = OMFITcgyro_nonlin.changeorder(self, E1_new)
        E2 = OMFITcgyro_nonlin.changeorder(self, E2_new)
        E3 = OMFITcgyro_nonlin.changeorder(self, E3_new)
        E1[0:ntheta_half - 1] = E1[0:ntheta_half - 1] - (E1[ntheta_half - 2] + E1[ntheta_half])
        E2[0:ntheta_half - 1] = E2[0:ntheta_half - 1] - (E2[ntheta_half - 2] + E2[ntheta_half])
        E3[0:ntheta_half - 1] = E3[0:ntheta_half - 1] - (E3[ntheta_half - 2] + E3[ntheta_half])
        fstar = 1 / E2_new[-1] * (
                    2 * np.pi * self.q * self.shear / self.rmin - 1 / self.rmin * E1_new[-1] + self.betastar *
                    E3_new[-1])
        THETA = Rs * BpoverBunit / BoverBunit * abs(gradr) * (1 / self.rmin * E1 + fstar * E2 - self.betastar * E3)
        Gq = 1 / self.q * (self.rmin / Rs * BoverBunit / BpoverBunit)
        Gtheta = BoverBunit * Rs / self.Rmaj / self.rmin / gradr * dldtheta
        # assuming partial/partial(theta_p)=-i k_theta, partial/partial(r)=-i k_r
        ktheta = self.ky[1]  # the ktheta is the kyrho_s specified in input.cgyro
        kr = 0
        vpar2v = 1 / 2  # assuming an isotropic distribution
        #  the sign here is consistent with the s-alpha geometry
        # the Rmaj exist in the denominator so that the output drift frequency has the unit of c_s/a,consistent with cgyro units
        wd1 = ktheta * Gq * 1 * (gcos1 + gcos2 + THETA * gsin) / self.Rmaj
        wd2 = -1 * ktheta * Gq * vpar2v * gcos2 / self.Rmaj
        wd3 = -1 * kr * 1 * gradr * gsin / self.Rmaj
        k_perp = (ktheta * Gq * THETA) ** 2 + (ktheta * Gq) ** 2
        k_perp = k_perp ** 0.5
        #  the local q and magnetic shear
        IoverBp = IoverBunit / BpoverBunit
        q_loc = IoverBp / Rs ** 2 * np.gradient(l) / np.gradient(theta_p)
        # solute for I' with given s
        D0_kernel = 1. / Rs * (2. / rc_theta / Rs - 2 * cosu / Rs ** 2) * IoverBp
        D1_kernel_part = 1. / Rs ** 2. * (BoverBunit / BpoverBunit) ** 2  # the dI/dr is unknow yet
        D2_kernel = -1. / 2 * 1. / Rs ** 2. * IoverBp * self.betastar / BpoverBunit ** 2
        D0 = integrate.cumtrapz(D0_kernel, l, initial=0)
        D1_part = integrate.cumtrapz(D1_kernel_part, l, initial=0)
        D2 = integrate.cumtrapz(D2_kernel, l, initial=0)
        dqdr = self.shear * self.q / self.rmin
        dIdr = (2 * np.pi * dqdr - D0[-1] - D2[-1]) / D1_part[-1]
        D1_kernel = 1. / Rs ** 2. * (BoverBunit / BpoverBunit) ** 2 * dIdr
        D1 = integrate.cumtrapz(D1_kernel, l, initial=0)
        # mu1=D0+D1+D2
        s_loc = self.rmin / q_loc * (D0_kernel + D1_kernel + D2_kernel) * np.gradient(l) / np.gradient(theta_p)
        #     get the output
        self.wd1 = wd1
        self.wd2 = wd2
        self.wd3 = wd3
        self.k_perp = k_perp
        self.Rs = Rs
        self.Zs = Zs
        self.l = l
        self.rc_theta = rc_theta
        self.BoverBunit = BoverBunit
        self.Gq = Gq
        self.Gtheta = Gtheta
        self.gcos1 = gcos1
        self.gcos2 = gcos2
        self.gsin = gsin
        self.gradr = gradr
        self.THETA = THETA
        self.q_loc = q_loc
        self.s_loc = s_loc

    def changeorder(self, arr):
        # called by miller_drffreq
        n_arr = len(arr)
        arr_new = np.zeros(n_arr)
        n_arr_half = np.int(np.round((n_arr + 1) / 2))
        arr_new[0:n_arr_half - 1] = arr[n_arr_half - 1:n_arr - 1]
        #    arr_new[n_arr_half-1:n_arr-1]=arr[0:n_arr_half-1];
        arr_new[n_arr_half - 1:n_arr] = arr[0:n_arr_half]
        return arr_new

    # some assist functions
    def fft_jian(self,tt,yt):
        """
        Usage: fft_jian(time.data)
        functionality: a good fourier transport of a complex series which is able to resolve the frequency sign
        :param data:
        :return:
        """
        ndata=len(yt)
        dt = tt[1]-tt[0]
        yw = np.fft.fft(yt)
        yw_shifted = np.fft.fftshift(yw)
        w=2.*np.pi*linspace(-1./(2.*dt),1/(2.*dt),ndata) # in unit of rad/s, not
        yw=yw_shifted/ndata
        return w, yw

# load cases that required to be plotted
t_ave=root['SETTINGS']['PLOTS']['nl']['t_ave']
t_end=root['SETTINGS']['PLOTS']['nl']['t_end'] # the t_end
case_plot=root['SETTINGS']['PLOTS']['nl']['case_plot']
ncase=len(case_plot)
ncount=0
for case in case_plot:
    if isinstance(t_ave,float):
        root['OUTPUTS'][case]=OMFITcgyro_nonlin(root['OUTPUTS'][case].filename,window=t_ave,t_end=t_end)
    else:
        root['OUTPUTS'][case] = OMFITcgyro_nonlin(root['OUTPUTS'][case].filename, window=t_ave[ncount], t_end=t_end[ncount])
    ncount=ncount+1

