import sys
sys.path.append('/home/vinu/scripts/BiasEstimator/bias')
import os
import numpy as np
import pylab as pl
import cosmolopy.distance as cd
import pyfits
from scipy import signal, optimize
import scipy.ndimage as nd
from convolution_mask import convolve_mask_fft, Gaussian
import kappa_utils as ku
import config as c
#import minuit

#from mayavi import mlab



class WeightedGalKappa:

    def __init__(self, ra, dec, z, opath, smooth, pixel_scale,
                 bin_ra, bin_dec, bin_z, mask,ipath,sourcefile,lensfile,
                 zs=0.8,pdf_zs=None, zmin_s=0.4, zmax_s=1.1,
                 zmin_l=0.1, zmax_l=1.1,rho_weight=None):

        self.sourcefile = os.path.join(ipath,sourcefile)
        self.lensfile = os.path.join(ipath,lensfile)
        self.smooth = smooth
        self.pixel_scale = pixel_scale
        self.bin_ra = bin_ra
        self.bin_dec = bin_dec
        self.bin_z = bin_z
        self.mask = mask

        self.zs = zs
        self.pdf_zs = pdf_zs
        self.zmin_l = zmin_l
        self.zmax_l = zmax_l
        self.zmin_s = zmin_s
        self.zmax_s = zmax_s
        self.initialize(ra, dec, z)
        #self.delta_rho_3d(bin_ra, bin_dec, bin_z)

    
    def initialize(self, ra, dec, z):
        self.cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 
                      'omega_k_0':0.0, 'h':0.72}

        self.z = z 
        if self.zmin_l is None:
            con = np.ones(self.z.shape).astype('bool')
        else: 
            con = (self.z >= self.zmin_l) & (self.z <= self.zmax_l)
        self.ra = ra
        self.dec = dec
        self.z = z

    def return_size(self, x, s=3):    
        """Return size of Gaussina kernal"""
        if np.ceil(2*s*x) % 2 == 0:
            size = np.ceil(2*s*x) + 1.
        else:
            size = np.ceil(2*s*x)
        return size
 
    def delta_rho_3d(self):

        if self.smooth == 0:
            self.sigma = 0.0
            self.kern_size = 1
            self.g_2d = np.array([[1]]) 
            self.g_3d = np.array([[[1]]]) 
        else:
            self.sigma = self.smooth/self.pixel_scale 
            self.kern_size = self.return_size(self.sigma, s=3)
            self.g_2d = Gaussian(self.sigma, size=self.kern_size, ndim=2)
            self.g_3d = Gaussian(self.sigma, size=self.kern_size, ndim=3)

        print 'Pix scale %2.2f arcmin'%self.pixel_scale
        print 'Sigma %2.2f pixels'%self.sigma

        self.N3d, edges = np.histogramdd(np.array([self.z, self.dec, 
                          self.ra]).T,
                          bins=(self.bin_z, self.bin_dec, self.bin_ra))

        self.raedges = edges[2]
        self.decedges = edges[1]
        self.zedges = edges[0]
        self.zavg = (self.zedges[:-1] + self.zedges[1:]) / 2.
        self.raavg = (self.raedges[:-1] + self.raedges[1:]) / 2.
        self.decavg = (self.decedges[:-1] + self.decedges[1:]) / 2.

        # The total galaxies per redshift slice
        N1d, zedge = np.histogram(self.z, bins=self.zedges) 

        # Average per redshift slices. Dividing it by the number of 
        # pixels in RA and DEC directions 
        self.n1d = N1d / (self.N3d.shape[1] * self.N3d.shape[2] * 1.0) 
        self.n1d[self.n1d == 0] = 1
 
        #print bin_ra, bin_dec, bin_z, self.N3d.shape, self.n1d.shape

        # subtracting the average number from each redshift slice.
        self.n3d = self.N3d - self.n1d[:,np.newaxis][:,np.newaxis] 

        #print self.n1d , N1d, self.edges[0]

        #delta = (rho - rho_m) / rho_m
        self.delta3d = self.n3d / self.n1d[:,np.newaxis][:,np.newaxis] 

        #print self.delta3d 

        np.savez(c.bigfilepath+'density.npz', raedge=self.raavg, decedge=self.decavg, 
                 zedge=self.zavg, N3d=self.N3d, n1d=self.n1d, 
                 n3d=self.n3d, delta3d=self.delta3d)

    def comoving_d(self):
        self.cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 
                      'omega_k_0':0.0, 'h':0.72}
        comoving_edges = cd.comoving_distance(self.zedges, **self.cosmo) #Mpc
        #comoving_edges /= (1. + self.zedges)
        self.d_c = cd.comoving_distance(self.zavg, **self.cosmo) #Mpc

        #There is some subtilities in this case. When using MICE, the answer
        #makes sense when commenting the following line. When using BCC
        #it requires the following line
        #self.d_c /= (1. + self.zavg)

        #self.d_s = comoving_edges[-1] #source distance
        self.d_s = cd.comoving_distance(self.zs, **self.cosmo) #source distance
        #self.d_s /= (1. + self.zs)
        self.delta_d = comoving_edges[1:] - comoving_edges[:-1]
        
        self.a = 1 / (1 + self.zavg)

    def kappa_predicted(self):
        self.comoving_d()
        c_light = 3e5 #km/s
 
        # Eq. 9 Amara et al.
        constant = ((100. * self.cosmo['h'])**2 * self.cosmo['omega_M_0']) * \
                   (3/2.) * (1/c_light**2)         

        if type(self.zs) is np.ndarray:
            if self.pdf_zs is None:
                self.pdf_zs = np.arange(len(self.d_c)*len(self.d_s)).reshape((len(self.d_c),len(self.d_s)))*0.0 + 1.0 #DEFAULT Flat Distribution
            else:
                self.pdf_zs = np.resize(self.pdf_zs,(len(self.d_c),len(self.d_s)))

            self.pdf_zs /= np.linalg.norm(self.pdf_zs[0,:],ord=1)#normalize probabilities to be used in integral
            self.pdf_zs = np.transpose(self.pdf_zs)
            
            twod_d_s = np.transpose(np.resize(self.d_s,(len(self.d_c),len(self.d_s))))
            twod_d_c = np.resize(self.d_c,(len(self.d_s),len(self.d_c)))
            
            integral_2 = (self.pdf_zs*(twod_d_s - twod_d_c) / twod_d_s)
            integral_2_summed = [integral_2[:,x].sum() for x in range(len(self.d_c))]

            integral_1 = ((self.d_c * integral_2_summed) * \
                          (self.delta_d / self.a))[:,np.newaxis][:,np.newaxis]
        else:
            integral_1 = ((self.d_c * (self.d_s - self.d_c) / self.d_s) * \
                          (self.delta_d / self.a))[:,np.newaxis][:,np.newaxis]#NOW 3D


        # Smooth the 3d density field and find kappa from that
        self.mask_3d = np.ones(self.delta3d.shape) * self.mask
        xxx, self.delta3d_sm, yyy = convolve_mask_fft(self.delta3d, \
                                        self.mask_3d, self.g_3d, ignore=0.0)
        self.kappa_pred_3d = constant * np.sum(integral_1 * self.delta3d_sm, \
                                               axis=0)

        # Use unsmoothed density field and generate kappa from that. Later
        # smooth the 2D kappa field
        self.kappa_pred = constant * np.sum(integral_1 * self.delta3d, axis=0)
        xxx, self.kappa_pred, yyy = convolve_mask_fft(self.kappa_pred, self.mask, \
                                        self.g_2d, ignore=0.0) 

        print integral_1.shape, self.delta3d.shape, self.kappa_pred.shape

        np.savez(c.opath+'kappa_predicted.npz', kappa=self.kappa_pred, \
                  kappa3d=self.kappa_pred_3d)


    def gamma_predicted(self):
        """Eq. 26 of Schenider"""
        @np.vectorize 
        def D_kernel(Dx, Dy, Dsq):
            if abs(Dsq)==0:
                return 0., 0.
            else:
                return (Dy**2 - Dx**2) / Dsq**2., (-2 * Dx * Dy) / Dsq**2.
        Dx, Dy = np.mgrid[-10:10:21j, -10:10:21j]
        Dsq = Dx**2 + Dy**2
        #D1 = (Dy**2 - Dx**2) / Dsq**2.
        #D2 = (-2 * Dx * Dy) / Dsq**2.
        D1, D2 = D_kernel(Dx, Dy, Dsq)
        D = D1 + 1j * D2

        #D = -1. / (Dx - 1j * Dy)**2.

        self.gamma_p = signal.convolve2d(self.kappa_pred, \
                                         D, mode='same') / np.pi
        self.gamma_tp = signal.convolve2d(self.kappa_true, \
                                         D, mode='same') / np.pi
   

    def plot():
        a = np.random.randint(0, 100, size=100)
        b = np.random.randint(0, 10, size=100)
        c = np.random.uniform(0, 1, size=100)
        h, e = np.histogramdd(array([a,b,c]).T, bins=(10, 10, 10))

        x,y,z = np.mgrid[0:100:10j, 0:100:10j, 0:1:10j]
        xx = x.ravel()
        yy = y.ravel()
        zz = z.ravel()
        hh = h.ravel()
        mlab.points3d(xx,yy,zz,hh)

        pl.show()

class BiasModeling:

    def __init__(self, g1t, g2t, g1p, g2p, bias_model='linear', bin_no=30, do_plot=False, sigma=1e10, boot_real=100, boot_sample=None):
        self.initialize(g1t, g2t, g1p, g2p, bin_no, sigma, do_plot)

        self.binning(valid=None)

        if bias_model == 'linear':
            self.linear_bias()
            self.linear_bias_cov_boot(boot_real=boot_real, 
                                  boot_sample=boot_sample)
            self.linear_bias_error()
        elif bias_model == 'linear_evolve':
            self.linear_evolve_bias()
        elif bias_model == 'nonlinear':
            self.nonlinear_bias()
        else:
            print 'Unknown bias model. Stopping'
        

    def initialize(self, g1t, g2t, g1p, g2p, bin_no, sigma, do_plot):
        self.g1t = g1t.ravel()
        self.g2t = g2t.ravel()
        self.g1p = g1p.ravel()
        self.g2p = g2p.ravel()
  
        self.gt = abs(self.g1t + 1j*self.g2t)
        self.gp = abs(self.g1p + 1j*self.g2p)

        self.bin1 = np.linspace(self.g1p.min(), self.g1p.max(), bin_no)
        self.bin2 = np.linspace(self.g2p.min(), self.g2p.max(), bin_no)
        self.bin = np.linspace(self.gp.min(), self.gp.max(), bin_no)
        
        self.sigma = sigma
        self.do_plot = do_plot
        self.rN = 5

    def binning(self, valid=None, boot=False):
        #Based on Amara et al. Not that the true value is binned
        #for a fixed predicted value.
        if valid is None:
            valid = np.arange(self.g1t.shape[0])

        self.g1p_b, g1p_be, self.g1t_b, self.g1t_be, N1, B1 = \
        MyF.AvgQ(self.g1p[valid], self.g1t[valid], self.bin1, sigma=self.sigma) 
        self.g2p_b, g2p_be, self.g2t_b, self.g2t_be, N2, B2 = \
        MyF.AvgQ(self.g2p[valid], self.g2t[valid], self.bin2, sigma=self.sigma) 
        self.gp_b, gp_be, self.gt_b, self.gt_be, N, B = \
        MyF.AvgQ(self.gp[valid], self.gt[valid], self.bin, sigma=self.sigma) 


        self.g1t_be[self.g1t_be == 0] = 9999.
        self.g2t_be[self.g2t_be == 0] = 9999.
        self.gt_be[self.gt_be == 0] = 9999.
        N1[N1 == 0] = 1
        N2[N2 == 0] = 1

        #print self.g1p_b, g1p_be, self.g1t_b, self.g1t_be
     
        if self.do_plot and boot is False: 
            gN = self.g1t_b.shape[0]
            N, E = np.histogramdd(np.array([self.g1t, self.g1p]).T, bins=(self.bin1, self.bin1))
            pl.subplot(121) 
            pl.contourf(N, origin='lower', extent=[self.bin1[0], self.bin1[-1], self.bin1[0], self.bin1[-1]])
            #pl.colorbar()
            pl.scatter(self.g1p, self.g1t, s=0.01) 
            pl.scatter(self.g1p[B1==1], self.g1t[B1==1], c='r', s=5.01, edgecolor='') 
            pl.scatter(self.g1p[B1==3], self.g1t[B1==3], c='b', s=5.01, edgecolor='') 
            pl.scatter(self.g1p[B1==7], self.g1t[B1==7], c='g', s=5.01, edgecolor='') 
            #for xx, yy in zip(self.g1p_b, self.g1t_b):
            # print xx, yy

            pl.errorbar(self.g1p_b[self.rN:gN-self.rN], self.g1t_b[self.rN:gN-self.rN], self.g1t_be[self.rN:gN-self.rN]/np.sqrt(N1[self.rN:gN-self.rN]), c='r')
            pl.xlabel(r'$\gamma_1^p$')
            pl.ylabel(r'$\gamma_1^t$')
            pl.xticks([-0.02,-0.015,-0.01,-0.005,0.0,0.005, 0.01, 0.015, 0.02])
            pl.yticks([-0.02,-0.015,-0.01,-0.005,0.0,0.005, 0.01, 0.015, 0.02])
            pl.axis([-0.015, 0.015, -0.015, 0.015])
            #pl.axis([-0.05, 0.05, -0.06, 0.06])
            pl.subplot(122)
            pl.scatter(self.g2p, self.g2t, s=1.01)
            pl.errorbar(self.g2p_b, self.g2t_b, self.g2t_be/np.sqrt(N2), c='r')
            pl.xlabel(r'$\gamma_2^p$')
            pl.ylabel(r'$\gamma_2^t$')
            pl.xticks([-0.02,-0.015,-0.01,-0.005,0.0,0.005, 0.01, 0.015, 0.02])
            pl.yticks([-0.02,-0.015,-0.01,-0.005,0.0,0.005, 0.01, 0.015, 0.02])
            #pl.axis([-0.05, 0.05, -0.06, 0.06])
            pl.axis([-0.015, 0.015, -0.015, 0.015])
            pl.show()
            #self.g1t_be /= np.sqrt(N1)
            #self.g2t_be /= np.sqrt(N2)


    def linear_bias(self, boot=False):
        """Predicted gamma = gamma_g * 1/b. The parameter b[0] used here
           is 1/b not b"""
        gN = self.g1t_b.shape[0] 
        b_init = [1]
        chi2_1 = lambda b: np.sum(((self.g1t_b[self.rN:gN-self.rN] - 
                                    self.g1p_b[self.rN:gN-self.rN] * b[0]) 
                                    / self.g1t_be[self.rN:gN-self.rN])**2)
        bias1 = optimize.fmin(chi2_1, b_init)
        chi2_2 = lambda b: np.sum(((self.g2t_b[self.rN:gN-self.rN] - 
                                    self.g2p_b[self.rN:gN-self.rN] * b[0]) 
                                   / self.g2t_be[self.rN:gN-self.rN])**2)
        bias2 = optimize.fmin(chi2_2, b_init)
        chi2_3 = lambda b: np.sum(((self.gt_b[self.rN:gN-self.rN] - 
                                    self.gp_b[self.rN:gN-self.rN] * b[0]) 
                                    / self.gt_be[self.rN:gN-self.rN])**2)
        bias = optimize.fmin(chi2_3, b_init)

        if boot is False:
            print 'Bias1 %2.2f'%(1 / bias1)
            print 'Bias2 %2.2f'%(1 / bias2)
            print 'Bias %2.2f'%(1 / bias)

        return 1/bias1, 1/bias2, 1/bias

    def linear_bias_cov_boot(self, boot_real=20, boot_sample=None):
        """Bias error using bootstrap"""
        trN = self.rN
        self.rN = 3
        if boot_sample is None:
            boot_sample = self.g1t.shape[0]

        b1_arr, b2_arr, b_arr = [], [], []
        for i in range(boot_real):
            #print 'Boot sample > %d'%i
            valid = np.random.randint(0, self.g1t.shape[0], boot_sample)
            s1 = np.std(self.g1p[valid])
            s2 = np.std(self.g2p[valid])
            s = np.std(self.gp[valid])
         
            b1 = np.cov(self.g1t[valid] / s1, self.g1p[valid] / s1)
            b2 = np.cov(self.g2t[valid] / s2, self.g2p[valid] / s2)
            b = np.cov(self.gt[valid] / s, self.gp[valid] / s)
           
            b1_arr.append(1/b1[0][1])
            b2_arr.append(1/b2[0][1])
            b_arr.append(1/b[0][1])
        b1_arr.sort()
        b2_arr.sort()
        b_arr.sort()

        #pl.hist(b1_arr, histtype='step', color='r', label='g1')
        #pl.hist(b2_arr, histtype='step', color='k', label='g2')
        #pl.legend() 
        #pl.show()

        larg = np.floor(0.16*boot_real - 1).astype(int)
        harg = np.floor(0.84*boot_real - 1).astype(int)

        self.b1_med = np.median(b1_arr)
        self.b1_l, self.b1_h =  self.b1_med - b1_arr[larg], b1_arr[harg] - self.b1_med
                         
        self.b2_med = np.median(b2_arr)
        self.b2_l, self.b2_h =  self.b2_med - b2_arr[larg], b2_arr[harg] - self.b2_med

        self.b_med = np.median(b_arr)
        self.b_l, self.b_h =  self.b_med - b_arr[larg], b_arr[harg] - self.b_med
  
        print 'b1 (boot) = %2.2f - %2.2f + %2.2f'%(self.b1_med, self.b1_l, self.b1_h) 
        print 'b2 (boot) = %2.2f - %2.2f + %2.2f'%(self.b2_med, self.b2_l, self.b2_h) 
        print 'b (boot) = %2.2f - %2.2f + %2.2f'%(self.b_med, self.b_l, self.b_h) 
        self.rN = trN


    def linear_bias_boot(self, boot_real=20, boot_sample=None):
        """Bias error using bootstrap"""
        trN = self.rN
        self.rN = 3
        if boot_sample is None:
            boot_sample = self.g1t.shape[0]

        b1_arr, b2_arr, b_arr = [], [], []
        for i in range(boot_real):
            print 'Boot sample > %d'%i
            boot_samples = np.random.randint(0, self.g1t.shape[0],
                                                      boot_sample)
            self.binning(valid=boot_samples)

            b1, b2, b = self.linear_bias(boot=True)    
            b1_arr.append(b1)
            b2_arr.append(b2)
            b_arr.append(b)
        b1_arr.sort()
        b2_arr.sort()
        b_arr.sort()

        larg = np.floor(0.16*boot_real - 1).astype(int)
        harg = np.floor(0.84*boot_real - 1).astype(int)

        self.b1_med = np.median(b1_arr)
        self.b1_l, self.b1_h =  self.b1_med - b1_arr[larg], b1_arr[harg] - self.b1_med
                         
        self.b2_med = np.median(b2_arr)
        self.b2_l, self.b2_h =  self.b2_med - b2_arr[larg], b2_arr[harg] - self.b2_med

        self.b_med = np.median(b_arr)
        self.b_l, self.b_h =  self.b_med - b_arr[larg], b_arr[harg] - self.b_med
  
        print 'b1 (boot) = %2.2f - %2.2f + %2.2f'%(self.b1_med, self.b1_l, self.b1_h) 
        print 'b2 (boot) = %2.2f - %2.2f + %2.2f'%(self.b2_med, self.b2_l, self.b2_h) 
        print 'b (boot) = %2.2f - %2.2f + %2.2f'%(self.b_med, self.b_l, self.b_h) 
        self.rN = trN


    def return_bias(self, chi2):
        m = minuit.Minuit(chi2)
        m.migrad()
        m.hesse()
        b = m.values['b']
        be = m.errors['b']
        bias = 1/b
        bias_e = be/b**2.
        print 'Bias = %2.2f \pm %2.2f'%(bias, bias_e)
        return bias, bias_e 

    def linear_bias_error(self):
        """Predicted gamma = gamma_g * 1/b. The parameter b[0] used here
           is 1/b not b"""
        gN = self.g1t_b.shape[0] 
        b_init = [1]
        chi2_1 = lambda b: np.sum(((self.g1t_b[self.rN:gN-self.rN] - 
                                    self.g1p_b[self.rN:gN-self.rN] * b) 
                                  / self.g1t_be[self.rN:gN-self.rN])**2)
        self.bias1, self.bias1_e = self.return_bias(chi2_1)
        chi2_2 = lambda b: np.sum(((self.g2t_b[self.rN:gN-self.rN] - 
                                    self.g2p_b[self.rN:gN-self.rN] * b) 
                               / self.g2t_be[self.rN:gN-self.rN])**2)
        self.bias2, self.bias2_e = self.return_bias(chi2_2)
        chi2_3 = lambda b: np.sum(((self.gt_b[self.rN:gN-self.rN] - 
                                    self.gp_b[self.rN:gN-self.rN] * b) 
                               / self.gt_be[self.rN:gN-self.rN])**2)
        self.bias3, self.bias3_e = self.return_bias(chi2_3)


    def linear_evolve_bias(self):
        print 'Not yet'
    def nonlinear_bias(self):
        print 'Not yet'


def linear_bias_kappa(kt, kp):
    kt = kt.ravel()
    kp = kp.ravel()
    kt = sigma_clip(kt, 8, 10)
    kp = kp[~kt.mask]
    kt = kt.compressed()
    print kp.shape, kt.shape
    bin = np.linspace(kp.min(), kp.max(), 30)
    kp_b, kp_be, kt_b, kt_be, N, B =  MyF.AvgQ(kp, kt, bin)
    N[N == 0] = 1
    kt_be /= np.sqrt(N)
    b_init = [1]
    kt_be[kt_be == 0] = 99.0
    chi2 = lambda b: np.sum(((kt_b - kp_b * b[0]) / kt_be)**2)
    bias = optimize.fmin(chi2, b_init)
    print 'Bias ', 1/bias

    
    chi2 = lambda b: np.sum(((kt_b - kp_b * b) / kt_be)**2)
    m = minuit.Minuit(chi2)
    m.migrad()
    m.hesse()
    b = m.values['b']
    be = m.errors['b']
    bias = 1/b
    bias_e = be/b**2.
    print 'b = %2.2f \pm %2.2f'%(bias, bias_e)
    return bias, bias_e

if __name__=='__main__':
   
    ipath = '.'
    ifile = 'kappamap_im3shape_r_1.0_0.5_g1.fits'
    opath = '.'
    k = KappaAmara(ipath, ifile, opath)    
    k.delta_rho_3d(50, 50, 10)
    k.kappa_predicted()
