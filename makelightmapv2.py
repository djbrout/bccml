#! /usr/global/paper/bin/python
import sys
sys.path.append('/home/dbrout/bccml/ks_mapping')
import pyfits as pf
import matplotlib.pyplot as plt
import os
import pyfits as pf
import numpy as np
import math
#import pylab as p
import sys
import config as c

class lightmap:
    def __init__(self,bin_ra,bin_dec, ipath='/home/dbrout/bccml/catalogs/',
                 lensfile='foreground.fits',z_lens=.33,
                 srcfile='background.fits',file_root='',trim=False):

        #self.bin_ra = [line.strip() for line in open('/home/dbrout/bccml/maps/kappa/kappa_ra_bins.txt')]
        #self.bin_dec = [line.strip() for line in open('/home/dbrout/bccml/maps/kappa/kappa_dec_bins.txt')]

        self.bin_ra = bin_ra
        self.bin_dec = bin_dec

        self.file_root = file_root
        self.lensfilepath = os.path.join(ipath,lensfile)
        self.lensfile = pf.open(self.lensfilepath)
        self.srcfilepath = os.path.join(ipath,srcfile)
        self.srcfile = pf.open(self.srcfilepath)
        #self.kmap = pf.open(kappafit)
        
        self.omega_M_0 = 0.3
        self.omega_lambda_0 = 0.7

        self.c = 3*10**5 #km/sec
        self.H_0 = 72.0
        self.omega_k = 1. - self.omega_M_0 - self.omega_lambda_0
        
        #self.z_lens = z_lens

        self.src_cols = self.srcfile[1].data
        self.z_s = self.src_cols.field('z')
        #print self.lensfile[1].columns
        self.cols = self.lensfile[1].data
        
        self.tmag_l = self.cols.field('tmagr')
        self.ra_l = self.cols.field('ra')
        self.dec_l = self.cols.field('dec')
        self.z_l = self.cols.field('z')
        self.kappa = self.cols.field('kappa')
        print 'step 1'
        self.dist_lookup()
        self.com_dist_lookup()
        self.weighting_fcn()
        self.plot_lum()


    def get_mass_from_kappa(self):
        G = 4.302 * 10**(-3)#pc*M^-1*(km/s)^2 
        pi = 3.14159
        z_l = .33#Fixing z_l based on weighting function peak
        d_l = self.com_distances[(np.round(z_l,3)*1000-100).astype(int)]
        z_s = np.mean(self.z_s)
        print z_s
        d_s = self.com_distances[(np.round(z_s,3)*1000-100).astype(int)]
        print d_s
        d_ls = d_l - d_s
        critical_mass_surface_density = self.c**2*d_s/(4*pi*G*d_ls*d_l)#Using comoving distances, c in km/s
        sigma = self.kappa2d * critical_mass_surface_density#Surface Mass Density! Need to multiply by area of pixel to get mass
        #NEed to multiply by area of pixel at z_l
        theta = c.pixel_scale# in arcmin
        #3437.74677078 arcmin per radian
        theta = theta/3437.746#theta in radians 
        size_of_pixel_in_pc = theta*d_l
        area_of_pix_in_pc = size_of_pixel_in_pc**2 
        mass = sigma*area_of_pix_in_pc
        return mass
                
    def weighting_fcn(self):
        G = 6.67*10**(-11)*10**(-9) #m^3kg^-1s^-2 convert to km *10^-9
        pi = 3.14159
        z_l = np.arange(0.1,.5,.01)
        d_l = self.distances[(np.round(z_l,3)*1000-100).astype(int)]
        z_s = np.mean(self.z_s)
        d_s = self.distances[(np.round(z_s,3)*1000-100).astype(int)]
        d_ls = d_l - d_s
        print len(z_l)
        print len(d_ls)
        critical_mass_surface_density = -self.c**2/(4*pi*G)*d_s/((1.+z_l)**2*d_l*d_ls)
        weights = 1./(critical_mass_surface_density)

        plt.figure(1)
        plt.plot(z_l,weights)
        plt.xlabel('z lens')
        plt.ylabel('weight')
        plt.savefig("/home/dbrout/bccml/figures/weights_test_zs"+str(z_s)+".pdf")
        #plt.figure(2)
        #plt.plot(z_l,critical_mass_surface_density)
        #plt.xlabel('z lens')
        #plt.ylabel('critical_mass_surface_density')
        #plt.savefig("./figures/critical_mass_surface_density_test.pdf")
        

    def get_mem(self):
            lines = open('/proc/%d/status' % pid).readlines()
            print '\n'.join([L for L in lines if L.find('VmSize') != -1])

    def com_dist_lookup(self):
        self.com_distances = []
        zrange = np.arange(0.1,0.9,.001)
        for zz in zrange:
            zi = np.arange(0.01,zz,.01)
            [self.com_distances.append(self.c/self.H_0 * self.integrate(zi,self.E(zi)) * 1000000)]#Lum Dist = Com Dist * (1+z)
        self.com_distances = np.asarray(self.com_distances)
            
            
    def dist_lookup(self):
        self.distances = []
        zrange = np.arange(0.1,0.9,.001)
        for zz in zrange:
            zi = np.arange(0.01,zz,.01)
            [self.distances.append(self.c/self.H_0 *self.integrate(zi,self.E(zi))* 1000000 * (1+zz))]#Lum Dist = Com Dist * (1+z)
        self.distances = np.asarray(self.distances)
        
    def integrate(self,x,y):
        delta = (x[1:] - x[:-1])
        yavgs = (y[:-1] + y[1:]) / 2.
        return np.sum(delta*yavgs)

    def E(self,z):
        return np.asarray([1./np.sqrt(self.omega_M_0*(1.+zi)**3+
                                      self.omega_k*(1.+zi)**2+
                                      self.omega_lambda_0) for zi in z])

    def plot_lum(self):                 
        #print 'reading file'
        #col = rdcol.read('./catalogs/fg_lum_cat/lum_catalog.dat',1,2)
        print 'rounding'
        self.z_li = np.round(self.z_l,3)*1000-100
        self.z_li = self.z_li.astype(int)
        print 'size of self.tmag_l '+str(self.tmag_l.shape)
        print 'size of self.z_l '+str(self.z_li.shape)
        print 'self.z_l'
        print self.z_l[0:50]
        print 'self.z_li'
        print self.z_li[0:50]
        
        print 'self.distances[z_l]'
        print self.distances[self.z_li[0:50]]
        
        self.absmag = 5.0 + self.tmag_l - 5.0*np.log10(self.distances[self.z_li])
        self.lum_l = 10**((4.65-self.absmag)/2.5)#AB system SDSS r-band Solar AbsMagnitude
        #self.lum_l = np.asarray(col['luminosity'])
        #self.dec_l = np.asarray(col['dec'])
        #self.ra_l = np.asarray(col['ra'])
        self.lum_l = self.lum_l - np.mean(self.lum_l.ravel()) #Subtract out mean
        self.kappa = self.kappa - np.mean(self.kappa.ravel())
        print 'pixelizing'
        self.lum2dw, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                           bins=(self.bin_dec, self.bin_ra),
                                           weights=self.lum_l)
        self.grid, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                           bins=(self.bin_dec, self.bin_ra))

        self.lum2d = self.lum2dw/self.grid
        
        #self.tmag2d, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
        #                                   bins=(self.bin_dec, self.bin_ra),
        #                                   weights=self.tmag_l)
        self.kappa2dw, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                            bins=(self.bin_dec, self.bin_ra),
                                            weights=self.kappa)
        self.kappa2d = self.kappa2dw/self.grid
        self.mass2d = self.get_mass_from_kappa()
        self.mass2d = self.mass2d - np.mean(self.mass2d.ravel())
        print 'saving'
        self.save_fits_image(self.lum2d,'/home/dbrout/bccml/maps/luminosity/lum_density'+str(self.file_root)+'.fits')
        #self.save_fits_image(self.tmag2d,'/home/dbrout/bccml/maps/luminosity/mag_density.fits')
        self.save_fits_image(self.kappa2d,'/home/dbrout/bccml/maps/luminosity/kappa_density'+str(self.file_root)+'.fits')
        self.save_fits_image(self.mass2d,'/home/dbrout/bccml/maps/luminosity/mass'+str(self.file_root)+'.fits')
        print 'plotting 3'
        
        #plt.figure()
        #n, bins, patches = plt.hist(self.lum2d[(self.lum2d < .4*10**10) & (self.lum2d > 1001.0)],100,log=True, histtype='bar')
        #plt.xlabel('Solar Luminosity')
        #plt.ylabel('# of Pixels')
        #print 'saving'
        #plt.savefig('/home/dbrout/bccml/maps/luminosity/lum_hist.png')
        return

    def save_fits_image(self,image,filename):
        hdu = pf.PrimaryHDU(image)
        if os.path.exists(filename):
            os.remove(filename)
        hdu.writeto(filename)
        return


if __name__ == '__main__':
    a = lightmap()
