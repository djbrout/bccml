#! /usr/global/paper/bin/python

import pyfits as pf
#import matplotlib.pyplot as plt
import os
import pyfits as pf
import numpy as np
import math
#import pylab as p
import sys

class lightmap:
    def __init__(self, ipath='',
                 lensfile='/home/dbrout/bccml/catalogs/foreground.fits',z_lens=.4):

        self.bin_ra = [line.strip() for line in open('/home/dbrout/bccml/maps/kappa/kappa_ra_bins.txt')]
        self.bin_dec = [line.strip() for line in open('/home/dbrout/bccml/maps/kappa/kappa_dec_bins.txt')]

        self.lensfilepath = os.path.join(ipath,lensfile)
        self.lensfile = pf.open(self.lensfilepath)
        #self.kmap = pf.open(kappafit)
        
        self.omega_M_0 = 0.3
        self.omega_lambda_0 = 0.7

        self.c = 3*10**5 #parsec/sec
        self.H_0 = 72.0
        self.omega_k = 1. - self.omega_M_0 - self.omega_lambda_0
        
        #self.z_lens = z_lens
        
        #print self.lensfile[1].columns
        self.cols = self.lensfile[1].data
        
        self.tmag_l = self.cols.field('tmagr')
        self.ra_l = self.cols.field('ra')
        self.dec_l = self.cols.field('dec')
        self.z_l = self.cols.field('z')
        self.kappa = self.cols.field('kappa')
        print 'step 1'
        self.dist_lookup()
        self.plot_lum()

    def get_mem(self):
            lines = open('/proc/%d/status' % pid).readlines()
            print '\n'.join([L for L in lines if L.find('VmSize') != -1])

    def dist_lookup(self):
        self.distances = []
        zrange = np.arange(0.1,0.9,.001)
        for zz in zrange:
            zi = np.arange(0.01,zz,.01)
            [self.distances.append(self.c/self.H_0 *
                                   self.integrate(zi,self.E(zi))
                                   * 1000000)]
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
        self.lum_l = 10**((4.83-self.absmag)/2.5)
        #self.lum_l = np.asarray(col['luminosity'])
        #self.dec_l = np.asarray(col['dec'])
        #self.ra_l = np.asarray(col['ra'])
        #self.lum_l = self.lum_l - np.mean(self.lum_l)
        print 'pixelizing'
        self.lum2d, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                           bins=(self.bin_dec, self.bin_ra),
                                           weights=self.lum_l)
        self.tmag2d, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                           bins=(self.bin_dec, self.bin_ra),
                                           weights=self.tmag_l)
        self.kappa2d, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                            bins=(self.bin_dec, self.bin_ra),
                                            weights=self.kappa)
        print 'saving'
        self.save_fits_image(self.lum2d,'/home/dbrout/bccml/maps/luminosity/lum_density.fits')
        self.save_fits_image(self.tmag2d,'/home/dbrout/bccml/maps/luminosity/mag_density.fits')
        self.save_fits_image(self.kappa2d,'/home/dbrout/bccml/maps/luminosity/kappa_density.fits')

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
