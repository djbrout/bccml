import pyfits as pf
import rdcol
import matplotlib.pyplot as plt
import os
import pyfits as pf
import numpy as np
import math

class lightmap:
    def __init__(self, ipath='',
                 lensfile='./catalogs/foreground.fits',z_lens=.4):

        self.bin_ra = [line.strip() for line in open('./maps/kappa/kappa_ra_bins.txt')]
        self.bin_dec = [line.strip() for line in open('./maps/kappa/kappa_dec_bins.txt')]

        self.lensfilepath = os.path.join(ipath,lensfile)
        self.lensfile = pf.open(self.lensfilepath)
        #self.kmap = pf.open(kappafit)
        
        self.omega_M_0 = 0.3
        self.omega_lambda_0 = 0.7
        self.h = 0.72

        self.z_lens = z_lens

        
        print self.lensfile[1].columns
        self.cols = self.lensfile[1].data
        
        self.tmag_l = self.cols['tmagr']
        self.ra_l = self.cols['ra']
        self.dec_l = self.cols['dec']
        self.z_l = self.cols['z']
        self.mag_map()
        self.lum_from_appmag()

    def mag_map(self):
        print self.dec_l
        print self.ra_l
        self.tmag2d, edges = np.histogramdd(np.array([self.dec_l,self.ra_l]).T,
                                            bins=(self.bin_dec, self.bin_ra),
                                            weights=self.tmag_l)
        self.raedges = edges[1]
        self.decedges = edges[0]
        self.raavg = (self.raedges[:-1] + self.raedges[1:]) / 2.
        self.decavg = (self.decedges[:-1] + self.decedges[1:]) / 2.

        self.save_fits_image(self.tmag2d,'./maps/luminosity/tmag_density.fits')
        
        #np.savez('./maps/luminosity/tmag_density.npz',
        #         raedge=self.raavg, decedge=self.decavg,
        #         tmag2d = self.tmag2d)
        #self.plot_2dmap(self.tmag2d,'True Apparent Mag Map','./maps/luminosity/tmag_2d_map.pdf')
        #NEED to test histogram and compare to kappamap then check numerical values
        #and then convert each bin(pixel) into luminosity
        return

    def integrate(self,x,y):
        delta = (x[1:] - x[:-1])
        yavgs = (y[:-1] + y[1:]) / 2.
        return np.sum(delta*yavgs)

    def E(self,z):
        self.omega_k = 1. - self.omega_M_0 - self.omega_lambda_0
        return np.asarray([1./np.sqrt(self.omega_M_0*(1.+zi)**3+
                                      self.omega_k*(1.+zi)**2+
                                      self.omega_lambda_0) for zi in z])

    def lum_from_appmag(self):
        c = 3*10**5 #parsec/sec
        H_0 = 72.0
        zi = np.arange(0.01,self.z_lens,.01)
        Dc = c/H_0 * self.integrate(zi,self.E(zi))
        parsec = Dc*1000000
        #parsec = 780000
        print 'parsec '+str(parsec)
        print 'tmag'
        print self.tmag2d
        absmag = 5.0 + self.tmag2d - 5.0*math.log(parsec,10)
        print 'absmag '
        print absmag
        self.lum2d = 10**((4.83-absmag)/2.5)
        print 'lum '
        print self.lum2d
        self.save_fits_image(self.lum2d,'./maps/luminosity/lum_density.fits')
               
        #self.plot_2dmap(self.lum2d,'Luminosity Map','./maps/luminosity/lum_2d.pdf')
        #self.save_fits_image(self.lum2d,'./maps/luminosity/lum_density.fits')
        return

    def plot_2dmap(self,image,title,filename):
        fig = plt.figure(figsize=(6, 3.2))

        ax = fig.add_subplot(111)
        ax.set_title(title)
        plt.imshow(self.tmag2d)
        ax.set_aspect('equal')

        #cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
        #cax.get_xaxis().set_visible(False)
        #cax.get_yaxis().set_visible(False)
        #cax.patch.set_alpha(0)
        #cax.set_frame_on(False)
        plt.colorbar(orientation='vertical')
        #plt.show()
        plt.savefig(filename)

    def save_fits_image(self,image,filename):
        hdu = pf.PrimaryHDU(image)
        if os.path.exists(filename):
            os.remove(filename)
        hdu.writeto(filename)
        return


if __name__ == '__main__':
    a = lightmap()
