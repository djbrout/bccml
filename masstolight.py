import numpy as np
import rdfits as r
import mytools
import sys
import pyfits as pf
import matplotlib.pyplot as plt
from pylab import *

class mlmap:
    def __init__(self,ks_massmap,true_massmap,lightmap,sigma_map,file_root=''):
        self.ks_massmap = pf.open(ks_massmap)[0].data
        self.massmap = pf.open(sigma_map)[0].data
        self.true_massmap = pf.open(true_massmap)[0].data
        self.lightmap = pf.open(lightmap)[0].data
        self.contour_plot(self.lightmap.ravel(),self.true_massmap.ravel(),bins=1000,
                          xlabel='Luminosity',ylabel='True Mass',title='True Mass to Light Contour for each pixel in 500sq deg BCC survey',
                          range=[[-1*10**11,1*10**11],[-1,1]],
                          file='./figures/truemass_to_light_2dhist.png')
        self.contour_plot(self.true_massmap.ravel(),self.ks_massmap.ravel(),bins=1000,
                          xlabel='True Mass',ylabel='KS Mass',title='Mass check for each pixel in 500sq deg BCC survey',
                          file='./figures/mass_check_2dhist.png')
        self.contour_plot(self.lightmap.ravel(),self.ks_massmap.ravel(),bins=1000,
                          xlabel='Luminosity',ylabel='KS Mass',title='KS Mass to Light Contour for each pixel in 500sq deg BCC survey',
                          range=[[-1*10**11,1*10**11],[-.02,.1]],
                          file='./figures/ksmass_to_light_2dhist.png')
        self.contour_plot(self.lightmap.ravel(),self.massmap.ravel(),bins=1000,
                          xlabel='Luminosity',ylabel='Sigma Mass',title='Mass to Light for each pixel in 500sq deg BCC survey',
                          file='./figures/sigmamass_to_light_2dhist.png',ylim=[-10**14,10**14])
        
        self.mass_hist(self.ks_massmap.ravel(),'./figures/ksmass_hist.png')
        self.mass_hist(self.true_massmap.ravel(),'./figures/truemass_hist.png')
        self.mass_hist(self.massmap.ravel(),'./figures/sigmamass_hist.png')
                
        self.crosscor()

    def mass_hist(self,vec,file):
        fig = plt.figure()
        plt.hist(vec,100, normed=0)
        plt.xlabel('mass')
        plt.ylabel('# of pixels')
        #plt.xlim(-1,1)
        savefig(file)
    def save_fits_image(self,image,filename):
        hdu = pf.PrimaryHDU(image)
        if os.path.exists(filename):
            os.remove(filename)
            hdu.writeto(filename)
        return

    def contour_plot(self,x,y,bins=100,range=None,xlabel='',ylabel='',file='test.png',title='',xlim=None,ylim=None):
        #y = self.true_massmap.ravel()
        #x = self.lightmap.ravel()
        fig = plt.figure()
        self.grid,xedges,yedges = np.histogram2d(x,y,bins=bins,range=range)
        print len(xedges)
        print len(yedges)
        plt.scatter(x,y,s=.05)
        plt.contour(xedges[0:-1],yedges[0:-1],self.grid.T)
        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.title(title)
        savefig(file)
                
    def crosscor(self):

        ksmass_light_corr = np.cov(self.ks_massmap.ravel(),self.lightmap.ravel())
        truemass_light_corr = np.cov(self.true_massmap.ravel(),self.lightmap.ravel())
        mass_light_corr = np.cov(self.massmap.ravel(),self.lightmap.ravel())
        
        print 'ksmass to light: ',(ksmass_light_corr)
        print 'truemass to light: ',(truemass_light_corr)
        print 'sigmamass to light: ',(mass_light_corr)
        

        

        fig = plt.figure()
        plt.scatter(self.lightmap.ravel(),self.true_massmap.ravel(),s=.1)
        plt.xlabel("luminosity (solar lum)")
        plt.ylabel("True Mass")
        fig.savefig("./figures/truemass_to_light.png")
        fig = plt.figure()
        plt.scatter(self.lightmap.ravel(),self.ks_massmap.ravel(),s=.1)
        plt.xlabel("luminosity (solar lum)")
        plt.ylabel("KS Mass")
        fig.savefig("./figures/ksmass_to_light.png")
        fig = plt.figure()
        plt.scatter(self.true_massmap.ravel(),self.ks_massmap.ravel(),s=.1)
        plt.plot([0,1],[0,1])
        plt.xlabel("True mass")
        plt.ylabel("KS Mass")
        fig.savefig("./figures/ksmass_check.png")        
        fig = plt.figure()
        plt.scatter(self.lightmap.ravel(),self.massmap.ravel(),s=.1)
        
        #plt.contour(self.lightmap.ravel(),self.massmap.ravel())
        #plt.plot([0,1],[0,1])
        plt.ylabel("Mass (Solar Mass)")
        plt.xlabel("Luminosity (Solar Lum)")
        plt.ylim(-1*10**14,1*10**14)
        fig.savefig("./figures/sigma_mass_to_light.png")
        
        
if __name__=='__main__':
    ks_massmap = '/home/dbrout/bccml/ks_mapping/kappamap__full500_unprojected.fits'
    true_massmap = '/home/dbrout/bccml/maps/luminosity/kappa_density_full500_unprojected_simpleavg.fits'
    lightmap = '/home/dbrout/bccml/maps/luminosity/lum_density_full500_unprojected_simpleavg.fits'
    mass_map = '/home/dbrout/bccml/maps/luminosity/mass_full500_unprojected_simpleavg.fits'
    print 'a'
    a = mlmap(ks_massmap,true_massmap,lightmap,mass_map)
    print 'b'
