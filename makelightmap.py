import pyfits as pf


class lightmap:
    def __init__(self, ra, dec, z, appmag, opath, pixel_scale,
                 bin_ra, bin_dec, bin_z, mask,ipath='catalogs/',
                 sourcefile='background.fits',
                 lensfile='foreground.fits',
                 zs=0.8,pdf_zs=None, zmin_s=0.6, zmax_s=1.5,
                 zmin_l=0.1, zmax_l=0.5):
        
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

        self.omega_M_0 = 0.3
        self.omega_lambda_0 = 0.7
        self.h = 0.72

        self.z = z
        self.ra = ra
        self.dec = dec

    
    def integrate(self,x,y):
        delta = (x[1:] - x[:-1])
        yavgs = (y[:-1] + y[1:]) / 2.
        return np.sum(delta*yavgs)

    def E(self,z):
        self.omega_k = 1. - self.omega_M_0 - self.omega_lambda_
        return np.asarray([1./np.sqrt(self.omega_M_0*(1.+zi)**3+self.omega_k*(1.+zi)**2+self.omega_lambda_0) for zi in z])

    def lum_from_appmag(self,appmag,z):
        c = 3*10**5 #parsec/sec                                                                                                                                                       
        zi = np.arange(0.01,z,.01)
        Dc = c/H_0 * integrate(zi,E(zi))
        parsec = Dc*1000000
        parsec = 780000
        print 'parsec '+str(parsec)
        absmag = -5.0 - appmag + 5.0*math.log(parsec,10)
        print 'absmag '+str(absmag)
        lum = 2.5118864315098**(4.83-absmag)
        print 'lum '+str(lum)



if __name__ == '__main__':
    #lum_from_appmag(3.44,1.2)
