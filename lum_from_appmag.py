import numpy as np
import math

def integrate(x,y):
    delta = (x[1:] - x[:-1])
    yavgs = (y[:-1] + y[1:]) / 2.
    return np.sum(delta*yavgs)

def E(z):
    omega_lambda_0 = .7
    omega_M_0 = .3
    omega_k = 1. - omega_M_0 - omega_lambda_0
    return np.asarray([1./np.sqrt(omega_M_0*(1.+zi)**3+
                                  omega_k*(1.+zi)**2+
                                  omega_lambda_0) for zi in z])
                        


c = 3*10**5 #parsec/sec
H_0 = 72.0
zi = np.arange(0.01,.4,.01)
Dc = c/H_0 * integrate(zi,E(zi))
parsec = Dc*1000000
parsec = 778000
print 'parsec '+str(parsec)
absmag = 5.0 +3.44  - 5.0*math.log(parsec,10)
print 'absmag '
print absmag
#lum2d = 2.5118864315098**(4.83-absmag)
lum2d = 10**((4.83-absmag)/2.5)
print 'lum '
print lum2d
#self.save_fits_image(self.lum2d,'./maps/luminosity/lum_density.fits')
