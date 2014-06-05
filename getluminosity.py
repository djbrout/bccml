import numpy as np
import math

def integrate(x,y):
    delta = (x[1:] - x[:-1])
    yavgs = (y[:-1] + y[1:]) / 2.
    return np.sum(delta*yavgs)

def E(z,omega_matter=.3,omega_lambda=.7):
    omega_k = 1. - omega_matter - omega_lambda
    return np.asarray([1./np.sqrt(omega_matter*(1.+zi)**3+omega_k*(1.+zi)**2+omega_lambda) for zi in z])

def lum_from_appmag(appmag,z):
    c = 3*10**5 #parsec/sec
    H_0 = 70.
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
    lum_from_appmag(3.44,1.2)
