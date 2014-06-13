import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import sys


def integrate(x,y):
    delta = (x[1:] - x[:-1])
    yavgs = (y[:-1] + y[1:]) / 2.
    return np.sum(delta*yavgs)

def E(z,omega_matter=.3,omega_lambda=.7):
    omega_k = 1. - omega_matter - omega_lambda
    return np.asarray([1./np.sqrt(omega_matter*(1.+zi)**3+omega_k*(1.+zi)**2+omega_lambda) for zi in z])

c = 3*10**5
H_0 = 70.
pi = 3.14159
G = 6.67*10**(-11)*10**(-9) #m^3kg^-1s^-2 convert to km *10^-9
print G
z_l = np.arange(0.7,1.5,.01)
z_s = z_l*0.0 + 0.5

z = np.asarray([np.arange(0.01,zl,.01) for zl in z_l])
#print z[1]
#print [zi[-1] for zi in z]
#print [integrate(zi,E(zi)) for zi in z]

Dc_l = np.asarray([c/H_0 * integrate(zi,E(zi)) for zi in z])
#print Dc_l
Dm_l = Dc_l #for zero curvature
#print Dm_l.shape
#print z_l.shape
D_l = Dm_l#/ (1.0+z_l)



z = np.asarray([np.arange(0.01,zs,.01) for zs in z_s])
Dc_s = np.asarray([c/H_0 * integrate(zi,E(zi)) for zi in z]) 
zt = np.asarray([.1,.2,.3,.4,.5])
#print 'EEEEE'
#print E(zt)
#print 'INTEGRATE'
#print integrate(zt,E(zt))
#print c/H_0
Dm_s = Dc_s #for zero curvature
D_s =  Dm_s#/ (1.0+z_s)
#print Dc_s
#print Dc_l
D_ls = D_l - D_s
print D_ls

shape_noise = z_l*0.0 + .01
e = shape_noise

sigma_e = z_l*0.0 + .01

#Mandelbaum et al. eq 2
critical_mass_surface_density = -c**2/(4*pi*G)*D_s/((1.+z_l)**2*D_l*D_ls)

#Mandelbaum et al. eq 4
weights = 1./(critical_mass_surface_density)
#weights = (1.+z_l)*D_l*D_s/D_s
#where e is the rms ellipticity and sigma_e is the ellipticity measurement error per component

#print D_s
#print D_l
#print weights

plt.figure(1)
plt.plot(z_l,weights)
plt.xlabel('z')
plt.ylabel('weight')
plt.savefig("weights_test_zs=0.5.pdf")
plt.figure(2)
plt.plot(z_l,critical_mass_surface_density)
plt.xlabel('z')
plt.ylabel('critical_mass_surface_density')
plt.savefig("critical_mass_surface_density_test_zs=0.5.pdf")
