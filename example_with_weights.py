import ks_mapping as ks
import os
import scipy.ndimage as nd
import numpy.ma as ma
import mytools as my
import config as c
import redmapper as rm
from weighted_galkappa_2d import WeightedGalKappa

def GenKappa(bg_f, fg_f):
    pixel_scale = c.pixel_scale #0.9375 # size of the pixel in arc min
    ipath = c.ipath
    opath = c.opath
    pipe = c.pipe
    filter = c.filter

    smooth_size = c.smooth_size #arcmin
    smooth_size_n = c.smooth_size_n #arcmin
    rotate = c.rotate #'_b' 
    esign = c.esign

    if esign[0] == 1 and  esign[1] == 1:
        sign = 'nsign'
    elif esign[0] == -1 and  esign[1] == -1:
        sign = 'sign'
    elif esign[0] == -1:
        sign = 'g1'
    elif esign[1] == -1:
        sign = 'g2'
      

    #root file name 
    file_root='%s_%s_%.1f_%.1f_%s'%(pipe, filter, pixel_scale, 
                                       smooth_size, sign)

    #generating kappamap
    k = ks.KappaMap(ipath, bg_f, opath, pixel_scale, skip=0, 
                 lens_quantity='shear', rotate=rotate, 
                 randomize=c.randomize,
                 constrain=c.constrain, coord=c.coord, 
                 bin_ra=None, bin_dec=None, 
                 project=c.project, reference_ra=c.reference_ra)

    xmin, xmax, ymin, ymax = k.ra_min, k.ra_max, k.dec_min, k.dec_max
    #saving kappamap
    k.savekappa_fits(smooth_size=smooth_size, file_root=file_root)

    #Uncomment following in the future to get bootstrap error
    ''' 
    ofile, mask = k.pixelize_shear(k.ra, k.dec, k.g1, k.g2, k.w, 
                                   savethis=False)

    k.gamma_to_kappa(k.gamma, pixel_scale)

    boot_f = 'bootstrap_kappa_%s.npz'%file_root
    if os.path.exists(boot_f):
        pass
    else:
        k.bootstrap(boot_realiz=c.boot_realiz, boot_sample=None, 
                    smooth_size=smooth_size, file_root=file_root)
    f = np.load(boot_f)
    '''

    #generating foreground galaxy number count

    n = ks.KappaMap(ipath, fg_f, opath, pixel_scale, skip=0, 
                 lens_quantity='count', rotate=rotate,
                 randomize=c.randomize, constrain=c.constrain, 
                 coord=c.coord, bin_ra=k.ra_b, bin_dec=k.dec_b,
                 project=c.project, reference_ra=c.reference_ra)

    ofile, mask = n.pixelize_galcount(n.ra, n.dec, 
                  smooth_size=smooth_size,
                  savethis=True, file_root=file_root, fg='fg')


    wn = WeightedGalKappa(n.ra, n.dec, n.z, '.', c.smooth_size, pixel_scale, n.ra_b, n.dec_b, 5, mask.T, zs=[.5,.6,.7,.8,.9,1.,1.1,1.2])

    wn.delta_rho_3d()
    wn.kappa_predicted()
    os.system('mv '+c.opath+'kappa_predicted.npz '+c.opath+'kappa_predicted_%s.npz'%file_root)

    return k.gamma, k.ra_b, k.dec_b, bg_f, fg_f

if __name__=='__main__':
    bg_f = 'background.fits' #background file. columsn should be RA, DEC, G1, G2, W (weight, just one for equal weight)

    fg_f = 'foreground.fits' #foreground file. columns should be RA DEC

    gamma, ra_b, dec_b, bg_f, fg_f = GenKappa(bg_f, fg_f)
