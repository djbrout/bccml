#! /usr/global/paper/bin/python
from math import pi
import numpy as np

from matplotlib.pyplot import *

import healpy as hp
import pyfits as pf
import os

pid = os.getpid()

def get_mem():
    lines = open('/proc/%d/status' % pid).readlines()
    print '\n'.join([L for L in lines if L.find('VmSize') != -1])

indir = 'bcc_v1.0_truth_orig/'
outdir = 'bcc_v1.0_hpix_photoz/'
#outdir = 'temp_output/'
form = '%.8e'

joe = '/data2/home/clampitt/bcc_v1.0/bcc_v1.0_truth_orig/'
local '/home/dbrout/bccml/corrected_healpix/'
nside = 8

# Read command line argument
#if len(sys.argv) != 2:
#    sys.exit('Must provide one value.')
#tpix = int(sys.argv[1]) - 1

ind = int(os.environ['SGE_TASK_ID']) - 1
#ind = 0
#tpix = np.loadtxt('potential_pix.txt')[ind]
tpix = np.loadtxt('des_pix.txt')[ind]
#tpix = np.loadtxt('potential_pix.txt')[0]

spixfile = joe+'source_files_into_pix%d.txt' % (tpix)
spix = np.loadtxt('source_pix_lists/' +spixfile)
print 'source pixels = ', spix

outfile = local+'aardvark_v1.0_hpix_truth.%d.fit' % (tpix)

#key = ['ra', 'dec', 'm200', 'central', 'amag', 'tmag', 'z', 'photoz_gaussian']
#form = ['E', 'E', 'E', 'I', '5E', '5E', 'E', 'E']
key = ['ra', 'dec','z', 'photoz',
       'tmag','omag','amag','gamma1','gamma2','kappa']
form = ['E', 'E', 'E', 'E',
        'E', 'E', 'E', 'E', 'E', 'E']

ct = 0
new_data = {}
for j in range(len(spix)):

    if (spix[j] == 5000): continue

    infile = joe+'Aardvark_v1.0_truth.%d.fit' % (spix[j])
    
    hdulist = pf.open(indir +infile)
    #print hdulist[1].columns, '\n\n'
    ra = hdulist[1].data.field('ra')
    dec = hdulist[1].data.field('dec')
    z = hdulist[1].data.field('z')

    abs_r = hdulist[1].data.field('AMAG')[:,1]
    m_r = hdulist[1].data.field('TMAG')[:,1]
    m_g = hdulist[1].data.field('TMAG')[:,0]

    # Obtain gaussian photoz with appropriate scatter for LRG
    zgauss = np.random.normal(z, 0.03*(1.+z), len(z))

    # Shift RA to continuous interval
    con = (ra > 200.)
    ra[con] = ra[con] - 360.

    theta = 90. - dec
    loc = hp.ang2pix(nside, theta * pi/180., ra * pi/180.)
    pcut = (loc == tpix)

    ct = ct + len(ra[pcut])
    print 'j, num = ', j, '\t', ct
    print get_mem()

    for k in range(len(key)):
        if (key[k] == 'ra'):
            if (j == 0): new_data[key[k]] = ra[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], ra[pcut]))
        elif (key[k] == 'abs_r'):
            if (j == 0): new_data[key[k]] = abs_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], abs_r[pcut]))
        elif (key[k] == 'm_r'):
            if (j == 0): new_data[key[k]] = m_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], m_r[pcut]))
        elif (key[k] == 'm_g'):
            if (j == 0): new_data[key[k]] = m_g[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], m_g[pcut]))
        elif (key[k] == 'photoz'):
            if (j == 0): new_data[key[k]] = zgauss[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], zgauss[pcut]))
        else:
            #print key[k]
            if (j == 0): new_data[key[k]] = hdulist[1].data.field(key[k])[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], hdulist[1].data.field(key[k])[pcut]))
                                                
    hdulist.close()

tmpcols = []
for i in range(len(key)):
    tmpcols.append(pf.Column(name=key[i], format=form[i],
                             array=new_data[key[i]]))

hdu = pf.PrimaryHDU()
tbhdu = pf.new_table(tmpcols)
thdulist = pf.HDUList([hdu, tbhdu])

thdulist.writeto(outdir +outfile)
thdulist.close()
