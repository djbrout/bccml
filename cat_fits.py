import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pylab as plt
import os
import pylab as p
import rdfits as r
import mytools
import sys
import pyfits as pf

mag_cut = 22.5

file_dir = "/data3/data2/home/clampitt/bcc_v1.0/bcc_v1.0_truth_orig/"
OUTDIR = "./catalogs"

if not os.path.exists('figures'):
	        os.makedirs('figures')

title_in = ""

if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
	title = str(title_in)
newOUTDIR = OUTDIR+"/"

fits = ['Aardvark_v1.0_truth.180.fit','Aardvark_v1.0_truth.307.fit'
       ]
tables = []

for fit in fits:
	tables.append(pf.open(file_dir+fit))

cols = []
for table in tables:
	cols.append(table[1].data)

z = []
photoz = []
TMAGr = []
AMAGr = []
OMAGr = []
RA = []
DEC = []
GAMMA1 = []
GAMMA2 = []
K = []
e1 = []
e2 = []
for col in cols:
	z.extend(col["Z"])
	photoz.extend(col["PHOTOZ_GAUSSIAN"])
	TMAGr.extend(col["TMAG"][:,2])
        AMAGr.extend(col['AMAG'][:,2])
        OMAGr.extend(col['OMAG'][:,2])
	RA.extend(col["RA"])
	DEC.extend(col["DEC"])
	GAMMA1.extend(col["GAMMA1"])
	GAMMA2.extend(col["GAMMA2"])
	K.extend(col["KAPPA"])
	e1.extend(col["EPSILON"][:,0])
	e2.extend(col["EPSILON"][:,1])



z = np.asarray(z)
photoz = np.asarray(photoz)
TMAGr = np.asarray(TMAGr)
AMAGr = np.asarray(AMAGr)
OMAGr = np.asarray(OMAGr)
RA = np.asarray(RA)
DEC = np.asarray(DEC)
GAMMA1 = np.asarray(GAMMA1)
GAMMA2 = np.asarray(GAMMA2)
K = np.asarray(K) 
e1 = np.asarray(e1)
e2 = np.asarray(e2)

bg = [(z >.5) & (z < .6) & (TMAGr < mag_cut)] # background galaxies are where z > zcut = .5
#print TMAGr[0:100]
print bg
RAbg = RA[bg]
DECbg = DEC[bg]
GAMMA1bg = GAMMA1[bg]
GAMMA2bg = GAMMA2[bg]
weightsbg = np.ones(GAMMA1[bg].shape)
zbg = z[bg]

fg = [(z < .5) & (z > .1) & (TMAGr < mag_cut)] # foreground galaxies are where z < zcut = .5
RAfg = RA[fg]
DECfg = DEC[fg]
GAMMA1fg = GAMMA1[fg]
GAMMA2fg = GAMMA2[fg]
weightsfg = np.ones(GAMMA1[fg].shape)
zfg = z[fg]
#print fg

fig = plt.figure()
plt.hist(zbg,30, normed=0)
plt.xlabel("redshift")
plt.ylabel("Source Distribution Counts")
plt.title("z (0.5 - 0.6)")
fig.savefig("./figures/source_distribution.png")
fig = plt.figure()
plt.hist(zfg,30, normed=0)
plt.xlabel("redshift")
plt.ylabel("Fg Distribution Counts")
plt.title("z (0.1 - 0.5)")
fig.savefig("./figures/fg_distribution.png")
        
if os.path.exists(newOUTDIR+'foreground.fits'):
    os.remove(newOUTDIR+'foreground.fits')
mytools.write_fits_table(newOUTDIR+'foreground.fits', ['z','RA','DEC'], [zfg,RAfg,DECfg])
if os.path.exists(newOUTDIR+'background.fits'):
    os.remove(newOUTDIR+'background.fits')
mytools.write_fits_table(newOUTDIR+'background.fits', ['RA','DEC','S1','S2','W','z'],
                         [RAbg,DECbg,GAMMA1bg,GAMMA2bg,weightsbg,zbg])
if os.path.exists(newOUTDIR+'catalogue.fits'):
    os.remove(newOUTDIR+'catalogue.fits')
mytools.write_fits_table(newOUTDIR+'catalogue.fits',
			 ['RA','DEC','z','S1','S2','TMAGr','OMAGr','AMAGr','KAPPA','PHOTOZ','E1','E2'],
			 [RA, DEC, z, GAMMA1, GAMMA2, TMAGr, OMAGr, AMAGr, K, photoz, e1, e2])

sys.exit()
