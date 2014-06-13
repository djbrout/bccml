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

pid = os.getpid()

def get_mem():
        lines = open('/proc/%d/status' % pid).readlines()
        print '\n'.join([L for L in lines if L.find('VmSize') != -1])

def getsquareregion(file_dir,topleft_secondrow,topright_firstrow,bottomleft_thirdrow):
	tl_sr = pf.open(file_dir+topleft_secondrow)
	tr_fr = pf.open(file_dir+topright_firstrow)
	bl_tr = pf.open(file_dir+bottomleft_thirdrow)
	tl_sr_cols = tl_sr[1].data
	tr_fr_cols = tr_fr[1].data
	bl_tr_cols = bl_tr[1].data
	minra = min(np.asarray(tl_sr_cols['RA']))
	maxra = max(np.asarray(tr_fr_cols['RA']))
	mindec = min(np.asarray(bl_tr_cols['DEC']))
	maxdec = max(np.asarray(tl_sr_cols['DEC']))

	return [minra,maxra,mindec,maxdec]
	
print 'MEM'
print get_mem()

#mag_cut = 22.5
mag_cut = 25.0

file_dir = "/home/dbrout/bccml/corrected_healpix/"
OUTDIR = "./catalogs"

if not os.path.exists('figures'):
	        os.makedirs('figures')

title_in = ""

if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
	title = str(title_in)
newOUTDIR = OUTDIR+"/"

#500 Contiguous Sq Deg Area #1
fits = ['aardvark_v1.0_hpix_truth.428.fit',
        'aardvark_v1.0_hpix_truth.429.fit',
        'aardvark_v1.0_hpix_truth.430.fit',
        'aardvark_v1.0_hpix_truth.431.fit',
        'aardvark_v1.0_hpix_truth.460.fit',
        'aardvark_v1.0_hpix_truth.461.fit',
        'aardvark_v1.0_hpix_truth.462.fit',
        'aardvark_v1.0_hpix_truth.463.fit',
        'aardvark_v1.0_hpix_truth.492.fit',
	'aardvark_v1.0_hpix_truth.493.fit',
	'aardvark_v1.0_hpix_truth.494.fit',
	'aardvark_v1.0_hpix_truth.495.fit',
	'aardvark_v1.0_hpix_truth.524.fit',
	'aardvark_v1.0_hpix_truth.525.fit',
	'aardvark_v1.0_hpix_truth.526.fit',
	'aardvark_v1.0_hpix_truth.527.fit'
       ]

esign = [1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 -1,
	 1,
	 1,
	 1,
	 -1,
	 1,
	 1,
	 -1,
	 -1
	 ]


topleft_secondrow = 'aardvark_v1.0_hpix_truth.460.fit'
topright_firstrow = 'aardvark_v1.0_hpix_truth.431.fit'
bottomleft_thirdrow = 'aardvark_v1.0_hpix_truth.493.fit'

#500 Contiguous Sq Deg Area #2
#fits = ['aardvark_v1.0_hpix_truth.400.fit',
#	'aardvark_v1.0_hpix_truth.401.fit',
#	'aardvark_v1.0_hpix_truth.402.fit',
#	'aardvark_v1.0_hpix_truth.403.fit',
#	'aardvark_v1.0_hpix_truth.432.fit',
#	'aardvark_v1.0_hpix_truth.433.fit',
#	'aardvark_v1.0_hpix_truth.434.fit',
#	'aardvark_v1.0_hpix_truth.435.fit',
#	'aardvark_v1.0_hpix_truth.464.fit',
#	'aardvark_v1.0_hpix_truth.465.fit',
#	'aardvark_v1.0_hpix_truth.466.fit',
#	'aardvark_v1.0_hpix_truth.467.fit',
#	'aardvark_v1.0_hpix_truth.496.fit',
#	'aardvark_v1.0_hpix_truth.497.fit',
#	'aardvark_v1.0_hpix_truth.498.fit',
#	'aardvark_v1.0_hpix_truth.499.fit'
#	]

#topleft_secondrow = 'aardvark_v1.0_hpix_truth.432.fit'
#topright_firstrow = 'aardvark_v1.0_hpix_truth.403.fit'
#bottomleft_thirdrow = 'aardvark_v1.0_hpix_truth.464.fit'

where = getsquareregion(file_dir,topleft_secondrow,topright_firstrow,bottomleft_thirdrow)
minra = where[0]
maxra = where[1]
mindec = where[2]
maxdec = where[3]


tables = []
indx=0
for fit in fits:
    indx+=1
    print indx
    tables.append(pf.open(file_dir+fit))
    print file_dir+fit



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
#e1 = []
#e2 = []

cols = []
index = -1
for table in tables:
    cc =table[1].columns
    index += 1
    #table[1].data['GAMMA1'] = table[1].data['GAMMA1']*esign[index]
    #table[1].data['GAMMA2'] = table[1].data['GAMMA2']*esign[index]
    print cc.names
    col = table[1].data
    #cols.append(table[1].data)
    #z.extend(col["Z"])
    z.extend(np.asarray(col.field('z')))
    #photoz.extend(col["PHOTOZ"])
    photoz.extend(np.asarray(col.field("PHOTOZ")))
    #TMAGr.extend(col["TMAG"])
    TMAGr.extend(np.asarray(col.field("TMAG")))
    #AMAGr.extend(col['AMAG'])
    AMAGr.extend(np.asarray(col.field('AMAG')))
    #OMAGr.extend(col['OMAG'])
    OMAGr.extend(np.asarray(col.field('OMAG')))
    #RA.extend(col["RA"])
    RA.extend(np.asarray(col.field("RA")))
    #DEC.extend(col["DEC"])
    DEC.extend(np.asarray(col.field("DEC")))
    #GAMMA1.extend(col["GAMMA1"])
    GAMMA1.extend(np.asarray(col.field("GAMMA1")).astype(float)*esign[index])
    #GAMMA2.extend(col["GAMMA2"])
    GAMMA2.extend(np.asarray(col["GAMMA2"]).astype(float)*esign[index])
    K.extend(np.asarray(col.field("KAPPA")))
    #K.extend(col["KAPPA"])

#z = []
#photoz = []
#TMAGr = []
#AMAGr = []
#OMAGr = []
#RA = []
#DEC = []
#GAMMA1 = []
#GAMMA2 = []
#K = []
##e1 = []
##e2 = []
#index = -1
#for col in cols:
#    index += 1
#    z.extend(col["Z"])
#    photoz.extend(col["PHOTOZ"])
#    TMAGr.extend(col["TMAG"])
#    AMAGr.extend(col['AMAG'])
#    OMAGr.extend(col['OMAG'])
#    RA.extend(col["RA"])
#    DEC.extend(col["DEC"])
#    GAMMA1.extend(col["GAMMA1"])
#    GAMMA2.extend(col["GAMMA2"])
#    K.extend(col["KAPPA"])
#    #e1.extend(col["EPSILON"][:,0])
#    #e2.extend(col["EPSILON"][:,1])
#    cols[index] = 0


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
#e1 = np.asarray(e1)
#e2 = np.asarray(e2)

ww = [(RA > minra) & (RA < maxra) & (DEC > mindec) & (DEC < maxdec)]
z = z[ww]
photoz = photoz[ww]
TMAGr = TMAGr[ww]
AMAGr = AMAGr[ww]
OMAGr = OMAGr[ww]
RA = RA[ww]
DEC = DEC[ww]
GAMMA1 = GAMMA1[ww]
GAMMA2 = GAMMA2[ww]
K = K[ww]

print 'step 1'
bg = [(z >.5) & (z < .6) & (TMAGr < mag_cut)] # background galaxies are where z > zcut = .5
#print TMAGr[0:100]
#print bg
RAbg = RA[bg]
DECbg = DEC[bg]
GAMMA1bg = GAMMA1[bg]
GAMMA2bg = GAMMA2[bg]
weightsbg = np.ones(GAMMA1[bg].shape)
zbg = z[bg]
TMAGrbg = TMAGr[bg]
AMAGrbg = AMAGr[bg]
OMAGrbg = OMAGr[bg]


print 'step 2'
fg = [(z < .5) & (z > .1) & (TMAGr < mag_cut)] # foreground galaxies are where z < zcut = .5
RAfg = RA[fg]
DECfg = DEC[fg]
GAMMA1fg = GAMMA1[fg]
GAMMA2fg = GAMMA2[fg]
weightsfg = np.ones(GAMMA1[fg].shape)
zfg = z[fg]
Kfg = K[fg]
TMAGrfg = TMAGr[fg]
AMAGrfg = AMAGr[fg]
OMAGrfg = OMAGr[fg]

#print fg

print 'plotting'
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

print 'saving'
if os.path.exists(newOUTDIR+'foreground.fits'):
    os.remove(newOUTDIR+'foreground.fits')
mytools.write_fits_table(newOUTDIR+'foreground.fits', ['z','RA','DEC','tmagr','omagr','amagr','kappa'], [zfg,RAfg,DECfg,TMAGrfg, OMAGrfg, AMAGrfg,Kfg])
if os.path.exists(newOUTDIR+'background.fits'):
    os.remove(newOUTDIR+'background.fits')
mytools.write_fits_table(newOUTDIR+'background.fits', ['RA','DEC','S1','S2','W','z'],
                         [RAbg,DECbg,GAMMA1bg,GAMMA2bg,weightsbg,zbg])
if os.path.exists(newOUTDIR+'catalogue.fits'):
    os.remove(newOUTDIR+'catalogue.fits')
mytools.write_fits_table(newOUTDIR+'catalogue.fits',
			 ['RA','DEC','z','S1','S2','TMAGr','OMAGr','AMAGr','KAPPA','PHOTOZ'],
			 [RA, DEC, z, GAMMA1, GAMMA2, TMAGr, OMAGr, AMAGr, K, photoz])

sys.exit()
