file_dir = "/data3/scratch/bcc_v1"
OUTDIR = "./out"
title_in = ""

import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pylab as plt
import os
import pylab as p
import rdfits as r
import mytools
import sys

if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
	title = str(title_in)
newOUTDIR = OUTDIR+"/"

import pyfits as pf
table1 = pf.open(file_dir+"/Aardvark_v1.0c_truth_des_rotated.114.fit")
table2 = pf.open(file_dir+"/Aardvark_v1.0c_truth_des_rotated.115.fit")
table3 = pf.open(file_dir+"/Aardvark_v1.0c_truth_des_rotated.86.fit")
table4 = pf.open(file_dir+"/Aardvark_v1.0c_truth_des_rotated.147.fit")

cols1 = table1[1].data	    
cols2 = table2[1].data
cols3 = table3[1].data
cols4 = table4[1].data

z = np.concatenate((cols1["Z"],cols2["Z"],cols3["Z"],cols4["Z"]))
photoz = np.concatenate((cols1["PHOTOZ_GAUSSIAN"],cols2["PHOTOZ_GAUSSIAN"],cols3["PHOTOZ_GAUSSIAN"],cols4["PHOTOZ_GAUSSIAN"]))
TMAGr = np.concatenate((cols1["TMAG"][:,2],cols2["TMAG"][:,2],cols3["TMAG"][:,2],cols4["TMAG"][:,2]))
RA = np.concatenate((cols1["RA"],cols2["RA"],cols3["RA"],cols4["RA"]))
DEC = np.concatenate((cols1["DEC"],cols2["DEC"],cols3["DEC"],cols4["DEC"]))
GAMMA1 = np.concatenate((cols1["GAMMA1"],cols2["GAMMA1"],cols3["GAMMA1"],cols4["GAMMA1"]))
GAMMA2 = np.concatenate((cols1["GAMMA2"],cols2["GAMMA2"],cols3["GAMMA2"],cols4["GAMMA2"]))
K = np.concatenate((cols1["KAPPA"],cols2["KAPPA"],cols3["KAPPA"],cols4["KAPPA"]))
e1 = np.concatenate((cols1["EPSILON"][:,0],cols2["EPSILON"][:,0],cols3["EPSILON"][:,0],cols4["EPSILON"][:,0]))
e2 = np.concatenate((cols1["EPSILON"][:,1],cols2["EPSILON"][:,1],cols3["EPSILON"][:,1],cols4["EPSILON"][:,1]))

mytools.write_fits_table('catalogue.fits', ['RA','DEC','Z','S1','S2','TMAGr','KAPPA','PHOTOZ','E1','E2'], [RA,DEC,z,GAMMA1,GAMMA2,TMAGr,K,photoz,e1,e2])

sys.exit()
