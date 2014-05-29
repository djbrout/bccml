
def run_catalogue(mag_cut,file_dir="",OUTDIR="./out"):
	#file_dir = "/data3/scratch/bcc_v1" 
	#file_dir = ""
	#OUTDIR = "./out"
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
	table1 = pf.open(file_dir+"catalogue.fits")
	cols = table1[1].data

	z=cols["Z"]
	RA=cols["RA"]
	DEC=cols["DEC"]
	GAMMA1=cols["S1"]
	GAMMA2=cols["S2"]
	TMAGr = cols["TMAGr"]
	
	bg = [(z >.5) & (z < 1.5) & (TMAGr < mag_cut)] # background galaxies are where z > zcut = .5
	RAbg = RA[bg]
	DECbg = DEC[bg]
	GAMMA1bg = GAMMA1[bg]
	GAMMA2bg = GAMMA2[bg]
	weightsbg = np.ones(GAMMA1[bg].shape)
	zbg = z[bg]

	fg = [(z < .5) & (TMAGr < mag_cut)] # foreground galaxies are where z < zcut = .5
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
	plt.title("Z_cut = 0.5")
	fig.savefig("source_distribution.png")
        fig = plt.figure()
	plt.hist(zfg,30, normed=0)
	plt.xlabel("redshift")
	plt.ylabel("Fg Distribution Counts")
	plt.title("Z_cut = 0.5")
	fig.savefig("fg_distribution.png")
	
	mytools.write_fits_table(OUTDIR+'foreground.fits', ['z','RA','DEC'], [zfg,RAfg,DECfg])
	mytools.write_fits_table(OUTDIR+'background.fits', ['RA','DEC','S1','S2','W','z'], [RAbg,DECbg,GAMMA1bg,GAMMA2bg,weightsbg,zbg])


run_catalogue(22.5,file_dir="/home/dbrout/",OUTDIR="./out")
