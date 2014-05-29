pixel_scale = 1.0 #size of the pixel in arc min 
ipath = '/data3/data2/home/dbrout' #input directory
opath = '.' #output directory

pipe = 'im3shape' # this is just for naming of output file
filter = 'r' #this is just for naming of output file

smooth_size = 0.5 #smoothing scale for background arcmin 
smooth_size_n = 0.5 #smoothing scale for background arcmin 

rotate = False #rotate=True will rotate galaxies 45 deg. good for bmode
randomize = False #randomize shear catalog
esign = [-1, 1] #sign of the g1 and g2

boot_realiz = 2 #bootstrap realizations

constrain = False
coord = [70., 12, -53., 12] #[ra, ra_width/2., dec, dec_width/2.]
#If contrain=True, then only galaxies within the coord will be used
project = False #project=False does not work with astropy.wcs
reference_ra = 70. #If project=True, then this will be the reference     
