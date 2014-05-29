import os
import pyfits as pf
import fnmatch
import numpy as np
import matplotlib.pylab as plt

filedir = '/data3/data2/home/clampitt/bcc_v1.0/bcc_v1.0_truth_orig/'
color_index = 0
colors = plt.cm.rainbow(np.linspace(0, 1, 15))
for file in os.listdir(filedir):
    if fnmatch.fnmatch(file,'*.*.fit'):
        if color_index < 10000:
            num = file.split('.')[2]
            color_index += 1
            print color_index
            table = pf.open(filedir+file)
            cols = table[1].data
            rai = np.interp(range(0,len(cols["RA"]),1000),range(len(cols["RA"])),cols["RA"])
            deci = np.interp(range(0,len(cols["DEC"]),1000),range(len(cols["DEC"])),cols["DEC"])
            ra, dec = np.meshgrid(rai,deci)
            
            clr_indx = ra*dec*0.0+color_index
            print 'plotting'
            plt.scatter(rai,deci,color=colors[color_index%9])
            a = np.mean(rai)
            b = np.mean(deci)
            print num
            plt.annotate(str(num),[a,b])

plt.xlabel("RA")
plt.ylabel("DEC")
plt.title(filedir.split('/')[-1]+' Map')
plt.savefig('./figures/bcc_map.pdf')
