import os
import pyfits as pf
import fnmatch
import numpy as np
import matplotlib.pylab as plt
import sys
#filedir = '/data3/data2/home/clampitt/bcc_v1.0/bcc_v1.0_hpix_truth/'
#f = open('/data3/data2/home/clampitt/bcc_v1.0/des_pix.txt','r')

filedir = '/data2/home/clampitt/bcc_v1.0/bcc_v1.0_truth_orig/'
#os.system('ls /data2/home/clampitt/bcc_v1.0/bcc_v1.0_truth_orig/*.fit > /home/dbrout/bccml/bcc_v1.0_truth_orig.txt')
#f = open('/home/dbrout/bccml/bcc_v1.0_truth_orig.txt','r')

#ff= set()
#for filenum in f:
#  #ff.add(filenum.replace("\n",""))
#  ff.add(filenum.replace("\n","").split('.')[3])

#print ff
#sys.exit()
color_index = 0
colors = plt.cm.rainbow(np.linspace(0, 1, 15))
nums = []
ays = []
bs = []
for file in os.listdir(filedir):
    if fnmatch.fnmatch(file,'*.*.fit'):
        if color_index < 10000:
            num = file.split('.')[2]
            #if num in ff:
            color_index += 1
            print color_index
            print 'reading'
            table = pf.open(filedir+file)
            #print 'donereading'
            cols = table[1].data
            #print 'donetabledata'
            #colra = table[1].data["RA"]
            #coldec = table[1].data["DEC"]
            #print 'testradone'
            #rai = np.interp(range(0,len(cols["RA"]),10),range(len(cols["RA"])),cols["RA"])
            #deci = np.interp(range(0,len(cols["DEC"]),10),range(len(cols["DEC"])),cols["DEC"])
            rand = [np.random.randint(0,len(cols['RA'])) for i in range(100)]
            print rand
            rai = [cols['RA'][i] for i in rand]
            deci = [cols['DEC'][i] for i in rand]
            #ra, dec = np.meshgrid(rai,deci)
            
            #clr_indx = ra*dec*0.0+color_index
            print 'plotting'
            plt.scatter(rai,deci,color=colors[color_index%9])
            a = np.mean(rai)
            b = np.mean(deci)
            print num
            nums.append(num)
            ays.append(a)
            bs.append(b)

#f.close()
for j in range(len(nums)):
    plt.annotate(str(nums[j]),[ays[j],bs[j]],size='xx-small')

plt.xlabel("RA")
plt.ylabel("DEC")
plt.title(filedir.split('/')[-1]+' Map')
plt.savefig('./figures/bcc_orig_v1.0_map.pdf')
