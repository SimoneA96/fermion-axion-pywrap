import numpy as np
import os, sys, glob 
import time 
import pylab as pl
import subprocess

path_mass = './'+sys.argv[1]

final = []

for file in glob.glob(path_mass+"/blackline_*.txt"):
    filename = os.path.basename(file) 
    mass = float(filename[10:-4]) 
    
    data=np.loadtxt(file)
    
    find_max=data[:,2]
    for i in range(1,len(find_max)) :
    	if (find_max[i]-find_max[i-1]) < 0. :
    		final.append([data[i,0],data[i,1]])
    		break

#new_final=sorted(final,key=lambda l:l[0])

f = open('blackline_'+sys.argv[1]+'.txt', "a")

for xx in final :
    f.write(str(xx[0])+'\t'+str(xx[1])+'\n')
    
f.close()
