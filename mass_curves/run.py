'''
 *    RUN: python3 run.py value_for_fa
 *
'''

import os, sys
import matplotlib 
import matplotlib.pyplot as plt
import numpy as np

home_path = os.getcwd()
path_exist = home_path+'/results/fa-'+sys.argv[1]+'/'

#path_plots = path_exist+'plots/'
   
try:
    os.makedirs(home_path+'/results/')
except OSError as error:
    print(error)
        
try:
    os.makedirs(path_exist)
except OSError as error:
    print(error)
    
try:
    os.makedirs(path_exist+'forPLOT/')
except OSError as error:
    print(error)


os.system('g++ -o LevelCurves MainLevelCurves.cpp LevelCurves.cpp root_shoot.cpp')
os.system('./LevelCurves ')


# Taking the data for the plot
data=np.loadtxt('out_equalmass.txt', skiprows=4)

plt.clf()
plt.xlabel(r'$\rho_c/\mu^2$')
plt.ylabel(r'$\phi_c$')
plt.grid(True)
plt.xlim(0.0,0.008)
plt.ylim(0.0,0.30)
plt.plot(data[:,1],data[:,0])
plt.savefig('plot_rho_vs_phi.png', format='png', dpi=400)


with open('out_equalmass.txt', 'r') as fin:
    data_file = fin.read().splitlines(True)
    mass=data_file[0][:-2]

with open('out_equalmass_'+mass+'.txt', 'w') as fout:
    fout.writelines(data_file[1:])
    
os.system('mv out_equalmass_* '+path_exist+'/.')
os.system('mv plot_rho_vs_phi.png '+path_exist+'/.')
