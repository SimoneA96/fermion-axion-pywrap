'''
 *    python scritp reproducing the 2-D existence plots in arXiv:2006.08583
 *    for fermion-axion star. 
 *
 *    RUN: python3 PlotsDensPhi.py value_for_fa
 *
'''

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
rc('font',family='serif')
import glob, os, sys, math

#=======================================================

def ReadFile(path):
    os.chdir(path)
    for file in glob.glob("existence*.txt"):
        existence_file = file
    data   = np.loadtxt(path+existence_file)
    lendata = len(data[:,0])
    
    j = 1
    for i in range(1,lendata):
        if data[i,0] < data[i-1,0]:
            break
        j = j+1
    nphi = j
    
    k = 1
    for i in range(1,lendata):
        if data[i,3] != data[i-1,3]:
            k = k+1
    nrho = k
    
    print('Number of points for central potential   : ', nphi)
    print('Number of points for central rest mass density   : ', nrho)
    
    phictmp   = np.array(np.array_split(data[:,0],nrho))
    phic = phictmp[0,:]
    rhotmp = data[:,3]
    rhoc   = np.zeros(nrho)
    
    for i in range(0,nrho):
        rhoc[i] = rhotmp[nphi*i]
    
    mass = np.array(np.array_split(data[:,2],nrho))
    
    return phic, rhoc, mass

#===========================================================

home_path = os.getcwd()
path_contour = home_path+'/out/fa-'+sys.argv[1]+'/'
path_plots = path_contour+'plots/'
    
try:
    os.makedirs(path_plots)
except OSError as error:
    print(error)

fa = float(sys.argv[1])
print('fa   = ',fa)
print('\n log(fa) = ', math.log(fa,10))

phic, rhoc, mass = ReadFile(path_contour)
mass = mass[:,1:]
mass = np.swapaxes(mass,0,1)
phic = phic[1:]

fig, ax = plt.subplots(figsize = (12,8)) 
c       = ax.pcolor(rhoc,phic, mass,cmap = 'coolwarm')
cbar    = fig.colorbar(c, ax=ax)
cbar.set_label(label = r'$M_{\rm T}\mu$', size = 16, weight = 'bold')

ax.set_xlabel(r'$\rho_c/\mu^2$', fontsize = 16)
ax.set_ylabel(r'$\phi_c$', fontsize = 16)
ax.set_title(r'$f_a = '+sys.argv[1]+'$', fontsize = 20, fontweight= 'bold')
#ax.set_ylim(top = 0.14)
ax.tick_params(axis = 'both', labelsize = 15)
cbar.ax.tick_params(labelsize=15)
    
fig.savefig(path_plots+'Count_fa_'+sys.argv[1]+'_rho_vs_phi.png', dpi=400)
