import numpy as np
import os 
import sys
import time 
import pylab as pl
import subprocess

#--------------------------------------------------------------
# define some functions
#--------------------------------------------------------------
def input_string_to_var(var_name, input_string):
    input_vars_strings = ['KK', 'polind', 'fa_axion', 'm', 'Lambda']
    input_vars_bool    = ['give_axion_exp', 'only_input_txt']
    input_vars_float   = ['phi0c_start', 'phi0c_end', 'rho0c_start', 'rho0c_end', 'r_max',     \
                          'r_max_low_phi0c', 'fa_axion_exp', 'zoom_phi_width','zoom_rho_width', \
                          'zoom_phi_points', 'zoom_rho_points']
    input_vars_int     = ['Ng', 'tol', 'phi0c_N', 'rho0c_N', 'zoom_phi_value', 'zoom_rho_value']
    if var_name in input_vars_float:
        var = float(input_string)
    elif var_name in input_vars_int:
        var = int(input_string)
    elif var_name in input_vars_bool:
        var = (input_string=='True' or input_string=='true')
    elif var_name in input_vars_strings:
        var = input_string
    else:
        print(var_name, 'not found in the input-vars list!')
        sys.exit()
    return var 

def append2file(fname, string):
    fp = open(fname, 'a')
    fp.write(string+'\n')
    fp.close()
    return

def create_1d_grid(x1,x2,N,zoom_points=[],zoom_width=0.01,zoom_value=2):
    xgrid = np.linspace(x1, x2, N)
    if len(zoom_points)>0:
        for zp in zoom_points:
            point_idx = np.argmin(np.abs(np.array(xgrid)-zp))
            if zp+zoom_width<=x2 and zp-zoom_width>=x1:
                new_subpoints = np.linspace(zp-zoom_width, zp+zoom_width, zoom_value)
                j1 = np.where(np.array(xgrid)<(zp-zoom_width))[0][-1]+1
                j2 = np.where(np.array(xgrid)>(zp+zoom_width))[0][0]
                subvec1 = xgrid[0:j1]
                subvec2 = xgrid[j2:]
                xgrid   = np.concatenate((subvec1, new_subpoints, subvec2))
    return xgrid

#--------------------------------------------------------------
# check input, check/create directories
#--------------------------------------------------------------
out_dir  = './out/'
cpp_dir  = './cpp/'
bin_dir  = './bin/'
if len(sys.argv)<3:
    print("Error usage of the wrapper!\nPass to the script the parfile's name",  \
          "and the directory where to save the output \n(the directory", \
          " will be created in ", out_dir, ")\nExample:\n   >> python3 ", \
          "wrapper.py par/test.par axionstar &", sep='')
    sys.exit()
else:
    parfile         = sys.argv[1]
    out_subdir_name = sys.argv[2]

if not os.path.exists(out_dir):
    os.system('mkdir '+out_dir)
if not os.path.exists(bin_dir):
    os.system('mkdir '+bin_dir)

#--------------------------------------------------------------
# read parfile
#--------------------------------------------------------------
fp = open(parfile, 'r')
input_dict = {}
for line in fp:
    line_wo_comments = line.partition('#')[0]
    part      = line_wo_comments.partition('=')
    var_name  = part[0].replace(' ' , '')
    var_value = (part[2].replace(' ' , '')).replace('\n', '')
    if '[' in var_value and ']' in var_value:
        list_in_dict = []
        var_value = var_value.replace(' ' , '')
        var_value = var_value.replace('\n', '')
        var_value = var_value[1:-1]
        elements  = var_value.split(',');
        for elem in elements:
            if elem!='':
                var = input_string_to_var(var_name, elem)
                list_in_dict.append(var)
        input_dict[var_name] = list_in_dict
    else:
        input_dict[var_name] = input_string_to_var(var_name, var_value)
fp.close()

Ng              = input_dict['Ng']
tol             = input_dict['tol']
phi0c_N         = input_dict['phi0c_N'] 
phi0c_start     = input_dict['phi0c_start']
phi0c_end       = input_dict['phi0c_end']
rho0c_N         = input_dict['rho0c_N'] 
rho0c_start     = input_dict['rho0c_start']
rho0c_end       = input_dict['rho0c_end']
r_max           = input_dict['r_max']
r_max_low_phi0c = input_dict['r_max_low_phi0c']

zoom_phi_points = input_dict['zoom_phi_points'] 
zoom_phi_value  = input_dict['zoom_phi_value']
zoom_phi_width  = input_dict['zoom_phi_width']
zoom_rho_points = input_dict['zoom_rho_points'] 
zoom_rho_value  = input_dict['zoom_rho_value']
zoom_rho_width  = input_dict['zoom_rho_width']

if zoom_phi_value%2!=0 or zoom_rho_value%2!=0:
    append2file(wrapper_log,'zoom_values must be even! Exit...')
    sys.exit()

only_input_txt   = input_dict['only_input_txt']

cv               = {}
cv['KK']         = input_dict['KK']
cv['polind']     = input_dict['polind']
cv['fa_axion']   = input_dict['fa_axion']
cv['m']          = input_dict['m']
cv['lambda']     = input_dict['Lambda']
compilation_variables_keys = cv.keys()

if input_dict['give_axion_exp']:
    fa_axion_exp   = input_dict['fa_axion_exp']
    fa_axion       = '{:.15f}'.format(10**fa_axion_exp)
    cv['fa_axion'] = fa_axion

#--------------------------------------------------------------
# create log file and write input used
#--------------------------------------------------------------
outfiles_dir = out_dir+out_subdir_name
ltime=time.strftime("h%Hm%Ms%S", time.localtime())
wrapper_log  = 'log_wrappy_'+ltime+'.txt'
dashes       = '--------------------------------------------------'
octothorps   = '##################################################'

if os.path.exists(wrapper_log):
    wrapper_log = wrapper_log.replace('.txt', '_'+str(time.monotonic())+'.txt')

append2file(wrapper_log, octothorps+'\nStarting at '+time.strftime("%H:%M:%S", time.localtime()))
t0            = time.perf_counter()
dict_keys     = list(input_dict.keys())
input_message = '';
excluded_keys = list(compilation_variables_keys)
excluded_keys.append('Lambda')
for i in range(0, len(dict_keys)):
    aa = dict_keys[i]
    if not aa[0]=='_' and aa not in excluded_keys:
        input_message += '{:<16} : {:>6}'.format(aa,str(input_dict[aa]))+'\n'   
append2file(wrapper_log, octothorps+'\n'+input_message+dashes)

#--------------------------------------------------------------
# generate the file Input.txt (read by Shoot)
#--------------------------------------------------------------
fname = 'Input.txt'
fp    = open(fname, 'w')
if phi0c_N>1:
    phi0c_vec = create_1d_grid(phi0c_start, phi0c_end, phi0c_N, zoom_points=zoom_phi_points, \
                               zoom_width=zoom_phi_width, zoom_value=zoom_phi_value)
    phi0c_N   = len(phi0c_vec)
else: 
    phi0c_vec = [phi0c_start]
    phi0c_N   = 1
if rho0c_N>1:
    rho0c_vec = create_1d_grid(rho0c_start, rho0c_end, rho0c_N, zoom_points=zoom_rho_points, \
                               zoom_width=zoom_rho_width, zoom_value=zoom_rho_value)
    rho0c_N   = len(rho0c_vec)
else:
    rho0c_vec = [rho0c_start]
    rho0c_N   = 1
Nmodels = phi0c_N*rho0c_N
print(tol, Ng, '\n', Nmodels, '\n', file=fp) 
for rho0c in rho0c_vec:
    for phi0c in phi0c_vec:
        if phi0c<0.01 and r_max<r_max_low_phi0c:
            r_max2write = r_max_low_phi0c
        else:
            r_max2write = r_max
        fp.write('{:<12.6f} {:<12.6f} {:<12f}\n'.format(phi0c, rho0c, r_max2write))
fp.close()

if only_input_txt:
    sys.exit()

#--------------------------------------------------------------
# compile the code 
#--------------------------------------------------------------
os.system('cp '+cpp_dir+'seeds/Shoot.cpp '+cpp_dir+'Shoot.cpp')
os.system('cp '+cpp_dir+'seeds/Solve.cpp '+cpp_dir+'Solve.cpp')
compile_message = '>> Compiling Shoot.cpp and Solve.cpp:\n'
for cvk in compilation_variables_keys:
    var_def = cvk+'='+str(cv[cvk])
    sed_opt = '|'+cvk+'=VALUE|/'+var_def
    os.system('sed -i.bak "s/'+sed_opt+'/" '+cpp_dir+'Shoot.cpp')
    os.system('sed -i.bak "s/'+sed_opt+'/" '+cpp_dir+'Solve.cpp')
    os.system('rm -f '+cpp_dir+'*.bak')
    compile_message += var_def.replace('=', ' '*(16-len(cvk))+' : ')+'\n'
os.system('g++ -o '+bin_dir+'Shoot '+cpp_dir+'Shoot.cpp -w')
os.system('g++ -o '+bin_dir+'Solve '+cpp_dir+'Solve.cpp -w')
append2file(wrapper_log, compile_message+dashes)

#--------------------------------------------------------------
# run Shoot and estimate the shooting-time
#--------------------------------------------------------------
fout_shoot   = 'out_shoot.txt'
if os.path.isdir(outfiles_dir):
    append2file(wrapper_log, outfiles_dir+' already exists! Exit')
    sys.exit()
os.popen('mkdir '+outfiles_dir).read()
append2file(wrapper_log, '>> Shooting...')
t1 = time.perf_counter()
shoot_command = bin_dir+'Shoot ' + outfiles_dir + ' > ' +fout_shoot + '  2>&1'
process = subprocess.Popen(shoot_command, shell=True)
lines = 0
sleep_step = 0.1
while not os.path.exists(outfiles_dir+'/Output.txt'):
    time.sleep(sleep_step)
while lines<3:
    lines_str = os.popen('wc -l '+outfiles_dir+'/Output.txt').read()
    lines = int(lines_str.partition(' ')[0])
    time.sleep(sleep_step)
t2 = time.perf_counter()
append2file(wrapper_log, 'Estimated time needed: {:.2f} s'.format((t2-t1)/lines*Nmodels))
p = process.communicate()
t3 = time.perf_counter()
append2file(wrapper_log, 'Elapsed time: '+str(t3-t1)+' s\n'+dashes)

#--------------------------------------------------------------
# run Solve
#--------------------------------------------------------------
fout_solve   = 'out_solve.txt'
append2file(wrapper_log, '>> Solving...')
os.popen(bin_dir+'Solve ' + outfiles_dir + ' > ' +fout_solve + '  2>&1 ').read()
os.popen('mv '+fout_shoot+' '+fout_solve+' '+outfiles_dir).read()
os.system('mv Input.txt '+outfiles_dir)

#--------------------------------------------------------------
# do a simple plot in the case of axion stars
#--------------------------------------------------------------
if rho0c_N==1:
    X    = np.loadtxt(outfiles_dir+'/existence_plot{:.6f}'.format(rho0c)+'.txt')
    phic = X[:,0]
    Mbar = X[:,2]
    f, ax = pl.subplots()
    pl.plot(phic, Mbar)
    pl.xlabel(r'$\Phi_c$')
    pl.ylabel(r'$M$')
    fa_exp_str = '{:.1f}'.format(np.log10(float(fa_axion)))
    pl.title(r'$f_a = 10^{'+fa_exp_str+'}$')
    ax.yaxis.get_ticklocs(minor=True)
    ax.minorticks_on()
    pl.savefig(outfiles_dir+'/existence_plot.png')

tend = time.perf_counter()
append2file(wrapper_log, dashes+'\nTotal elapsed time: '+str(tend-t0)+' s')
