'''
 *  The aim of this code is to compute the equilibrium configurations
 *  for a Fermion-Axion star system in a 2 parameters space of solutions. 
 *  
 *
 *  LUNCH:
 *     python3 wrapper.py path_directory_with_equal_mass_files fa-number &
 *
 *
 *  OUTPUT: various txt file devided in different directories 
 *      + the existence plot of the pure axion star to be verified with
 *      the solutions in arXiv:1909.05515 
 *		 
 *
 
 
 *
 *   Copyright (c) 2021-2022 Saeed Fakhry (the main author of this code),
 *   Simone Albanesi, Fabrizio di Giovanni, Davide Guerra, Miquel Miravet
 *
 *   This file is part of FAS (Fermion-Axion Star program).
 *
 *   FAS is free software (we accept a coffee or a beer if you want); 
 *   you can redistribute it and/or modify it under the terms of the 
 *   people that develope it.
 *
 *   This script will use the codes Shoot.cpp and Solve.cpp in order
 *   to build the 2-D existence plot of a fermion-axion star. 
 *
 *


 *
 *
 *   Revision 1.1  2021/12/02  21:00:00  Davide Guerra
 *       *** first patch ***
 *
 *
 *
'''

import numpy as np
import os, sys, glob 
import time 
import pylab as pl
import subprocess

#--------------------------------------------------------------
# define some functions
#--------------------------------------------------------------
def append2file(fname, string):
    fp = open(fname, 'a')
    fp.write(string+'\n')
    fp.close()
    return


#--------------------------------------------------------------
# check input, check/create directories
#--------------------------------------------------------------
'''
if len(sys.argv)<3:
    print("Error usage of the wrapper!\nPass to the script the parfile's name",  \
          "and the directory where to save the output \n(the directory", \
          " will be created in ", out_dir, ")\nExample:\n   >> python3 ", \
          "wrapper.py path_directory_with_equal_mass_files axionstar &", sep='')
    sys.exit()
else:
    


if not os.path.exists(out_dir):
    os.system('mkdir '+out_dir)
if not os.path.exists(out_final_dir):
    os.system('mkdir '+out_final_dir)
if not os.path.exists(bin_dir):
    os.system('mkdir '+bin_dir)
'''
#--------------------------------------------------------------
# read the file Input.txt (read by Shoot) and run
#--------------------------------------------------------------

path_mass       = sys.argv[1]    
    
os.system('g++ -o Shoot Shoot.cpp root_shoot.cpp -w')
os.system('g++ -o Solve Solve.cpp root_shoot.cpp -w')
print('\n Compiling Shoot.cpp and Solve.cpp:\n')

for file in glob.glob(path_mass+"/out_equalmass_*.txt"):
    filename = os.path.basename(file) 
    mass = float(filename[14:-4]) 
    print('Fixed mass value : ', mass)

    os.system('cp '+file+' Input.txt ')
    print("file copied\n")
    
    time.sleep(1.0)
    
    os.system('./Shoot dummy >/dev/null')
    print('shoot done\n')
    os.system('./Solve dummy >/dev/null')
    print('solve done\n')
    time.sleep(1.0)
    os.system('mv dummy/blackline* /home/davide/fermion-axion-pywrap/blackline/out/fa-0.02/blackline_'+str(mass)+'.txt')
    print('move blacline.txt\n')
    time.sleep(1.0)
    os.system('rm -r dummy')
    print('remove dummy directory\n\n')


#--------------------------------------------------------------

