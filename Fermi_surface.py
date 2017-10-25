# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 00:26:05 2016

@author: Yongjin
How to use:
    $ python Fermi_surface.py [filepath, filename]
Input:
    File name of EIGENVAL files. If not stated, 'EIGENVAL' is default.
Output:
    Filename of output. Defalut is 'Xcrysden.bxsf'

    This script extract spatial eigen values from EIGENVAL files for drawing
fermi-surface. Mind that the calculation should be done with pecular k-point mesh,
which is without ISYM-tag. Output file is the body for bxcf file, which is used for Xcrysden file.
Spin down elements are listed after up elements.0
    
In case of less than three argument is given, priority coms to input file.
For example,
1. $ python Fermi_surface.py
    Read 'OUTCAR' 'EIGENVAL' file, and write in 'Xcrysden.bxsf' file
2. $ python Fermi_surface.py OUTCAR.fermi EIGENVAL_SLAO Fermi_SLAO
    Read from 'OUTCAR.fermi', 'EIGENVAL_SLAO', and write in 'Fermi_SLAO' file
3. $ python Fermi_surface.py MY_FILE
    Read OUTCAR from 'MY_FILE', 'EIGENVAL', and write in 'Xcrysden.bxsf' file.
"""

import re
import numpy as np
import sys
import os


print ("\n### Conversion script: from VASP EIGENVAL file to Xcrysden input format ")
print ("### Author: Yongjin Shin\n")

#Case of different arguments are given.
if len(sys.argv) == 1:
    print ("Defalut Setting:\n")
    print ("Input OUTCAR file: OUTCAR")
    print ("Input EIGENVAL file: EIGENVAL\n")
    print ("Output file: Xcrysden.bxsf")
    if not os.path.isfile('OUTCAR'):
        print ("\n** ERROR: Input file 'OUTCAR' was not found.")
        sys.exit(0)    
    if not os.path.isfile('EIGENVAL'):
        print ("\n** ERROR: Input file 'EIGENVAL' was not found.")
        sys.exit(0)    
    outcar_filename='OUTCAR'
    in_filename='EIGENVAL'    
    out_filename='Xcrysden.bxsf'
elif len(sys.argv) == 2:
    print ("One external argument is given: \n")
    print ("Input OUTCAR file: {0}".format(sys.argv[1]))
    print ("Input EIGENVAL file: EIGENVAL\n")
    print ("Output file: Xcrysden.bxsf")
    if not os.path.isfile(sys.argv[1]):
        print ("\n** ERROR: Input file {0} was not found.".format(sys.argv[1]))
        sys.exit(0)
    if not os.path.isfile('EIGENVAL'):
        print ("\n** ERROR: Input file 'EIGENVAL' was not found.")
        sys.exit(0)    
    outcar_filename=sys.argv[1]
    in_filename='EIGENVAL'
    out_filename='Xcrysden.bxsf'
elif len(sys.argv) == 3:
    print ("Two external argument is given: \n")
    print ("Input OUTCAR file: {0}".format(sys.argv[1]))
    print ("Input EIGENVAL file: {0}\n".format(sys.argv[2]))
    print ("Output file: Xcrysden.bxsf")
    if not os.path.isfile(sys.argv[1]):
        print ("\n** ERROR: Input file {0} was not found.".format(sys.argv[1]))
        sys.exit(0)
    if not os.path.isfile(sys.argv[2]):
        print ("\n** ERROR: Input file {0} was not found.".format(sys.argv[2]))
        sys.exit(0)
    outcar_filename=sys.argv[1]
    in_filename=sys.argv[2]
    out_filename='Xcrysden.bxsf'
elif len(sys.argv) == 4:
    print ("Three external argument is given: \n")
    print ("Input OUTCAR file: {0}".format(sys.argv[1]))
    print ("Input EIGENVAL file: {0}\n".format(sys.argv[2]))
    print ("Output file: {0}".format(sys.argv[3]))
    if not os.path.isfile(sys.argv[1]):
        print ("\n** ERROR: Input file {0} was not found.".format(sys.argv[1]))
        sys.exit(0)
    if not os.path.isfile(sys.argv[2]):
        print ("\n** ERROR: Input file {0} was not found.".format(sys.argv[2]))
        sys.exit(0)
    outcar_filename=sys.argv[1]
    in_filename=sys.argv[2]
    out_filename=sys.argv[3]
    
############  READ OUTCAR FILE  ############################
outcar_file=open(outcar_filename,'r')
lines=outcar_file.readlines()
#Initialize
fermi_energy=0.0
reci_a=[]
reci_b=[]
reci_c=[]

for i,line in enumerate(lines): #starts from 8th line
    if 'E-fermi' in line:
        temp=re.findall('[-a-zA-Z0-9.\+]+',line)
        fermi_energy=float(temp[1])
    if 'reciprocal lattice vectors' in line:
        temp=re.findall('[-a-zA-Z0-9.\+]+',lines[i+1])
        reci_a=np.array([float(temp[-3]),float(temp[-2]),float(temp[-1])])
        temp=re.findall('[-a-zA-Z0-9.\+]+',lines[i+2])
        reci_b=np.array([float(temp[-3]),float(temp[-2]),float(temp[-1])])
        temp=re.findall('[-a-zA-Z0-9.\+]+',lines[i+3])
        reci_c=np.array([float(temp[-3]),float(temp[-2]),float(temp[-1])])
       
if fermi_energy*len(reci_a)*len(reci_b)*len(reci_c) == 0:
    print ("** ERROR: OUTCAR file is not properly read. Please check the file.")
    sys.exit(0)
        
outcar_file.close()
###########################################################

###########  READ EIGENVAL FILE ###########################
#filename='EIGENVAL'
in_file=open(in_filename, 'r')
lines=in_file.readlines()
try:    
    # Information: ispin
    [dump1,dump2,dump3,ispin]=re.findall('[0-9]+',lines[0])
    name_of_system=lines[4].strip()
    # Information : Number of electrons, number of k-points, number of bands
    [nelect,nkpoints,nbands]=re.findall('[0-9]+',lines[5])
    # Convert to integers
    ispin = int(ispin)
    nelect=int(nelect)
    nkpoints=int(nkpoints)
    nbands=int(nbands)
except ValueError:
        print ("** ERROR: EIGENVAL file is not properly read. Please check the file.")
        sys.exit(0)

#Data arrays:
# kpoints: col1--k-point of a-axis, col2--k-point of b-axis, col3--k-point of c-axis
# col4-dump
# eigenvalues: 3-dimensional array. eigenvaluse[i-th kpoint,j-th band,k-th spin]
# eigenvalues[:,j,k] will have 1-D array matched to kpoints array with j-th band with k-spin  
kpoints=np.zeros(shape=(nkpoints,4))
eigenvalues=np.zeros(shape=(nkpoints,nbands,ispin))

# EIGENVAL file is format blocked by k-point. To count the block,
# kp_counter is introduced.
kp_counter=0


# For loop goes to line by line of EIGENVAL file (afer 8th line), and record kpoints and eigenvalues
# depending on how many values are given. If # values are 4, it is for k_point,
# otherwise it is eigen value file.
for i,line in enumerate(lines[7::]): #starts from 8th line
    # If the line is for k-point mash
    if len(re.findall("[-a-zA-Z0-9.\+]+",line))==4:
        temp=re.findall('[-a-zA-Z0-9.\+]+',line)
        num_temp=[float(j) for j in temp]
        kpoints[kp_counter,:]=num_temp
        kp_counter+=1
        
    # If the line is for Eigenvalues
    elif len(re.findall("[-a-zA-Z0-9.\+]+",line)) > 0:
        temp=re.findall('[-a-zA-Z0-9.\+]+',line)
        band_num=int(temp[0])        
        ei_val_string=temp[1:1+ispin:]
        ei_val=[float(k) for k in ei_val_string] #Eigenvlue of current Kpoint
        eigenvalues[kp_counter-1,band_num-1,:]=ei_val
                
in_file.close()
###############################################

#############################################
####### Making Xcrysden file ################
nkx=len(set(kpoints[:,0]))
nky=len(set(kpoints[:,1]))
nkz=len(set(kpoints[:,2]))
#### PART 1: HEADER##################
out_file=open(out_filename,'w')
out_file.write('BEGIN_INFO\n')
out_file.write('   #\n')
out_file.write('   # Case:  {0}\n'.format(name_of_system))
out_file.write('   #\n')
out_file.write('   # Launch as: xcrysden --bxsf example.bxsf\n')
out_file.write('   #\n')
out_file.write('   Fermi Energy: {0}\n'.format(fermi_energy))
out_file.write(' END_INFO\n')
out_file.write('\n')
out_file.write(' BEGIN_BLOCK_BANDGRID_3D\n')
out_file.write(' Num_bands_are_sum_of_spin_up/down._Change_reci_dimension_based_on_your_unit_cell\n')
out_file.write('   BEGIN_BANDGRID_3D\n')
num_total_bands=eigenvalues.shape[1]*eigenvalues.shape[2]
out_file.write('       {0}\n     {1} {2} {3}\n'.format(num_total_bands,nkx,nky,nkz))
out_file.write('     0.0000 0.0000 0.0000\n')
out_file.write('     {0:.4f} {1:.4f} {2:.4f}\n'.format(reci_a[0],reci_a[1],reci_a[2]))
out_file.write('     {0:.4f} {1:.4f} {2:.4f}\n'.format(reci_b[0],reci_b[1],reci_b[2]))
out_file.write('     {0:.4f} {1:.4f} {2:.4f}\n'.format(reci_c[0],reci_c[1],reci_c[2]))

###### PART 2: BAND BLOCKS ###########
for j in range(eigenvalues.shape[1]): # Axis1: For each band
    out_file.write('   BAND:  {0}\n'.format(j+1)) # One band block
    # Convert eigenvalues of one band into array of (nkx*nky, nkz) matrix
    one_band=np.reshape(eigenvalues[:,j,0], (nkx*nky,nkz)) # index 0 is for first spin element
    # block_counter is index of y-z block.
    block_counter=0
    # there are 'nkx' number of blocks within one BAND
    for i in range(nkx):
        # there are 'nky' number of lines within one block
        for k in range(nky):
            # there are 'nkz' number of columns within one line
            # block_counter*nkx+k row will be line to be printed
            out_file.write('       '+'  '.join(map(str,one_band[block_counter*nkx+k,:]))+'\n')
        # each blocks are separated by empty line
        out_file.write('\n')
        block_counter+=1

#For spin_down
if np.shape(eigenvalues)[2]==2: # Axis2: When more than one spin element exists
    for j in range(eigenvalues.shape[1]) :# Axis1: For each band
        out_file.write('   BAND:  {0}\n'.format(j+1+eigenvalues.shape[1])) # One band block
        
        one_band=np.reshape(eigenvalues[:,j,1], (nkx*nky,nkz)) # index 1 is for second spin element
        # block_counter is index of y-z block.
        block_counter=0
        # there are 'nkx' number of blocks within one BAND
        for i in range(nkx):
            # there are 'nky' number of lines within one block
            for k in range(nky):
                # there are 'nkz' number of columns within one line
                # block_counter*nkx+k row will be line to be printed
                out_file.write('       '+'  '.join(map(str,one_band[block_counter*nkx+k,:]))+'\n')
            # each blocks are separated by empty line
            out_file.write('\n')
            block_counter+=1

###### PART 3: TAILS  ####################
out_file.write('   END_BANDGRID_3D\n')
out_file.write(' END_BLOCK_BANDGRID_3D')
out_file.close()

print ("Done.")

