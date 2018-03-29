# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 23:07:02 2017

@author: Yongjin

How to use:
    $ python Fermi_kpoints.py KX KY KZ
    KX, KY, and KZ are number of points along each reciprocal axis
# -*- coding: utf-8 -*-
"""


import numpy as np
import sys

print ("\n### KPOINTS script: Generates KPOINTS file for Fermi-surface calcaultion ")
print ("### Author: Yongjin Shin\n")

#Case of different arguments are given.
if len(sys.argv) == 1:
    kx=9
    ky=9
    kz=9
    print("No input argument is given. KX,KY,KZ=9")
elif len(sys.argv) == 2:
    kx=int(sys.argv[1])
    ky=int(sys.argv[1])
    kz=int(sys.argv[1])
    print("One input argument is given. KX,KY,KZ={0}".format(kx))
elif len(sys.argv) == 3:
    kx=int(sys.argv[1])
    ky=int(sys.argv[2])
    kz=int(sys.argv[2])
    print("Two input argument is given. KX={0}, KY={1}, KZ={1} ".format(kx,ky))
elif len(sys.argv) == 4:
    kx=int(sys.argv[1])
    ky=int(sys.argv[2])
    kz=int(sys.argv[3])
    print("Three input argument is given. KX={0}, KY={1}, KZ={2} ".format(kx,ky,kz))



# Reciprocal points are described as (kx+1) style as Xcrysden does not accept periodic
# boundary conditions
data_array=np.ones(shape=((kx+1)*(ky+1)*(kz+1),4));
for i in range(kx+1):
    for j in range(ky+1):
        for k in range(kz+1):
            data_array[i*(ky+1)*(kz+1)+j*(kz+1)+k,:-1]=[i/kx, j/ky, k/kz]

out_file=open('KPOINTS','w')
out_file.write("k-points for fermi-surface. RP-phase {0}x{1}x{2}\n".format(kx,ky,kz))
out_file.write("{0}\n".format((kx+1)*(ky+1)*(kz+1)))
out_file.write("Reciprocal\n")
for i in range(np.shape(data_array)[0]):
    out_file.write("{0:.4f} {1:.4f} {2:.4f} {3:.0f}\n".format(data_array[i,0],data_array[i,1],data_array[i,2],data_array[i,3]))
out_file.close()

print ("Done.")

