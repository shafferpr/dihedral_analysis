#! /usr/bin/env python

import os
import sys
import math
import string
import numpy as np
from dihedral_analysis_tools import Readcolvar
from dihedral_analysis_tools import Readcrystal

N_CVs=int(sys.argv[1])
cv_phi_file=str(sys.argv[2])
cv_psi_file=str(sys.argv[3])
phi_crystal_file=str(sys.argv[4])
psi_crystal_file=str(sys.argv[5])
histofile1=str(sys.argv[6])
histofile2=str(sys.argv[7])

dihedrals_phi=Readcolvar(cv_phi_file, N_CVs)
dihedrals_psi=Readcolvar(cv_psi_file, N_CVs)

phi_crystal=Readcrystal(phi_crystal_file, N_CVs)
psi_crystal=Readcrystal(psi_crystal_file, N_CVs)


print dihedrals_phi.shape[0]
print dihedrals_psi.shape
for i in range(dihedrals_phi.shape[0]):
    for j in range(dihedrals_phi.shape[1]):
        x=dihedrals_phi[i][j]-phi_crystal[j]
        if(x>3.14159):
            x=x-2*3.14159
        if(x<-3.14159):
            x=2*3.14159+x
        dihedrals_phi[i][j]=x
        
        y=dihedrals_psi[i][j]-psi_crystal[j]
        if(y>3.14159):
            y=y-2*3.14159
        if(y<-3.14159):
            y=2*3.14159+y
        dihedrals_psi[i][j]=y

filename="deltapsi.data"
outputfile=open(filename,"a")
for i in range(dihedrals_psi.shape[0]):
    for j in range(dihedrals_psi.shape[1]):
        outputfile.write('%f ' %(dihedrals_psi[i][j]))
    outputfile.write('\n')

rmsd1=0
rmsd2=0

histo1=np.zeros(101);
histo2=np.zeros(101);
histoboth=np.zeros([101,101]);

for i in range(dihedrals_phi.shape[0]):
    rmsd1=0
    for j in range(0, 4, 1):
        rmsd1+=dihedrals_phi[i][j]*dihedrals_phi[i][j]+dihedrals_psi[i][j]*dihedrals_psi[i][j]
    rmsd1=math.sqrt(rmsd1/8)
    rmsd2=0
    for j in range(4, 9, 1):
        rmsd2+=dihedrals_phi[i][j]*dihedrals_phi[i][j]+dihedrals_psi[i][j]*dihedrals_psi[i][j]

    rmsd2=math.sqrt(rmsd2/10)
    histo1[rmsd1*100/3.14159]+=1
    histo2[rmsd2*100/3.14159]+=1
    histoboth[rmsd1*100/3.14159][rmsd2*100/3.14159]+=1

filename=histofile1
outputfile=open(filename, "w")
for i in range(histo1.size):
    outputfile.write('%f %f %f\n' %(i*3.14159/100, histo1[i], histo2[i]))
outputfile.close

filename=histofile2
outputfile=open(filename, "w")
for i in range(histoboth.shape[0]):
    for j in range(histoboth.shape[1]):
        outputfile.write('%f %f %f\n' %(i*3.14159/100, j*3.14159/100, histoboth[i][j]))
outputfile.close
