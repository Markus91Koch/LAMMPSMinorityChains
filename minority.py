import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
plt.switch_backend('Qt4Agg')
import seaborn as sns; #sns.set()
import pandas as pd
import re
from collections import Counter
import glob, os
from itertools import groupby
from operator import itemgetter
import sys

#print('Number of arguments:', len(sys.argv), 'arguments.')
#print('Argument List:', str(sys.argv))
#print(int(sys.argv[1]))

os.chdir(".")

fdata=np.loadtxt("atoms.txt", skiprows=0)
gdata=np.loadtxt("footer.txt", skiprows=0)

N=64
Nm=int(sys.argv[1])
Nmiddle=int(Nm/2)
print("Nmiddle:", Nmiddle)


Ntot=len(fdata)
Nbond=len(gdata)

#print(Ntot)
#exit()

brush=fdata[fdata[:,2]!=4] # exclude wall particles (type 4)
wall = fdata[fdata[:,2]==4]


#print(brush)
#print(np.abs(brush[:,3:5]))
#print(np.min(np.abs(brush[:,3])))
#print(np.min(np.abs(brush[:,4])))

#brush[:,3].index(np.min(np.abs(brush[:,3])))
xelem = np.where(np.abs(brush[:,3]) == np.min(np.abs(brush[:,3]))) # find most central chain in x and y
yelem = np.where(np.abs(brush[:,4]) == np.min(np.abs(brush[:,4])))

#print(xelem, yelem)


#exit()

common = np.intersect1d(xelem, yelem)   #
#print(list(common))
#print(fdata[common[:],:])
print(len(common))
print(fdata[common[0], 1])
#exit()

Nm_start = common[0]            # start index of chain
Nm_eq = common[-1]              # end index of chain

start_index = int(brush[common[0],0])
end_index = int(brush[common[-1],0])

brush[common[0],2]=5                    # change type of head
#brush[common[-1],2]=7                   # change type of end
brush[common[-1],2]=6
brush[common[0]+1:common[-1],2]=6       # change type of mid

if Nmiddle<=N:
    brush[common[0]+Nmiddle-1,2]=8

x_last_mon = brush[common[-1],3]
y_last_mon = brush[common[-1],4] 
z_last_mon = brush[common[-1],5]

molnr = int(brush[common[0],1])         # find molecule nr
#print(molnr)

#### add monomers at the end if needed!!!!

if Nm>N:
    Nadd=Nm-N

    newmon=[]
    for i in range(Nadd):
        newmon.append([Ntot+i+1, molnr, 6, x_last_mon, y_last_mon, z_last_mon+i+1])

    new = np.array(newmon)
    new[-1,2]=7

    if Nmiddle>N:
        #print("Nmiddle-N=", Nmiddle-N)
        new[Nmiddle-N-1,2]=8

    #print(new)
    old=np.concatenate((brush, wall), axis=0)
    full=np.concatenate((old, new), axis=0)

    #print(full)

    #### now adjust bonds
    lastbond=int(gdata[-1,0])    #last bond index
    newb=[]
    for i in range(Nadd):
        if i==0:
            newb.append([lastbond+i+1, 1, end_index, Ntot+i+1])  
        else:
            newb.append([lastbond+i+1, 1, Ntot+i, Ntot+i+1])
    newbonds=np.array(newb)
    fullbonds=np.concatenate((gdata, newb), axis=0)
    #print(newb)
else:
    brush[common[-1],2]=7
    full=np.concatenate((brush, wall), axis=0)
    fullbonds=gdata    

## WRITE OUT FILES

#print(full.shape)
#print(int(full[5,0]))
### first of all atom coords
file = open("newatoms.txt", "w")
#file.write("# Nm\tEpsilon_crit\tUbarrier\tFmin\tFmax\tFslope\n")
for i in range(len(full)):
    file.write("%d\t%d\t%d\t%f\t%f\t%f\n" % (int(full[i,0]), int(full[i,1]), int(full[i,2]), full[i,3], full[i,4], full[i,5]))
file.close()

file = open("newfooter.txt", "w")
#file.write("# Nm\tEpsilon_crit\tUbarrier\tFmin\tFmax\tFslope\n")
for i in range(len(fullbonds)):
    file.write("%d\t%d\t%d\t%d\n" % (int(fullbonds[i,0]), int(fullbonds[i,1]), int(fullbonds[i,2]), int(fullbonds[i,3])))
file.close()

file = open("newinfo.txt", "w")
file.write("%d\n" % len(full))
file.write("%d\n" % len(fullbonds))
file.close()


