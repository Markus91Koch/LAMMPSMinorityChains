import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import seaborn as sns; #sns.set()
import pandas as pd
import re
from collections import Counter
import glob, os
from itertools import groupby
from operator import itemgetter
import sys

plt.switch_backend('Qt4Agg')

os.chdir(".")

##### read in atoms section and bonds section ######

fdata=np.loadtxt("atoms.txt", skiprows=0)
gdata=np.loadtxt("footer.txt", skiprows=0)


N = 128 # length of brush chains
Nm = int(sys.argv[1]) # final length of minority chain (1st command line input)
Nmiddle = int(Nm/2) # middle monomer position of minority chain

# total number of monomers and boonds so far
Ntot = len(fdata)
Nbond = len(gdata)

# take only atoms of types unequal 4 (brush) or equal 4 (wall)
brush = fdata[fdata[:,2]!=4]
wall = fdata[fdata[:,2]==4]

# find most central chain in x and y
xelem = np.where(np.abs(brush[:,3]) == np.min(np.abs(brush[:,3])))
yelem = np.where(np.abs(brush[:,4]) == np.min(np.abs(brush[:,4])))
common = np.intersect1d(xelem, yelem) 

Nm_start = common[0]            # start index of chain to be modified
Nm_eq = common[-1]              # end index of that chain
start_index = int(brush[common[0],0])
end_index = int(brush[common[-1],0])

brush[common[0],2] = 5                    # change atom type of previous head
brush[common[-1],2] = 6                   # change atom type of previous end
brush[common[0]+1:common[-1],2] = 6       # change atom type of prev. midsection

# if middle monomer is below N, edit in old data 
if Nmiddle <= N:
    brush[common[0]+Nmiddle-1,2] = 8

# get cartesian coords of last monomer of that chain
x_last_mon = brush[common[-1],3]
y_last_mon = brush[common[-1],4] 
z_last_mon = brush[common[-1],5]
molnr = int(brush[common[0],1])         # find its molecule index

#### add monomers at the end of file if minority chain is longer than brush chain

# if minority chain is longer than bulk chains
if Nm > N:

    # number of protruding monomers
    Nadd = Nm - N

    # collect data (coords, types, ..) of newly added monomers
    newmon = []
    for i in range(Nadd):
        newmon.append([Ntot+i+1, molnr, 6, x_last_mon, y_last_mon, z_last_mon+i+1])

    # change atom type of new end monomer
    new = np.array(newmon)
    new[-1,2] = 7

    # atom type of middle monomer
    if Nmiddle > N:
        new[Nmiddle-N-1,2] = 8

    # add new atom coords to old ones
    old = np.concatenate((brush, wall), axis=0)
    full = np.concatenate((old, new), axis=0)

    # generate new bonds info

    lastbond=int(gdata[-1,0])    # index of last bond of "old" chain
    newb = []
    for i in range(Nadd):
        if i == 0:
            newb.append([lastbond+i+1, 1, end_index, Ntot+i+1])  
        else:
            newb.append([lastbond+i+1, 1, Ntot+i, Ntot+i+1])

    #add new bond info to old one
    newbonds = np.array(newb)
    fullbonds = np.concatenate((gdata, newb), axis=0)

# if minority chain is not longer than bulk chains
# just change end monomer and recombine brush and wall data
else:
    brush[common[-1],2] = 7
    full = np.concatenate((brush, wall), axis=0)
    fullbonds = gdata    

#### WRITE OUT FILES ######

# first of all atom coords
file = open("newatoms.txt", "w")
#file.write("# Nm\tEpsilon_crit\tUbarrier\tFmin\tFmax\tFslope\n")
for i in range(len(full)):
    file.write("%d\t%d\t%d\t%f\t%f\t%f\n" % (int(full[i,0]), int(full[i,1]), int(full[i,2]), full[i,3], full[i,4], full[i,5]))
file.close()

# secondly bond section
file = open("newfooter.txt", "w")
#file.write("# Nm\tEpsilon_crit\tUbarrier\tFmin\tFmax\tFslope\n")
for i in range(len(fullbonds)):
    file.write("%d\t%d\t%d\t%d\n" % (int(fullbonds[i,0]), int(fullbonds[i,1]), int(fullbonds[i,2]), int(fullbonds[i,3])))
file.close()

# generate file with info about new total no. of atoms and bonds
file = open("newinfo.txt", "w")
file.write("%d\n" % len(full))
file.write("%d\n" % len(fullbonds))
file.close()


