#!/usr/bin/env python2.7
import sys
import numpy as np
from MDAnalysis import *
import MDAnalysis
import MDAnalysis.lib.distances # Previously, this was MDAnalysis.core.distances
import numpy.linalg
import scipy.stats
import matplotlib.pyplot as plt
import math
import MDAnalysis.lib.NeighborSearch as NS
from scipy.spatial import Voronoi, voronoi_plot_2d
import multiprocessing as mp


print 'Initiating Voroni Tesselation'
#Lpd1 is Phospolipid1 
#Lpd2 is Phospolipid 2 
#Lpd3 is Phospolipid 3 or Cholesterol
 
###############INPUTS######################################

#Indexs for atom selection are according to the psf file

side     = sys.argv[1] # up for upper leaflet, down for lower leaflet

#Voronoi Snapshots
PRINT_VOR    = 'NO'
PRINT_FREQ   = 1
PROT_color   = 'yellow'


# MDAnalysisInputs
print "Loading inputs..."

psf  = '../../step5_assembly.psf'
trr  = '../prod.run.trr'

u = MDAnalysis.Universe(psf,trr)

#Lipid Residue names
lipid1 ='DPPC'
lipid2 ='DIPC'
lipid3 ='CHOL'

#Specify the qasi-2D plne for the membrane tesselation
plane_of_interest = 'TAIL' #HEAD, TAIL

#Selcts proteins in the analysis
protein = 'YES' # YES, NO
selprot = 'segid PRO* and name BB'  
prot_res = np.array(('ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC'))
PROT_seg  = np.array(('PROA','PROB'))

# Identify number of residues in each lipid and extract only he top residues (for now)

#Frames to be calculated
end_f = u.trajectory.n_frames
print end_f
start_f = 1
skip    = 1

# Number of processors to use in multiprocessing
nprocs = 28

frames = np.arange(start_f, end_f)[::skip]
n_frames = len(frames)

######################################################################################

# Make empty arrays to hold contact counts
ens_Lpd1_Lpd1 = np.zeros(n_frames)
ens_Lpd2_Lpd2 = np.zeros(n_frames)
ens_Lpd3_Lpd3 = np.zeros(n_frames)
ens_Lpd1_Lpd2 = np.zeros(n_frames)
ens_Lpd1_Lpd3 = np.zeros(n_frames)
ens_Lpd2_Lpd3 = np.zeros(n_frames)
ens_Prot_Lpd1 = np.zeros(n_frames)
ens_Prot_Lpd2 = np.zeros(n_frames)
ens_Prot_Lpd3 = np.zeros(n_frames)
ens_Prot_Prot = np.zeros(n_frames)

######################################################################################
def voronoi_tessel(frame):
    # set the time step
    print 'Frame %i in %i'%(frame, end_f)
    u = MDAnalysis.Universe(psf,trr)
    u.trajectory[int(frame)]
    #u.trajectory[frame-1] # What the hell is MDAnalysis doing...? This changes the frame to frame "ts"

# Select atoms within this particular frame
    lpd1_atoms = u.select_atoms('resname %s and name PO4'%lipid1) # previously, this was u.selectAtoms
    lpd2_atoms = u.select_atoms('resname %s and name PO4'%lipid2)
    lpd3_atoms = u.select_atoms('resname %s and name ROH'%lipid3)
	
    prot_atoms     = u.select_atoms(selprot) 
    
    num_lpd1 = lpd1_atoms.n_atoms
    num_lpd2 = lpd2_atoms.n_atoms

    # atoms in the upper leaflet as defined by insane.py or the CHARMM-GUI membrane builders
    # select cholesterol headgroups within 1.5 nm of lipid headgroups in the selected leaflet
    if side == 'up':
        lpd1i = lpd1_atoms[:((num_lpd1)/2)]
        lpd2i = lpd2_atoms[:((num_lpd2)/2)]

        lipids = lpd1i + lpd2i
        ns_lipids = NS.AtomNeighborSearch(lpd3_atoms)

        lpd3i = ns_lipids.search(lipids,15.0)
    elif side == 'down':
        lpd1i = lpd1_atoms[((num_lpd1)/2):]
        lpd2i = lpd2_atoms[((num_lpd2)/2):]

        lipids = lpd1i + lpd2i
        ns_lipids = NS.AtomNeighborSearch(lpd3_atoms)

        lpd3i = ns_lipids.search(lipids,15.0)

    if plane_of_interest == 'TAIL' :
        #Using tail definitions
        groups1 = []
        for i in np.arange(len(lpd1i.resnums)):
            resnum = lpd1i.resnums[i]
            group = u.select_atoms('resname %s and resnum %i and (name C2A or name C2B)'%(lipid1,resnum))
            groups1.append(group)
        lpd1i = np.sum(groups1)
        
        groups2 = []
        for i in np.arange(len(lpd2i.resnums)):
            resnum = lpd2i.resnums[i]
            group = u.select_atoms('resname %s and resnum %i and (name D2A or name D2B)'%(lipid2,resnum))
            groups2.append(group)
        lpd2i = np.sum(groups2)

        groups3 = []
        for i in np.arange(len(lpd3i.resnums)):
            resnum = lpd3i.resnums[i]
            group = u.select_atoms('resname CHOL and resnum %i and (name R3)'%resnum)
            groups3.append(group)
        lpd3i = np.sum(groups3)

    if protein == 'YES' :    
    #Selecting Protein atoms that are on the membrane plane.
        phospholipidgrid  = lpd1i + lpd2i #Phospholipid tail atoms that defines the plane of interest.
        ns_PROT           = NS.AtomNeighborSearch(prot_atoms) #Protein
        proti             = ns_PROT.search(phospholipidgrid,10.0) #Protein atoms that are within 1.0nm of the membrane   
        u_leaflet     = lpd1i + lpd2i + lpd3i + proti
		
    else :
        u_leaflet     = lpd1i + lpd2i + lpd3i	


#Extracting the coordinates
    lpd_atms = u_leaflet
    Pxyz     = lpd_atms.positions
    Pxy      = []
    for l in range(0,len(Pxyz)) :
        Pxy.append([Pxyz[l][0],Pxyz[l][1]])
#Extracting xy coordinates and residue names

    atm_list = []
    for a in range(0, len(Pxyz)):
        #print lpd_atms[a].resname 
        atm_list.append([Pxyz[a][0],Pxyz[a][1],lpd_atms[a].resname,lpd_atms[a].segid])

#Introducing PBC
    x_box = u.dimensions[0]
    y_box = u.dimensions[1]

    xplus   = []
    xminus  = []
    xyplus  = []
    xyminus = []

    for atm in range(0 ,len(atm_list)):
        xplus.append([atm_list[atm][0]+x_box,atm_list[atm][1],atm_list[atm][2],atm_list[atm][3]])
        xminus.append([atm_list[atm][0]-x_box,atm_list[atm][1],atm_list[atm][2],atm_list[atm][3]])

    atm_list_px = atm_list + xplus + xminus

    for atm in range(0 ,len(atm_list_px)):
        xyplus.append([atm_list_px[atm][0],atm_list_px[atm][1]+y_box,atm_list_px[atm][2],atm_list_px[atm][3]])
        xyminus.append([atm_list_px[atm][0],atm_list_px[atm][1]-y_box,atm_list_px[atm][2],atm_list_px[atm][3]])

    atm_list_p = atm_list_px + xyplus + xyminus


    atm_xy = []
    for i in range(0,len(atm_list_p)) :
        atm_xy.append([atm_list_p[i][0],atm_list_p[i][1]])


    vor = Voronoi(atm_xy)
    vor_s = Voronoi(Pxy)
    vertices       = vor.vertices

    ridge_points = vor.ridge_points
    
    #Plotting Voroni Diagrams
    vor_rgns     = vor.regions
    l_vor_rgns   = len(vor_rgns)

    vor_points   = vor.point_region
    l_vor_points = len(vor_points)    
    if PRINT_VOR == 'YES' and np.mod(frame,PRINT_FREQ) == 0 :
        plt.clf()
        #voronoi_plot_2d(vor_s)
    
        for p in range(0, l_vor_points):
            rgn = vor_rgns[vor_points[p]]
            L = atm_list_p[p]
            if not -1 in rgn and 0 <= L[0] < x_box and 0 <= L[1] < y_box :
                if L[2] == lipid1:
                    polygon = [vor.vertices[i] for i in rgn]
                    plt.plot(*zip(*polygon), color='black')
                    plt.fill(*zip(*polygon), color='blue')
                elif L[2] == lipid2:
                    polygon = [vor.vertices[i] for i in rgn]
                    plt.plot(*zip(*polygon), color='black')
                    plt.fill(*zip(*polygon), color='red')
                elif L[2] == lipid3:
                    polygon = [vor.vertices[i] for i in rgn]
                    plt.plot(*zip(*polygon), color='black')
                    plt.fill(*zip(*polygon), color='green')				
                else :
                    polygon = [vor.vertices[i] for i in rgn]
                    plt.plot(*zip(*polygon), color='black')
                    plt.fill(*zip(*polygon), color=PROT_color)
        plt.axis('equal')
        if side == 'up':
            plt.savefig('VT/img' + str('%03d' %frame) + '.up.png')
        if side == 'down':
            plt.savefig('VT/img' + str('%03d' %frame) + '.down.png')    
    
    
    Lpd1_Lpd1_I = 0
    Lpd2_Lpd2_I = 0
    Lpd3_Lpd3_I = 0
    Lpd1_Lpd2_I = 0
    Lpd1_Lpd3_I = 0
    Lpd2_Lpd3_I = 0

    Lpd1_Lpd1_E = 0
    Lpd2_Lpd2_E = 0
    Lpd3_Lpd3_E = 0
    Lpd1_Lpd2_E = 0
    Lpd1_Lpd3_E = 0
    Lpd2_Lpd3_E = 0
    
    Prot_Lpd1_I = 0
    Prot_Lpd2_I = 0
    Prot_Lpd3_I = 0
    Prot_Prot_I = 0    

    Prot_Lpd1_E = 0
    Prot_Lpd2_E = 0
    Prot_Lpd3_E = 0
    Prot_Prot_E = 0   
    
    r_length  = len(ridge_points)


  
    for k in range (0,r_length) :
        ridge_k = ridge_points[k]
        Li = atm_list_p[int(ridge_k[0])]
        Lj = atm_list_p[int(ridge_k[1])]
        
        #print 'Li', Li[3], Li[2]
        #print 'Lj', Lj[3], Lj[2] 
                
#Lipids INSIDE the box 

        if 0 < Li[0] < x_box and 0 < Li[1] < y_box and 0 < Lj[0] < x_box and 0 < Lj[1] < y_box :
#Lipid Lipid contacts 
    
            if Li[2] == lipid1 and Lj[2] == lipid1:
                Lpd1_Lpd1_I = Lpd1_Lpd1_I + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid2:
                Lpd2_Lpd2_I = Lpd2_Lpd2_I + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid3:
                Lpd3_Lpd3_I = Lpd3_Lpd3_I + 1
                
            if Li[2] == lipid1 and Lj[2] == lipid2:
                Lpd1_Lpd2_I  = Lpd1_Lpd2_I + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid1:
                Lpd1_Lpd2_I = Lpd1_Lpd2_I + 1

            if Li[2] == lipid1 and Lj[2] == lipid3:
                Lpd1_Lpd3_I  = Lpd1_Lpd3_I + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid1:
                Lpd1_Lpd3_I = Lpd1_Lpd3_I + 1

            if Li[2] == lipid2 and Lj[2] == lipid3:
                Lpd2_Lpd3_I  = Lpd2_Lpd3_I + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid2:
                Lpd2_Lpd3_I = Lpd2_Lpd3_I + 1 
                
#Protein lipid contacts
            if Li[3] in PROT_seg and Lj[2] == lipid1:
                Prot_Lpd1_I = Prot_Lpd1_I + 1

            if Li[2] == lipid1 and Lj[3] in PROT_seg:
                Prot_Lpd1_I = Prot_Lpd1_I + 1
                
            if Li[3] in PROT_seg and Lj[2] == lipid2:
                Prot_Lpd2_I = Prot_Lpd2_I + 1
                
            if Li[2] == lipid2 and Lj[3] in PROT_seg:
                Prot_Lpd2_I = Prot_Lpd2_I + 1                
                
            if Li[3] in PROT_seg and Lj[2] == lipid3:
                Prot_Lpd3_I = Prot_Lpd3_I + 1
                
            if Li[2] == lipid3 and Lj[3] in PROT_seg :
                Prot_Lpd3_I = Prot_Lpd3_I + 1                
                
            if Li[3] in PROT_seg and Lj[3] in PROT_seg and Li[3] != Lj[3]: # Avoid self counts
                Prot_Prot_I = Prot_Prot_I + 1                                               
                
#Lipids at the EDGE of the box                
#Lipid Lipid contacts                
        if 0 <= Li[0] < x_box and 0 <= Li[1] < y_box or 0 <= Lj[0] < x_box and 0 <= Lj[1] < y_box :

            if Li[2] == lipid1 and Lj[2] == lipid1:
                Lpd1_Lpd1_E = Lpd1_Lpd1_E + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid2:
                Lpd2_Lpd2_E = Lpd2_Lpd2_E + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid3:
                Lpd3_Lpd3_E = Lpd3_Lpd3_E + 1
                
            if Li[2] == lipid1 and Lj[2] == lipid2:
                Lpd1_Lpd2_E  = Lpd1_Lpd2_E + 1
                
            if Li[2] == lipid2 and Lj[2] == lipid1:
                Lpd1_Lpd2_E = Lpd1_Lpd2_E + 1

            if Li[2] == lipid1 and Lj[2] == lipid3:
                Lpd1_Lpd3_E  = Lpd1_Lpd3_E + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid1:
                Lpd1_Lpd3_E = Lpd1_Lpd3_E + 1

            if Li[2] == lipid2 and Lj[2] == lipid3:
                Lpd2_Lpd3_E  = Lpd2_Lpd3_E + 1
                
            if Li[2] == lipid3 and Lj[2] == lipid2:
                Lpd2_Lpd3_E = Lpd2_Lpd3_E + 1    

#Protein lipid contacts
            if Li[3] in PROT_seg and Lj[2] == lipid1:
                Prot_Lpd1_E = Prot_Lpd1_E + 1

            if Li[2] == lipid1 and Lj[3] in PROT_seg:
                Prot_Lpd1_E = Prot_Lpd1_E + 1
                
            if Li[3] in PROT_seg and Lj[2] == lipid2:
                Prot_Lpd2_E = Prot_Lpd2_E + 1
                
            if Li[2] == lipid2 and Lj[3] in PROT_seg:
                Prot_Lpd2_E = Prot_Lpd2_E + 1                
                
            if Li[3] in PROT_seg and Lj[2] == lipid3:
                Prot_Lpd3_E = Prot_Lpd3_E + 1
                
            if Li[2] == lipid3 and Lj[3] in PROT_seg :
                Prot_Lpd3_E = Prot_Lpd3_E + 1                
                
            if Li[3] in PROT_seg and Lj[3] in PROT_seg and Li[3] != Lj[3]: # Avoid self counts
                Prot_Prot_E = Prot_Prot_E + 1     
                
#Total = LipidsInside + (Lipids including EDGES - Lipids Inside)/2 -----> Correction for over counting the lipids in periodic images
    Lpd1_Lpd1 = Lpd1_Lpd1_I + (Lpd1_Lpd1_E - Lpd1_Lpd1_I)/2
    Lpd2_Lpd2 = Lpd2_Lpd2_I + (Lpd2_Lpd2_E - Lpd2_Lpd2_I)/2
    Lpd3_Lpd3 = Lpd3_Lpd3_I + (Lpd3_Lpd3_E - Lpd3_Lpd3_I)/2
    Lpd1_Lpd2 = Lpd1_Lpd2_I + (Lpd1_Lpd2_E - Lpd1_Lpd2_I)/2
    Lpd1_Lpd3 = Lpd1_Lpd3_I + (Lpd1_Lpd3_E - Lpd1_Lpd3_I)/2
    Lpd2_Lpd3 = Lpd2_Lpd3_I + (Lpd2_Lpd3_E - Lpd2_Lpd3_I)/2
    
    Prot_Lpd1 = Prot_Lpd1_I + (Prot_Lpd1_E - Prot_Lpd1_I)/2
    Prot_Lpd2 = Prot_Lpd2_I + (Prot_Lpd2_E - Prot_Lpd2_I)/2
    Prot_Lpd3 = Prot_Lpd3_I + (Prot_Lpd3_E - Prot_Lpd3_I)/2
    Prot_Prot = Prot_Prot_I + (Prot_Prot_E - Prot_Prot_I)/2       

    return Lpd1_Lpd1, Lpd2_Lpd2, Lpd3_Lpd3, Lpd1_Lpd2, Lpd1_Lpd3, Lpd2_Lpd3, Prot_Lpd1, Prot_Lpd2, Prot_Lpd3, Prot_Prot

pool = mp.Pool(processes=nprocs)
print 'Initiating multiprocessing with %i processors'%nprocs
results = pool.map(voronoi_tessel, frames)

ens_Lpd1_Lpd1_np = []
ens_Lpd2_Lpd2_np = []
ens_Lpd3_Lpd3_np = []
ens_Lpd1_Lpd2_np = []
ens_Lpd1_Lpd3_np = []
ens_Lpd2_Lpd3_np = []
ens_Prot_Lpd1_np = []
ens_Prot_Lpd2_np = []
ens_Prot_Lpd3_np = []
ens_Prot_Prot_np = []


for i in range(n_frames):
    ens_Lpd1_Lpd1_np.append(results[i][0])
    ens_Lpd2_Lpd2_np.append(results[i][1])
    ens_Lpd3_Lpd3_np.append(results[i][2])
    ens_Lpd1_Lpd2_np.append(results[i][3])
    ens_Lpd1_Lpd3_np.append(results[i][4])
    ens_Lpd2_Lpd3_np.append(results[i][5])
    ens_Prot_Lpd1_np.append(results[i][6])
    ens_Prot_Lpd2_np.append(results[i][7])
    ens_Prot_Lpd3_np.append(results[i][8])
    ens_Prot_Prot_np.append(results[i][9])    

ens_Lpd1_Lpd1_np = np.asarray(ens_Lpd1_Lpd1_np)
ens_Lpd2_Lpd2_np = np.asarray(ens_Lpd2_Lpd2_np)
ens_Lpd3_Lpd3_np = np.asarray(ens_Lpd3_Lpd3_np)
ens_Lpd1_Lpd2_np = np.asarray(ens_Lpd1_Lpd2_np)
ens_Lpd1_Lpd3_np = np.asarray(ens_Lpd1_Lpd3_np)
ens_Lpd2_Lpd3_np = np.asarray(ens_Lpd2_Lpd3_np)
ens_Prot_lpd1_np = np.asarray(ens_Prot_Lpd1_np)
ens_Prot_Lpd2_np = np.asarray(ens_Prot_Lpd2_np)
ens_Prot_Lpd3_np = np.asarray(ens_Prot_Lpd3_np)
ens_Prot_Prot_np = np.asarray(ens_Prot_Prot_np)


# Define output file names
if side == "up":
    Lpd1_Lpd1_fn   = 'upper_Lpd1_Lpd1.dat'
    Lpd2_Lpd2_fn   = 'upper_Lpd2_Lpd2.dat'
    Lpd3_Lpd3_fn   = 'upper_Lpd3_Lpd3.dat'
    Lpd1_Lpd2_fn   = 'upper_Lpd1_Lpd2.dat'
    Lpd1_Lpd3_fn   = 'upper_Lpd1_Lpd3.dat'
    Lpd2_Lpd3_fn   = 'upper_Lpd2_Lpd3.dat'
    Prot_Lpd1_fn   = 'upper_Prot_Lpd1.dat'
    Prot_Lpd2_fn   = 'upper_Prot_Lpd2.dat'
    Prot_Lpd3_fn   = 'upper_Prot_Lpd3.dat'
    Prot_Prot_fn   = 'upper_Prot_Prot.dat'    

elif side == "down":
    Lpd1_Lpd1_fn   = 'lower_Lpd1_Lpd1.dat'
    Lpd2_Lpd2_fn   = 'lower_Lpd2_Lpd2.dat'
    Lpd3_Lpd3_fn   = 'lower_Lpd3_Lpd3.dat'
    Lpd1_Lpd2_fn   = 'lower_Lpd1_Lpd2.dat'
    Lpd1_Lpd3_fn   = 'lower_Lpd1_Lpd3.dat'
    Lpd2_Lpd3_fn   = 'lower_Lpd2_Lpd3.dat'
    Prot_Lpd1_fn   = 'lower_Prot_Lpd1.dat'
    Prot_Lpd2_fn   = 'lower_Prot_Lpd2.dat'
    Prot_Lpd3_fn   = 'lower_Prot_Lpd3.dat'
    Prot_Prot_fn   = 'lower_Prot_Prot.dat' 

#Writing Outputs
np.savetxt(Lpd1_Lpd1_fn,ens_Lpd1_Lpd1_np)
np.savetxt(Lpd2_Lpd2_fn,ens_Lpd2_Lpd2_np)
np.savetxt(Lpd3_Lpd3_fn,ens_Lpd3_Lpd3_np)
np.savetxt(Lpd1_Lpd2_fn,ens_Lpd1_Lpd2_np)
np.savetxt(Lpd1_Lpd3_fn,ens_Lpd1_Lpd3_np)
np.savetxt(Lpd2_Lpd3_fn,ens_Lpd2_Lpd3_np)
np.savetxt(Prot_Lpd1_fn,ens_Prot_Lpd1_np)
np.savetxt(Prot_Lpd2_fn,ens_Prot_Lpd2_np)
np.savetxt(Prot_Lpd3_fn,ens_Prot_Lpd3_np)
np.savetxt(Prot_Prot_fn,ens_Prot_Prot_np)

print 'Calculation Complete'
