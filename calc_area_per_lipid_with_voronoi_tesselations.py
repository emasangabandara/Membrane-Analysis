#!/usr/bin/env python2.7
import string
import sys
import numpy as np
from MDAnalysis import *
import MDAnalysis
from scipy.spatial import Voronoi
#from tadlib.polygon import shoelace
import MDAnalysis.lib.NeighborSearch as NS
import multiprocessing as mp
from glob import glob
import sys

#Lpd1 is Phospolipid1 
#Lpd2 is Phospolipid 2 
#Lpd3 is Phospolipid 3 or Cholesterol
side     = sys.argv[1] # up for upper leaflet, down for lower leaflet

###############INPUTS######################################
# MDAnalysisInputs
print "Loading inputs..."

psf  = '../../step5_assembly.psf'
trr  = '../prod.run.trr'

u        = MDAnalysis.Universe(psf,trr)

nprocs   = 28 # number of stupidly parallel processors to use

start_f  = 1
end_f    = u.trajectory.n_frames
skip     = 1
frames   = np.arange(start_f,end_f,skip)

#Lipid Residue names
lipid1 ='DPPC'
lipid2 ='DIPC'
lipid3 ='CHOL'

def aplHist(data) :
    y,binEdges=np.histogram(data,normed=True,bins=25,range=(0,100))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    
    return bincenters , y
    
#Introducing PBC in x and y dimensions
def introducePBC(x_box,y_box,atm_list):
    xplus      = []
    xminus     = []
    xyplus     = []
    xyminus    = []

    for atm in range(0 ,len(atm_list)):
        xplus.append([atm_list[atm][0]+x_box,atm_list[atm][1],atm_list[atm][2]])
        xminus.append([atm_list[atm][0]-x_box,atm_list[atm][1],atm_list[atm][2]])

    atm_list_pbcx = atm_list + xplus + xminus

    for atm in range(0 ,len(atm_list_pbcx)):
        xyplus.append([atm_list_pbcx[atm][0],atm_list_pbcx[atm][1]+y_box,atm_list_pbcx[atm][2]])
        xyminus.append([atm_list_pbcx[atm][0],atm_list_pbcx[atm][1]-y_box,atm_list_pbcx[atm][2]])

    atm_list_pbc = atm_list_pbcx + xyplus + xyminus

    return atm_list_pbc
  
#Polygon area calculation using shoelace algorithm
#def calcAreaofPolygon(vert) :
#    area = np.abs(0.5*shoelace(vert))
#    return area

def calcAreaofPolygon(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

    
#Calculate the coordinates of the Voronoi Vertices
def calcVerticeCoordinates(vert_arr,vor_struct) :
    vert_coor = []
    for vt in vert_arr :
        vert_coor.append(vor_struct.vertices[vt])
    return vert_coor
    
#Calculate Area per lipid
def calcAPL(frame):
    print frame
    u        = MDAnalysis.Universe(psf,trr)
    u.trajectory[frame]

    #Lipid Residue names
    lpd1_atoms = u.select_atoms('resname %s and name PO4'%lipid1) 
    lpd2_atoms = u.select_atoms('resname %s and name PO4'%lipid2)
    lpd3_atoms = u.select_atoms('resname %s and name ROH'%lipid3)


    num_lpd1 = lpd1_atoms.n_atoms
    num_lpd2 = lpd2_atoms.n_atoms

    #Separating upper and lower leaflets
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

    lpd_atms = lpd1i + lpd2i + lpd3i

    nlipid = lpd_atms.n_atoms
    
    frm_l1_lpd_area   = []
    frm_l2_lpd_area   = []
    frm_l3_lpd_area   = []

    #Extracting the coordinates
    Pxyz = lpd_atms.positions
    Pxy = []
    for l in range(0,len(Pxyz)) :
        Pxy.append([Pxyz[l][0],Pxyz[l][1]])
    #Extracting xy coordinates and residue names

    atmlist = []
    for a in range(0, len(Pxyz)):
        atmlist.append([Pxyz[a][0],Pxyz[a][1],lpd_atms[a].resname])

    xbox = u.dimensions[0]
    ybox = u.dimensions[1]

    atm_list_pbc = introducePBC(xbox,ybox,atmlist)

    atm_xy = []
    for i in range(0,len(atm_list_pbc)) :
        atm_xy.append([atm_list_pbc[i][0],atm_list_pbc[i][1]])

    vor       = Voronoi(atm_xy)
    vor_s     = Voronoi(Pxy)
    vertices  = vor.vertices

    #Selecting the Index of the Voronoi region for each input point in siede the box
    point_region_box = vor.point_region[:nlipid]
    point_region_pbc = vor.point_region
    
    for k in range(len(point_region_box)) :
        point_region_k = point_region_box[k]      
        vor_vert       = vor.regions[point_region_k]
        vor_vert_coor  = calcVerticeCoordinates(vor_vert,vor)
        lpd_area       = calcAreaofPolygon(vor_vert_coor)

        if atm_list_pbc[k][2]== lipid1 :
            frm_l1_lpd_area.append(lpd_area)
        elif atm_list_pbc[k][2]== lipid2 :
            frm_l2_lpd_area.append(lpd_area)
        elif atm_list_pbc[k][2]== lipid3 :
            frm_l3_lpd_area.append(lpd_area)    
        else:
            frm_GM_lpd_area.append(lpd_area)   

    return frm_l1_lpd_area , frm_l2_lpd_area , frm_l3_lpd_area , frm_GM_lpd_area

pool = mp.Pool(processes=nprocs)
print 'Initiating multiprocessing with %i processors'%nprocs
results = pool.map(calcAPL, frames)

#Collecting data and saving results
l1_apl    = []
l2_apl    = []
l3_apl    = []

for i in range(len(results)):
    l1_apl.append(results[i][0])
    l2_apl.append(results[i][1])
    l3_apl.append(results[i][2])

l1_apl   = np.array(l1_apl)
l2_apl   = np.array(l2_apl)
l3_apl   = np.array(l3_apl)

if side == 'up':
    np.save('%s_APL_upper.npy'%lipid1,l1_apl)
    np.save('%s_APL_upper.npy'%lipid2,l2_apl)
    np.save('%s_APL_upper.npy'%lipid3,l3_apl)
if side == 'down':
    np.save('%s_APL_lower.npy'%lipid1,l1_apl)
    np.save('%s_APL_lower.npy'%lipid2,l2_apl)
    np.save('%s_APL_lower.npy'%lipid3,l3_apl)
