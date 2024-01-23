### IMPORTS ###


import re
import math
import time
import pickle
import os.path
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.path as mpath
import matplotlib.collections as mcoll
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import axes3d
import pytraj as pt
import mdtraj as md


### CORE ###


def orient(p, p0, pZ, pX):
    '''
    Cambia las coordenadas de "punto" a un nuevo sistema de referencia en el que "p0" es el 0,
    el eje Z está orientado de "p0" a "pZ", y el eje X está orientado de "p0" al punto "pX" proyectado
    sobre el plano perpendicular al nuevo eje Z
    '''
    vZ = pZ - p0 # Vector from p0 to pZ
    uZ = vZ/np.linalg.norm(vZ) # Vector del nuevo eje Z
    
    if pX is None:
        uX = np.array([1.0, 0.0, 0.0])
    else:
        wX = pX - p0 # Vector de p0 a pX
        vX = wX - np.dot(wX, uZ)*uZ # wX proyectado sobre el plano perpendicular a uZ
        uX = vX/np.linalg.norm(vX) # Vector del nuevo eje X
    
    uY = np.cross(uZ, uX) # Vector del nuevo eje Y
    
    vp = p - p0 # Vector de p0 a punto
    orientado = np.array([np.dot(vp, uX), np.dot(vp, uY), np.dot(vp, uZ)])
    return orientado
#end


def recenter_traj_RMSD(run_name, N_tubes, N_res):
    '''
    Vamos a tomar:
    como p0 el centro de masa de los carbonos alfas de los 4 tubos
    como pZ el centro de masa de los carbonos alfas de los primeros anillos de los nanotubos
    como pX el centro de masa de los carbonos alfas del tubo 1
    '''
    traj = md.load(run_name+"_MD.nc", top=run_name+".parm7")
    
    CAs = traj.top.select("name==CA")
    CAs_tube1 = CAs[0:int(CAs.size/N_tubes)]
    CAs_tube2 = CAs[int(CAs.size/N_tubes):int(2*CAs.size/N_tubes)]
    CAs_tube3 = CAs[int(2*CAs.size/N_tubes):int(3*CAs.size/N_tubes)]
    CAs_tube4 = CAs[int(3*CAs.size/N_tubes):int(4*CAs.size/N_tubes)]
    CAs_top = np.concatenate((CAs_tube1[0:N_res], CAs_tube2[0:N_res], CAs_tube3[0:N_res], CAs_tube4[0:N_res]))
    CAs_bot = np.concatenate((CAs_tube1[-N_res:], CAs_tube2[-N_res:], CAs_tube3[-N_res:], CAs_tube4[-N_res:]))
    
    step = len(traj)-1
    p0 = np.sum(traj.xyz[step][CAs], axis=0)/CAs.size
    pZ = np.sum(traj.xyz[step][CAs_top], axis=0)/CAs_top.size
    pX = np.sum(traj.xyz[step][CAs_tube1], axis=0)/CAs_tube1.size
    
    selection = CAs
    xyz = []
    for atom in selection:
        xyz.append(orient(traj.xyz[step][atom], p0, pZ, pX))

    oriented = md.Trajectory(xyz, traj.top.subset(selection))
    
    traj.superpose(oriented, 0, atom_indices=selection, ref_atom_indices=range(oriented.n_atoms))
    traj.save(run_name+"_RMSD.nc")
#end


class MyAtom():
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    
    def __init__(self, top, N_rings, N_res, index):
        atom = top.atom(index)
        self.index = index
        self.name = atom.name
        self.resid = atom.residue.index
        self.resname = atom.residue.name
        if atom.residue.is_protein:
            self.tube = atom.residue.index//(N_rings*N_res)
            self.layer = (atom.residue.index//N_res)%N_rings
        else:
            self.tube = None
            self.layer = None
            
    @classmethod
    def from_string(cls, string):
        regex = r"MyAtom\(index=(\d+), name=(\w+), resid=(\d+), resname=(\w+), tube=(\w+), layer=(\w+)\)"
        result = re.search(regex, string)
        myatom = cls.__new__(cls)
        myatom.index = int(result.groups()[0])
        myatom.name = result.groups()[1]
        myatom.resid = int(result.groups()[2])
        myatom.resname = result.groups()[3]
        if result.groups()[4] == "None":
            myatom.tube = None
            myatom.layer = None
        else:
            myatom.tube = int(result.groups()[4])
            myatom.layer = int(result.groups()[5])
        return myatom
    
    def __str__(self): # print()
        return "MyAtom(index=" + str(self.index) + ", name=" + str(self.name) + ", resid=" + str(self.resid) + ", resname=" + str(self.resname) + ", tube=" + str(self.tube) + ", layer=" + str(self.layer) + ")"
    
    def __repr__(self):
        return "MyAtom(index=" + str(self.index) + ", name=" + str(self.name) + ", resid=" + str(self.resid) + ", resname=" + str(self.resname) + ", tube=" + str(self.tube) + ", layer=" + str(self.layer) + ")"
#end


def select_atoms(top, N_rings, N_res, selection):
    return np.array([
        MyAtom(top, N_rings, N_res, index)
        for index in top.select(selection)
    ])
#end


class MyParams:
    def __init__(self, traj, N_tubes, N_rings, N_res, selections):
        
        self.N_tubes = N_tubes # Número de tubos en el sistema
        self.N_rings = N_rings # Número de anillos en un tubo
        self.N_res = N_res # Número de residuos en un anillo
        N_allres = N_tubes*N_rings*N_res
        self.N_allres = N_allres
        
        top = traj.top
        self.CAs = select_atoms(top, N_rings, N_res, "name CA")
        self.bbNs = select_atoms(top, N_rings, N_res, "name N and resid 0 to " + str(N_allres-1))
        self.bbOs = select_atoms(top, N_rings, N_res, "name O and resid 0 to " + str(N_allres-1))
        
        bondable = np.array([], dtype=int)
        for selection in selections:
            bondable = np.concatenate((bondable, select_atoms(top, N_rings, N_res, selection)))
        self.bondable = bondable
        
        self.WATs = traj.top.select("water and name O")
        self.IONs = traj.top.select("element Cl")
#end


def get_reslist(N_rings, N_res, tube, residues):
    '''
    Uso: get_reslist(6, 8, 0, [1, 3])
    Devuelve: [1, 11, 17, 27, 33, 43]
    Cada ring tiene 8 residuos, que numeramos del 0 al 7.
    Para el primer ring del tubo 0 devuelve el resid del residuo 1 del ring,
    para el segundo ring del tubo 0 devuelve el resid del residuo 3 del ring.
    Para el tercer ring de nuevo el residuo 1, para el cuarto ring el residuo 3 y así...
    '''
    reslist = []
    for layer in range(N_rings):
        reslist.append(tube*N_rings*N_res + layer*N_res + residues[layer%len(residues)])
    return reslist
#end


def get_channel_reslist(N_rings, N_res, tubes, tuberesidues):
    '''
    Uso: get_channel_reslist(6, 8, [0, 1], [[1, 3], [5, 5])
    Devuelve concatenadas get_reslist(6, 8, 0, [1, 3]) y get_reslist(6, 8, 1, [5, 5])
    '''
    reslist = []
    for i, tube in enumerate(tubes):
        reslist += get_reslist(N_rings, N_res, tube, tuberesidues[i])
    return reslist
#end


def get_atoms_in_reslist(myatoms, reslist):
    return np.array([atom for atom in myatoms if atom.resid in reslist])
#end


def get_indices(myatoms):
    return np.array([atom.index for atom in myatoms])
#end


def get_indices_in_layer(myatoms, layer):
    return np.array([atom.index for atom in myatoms if atom.layer == layer])
#end


def get_indices_between_layers(myatoms, firstlayer, lastlayer):
    indices = np.array([], dtype=int)
    for layer in range(firstlayer, lastlayer+1):
        indices = np.concatenate((indices, get_indices_in_layer(myatoms, layer)))
    return indices
#end


def wrap_coordinates(point, box):
    # De momento, solo para celdas ortoédricas...
    wrapped_point = np.copy(point)
    
    for dim in range(3):
        wrapped_point[dim] = ((point[dim] + box[dim]/2) % box[dim]) - box[dim]/2
    
    return wrapped_point
#end


def periodic_pdist(positions, box):
    dist = 0
    
    for dim in range(3):
        pddim = pdist(positions[:, dim].reshape(positions.shape[0], 1))
        pddim[pddim > 0.5*box[dim]] -= box[dim] # ? apply boundary conditions
        dist += pddim**2
        
    return np.sqrt(dist)
#end


def decorate_ax(ax, title, titlesize, xlabel, ylabel, labelsize, ticksize, linewidth, length, legend):
    ax.set_title(title, fontsize=titlesize, pad=titlesize)
    ax.set_xlabel(xlabel, fontsize=labelsize)#, labelpad=15)
    ax.set_ylabel(ylabel, fontsize=labelsize)#, labelpad=15)
    ax.tick_params(top=True, right=True, labelsize=ticksize, width=linewidth, length=length)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    if legend:
        ax.legend(fontsize=ticksize, edgecolor='0.75')
#end


### EOF ###