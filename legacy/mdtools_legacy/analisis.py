### IMPORTS ###


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import axes3d
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import pytraj as pt
import mdtraj as md
import pandas as pd
import math
import time

from main import *


### ANÁLISIS ###


def get_indices(traj, WATs, CLs, CAs, N_res, layer=0, boundary=None,
                delta=1.0, deltaR=None, deltaZ=None, offset=None,
                preselected=False, save=True, savefileWATs="iWATs", savefileCLs="iCLs", first=None, last=None):
    
    # layer y boundary son el número de anillos que se excluyen de sus respectivas zonas:
    #               si es 0, se selecciona todo el nanotubo
    #               si es n, se selecciona el nanotubo menos los n primeros y los n últimos anillos
    
    # CAs son los carbonos alfa de interés:
    # tiene que ser CAs_tubej si es un tubo
    #               CAs_canal si es el canal
    #               CAs_topbot si es todo el bundle
    # N_res es el número de residuos que hay en una "capa" por lo que:
    # tiene que ser N_res si es un tubo
    #               N_tubes si es el canal
    #               N_res*N_tubes si es todo el bundle
    # layer y boundary tienen que ser 0 si es todo el bundle
    
    if boundary is None:
        boundary = layer
    
    if deltaR is None:
        deltaR = delta
    if deltaZ is None:
        deltaZ = delta
    
    if offset is None:
        offset = np.array([0.0, 0.0, 0.0])
    
    if first is None:
        first = 0
    if last is None:
        last = len(traj)
    
    iWATs = []
    iCLs = []
    if layer != boundary:
        iWATs_b = []
        iCLs_b = []
    
    atoms_top = CAs[layer*N_res:(layer+1)*N_res]
    atoms_bot = CAs[-(layer+1)*N_res:len(CAs)-layer*N_res]
    if layer != boundary:
        atoms_top_b = CAs[boundary*N_res:(boundary+1)*N_res]
        atoms_bot_b = CAs[-(boundary+1)*N_res:len(CAs)-boundary*N_res]
    
    if preselected:
        auxWATs = WATs
        auxCLs = CLs
    
    for step in range(first, last):
        frame = traj[step]
        
        if preselected:
            WATs = auxWATs[step]
            CLs = auxCLs[step]
        
        # Centro de la región
        centertop = np.sum(frame[atoms_top], axis=0)/atoms_top.size
        centerbot = np.sum(frame[atoms_bot], axis=0)/atoms_bot.size
        center = (centertop + centerbot)/2
        # Radio de la región
        Rtop = np.max(distance_matrix(frame[atoms_top], frame[atoms_top]))
        Rbot = np.max(distance_matrix(frame[atoms_bot], frame[atoms_bot]))
        R = max(Rtop, Rbot)/2 + deltaR
        # Alturas máxima y mínima de la región
        Zmax = np.sum(frame[atoms_top][:,2])/atoms_top.size - center[2] + deltaZ
        Zmin = np.sum(frame[atoms_bot][:,2])/atoms_bot.size - center[2] - deltaZ
        if layer != boundary:
            Zmax_b = np.sum(frame[atoms_top_b][:,2])/atoms_top_b.size - center[2] + deltaZ
            Zmin_b = np.sum(frame[atoms_bot_b][:,2])/atoms_bot_b.size - center[2] - deltaZ
        
        # Aguas (solo los oxígenos) en la región
        aux = []
        for atom in WATs:
            xyz = frame[atom] - (center + offset)
            if (Zmin < xyz[2]) and (xyz[2] < Zmax) and (xyz[0]**2+xyz[1]**2 < R**2):
                aux.append(atom)
        aux = np.array(aux)
        iWATs.append(aux)
        
        if layer != boundary:
            aux = []
            for atom in WATs:
                xyz = frame[atom] - (center + offset)
                if (Zmax < xyz[2]) and (xyz[2] < Zmax_b) and (xyz[0]**2+xyz[1]**2 < R**2):
                    aux.append(atom)
                elif (Zmin_b < xyz[2]) and (xyz[2] < Zmin) and (xyz[0]**2+xyz[1]**2 < R**2):
                    aux.append(atom)
            aux = np.array(aux)
            iWATs_b.append(aux)
        
        # Cloros en la región
        aux = []
        for atom in CLs:
            xyz = frame[atom] - (center + offset)
            if (Zmin < xyz[2]) and (xyz[2] < Zmax) and (xyz[0]**2+xyz[1]**2 < R**2):
                aux.append(atom)
        aux = np.array(aux)
        iCLs.append(aux)
        
        if layer != boundary:
            aux = []
            for atom in CLs:
                xyz = frame[atom] - (center + offset)
                if (Zmax < xyz[2]) and (xyz[2] < Zmax_b) and (xyz[0]**2+xyz[1]**2 < R**2):
                    aux.append(atom)
                elif (Zmin_b < xyz[2]) and (xyz[2] < Zmin) and (xyz[0]**2+xyz[1]**2 < R**2):
                    aux.append(atom)
            aux = np.array(aux)
            iCLs_b.append(aux)
    
    iWATs = np.array(iWATs, dtype=object)
    iCLs = np.array(iCLs, dtype=object)
    if layer != boundary:
        iWATs_b = np.array(iWATs_b, dtype=object)
        iCLs_b = np.array(iCLs_b, dtype=object)
    
    if save:
        np.save(savefileWATs+".npy", iWATs)
        np.save(savefileCLs+".npy", iCLs)
        if layer != boundary:
            np.save(savefileWATs+"_b.npy", iWATs_b)
            np.save(savefileCLs+"_b.npy", iCLs_b) 
#end


def get_indices_xtal(traj, WATs, CLs, CAs, N_res, deltaR=0.0, offsets=None, # En unidades de los lattice vectors
                     save=True, savefileWATs="iWATs", savefileCLs="iCLs", first=None, last=None):

    # CAs son los carbonos alfa de interés:
    # tiene que ser CAs_tubej si es un tubo
    #               CAs_canal si es el canal
    #               CAs_topbot si es todo el bundle
    # N_res es el número de residuos que hay en una "capa" por lo que:
    # tiene que ser N_res si es un tubo
    #               N_tubes si es el canal
    #               N_res*N_tubes si es todo el bundle
    
    if offsets is None:
        offsets = [np.array([0.0, 0.0, 0.0])]
    
    if first is None:
        first = 0
    if last is None:
        last = len(traj)
    
    iWATs = []
    iCLs = []
    
    layer = 0
    atoms_top = CAs[layer*N_res:(layer+1)*N_res]
    atoms_bot = CAs[-(layer+1)*N_res:len(CAs)-layer*N_res]
    
    for step in range(first, last):
        frame = traj[step]
        lvs = np.array([frame.box[0], frame.box[1], frame.box[2]])
        
        # Centro de la región
        centertop = np.sum(frame[atoms_top], axis=0)/atoms_top.size
        centerbot = np.sum(frame[atoms_bot], axis=0)/atoms_bot.size
        center = (centertop + centerbot)/2
        # Radio de la región
        Rtop = np.max(distance_matrix(frame[atoms_top], frame[atoms_top]))
        Rbot = np.max(distance_matrix(frame[atoms_bot], frame[atoms_bot]))
        R = max(Rtop, Rbot)/2 + deltaR
        
        # Aguas (solo los oxígenos) en la región
        aux = []
        for atom in WATs:
            for offset in offsets:
                xyz = frame[atom] - (center + offset*lvs)
                if (xyz[0]**2+xyz[1]**2 < R**2):
                    aux.append(atom)
        aux = np.array(aux)
        iWATs.append(aux)
        
        # Cloros en la región
        aux = []
        for atom in CLs:
            for offset in offsets:
                xyz = frame[atom] - (center + offset*lvs)
                if (xyz[0]**2+xyz[1]**2 < R**2):
                    aux.append(atom)
        aux = np.array(aux)
        iCLs.append(aux)
    
    iWATs = np.array(iWATs, dtype=object)
    iCLs = np.array(iCLs, dtype=object)
    
    if save:
        np.save(savefileWATs+".npy", iWATs)
        np.save(savefileCLs+".npy", iCLs)
#end


def compute_distance(frame, atom_a, atom_b, lvs):
    
    a = frame[atom_a]
    b = frame[atom_b]
    ba = a - b
    
    if lvs is not None:
#         lvs = np.array([frame.box[0], frame.box[1], frame.box[2]])
        for dim in range(3):
            if np.abs(ba[dim]) > lvs[dim]/2:
                ba[dim] = lvs[dim] - np.abs(ba[dim])
    
    distance = np.linalg.norm(ba)
    
    return distance
#end


def compute_angle(frame, atom_a, atom_b, atom_c, lvs):
    
    a = frame[atom_a]
    b = frame[atom_b]
    c = frame[atom_c]
    
    if lvs is not None:
#         lvs = np.array([frame.box[0], frame.box[1], frame.box[2]])
        a = wrap_coordinates(a-b, lvs)
        c = wrap_coordinates(c-b, lvs)
        b = wrap_coordinates(b-b, lvs)
    
    ba = a - b
    bc = c - b

    cosine = np.dot(ba, bc)/(np.linalg.norm(ba)*np.linalg.norm(bc))
    angle = np.arccos(cosine)
    
    return angle
#end


def check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, ihs, ilabel, jatom, jhs, jlabel, lvs,
                 hbondsdf, N_hbonds, d_ave):
    
    if ihs != 0 and jhs != 0:
        dim = max(ihs, jhs) + 1
        ds = 999.9*np.ones((dim, dim))
        
        for h in range(1, ihs+1):
            ds[h, 0] = compute_distance(frame, iatom+h, jatom, lvs)
            # ds[h, 0] = np.linalg.norm(frame[iatom+h] - frame[jatom]) # del h de iatom a jatom
        for h in range(1, jhs+1):
            ds[0, h] = compute_distance(frame, iatom, jatom+h, lvs)
            # ds[0, h] = np.linalg.norm(frame[iatom] - frame[jatom+h]) # de iatom al h de jatom
        for index, d in enumerate(ds.flatten()):
            if d < distance_cutoff:
                bond = np.unravel_index(index, ds.shape)
                if bond[0] > bond[1]: # iatom es el donor
                    if compute_angle(frame, iatom, iatom+bond[0], jatom, lvs) > angle_cutoff:
                        N_hbonds += 1
                        d_ave += d
                        hbondsdf.loc[hbondsdf.shape[0]] = [step, iatom, iatom+bond[0], jatom, ilabel+"-"+jlabel, d]
                else: # jatom es el donor
                    if compute_angle(frame, jatom, jatom+bond[1], iatom, lvs) > angle_cutoff:
                        N_hbonds += 1
                        d_ave += d
                        hbondsdf.loc[hbondsdf.shape[0]] = [step, jatom, jatom+bond[1], iatom, jlabel+"-"+ilabel, d]
                            
    else:
        dim = max(ihs, jhs) + 1
        ds = 999.9*np.ones(dim)
        
        if jhs == 0:
            for h in range(1, ihs+1):
                ds[h] = compute_distance(frame, iatom+h, jatom, lvs)
                # ds[h] = np.linalg.norm(frame[iatom+h] - frame[jatom]) # del h de iatom a jatom
            for bond, d in enumerate(ds):
                if d < distance_cutoff:
                    if compute_angle(frame, iatom, iatom+bond, jatom, lvs) > angle_cutoff:
                        N_hbonds += 1
                        d_ave += d
                        hbondsdf.loc[hbondsdf.shape[0]] = [step, iatom, iatom+bond, jatom, ilabel+"-"+jlabel, d]
        else:
            for h in range(1, jhs+1):
                ds[h] = compute_distance(frame, iatom, jatom+h, lvs)
                # ds[h] = np.linalg.norm(frame[iatom] - frame[jatom+h]) # de iatom al h de jatom
            for bond, d in enumerate(ds):
                if d < distance_cutoff:
                    if compute_angle(frame, jatom, jatom+bond, iatom, lvs) > angle_cutoff:
                        N_hbonds += 1
                        d_ave += d
                        hbondsdf.loc[hbondsdf.shape[0]] = [step, jatom, jatom+bond, iatom, jlabel+"-"+ilabel, d]
    
    return N_hbonds, d_ave
#end


def analyse(p, traj, label, distance_cutoff=2.5, angle_cutoff=2*np.pi/3, canal=False, layer=0, boundary=None,
            KY=False, first=None, last=None, bblabel=None, xtal=False):
    
    if boundary is None:
        boundary = layer
        
    if bblabel is None:
        bblabel = label
    
    if not xtal:
        lvs = None
    
    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        if layer != boundary:
            NZs_LYS_b, NZs_LYN_b = select_atoms(p, boundary) # incluye NZs_LYS y NZs_LYN
            NZs_LYS_b = np.setdiff1d(NZs_LYS_b, NZs_LYS)
            NZs_LYN_b = np.setdiff1d(NZs_LYN_b, NZs_LYN)
        else:
            NZs_LYS_b = np.array([], dtype=int)
            NZs_LYN_b = np.array([], dtype=int)
            
        if KY:
            TYDs = select_atoms(p, layer, KY)
            if layer != boundary:
                TYDs_b = select_atoms(p, boundary, KY)
                TYDs_b = np.setdiff1d(TYDs_b, TYDs)
            else:
                TYDs_b = np.array([], dtype=int)
    
    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)
    if layer != boundary:
        iWATs_b = np.load("iWATs_"+label+"_b.npy", allow_pickle=True)
        iCLs_b = np.load("iCLs_"+label+"_b.npy", allow_pickle=True)
    
    if first is None:
        first = 0
    if last is None:
        last = len(traj)
    
    # distance_cutoff = 2.5 # Solo contamos los puentes con una distancia menor a una dada
    # angle_cutoff = 2*np.pi/3
    
    # Base de datos de las estadísticas de las aguas
    # Step | Número de aguas | Número de cloros | Número de puentes| Distancia media puente |
    statdf = pd.DataFrame(columns=['istep', 'N_wats', 'N_cls', 'N_hbonds', 'ave_dist'])
    
    # Base de datos de los puentes de H
    # Step | Índice del "donor" | Índice del H (del donor) | Índice del "acceptor" | Residuos | Distancia
    hbondsdf = pd.DataFrame(columns=['istep', 'donor', 'H', 'acceptor', 'residues', 'dist'])
    
    for step in range(first, last):
        frame = traj[step]
        
        WATs = iWATs[step]
        CLs = iCLs[step]
        if layer != boundary:
            WATs_b = iWATs_b[step]
            CLs_b = iCLs_b[step]
        else:
            WATs_b = np.array([], dtype=int)
            CLs_b = np.array([], dtype=int)
        
        if xtal:
            lvs = np.array([frame.box[0], frame.box[1], frame.box[2]])
        
        N_hbonds = 0
        d_ave = 0.0
        
        ### Bucle sobre las aguas ###
        for i in range(0, len(WATs)):
            iatom = WATs[i]
            
            # Puentes de H agua-agua
            for jatom in np.concatenate((WATs[i+1:], WATs_b)).astype(int):
#             for j in range(i+1, len(WATs)):
#                 jatom = WATs[j]
                N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "WAT", jatom, 2, "WAT", lvs, hbondsdf, N_hbonds, d_ave)
            
            # Puentes de H agua(H)-Cl
            for jatom in np.concatenate((CLs, CLs_b)).astype(int):
                N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "WAT", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
            
            # Puentes de H agua-backbone ### Creo que se puede optimizar... ###
            for jatom in p.bbNs[bblabel]:
                N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 0, "WAT", jatom, 1, "N", lvs, hbondsdf, N_hbonds, d_ave)
            
            for jatom in p.bbOs[bblabel]:
                N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "WAT", jatom, 0, "O", lvs, hbondsdf, N_hbonds, d_ave)
        
        if canal:
            
            if KY:
                ### Bucle sobre las TYD ###
                for i in range(0, len(TYDs)):
                    iatom = TYDs[i]
                    
                    # Puentes de H TYD-agua
                    for jatom in np.concatenate((WATs, WATs_b)).astype(int):
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 2, "WAT", lvs, hbondsdf, N_hbonds, d_ave)

                    # Puentes de H TYD-TYD
                    for jatom in np.concatenate((TYDs[i+1:], TYDs_b)).astype(int):
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 1, "TYD", lvs, hbondsdf, N_hbonds, d_ave)

                    # Puentes de H TYD-LYN
                    for jatom in np.concatenate((NZs_LYN, NZs_LYN_b)).astype(int):
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 2, "LYN", lvs, hbondsdf, N_hbonds, d_ave)
                    
                    # Puentes de H TYD-LYS
                    for jatom in np.concatenate((NZs_LYS, NZs_LYS_b)).astype(int):
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 0, "TYD", jatom, 3, "LYS", lvs, hbondsdf, N_hbonds, d_ave)

                    # Puentes de H TYD-CL
                    for jatom in np.concatenate((CLs, CLs_b)).astype(int):
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
                            
            ### Bucle sobre las LYN ###
            for i in range(0, len(NZs_LYN)):
                iatom = NZs_LYN[i]
                
                # Puentes de H LYN-agua
                for jatom in np.concatenate((WATs, WATs_b)).astype(int):
                    N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "LYN", jatom, 2, "WAT", lvs, hbondsdf, N_hbonds, d_ave)
                
                # Puentes de H LYN-LYN
                # for j in range(i+1, len(NZs_LYN)):
                for jatom in np.concatenate((NZs_LYN[i+1:], NZs_LYN_b)).astype(int):
                    # jatom = NZs_LYN[j]
                    N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "LYN", jatom, 2, "LYN", lvs, hbondsdf, N_hbonds, d_ave)
                
                # Puentes de H LYN(H)-Cl
                for jatom in np.concatenate((CLs, CLs_b)).astype(int):
                    N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "LYN", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
                
            ### Bucle sobre las LYS ###
            for iatom in NZs_LYS:
                
                # Puentes de H LYS(H)-agua(O)
                for jatom in np.concatenate((WATs, WATs_b)).astype(int):
                    N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 3, "LYS", jatom, 0, "WAT", lvs, hbondsdf, N_hbonds, d_ave)
                
                # Puentes de H LYS(H)-LYN(N)
                for jatom in np.concatenate((NZs_LYN, NZs_LYN_b)).astype(int):
                    N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 3, "LYS", jatom, 0, "LYN", lvs, hbondsdf, N_hbonds, d_ave)
                
                # Puentes de H LYS(H)-Cl
                for jatom in np.concatenate((CLs, CLs_b)).astype(int):
                    N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 3, "LYS", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
                                
            # Puentes que nos pueden faltar en el caso con boundary
            if layer != boundary:
                for iatom in WATs_b:
                    for jatom in CLs:
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "WAT", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
                
                for iatom in NZs_LYS_b:
                    for jatom in WATs:
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 3, "LYS", jatom, 0, "WAT", lvs, hbondsdf, N_hbonds, d_ave)

                    for jatom in CLs:
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 3, "LYS", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)

                    for jatom in NZs_LYN:
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 3, "LYS", jatom, 0, "LYN", lvs, hbondsdf, N_hbonds, d_ave)
                        
                for iatom in NZs_LYN_b:
                    for jatom in WATs:
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "LYN", jatom, 2, "WAT", lvs, hbondsdf, N_hbonds, d_ave)
                        
                    for jatom in CLs:
                        N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 2, "LYN", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
                        
                if KY:
                    for iatom in TYDs_b:
                        for jatom in WATs:
                            N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 2, "WAT", lvs, hbondsdf, N_hbonds, d_ave)
                        for jatom in CLs:
                            N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 0, "CL", lvs, hbondsdf, N_hbonds, d_ave)
                        for jatom in NZs_LYN:
                            N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 1, "TYD", jatom, 2, "LYN", lvs, hbondsdf, N_hbonds, d_ave)
                        for jatom in NZs_LYS:
                            N_hbonds, d_ave = check_hbonds(distance_cutoff, angle_cutoff, step, frame, iatom, 0, "TYD", jatom, 3, "LYS", lvs, hbondsdf, N_hbonds, d_ave)
            
        if N_hbonds != 0:
            d_ave = d_ave/N_hbonds
        statdf.loc[statdf.shape[0]] = [step, len(WATs), len(CLs), N_hbonds, d_ave]
    
    hbondsdf.to_csv(label+"_hbonds.csv")
    statdf.to_csv(label+"_stats.csv")
#end


def detail_hbonds(label):
    hbonds = pd.read_csv(label+"_hbonds.csv")
    Nsteps = hbonds['istep'].max()+1
    # Base de datos de los puentes de H del canal
    # Número de frame | Residuos | Número de puentes | Distancia media
    detail = pd.DataFrame(columns=['istep', 'residues', 'N_hbonds', 'ave_dist'])
    
    for step in range(Nsteps):
        aux = hbonds[hbonds["istep"] == step]
        for pair in hbonds["residues"].unique():
            N = aux[aux["residues"] == pair].shape[0]
            ave_dist = aux[aux["residues"] == pair]["dist"].mean()
            detail.loc[detail.shape[0]] = [step, pair, N, ave_dist]

    detail.to_csv(label+"_hbonds_detail.csv")
#end


def stability_waters(label, boundary=False):
    
    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    if boundary:
        iWATs_b = np.load("iWATs_"+label+"_b.npy", allow_pickle=True)
    aux = np.array([], dtype=int)
    for WATs in iWATs:
        aux = np.append(aux, WATs)
    inWATs = np.unique(aux)
    
    stability = pd.DataFrame(columns=['atom', 'steps', 'consecutives'])
    
    aux = [[] for water in inWATs]
    aux_b = [[] for water in inWATs]
    for step in range(len(iWATs)):
        for i, water in enumerate(inWATs):
            if water in iWATs[step]:
                aux[i].append(step)
            elif boundary and (water in iWATs_b[step]):
                aux_b[i].append(step)
            
            if step == len(iWATs)-1:
                stability.loc[stability.shape[0]] = [water, aux[i], '']
    steps_b = aux_b
    
    for index, water in stability.iterrows():
        steps = water["steps"]
        consecutives = 1
        aux = []

        for i in range(0, len(steps)-1):
            if steps[i]+1 == steps[i+1]:
                consecutives += 1
            else:
                aux.append(consecutives)
                consecutives = 1
        
        aux.append(consecutives)
        stability.at[index, "consecutives"] = aux

    stability.to_csv(label+"_waters_stability.csv")
#end


def stability_hbonds(label):
    hbonds = pd.read_csv(label+"_hbonds.csv")
    stability = pd.DataFrame(columns=['donor', 'H', 'acceptor', 'residues', 'steps', 'consecutives'])
    
    for donor in hbonds["donor"].unique():
        aux = hbonds[hbonds["donor"] == donor]

        for acceptor in aux["acceptor"].unique():
            aux2 = aux[aux["acceptor"] == acceptor]
            residues = aux2["residues"].iloc[0]

            for H in aux2["H"].unique():
                aux3 = aux2[aux2["H"] == H]
                steps = aux3["istep"].unique().tolist()
                stability.loc[stability.shape[0]] = [donor, H, acceptor, residues, steps, '']
    
    for index, bond in stability.iterrows():
        steps = bond["steps"]
        consecutives = 1
        aux = []

        for i in range(0, len(steps)-1):
            if steps[i]+1 == steps[i+1]:
                consecutives += 1
            else:
                aux.append(consecutives)
                consecutives = 1
        
        aux.append(consecutives)
        stability.at[index, "consecutives"] = aux

    stability.to_csv(label+"_hbonds_stability.csv")
#end


def search_path(frame, graph_df, start_node, direction, xtal, temp_path=[]):
    '''    
    # ARGUMENTOS #
    - graph_df es un dataframe con los nodos y conexiones (el dataframe de hbonds)
    
    - start_node es el nodo inicial por el que se empieza a buscar un camino
    
    - direction es el sentido en el que estamos buscando caminos (hacia arriba o hacia abajo)
    
    - xtal es si tenemos que tener en cuenta PBCs (en concreto en el eje z)
    
    - temp_path es el path que llevamos construido hasta el momento
    
    # SALIDA #
    - path es la lista de nodos que empieza en el nodo inicial start_node y termina en el nodo "más alejado" de
      start_node al que puede se puede llegar usando las conexiones entre nodos
    
    - path_length es la longitud de path, es decir, la suma de las distancias de las conexiones entre nodos que
      hay que usar para llegar del nodo inicial al final
      
    - dz es la diferencia de las coordenadas z del nodo inicial y final
    
    En caso de que haya dos path que lleguen al mismo nodo final, se toma el path con menor path_length
     (en nuestro caso un menor path_length significa menores distancias de puentes de H / menos puentes de H)
    '''
    
    path = [start_node]
    residues = []
    path_length = 0
    dz = 0
    
    if xtal and (-start_node in temp_path):
        return path, residues, path_length, dz
    aux_temp_path = temp_path + [start_node]

    # Dataframe auxiliar solo con los puentes que forma el start_node como acceptor
    aux_df = graph_df[graph_df["acceptor"]==abs(start_node)]
    if aux_df.empty:
        # En el caso en que el nodo no sea un acceptor y por tanto no tenga nodos conectados,
        # devolvemos los valores iniciales
        return path, residues, path_length, dz
    
    # Iteramos sobre las filas
    for i, row in aux_df.iterrows():
        donor = row["donor"]

        # Añadimos la distancia de start_node a neighbour_node a aux_path_length,
        # pero de momento no agregamos neighbour_node a aux_path, ya que esto se hace en el paso siguiente
        aux_path = [start_node]
        aux_residues = [row["residues"]]
        aux_path_length = row["dist"]
        aux_direction = np.sign(frame[abs(start_node)][2] - frame[donor][2])
        
        wrap = np.sign(start_node)
        if xtal:
            lvz = frame.box[2]
            if np.abs(frame[abs(start_node)][2] - frame[donor][2]) > lvz/2:
                if wrap < 0:
                    return path, residues, path_length, dz
                wrap = -1
                aux_direction *= -1
        
        if direction != 0:
            if direction != aux_direction:
                continue

        # Aplicamos la función al donor: result[0] es el camino que devuelve y comienza por sí mismo,
        # por eso no lo agregamos a aux_path antes
        result = search_path(frame, graph_df, wrap*donor, aux_direction, xtal, aux_temp_path)

        aux_path += result[0] # ¡No es una suma! Concatenamos las listas de nodos de los caminos
        aux_residues += result[1] # ¡No es una suma! Concatenamos las listas de los tipos de hbonds de los caminos
        aux_path_length += result[2] # Sumamos las distancias de los caminos
        if aux_path[-1] < 0:
            aux_dz = frame[start_node][2] - (frame[abs(aux_path[-1])][2] - aux_direction*lvz)
        else:
            aux_dz = frame[start_node][2] - frame[aux_path[-1]][2]
        
        if (abs(aux_dz) > abs(dz)) or (aux_path[-1] == path[-1] and aux_path_length < path_length):
            # Nos quedamos con aux_path si llega a un nodo "más alejado" que path o, en caso de llegar al mismo nodo,
            # si aux_path_length es menor que path_length
            path = aux_path
            residues = aux_residues
            path_length = aux_path_length
            dz = aux_dz
    
    return path, residues, path_length, dz
#end



def save_paths(traj, label, xtal=False):
    hbonds = pd.read_csv(label+"_hbonds.csv") # base de datos con los puentes de H
    paths = pd.DataFrame(columns=['istep', 'path', 'residues', 'path_length', 'dz'])

    for step in range(len(traj)):
        frame = traj[step]
        aux = hbonds[hbonds["istep"] == step]

        path = []
        residues = []
        path_length = 0
        dz = 0
        
        for acceptor in aux["acceptor"].unique():
            result = search_path(frame, aux, acceptor, 0, xtal)

            if abs(result[3]) > abs(dz):
                path = result[0]
                residues = result[1]
                path_length = result[2]
                dz = result[3]

        paths.loc[paths.shape[0]] = [step, path, residues, path_length, dz]

    paths.to_csv(label+"_paths.csv")
#end


def search_path_res(step, graph_df, res_wanted, iteration, path_df, path=[], residues=[], path_length=0):
    '''    
    # ARGUMENTOS #
    - graph_df es un dataframe con los nodos y conexiones (el dataframe de hbonds)
    
    - res_wanted son los residuos que queremos que formen el path
    
    - path_df es un dataframe en el que escribimos los paths encontrados
    
    # SALIDA #
    - path es la lista de nodos
    
    - path_length es la longitud de path, es decir, la suma de las distancias de las conexiones entre nodos que
      hay que usar para llegar del nodo inicial al final
      
    - dz es la diferencia de las coordenadas z del nodo inicial y final
    '''
    
    # Dataframe auxiliar solo con los puentes que forma el start_node como acceptor
    aux_df = graph_df[graph_df["residues"] == res_wanted[iteration]+"-"+res_wanted[iteration+1]]
    
    if aux_df.empty:
        # En el caso en el que no haya ningún hbond del tipo que buscamos
        return
    
    # Iteramos sobre las filas (cada enlace)
    for i, row in aux_df.iterrows():
        donor = row["donor"]
        hydrogen = row["H"]
        acceptor = row["acceptor"]
        
        # Añadimos la distancia de start_node a neighbour_node a aux_path_length,
        # pero de momento no agregamos neighbour_node a aux_path, ya que esto se hace en el paso siguiente
        aux_path = path + [donor, hydrogen]
        aux_residues = residues + [row["residues"]]
        aux_path_length = path_length + row["dist"]
        
        if iteration == len(res_wanted)-2:
            # Si es la última iteración, concatenamos el donor al path y devolvemos el resultado
            aux_path += [acceptor] # ¡No es una suma! Concatenamos las listas de nodos de los caminos
            path_df.loc[path_df.shape[0]] = [step, aux_path, aux_residues, aux_path_length]
        else:
            # Aplicamos la función a los enlaces del donor
            donor_df = graph_df[graph_df["donor"] == acceptor]
            search_path_res(step, donor_df, res_wanted, iteration+1, path_df, aux_path, aux_residues, aux_path_length)
    
    return
#end


def save_paths_res(traj, label, res):
    hbonds = pd.read_csv(label+"_hbonds.csv") # base de datos con los puentes de H
    paths = pd.DataFrame(columns=['istep', 'path', 'residues', 'path_length'])
    res_wanted = res.split("-")
    
    for step in range(len(traj)):
        aux = hbonds[hbonds["istep"] == step]
        search_path_res(step, aux, res_wanted, 0, paths)
    
    paths.to_csv(label+"_paths_"+res+".csv")
#end


def hbonds_correlation(label, bondtype=False):
    hbonds = pd.read_csv(label+"_hbonds.csv")
    Nsteps = hbonds['istep'].max()+1

    if bondtype:
        hbonds = hbonds[hbonds['residues']==bondtype]
    ks = hbonds['H'].unique()
    NHs = len(ks)
    Hk0 = np.zeros((NHs))

    step = 0
    aux = hbonds[hbonds['istep'] == step]

    for index, H in enumerate(ks):
        if H in aux['H'].unique():
            Hk0[index] = aux[aux['H'] == H]['acceptor'].iloc[0]

    eta = []

    for step in range(0, Nsteps):
        aux = hbonds[hbonds['istep'] == step]

        delta = 0
        for index, H in enumerate(ks):
            if H in aux['H'].unique():
                acceptor = aux[aux['H'] == H]['acceptor'].iloc[0]
                if Hk0[index] == acceptor:
                    delta += 1

        eta.append(delta)

    eta = np.asarray(eta)#/NHs

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(range(Nsteps), eta)
    return fig, ax
#end


def hbonds_orientation(traj, label, Nbins=100, bondtype=False):
    hbonds = pd.read_csv(label+"_hbonds.csv")
    
    if bondtype:
        hbonds = hbonds[hbonds['residues'] == bondtype]
    
    thetas = []
    for step in range(len(traj)):
        frame = traj[step]
        aux = hbonds[hbonds['istep'] == step]
        
        for index, bond in aux.iterrows():
            acceptor = bond['acceptor']
            donor = bond['donor']
            vector = frame[acceptor] - frame[donor] # vector del donor al acceptor
            theta = np.arccos(np.dot(vector, [0, 0, 1])/np.linalg.norm(vector)) * 180/np.pi
            thetas.append(theta)
    
    thetas = np.asarray(thetas)
    hist, bin_edges = np.histogram(thetas, bins=Nbins)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, where='mid')
    return fig, ax
#end


def average_positions(p, traj, label, canal=False, layer=0, KY=False):
    
    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        if KY:
            TYDs = select_atoms(p, layer, KY)
    
    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)
    average_pos = pd.DataFrame(columns=['atom', 'name', 'count', 'ave_pos', 'var_pos'])
    
    for step in range(len(traj)):
        frame = traj[step]
        
        if canal:
            for atom in NZs_LYS:
                xyz = frame[atom]
                if step == 0:
                    average_pos.loc[average_pos.shape[0]] = np.asarray([atom, "NZ_LYS", 1, xyz.tolist(), np.array([0.0,0.0,0.0]).tolist()], dtype=object)
                else:
                    index = average_pos.index[average_pos['atom'] == atom].tolist()
                    count = average_pos[average_pos['atom'] == atom]['count'].iloc[0] + 1
                    ave = average_pos[average_pos['atom'] == atom]['ave_pos'].iloc[0]
                    delta = xyz - ave
                    ave += delta/count
                    M2 = average_pos[average_pos['atom'] == atom]['var_pos'].iloc[0]
                    delta2 = xyz - ave
                    M2 += delta*delta2
                    if step == len(traj)-1:
                        average_pos.loc[index] = [atom, "NZ_LYS", count, ave.tolist(), (M2/count).tolist()]
                    else:
                        average_pos.loc[index] = np.asarray([atom, "NZ_LYS", count, ave.tolist(), M2.tolist()], dtype=object)

            for atom in NZs_LYN:
                xyz = frame[atom]
                if step == 0:
                    average_pos.loc[average_pos.shape[0]] = np.asarray([atom, "NZ_LYN", 1, xyz.tolist(), np.array([0.0,0.0,0.0]).tolist()], dtype=object)
                else:
                    index = average_pos.index[average_pos['atom'] == atom].tolist()
                    count = average_pos[average_pos['atom'] == atom]['count'].iloc[0] + 1
                    ave = average_pos[average_pos['atom'] == atom]['ave_pos'].iloc[0]
                    delta = xyz - ave
                    ave += delta/count
                    M2 = average_pos[average_pos['atom'] == atom]['var_pos'].iloc[0]
                    delta2 = xyz - ave
                    M2 += delta*delta2
                    if step == len(traj)-1:
                        average_pos.iloc[index] = [atom, "NZ_LYN", count, ave.tolist(), (M2/count).tolist()]
                    else:
                        average_pos.iloc[index] = [atom, "NZ_LYN", count, ave.tolist(), M2.tolist()]
                        
            if KY:
                for atom in TYDs:
                    xyz = frame[atom]
                    if step == 0:
                        average_pos.loc[average_pos.shape[0]] = np.asarray([atom, "TYD", 1, xyz.tolist(), np.array([0.0,0.0,0.0]).tolist()], dtype=object)
                    else:
                        index = average_pos.index[average_pos['atom'] == atom].tolist()
                        count = average_pos[average_pos['atom'] == atom]['count'].iloc[0] + 1
                        ave = average_pos[average_pos['atom'] == atom]['ave_pos'].iloc[0]
                        delta = xyz - ave
                        ave += delta/count
                        M2 = average_pos[average_pos['atom'] == atom]['var_pos'].iloc[0]
                        delta2 = xyz - ave
                        M2 += delta*delta2
                        if step == len(traj)-1:
                            average_pos.iloc[index] = [atom, "TYD", count, ave.tolist(), (M2/count).tolist()]
                        else:
                            average_pos.iloc[index] = [atom, "TYD", count, ave.tolist(), M2.tolist()]

        WATs = iWATs[step]
        for atom in WATs:
            xyz = frame[atom]
            if atom not in average_pos['atom'].unique():
                average_pos.loc[average_pos.shape[0]] = np.asarray([atom, "WAT", 1, xyz.tolist(), np.array([0.0,0.0,0.0]).tolist()], dtype=object)
            else:
                index = average_pos.index[average_pos['atom'] == atom].tolist()
                count = average_pos[average_pos['atom'] == atom]['count'].iloc[0] + 1
                ave = average_pos[average_pos['atom'] == atom]['ave_pos'].iloc[0]
                delta = xyz - ave
                ave += delta/count
                M2 = average_pos[average_pos['atom'] == atom]['var_pos'].iloc[0]
                delta2 = xyz - ave
                M2 += delta*delta2
                if step == len(traj)-1:
                    average_pos.iloc[index] = [atom, "WAT", count, ave.tolist(), (M2/count).tolist()]
                else:
                    average_pos.iloc[index] = [atom, "WAT", count, ave.tolist(), M2.tolist()]

        CLs = iCLs[step]
        for atom in CLs:
            xyz = frame[atom]
            if atom not in average_pos['atom'].unique():
                average_pos.loc[average_pos.shape[0]] = np.asarray([atom, "CL", 1, xyz.tolist(), np.array([0.0,0.0,0.0]).tolist()], dtype=object)
            else:
                index = average_pos.index[average_pos['atom'] == atom].tolist()
                count = average_pos[average_pos['atom'] == atom]['count'].iloc[0] + 1
                ave = average_pos[average_pos['atom'] == atom]['ave_pos'].iloc[0]
                delta = xyz - ave
                ave += delta/count
                M2 = average_pos[average_pos['atom'] == atom]['var_pos'].iloc[0]
                delta2 = xyz - ave
                M2 += delta*delta2
                if step == len(traj)-1:
                    average_pos.iloc[index] = [atom, "CL", count, ave.tolist(), (M2/count).tolist()]
                else:
                    average_pos.iloc[index] = [atom, "CL", count, ave.tolist(), M2.tolist()]

    average_pos.to_csv(label+"_averages.csv")
#end


def RDF(p, traj, label, layer=0): # ¿meter boundary también?
    
    canal = True
    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
    
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)
    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)

    fig, axs = plt.subplots(4, 4, figsize=(12, 8), constrained_layout=True)

    # NZs_LYS, NZs_LYN, iWATs, iCLs
    selections = [NZs_LYS, NZs_LYN, iWATs, iCLs]
    iterselection = [False, False, True, True]

    for a in range(4):
        for b in range(4):

            # Inicilizar histograma
            lower = 0
            upper = 20
            dbin = 0.2
            Nbins = int(np.floor((upper-lower)/dbin))
            hist_r = np.zeros((Nbins))
            hist_z = np.zeros((Nbins))
            out_r = 0
            out_z = 0

            # Escoger selecciones de átomos
            selection1 = selections[a]
            iter1 = iterselection[a]
            selection2 = selections[b]
            iter2 = iterselection[b]
            sameselection = False
            if a == b:
                sameselection = True

            aux1 = selection1
            aux2 = selection2

            # Computar los histogramas y graficarlos
            for step in range(0, len(traj)):
                frame = traj[step]

                if iter1:
                    selection1 = aux1[step]
                if iter2:
                    selection2 = aux2[step]

                for i in range(0, len(selection1)):
                    first = 0
                    if sameselection:
                        first = i
                    for j in range(first+1, len(selection2)):
                        dist_r = np.linalg.norm(frame[selection1[i]] - frame[selection2[j]])
                        ibin = int(np.floor(Nbins*(dist_r - lower)/(upper - lower)))
                        if (ibin >= 0) & (ibin < Nbins):
                            hist_r[ibin] += 1
                        else:
                            out_r += 1
                        dist_z = np.abs(frame[selection1[i]][2] - frame[selection2[j]][2])
                        ibin = int(np.floor(Nbins*(dist_z - lower)/(upper - lower)))
                        if (ibin >= 0) & (ibin < Nbins):
                            hist_z[ibin] += 1
                        else:
                            out_z += 1

            axs[a, b].plot(np.linspace(lower+dbin/2, upper-dbin/2, Nbins), hist_r, label='r')
            axs[a, b].plot(np.linspace(lower+dbin/2, upper-dbin/2, Nbins), hist_z, label='z')
            axs[a, b].plot([3.7, 3.7], [0, 5000], 'k--')
            axs[a, b].legend(edgecolor='0.75')
            #print(out_r)
            #print(out_z)

    axs[0, 0].set_ylabel("LYS(NZ)")
    axs[1, 0].set_ylabel("LYN(NZ)")
    axs[2, 0].set_ylabel("WAT(O)")
    axs[3, 0].set_ylabel("Cl")
    axs[0, 0].set_title("LYS(NZ)")
    axs[0, 1].set_title("LYN(NZ)")
    axs[0, 2].set_title("WAT(O)")
    axs[0, 3].set_title("Cl")
    
    return fig, axs
#end


def coordination(p, traj, label, cutoff=3.6, layer=0, boundary=None):
    
    if boundary is None:
        boundary = layer
    
    canal = True
    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        if layer != boundary:
            NZs_LYS_b, NZs_LYN_b = select_atoms(p, boundary) # incluye NZs_LYS y NZs_LYN
            NZs_LYS_b = np.setdiff1d(NZs_LYS_b, NZs_LYS)
            NZs_LYN_b = np.setdiff1d(NZs_LYN_b, NZs_LYN)
        else:
            NZs_LYS_b = np.array([], dtype=int)
            NZs_LYN_b = np.array([], dtype=int)
    else:
        NZs_LYS = np.array([], dtype=int)
        NZs_LYN = np.array([], dtype=int)
        NZs_LYS_b = np.array([], dtype=int)
        NZs_LYN_b = np.array([], dtype=int)
    
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)
    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    if layer != boundary:
        iWATs_b = np.load("iWATs_"+label+"_b.npy", allow_pickle=True)
        iCLs_b = np.load("iCLs_"+label+"_b.npy", allow_pickle=True)

    # Base de datos
    # CN = Coordinaton Number
    # | Número de frame | Índice del átomo | Tipo de átomo | CN Cl | CN WAT(O) | CN LYS(NZ) | CN LYN(NZ) |
    nc_df = pd.DataFrame(columns=['step', 'atom', 'name', 'cn_cl', 'cn_wat', 'cn_lys', 'cn_lyn'])

    for step in range(0, len(traj)):
        frame = traj[step]
        CLs = iCLs[step]
        WATs = iWATs[step]
        if layer != boundary:
            CLs_b = iCLs_b[step]
            WATs_b = iWATs_b[step]
        else:
            CLs_b = np.array([], dtype=int)
            WATs_b = np.array([], dtype=int)

        for i, atom_i in enumerate(np.concatenate((CLs, WATs, NZs_LYS, NZs_LYN)).astype(int)):
            cn_cl = 0
            cn_wat = 0
            cn_lys = 0
            cn_lyn = 0
            if i < len(CLs):
                name = "CL"
            elif i < (len(CLs)+len(WATs)):
                name = "WAT"
            elif i < (len(CLs)+len(WATs)+len(NZs_LYS)):
                name = "LYS"
            else:
                name = "LYN"
            for j, atom_j in enumerate(np.concatenate((CLs, CLs_b, WATs, WATs_b, NZs_LYS, NZs_LYS_b, NZs_LYN, NZs_LYN_b)).astype(int)):
                if atom_i != atom_j:
                    distance = np.linalg.norm(frame[atom_i] - frame[atom_j])
                    if distance < cutoff:
                        if atom_j in np.concatenate((CLs, CLs_b)).astype(int):
                            cn_cl += 1
                        elif atom_j in np.concatenate((WATs, WATs_b)).astype(int):
                            cn_wat += 1
                        elif atom_j in np.concatenate((NZs_LYS, NZs_LYS_b)).astype(int):
                            cn_lys += 1
                        else:
                            cn_lyn += 1

            nc_df.loc[nc_df.shape[0]] = [step, atom_i, name, cn_cl, cn_wat, cn_lys, cn_lyn]

    nc_df.to_csv(label+"_coordination.csv")
#end


def coordination_hbonds(p, traj, label, layer=0, boundary=None):
    
    if boundary is None:
        boundary = layer
    
    canal = True
    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        if layer != boundary:
            NZs_LYS_b, NZs_LYN_b = select_atoms(p, boundary) # incluye NZs_LYS y NZs_LYN
            NZs_LYS_b = np.setdiff1d(NZs_LYS_b, NZs_LYS)
            NZs_LYN_b = np.setdiff1d(NZs_LYN_b, NZs_LYN)
        else:
            NZs_LYS_b = np.array([], dtype=int)
            NZs_LYN_b = np.array([], dtype=int)
    else:
        NZs_LYS = np.array([], dtype=int)
        NZs_LYN = np.array([], dtype=int)
        NZs_LYS_b = np.array([], dtype=int)
        NZs_LYN_b = np.array([], dtype=int)
    
    hbonds = pd.read_csv(label+"_hbonds.csv")
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)
    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    if layer != boundary:
        iWATs_b = np.load("iWATs_"+label+"_b.npy", allow_pickle=True)
        iCLs_b = np.load("iCLs_"+label+"_b.npy", allow_pickle=True)

    # Base de datos
    # HB = Hydrogen Bond
    # hb_d = formando un HB como donor
    # hb_a = formando un HB como acceptor
    # | Número de frame | Índice del átomo | Tipo de átomo | HB Cl | HB WAT(O) | HB LYS(NZ) | HB LYN(NZ) |
    data = pd.DataFrame(columns=['step', 'atom', 'name', 'hb_d_cl', 'hb_d_wat', 'hb_d_lyn', 'hb_a_wat', 'hb_a_lys', 'hb_a_lyn'])

    for step in range(0, len(traj)):
        frame = traj[step]
        aux = hbonds[hbonds['istep']==step]
        CLs = iCLs[step]
        WATs = iWATs[step]
        if layer != boundary:
            CLs_b = iCLs_b[step]
            WATs_b = iWATs_b[step]
        else:
            CLs_b = np.array([], dtype=int)
            WATs_b = np.array([], dtype=int)

        for i, atom_i in enumerate(np.concatenate((CLs, WATs, NZs_LYS, NZs_LYN)).astype(int)):
            hb_d_cl = 0
            hb_d_wat = 0
            hb_d_lyn = 0
            hb_a_wat = 0
            hb_a_lys = 0
            hb_a_lyn = 0
            if i < len(CLs):
                name = "CL"
            elif i < (len(CLs)+len(WATs)):
                name = "WAT"
            elif i < (len(CLs)+len(WATs)+len(NZs_LYS)):
                name = "LYS"
            else:
                name = "LYN"

            if atom_i in aux['donor'].values:
                aux2 = aux[aux['donor'] == atom_i]
                for index, bond in aux2.iterrows():
                    atom_j = bond['acceptor']
                    if atom_j in np.concatenate((CLs, CLs_b)).astype(int):
                        hb_d_cl += 1
                    elif atom_j in np.concatenate((WATs, WATs_b)).astype(int):
                        hb_d_wat += 1
                    else: # Por descarte, atom_j es una LYN o LYN_b
                        hb_d_lyn +=1

            if atom_i in aux['acceptor'].values:
                aux2 = aux[aux['acceptor'] == atom_i]
                for index, bond in aux2.iterrows():
                    atom_j = bond['donor']
                    if atom_j in np.concatenate((WATs, WATs_b)).astype(int):
                        hb_a_wat +=1
                    elif atom_j in np.concatenate((NZs_LYS, NZs_LYS_b)).astype(int):
                        hb_a_lys +=1
                    else:
                        hb_a_lyn +=1

            data.loc[data.shape[0]] = [step, atom_i, name, hb_d_cl, hb_d_wat, hb_d_lyn, hb_a_wat, hb_a_lys, hb_a_lyn]

    data.to_csv(label+"_coordination_hbonds.csv")
#end


def closest_atoms(traj, selection1, iterselection1, selection2):
    distances = np.zeros((len(selection1), len(traj)))
    bonds = np.zeros((len(selection1), len(traj)))

    for step in range(0, len(traj)):
        frame = traj[step]

        for index, atom1 in enumerate(selection1):
            if not (atom1 in iterselection1[step]):
                distances[index][step] = -1
                bonds[index][step] = -1
            else:
                d0 = 999.9
                b0 = -1
                for atom2 in selection2:
                    d = np.linalg.norm(frame[atom1] - frame[atom2])
                    if d < d0:
                        d0 = d
                        b0 = atom2
                distances[index][step] = d0
                bonds[index][step] = b0
    
    return distances, bonds
#end


def plot_closest_atoms(distances, bonds, index):
    Nsteps = np.shape(distances)[1]
    cmap = plt.get_cmap('jet', len(np.unique(bonds[index][bonds[index]>0])))
    cmap.set_under('gray')
    midpoints = (np.unique(bonds[index][bonds[index]>0])[1:] + np.unique(bonds[index][bonds[index]>0])[:-1]) / 2
    spacing = 100
    VariableLimits = np.concatenate(([bonds[index][bonds[index]>0].min()-spacing], midpoints, [bonds[index].max()+spacing]))
    norm = clr.BoundaryNorm(VariableLimits, ncolors=len(np.unique(bonds[index][bonds[index]>0])))

    fig, ax = plt.subplots(figsize=(13, 7))
    cax = ax.scatter(range(Nsteps), distances[index], c=bonds[index], cmap=cmap, norm=norm, s=3)
    cbar = fig.colorbar(cax, extend='min', ticks=np.unique(bonds[index][bonds[index]>0]))
    
    return fig, ax
#end


def plot_distance(traj, selection1, index1, atom2, hbond=True):
    fig, ax = plt.subplots(figsize=(10, 7))
    mask = '@'+str(selection1[index1]+1)+' @'+str(atom2+1)
    distance = pt.distance(traj, mask)
    ax.plot(range(len(traj)), distance, lw=1)
    if hbond:
        ax.plot([0, len(traj)], [3, 3], 'k--', lw=2)
        
    return fig, ax
#end


def save_distances(traj, atoms):
    m = len(atoms)
    N = (m**2-m)//2
    dmatrix = np.empty((0, N), float)
    
    for step in range(len(traj)):
        frame = traj[step]
        aux = pdist(frame[atoms])
        dmatrix = np.vstack((dmatrix, aux))
    
    np.save("distances.npy", dmatrix)
#end


###


# datadf = pd.read_csv('data.csv')
# datadf['pos'] = datadf['pos'].apply(lambda x: np.array([float(i) for i in x[1:-1].split()]))