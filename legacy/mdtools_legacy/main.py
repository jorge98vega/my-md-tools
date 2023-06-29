### IMPORTS ###


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d   
import numpy as np
from scipy.spatial import distance_matrix
import pytraj as pt
import mdtraj as md
import pandas as pd
import math
import time


### MAIN ###


def orientar(punto, p0, pZ, pX):
    '''
    Cambia las coordenadas de "punto" a un nuevo sistema de referencia en el que "p0" es el 0,
    el eje Z está orientado de "p0" a "pZ", y el eje X está orientado de "p0" al punto "pX" proyectado
    sobre el plano perpendicular al nuevo eje Z
    '''
    
    vZ = pZ - p0 # Vector de p0 a pZ
    uZ = vZ/np.linalg.norm(vZ) # Vector del nuevo eje Z
    
    wX = pX - p0 # Vector de p0 a pX
    vX = wX - np.dot(wX, uZ)*uZ # wX proyectado sobre el plano perpendicular a uZ
    uX = vX/np.linalg.norm(vX) # Vector del nuevo eje X
    
    uY = np.cross(uZ, uX) # Vector del nuevo eje Y
    
    vp = punto - p0 # Vector de p0 a punto
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
        xyz.append(orientar(traj.xyz[step][atom], p0, pZ, pX))

    oriented = md.Trajectory(xyz, traj.top.subset(selection))
    
    traj.superpose(oriented, 0, atom_indices=selection, ref_atom_indices=range(oriented.n_atoms))
    traj.save(run_name+"_RMSD.nc")
#end


class Params:
    def __init__(self, traj, N_tubes, N_rings, N_res, channel_res_1, channel_res_2, channel_res_3, channel_res_4, channel_res_list=None):
        
        self.N_tubes = N_tubes
        self.N_rings = N_rings
        self.N_res = N_res
        N_allres = N_tubes*N_rings*N_res
        self.N_allres = N_allres

        CAs = traj.top.select("@CA")
        self.CAs = CAs
        CAs_tube1 = CAs[0:int(CAs.size/N_tubes)]
        self.CAs_tube1 = CAs_tube1
        CAs_tube2 = CAs[int(CAs.size/N_tubes):int(2*CAs.size/N_tubes)]
        self.CAs_tube2 = CAs_tube2
        CAs_tube3 = CAs[int(2*CAs.size/N_tubes):int(3*CAs.size/N_tubes)]
        self.CAs_tube3 = CAs_tube3
        CAs_tube4 = CAs[int(3*CAs.size/N_tubes):int(4*CAs.size/N_tubes)]
        self.CAs_tube4 = CAs_tube4
        CAs_top = np.concatenate((CAs_tube1[0:N_res], CAs_tube2[0:N_res], CAs_tube3[0:N_res], CAs_tube4[0:N_res]))
        CAs_bot = np.concatenate((CAs_tube1[-N_res:], CAs_tube2[-N_res:], CAs_tube3[-N_res:], CAs_tube4[-N_res:]))
        CAs_topbot = np.concatenate((CAs_top, CAs_bot))
        self.CAs_topbot = CAs_topbot

        bbNs = {}
        bbNs["canal"] = traj.top.select(":1-"+str(N_allres)+" & @N")
        bbNs["tube1"] = bbNs["canal"][0:int(bbNs["canal"].size/N_tubes)]
        bbNs["tube2"] = bbNs["canal"][int(bbNs["canal"].size/N_tubes):int(2*bbNs["canal"].size/N_tubes)]
        bbNs["tube3"] = bbNs["canal"][int(2*bbNs["canal"].size/N_tubes):int(3*bbNs["canal"].size/N_tubes)]
        bbNs["tube4"] = bbNs["canal"][int(3*bbNs["canal"].size/N_tubes):int(4*bbNs["canal"].size/N_tubes)]
        self.bbNs = bbNs
        bbOs = {}
        bbOs["canal"] = traj.top.select(":1-"+str(N_allres)+" & @O")
        bbOs["tube1"] = bbOs["canal"][0:int(bbOs["canal"].size/N_tubes)]
        bbOs["tube2"] = bbOs["canal"][int(bbOs["canal"].size/N_tubes):int(2*bbOs["canal"].size/N_tubes)]
        bbOs["tube3"] = bbOs["canal"][int(2*bbOs["canal"].size/N_tubes):int(3*bbOs["canal"].size/N_tubes)]
        bbOs["tube4"] = bbOs["canal"][int(3*bbOs["canal"].size/N_tubes):int(4*bbOs["canal"].size/N_tubes)]
        self.bbOs = bbOs
        
        channel_res = {}
        if channel_res_list is None:
            channel_res["0"] = channel_res_1
            channel_res["1"] = channel_res_2
            channel_res["2"] = channel_res_3
            channel_res["3"] = channel_res_4
        else: # channel_res_list = [[5, 5], [3, 7], [1, 1], [7, 3]]
            for tube in range(4):
                channel_res[str(tube)] = ""
                while len(channel_res_list[tube]) >= 2:
                    for ring in range(N_rings):
                        channel_res[str(tube)] += str(tube*N_rings*N_res + ring*N_res + channel_res_list[tube][ring%2])
                        channel_res[str(tube)] += ","
                    channel_res_list[tube] = channel_res_list[tube][2:]
        
        CAs_canal1 = traj.top.select(":"+channel_res["0"]+" & :LYS,LYN & @CA") # Solo los CAs de LYS/LYN
        LYS1 = traj.top.select(":"+channel_res["0"]+" & :LYS & @NZ")
        self.LYS1 = LYS1
        LYN1 = traj.top.select(":"+channel_res["0"]+" & :LYN & @NZ")
        self.LYN1 = LYN1
        TYD1 = traj.top.select(":"+channel_res["0"]+" & :TYD & @OH")
        self.TYD1 = TYD1
        CAs_canal2 = traj.top.select(":"+channel_res["1"]+" & :LYS,LYN & @CA")
        LYS2 = traj.top.select(":"+channel_res["1"]+" & :LYS & @NZ")
        self.LYS2 = LYS2
        LYN2 = traj.top.select(":"+channel_res["1"]+" & :LYN & @NZ")
        self.LYN2 = LYN2
        TYD2 = traj.top.select(":"+channel_res["1"]+" & :TYD & @OH")
        self.TYD2 = TYD2
        CAs_canal3 = traj.top.select(":"+channel_res["2"]+" & :LYS,LYN & @CA")
        LYS3 = traj.top.select(":"+channel_res["2"]+" & :LYS & @NZ")
        self.LYS3 = LYS3
        LYN3 = traj.top.select(":"+channel_res["2"]+" & :LYN & @NZ")
        self.LYN3 = LYN3
        TYD3 = traj.top.select(":"+channel_res["2"]+" & :TYD & @OH")
        self.TYD3 = TYD3
        CAs_canal4 = traj.top.select(":"+channel_res["3"]+" & :LYS,LYN & @CA")
        LYS4 = traj.top.select(":"+channel_res["3"]+" & :LYS & @NZ")
        self.LYS4 = LYS4
        LYN4 = traj.top.select(":"+channel_res["3"]+" & :LYN & @NZ")
        self.LYN4 = LYN4
        TYD4 = traj.top.select(":"+channel_res["3"]+" & :TYD & @OH")
        self.TYD4 = TYD4
        
        # El canal es como un nanotubo en el que N_res = N_tubes
        CAs_canal = np.array([], dtype=int)
        for i in range(N_rings):
            for selection in (CAs_canal1, CAs_canal2, CAs_canal3, CAs_canal4):
                if len(selection) != 0:
                    CAs_canal = np.append(CAs_canal, selection[i])
        self.CAs_canal = CAs_canal
        
        # 1er CA de la primera capa de un nanotubo: CAs_tubej[0]
        # 1er CA de la segunda capa de un nanotubo: CAs_tubej[N_res] ([N_tubes] en el caso tubej = canal)
        # ...
        # 1er CA de la última capa de un nanotubo: CAs_tubej[-Nres]
        # Los átomos (menos los N y H del backbone) cuyos índices están entre estos dos valores, perteneces a las capas intermedias (2-penúltima)
#         NZs_LYS = np.array([], dtype=int)
#         for NZs, CAs in zip((LYS1, LYS2, LYS3, LYS4), (CAs_tube1, CAs_tube2, CAs_tube3, CAs_tube4)):
#             for atom in NZs:
#                 if (CAs[N_res] < atom) and (atom < CAs[-N_res]):
#                     NZs_LYS = np.append(NZs_LYS, atom)
#         self.NZs_LYS = NZs_LYS
#         NZs_LYN = np.array([], dtype=int)
#         for NZs, CAs in zip((LYN1, LYN2, LYN3, LYN4), (CAs_tube1, CAs_tube2, CAs_tube3, CAs_tube4)):
#             for atom in NZs:
#                 if (CAs[N_res] < atom) and (atom < CAs[-N_res]):
#                     NZs_LYN = np.append(NZs_LYN, atom)
#         self.NZs_LYN = NZs_LYN
        
        WATs = traj.top.select(":WAT&@O")
        self.WATs = WATs
        CLs = traj.top.select("@Cl-")
        self.CLs = CLs
    #end
#end


def select_atoms(p, layer, KY=False):
    if KY:
        TYDs = np.array([], dtype=int)
        for OHs, CAs in zip((p.TYD1, p.TYD2, p.TYD3, p.TYD4), (p.CAs_tube1, p.CAs_tube2, p.CAs_tube3, p.CAs_tube4)):
            for atom in OHs:
                lastCA = CAs[-layer*p.N_res] if layer != 0 else atom+1
                if (CAs[layer*p.N_res] < atom) and (atom < lastCA):
                    TYDs = np.append(TYDs, atom)
        return TYDs
    
    NZs_LYS = np.array([], dtype=int)
    for NZs, CAs in zip((p.LYS1, p.LYS2, p.LYS3, p.LYS4), (p.CAs_tube1, p.CAs_tube2, p.CAs_tube3, p.CAs_tube4)):
        for atom in NZs:
            lastCA = CAs[-layer*p.N_res] if layer != 0 else atom+1
            if (CAs[layer*p.N_res] < atom) and (atom < lastCA):
                NZs_LYS = np.append(NZs_LYS, atom)

    NZs_LYN = np.array([], dtype=int)
    for NZs, CAs in zip((p.LYN1, p.LYN2, p.LYN3, p.LYN4), (p.CAs_tube1, p.CAs_tube2, p.CAs_tube3, p.CAs_tube4)):
        for atom in NZs:
            lastCA = CAs[-layer*p.N_res] if layer != 0 else atom+1
            if (CAs[layer*p.N_res] < atom) and (atom < lastCA):
                NZs_LYN = np.append(NZs_LYN, atom)
    
    return NZs_LYS, NZs_LYN
#end


def layer_of(p, atom):
    for CAs in (p.CAs_tube1, p.CAs_tube2, p.CAs_tube3, p.CAs_tube4):
        for layer in range(p.N_rings):
            lastCA = CAs[(layer+1)*p.N_res] if (layer+1) != p.N_rings else (2*CAs[layer*p.N_res] - CAs[(layer-1)*p.N_res])
            if (CAs[layer*p.N_res] < atom) and (atom < lastCA):
                return layer
#end


def wrap_coordinates(point, box):
    
    for dim in range(3):
        point[dim] = ((point[dim] + box[dim]/2) % box[dim]) - box[dim]/2
    
    return point
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