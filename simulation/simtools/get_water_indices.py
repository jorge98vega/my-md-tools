#!/home/jorge/AMBER/ambertools21/amber20_src/miniconda/bin/python3

import sys
import os.path
import numpy as np
import pytraj as pt
from scipy.spatial import distance_matrix

# Uso: ./get_water_indices.py 4t10s_run03.top 4t10s_run03_MD.nc 4t10s_run03_RMSD.nc [4990] [5000]
input_top = str(sys.argv[1]) # Topología del sistema
input_traj = str(sys.argv[2]) # Trayectoria previa al fit RMSD
output_traj = str(sys.argv[3]) # Trayectoria post fit RMSD
if len(sys.argv)>4:
    first = int(sys.argv[4]) # Primer frame a analizar (incluído)
else:
    first = None
if len(sys.argv)>5:
    last = int(sys.argv[5]) # Último frame a analizar (no incluído)
else:
    last = None

N_tubes = 4 # Número de tubos en el sistema
N_rings = 10 # Número de anillos en un tubo
N_res = 8 # Número de residuos en un anillo
channel_res_1 = "5, 13, 21, 29, 37, 45, 53, 61, 69, 77"
channel_res_2 = "83, 95, 99, 111, 115, 127, 131, 143, 147, 159"
channel_res_3 = "161, 169, 177, 185, 193, 201, 209, 217, 225, 233"
channel_res_4 = "247, 251, 263, 267, 279, 283, 295, 299, 311, 315"


class AlphaCarbons:
    def __init__(self, traj, N_tubes, N_rings, N_res, channel_res_1, channel_res_2, channel_res_3, channel_res_4):
        CAs = traj.top.select("@CA")
        self.all = CAs
        CAs_tube1 = CAs[0:int(CAs.size/N_tubes)]
        self.tube1 = CAs_tube1
        CAs_tube2 = CAs[int(CAs.size/N_tubes):int(2*CAs.size/N_tubes)]
        self.tube2 = CAs_tube2
        CAs_tube3 = CAs[int(2*CAs.size/N_tubes):int(3*CAs.size/N_tubes)]
        self.tube3 = CAs_tube3
        CAs_tube4 = CAs[int(3*CAs.size/N_tubes):int(4*CAs.size/N_tubes)]
        self.tube4 = CAs_tube4
        CAs_top = np.concatenate((CAs_tube1[0:N_res], CAs_tube2[0:N_res], CAs_tube3[0:N_res], CAs_tube4[0:N_res]))
        self.top = CAs_top
        CAs_bot = np.concatenate((CAs_tube1[-N_res:], CAs_tube2[-N_res:], CAs_tube3[-N_res:], CAs_tube4[-N_res:]))
        self.bot = CAs_bot
        CAs_topbot = np.concatenate((CAs_top, CAs_bot))
        self.topbot = CAs_topbot
        
        # Solo los CAs de LYS/LYN
        CAs_canal1 = traj.top.select(":"+channel_res_1+" & :LYS,LYN & @CA")
        CAs_canal2 = traj.top.select(":"+channel_res_2+" & :LYS,LYN & @CA")
        CAs_canal3 = traj.top.select(":"+channel_res_3+" & :LYS,LYN & @CA")
        CAs_canal4 = traj.top.select(":"+channel_res_4+" & :LYS,LYN & @CA")
        # El canal es como un nanotubo en el que N_res = N_tubes
        CAs_canal = np.array([], dtype=int)
        for i in range(N_rings):
            for selection in (CAs_canal1, CAs_canal2, CAs_canal3, CAs_canal4):
                CAs_canal = np.append(CAs_canal, selection[i])
        self.canal = CAs_canal
#end


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


def reorientate_traj_RMSD(traj, CAs, output_trajname):
    '''
    Vamos a tomar:
    como p0 el centro de masa de los carbonos alfas de los 4 tubos
    como pZ el centro de masa de los carbonos alfas de los primeros anillos de los nanotubos
    como pX el centro de masa de los carbonos alfas del tubo 1
    '''
    
    step = len(traj)-1
    p0 = np.sum(traj[step][CAs.all], axis=0)/CAs.all.size
    pZ = np.sum(traj[step][CAs.top], axis=0)/CAs.top.size
    pX = np.sum(traj[step][CAs.tube1], axis=0)/CAs.tube1.size
    
    selection = CAs.all
    xyz = []
    for atom in range(traj.n_atoms):
        if atom in selection:
            xyz.append(orientar(traj[step][atom], p0, pZ, pX))
        else:
            xyz.append(traj[step][atom])
    xyz = np.array([xyz])
    
    oriented = pt.Trajectory(top=traj.top, xyz=xyz)
    
    traj = traj.superpose(mask="@CA", ref=oriented[0], ref_mask="@CA")
    traj.save(filename=output_trajname)
    
    return traj
#end


def get_indices(traj, WATs, CAs, N_res, layer=0, first=None, last=None,
                delta=0.0, deltaR=None, deltaZ=None, offset=None,
                preselected=False, save=True, savefileWATs="iWATs"):
    
    '''
    layer es el número de anillos que se excluyen de sus respectivas zonas:
        si es 0, se selecciona todo el nanotubo
        si es n, se selecciona el nanotubo menos los n primeros y los n últimos anillos
    CAs son los carbonos alfa de interés:
    tiene que ser CAs_tubej si es un tubo
    _____________ CAs_canal si es el canal
    _____________ CAs_topbot si es todo el bundle
    N_res es el número de residuos que hay en una "capa" por lo que:
    tiene que ser N_res si es un tubo
    _____________ N_tubes si es el canal
    _____________ N_res*N_tubes si es todo el bundle
    layer tiene que ser 0 si es todo el bundle
    '''
    
    if first is None:
        first = 0
    if last is None:
        last = len(traj)
    Nsteps = last - first
    
    if deltaR is None:
        deltaR = delta
    if deltaZ is None:
        deltaZ = delta
    
    if offset is None:
        offset = np.array([0.0, 0.0, 0.0])
    
    iWATs = []
    
    atoms_top = CAs[layer*N_res:(layer+1)*N_res]
    atoms_bot = CAs[-(layer+1)*N_res:len(CAs)-layer*N_res]
    
    if preselected:
        auxWATs = WATs
    
    for step in range(0, Nsteps):
        frame = traj[first+step]
        
        if preselected:
            WATs = auxWATs[step]
        
        # Centro de la región
        centertop = np.sum(frame[atoms_top], axis=0)/atoms_top.size
        centerbot = np.sum(frame[atoms_bot], axis=0)/atoms_bot.size
        center = (centertop + centerbot)/2 + offset
        # Radio de la región
        Rtop = np.max(distance_matrix(frame[atoms_top], frame[atoms_top]))
        Rbot = np.max(distance_matrix(frame[atoms_bot], frame[atoms_bot]))
        R = max(Rtop, Rbot)/2 + deltaR
        # Alturas máxima y mínima de la región
        Zmax = np.sum(frame[atoms_top][:,2])/atoms_top.size - center[2] + deltaZ
        Zmin = np.sum(frame[atoms_bot][:,2])/atoms_bot.size - center[2] - deltaZ
        
        # Aguas (solo los oxígenos) en la región
        aux = []
        for atom in WATs:
            xyz = frame[atom] - center
            if (Zmin < xyz[2]) and (xyz[2] < Zmax) and (xyz[0]**2+xyz[1]**2 < R**2):
                aux.append(atom)
        aux = np.array(aux)
        iWATs.append(aux)
    
    if save:
        with open(savefileWATs+'.dat', 'w') as f:
            for lista in iWATs:
                for element in lista:
                    f.write(str(element)+' ')
                f.write('\n')
        with open(savefileWATs+'_res.dat', 'w') as f:
            for lista in iWATs:
                for element in lista:
                    f.write(str(traj.top[element].resid+1)+' ')
                f.write('\n')
        if preselected:
            with open(savefileWATs+'_indexres.dat', 'w') as f:
                for step, lista in enumerate(iWATs):
                    auxres = []
                    for atom in auxWATs[step]:
                        auxres.append(traj.top[atom].resid+1)
                    for element in lista:
                        f.write(str(auxres.index(traj.top[element].resid+1)+321)+' ')
                    f.write('\n')
        iWATs = np.array(iWATs, dtype=object)
        np.save(savefileWATs+'.npy', iWATs)
#end


if not os.path.isfile(output_traj):
    print("Orientando trayectoria")
    traj = pt.iterload(input_traj, input_top)
    CAs = AlphaCarbons(traj, N_tubes, N_rings, N_res, channel_res_1, channel_res_2, channel_res_3, channel_res_4)
    traj = reorientate_traj_RMSD(traj, CAs, output_traj)
else:
    traj = pt.iterload(output_traj, input_top)
    CAs = AlphaCarbons(traj, N_tubes, N_rings, N_res, channel_res_1, channel_res_2, channel_res_3, channel_res_4)
WATs = traj.top.select(":WAT&@O")

# Para todo el bundle
print("Analizando el bundle")
get_indices(traj, WATs, CAs.topbot, N_res*N_tubes, delta=5.0, first=first, last=last)
# Para el canal
print("Analizando el canal")
iWATs = np.load("iWATs.npy", allow_pickle=True)
get_indices(traj, iWATs, CAs.canal, N_tubes, delta=-1.0, preselected=True, savefileWATs="iWATs_canal", first=first, last=last)
# Para un tubo
print("Analizando el tubo 1")
iWATs = np.load("iWATs.npy", allow_pickle=True)
get_indices(traj, iWATs, CAs.tube1, N_res, preselected=True, savefileWATs="iWATs_tube1", first=first, last=last)
