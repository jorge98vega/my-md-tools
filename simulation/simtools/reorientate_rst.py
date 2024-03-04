#!/home/jorge/AMBER/ambertools21/amber20_src/miniconda/bin/python3

import sys
import os
import numpy as np
import mdtraj as md


'''
Para que el .rst que se genera funcione en vmd, hay que borrar
el segundo número de la segunda fila (creo)
'''


N_tubes = 9
N_res = 8
formato = str(sys.argv[1]) # "rst" o "pdb"
folder = str(sys.argv[2])+"/" # directorio donde está la topología y el rst/pdb
traj_name = str(sys.argv[3]) # nombre de la topología sin .top/.pdb
if len(sys.argv) > 4: # se puede dejar el sufijo en blanco
    rst_name = "_" + str(sys.argv[4]) # sufijo (entre traj_name_ y .rst) del rst/pdb
else:
    rst_name = ""


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


def reorientate_rst(formato, folder, traj_name, rst_name, N_tubes, N_res):
    '''
    N_tubes es el número de tubos en el sistema
    N_res es el número de residuos en un anillo
    
    Vamos a tomar:
    como p0 el centro de masa de los carbonos alfas de los 4 tubos
    como pZ el centro de masa de los carbonos alfas de los primeros anillos de los nanotubos
    como pX el centro de masa de los carbonos alfas del tubo 1
    '''
    if formato == "pdb":
        traj = md.load(folder+traj_name+rst_name+".pdb")
    elif formato == "rst":
        traj = md.load(folder+traj_name+rst_name+".rst7", top=folder+traj_name+".parm7")
    
    CAs = traj.top.select("name==CA")
    CAs_tubes = []
    CAs_top = np.array([], dtype=int)
    CAs_bot = np.array([], dtype=int)
    for i in range(N_tubes):
        CAs_tube_i = CAs[int(i*CAs.size/N_tubes):int((i+1)*CAs.size/N_tubes)]
        CAs_tubes.append(CAs_tube_i)
        CAs_top = np.concatenate((CAs_top, CAs_tube_i[0:N_res]))
        CAs_bot = np.concatenate((CAs_bot, CAs_tube_i[-N_res:]))
    
    step = 0 # len(traj)-1
    p0 = np.sum(traj.xyz[step][CAs], axis=0)/CAs.size
    pZ = np.sum(traj.xyz[step][CAs_top], axis=0)/CAs_top.size
    pX = np.sum(traj.xyz[step][CAs_tubes[0]], axis=0)/CAs_tubes[0].size
    
    selection = CAs
    xyz = []
    for atom in selection:
        xyz.append(orientar(traj.xyz[step][atom], p0, pZ, pX))

    oriented = md.Trajectory(xyz, traj.top.subset(selection))
    traj.superpose(oriented, 0, atom_indices=selection, ref_atom_indices=range(oriented.n_atoms))
    if formato == "pdb":
        traj.save_pdb(folder+traj_name+rst_name+"_oriented.pdb")
    elif formato == "rst":
        traj.save_amberrst7(folder+traj_name+rst_name+"_oriented.rst")
#end


if formato == "rst":
    os.system("cp "+folder+traj_name+".top "+folder+traj_name+".parm7")
    os.system("cp "+folder+traj_name+rst_name+".rst "+folder+traj_name+rst_name+".rst7")
reorientate_rst(formato, folder, traj_name, rst_name, N_tubes, N_res)
if formato == "rst":
    os.remove(folder+traj_name+".parm7")
    os.remove(folder+traj_name+rst_name+".rst7")
