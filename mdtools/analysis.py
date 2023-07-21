### IMPORTS ###


from mdtools.core import *
from mdtools.hbond import *


### ANALYSIS ###


def get_indices(traj, WATs, IONs, CAs, N_rings, layer=0, boundary=None,
                delta=1.0, delta_r=None, delta_z=None, offset=None,
                preselected=False, save=True, savefileWATs="iterWATs", savefileIONs="iterIONs", first=None, last=None):
    
    # layer y boundary son el número de anillos que se excluyen de sus respectivas zonas:
    #               si es 0, se selecciona todo el nanotubo
    #               si es n, se selecciona el nanotubo menos los n primeros y los n últimos anillos
    
    # CAs son los carbonos alfa de interés:
    # tiene que ser CAs del tubej si es un tubo
    #               CAs del canal si es el canal
    #               todos los CAs si es todo el bundle
    
    if boundary is None: boundary = layer
    if delta_r is None: delta_r = delta
    if delta_z is None: delta_z = delta
    if offset is None: offset = np.array([0.0, 0.0, 0.0])
    
    if first is None: first = 0
    if last is None: last = len(traj)
    
    iterWATs = []
    iterIONs = []
    if layer != boundary:
        iterWATs_b = []
        iterIONs_b = []
    
    atoms_top = get_indices_in_layer(CAs, layer)
    atoms_bot = get_indices_in_layer(CAs, Nrings-layer)
    if layer != boundary:
        atoms_top_b = get_atoms_in_layer(CAs, layer-boundary)
        atoms_bot_b = get_atoms_in_layer(CAs, Nrings-(layer-boundary))
    
    if preselected:
        auxWATs = WATs
        auxIONs = IONs
    
    for step in range(first, last):
        frame = traj[step]
        
        if preselected:
            WATs = auxWATs[step]
            IONs = auxIONs[step]
        
        # Centro de la región
        centertop = np.sum(frame[atoms_top], axis=0)/atoms_top.size
        centerbot = np.sum(frame[atoms_bot], axis=0)/atoms_bot.size
        center = (centertop + centerbot)/2
        # Radio de la región
        rtop = np.max(distance_matrix(frame[atoms_top], frame[atoms_top]))
        rbot = np.max(distance_matrix(frame[atoms_bot], frame[atoms_bot]))
        r = max(rtop, rbot)/2 + delta_r
        # Alturas máxima y mínima de la región
        zmax = np.sum(frame[atoms_top][:,2])/atoms_top.size - center[2] + delta_z
        zmin = np.sum(frame[atoms_bot][:,2])/atoms_bot.size - center[2] - delta_z
        if layer != boundary:
            zmax_b = np.sum(frame[atoms_top_b][:,2])/atoms_top_b.size - center[2] + delta_z
            zmin_b = np.sum(frame[atoms_bot_b][:,2])/atoms_bot_b.size - center[2] - delta_z
        
        # Aguas (solo los oxígenos) en la región
        aux = []
        aux_b = []
        for atom in WATs:
            xyz = frame[atom] - (center + offset)
            if (zmin < xyz[2]) and (xyz[2] < zmax) and (xyz[0]**2+xyz[1]**2 < r**2):
                aux.append(atom)
            elif layer != boundary:
                if (zmin_b < xyz[2]) and (xyz[2] < zmax_b) and (xyz[0]**2+xyz[1]**2 < r**2):
                    aux_b.append(atom)
                
        aux = np.array(aux)
        iterWATs.append(aux)
        aux_b = np.array(aux_b)
        iterWATs_b.append(aux_b)
        
        # Iones en la región
        aux = []
        aux_b = []
        for atom in IONs:
            xyz = frame[atom] - (center + offset)
            if (zmin < xyz[2]) and (xyz[2] < zmax) and (xyz[0]**2+xyz[1]**2 < r**2):
                aux.append(atom)
            elif layer != boundary:
                if (zmin_b < xyz[2]) and (xyz[2] < zmax_b) and (xyz[0]**2+xyz[1]**2 < r**2):
                    aux_b.append(atom)
                
        aux = np.array(aux)
        iterIONs.append(aux)
        aux_b = np.array(aux_b)
        iterIONs_b.append(aux_b)
    
    iterWATs = np.array(iterWATs, dtype=object)
    iterIONs = np.array(iterIONs, dtype=object)
    if layer != boundary:
        iterWATs_b = np.array(iterWATs_b, dtype=object)
        iterIONs_b = np.array(iterIONs_b, dtype=object)
    
    if save:
        np.save(savefileWATs+".npy", iterWATs)
        np.save(savefileIONs+".npy", iterIONs)
        if layer != boundary:
            np.save(savefileWATs+"_b.npy", iterWATs_b)
            np.save(savefileIONs+"_b.npy", iterIONs_b) 
#end


def get_indices_xtal(traj, WATs, IONs, CAs, N_rings, delta_r=0.0, offsets=None, # En unidades de los lattice vectors
                     save=True, savefileWATs="iterWATs", savefileIONs="iterIONs", first=None, last=None):

    # CAs son los carbonos alfa de interés:
    # tiene que ser CAs_tubej si es un tubo
    #               CAs_canal si es el canal
    
    if offsets is None: offsets = [np.array([0.0, 0.0, 0.0])]
    
    if first is None: first = 0
    if last is None: last = len(traj)
    
    iterWATs = []
    iterIONs = []
    
    layer = 0
    atoms_top = get_atoms_in_layer(CAs, layer)
    atoms_bot = get_atoms_in_layer(CAs, Nrings-layer)
    
    for step in range(first, last):
        frame = traj[step]
        lvs = np.array([frame.box[0], frame.box[1], frame.box[2]])
        
        # Centro de la región
        centertop = np.sum(frame[atoms_top], axis=0)/atoms_top.size
        centerbot = np.sum(frame[atoms_bot], axis=0)/atoms_bot.size
        center = (centertop + centerbot)/2
        # Radio de la región
        rtop = np.max(distance_matrix(frame[atoms_top], frame[atoms_top]))
        rbot = np.max(distance_matrix(frame[atoms_bot], frame[atoms_bot]))
        r = max(rtop, rbot)/2 + delta_r
        
        # Aguas (solo los oxígenos) en la región
        aux = []
        for atom in WATs:
            for offset in offsets:
                xyz = frame[atom] - (center + offset*lvs)
                if (xyz[0]**2+xyz[1]**2 < r**2):
                    aux.append(atom)
        aux = np.array(aux)
        iterWATs.append(aux)
        
        # Cloros en la región
        aux = []
        for atom in IONs:
            for offset in offsets:
                xyz = frame[atom] - (center + offset*lvs)
                if (xyz[0]**2+xyz[1]**2 < r**2):
                    aux.append(atom)
        aux = np.array(aux)
        iterIONs.append(aux)
    
    iterWATs = np.array(iterWATs, dtype=object)
    iterIONs = np.array(iterIONs, dtype=object)
    
    if save:
        np.save(savefileWATs+".npy", iterWATs)
        np.save(savefileIONs+".npy", iterIONs)
#end


def compute_distance(frame, atom_a, atom_b, lvs):
    
    a = frame[atom_a]
    b = frame[atom_b]
    ba = a - b
    
    if lvs is not None:
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
    
    # Baker-Hubbard Hydrogen Bond Identification
    
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


def analyse(p, traj, label, res_list=[], layer=0, boundary=None,
            distance_cutoff=2.5, angle_cutoff=120, first=None, last=None, xtal=False):
    
    if boundary is None: boundary = layer
    
    if first is None: first = 0
    if last is None: last = len(traj)
    
    iterWATs = np.load("iterWATs_"+label+".npy", allow_pickle=True)
    iterIONs = np.load("iterIONs_"+label+".npy", allow_pickle=True)
    bondable = get_indices_between_layers(p.bondable, layer, p.N_rings-layer-1)
    if layer != boundary:
        iterWATs_b = np.load("iterWATs_"+label+"_b.npy", allow_pickle=True)
        iterIONs_b = np.load("iterIONs_"+label+"_b.npy", allow_pickle=True)
        bondable_b = get_indices_between_layers(p.bondable, boundary, layer-1)
        bondable_b = np.concatenate(bondable_b, get_indices_between_layers(p.bondable, p.N_rings-layer, p.N_rings-boundary-1))
    
    # Base de datos de las estadísticas como una lista de diccionarios
    # | Step | Número de aguas | Número de iones | Número de puentes | Distancia media puentes |
    stats_dicts = []
    
    # Base de datos de los puentes de hidrógeno como una lista de diccionarios - complementario al grafo (abajo)
    # | Step | Donor | Hydrogen | Acceptor | Distance |
    hbonds_dicts = []
    
    # Grafo los puentes de H
    hbonds_G = nx.MultiDiGraph()
    
    for step in range(first, last):
        frame = traj[step]
        
        WATs = iterWATs[step]
        IONs = iterIONs[step]
        if layer != boundary:
            WATs_b = iterWATs_b[step]
            IONs_b = iterIONs_b[step]
        else:
            WATs_b = np.array([], dtype=int)
            IONS_b = np.array([], dtype=int)
        b = np.concatenate((WATs_b, IONS_b, bondable_b))
        
        N_hbonds = 0
        d_ave = 0.0
        
        # Buscamos los puentes de H del frame
        
        interesting_atoms = np.concatenate((WATs, IONs, bondable, b))
        triplets, distances, presence = baker_hubbard(frame, periodic=xtal,
                                                      interesting_atoms=interesting_atoms, return_distances=True,
                                                      distance_cutoff=0.1*distance_cutoff, angle_cutoff=angle_cutoff)
        for (donor, h, acceptor), d in zip(triplets[presence[0]], distances[0][presence[0]]):
            if donor in b and acceptor in b:
                continue
            mydonor = MyAtom(traj.top, p.N_rings, p.N_res, donor)
            myacceptor = MyAtom(traj.top, p.N_rings, p.N_res, acceptor)
            hbonds_G.add_edge(mydonor, myacceptor, step=step, h=h, d=10.0*d)
            hbonds_dicts.append({'step': step, 'donor': mydonor, 'h': h, 'acceptor': myacceptor, 'd': 10.0*d})
            
        # Guardamos las estadísticas
        
        if N_hbonds != 0: d_ave = d_ave/N_hbonds
        stats_dicts.append({'step': step, 'N_WATs': len(WATs), 'N_IONs': len(IONs), 'N_HBonds': N_hbonds, 'ave_dist': d_ave})
    
    # Guardamos la información
    
    pickle.dump(hbonds_G, open(label+'hbondsG.txt', 'wb'))
    hbonds_df = pd.DataFrame(hbonds_dicts)
    hbonds_df.to_csv(label+"_hbonds.csv")
    stats_df = pd.DataFrame(stats_dicts)
    stats_df.to_csv(label+"_stats.csv")
#end


### EOF ###