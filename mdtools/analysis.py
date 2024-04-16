### IMPORTS ###


from mdtools.core import *


### ANALYSIS ###


def get_indices(traj, WATs, IONs, CAs, N_rings, layer=0, boundary=None,
                delta=0.1, delta_r=None, delta_z=None, offset=None,
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
    atoms_bot = get_indices_in_layer(CAs, N_rings-layer-1)
    if layer != boundary:
        atoms_top_b = get_indices_in_layer(CAs, boundary)
        atoms_bot_b = get_indices_in_layer(CAs, N_rings-boundary-1)
    
    if preselected:
        auxWATs = WATs
        auxIONs = IONs
    
    for step in range(first, last):
        frame = traj.slice(step, copy=False).xyz[0]
        
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
        if layer != boundary:
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
        if layer != boundary:
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
    atoms_top = get_indices_in_layer(CAs, layer)
    atoms_bot = get_indices_in_layer(CAs, N_rings-layer-1)
    
    for step in range(first, last):
        frame = traj.slice(step, copy=False).xyz[0]
        lvs = traj.slice(0, copy=False).unitcell_lengths[0]
        
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
                xyz = wrap_coordinates(frame[atom], lvs) - (center + offset*lvs)
                if (xyz[0]**2+xyz[1]**2 < r**2):
                    aux.append(atom)
                    break
        aux = np.array(aux)
        iterWATs.append(aux)
        
        # Iones en la región
        aux = []
        for atom in IONs:
            for offset in offsets:
                xyz = wrap_coordinates(frame[atom], lvs) - (center + offset*lvs)
                if (xyz[0]**2+xyz[1]**2 < r**2):
                    aux.append(atom)
                    break
        aux = np.array(aux)
        iterIONs.append(aux)
    
    iterWATs = np.array(iterWATs, dtype=object)
    iterIONs = np.array(iterIONs, dtype=object)
    
    if save:
        np.save(savefileWATs+".npy", iterWATs)
        np.save(savefileIONs+".npy", iterIONs)
#end


def analyse(p, traj, label, reslist=[], layer=0, boundary=None,
            distance_cutoff=2.5, angle_cutoff=120, first=None, last=None, xtal=False):
    
    if boundary is None: boundary = layer
    
    if first is None: first = 0
    if last is None: last = len(traj)
    
    iterWATs = np.load("iterWATs_"+label+".npy", allow_pickle=True)
    iterIONs = np.load("iterIONs_"+label+".npy", allow_pickle=True)
    bondable_atoms = get_atoms_in_reslist(p.bondable, reslist)
    bondable = get_indices_between_layers(bondable_atoms, layer, p.N_rings-layer-1)
    backbone = get_indices_between_layers(np.concatenate((p.bbNs, p.bbOs)), layer, p.N_rings-layer-1)
    if layer != boundary:
        iterWATs_b = np.load("iterWATs_"+label+"_b.npy", allow_pickle=True)
        iterIONs_b = np.load("iterIONs_"+label+"_b.npy", allow_pickle=True)
        bondable_b = get_indices_between_layers(bondable_atoms, boundary, layer-1)
        bondable_b = np.concatenate((bondable_b, get_indices_between_layers(bondable_atoms, p.N_rings-layer, p.N_rings-boundary-1)))
        backbone_b = get_indices_between_layers(np.concatenate((p.bbNs, p.bbOs)), boundary, layer-1)
        backbone_b = np.concatenate((backbone_b, get_indices_between_layers(np.concatenate((p.bbNs, p.bbOs)), p.N_rings-layer, p.N_rings-boundary-1)))
    else:
        bondable_b = np.array([], dtype=int)
        backbone_b = np.array([], dtype=int)
    
    # Base de datos de las estadísticas como una lista de diccionarios
    # | Step | Número de aguas | Número de iones | Número de puentes | Distancia media puentes |
    stats_dicts = []
    
    # Base de datos de los puentes de hidrógeno como una lista de diccionarios - complementario al grafo (abajo)
    # | Step | Donor | Hydrogen | Acceptor | Distance |
    hbonds_dicts = []
    
    # Grafo los puentes de H
    hbonds_G = nx.MultiDiGraph()
    
    for step in range(first, last):
        frame = traj.slice(step, copy=False)
        
        WATs = iterWATs[step]
        IONs = iterIONs[step]
        if layer != boundary:
            WATs_b = iterWATs_b[step]
            IONs_b = iterIONs_b[step]
        else:
            WATs_b = np.array([], dtype=int)
            IONs_b = np.array([], dtype=int)
        b = np.concatenate((WATs_b, IONs_b, bondable_b, backbone_b))
        
        N_hbonds = 0
        d_ave = 0.0
        
        # Buscamos los puentes de H del frame
        
        interesting_atoms = np.concatenate((WATs, IONs, bondable, backbone, b))
        triplets, distances, presence = md.baker_hubbard(frame, periodic=xtal,
                                                         interesting_atoms=interesting_atoms, return_distances=True,
                                                         distance_cutoff=0.1*distance_cutoff, angle_cutoff=angle_cutoff)
        
        # Evitar el conteo doble
        
        ignore_indices = []
        u, c = np.unique(triplets[presence[0]][:, 1], return_counts=True) 
        for duplicate in u[c > 1]:
            indices, = np.where(triplets[presence[0]][:, 1] == duplicate) # índices en "triplets[presence[0]][:,1]", "distances[0]"
            dmin_index = np.argmin(distances[0][presence[0]][indices]) # índice en "indices"
            ignore_indices += [index for index in np.delete(indices, dmin_index)] # índices en "triplets[presence[0]][:,1]", "distances[0]"
            
        # Guardamos los hbonds
        
        for index, ((donor, h, acceptor), d) in enumerate(zip(triplets[presence[0]], distances[0][presence[0]])):
            if index in ignore_indices:
                continue
            if donor in b and acceptor in b:
                continue
            if donor in backbone and acceptor in backbone:
                continue
            mydonor = MyAtom(traj.top, p.N_rings, p.N_res, donor)
            myacceptor = MyAtom(traj.top, p.N_rings, p.N_res, acceptor)
            hbonds_G.add_edge(donor, acceptor, step=step, h=h, d=10.0*d)
            hbonds_dicts.append({'step': step, 'donor': mydonor, 'h': h, 'acceptor': myacceptor, 'd': 10.0*d})
            N_hbonds += 1
            d_ave += d
            
        # Guardamos las estadísticas
        
        if N_hbonds != 0: d_ave = d_ave/N_hbonds
        stats_dicts.append({'step': step, 'N_WATs': len(WATs), 'N_IONs': len(IONs), 'N_HBonds': N_hbonds, 'ave_dist': 10.0*d_ave})
    
    # Guardamos la información
    
    pickle.dump(hbonds_G, open(label+'_hbondsG.dat', 'wb'))
    hbonds_df = pd.DataFrame(hbonds_dicts)
    hbonds_df.to_csv(label+"_hbonds.csv")
    stats_df = pd.DataFrame(stats_dicts)
    stats_df.to_csv(label+"_stats.csv")
#end


def detail_hbonds(label):
    hbonds_df = pd.read_csv(label+"_hbonds.csv")
    Nsteps = hbonds_df['step'].max()+1
    
    # Base de datos de los puentes de H del canal como una lista de diccionarios
    # | Step | Donor | Acceptor | Número de puentes | Distancia media |
    detail_dicts = []
    
    for step in range(Nsteps):
        aux_df = hbonds_df[hbonds_df["step"] == step]
        atoms = []
        nhbonds = []
        dists = []
        for index, hbond in aux_df.iterrows():
            donor = MyAtom.from_string(hbond['donor'])
            acceptor = MyAtom.from_string(hbond['acceptor'])
            atoms_dict = {'donor': donor.resname + "-" + re.sub(r'\d+', '', donor.name),
                          'acceptor': acceptor.resname + "-" + re.sub(r'\d+', '', acceptor.name)}
            if atoms_dict not in atoms:
                atoms.append(atoms_dict)
                nhbonds.append(1)
                dists.append(hbond['d'])
            else:
                index = atoms.index(atoms_dict)
                nhbonds[index] += 1
                dists[index] += hbond['d']
        for (pair, nhb, d) in zip(atoms, nhbonds, dists):
            aux_dict = {'step': step, 'donor': pair['donor'], 'acceptor': pair['acceptor'],
                        'N_HBonds': nhb, 'd': d/nhb}
            detail_dicts.append(aux_dict)
                
    
    detail_df = pd.DataFrame(detail_dicts)
    detail_df.to_csv(label+"_detail.csv")
#end


def search_longestpaths(traj, label, xtal=False, first=None, last=None):
    hbondsG = pickle.load(open(label + '_hbondsG.dat', 'rb'))
    if first is None: first = 0
    if last is None: last = len(traj)

    paths_dicts = []
    for step in range(first, last):
        frame = traj.slice(step, copy=False).xyz[0]
        if xtal: lvs = traj.slice(0, copy=False).unitcell_lengths[0] # nm
        auxG = nx.MultiDiGraph(((u,v,d) for u,v,d in hbondsG.edges(data=True) if d['step'] == step))
        longest_path = {'step': step, 'path': [], 'residues': [], 'dz': 0.0}

        # Paths
        paths = dict(nx.all_pairs_shortest_path(auxG))
        for node1 in paths:
            for node2 in paths[node1]:
                path = paths[node1][node2]
                if xtal:
                    totaldz = 0.0
                    for i in range(len(path)-1):
                        dz = 10.0*(frame[path[i+1]][2] - frame[path[i]][2]) # Å
                        if abs(dz) > 10.0*lvs[2]/2: dz = dz - np.sign(dz)*10.0*lvs[2]
                        totaldz += dz
                    dz = totaldz
                else:
                    dz = 10.0*(frame[node2][2] - frame[node1][2])
                if abs(dz) > abs(longest_path['dz']):
                    longest_path['path'] = path
                    longest_path['residues'] = [traj.top.atom(node).residue.name for node in path]
                    longest_path['dz'] = dz

        # Cycles
        if xtal:
            cycles = sorted(nx.simple_cycles(auxG))
            for cycle in cycles:
                cycle.append(cycle[0])
                totaldz = 0.0
                for i in range(len(cycle)-1):
                    dz = 10.0*(frame[cycle[i+1]][2] - frame[cycle[i]][2]) # Å
                    if abs(dz) > 10.0*lvs[2]/2: dz = dz - np.sign(dz)*10.0*lvs[2]
                    totaldz += dz
                dz = totaldz
                if abs(dz) > abs(longest_path['dz']):
                    longest_path['path'] = cycle
                    longest_path['residues'] = [traj.top.atom(node).residue.name for node in cycle]
                    longest_path['dz'] = dz

        paths_dicts.append(longest_path)
        
    paths_df = pd.DataFrame(paths_dicts)
    paths_df.to_csv(label+"_longestpaths.csv")
#end


def search_paths_res(traj, label, resnamelist, first=None, last=None):
    
    def search_neighbors_res(G, nodes, resnamelist, path, paths_dicts):
        if len(resnamelist) == 0:
            paths_dicts.append({'step': step, 'path': path})
            return
        for node in nodes:
            if traj.top.atom(node).residue.name != resnamelist[0]: continue
            aux_path = path.copy()
            aux_path.append(node)
            neighbors = G.neighbors(node)
            search_neighbors_res(G, neighbors, resnamelist[1:], aux_path, paths_dicts)
    
    hbondsG = pickle.load(open(label + '_hbondsG.dat', 'rb'))
    if first is None: first = 0
    if last is None: last = len(traj)
    reslabel = ""
    for resname in resnamelist: reslabel += resname + "-"
    reslabel = reslabel[:-1]
    
    paths_dicts = []
    for step in range(first, last):
        frame = traj.slice(step, copy=False).xyz[0]
        auxG = nx.MultiDiGraph(((u,v,d) for u,v,d in hbondsG.edges(data=True) if d['step'] == step))
        nodes = list(auxG.nodes)
        search_neighbors_res(auxG, nodes, resnamelist, [], paths_dicts)
    
    paths_df = pd.DataFrame(paths_dicts)
    paths_df.to_csv(label+"_paths_"+reslabel+".csv")
#end


### EOF ###
