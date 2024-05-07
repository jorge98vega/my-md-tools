### IMPORTS ###


from mdtools.core import *


### VISUALIZATION ###


def plot_CAs(p, traj, step):
    frame = traj.slice(step, copy=False).xyz[0]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 30))

    for atom in p.CAs:
        xyz = frame[atom.index]
        ax.scatter(xyz[0], xyz[1], xyz[2], c='grey')
    for tube in range(p.N_tubes):
        for ring in range(p.N_rings):
            atoms = [atom.index for atom in p.CAs if atom.tube == tube and atom.layer == ring]
            atoms.append(atoms[0])
            for (atom1, atom2) in zip(atoms[:-1], atoms[1:]):
                xyz1 = frame[atom1]
                xyz2 = frame[atom2]
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], color='grey', linewidth=1)
#end


def plot_region(p, traj, step, CAs, WATfile=None, layer=0, delta_r=0.0, delta_z=0.0, offsets=None, lvsunits=False, wrap=False):
    frame = traj.slice(step, copy=False).xyz[0]
    lvs = traj.slice(0, copy=False).unitcell_lengths[0]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 30))

    for atom in p.CAs:
        xyz = frame[atom.index]
        ax.scatter(xyz[0], xyz[1], xyz[2], c='grey')
    for tube in range(p.N_tubes):
        for ring in range(p.N_rings):
            atoms = [atom.index for atom in p.CAs if atom.tube == tube and atom.layer == ring]
            atoms.append(atoms[0])
            for (atom1, atom2) in zip(atoms[:-1], atoms[1:]):
                xyz1 = frame[atom1]
                xyz2 = frame[atom2]
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], color='grey', linewidth=1)
    
    if WATfile is None:
        WATs = p.WATs
    else:
        iterWATs = np.load(WATfile+".npy", allow_pickle=True)
        WATs = iterWATs[step]
    for atom in WATs:
        if wrap:
            xyz = wrap_coordinates(frame[atom], lvs)
        else:
            xyz = frame[atom]
        ax.scatter(xyz[0], xyz[1], xyz[2], c='r')
    
    atoms_top = get_indices_in_layer(CAs, layer)
    atoms_bot = get_indices_in_layer(CAs, p.N_rings-layer-1)
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
    
    # Plot cylinders
    if offsets is None: offsets = [np.array([0.0, 0.0, 0.0])]
    if not lvsunits: lvs = np.array([1.0, 1.0, 1.0])
    for offset in offsets:
        z = np.linspace(zmin + center[2] + offset[2]*lvs[2], zmax + center[2] + offset[2]*lvs[2], 50)
        phi = np.linspace(0, 2*np.pi, 50)
        phi_grid, z_grid = np.meshgrid(phi, z)
        x_grid = r*np.cos(phi_grid) + center[0] + offset[0]*lvs[0]
        y_grid = r*np.sin(phi_grid) + center[1] + offset[1]*lvs[1]
        ax.plot_surface(x_grid, y_grid, z_grid, color='g', alpha=0.5)
#end


def plot_network_3D(p, traj, step, label, reslist=[], ifpath=True, layer=0, xtalcenter=None, colordict=None):
    # Aguas y cloros del canal

    iterWATs = np.load("iterWATs_"+label+".npy", allow_pickle=True)
    iterIONs = np.load("iterIONs_"+label+".npy", allow_pickle=True)

    # Puentes de hidrógeno

    hbonds = pd.read_csv(label+"_hbonds.csv")
    if ifpath:
        if os.path.isfile(label+"_longestpaths.csv"):
            paths = pd.read_csv(label+"_longestpaths.csv")
            paths['path'] = paths['path'].apply(lambda x: np.array([int(i) for i in x[1:-1].split(",")]))
        else:
            ifpath = False

    ###
    
    if colordict is None: colordict = {'HOH-O': "r", 'LYS-N': "b", 'LYN-N': "y", 'TFA-O': "m", 'TFA-F': "k"}
    
    frame = traj.slice(step, copy=False).xyz[0]
    lvs = traj.slice(0, copy=False).unitcell_lengths[0]
    
    if xtalcenter is None:
        xtal = False
    else:
        xtal = True

    WATs = iterWATs[step]
    IONs = iterIONs[step]
    if ifpath: path = paths['path'].iloc[step]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 30))

    for index in WATs:
        atom = traj.top.atom(index)
        key = str(atom.residue.name) + '-' + str(atom.element.symbol)
        xyz = frame[index]
        if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
        ax.scatter(xyz[0], xyz[1], xyz[2], c=colordict[key])
    
    resids = []
    for index in IONs:
        resids.append(traj.top.atom(index).residue.index)
    coms = []
    for resid in set(resids):
        indices = [index for index, element in enumerate(resids) if element == resid]
        com = np.array([0.0,0.0,0.0])
        for index in indices:
            com = com + frame[IONs[index]]
        coms.append(com/len(indices))
    for index in IONs:
        atom = traj.top.atom(index)
        key = str(atom.residue.name) + '-' + str(atom.element.symbol)
        xyz = frame[index]
        if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
        ax.scatter(xyz[0], xyz[1], xyz[2], c=colordict[key])
        
        com = coms[list(set(resids)).index(traj.top.atom(index).residue.index)]
        if xtal: com = wrap_coordinates(com - xtalcenter*lvs, lvs)
        if xtal and (np.abs(com[2] - xyz[2]) > lvs[2]/2):
            case = np.sign(com[2] - xyz[2])
            ax.plot3D([xyz[0], com[0]], [xyz[1], com[1]], [xyz[2]+case*lvs[2], com[2]], color="grey", linestyle=":", linewidth=1)
            ax.plot3D([xyz[0], com[0]], [xyz[1], com[1]], [xyz[2], com[2]-case*lvs[2]], color="grey", linestyle=":", linewidth=1)
        else:
            ax.plot3D([xyz[0], com[0]], [xyz[1], com[1]], [xyz[2], com[2]], color="grey", linestyle=":", linewidth=1)
        
    bondable_atoms = get_atoms_in_reslist(p.bondable, reslist)
    bondable = get_indices_between_layers(bondable_atoms, layer, p.N_rings-layer-1)
        
    for index in bondable:
        atom = traj.top.atom(index)
        key = str(atom.residue.name) + '-' + str(atom.element.symbol)
        xyz = frame[index]
        if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
        ax.scatter(xyz[0], xyz[1], xyz[2], c=colordict[key])

    for index, bond in hbonds[hbonds["step"] == step].iterrows():
        xyz1 = frame[MyAtom.from_string(bond["donor"]).index]
        if xtal: xyz1 = wrap_coordinates(xyz1 - xtalcenter*lvs, lvs)
        xyz2 = frame[MyAtom.from_string(bond["acceptor"]).index]
        if xtal: xyz2 = wrap_coordinates(xyz2 - xtalcenter*lvs, lvs)
        if xtal and (np.abs(xyz2[2] - xyz1[2]) > lvs[2]/2):
            case = np.sign(xyz2[2] - xyz1[2])
            ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2]+case*lvs[2], xyz2[2]], color='k', linewidth=1)
            ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]-case*lvs[2]], color='k', linewidth=1)
        else:
            ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], color='k', linewidth=1)
    
    if ifpath:
        for i in range(len(path)-1):
            xyz1 = frame[abs(path[i])]
            if xtal: xyz1 = wrap_coordinates(xyz1 - xtalcenter*lvs, lvs)
            xyz2 = frame[abs(path[i+1])]
            if xtal: xyz2 = wrap_coordinates(xyz2 - xtalcenter*lvs, lvs)
            if xtal and (np.abs(xyz2[2] - xyz1[2]) > lvs[2]/2):
                case = np.sign(xyz2[2] - xyz1[2])
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2]+case*lvs[2], xyz2[2]], 'k--', linewidth=2.5)
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]-case*lvs[2]], 'k--', linewidth=2.5)
            else:
                xyz2 = frame[abs(path[i+1])]
                if xtal: xyz2 = wrap_coordinates(xyz2 - xtalcenter*lvs, lvs)
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)
#end


def plot_network_2D(p, traj, step, label, reslist=[], phimin=-4, ifpath=True, tube=0, layer=0, colordict=None):
    # Aguas y cloros del canal

    iterWATs = np.load("iterWATs_"+label+".npy", allow_pickle=True)
    iterIONs = np.load("iterIONs_"+label+".npy", allow_pickle=True)

    # Puentes de hidrógeno

    hbonds = pd.read_csv(label+"_hbonds.csv")
    if ifpath:
        if os.path.isfile(label+"_paths.csv"):
            paths = pd.read_csv(label+"_paths.csv")
            paths['path'] = paths['path'].apply(lambda x: np.array([int(i) for i in x[1:-1].split(",")]))
        else:
            ifpath = False

    ###
    
    frame = traj[step]

    WATs = iterWATs[step]
    IONs = iterIONs[step]
    if ifpath: path = paths['path'].iloc[step]

    fig = plt.figure()
    ax = fig.add_subplot()
    
    bondable_atoms = get_atoms_in_reslist(p.bondable, reslist)
    bondable = get_indices_between_layers(bondable_atoms, layer, p.N_rings-layer-1)
        
    for index in bondable:
        atom = traj.top.atom(index)
        key = str(atom.residue.name) + '-' + str(atom.element)
        xyz = frame[atom]
        phi = np.arctan2(xyz[1], xyz[0])
        if phi < phimin:
            phi += 2*np.pi
        ax.scatter(phi, xyz[2], c=colordict[key])
    
    center = np.array([0.0, 0.0, 0.0])
    if tube > 0:
        atoms_center = p.CAs[int((tube-1)*p.CAs.size/p.N_tubes):int(tube*p.CAs.size/p.N_tubes)]
        center = np.sum(frame[atoms_center], axis=0)/atoms_center.size
    
    for atom in WATs:
        xyz = frame[atom]
        if tube > 0:
            xyz = xyz - center
        phi = np.arctan2(xyz[1], xyz[0])
        if phi < phimin:
            phi += 2*np.pi
        ax.scatter(phi, xyz[2], c='r')

    for atom in IONs:
        xyz = frame[atom]
        phi = np.arctan2(xyz[1], xyz[0])
        if phi < phimin:
            phi += 2*np.pi
        ax.scatter(phi, xyz[2], c='k')

    xmin = 999
    xmax = -999
    ymin = 999
    ymax = -999
    for index, bond in hbonds[hbonds["istep"] == step].iterrows():
        xyzd = frame[bond["donor"]]
        if tube > 0:
            xyzd = xyzd - center
        phid = np.arctan2(xyzd[1], xyzd[0])
        if phid < phimin:
            phid += 2*np.pi
        xyza = frame[bond["acceptor"]]
        if tube > 0:
            xyza = xyza - center
        phia = np.arctan2(xyza[1], xyza[0])
        if phia < phimin:
            phia += 2*np.pi
        if phia - phid > np.pi:
            ax.annotate("", xy=(phid, xyzd[2]), xytext=(phia-2*np.pi, xyza[2]), arrowprops=dict(arrowstyle="<|-", lw=1.2, fc='y'))
            ax.annotate("", xy=(phid+2*np.pi, xyzd[2]), xytext=(phia, xyza[2]), arrowprops=dict(arrowstyle="<|-", lw=1.2, fc='y'))
            if phia-2*np.pi < xmin: xmin = phia-2*np.pi
            if phid+2*np.pi > xmax: xmax = phid+2*np.pi
        elif phid - phia > np.pi:
            ax.annotate("", xy=(phid, xyzd[2]), xytext=(phia+2*np.pi, xyza[2]), arrowprops=dict(arrowstyle="<|-", lw=1.2, fc='y'))
            ax.annotate("", xy=(phid-2*np.pi, xyzd[2]), xytext=(phia, xyza[2]), arrowprops=dict(arrowstyle="<|-", lw=1.2, fc='y'))
            if phid-2*np.pi < xmin: xmin = phid-2*np.pi
            if phia+2*np.pi > xmax: xmax = phia+2*np.pi
        else:
            ax.annotate("", xy=(phid, xyzd[2]), xytext=(phia, xyza[2]), arrowprops=dict(arrowstyle="<|-", lw=1.2, fc='y'))
            if phid < xmin: xmin = phid
            if phia < xmin: xmin = phia
            if phid > xmax: xmax = phid
            if phia > xmax: xmax = phia
            if xyzd[2] < ymin: ymin = xyzd[2]
            if xyza[2] < ymin: ymin = xyza[2]
            if xyzd[2] > ymax: ymax = xyzd[2]
            if xyza[2] > ymax: ymax = xyza[2]
    
    if ifpath:
        for i in range(len(path)-1):
            xyz1 = frame[path[i]]
            if tube > 0:
                xyz1 = xyz1 - center
            phi1 = np.arctan2(xyz1[1], xyz1[0])
            if phi1 < phimin:
                phi1 += 2*np.pi
            xyz2 = frame[path[i+1]]
            if tube > 0:
                xyz2 = xyz2 - center
            phi2 = np.arctan2(xyz2[1], xyz2[0])
            if phi2 < phimin:
                phi2 += 2*np.pi

            if phi1 - phi2 > np.pi:
                ax.plot([phi1, phi2+2*np.pi], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)
                ax.plot([phi1-2*np.pi, phi2], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)
            elif phi2 - phi1 > np.pi:
                ax.plot([phi1, phi2-2*np.pi], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)
                ax.plot([phi1+2*np.pi, phi2], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)
            else:
                ax.plot([phi1, phi2], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)


    xticks = np.arange(-2*np.pi, 2.1*np.pi, 0.5*np.pi)
    xlabels = ["-2$\pi$", "-1.5$\pi$", "-$\pi$", "-0.5$\pi$", "0", "$0.5\pi$", "$\pi$", "$1.5\pi$", "2$\pi$"]
    plt.xticks(xticks, xlabels)
    plt.xlim(xmin-0.01, xmax+0.01)
    plt.ylim(ymin-0.01, ymax+0.01)
    ax.set_xlabel("Polar angle $\Phi$ (rad)")
    ax.set_ylabel("z ($\AA$)")
    
    return fig, ax
#end


def compute_xyz(p, traj, atoms, iterate, first, last, tube=0):
    xyz = np.array([]).reshape(0, 3)
    
    if iterate:
        aux = atoms

    for step in range(first, last):
        frame = traj[step]
        if iterate:
            atoms = aux[step]
        for atom in atoms:
            coords = frame[atom]
            if tube > 0:
                atoms_center = p.CAs[int((tube-1)*p.CAs.size/p.N_tubes):int(tube*p.CAs.size/p.N_tubes)]
                center = np.sum(frame[atoms_center], axis=0)/atoms_center.size
                coords = coords - center
            xyz = np.vstack((xyz, coords))
    
    return xyz
#end


def compute_density_profile(xyz, Nbins, axs, color, name):
    xs = xyz[:,0]
    ys = xyz[:,1]
    phis = np.arctan2(ys, xs)
    rs = np.sqrt(xs**2 + ys**2)
    zs = xyz[:,2]
    hist, bin_edges = np.histogram(phis, bins=Nbins)
    axs[0].plot((bin_edges[1:]+bin_edges[:-1])/2, hist, color=color, label=name)
    hist, bin_edges = np.histogram(rs, bins=Nbins)
    axs[1].plot((bin_edges[1:]+bin_edges[:-1])/2, hist, color=color, label=name)
    hist, bin_edges = np.histogram(zs, bins=Nbins)
    axs[2].plot((bin_edges[1:]+bin_edges[:-1])/2, hist, color=color, label=name)
#end


def plot_density_profiles(p, traj, label, first, last, canal=True, tube=0, Nbins=100, layer=0, KY=False):
    iWATs_canal = np.load("iWATs_"+label+".npy", allow_pickle=True)
    WATs_xyz = compute_xyz(p, traj, iWATs_canal, True, first, last, tube)
    #print(WATs_xyz.shape)
    if canal:
        iCLs_canal = np.load("iCLs_"+label+".npy", allow_pickle=True)
        CLs_xyz = compute_xyz(p, traj, iCLs_canal, True, first, last, tube)
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        NZs_LYN_xyz = compute_xyz(p, traj, NZs_LYN, False, first, last, tube)
        NZs_LYS_xyz = compute_xyz(p, traj, NZs_LYS, False, first, last, tube)
        if KY:
            TYDs = select_atoms(p, layer, KY)
            TYDs_xyz = compute_xyz(p, traj, TYDs, False, first, last, tube)

    fig, axs = plt.subplots(3, 1, figsize=(8, 12))

    compute_density_profile(WATs_xyz, Nbins, axs, 'r', "WAT(O)")
    if canal:
        compute_density_profile(CLs_xyz, Nbins, axs, 'k', "Cl")
        compute_density_profile(NZs_LYN_xyz, Nbins, axs, 'y', "LYN(NZ)")
        compute_density_profile(NZs_LYS_xyz, Nbins, axs, 'b', "LYS(NZ)")
        if KY:
            compute_density_profile(TYDs_xyz, Nbins, axs, 'm', "TYD(O)")

    axs[0].set_xlabel("$\phi (rad)$")
    axs[1].set_xlabel("$r (\AA)$")
    axs[2].set_xlabel("$z (\AA)$")
    axs[0].legend(edgecolor='0.75')
    axs[1].legend(edgecolor='0.75')
    axs[2].legend(edgecolor='0.75')
    
    return fig, axs
#end


def plot_water_stability(traj, file="water_stability.csv", index=None, figsize=(12, 6)):
    stab_df = pd.read_csv(file)

    indices = stab_df['index'].tolist()
    if index is None:
        index = indices[np.random.randint(len(indices))]
    else:
        if not index in indices:
            print("Error: index not present in file")
            return

    row = stab_df[stab_df['index'] == index]
    field = row['bkps']
    data = field[field.index[0]]
    bkps = [0] + json.loads(data)

    zs = np.array([traj.slice(step, copy=False).xyz[0][index][2] for step in range(len(traj))])
    window_size = 100 # del rolling average
    zseries = pd.Series(zs)
    windows = zseries.rolling(window_size, center=True)
    rolling_averages = windows.mean()
    zra = np.array(rolling_averages.tolist())

    fig, axs = rpt.display(zs, bkps, figsize=figsize)
    axs[0].plot(zra, 'r', lw=3)

    for ibkp in range(len(bkps)-1):
        i = bkps[ibkp]
        j = bkps[ibkp+1]
        axs[0].plot([i, j], [zs[i:j].mean(), zs[i:j].mean()], 'k--', lw=3)

    axs[0].set_xlim([bkps[0], bkps[-1]])
    axs[0].set_title("index = " + str(index))
    axs[0].set_xlabel("Step")
    axs[0].set_ylabel("z (nm)")
    return fig, axs, index
#end


### EOF ###