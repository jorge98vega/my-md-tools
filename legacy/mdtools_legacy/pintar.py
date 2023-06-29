### IMPORTS ###


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import numpy as np
import pytraj as pt

from main import *


### PINTAR ###


# def select_atoms(p, layer):
#     NZs_LYS = np.array([], dtype=int)
#     for NZs, CAs in zip((p.LYS1, p.LYS2, p.LYS3, p.LYS4), (p.CAs_tube1, p.CAs_tube2, p.CAs_tube3, p.CAs_tube4)):
#         for atom in NZs:
#             lastCA = CAs[-layer*p.N_res] if layer != 0 else atom+1
#             if (CAs[layer*p.N_res] < atom) and (atom < lastCA):
#                 NZs_LYS = np.append(NZs_LYS, atom)

#     NZs_LYN = np.array([], dtype=int)
#     for NZs, CAs in zip((p.LYN1, p.LYN2, p.LYN3, p.LYN4), (p.CAs_tube1, p.CAs_tube2, p.CAs_tube3, p.CAs_tube4)):
#         for atom in NZs:
#             lastCA = CAs[-layer*p.N_res] if layer != 0 else atom+1
#             if (CAs[layer*p.N_res] < atom) and (atom < lastCA):
#                 NZs_LYN = np.append(NZs_LYN, atom)
    
#     return NZs_LYS, NZs_LYN
# #end


def plot_CAs(traj, step):
    frame = traj[step]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 30))

    CAs = traj.top.select("@CA")
    for atom in CAs:
        xyz = frame[atom]
        ax.scatter(xyz[0], xyz[1], xyz[2], c='b')
#end


def plot_network_3D(p, traj, step, label, canal=False, ifpath=True, layer=0, KY=False, xtalcenter=None, TFA=False):
    # Aguas y cloros del canal

    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)

    # Puentes de hidrógeno

    hbonds = pd.read_csv(label+"_hbonds.csv")
    if ifpath:
        paths = pd.read_csv(label+"_paths.csv")
        paths['path'] = paths['path'].apply(lambda x: np.array([int(i) for i in x[1:-1].split(",")]))

    ###
    
    frame = traj[step]
    
    if xtalcenter is None:
        xtal = False
    else:
        xtal = True
    if xtal:
        lvs = np.array([frame.box[0], frame.box[1], frame.box[2]])

    WATs = iWATs[step]
    CLs = iCLs[step]
    if ifpath:
        path = paths['path'].iloc[step]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 30))

    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        
        for atom in NZs_LYS:
            xyz = frame[atom]
            if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
            ax.scatter(xyz[0], xyz[1], xyz[2], c='b')

        for atom in NZs_LYN:
            xyz = frame[atom]
            if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
            ax.scatter(xyz[0], xyz[1], xyz[2], c='y')
            
        if KY:
            TYDs = select_atoms(p, layer, KY)
            for atom in TYDs:
                xyz = frame[atom]
                if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
                ax.scatter(xyz[0], xyz[1], xyz[2], c='m')

    for atom in WATs:
        xyz = frame[atom]
        if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
        ax.scatter(xyz[0], xyz[1], xyz[2], c='r')

    for atom in CLs:
        xyz = frame[atom]
        if xtal: xyz = wrap_coordinates(xyz - xtalcenter*lvs, lvs)
        ax.scatter(xyz[0], xyz[1], xyz[2], c='k')
        
    if TFA:
        resids = []
        for atom in CLs:
            resids.append(traj.top[atom].resid)
        coms = []
        for resid in set(resids):
            indices = [index for index, element in enumerate(resids) if element == resid]
            com = np.array([0.0,0.0,0.0])
            for index in indices:
                com += frame[CLs[index]]
            coms.append(com/len(indices))
        for atom in CLs:
            xyz1 = frame[atom]
            if xtal: xyz1 = wrap_coordinates(xyz1 - xtalcenter*lvs, lvs)
            xyz2 = coms[list(set(resids)).index(traj.top[atom].resid)]
            if xtal: xyz2 = wrap_coordinates(xyz2 - xtalcenter*lvs, lvs)
            if xtal and (np.abs(xyz2[2] - xyz1[2]) > lvs[2]/2):
                case = np.sign(xyz2[2] - xyz1[2])
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2]+case*lvs[2], xyz2[2]], color='m', linewidth=1)
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]-case*lvs[2]], color='m', linewidth=1)
            else:
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], color='m', linewidth=1)

    for index, bond in hbonds[hbonds["istep"] == step].iterrows():
        xyz1 = frame[bond["donor"]]
        if xtal: xyz1 = wrap_coordinates(xyz1 - xtalcenter*lvs, lvs)
        xyz2 = frame[bond["acceptor"]]
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
            if xtal and (np.sign(path[i]) != np.sign(path[i+1])):
                xyz2 = frame[abs(path[i+1])]
                xyz2 = wrap_coordinates(xyz2 - xtalcenter*lvs, lvs)
                case = np.sign(xyz2[2] - xyz1[2])
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2]+case*lvs[2], xyz2[2]], 'k--', linewidth=2.5)
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]-case*lvs[2]], 'k--', linewidth=2.5)
            else:
                xyz2 = frame[abs(path[i+1])]
                if xtal: xyz2 = wrap_coordinates(xyz2 - xtalcenter*lvs, lvs)
                ax.plot3D([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], 'k--', linewidth=2.5)
#end


def plot_network_2D(p, traj, step, label, canal=False, phimin=-4, ifpath=True, tube=0, layer=0, KY=False):
    # Aguas y cloros del canal

    iWATs = np.load("iWATs_"+label+".npy", allow_pickle=True)
    iCLs = np.load("iCLs_"+label+".npy", allow_pickle=True)

    # Puentes de hidrógeno

    hbonds = pd.read_csv(label+"_hbonds.csv")
    if ifpath:
        paths = pd.read_csv(label+"_paths.csv")
        paths['path'] = paths['path'].apply(lambda x: np.array([int(i) for i in x[1:-1].split(",")]))

    ###
    
    frame = traj[step]

    WATs = iWATs[step]
    CLs = iCLs[step]
    if ifpath:
        path = paths['path'].iloc[step]

    fig = plt.figure()
    ax = fig.add_subplot()

    if canal:
        NZs_LYS, NZs_LYN = select_atoms(p, layer)
        
        for atom in NZs_LYS:
            xyz = frame[atom]
            phi = np.arctan2(xyz[1], xyz[0])
            if phi < phimin:
                phi += 2*np.pi
            ax.scatter(phi, xyz[2], c='b')

        for atom in NZs_LYN:
            xyz = frame[atom]
            phi = np.arctan2(xyz[1], xyz[0])
            if phi < phimin:
                phi += 2*np.pi
            ax.scatter(phi, xyz[2], c='y')
            
        if KY:
            TYDs = select_atoms(p, layer, KY)
            for atom in TYDs:
                xyz = frame[atom]
                phi = np.arctan2(xyz[1], xyz[0])
                if phi < phimin:
                    phi += 2*np.pi
    
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

    for atom in CLs:
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


def plot_average_3D(label, frame_cutoff):
    positions = pd.read_csv(label+"_averages.csv")
    positions['ave_pos'] = positions['ave_pos'].apply(lambda x: np.array([float(i) for i in x[1:-1].split()]))
    positions['var_pos'] = positions['var_pos'].apply(lambda x: np.array([float(i) for i in x[1:-1].split()]))

    ###

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 30))

    for atom in positions['atom'].unique():
        if positions[positions['atom'] == atom]['count'].iloc[0] >= frame_cutoff:
            xyz = positions[positions['atom'] == atom]['ave_pos'].iloc[0]
            var = positions[positions['atom'] == atom]['var_pos'].iloc[0]
            if positions[positions['atom'] == atom]['name'].iloc[0] == "NZ_LYS": color = 'b'
            elif positions[positions['atom'] == atom]['name'].iloc[0] == "NZ_LYN": color = 'y'
            elif positions[positions['atom'] == atom]['name'].iloc[0] == "CL": color = 'k'
            elif positions[positions['atom'] == atom]['name'].iloc[0] == "WAT": color = 'r'
            ax.scatter(xyz[0], xyz[1], xyz[2], c=color, s=500*np.linalg.norm(var))
#end


def plot_average_2D(label, frame_cutoff, phimin=-4):
    positions = pd.read_csv(label+"_averages.csv")
    positions['ave_pos'] = positions['ave_pos'].apply(lambda x: np.array([float(i) for i in x[1:-1].split()]))
    positions['var_pos'] = positions['var_pos'].apply(lambda x: np.array([float(i) for i in x[1:-1].split()]))

    ###

    fig = plt.figure()
    ax = fig.add_subplot()

    for atom in positions['atom'].unique():
        if positions[positions['atom'] == atom]['count'].iloc[0] >= frame_cutoff:
            xyz = positions[positions['atom'] == atom]['ave_pos'].iloc[0]
            var = positions[positions['atom'] == atom]['var_pos'].iloc[0]
            phi = np.arctan2(xyz[1], xyz[0])
            if phi < phimin:
                phi += 2*np.pi
            if positions[positions['atom'] == atom]['name'].iloc[0] == "NZ_LYS": color = 'b'
            elif positions[positions['atom'] == atom]['name'].iloc[0] == "NZ_LYN": color = 'y'
            elif positions[positions['atom'] == atom]['name'].iloc[0] == "CL": color = 'k'
            elif positions[positions['atom'] == atom]['name'].iloc[0] == "WAT": color = 'r'
            ax.scatter(phi, xyz[2], c=color, s=500*np.linalg.norm(var))


    ax.set_xlabel("Polar angle $\Phi$")
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


def probability_cloud(xyz, layers, color, scale=5):
    xs = xyz[:,0]
    ys = xyz[:,1]
    zs = xyz[:,2]
    phis = np.arctan2(ys, xs)
    
    phi_scale = round(max(phis) - min(phis))*scale
    z_scale = round(max(zs) - min(zs))*scale
    bins = [phi_scale, z_scale]
    
    fig = plt.figure(figsize = (8, 12), constrained_layout=True)
    gs = GridSpec(layers, 2, figure=fig)
    ax0 = fig.add_subplot(gs[:, 0])
    cmap = mpl.colors.LinearSegmentedColormap.from_list('', ['white', color])
    h = ax0.hist2d(phis, zs, bins=bins, cmap=cmap)#, norm=mpl.colors.LogNorm())
    ax0.set_aspect('equal')
    ax0.set_xlabel('$\Phi$')
    ax0.set_ylabel('z ($\AA$)')
    # cb = fig.colorbar(h[3], ax=ax0)
    for i in range(1, layers):
        ax0.plot([min(phis), max(phis)], [i*(max(zs) - min(zs))/layers + min(zs), i*(max(zs) - min(zs))/layers + min(zs)], 'k--')

    h, edges = np.histogramdd(xyz, np.array([20, 20, layers]))
    for i in range(0, layers):
        _ax = fig.add_subplot(gs[layers-1-i, -1])
        _ax.set_aspect('equal')
        _ax.set_xlabel('x ($\AA$)')
        _ax.set_ylabel('y ($\AA$)')
        _ax.imshow(np.flipud(np.transpose(h[:, :, i])), cmap=cmap, extent=[edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
        
    return fig
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



