### IMPORTS ###


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
import pandas as pd
import numpy as np
import pytraj as pt
from scipy.stats import sem
from icosphere import icosphere
from tabulate import tabulate

from main import *


### ESTADÍSTICAS ###


def plot_coordination_atom(label, atom):
    cn = pd.read_csv(label+"_coordination.csv")
    
    Nbins = np.linspace(-0.5, 10.5, 12)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    aux = cn[cn['atom'] == atom]
    print(aux['name'].iloc[0])

    cn_cl = aux['cn_cl'].iloc[:].values
    cn_wat = aux['cn_wat'].iloc[:].values
    cn_lys = aux['cn_lys'].iloc[:].values
    cn_lyn = aux['cn_lyn'].iloc[:].values
    cn_tot = cn_cl + cn_wat + cn_lys + cn_lyn

    hist, bin_edges = np.histogram(cn_cl, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'k', where='mid', label="CL", lw=2)
    hist, bin_edges = np.histogram(cn_wat, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'r', where='mid', label="WAT(O)", lw=2)
    hist, bin_edges = np.histogram(cn_lys, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'b', where='mid', label="LYS(NZ)", lw=2)
    hist, bin_edges = np.histogram(cn_lyn, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'y', where='mid', label="LYN(NZ)", lw=2)
    hist, bin_edges = np.histogram(cn_tot, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'k--', where='mid', label="TOT", lw=2)

    ax.set_xticks(range(11))
    ax.legend()
    return fig, ax
#end


def plot_coordinations(label):
    cn = pd.read_csv(label+"_coordination.csv")
    Nbins = np.linspace(-0.5, 10.5, 12)
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))

    for index, name in enumerate(["CL", "WAT", "LYS", "LYN"]):
        aux = cn[cn['name'] == name]
        cn_cl = aux['cn_cl'].iloc[:].values
        cn_wat = aux['cn_wat'].iloc[:].values
        cn_lys = aux['cn_lys'].iloc[:].values
        cn_lyn = aux['cn_lyn'].iloc[:].values
        cn_tot = cn_cl + cn_wat + cn_lys + cn_lyn

        i = int(np.floor(index/2))
        j = index%2

        hist, bin_edges = np.histogram(cn_cl, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'k', where='mid', label="CL", lw=2)
        hist, bin_edges = np.histogram(cn_wat, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'r', where='mid', label="WAT(O)", lw=2)
        hist, bin_edges = np.histogram(cn_lys, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'b', where='mid', label="LYS(NZ)", lw=2)
        hist, bin_edges = np.histogram(cn_lyn, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'y', where='mid', label="LYN(NZ)", lw=2)
        hist, bin_edges = np.histogram(cn_tot, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'k--', where='mid', label="TOT", lw=2)
        axs[i, j].set_title(name)
        axs[i, j].set_xticks(range(11))
        axs[i, j].legend()
        
    return fig, axs
#end


def plot_coordination_averages(label):
    cn = pd.read_csv(label+"_coordination.csv")
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))

    for index, name in enumerate(["CL", "WAT", "LYS", "LYN"]):
        aux = cn[cn['name'] == name]
        cn_cl = aux['cn_cl'].iloc[:].values
        cn_wat = aux['cn_wat'].iloc[:].values
        cn_lys = aux['cn_lys'].iloc[:].values
        cn_lyn = aux['cn_lyn'].iloc[:].values
        cn_tot = cn_cl + cn_wat + cn_lys + cn_lyn

        i = int(np.floor(index/2))
        j = index%2

        xs = range(5)
        means = [cn_cl.mean(), cn_wat.mean(), cn_lys.mean(), cn_lyn.mean(), cn_tot.mean()]
        errors = [cn_cl.std(), cn_wat.std(), cn_lys.std(), cn_lyn.std(), cn_tot.std()]
        colors = ['k', 'r', 'b', 'y', 'gray']
        axs[i, j].scatter(xs, means, c=colors)
        for bar in range(5):
            axs[i, j].errorbar(xs[bar], means[bar], yerr=errors[bar], fmt="o", capsize=4, capthick=2, color=colors[bar])
        axs[i, j].set_xticks(range(5), ["CL", "WAT", "LYS", "LYN", "TOT"])
        axs[i, j].set_title(name)
        
    return fig, axs
#end


def plot_coordinationbonds_atom(label, atom):
    cn = pd.read_csv(label+"_coordination_hbonds.csv")
    
    Nbins = np.linspace(-0.5, 8.5, 10)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    aux = cn[cn['atom'] == atom]
    print(aux['name'].iloc[0])
    
    hb_d_cl = aux['hb_d_cl'].iloc[:].values
    hb_d_lyn = aux['hb_d_lyn'].iloc[:].values
    hb_d_wat = aux['hb_d_wat'].iloc[:].values
    hb_a_lys = aux['hb_a_lys'].iloc[:].values
    hb_a_lyn = aux['hb_a_lyn'].iloc[:].values
    hb_a_wat = aux['hb_a_wat'].iloc[:].values
    hb_d_tot = hb_d_cl + hb_d_wat + hb_d_lyn
    hb_a_tot = hb_a_wat + hb_a_lys + hb_a_lyn
    hb_tot = hb_d_tot + hb_a_tot

    hist, bin_edges = np.histogram(hb_d_cl, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'k', where='mid', label="donor with CL", lw=2)
    hist, bin_edges = np.histogram(hb_d_lyn, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'y', where='mid', label="donor with LYN", lw=2)
    hist, bin_edges = np.histogram(hb_d_wat, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'r', where='mid', label="donor with WAT", lw=2)
    hist, bin_edges = np.histogram(hb_a_lys, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'b--', where='mid', label="acceptor with LYS", lw=2)
    hist, bin_edges = np.histogram(hb_a_lyn, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'y--', where='mid', label="acceptor with LYN", lw=2)
    hist, bin_edges = np.histogram(hb_a_wat, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'r--', where='mid', label="acceptor with WAT", lw=2)
    hist, bin_edges = np.histogram(hb_d_tot, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', where='mid', label="TOTAL as donor", lw=2)
    hist, bin_edges = np.histogram(hb_a_tot, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', linestyle='--', where='mid', label="TOTAL as acceptor", lw=2)
    hist, bin_edges = np.histogram(hb_tot, bins=Nbins)
    ax.step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', linestyle=':', where='mid', label="TOTAL", lw=2)

    ax.set_xticks(range(9))
    ax.legend()
    return fig, ax
#end


def plot_coordinationsbonds(label):
    cn = pd.read_csv(label+"_coordination_hbonds.csv")
    Nbins = np.linspace(-0.5, 8.5, 10)
    fig, axs = plt.subplots(4, 2, figsize=(14, 18))

    for index, name in enumerate(["CL", "WAT", "LYS", "LYN"]):
        aux = cn[cn['name'] == name]
        hb_d_cl = aux['hb_d_cl'].iloc[:].values
        hb_d_lyn = aux['hb_d_lyn'].iloc[:].values
        hb_d_wat = aux['hb_d_wat'].iloc[:].values
        hb_a_lys = aux['hb_a_lys'].iloc[:].values
        hb_a_lyn = aux['hb_a_lyn'].iloc[:].values
        hb_a_wat = aux['hb_a_wat'].iloc[:].values
        hb_d_tot = hb_d_cl + hb_d_wat + hb_d_lyn
        hb_a_tot = hb_a_wat + hb_a_lys + hb_a_lyn
        hb_tot = hb_d_tot + hb_a_tot

        i = 2*int(np.floor(index/2))
        j = index%2

        hist, bin_edges = np.histogram(hb_d_cl, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'k', where='mid', label="donor with CL", lw=2)
        hist, bin_edges = np.histogram(hb_d_lyn, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'y', where='mid', label="donor with LYN", lw=2)
        hist, bin_edges = np.histogram(hb_d_wat, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'r', where='mid', label="donor with WAT", lw=2)
        hist, bin_edges = np.histogram(hb_a_lys, bins=Nbins)
        axs[i+1, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'b--', where='mid', label="acceptor with LYS", lw=2)
        hist, bin_edges = np.histogram(hb_a_lyn, bins=Nbins)
        axs[i+1, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'y--', where='mid', label="acceptor with LYN", lw=2)
        hist, bin_edges = np.histogram(hb_a_wat, bins=Nbins)
        axs[i+1, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'r--', where='mid', label="acceptor with WAT", lw=2)
        hist, bin_edges = np.histogram(hb_d_tot, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', where='mid', label="TOTAL as donor", lw=2)
        hist, bin_edges = np.histogram(hb_a_tot, bins=Nbins)
        axs[i+1, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', linestyle='--', where='mid', label="TOTAL as acceptor", lw=2)
        hist, bin_edges = np.histogram(hb_tot, bins=Nbins)
        axs[i, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', linestyle=':', where='mid', label="TOTAL", lw=2)
        axs[i+1, j].step((bin_edges[1:]+bin_edges[:-1])/2, hist, 'gray', linestyle=':', where='mid', label="TOTAL", lw=2)
        axs[i, j].set_title(name)
        axs[i, j].set_xticks(range(9))
        axs[i, j].legend()
        axs[i+1, j].set_xticks(range(9))
        axs[i+1, j].legend()
    
    return fig, axs
#end


def plot_coordinationbonds_averages(label):
    cn = pd.read_csv(label+"_coordination_hbonds.csv")
    fig, axs = plt.subplots(2, 2, figsize=(18, 12))

    for index, name in enumerate(["CL", "WAT", "LYS", "LYN"]):
        aux = cn[cn['name'] == name]
        hb_d_cl = aux['hb_d_cl'].iloc[:].values
        hb_d_lyn = aux['hb_d_lyn'].iloc[:].values
        hb_d_wat = aux['hb_d_wat'].iloc[:].values
        hb_a_lys = aux['hb_a_lys'].iloc[:].values
        hb_a_lyn = aux['hb_a_lyn'].iloc[:].values
        hb_a_wat = aux['hb_a_wat'].iloc[:].values
        hb_d_tot = hb_d_cl + hb_d_wat + hb_d_lyn
        hb_a_tot = hb_a_wat + hb_a_lys + hb_a_lyn
        hb_tot = hb_d_tot + hb_a_tot

        i = int(np.floor(index/2))
        j = index%2

        xs = range(9)
        means = [hb_d_cl.mean(), hb_d_lyn.mean(), hb_d_wat.mean(), hb_d_tot.mean(), hb_a_lys.mean(), hb_a_lyn.mean(), hb_a_wat.mean(), hb_a_tot.mean(), hb_tot.mean()]
        errors = [hb_d_cl.std(), hb_d_lyn.std(), hb_d_wat.std(), hb_d_tot.std(), hb_a_lys.std(), hb_a_lyn.std(), hb_a_wat.std(), hb_a_tot.std(), hb_tot.std()]
        colors = ['k', 'y', 'r', 'gray', 'b', 'y', 'r', 'gray', 'gray']
        axs[i, j].scatter(xs, means, c=colors)
        for bar in range(9):
            axs[i, j].errorbar(xs[bar], means[bar], yerr=errors[bar], fmt="o", capsize=4, capthick=2, color=colors[bar])
        axs[i, j].set_xticks(range(9), ["d w/ CL", "d w/ LYN", "d w/ WAT", "TOT as d", "a w/ LYS", "a w/ LYN", "a w/ WAT", "TOT as a", "TOT"])
        axs[i, j].set_title(name)
        
    return fig, axs
#end


def print_configurations(label, atomtype):
    cnb = pd.read_csv(label+"_coordination_hbonds.csv")
    aux = cnb[cnb["name"]==atomtype]
    new = pd.Series(dtype=int)
    for index, row in aux.iterrows():
        new.loc[index] = 100000*row['hb_d_cl'] + 10000*row['hb_d_lyn'] + 1000*row['hb_d_wat'] + \
        100*row['hb_a_lys'] + 10*row['hb_a_lyn'] + row['hb_a_wat']
        
    print("Todas las configuraciones de " + atomtype + " a lo largo de la simulación:")
    print(new.shape[0])
    print("\nConfiguraciones únicas de " + atomtype + ":")
    print(new.value_counts().shape[0])
        
    table = [['Conf count', 'Dona a CL', 'Dona a LYN', 'Dona a WAT',
          'Acepta de LYS', 'Acepta de LYN', 'Acepta de WAT']]
    for i in range(new.value_counts().shape[0]):
        table.append([new.value_counts().iloc(0)[i],
                      new.value_counts().index[i]//100000,
                      new.value_counts().index[i]//10000-new.value_counts().index[i]//100000*10,
                      new.value_counts().index[i]//1000-new.value_counts().index[i]//10000*10,
                      new.value_counts().index[i]//100-new.value_counts().index[i]//1000*10,
                      new.value_counts().index[i]//10-new.value_counts().index[i]//100*10,
                      new.value_counts().index[i]//1-new.value_counts().index[i]//10*10])
    
    print(tabulate(table, headers='firstrow'))
#end


def TSV(a, b, c, d): # puntos en 3d
    # "Tetahedron Signed Volume
    return (1.0/6.0)*np.dot(np.cross(b-a, c-a), d-a)
#end


def LTI(p, t):
    # "Line - Triangle Intersection
    if np.sign(TSV(np.array([0, 0, 0]), t[0], t[1], t[2])) != np.sign(TSV(p, t[0], t[1], t[2])):
        if (np.sign(TSV(np.array([0, 0, 0]), p, t[0], t[1])) == np.sign(TSV(np.array([0, 0, 0]), p, t[1], t[2])) and
                np.sign(TSV(np.array([0, 0, 0]), p, t[1], t[2])) == np.sign(TSV(np.array([0, 0, 0]), p, t[2], t[0]))):
            return True
    return False
#end


def compute_SphereHeatMap(traj, label, atomid, bondname, nu=5):
    vertices, faces = icosphere(nu)
    
    hbonds = pd.read_csv(label+"_hbonds.csv")
    aux = hbonds[hbonds["acceptor"] == atomid]
    aux = aux[aux["residues"] == bondname]
    
    counts = np.zeros(len(faces), dtype=int)
    for step in aux["istep"].unique():
        frame = traj[step]
        aux2 = aux[aux["istep"] == step]

        for i, row in aux2.iterrows():
            donor = row["donor"]
            relpos = frame[donor] - frame[atomid]

            for index, face in enumerate(faces):
                if LTI(100*relpos, vertices[face]):
                    counts[index] += 1
    
    return vertices, faces, counts
#end


def draw_SphereHeatMap(vertices, faces, counts):
                    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((15, 15, 15))
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-1.2, 1.2)
    cmap = plt.cm.get_cmap('YlOrBr')
    for i, face in enumerate(faces):
        tri = Poly3DCollection([vertices[face]])
        tri.set_color(cmap(counts[i]/counts.max()))
        tri.set_edgecolor('k')
        ax.add_collection3d(tri)
    norm = mpl.colors.Normalize(vmin=0, vmax=counts.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ticks=np.linspace(0, counts.max(), 10), 
    boundaries=np.arange(0, counts.max()+1, 1))
    
#end
