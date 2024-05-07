### IMPORTS ###


from mdtools.core import *


### STATISTICS ###


def stability_hist(df, nbins=100, color='b', label=None, fig=None, ax=None):
    data = []
    for index, row in df.iterrows():
        data += json.loads(row['intervals'])

    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    ax.hist(data, nbins, color=color, histtype='step')
    h = ax.hist(data, nbins, color=color, alpha=0.5, label=label)
    ax.plot([np.array(data).mean(), np.array(data).mean()], [0, np.max(h[0])], color=color, ls='--', lw=2)
    return fig, ax, data, h
#end


### EOF ###
