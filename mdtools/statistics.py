### IMPORTS ###


from mdtools.core import *


### STATISTICS ###


def N_hist(data, x=None, bins='auto', hue=None, element='bars', alpha=1, title=None, xlabel=None):
    N_data = data[x]
    if bins=='int':
        bins = np.arange(min(N_data)-0.5, max(N_data)+1, 1)
    
    fig, ax = plt.subplots()
    sns.histplot(data=data, x=x, hue=hue, bins=bins, lw=2, element=element, alpha=alpha)
    ylabel = 'Count (number of frames)'
    decorate_ax(ax, title, 16, xlabel, ylabel, 14, 12, 2, 4, False)
#end


def evolution(data, x=None, y=None, hue=None, label=None, title=None, ylabel=None):
    fig, ax = plt.subplots()
    sns.lineplot(data=data, x=x, y=y, hue=hue)
    xlabel = 'Step'
    decorate_ax(ax, title, 16, xlabel, ylabel, 14, 12, 2, 4, False)
#end


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
