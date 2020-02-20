from astropy.table import Table
from lifelines import KaplanMeierFitter

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, savepath='/home/jotter/nrao/summer_research_2018/plots/KM_comparisons.png'):

    kmf = KaplanMeierFitter()

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    for ind in range(len(arrays)):
        kmf.fit(arrays[ind], label=labels[ind])
        kmf.plot(ax=ax)

    plt.legend()
    ax.set_xlabel(r'$\log(M_{dust}/M_{\oplus})$')
    ax.set_ylabel(r'$P \geq M_{dust}$')
    plt.savefig(savepath)


