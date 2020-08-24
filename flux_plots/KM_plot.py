from lifelines import KaplanMeierFitter
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, upper_lim_flags, savepath):
    kmf = KaplanMeierFitter()

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    colors = ['tab:red','tab:blue','tab:green', 'tab:orange', 'yellow', 'tab:purple']
    for ind in range(len(arrays)):
        col = colors[ind]
        print(labels[ind])
        if upper_lim_flags[ind] is not None:
            #kmf.fit_left_censoring(arrays[ind], upper_lim_flags[ind], label=labels[ind])                                                                                                     
            kmf.fit(arrays[ind], event_observed=upper_lim_flags[ind], label=labels[ind]) #right censoring                                                                                     
        else:
            kmf.fit(arrays[ind], upper_lim_flags[ind], label=labels[ind])

        #kmf.confidence_interval_survival_function_.plot(ax=ax)                                                                                                                               
        #kmf.survival_function_.plot(ax=ax)                                                                                                                                                   

        size = np.array(kmf.survival_function_[labels[ind]].axes).flatten()[1:]
        prob = kmf.survival_function_[labels[ind]].values[1:]

        lower = np.array(kmf.confidence_interval_survival_function_[f'{labels[ind]}_lower_0.95'].values[1:])
        upper = np.array(kmf.confidence_interval_survival_function_[f'{labels[ind]}_upper_0.95'].values[1:])

        ax.plot(size, prob, label=labels[ind], color=col)
        ax.fill_between(size, lower, upper, color=col, alpha=0.25)

    plt.legend()
    ax.set_xlabel(r'FWHM (AU)')
    ax.set_ylabel(r'$P(R \geq$ FWHM)')
    ax.set_xscale('log')
    ax.set_xlim(10, 150)
    plt.savefig(savepath)
    print(f'saved figure at {savepath}')
