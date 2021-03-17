from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, upper_lim_flags, savepath, left_censor=True, cdf=False, plot_quantity='Mdust'):
    kmf = KaplanMeierFitter()

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    colors = ['tab:red','tab:blue','tab:green', 'tab:orange', 'tab:purple', 'gray', 'brown']
    for ind in range(len(arrays)):
        col = colors[ind]
        print(labels[ind])
        if upper_lim_flags[ind] is not None:
            if left_censor == True:
                kmf.fit_left_censoring(arrays[ind], upper_lim_flags[ind], label=labels[ind])
            else:
                kmf.fit(arrays[ind], event_observed=upper_lim_flags[ind], label=labels[ind]) #right censoring                                                                                     
        else:
            kmf.fit(arrays[ind], upper_lim_flags[ind], label=labels[ind])

        #kmf.confidence_interval_survival_function_.plot(ax=ax)                                                                                                                               
        #kmf.survival_function_.plot(ax=ax)                                                                                                                                                   
        if cdf == True:
            kmf.plot(ax=ax)
        else:
            size = np.array(kmf.survival_function_[labels[ind]].axes).flatten()[1:]
            prob = kmf.survival_function_[labels[ind]].values[1:]

            lower = np.array(kmf.confidence_interval_survival_function_[f'{labels[ind]}_lower_0.95'].values[1:])
            upper = np.array(kmf.confidence_interval_survival_function_[f'{labels[ind]}_upper_0.95'].values[1:])

            ax.plot(size, prob, label=labels[ind], color=col)
            ax.fill_between(size, lower, upper, color=col, alpha=0.25)

    plt.legend(fontsize=16)
    
    ax.set_xscale('log')
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    
    if plot_quantity == 'Mdust':
        ax.set_xlim(0.03, 1000)   
        ax.set_xlabel(r'$M_{dust}$ ($M_{\oplus}$)', fontsize=18)
        ax.set_ylabel(r'$P \geq M_{dust}$', fontsize=18)

    if plot_quantity == 'Rdisk':
        ax.set_xlim(3, 200)
        ax.set_xlabel(r'$R_{disk}$ (AU)', fontsize=18)
        ax.set_ylabel(r'$P \geq R_{disk}$', fontsize=18)
    
    plt.savefig(savepath)
    print(f'saved figure at {savepath}')


def KM_median(array, upper_lim_flags, left_censor=True, percentile=True):
    kmf = KaplanMeierFitter()

    if upper_lim_flags is not None:
        if left_censor == True:
            kmf.fit_left_censoring(array, upper_lim_flags)
        else:
            kmf.fit(array, event_observed=upper_lim_flags) #right censoring                                                                                    
    else:
        kmf.fit(array, upper_lim_flags)

    median = median_survival_times(kmf.survival_function_)

    median_ci = median_survival_times(kmf.confidence_interval_).values 
    
    upper_perc = kmf.percentile(0.25)
    lower_perc = kmf.percentile(0.75)

    if percentile == True:
        print(f'median and 1st/3rd quartiles: {median}, {lower_perc}, {upper_perc}')
        return median, upper_perc, lower_perc

    else:
        print(f'median and CI: {median}, {median_ci}')
        return median, median_ci[0][0], median_ci[0][1]
