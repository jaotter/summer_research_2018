from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from scipy import stats
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_KM(arrays, labels, upper_lim_flags, savepath, left_censor=True, cdf=False, plot_quantity='Mdust', noerr_inds=[]):
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
            if ind > 0:
                prev_prob = prob
                prev_size = size
            size = np.array(kmf.survival_function_[labels[ind]].axes).flatten()[1:]
            prob = kmf.survival_function_[labels[ind]].values[1:]

            lower = np.array(kmf.confidence_interval_survival_function_[f'{labels[ind]}_lower_0.95'].values[1:])
            upper = np.array(kmf.confidence_interval_survival_function_[f'{labels[ind]}_upper_0.95'].values[1:])

            ax.plot(size, prob, label=labels[ind], color=col)
            if ind not in noerr_inds:
                ax.fill_between(size, lower, upper, color=col, alpha=0.25)

            elif ind == noerr_inds[1]:
                interp_new = interp1d(size, prob)
                interp_prob = interp_new(prev_size)
                ax.fill_between(prev_size, prev_prob, interp_prob, color=col, alpha=0.25)
                #ax.fill(np.concatenate((prev_size,size)), np.concatenate((prev_prob,prob)), color=col, alpha=0.25)
                
    plt.legend(fontsize=16)
    
    ax.set_xscale('log')
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    
    if plot_quantity == 'Mdust':
        ax.set_xlim(0.03, 2000)   
        ax.set_xlabel(r'$M_{dust}$ ($M_{\oplus}$)', fontsize=18)
        ax.set_ylabel(r'$P \geq M_{dust}$', fontsize=18)

    if plot_quantity == 'Rdisk':
        ax.set_xlim(3, 200)
        ax.set_xlabel(r'$R_{disk}$ (AU)', fontsize=18)
        ax.set_ylabel(r'$P \geq R_{disk}$', fontsize=18)
    
    plt.savefig(savepath)
    print(f'saved figure at {savepath}')


def KM_median(array, upper_lim_flags, left_censor=True, return_type='percentile'):
    kmf = KaplanMeierFitter()

    if upper_lim_flags is not None:
        if left_censor == True:
            kmf.fit_left_censoring(array, upper_lim_flags)
        else:
            kmf.fit(array, event_observed=upper_lim_flags) #right censoring                                                                                    
    else:
        kmf.fit(array, upper_lim_flags)

    median = median_survival_times(kmf.survival_function_)


    if return_type == 'percentile':
        upper_perc = kmf.percentile(0.25)
        lower_perc = kmf.percentile(0.75)
    
        print(f'median and 1st/3rd quartiles: {median}, {lower_perc}, {upper_perc}')
        return median, upper_perc, lower_perc

    elif return_type == 'ci':
        median_ci = median_survival_times(kmf.confidence_interval_).values 
        print(f'median and CI: {median}, {median_ci}')
        return median, median_ci[0][0], median_ci[0][1]

    elif return_type == 'median':
        return median


def bootstrap_ci(num_samples, data, ulim_flag, alpha=0.99):

    bootstrap_vals = []
    for i in range(num_samples):
        resample_ind = np.random.randint(0,high=len(data),size=len(data))
        resamp_data = data[resample_ind]
        resamp_flag = ulim_flag[resample_ind]

        median = KM_median(resamp_data, resamp_flag, left_censor=True, return_type='median')
        if median == 0:
            print('zero median')
            print(len(np.where(resamp_flag == True)[0]))
        bootstrap_vals.append(median)

    se = np.std(bootstrap_vals)
    Zval = stats.norm.ppf(alpha)

    CI_range = se*Zval
    print(f'CI range: {CI_range}')
    
    return CI_range
