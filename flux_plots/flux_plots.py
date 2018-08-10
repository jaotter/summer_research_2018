from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import radio_beam
import astropy.units as u
from deconv_sources import deconv_srcs
from functools import reduce
import fnmatch


#names = ['B3', 'B6', 'B7_hr', 'B7_lr']
#imgs = ['/lustre/aoc/students/jotter/directory/Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/Orion_SourceI_B7_continuum_r-2.mask5mJy.clean4mJy.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/member.uid___A001_X88e_X1dd.Orion_BNKL_source_I_sci.spw25_27_29_31.cont.I.pbcor.fits']

def get_ind(names): #returns indices where sources are detected in all bands with good gaussian fits
    data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits')
    flux_inds = [np.where(np.isnan(data['ap_flux_'+n]) == False)[0] for n in names]
    fit_inds = [np.where(data['fit_goodness_'+n] == 'y')[0] for n in names]
    ind1 = reduce(np.intersect1d, flux_inds)
    ind2 = reduce(np.intersect1d, fit_inds)
    ind = np.intersect1d(ind1, ind2)
    return ind

def plot_alpha(names, imgs, only_deconv=False, flux_type='aperture'): 
    #plot flux-flux alpha scatterplot, and alpha histogram
    #accepts 2 names/imgs, first will be plotted on x axis, second on y
    #if only_deconv = True, then only plot sources which can be deconvolved
    #flux_type = 'aperture' uses aperture fluxes, flux_type = 'gaussian' uses corrected gaussian amplitude

    if len(names) != 2:
        raise ValueError('Must have exactly two datasets for alpha plots')
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7_hr':339.7672758867*u.GHz, 'B7_lr':339.7672758867*u.GHz}  #first two from Adam, third from member.uid header
    data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits')
    
    
    if only_deconv == True:
        ind = deconv_srcs(names,imgs)
        lab = 'deconvolved sources'
    else:
        flux_ind = np.intersect1d(np.where(np.isnan(data['ap_flux_'+names[0]])== False)[0],np.where(np.isnan(data['ap_flux_'+names[1]])==False)[0]) #sources with measured aperture fluxes
        fit_ind = np.intersect1d(np.where(data['fit_goodness_'+names[0]] =='y'), np.where(data['fit_goodness_'+names[1]] =='y')) #sources w good gaussian fits in both bands
        ind = np.intersect1d(flux_ind, fit_ind)
        lab = 'all sources'
        
    if flux_type == 'gaussian':
        amp_corr = [[],[]] #corrected gaussian amplitude
        amp_err_corr = [[],[]]
        for i in range(len(names)):
            fl = fits.open(imgs[i])
            beam = radio_beam.Beam.from_fits_header(fl[0].header)
            corr_fac = (data['FWHM_major_'+names[i]]*data['FWHM_minor_'+names[i]])/(beam.major.to(u.arcsec)*beam.minor.to(u.arcsec))[ind] #corrective factor to amplitude: FWHM maj*min in arcsec / beam maj*min in arcsec
            amp_corr[i].append(data['gauss_amplitude_'+names[i]][ind]*corr_fac)
            amp_err_corr[i].append(data['amplitude_err_'+names[i]][ind]*corr_fac)
        plot1 = amp_corr[0]
        plot2 = amp_corr[1]
        plot1_err = amp_err_corr[0]
        plot2_err = amp_err_corr[1]
        fl_type = 'corrected gaussian amplitude'

    elif flux_type == 'aperture':
        plot1 = data['ap_flux_'+names[0]][ind]
        plot2 = data['ap_flux_'+names[1]][ind]
        plot1_err = data['ap_flux_err_'+names[0]][ind]
        plot2_err = data['ap_flux_err_'+names[1]][ind]
        fl_type = 'aperture flux'

    freq1 = freqs[names[0]]
    freq2 = freqs[names[1]]

    

    fig1 = plt.figure() #figure 1: flux-flux plot
    
    plt.xlabel('{name} {flux}'.format(name=names[0], flux = fl_type))
    plt.ylabel('{name} {flux}'.format(name=names[1], flux = fl_type))

    plt.errorbar(plot1, plot2, xerr=plot1_err, yerr=plot2_err, linestyle='', marker='*',color='blue',label=lab)
    
    F1 = np.linspace(np.min(plot1), np.max(plot1), 10) #plot lines of different alpha values
    alphas = np.arange(0.5, 3, 0.5)
    for a in alphas:
        alpha_F2 = F1*((freq2/freq1)**a)
        plt.loglog(F1, alpha_F2, linestyle=':', label='alpha={a}'.format(a=a))
    plt.legend()

    plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/'+names[0]+names[1]+'_fluxflux.png')


    fig2 = plt.figure() #figure 2: alpha histogram

    alpha_gauss = np.log((plot1/plot2))/np.log(freq1/freq2)
    plt.hist(alpha_gauss)

    plt.xlabel('alpha')
    plt.ylabel('number of sources')
    plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/'+names[0]+names[1]+'_alpha_hist.png')

def colorcolor_plot(x_name, y_name, num_name): 
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7_hr':339.7672758867*u.GHz, 'B7_lr':339.7672758867*u.GHz}  #first two from Adam, third from member.uid header
    data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits')

    flux_ind = reduce(np.intersect1d, (np.where(np.isnan(data['ap_flux_'+x_name])== False)[0],np.where(np.isnan(data['ap_flux_'+y_name])==False)[0],np.where(np.isnan(data['ap_flux_'+num_name])==False)[0])) #sources with measured aperture fluxes
    fit_ind = reduce(np.intersect1d, (np.where(data['fit_goodness_'+x_name]=='y')[0],np.where(data['fit_goodness_'+y_name]=='y')[0],np.where(data['fit_goodness_'+num_name]=='y')[0])) #sources w good gaussian fits in both bands
    ind = np.intersect1d(flux_ind, fit_ind)
    

    plot1 = data['ap_flux_'+x_name][ind]
    plot2 = data['ap_flux_'+y_name][ind]
    plot1_err = data['ap_flux_err_'+x_name][ind]
    plot2_err = data['ap_flux_err_'+y_name][ind]
    plot_num = data['ap_flux_'+num_name][ind]
    plot_num_err = data['ap_flux_err_'+num_name][ind]

    plotx = plot_num/plot1
    ploty = plot_num/plot2
    plotx_err = plot1*np.sqrt((plot1_err/plot1)**2 + (plot_num_err/plot_num)**2) #error is fractional error added in quad
    ploty_err = plot1*np.sqrt((plot2_err/plot2)**2 + (plot_num_err/plot_num)**2)

    plt.errorbar(plotx, ploty, xerr=plotx_err, yerr=ploty_err, linestyle='', marker='*')

    plt.xlabel('{num}/{x}'.format(num=num_name,x=x_name))
    plt.ylabel('{num}/{y}'.format(num=num_name,y=y_name))

    totmin = np.amin((plotx, ploty))
    totmax = np.amax((plotx, ploty))
    plt.xlim(totmin-0.05, totmax+0.05)
    plt.ylim(totmin-0.05, totmax+0.05)
    
    #plot alphas:
    alpha = [1.5,2,2.5]
    cols = ['r', 'g', 'k']
    for i,a in enumerate(alpha):
        plt.axvline((freqs[num_name]/freqs[x_name])**a, linestyle='--', label=r'$\alpha={a}$'.format(a=a), color=cols[i])
        plt.axhline((freqs[num_name]/freqs[y_name])**a, linestyle='--', label=r'$\alpha={a}$'.format(a=a), color=cols[i])

    plt.legend()
    plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/'+num_name+'div'+x_name+y_name+'_colcol.png', dpi=300)


def plot_SEDs(names):
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7_hr':339.7672758867*u.GHz, 'B7_lr':339.7672758867*u.GHz}
    data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog_B7.fits')
    ind = get_ind(names)

    alpha = [1.5,2,2.5]
    freq_x = np.array([freqs[n].value for n in names])
    for i in ind:
        fluxes = [data['ap_flux_'+n][i] for n in names]
        fluxes = np.array(fluxes) - np.array([data['bg_median_'+n][i] for n in names])
        flux_err = [data['ap_flux_err_'+n][i] for n in names]

        F2 = fluxes[1]
        nu2 = freq_x[1]
        
        plt.figure()
        plt.errorbar(freq_x, fluxes, yerr=flux_err, linestyle='', marker='o')
        for j in range(len(alpha)):
            F1 = F2*((freq_x/nu2)**alpha[j])
            plt.plot(freq_x, F1, linestyle='-', label=r'$\alpha = $'+str(alpha[j]))
		
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('aperture flux (mJy)')
        plt.xlabel('frequency (GHz)')
        plt.legend()
        plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/SEDs/SED_'+str(i)+'_'+names[-1]+'_sub_medbg.png', dpi=300)


def image_fluxes(srcID):
    ref_names = ['B3', 'B6', 'B7_hr']
    
    for ref_name in ref_names:
        if ref_name == 'B7_hr':
            data = Table.read('/lustre/aoc/students/jotter/dendro_catalogs/'+ref_name+'_diff_imgs_catalog.fits')
        else:
            data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/'+ref_name+'_diff_imgs_catalog.txt')
        srcind = np.where(data['D_ID'] == srcID)
        flux_names = fnmatch.filter(data.colnames, 'ap_flux_r*')
        bg_names = fnmatch.filter(data.colnames, 'bg_median_r*')
        sub_fluxes = []
        for n in range(len(flux_names)):
            sub_fluxes.append((data[flux_names[n]][srcind] - data[bg_names[n]][srcind], flux_names[n][8:]))
        sub_fluxes.sort()
        
        x = np.arange(len(sub_fluxes))

        plt.figure(figsize=(10,10))
        plt.semilogy(x, [sf[0] for sf in sub_fluxes], marker='o')
        plt.xticks(x, [sf[1] for sf in sub_fluxes], rotation='vertical')
        plt.subplots_adjust(bottom=0.5) 
        plt.grid()
    plt.show()


