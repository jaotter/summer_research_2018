from scipy.stats import gaussian_kde
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import radio_beam
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from deconv_sources import deconv_srcs
import regions
from scipy.stats import ks_2samp
from functools import reduce
import fnmatch
from astropy.coordinates import Angle, SkyCoord

#names = ['B3', 'B6', 'B7_hr', 'B7_lr']
#imgs = ['/lustre/aoc/students/jotter/directory/Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/Orion_SourceI_B7_continuum_r-2.mask5mJy.clean4mJy.image.tt0.pbcor.fits', '/lustre/aoc/students/jotter/directory/member.uid___A001_X88e_X1dd.Orion_BNKL_source_I_sci.spw25_27_29_31.cont.I.pbcor.fits']

def dist_flux_plot():
    data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_nov20_ulim.fits')
    theta1c = SkyCoord('05:35:16.46375', '-05:23:22.8486', unit=(u.hourangle, u.degree))
    b3_coord = SkyCoord(ra=data['RA_B3'], dec=data['DEC_B3'], unit=(u.degree, u.degree))

    
    dist = theta1c.separation(b3_coord).to(u.arcsecond)
    dist_pc = dist.to(u.radian)*400*u.pc
    
    plt.plot(dist_pc, data['fwhm_maj_deconv_B3'], linestyle='', marker='o')
    plt.xlabel('Distance to Theta 1C (pc)')
    plt.ylabel('Disk Deconvolved FWHM (arcseconds)')
    plt.savefig('/home/jotter/nrao/plots/dist_size_plot.png')

    xlim = plt.xlim()
    
    plt.clf()

    fb = Table.read('/home/jotter/nrao/tables/Forbrich_2016.fits')
    fb_coord = SkyCoord(ra=fb['RAJ2000'], dec=fb['DEJ2000'], unit=(u.degree, u.degree))

    theta1a = SkyCoord('05:35:15.8252169327', '-05:23:14.332546552', unit=(u.hourangle, u.degree))
    theta1b = SkyCoord('05:35:16.112', '-05:23:06.89', unit=(u.hourangle, u.degree))
    theta1d = SkyCoord('05:35:17.2574309230', '-05:23:16.567215876', unit=(u.hourangle, u.degree))
    
    dist_1c = (theta1c.separation(fb_coord).to(u.radian)*400*u.pc).value
    dist_1a = (theta1a.separation(fb_coord).to(u.radian)*400*u.pc).value
    dist_1b = (theta1b.separation(fb_coord).to(u.radian)*400*u.pc).value
    dist_1d = (theta1d.separation(fb_coord).to(u.radian)*400*u.pc).value

    print(dist_1c, dist_1a)
    
    dist_fb = []
    for i in range(len(fb_coord)):
        dist_fb.append(np.min([dist_1a[i],dist_1b[i],dist_1c[i],dist_1d[i]]))
    
    plt.plot(dist_fb, np.log10(fb['Spk']), linestyle='', marker='*', color='tab:orange', label='Forbrich+2016')
    plt.xlabel('Distance to Trapezium star (pc)')
    plt.ylabel(r'$\log(S_{5GHz})$ mJy')
    #plt.xlim(xlim)
    plt.savefig('/home/jotter/nrao/plots/dist_flux_forbrich_plot_full.png')

    
dist_flux_plot()

def get_ind(names): #returns indices where sources are detected in all bands with good gaussian fits
    data = fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted.fits')
    flux_inds = [np.where(np.isnan(data['ap_flux_'+n]) == False)[0] for n in names]
    #fit_inds = [np.where(data['fit_goodness_'+n] == 'y')[0] for n in names]
    fit_ind =  np.where(data['good_fit_flag'] == True)[0]
    ind1 = reduce(np.intersect1d, flux_inds)
    #ind2 = reduce(np.intersect1d, fit_inds)
    ind = np.intersect1d(ind1, fit_ind)
    return ind

def flux_hist(sources='all', exclude_BN=False, KS_test=False):
    data = Table.read('/home/jotter/nrao/summer_research_2018/tables/r0.5_catalog_bgfit_jun20.fits')
    eis = Table.read('/home/jotter/nrao/tables/eisner_tbl.txt', format='ascii')
        
    if sources == 'ONC':
        IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_apr20_full_edit.fits')
        IR_src = IR_tab['D_ID']
        IR_ind = [np.where(data['D_ID']==d_id)[0][0] for d_id in IR_src]
        data = data[IR_ind]

    if sources == 'OMC1':
        IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_apr20_full_edit.fits')
        nonIR_src = np.setdiff1d(data['D_ID'], IR_tab['D_ID'])
        nonIR_ind = [np.where(data['D_ID']==d_id)[0][0] for d_id in nonIR_src]
        data = data[nonIR_ind]

    if KS_test == True:
        IR_tab = Table.read('/home/jotter/nrao/summer_research_2018/tables/IR_matches_MLLA_apr20_full_edit.fits')
        nonIR_src = np.setdiff1d(data['D_ID'], IR_tab['D_ID'])
        print(nonIR_src)
        nonIR_ind = [np.where(data['D_ID']==d_id)[0][0] for d_id in nonIR_src]
        IR_src = IR_tab['D_ID']
        IR_ind = [np.where(data['D_ID']==d_id)[0][0] for d_id in IR_src]

        B7flux = np.log10(data['ap_flux_B7']*1000)
        
        B7flux_IR = B7flux[IR_ind]
        B7flux_IR = B7flux_IR[np.where(np.isnan(B7flux_IR) == False)[0]]

        print(nonIR_ind)
        B7flux_nonIR = B7flux[nonIR_ind]
        #B7flux_nonIR = B7flux_nonIR[np.where(np.isnan(B7flux_nonIR) == False)[0]]

        print(B7flux_IR, B7flux_nonIR)
        
        Dval, pval = ks_2samp(B7flux_IR, B7flux_nonIR)
        print(f'{pval} p value for ks test of IR/nonIR')
    
        
    if exclude_BN == True:
        BN_id = 43
        BN_ind = np.where(data['D_ID'] == 43)[0]
        data.remove_rows([BN_ind[0]])

        
    B7flux = np.log10(data['ap_flux_B7']*1000)
    B7flux = B7flux[np.where(np.isnan(B7flux) == False)[0]]

    eisflux_str = eis['F_lambda 850mum']
    eisflux = []
    for es in eisflux_str:
        eisflux.append(np.log10(float(es.split()[0])))
    
    eisflux = np.array(eisflux)[np.where(np.isinf(eisflux) == False)[0]]
        
    allhist, abins = np.histogram(np.concatenate((B7flux, eisflux)))
    B7hist, b = np.histogram(B7flux, bins=abins)
    eishist, b = np.histogram(eisflux, bins=abins)

    
    plotpts = []
    widths = []
    for b in range(len(abins[:-1])): #creating points to plot - midpoints of bins
        plotpts.append(abins[b] + (abins[b+1]-abins[b])/2)
        widths.append((abins[b+1]-abins[b]))

    flux_grid = np.linspace(-1,3,100)
    flux_kde = gaussian_kde(B7flux)
    flux_pdf = flux_kde.evaluate(flux_grid)
    norm_flux_pdf = flux_pdf*widths[0]*np.sum(B7hist)
    plt.plot(flux_grid, norm_flux_pdf, label='Band 7 KDE')

    eis_kde = gaussian_kde(eisflux)
    eis_pdf = eis_kde.evaluate(flux_grid)
    norm_eis_pdf = eis_pdf*widths[0]*np.sum(eishist)
    plt.plot(flux_grid, norm_eis_pdf, label='E18 KDE')

    print(flux_kde.n, eis_kde.n)
    
    Dval, pval = ks_2samp(B7flux, eisflux)
    print(f'{pval} p value for ks test of b7/eisner')
    
    plt.bar(plotpts, B7hist, widths, alpha=0.5, edgecolor='black', label=f'Band 7 {sources} sources')
    plt.bar(plotpts, eishist, widths, alpha=0.5, edgecolor='black', label='E18')
    plt.legend()
    plt.ylabel('number')
    plt.xlabel(r'$\log(F_{\nu=350 GHz} (mJy))$')
    plt.savefig(f'/home/jotter/nrao/plots/E18_comp/flux_hist_eis_B7_{sources}.pdf', dpi=400)
    plt.close()
    
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
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7':339.7672758867*u.GHz}  #first two from Adam, third from member.uid header
    data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/r0.5_catalog_best.fits')

    flux_ind = reduce(np.intersect1d, (np.where(np.isnan(data['ap_flux_'+x_name])== False)[0],np.where(np.isnan(data['ap_flux_'+y_name])==False)[0],np.where(np.isnan(data['ap_flux_'+num_name])==False)[0])) #sources with measured aperture fluxes
    fit_ind = reduce(np.intersect1d, (np.where(data['fit_flag_'+x_name]=='y')[0],np.where(data['fit_flag_'+y_name]=='y')[0],np.where(data['fit_flag_'+num_name]=='y')[0])) #sources w good gaussian fits in both bands
    ind = np.intersect1d(flux_ind, fit_ind)
    

    plot1 = data['ap_flux_'+x_name][ind]
    plot2 = data['ap_flux_'+y_name][ind]
    plot1_err = data['ap_flux_err_'+x_name][ind]
    plot2_err = data['ap_flux_err_'+y_name][ind]
    plot_num = data['ap_flux_'+num_name][ind]
    plot_num_err = data['ap_flux_err_'+num_name][ind]

    plotx = plot1/plot_num
    ploty = plot_num/plot2
    plotx_err = (plot1/plot_num)*np.sqrt((plot1_err/plot1)**2 + (plot_num_err/plot_num)**2) #error is fractional error added in quad
    ploty_err = (plot_num/plot2)*np.sqrt((plot2_err/plot2)**2 + (plot_num_err/plot_num)**2)

    plt.errorbar(plotx, ploty, xerr=plotx_err, yerr=ploty_err, linestyle='', marker='*')

    plt.xlabel('{x}/{num}'.format(num=num_name,x=x_name))
    plt.ylabel('{num}/{y}'.format(num=num_name,y=y_name))

    totmin = np.amin((plotx, ploty))
    totmax = np.amax((plotx, ploty))
    plt.xlim(totmin-0.05, totmax+0.05)
    plt.ylim(totmin-0.05, totmax+0.05)
    
    #plot alphas:
    alpha = [1.5,2,2.5]
    cols = ['r', 'g', 'k']
    for i,a in enumerate(alpha):
        plt.axvline((freqs[x_name]/freqs[num_name])**a, linestyle='--', label=r'$\alpha={a}$'.format(a=a), color=cols[i])
        plt.axhline((freqs[num_name]/freqs[y_name])**a, linestyle='--', color=cols[i])
    plt.xlim(1,4)
    plt.ylim(1,9)
    plt.legend()
    plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/colorcolor_plot.png', dpi=300)


def plot_SEDs(names):
    freqs = {'B3':98*u.GHz, 'B6':223.5*u.GHz, 'B7':339.7672758867*u.GHz}
    data = fits.getdata('/users/jotter/summer_research_2018/tables/r0.5_catalog_conv_bgfitted.fits')
    ind = get_ind(names)

    alpha = [1.5,2,2.5]
    freq_x = np.array([freqs[n].value for n in names])
    for i in ind:
        fluxes = [data['ap_flux_'+n][i] for n in names]
        #fluxes = np.array(fluxes) - np.array([data['bg_median_'+n][i] for n in names])
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
        plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/SEDs/SED_'+str(i)+'_'+names[-1]+'_updt_SED.png', dpi=300)

def multi_img_SED(srcID, B3_img, B6_img, B7_img, name, flux_type = 'amp'):

    #srcID = 8
    #B3_img = 'r-2.clean0.1mJy.500klplus.deepmask'
    #B6_img = 'r-2.clean0.1mJy.500klplus.deepmask'
    #B7_img = 'r-2.clean0.1mJy.500klplus.deepmask'
    #name = 'r-2.0.1.500klplusdeepmask'
    #flux_type='ap'

    imgs = [B3_img, B6_img, B7_img]

    freqs = {'B3':98, 'B6':223.5, 'B7':339.7672758867, 'Fb':6.099018346128}
    B3data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B3_500klplus_img_catalog_B6_ref.txt')
    B6data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B6_500klplus_img_catalog_B6_ref.txt')
    B7data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/B7_500klplus_img_catalog_B6_ref.txt')
    catalogs = [B3data, B6data, B7data]

    B3ind = np.where(B3data['D_ID'] == srcID)
    B6ind = np.where(B6data['D_ID'] == srcID)
    B7ind = np.where(B7data['D_ID'] == srcID)
    inds = [B3ind, B6ind, B7ind]

    alpha = [1,1.5,2,2.5,3]
    names = ['B3', 'B6', 'B7']

    freq_x = []
    fluxes = []
    flux_err = []
    bgs = []
    
    for i in range(len(names)):
        freq_x.append(freqs[names[i]])
        if flux_type == 'ap':
            fluxes.append(catalogs[i]['ap_flux_'+imgs[i]][inds[i]][0]*1000)
            flux_err.append(catalogs[i]['ap_flux_err_'+imgs[i]][inds[i]]*1000)
            bgs.append(catalogs[i]['bg_ap_'+imgs[i]][inds[i]].data.data[0]*1000)
        elif flux_type == 'amp':
            fluxes.append(catalogs[i]['g_amplitude_'+imgs[i]][inds[i]][0]*1000)
            flux_err.append(catalogs[i]['amp_err_'+imgs[i]][inds[i]]*1000)
            bgs.append(catalogs[i]['bg_median_'+imgs[i]][inds[i]].data.data[0]*1000)
        elif flux_type == 'circ':
             fluxes.append(catalogs[i]['circ_flux_'+imgs[i]][inds[i]][0]*1000)
             flux_err.append(catalogs[i]['circ_flux_err_'+imgs[i]][inds[i]]*1000)
             bgs.append(catalogs[i]['bg_circ_'+imgs[i]][inds[i]].data.data[0]*1000)

    plt.figure()

    Fb_noise_lim = 0.01 #mJy
    if flux_type == 'ap':
        Fb_flux = catalogs[0]['Speak_Fb'][inds[0][0]]
        if np.isnan(Fb_flux) == True:
            upper_lim = Fb_noise_lim
            plt.errorbar(freqs['Fb'], upper_lim, marker='v', uplims=True)
        else:
            fluxes.append(Fb_flux)
            flux_err.append(catalogs[0]['e_Speak_Fb'][inds[0][0]])
            bgs.append(0)
            freq_x.append(freqs['Fb'])
            Fb_alpha = catalogs[0]['alpha'][inds[0][0]]
            b = np.log10(Fb_flux/(freqs['Fb']**Fb_alpha))
            nuvals = np.array([5,8])
            Fvals = 10**b * nuvals**Fb_alpha
            plt.plot(nuvals, Fvals)
            

    F2 = fluxes[1] - bgs[1]
    nu2 = freq_x[1]
        
    
    freq_x = np.array(freq_x)
    fluxes = np.array(fluxes) 
    bgs = np.array(bgs)

    plt.errorbar(freq_x, fluxes-bgs, yerr=flux_err, linestyle='', marker='o')

    freq_x = np.array([freqs[fr] for fr in freqs])
    for j in range(len(alpha)):
        F1 = F2*((freq_x/nu2)**alpha[j])
        plt.plot(freq_x, F1, linestyle='-', label=r'$\alpha = $'+str(alpha[j]))
		
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('aperture flux (mJy)')
    plt.xlabel('frequency (GHz)')
    plt.legend(loc='best')
    plt.title(str(srcID)+'SED')

    locs_x = [0.3, 0.5, 0.7]
    locs_y = [0.15, 0.15, 0.15]
	
    for i in range(len(imgs)):
        nm = names[i]
        img_file = '/lustre/aoc/students/jotter/directory/Orion'+nm+'/Orion_SourceI_'+nm+'_continuum_'+imgs[i]+'.image.tt0.pbcor.fits'
        img_fl = fits.open(img_file)
        img_data = img_fl[0].data.squeeze()
        img_header = img_fl[0].header
        mywcs = WCS(img_header).celestial
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg 
        if nm == 'B7':
            nm = 'B7_hr'
        center_coord = SkyCoord(catalogs[i]['gauss_x_'+nm][inds[i]], catalogs[i]['gauss_y_'+nm][inds[i]], frame='icrs', unit=(u.deg, u.deg))
        center_coord_pix = center_coord.to_pixel(mywcs)
        if i == 0:
            pix_major_fwhm = ((catalogs[i]['FWHM_major_'+nm][inds[i]]*u.arcsec).to(u.degree)/pixel_scale).decompose()
            size = 2.3*pix_major_fwhm.value
        cutout = Cutout2D(img_data, center_coord_pix, size[0], mywcs, mode='partial')
        a = plt.axes([locs_x[i], locs_y[i], .2, .2])
        plt.imshow(cutout.data, origin='lower')
        plt.title(names[i])
        a.tick_params(labelleft=False, labelbottom=False)
    if flux_type == 'circ':
        plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/SEDs/500klplus/SED_circ_'+str(srcID)+'_'+name+'_SED.png', dpi=300)
    else:
        plt.savefig('/users/jotter/summer_research_2018/flux_plots/plots/SEDs/500klplus/SED_'+str(srcID)+'_'+name+'_SED.png', dpi=300)
    

def image_fluxes(srcID, flux_type='ap'):
    ref_names = ['B3', 'B6', 'B7']
    flagged_B3 = ['r2', 'r2_dirty','r0.5_dirty', 'r2.clean2mJy.200mplus', 'r2.clean2mJy.150mplus', 'r2.clean2mJy.50mplus', 'r2.clean1mJy.200mplus.deepmask', 'r2.clean1mJy.150mplus.deepmask', 'r2.clean1mJy.50mplus.deepmask', 'r2.clean2mJy.allbaselines', 'r2.clean1mJy.allbaselines.deepmask']
    flagged_B6 = ['r2.clean0.5mJy.selfcal.phase4.allbaselines.highressta', 'r2.clean1mJy.deepmask.allbaselines', 'r2.clean2mJy.allbaselines', 'r2.clean2mJy.selfcal.phase4.allbaselines', 'r2.clean1mJy.selfcal.phase4.deepmask.allbaselines', 'r2_dirty', 'r2.clean.1mJy.selfcal.phase4.deepmask.50mplus', 'r2.clean1mJy.1500kplus.deepmask', 'r-2.clean15.0mJy.pmcheck_100to3200m.2016.deepmask', 'r-2.clean35.0mJy.pmcheck_100to3200m.2016', 'r2', 'r0.5_dirty', 'r2.clean2mJy.selfcal.phase4.200mplus', 'r0.5', 'r2.clean2mJy.200mplus', 'r0.5.clean10mJy.pmcheck_100-3200m.2016.deepmask', 'r2.clean1mJy.selfcal.phase4.deepmask.200mplus', 'r2.clean1mJy.200mplus.deepmask', 'r2.clean2mJy.selfcal.phase4.50mplus', 'r0.5.clean1mJy.selfcal.phase4.200mplus', 'r0.5.clean1mJy.selfcal.phase4.150mplus', 'r2.clean2mJy.selfcal.phase4.50mplus', 'r0.5.clean0.5mJy.selfcal.phase4.deepmask.200mplus', 'r0.5.clean1mJy.200mplus', 'r0.5.clean0.5mJy.selfcal.phase4.deepmask.150mplus', 'r0.5.clean0.5mJy.150mplus.deepmask', 'r2.clean2mJy.selfcal.phase4.150mplus', 'r2.clean1mJy.selfcal.phase4.deepmask.150mplus', 'r2.clean2mJy.50mplus', 'r2.clean1mJy.50mplus.deepmask', 'r-2.clean0.4mJy.selfcal.ampphase6', 'r0.5.clean5.0mJy.pmcheck_100to3200m.2017.deepmask', 'r0.5.clean7.5mJy.pmcheck_100to3200m.2017', 'r0.5.clean25.0mJy.pmcheck_100to3200m.2016', 'r0.5.clean10.0mJy.pmcheck_100to3200m.2016.deepmask', 'r2.clean1mJy.selfcal.phase4.deepmask.50mplus']
    flagged_B7 = ['r2', 'r0.5', 'r2_dirty', 'r0.5_dirty', 'r0.5.clean0.5mJy.50mplus.deepmask']
    flagged_imgs = [flagged_B3, flagged_B6, flagged_B7]


    for ind,ref_name in enumerate(ref_names):
        data = ascii.read('/lustre/aoc/students/jotter/dendro_catalogs/'+ref_name+'_diff_imgs_catalog_B6_ref.txt')
        srcind = np.where(data['D_ID'] == srcID)
        flux_names = fnmatch.filter(data.colnames, 'ap_flux_r*')
        charnum = 8
        bg_names = fnmatch.filter(data.colnames, 'bg_ap_r*')

        if flux_type == 'circ':
            flux_names = fnmatch.filter(data.colnames, 'circ_flux_r*')
            charnum = 10
            bg_names = fnmatch.filter(data.colnames, 'bg_circ_r*')


        sub_fluxes = []
        for n in range(len(flux_names)):
            if flux_names[n][8:] not in flagged_imgs[ind]:
                sub_fluxes.append([data[flux_names[n]][srcind] - data[bg_names[n]][srcind], flux_names[n][charnum:], data[bg_names[n]][srcind]])
        sub_fluxes.sort()
        
        x = np.arange(len(sub_fluxes))
        fluxes = np.array([sf[0] for sf in sub_fluxes])
        bgs = np.array([sf[2] for sf in sub_fluxes])
        plt.figure(figsize=(10,10))
        plt.semilogy(x, fluxes, marker='o', label='bg subtracted')
        plt.semilogy(x, fluxes + bgs, marker='o', label=flux_type+' fluxes')
        plt.xticks(x, [sf[1] for sf in sub_fluxes], rotation='vertical')
        plt.subplots_adjust(bottom=0.5) 
        plt.grid()
        plt.legend()
        plt.title(ref_name+' source '+str(srcID)+' fluxes')  
    plt.show()
    

#flux_hist(sources='OMC1')
#flux_hist(sources='ONC')
#flux_hist(sources='all')
