from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import radio_beam
import astropy.units as u

data = fits.getdata('/lustre/aoc/students/jotter/dendro_catalogs/IR_matched_catalog.fits')
B3nu = 93
B6nu = 223.5

B6_img = fits.open('/lustre/aoc/students/jotter/directory/Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines.image.tt0.pbcor.fits')
B6_head = B6_img[0].header
beam_B6 = radio_beam.Beam.from_fits_header(B6_head)

B3_img = fits.open('/lustre/aoc/students/jotter/directory/Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits')
B3_head = B3_img[0].header
beam_B3 = radio_beam.Beam.from_fits_header(B3_head)

#B3/B6 plot:
ind = np.intersect1d(np.where(np.isnan(data['ap_flux_B6'])== False)[0],np.where(np.isnan(data['ap_flux_B3'])==False)[0])
size_inds = [ 3,  6,  8, 10, 11, 14, 16, 17] #sources which could be deconvolved and have good B3/B6 gaussian fits
ind2 = np.intersect1d(np.where(data['fit_goodness_B6'] =='y'), np.where(data['fit_goodness_B3'] =='y')) #sources w good gaussian fits in both bands
all_ind = np.intersect1d(ind, ind2)

del_ind = np.delete(ind, size_inds)
all_del_ind = np.intersect1d(all_ind, del_ind)

B6_fwhm_major = data['FWHM_major_B6']
B3_fwhm_major = data['FWHM_major_B3']

#B6_fwhm_major_corrected = np.sqrt(data['FWHM_major_B6']*data['FWHM_minor_B6'] - beam_B6.major.to(u.arcsec).value*beam_B6.minor.to(u.arcsec).value)
#B3_fwhm_major_corrected = np.sqrt(data['FWHM_major_B3']*data['FWHM_minor_B3'] - beam_B6.major.to(u.arcsec).value*beam_B6.minor.to(u.arcsec).value)

#B6_amp_corrected = (data['gauss_amplitude_B6']*(B6_fwhm_major_corrected**2)/(beam_B6.major.to(u.arcsec)*beam_B6.minor.to(u.arcsec)))

B6_amp_corrected = (data['gauss_amplitude_B6']*data['FWHM_major_B6']*data['FWHM_minor_B6']/(beam_B6.major.to(u.arcsec)*beam_B6.minor.to(u.arcsec)))
B3_amp_corrected = (data['gauss_amplitude_B3']*data['FWHM_major_B3']*data['FWHM_minor_B3']/(beam_B3.major.to(u.arcsec)*beam_B3.minor.to(u.arcsec)))

B6_amp_err_corrected = (data['amplitude_err_B6']*data['FWHM_major_B6']*data['FWHM_minor_B6']/(beam_B6.major.to(u.arcsec)*beam_B6.minor.to(u.arcsec)))
B3_amp_err_corrected = (data['amplitude_err_B3']*data['FWHM_major_B3']*data['FWHM_minor_B3']/(beam_B3.major.to(u.arcsec)*beam_B3.minor.to(u.arcsec)))

B3_amp = data['gauss_amplitude_B3']
B3_amp_err = data['amplitude_err_B3']
B6_amp = data['gauss_amplitude_B6']
B6_amp_err = data['amplitude_err_B6']

B3_plot = data['ap_flux_B3']
B3_err = data['ap_flux_err_B3']
B6_plot = data['ap_flux_B6']
B6_err = data['ap_flux_err_B6']

F1 = np.linspace(np.min(B3_plot[all_ind]), np.max(B3_plot[all_ind]), 10)
alpha05_F2 = F1*((B6nu/B3nu)**0.5)
alpha1_F2 = F1*(B6nu/B3nu)
alpha15_F2 = F1*((B6nu/B3nu)**1.5)
alpha2_F2 = F1*((B6nu/B3nu)**2)
alpha25_F2 = F1*((B6nu/B3nu)**2.5)

fig1 = plt.figure()
plt.xlabel('B3 aperture flux')
plt.ylabel('B6 aperture flux')
plt.errorbar(B3_plot[all_ind], B6_plot[all_ind], xerr=B3_err[all_ind], yerr=B6_err[all_ind], linestyle='', marker='*',color='blue',label='aperture flux')
#plt.errorbar(B3_plot[ind[size_inds]], B6_plot_deconv[ind[size_inds]], xerr=B3_err_deconv[ind[size_inds]], yerr=B6_err_deconv[ind[size_inds]], linestyle='', marker='o',color='blue')
plt.errorbar(B3_amp_corrected[all_ind].value, B6_amp_corrected[all_ind].value, xerr=B3_amp_err_corrected[all_ind].value, yerr=B6_amp_err_corrected[all_ind].value, linestyle='', marker='*', color='red',label='corrected gaussian amplitude')
#plt.errorbar(B3_amp_corrected, B6_amp_corrected, xerr=B3_amp_err, yerr=B6_amp_err, linestyle='', marker='*',color='red')

plt.loglog(F1, alpha05_F2, linestyle=':',label='alpha=0.5')
plt.loglog(F1, alpha1_F2, linestyle=':',label='alpha=1')
plt.loglog(F1, alpha15_F2, linestyle=':',label='alpha=1.5')
plt.loglog(F1, alpha2_F2, linestyle=':',label='alpha=2')
plt.loglog(F1, alpha25_F2, linestyle=':',label='alpha=2.5')

plt.legend()

plt.savefig('/users/jotter/summer_research_2018/flux_plots/B3B6_spectral_indices.png')
fig4 = plt.figure()

alpha_gauss = np.log((B3_amp_corrected[all_ind]/B6_amp_corrected[all_ind]).decompose())/np.log(B3nu/B6nu)
g_bins = plt.hist(alpha_gauss, alpha=0.3, label='corrected gaussian amplitude')


alpha_apflux = np.log(B3_plot[all_ind]/B6_plot[all_ind])/np.log(B3nu/B6nu)
plt.hist(alpha_apflux, alpha=0.3, label='aperture flux')#, bins=g_bins[1])

plt.legend()
plt.xlabel('alpha')
plt.savefig('/users/jotter/summer_research_2018/flux_plots/B3B6_specind_hist.png')
plt.show()
#Diagnostic plots below
'''
fig2 = plt.figure()
plt.ylabel('B6 aperture flux/amplitude')
plt.xlabel('B6 FWHM major')

plt.scatter(B6_fwhm_major[ind[size_inds]], B6_plot[ind[size_inds]]/B6_amp_corrected[ind[size_inds]])
plt.scatter(B6_fwhm_major[all_del_ind], B6_plot[all_del_ind]/B6_amp[all_del_ind])
plt.scatter(B6_fwhm_major[all_del_ind], B6_plot[all_del_ind]/B6_amp_corrected[all_del_ind])

plt.hlines(1,0.03,0.08)

fig3 = plt.figure()
plt.ylabel('B3 aperture flux/amplitude')
plt.xlabel('B3 FWHM major')

plt.scatter(B3_fwhm_major[ind[size_inds]], B3_plot[ind[size_inds]]/B3_amp_corrected[ind[size_inds]])
plt.scatter(B3_fwhm_major[all_del_ind], B3_plot[all_del_ind]/B3_amp[all_del_ind])
plt.scatter(B3_fwhm_major[all_del_ind], B3_plot[all_del_ind]/B3_amp_corrected[all_del_ind])

plt.hlines(1,0.05,0.12)

'''
'''
#B3 Kmag plot
ind = np.concatenate((np.where(np.isnan(data['flux_B3'])== False)[0],np.where(np.isnan(data['Kmag'])==False)[0]))
B3_plot_kmag = data['flux_B3'][ind]
kmag_plot_B3 = data['Kmag'][ind]

fig2 = plt.figure()
plt.xlabel('B3 flux')
plt.ylabel('K-band magnitude')
plt.semilogx(B3_plot_kmag, kmag_plot_B3, linestyle='', marker='*')
plt.ylim(15,9)

#B3 Hmag plot
ind = np.concatenate((np.where(np.isnan(data['flux_B3'])== False)[0],np.where(np.isnan(data['Hmag'])==False)[0]))
B3_plot_hmag = data['flux_B3'][ind]
hmag_plot_B3 = data['Hmag'][ind]

fig3 = plt.figure()
plt.xlabel('B3 flux')
plt.ylabel('H-band magnitude')
plt.semilogx(B3_plot_hmag, hmag_plot_B3, linestyle='', marker='*')
plt.ylim(17,9)
plt.show()'''
