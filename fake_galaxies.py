#import warnings
#warnings.filterwarnings("ignore", message="using a non-integer number instead of an integer will result in an error in the future")

def template_names():
	import glob as glob

	template_files = glob.glob('miles_models/Mun*.fits')
	return template_files

def choose_templates(templates, age_lim = 20.0, max_nonzero = 5):
	#start out by loading in the template files as a table

	import numpy as np
	import astropy.table as table
	import SPaCT

	ssp_rows = []
	for template in templates:
	    template = template.rstrip('.fits').split('/')[1]
	    
	    spectral_range = template[0]
	    IMF_type = template[1:3]
	    IMF_slope = float(template[3:7])
	    Z = SPaCT.plusminus(template[8])*float(template[9:13])
	    T = float(template[14:])
	    
	    #print template + ':', spectral_range, IMF_type, IMF_slope, Z, T
	    ssp_i = [template, spectral_range, IMF_type, IMF_slope, Z, T]
	    ssp_rows.append(ssp_i)
	    
	ssps = table.Table(map(list, zip(*ssp_rows)), names = ['name', 'spectral range', 'IMF type', 'IMF slope', 'Z', 't'])
	ssps = ssps[ssps['t'] <= age_lim]
	#then pick up to `max_nonzero` number of templates to be nonzero
	nonzero_templates = np.random.choice(ssps['name'], np.random.randint(1, max_nonzero + 1), replace = False)
	template_weights = np.random.rand(len(ssps['name'])) * [1. if i in nonzero_templates else 0. for i in ssps['name']]
	template_weights /= template_weights.sum()
	
	ssps.add_column(table.Column(name = 'weight', data = template_weights))

	return ssps

def generate_spectrum(ssps):
	'''
	generate a pristine spectrum based on weights given in an astropy table of templates
	'''

	import numpy as np
	from astropy.io import fits
	import astropy.table as table

	#now load in each template as a row in an array
	all_templates = np.empty([len(ssps['name']), fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['NAXIS1']])

	for i, row in enumerate(all_templates):
		all_templates[i] = ssps['weight'][i] * fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].data
		#if ssps['weight'][i] != 0: print all_templates[i]

	clean_spectrum = all_templates.sum(axis = 0)
	CRVAL1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['CRVAL1']
	CDELT1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['CDELT1']
	NAXIS1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['NAXIS1']
	l_full = CRVAL1 + np.linspace(0., NAXIS1 * CDELT1, NAXIS1)
	clean_spectrum /= np.median(clean_spectrum)
	return clean_spectrum, l_full

def generate_LOSVD(spectrum, v_res, moments, plots = False):
	'''
	Convolve `spectrum` with a Gaussian-like filter, except with nonzero higher-order moments.

	This reproduces a velocity field that pPXF will fit

	NOTE: nonzero higher-order moments not supported at this time
	NOTE: nonzero m1 is not supported (and is a very bad idea) - always use redshift routine to apply this!
	'''
	import numpy as np
	import scipy.ndimage as ndimage
	import matplotlib.pyplot as plt

	#generate a kernel with the given moments
	m1, m2, m3, m4, m5, m6 = moments

	if m1 != 0.: 
		while a not in ['y', 'n']:
			a = raw_input('Warning! non-zero-centered LOSVDs are not recommended! Proceed? (y/n)')
			if a == 'y': break
			elif a == 'n': exit()

	if moments[2:] != [0., 0., 0., 0.]:
		raise ValueError('only nonzero higher-order G-H moments are supported!')
	else:
		spectrum_LOSVD = ndimage.gaussian_filter1d(spectrum, m2/v_res)

	if plots == True:
		plt.figure(figsize = (6, 4))
		plt.plot(spectrum, c = 'b', label = 'rest-frame')
		plt.plot(spectrum_LOSVD, c = 'g', label = 'LOSVD spectrum')
		plt.plot(np.abs(spectrum - spectrum_LOSVD), label = 'residual')
		plt.legend(loc = 'best')
		plt.show()

	return spectrum_LOSVD

def redshift_spectrum(l_0, z = None, dz = None):
	#redshift a spectrum randomly, and return the new wavelength array, a "real" redshift, and a redshift measurement error

	import numpy as np

	if z == None:
		z = np.random.uniform(0.01, 0.025)
	if dz == None:
		dz = np.sign(np.random.random() - 0.5) * (10**(np.random.uniform(-1.0, -0.5))) * z #random error beween 1% and 10%, equal probabilities of + and -

	#print z, dz
	l_1 = l_0 * (1. + z + dz)
	return z, dz, l_1

def adjust_FWHM(sharp_spectrum, res_old, res_new, FWHM_old, FWHM_new):
	#convolve the spectrum with a Gaussian with a width of the square root of the difference of the squares of the intrument FWHMs

	import numpy as np
	import scipy.ndimage as ndimage

	assert FWHM_new >= FWHM_old
	FWHM_dif = np.sqrt(FWHM_new**2. - FWHM_old**2.)
	sigma_diff = FWHM_dif/2.355/res_old # Sigma difference in pixels
	blurred_spectrum = ndimage.gaussian_filter1d(sharp_spectrum, sigma_diff)
	return blurred_spectrum

def downsample_spectrum(l_dense, dense_spectrum, l_sparse):
	#linearly interpolate the input dense_spectrum (which has values at all the locations in l_0), to the values in l_1

	import numpy as np
	import scipy.interpolate as interp

	sparse_spectrum = interp.interp1d(l_dense, dense_spectrum, kind = 'linear')(l_sparse)

	return l_sparse, sparse_spectrum

def noisify_ifu(spectrum, n, SNR):
	#make some number `n` of rows with pure noise, and add similar noise profile to first row
	import numpy as np

	NAXIS1 = len(spectrum)

	raw_noise_IFU = np.random.normal(loc = 0.0, scale = 1./SNR, size = (n, NAXIS1))

	empty_fibers = raw_noise_IFU * np.tile(spectrum, reps = (n, 1))
	IFU = np.vstack((spectrum, empty_fibers))
	galaxy_noise = np.random.normal(loc = 0.0, scale = 1./SNR, size = NAXIS1) * spectrum
	IFU[0] = spectrum + galaxy_noise
	return IFU, galaxy_noise

def population_sum_models(ssps):
	#take a table of templates and weights, and return a spectrum in the specified range
	import numpy as np
	import scipy.ndimage as ndimage
	import astropy.table as table
	from astropy.io import fits

	#now load in each template as a row in an array
	all_templates = np.empty([len(ssps['name']), fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['NAXIS1']])

	for i, row in enumerate(all_templates):
		all_templates[i] = ssps['weight'][i] * fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].data
		#if ssps['weight'][i] != 0: print all_templates[i]

	real_spectrum = all_templates.sum(axis = 0)
	CRVAL1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['CRVAL1']
	CDELT1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['CDELT1']
	NAXIS1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['NAXIS1']
	l_full = CRVAL1 + np.linspace(0., NAXIS1 * CDELT1, NAXIS1)
	real_spectrum /= np.median(real_spectrum)

	return real_spectrum, l_full

def population_sum_fit(ssps):
	import numpy as np
	import scipy.ndimage as ndimage
	import astropy.table as table
	from astropy.io import fits

	#now load in each template as a row in an array
	all_templates = np.empty([len(ssps['name']), fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['NAXIS1']])

	for i, row in enumerate(all_templates):
		all_templates[i] = ssps['best-fit weights'][i] * fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].data
		#if ssps['weight'][i] != 0: print all_templates[i]

	derived_spectrum = all_templates.sum(axis = 0)
	CRVAL1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['CRVAL1']
	CDELT1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['CDELT1']
	NAXIS1 = fits.open('miles_models/' + ssps['name'][0] + '.fits')[0].header['NAXIS1']
	l_full = CRVAL1 + np.linspace(0., NAXIS1 * CDELT1, NAXIS1)
	derived_spectrum /= np.median(derived_spectrum)

	return derived_spectrum, l_full

def pPXF_summary_plots(ssps, instrument_info, pp, lam_sparse, vel, verbose = False):
	#make sure `vel` is the sum of the redshift and the kinematic velocity fit
	import numpy as np
	import matplotlib.pyplot as plt
	import astropy.table as table
	import colorpy.ciexyz as ciexyz
	import colorpy.colormodels as cmodels
	import warnings

	c = 299792.458

	if verbose == True:
		print 'non-zero fit templates:'
		print ssps[ssps['best-fit weights'] != 0.]['Z', 't', 'best-fit weights']
		print 'non-zero real solution templates:'
		print ssps[ssps['weight'] != 0.]['Z', 't', 'weight']

	#first plot the original and resultant populations
	f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex = True, figsize = (8, 6))
	ax1.set_title('fit')
	a = ax1.scatter(ssps['Z'], ssps['t'], c = pp.weights, cmap = 'gnuplot', s = 40, vmin = 0.0, vmax = 1.0, edgecolor = 'grey')
	ax2.set_title('reality')
	ax2.scatter(ssps['Z'], ssps['t'], c = ssps['weight'], cmap = 'gnuplot', s = 40, vmin = 0.0, vmax = 1.0, edgecolor = 'grey')
	plt.colorbar(a)
	plt.suptitle('population fit comparison', size = 16)
	plt.show()

	#now plot the result with the input

	instrument_lam_lims = (instrument_info['CRVAL1'], instrument_info['CRVAL1'] + instrument_info['NAXIS1'] * instrument_info['CDELT1'])

	lines = [
		['Ca H', 3968.5], ['Ca K', 3933.7], ['H-alpha', 6562.8], ['H-beta', 4861.], ['Mg I', 5175.], ['Ca I', 4307.]
		]

	#plt.figure(figsize = (10, 6))
	#ax = plt.subplot(111)
	#ax.plot(lam_sparse, pp.bestfit)
	print  'vel:', vel, 'km/s'

	#now plot relevant spectral lines
	'''
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", category = DeprecationWarning)
		for i, line in enumerate(lines):
			line_c = cmodels.irgb_string_from_xyz(ciexyz.xyz_from_wavelength(line[1]/10.))
			#print line_c
			ax.axvline(line[1] * (1. + vel / c), color = line_c)
			ax.annotate(line[0], xy = (line[1], 1.2), xytext = (line[1]+10., 1.1 - 0.1 * i%2), size = 14)
	'''

	#plt.show()

	'''
	real_spectrum, l_full = population_sum_models(ssps = ssps)
	derived_spectrum, l_full = population_sum_fit(ssps = ssps)
	plt.figure(figsize = (10, 6))
	ax1 = plt.subplot(211)
	#first at full resolution
	ax1real = ax1.plot(l_full, real_spectrum, c = 'g', label = 'Reality', linewidth = 0.25)
	ax1der = ax1.plot(l_full, derived_spectrum, c = 'b', linestyle = '--', label = 'Fit', linewidth = 0.25)
	for val in instrument_lam_lims:
		ax1.axvline(val, c = 'r', linestyle = ':')
	ax1.set_title('Full-Resolution spectra', size = 16)

	ax1_1 = ax1.twinx()
	ax1err = ax1_1.plot(l_full, np.abs(real_spectrum - derived_spectrum), linewidth = 0.25, c = 'tomato', label = 'Error')
	for tl in ax1_1.get_yticklabels(): tl.set_color('tomato')

	ax1_l = ax1_1.legend(ax1real + ax1der + ax1err, [l.get_label() for l in (ax1real + ax1der + ax1err)], loc = 'best')
	ax1_l.set_zorder(5)

	ax2 = plt.subplot(212, sharex = ax1)
	ax2.set_title('Downsampled spectra', size = 16)
	ax2.set_xlabel(r'$\lambda[\AA]$', size = 16)
	
	#now after blurring and downsampling
	l_sparse = np.linspace(instrument_info['CRVAL1'], instrument_info['CRVAL1'] + instrument_info['NAXIS1'] * instrument_info['CDELT1'], instrument_info['NAXIS1'])
	l_sparse, sparse_spectrum_real = downsample_spectrum(l_dense = l_full, dense_spectrum = real_spectrum, l_sparse = l_sparse) #this accomplishes both downsampling and paring!!
	l_sparse, sparse_spectrum_der = downsample_spectrum(l_dense = l_full, dense_spectrum = derived_spectrum, l_sparse = l_sparse) #this accomplishes both downsampling and paring!!

	ax2real = ax2.plot(l_sparse, sparse_spectrum_real, c = 'g', label = 'Reality', linewidth = 0.25)
	ax2der = ax2.plot(l_sparse, sparse_spectrum_der, c = 'b', label = 'Fit', linewidth = 0.25, linestyle = '--')
	for val in instrument_lam_lims:
		ax2.axvline(val, c = 'r', linestyle = ':')
	ax2.set_title('Downsampled spectra', size = 16)
	ax2.set_xlabel(r'$\lambda[\AA]$', size = 16)
	ax2_1 = ax2.twinx()
	
	ax2err = ax2_1.plot(l_sparse, np.abs(sparse_spectrum_real - sparse_spectrum_der), linewidth = 0.25, c = 'tomato', label = 'Error')
	for tl in ax2_1.get_yticklabels(): tl.set_color('tomato')

	ax2_l = ax2_1.legend(ax2real + ax2der + ax2err, [l.get_label() for l in (ax2real + ax2der + ax2err)], loc = 'best')
	ax2_l.set_zorder(5)

	plt.tight_layout()
	plt.show()
	'''

def simulate_noise(sparse_spectrum, SNR, n_skyfiber_range = [1, 20, 3]):
	'''
	generate synthetic noise spectra for a given input spectrum, and test the required number of sky fibers (with similar noise profiles) to accurately get the SNR
	'''

	import numpy as np
	import SPaCT
	import matplotlib.pyplot as plt

	plt.figure(figsize = (6, 4))

	for n_skyfibers in range(n_skyfiber_range[0], n_skyfiber_range[1] + 1, n_skyfiber_range[2]):
		ifu, galaxy_noise = noisify_ifu(sparse_spectrum, n = n_skyfibers, SNR = SNR)
		fiberlist = range(1, n_skyfibers + 1)
		SNR_calc = ifu[0] / SPaCT.noise_edgefibers(ifu, width = 3, fiberlist = fiberlist, verbose = False)
		bins, edges = np.histogram(SNR_calc, 50, normed = 1)
		left, right = edges[:-1],edges[1:]
		X = np.array([left,right]).T.flatten()
		Y = np.array([bins,bins]).T.flatten()
		plt.plot(X, Y/Y.max(), label = str(n_skyfibers) + ' fibers')

	plt.axvline(SNR, c = 'k', linestyle = ':')
	SNR_annotation = plt.text(SNR, 0.35, '$S/N=' + str(SNR) + '$')
	SNR_annotation.set_rotation('vertical')
	plt.title('Effect of # of sky fibers on SNR', size = 18)
	plt.xscale('log')
	plt.ylim([-0.05, 1.05])
	plt.xlabel('SNR', size = 18)
	plt.ylabel('normed fraction', size = 18)
	plt.legend(loc = 'best', prop = {'size':6})
	plt.tight_layout()
	plt.show()

def simulate_single_spectrum():
	'''
	STEPS:
		1. choose templates
		2. make spectrum
			2.1. convolve with a LOSVD
		3. blur to correct FWHM
		4. redshift to within some error
		5. downsample to correct wavelengths
		6. noisify and create an IFU with the same noise characteristics
		7. run pPXF
	'''

	from astropy.io import fits
	from astropy import table
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	import SPaCT
	import scipy.stats as stats
	from ppxf import robust_sigma

	import warnings

	SPSPK_info = fits.open('NGC2558.msobj.fits')[0].header

	template_files = template_names()

	ssps = choose_templates(templates = template_files, max_nonzero = 4)
	clean_spectrum, l_full = generate_spectrum(ssps = ssps)
	#now redshift the spectrum
	MILES_res = fits.open(template_files[0])[0].header['CDELT1']
	SPSPK_res = 1.4
	#FWHMs should be in Angstroms
	FWHM_MILES = 1.36
	FWHM_SPSPK = 4.877 #this is specific to one particular configuration, so handle with care!
	SNR = 100.
	n_skyfibers = 8
	n_moments = 4 #how many moments to fit

	'''
	This is a temporary solution to the problem of generating moments.
	Basically, just set the first one equal to zero (since that rolls out of redshift)
	and set 2 - 4 equal to some reasonable values
	'''
	moments = [0., 45.]
	moments += [0. for _ in range(6 - len(moments))] #pad moments out to the length that the LOSVD function accepts

	c = 299792.458

	l_sparse = np.linspace(SPSPK_info['CRVAL1'], SPSPK_info['CRVAL1'] + SPSPK_info['NAXIS1'] * SPSPK_info['CDELT1'], SPSPK_info['NAXIS1'])
	v_res = np.mean(c / (l_sparse / FWHM_SPSPK))
	#print 'Instrument velocity resolution:', v_res
	
	generate_LOSVD(spectrum = clean_spectrum, v_res = v_res, moments = moments, plots = False)

	blurred_spectrum = adjust_FWHM(sharp_spectrum = clean_spectrum, res_old = MILES_res, res_new = SPSPK_res, FWHM_old = FWHM_MILES, FWHM_new = FWHM_SPSPK)

	#now redshift the new blurred (but still full-resolution) spectrum into the observer frame
	z, dz, l_full = redshift_spectrum(l_0 = l_full, dz = 0.)
	
	l_sparse, sparse_spectrum = downsample_spectrum(l_dense = l_full, dense_spectrum = blurred_spectrum, l_sparse = l_sparse) #this accomplishes both downsampling and paring!!

	#now construct a fake IFU with 8 rows of pure noise at some SNR
	ifu, galaxy_noise = noisify_ifu(sparse_spectrum, n = 2, SNR = SNR)

	#simulate_noise(sparse_spectrum, SNR = SNR)

	#more debugs
	'''plt.plot(l_sparse, sparse_spectrum, linewidth = 0.25, label = 'original')
	plt.plot(l_sparse, ifu[0], linewidth = 0.25, label = 'noisy')
	plt.plot(l_sparse, galaxy_noise, linewidth = 0.25, label = 'sample noise')
	plt.legend(loc = 'best')
	plt.show()
	'''
	
	edgefibers = range(1, len(ifu))
	pp = SPaCT.SP_pPXF(ifu/np.median(ifu[0]), fiber = 0, l_summ = (3907., 1.4, 1934), 
		z = z + dz, verbose = False, noise_plots = False, fit_plots = False, 
		edgefibers = edgefibers, age_lim = 20., n_moments = n_moments, bias = 100.)

	#now compare the resulting redshift
	print 'Best-fitting redshift:\t\t', z + dz + pp.sol[0]/c
	print 'Real redshift:\t\t\t\t', z
	print 'Guess redshift:\t\t\t\t', z + dz
	print 'Reduced chi2:', pp.chi2

	print ' #      |  guess  |  real'
	for (i, fit_guess, real_value) in zip(range(1, n_moments + 1), pp.sol, moments[:n_moments]):
		print 'moment', str(i), ':', str(np.round(fit_guess, 2)), ':', str(np.round(real_value, 2))

	#compare the resulting population fits
	#print pp.weights
	ssps.add_column(table.Column(name = 'best-fit weights', data = pp.weights/pp.weights.sum()))

	pPXF_summary_plots(ssps = ssps, instrument_info = SPSPK_info, pp = pp, lam_sparse = l_sparse, vel = (z + dz) * c + pp.sol[0], verbose = True)

	#now return the chi2 parameter for the best-fit, as opposed to the "reality"
	print 'Chi-square test'



simulate_single_spectrum()