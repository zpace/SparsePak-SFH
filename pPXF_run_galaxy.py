def pPXF_run_galaxy(objname, save = False):
	import SPaCT
	import ppxf
	from astroML.datasets import fetch_sdss_spectrum
	import warnings
	import matplotlib.pyplot as plt
	import astropy.io.fits as fits
	import astropy.io.ascii as ascii
	import astropy.table as table
	import numpy as np

	warnings.filterwarnings("ignore", message="using a non-integer number instead of an integer will result in an error in the future")

	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.rcParams['figure.figsize'] = 16, 12

	objects = {
		'NGC2558': {'z': 0.0166076, 'EBV': .0417, 'sdss_spec': (1927, 53321, 290)},
		'NGC2557': {'z': 0.0161257, 'EBV': .0366, 'sdss_spec': (1927, 53321, 387)},
		'NGC2562': {'z': 0.0168859, 'EBV': .0376, 'sdss_spec': (1927, 53321, 186)},
		'NGC2563': {'z': 0.0150407, 'EBV': .0381, 'sdss_spec': (1927, 53321, 194)},
		'NGC4061': {'z': 0.0244529, 'EBV': .0290, 'sdss_spec': (2608, 54474, 514)},
		'NGC6109': {'z': 0.0298147, 'EBV': .0196, 'sdss_spec': (1057, 52522, 391)},
		'NGC5846': {'z': None, 'EBV': .0477, 'sdss_spec': (None, None, None)},
		'UGC04344': {'z': 0.0167887, 'EBV': .0439, 'sdss_spec': (1927, 53321, 230)},
		'UGC04329': {'z': 0.0136599, 'EBV': .0440, 'sdss_spec': (4479, 55592, 204)}
	}

	#reddenings from IRSA (Schlafly & Finkbeiner 2011), using 5' radius from center

	z = objects[objname]['z']
	EBV = objects[objname]['EBV']
	dz = 0.

	fiberdata = ascii.read('fiberdata.dat')
	#print fiberdata

	im, fiberflat = SPaCT.load_ifus_precorrection(objname)

	plate = objects[objname]['sdss_spec'][0]
	mjd = objects[objname]['sdss_spec'][1]
	sdss_fiber = objects[objname]['sdss_spec'][2]
	sdss = fetch_sdss_spectrum(plate, mjd, sdss_fiber)

	dz, ifu_corr, corr_poly = SPaCT.sdss_cal(im, fiberflat, sdss, dz = 0., z = z, verbose = False, blur = 75, full_output = True)
	SPaCT.write_corr_frame(ifu_corr, im, z, dz, objname, verbose = False)

	edgefibers = [4, 51, 65, 10, 53, 1, 71, 8]
	fiberdata.add_column(table.Column(name = 'Z', data = np.nan*np.ones(len(fiberdata['row']))))
	fiberdata.add_column(table.Column(name = 't', data = np.nan*np.ones(len(fiberdata['row']))))
	fiberdata.add_column(table.Column(name = 'V', data = np.nan*np.ones(len(fiberdata['row']))))
	fiberdata.add_column(table.Column(name = 'sigma', data = np.nan*np.ones(len(fiberdata['row']))))

	for i, fiber_entry in enumerate(fiberdata):
		#print fiber_entry
		fiber = fiber_entry['row']
		if (fiber not in edgefibers) and ([fiber_entry['sky'] != 1]):
			print 'Running fiber', fiber

			n_moments = 2

			ifu = fits.open(objname + '_fluxcal.fits')[0].data
			try:
				pp, ssps = SPaCT.SP_pPXF((ifu.T/np.median(ifu, axis = 1)).T, fiber = fiber, l_summ = (3907., 1.4, 1934), 
					z = z + dz, verbose = False, noise_plots = False, fit_plots = False, clean = False, quiet = True, 
					edgefibers = edgefibers, age_lim = 13.5, n_moments = n_moments, bias = None, objname = objname, 
					oversample = False, reddening = EBV)

				fiber_Z = np.log10(np.average(10.**ssps['Z'], weights = ssps['fits']))
				fiber_age = np.average(ssps['t'], weights = ssps['fits'])
				fiber_V = pp.sol[0]
				fiber_sigma = pp.sol[1]

				fiberdata['Z'][np.where(fiberdata['row'] == fiber)] = fiber_Z
				fiberdata['t'][np.where(fiberdata['row'] == fiber)] = fiber_age
				fiberdata['V'][np.where(fiberdata['row'] == fiber)] = fiber_V
				fiberdata['sigma'][np.where(fiberdata['row'] == fiber)] = fiber_sigma
				#print fit_data[fiber]
				#print fit_data
			except:
				pass
			#if i == 3: break

	#fiberdata.pprint(max_lines = -1)

	SPaCT.gal_im_fiber_plot(objname = objname, fibers = fiberdata, quantity = 't', qty_dets = '[Gyr]', save = save)
	SPaCT.gal_im_fiber_plot(objname = objname, fibers = fiberdata, quantity = 'Z', qty_dets = '[M/H]', save = save)
	SPaCT.gal_im_fiber_plot(objname = objname, fibers = fiberdata, quantity = 'V', qty_dets = 'km/s', save = save)
	SPaCT.gal_im_fiber_plot(objname = objname, fibers = fiberdata, quantity = 'sigma', qty_dets = 'km/s', save = save)

	SPaCT.gal_rad_dep_plot(objname = objname, fibers = fiberdata, quantity = 't', qty_dets = '[Gyr]', save = save)
	SPaCT.gal_rad_dep_plot(objname = objname, fibers = fiberdata, quantity = 'Z', qty_dets = '[M/H]', save = save)
	SPaCT.gal_rad_dep_plot(objname = objname, fibers = fiberdata, quantity = 'Z', qty_dets = '[M/H]', save = save)
	SPaCT.gal_rad_dep_plot(objname = objname, fibers = fiberdata, quantity = 'Z', qty_dets = '[M/H]', save = save)