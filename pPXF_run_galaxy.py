'''
Run a single fiber of a single galaxy
'''

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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['figure.figsize'] = 16, 12

objname = 'NGC2558'
z = 0.0166076
dz = 0.

fiberdata = ascii.read('fiberdata.dat')
#print fiberdata

im, fiberflat = SPaCT.load_ifus_precorrection(objname)

plate = 1927
mjd = 53321
sdss_fiber = 290
sdss = fetch_sdss_spectrum(plate, mjd, sdss_fiber)

dz, ifu_corr, corr_poly = SPaCT.sdss_cal(im, fiberflat, sdss, dz = 0., z = z, verbose = False, blur = 75, full_output = True)
SPaCT.write_corr_frame(ifu_corr, im, z, dz, objname, verbose = False)

edgefibers = [4, 51, 65, 10, 53, 1, 71, 8]
fiberdata.add_column(table.Column(name = 'Z', data = np.nan*np.ones(len(fiberdata['row']))))
fiberdata.add_column(table.Column(name = 't', data = np.nan*np.ones(len(fiberdata['row']))))

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
				edgefibers = edgefibers, age_lim = 13.5, n_moments = n_moments, bias = None, objname = objname)

			fiber_Z = np.log10(np.average(10.**ssps['Z'], weights = ssps['fits']))
			fiber_age = np.average(ssps['t'], weights = ssps['fits'])

			fiberdata['Z'][np.where(fiberdata['row'] == fiber)] = fiber_Z
			fiberdata['t'][np.where(fiberdata['row'] == fiber)] = fiber_age
			#print fit_data[fiber]
			#print fit_data
		except:
			pass
		#if i == 3: break

#fiberdata.pprint(max_lines = -1)

#SPaCT.gal_im_fiber_plot(objname = objname, fibers = fiberdata, quantity = 't', qty_dets = '[Gyr]')
#SPaCT.gal_im_fiber_plot(objname = objname, fibers = fiberdata, quantity = 'Z', qty_dets = '[M/H]')

SPaCT.gal_rad_dep_plot(objname = objname, fibers = fiberdata, quantity = 't', qty_dets = '[Gyr]')
SPaCT.gal_rad_dep_plot(objname = objname, fibers = fiberdata, quantity = 'Z', qty_dets = '[M/H]')