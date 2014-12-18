'''
SparsePak Correction Tools (SPaCT)
Author: Zach Pace (U Wisc-Madison)
License: GNU GPLv2

Some tools for correcting un-flux-calibrated SparsePak data.

Performance not guaranteed. Success dependent on correct redshift, correct centering of SparsePak.

Changelog:
Version 0.1 (Nov 2014) 
	- correct using interpolated fiberflat
Version 0.2 (Dec 2014) 
	- move from NGC2558 ipynb into standalone file
	- add correction with SDSS spectrum
	- add output to <OBJNAME>_fluxcal.fits
'''

import numpy as np
import astroML
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def interp_corr(ifu, fiberflat, verbose = False, full = False):
	'''
	Use fiberflat image to provide a small correction to the detector counts. Applied to all fibers at once.

	Arguments:
	 - ifu: .fits HDU for the full science IFU frame
	 - fiberflat: .fits HDU for the dome flat frame
	 - verbose: decides whether plots and supplementary information 
		are displayed, defaults to False (no plots)
	 - full: decides whether to output intermediate steps like polynomial fits

	Returns: [] denotes only for full = True mode
	 - ifu_interp_corr: interpolation-corrected science frame
	[- ifu_div_cor: division-corrected science frame]
	[- ifu_div_corr_m: division-correction (unsmoothed) matrix]
	[- ifu_interp_corr_m: interpolation-correction (smoothed) matrix]
	'''
	central_gal_fiber = 52
	central_det_fiber = int(np.floor(np.median(range(len(fiberflat[0].data)))))
	goodxllim = 200
	goodxhlim = np.shape(fiberflat[0].data)[1] - 100
	pc = polyfit(range(goodxllim, goodxhlim), fiberflat[0].data[central_det_fiber][goodxllim:goodxhlim], 4)

	if verbose == True:
		print np.shape(fiberflat[0].data)

		# display the un-corrected fiberflat
		fig = plt.figure(figsize=(10,20))
		from matplotlib.colors import LogNorm
		plt.imshow(fiberflat[0].data, origin = 'lower', aspect = 50, interpolation = 'None', cmap = 'binary')
		plt.axhline(central_gal_fiber+0.5, c ='r', linestyle = '--')
		plt.axhline(central_gal_fiber-0.5, c ='r', linestyle = '--')
		plt.text(250, central_gal_fiber+1, 'Bundle Center', size = 14, color = 'r')
		plt.axhline(central_det_fiber+0.5, c ='k', linestyle = '--')
		plt.axhline(central_det_fiber-0.5, c ='k', linestyle = '--')
		plt.text(250, central_det_fiber+1, 'Detector Center', size = 14, color = 'k')
		plt.axvline(goodxllim, c = 'b')
		plt.axvline(goodxhlim, c = 'b')
		plt.xlabel('detector position', size = 16)
		plt.ylabel('fiber num', size = 16)
		plt.colorbar(shrink = 0.8)
		plt.tight_layout()
		plt.title('Un-corrected fiberflat frame', size = 18)
		plt.show()

		# now overplot all the fiberflat spectra
		fig = plt.figure(figsize=(8,6))
		fibernums = range(0, len(fiberflat[0].data))
		for fiber in fibernums:
		    plt.plot(fiberflat[0].data[fiber], label = 'Fiber ' + str(fiber), linewidth = 0.5)
		plt.axvline(goodxllim, c = 'gray')
		plt.axvline(goodxhlim, c = 'gray')
		plt.ylim([-0.1, np.percentile(fiberflat[0].data, 99.9)])
		plt.title('fiberflat response, by fiber', size = 18)
		plt.xlabel('detector position', size = 14)
		plt.show()

		#finally, plot the response of the central fiber
		fig = plt.figure(figsize=(8,6))
		plt.plot(fiberflat[0].data[central_det_fiber], label = 'Detector counts')
		plt.plot(range(0, np.shape(fiberflat[0].data)[1]), np.polyval(pc, range(0, np.shape(fiberflat[0].data)[1])), label = 'Fit')
		plt.plot( np.abs( np.polyval(pc, range(0, np.shape(fiberflat[0].data)[1]))- fiberflat[0].data[central_det_fiber] ) , label = 'Residual' )
		plt.ylim([-.1, 1.5])
		plt.axvline(goodxllim, c = 'gray', label = 'lower lim')
		plt.axvline(goodxhlim, c = 'gray', label = 'upper lim')
		rowcenter = len(fiberflat[0].data[central_det_fiber])/2
		plt.axvline(rowcenter, c = 'orange', label = 'central pixel')
		plt.legend(loc = 'best')
		plt.xlabel('detector position', size = 14)
		plt.title('central fiber flat-field spectrum')
		plt.show()

	skyfibers = [16, 80, 2, 22, 70, 54, 37]

	# now  make a list of interp coefficients for each row
	row_interp_coeffs = []
	row_interp = []
	for row in fiberflat[0].data:
	    p = polyfit(range(goodxllim, goodxhlim), row[goodxllim:goodxhlim], 4)
	    row_interp_coeffs.append(p)
	    row_interp.append( np.polyval(p, np.arange(0., np.shape(fiberflat[0].data)[1], 1.) ) )

	row_mult = row_interp[central_det_fiber][rowcenter] / row_interp

	# neglect the sky fibers, and return 
	ifu_interp_corr = ifu * np.delete(row_mult, skyfibers, axis = 0)

	if full == True:
		ifu_interp_corr_m = np.delete(row_mult, skyfibers, axis = 0)
		ifu_div_corr_m = np.delete(np.asarray(fiberflat[0].data))
		ifu_div_cor = ifu * ifu_div_corr_m
		return ifu_interp_corr, ifu_div_corr, ifu_interp_corr_m, ifu_div_corr_m
	else:
		return ifu_interp_corr