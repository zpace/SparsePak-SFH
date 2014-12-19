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

Version 0.2.1 (Dec 2014)
	- add z_test function
'''

import numpy as np
import astroML
import matplotlib.pyplot as plt
import pyfits as fits
from astroML.datasets import fetch_sdss_spectrum

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def interp_corr(im, fiberflat, verbose = False, full = False):
	'''
	Use fiberflat image to provide a small correction to the detector counts. Applied to all fibers at once.

	Arguments:
	 - im: .fits HDU for the full science IFU frame
	 - fiberflat: .fits HDU for the dome flat frame
	 - verbose (False): decides whether plots and supplementary information 
		are displayed, defaults to False (no plots)
	 - full (False): decides whether to output intermediate steps like polynomial fits

	Returns: [] denotes only for full = True mode
	 - ifu_interp_corr: interpolation-corrected science frame
	[- ifu_div_cor: division-corrected science frame]
	[- ifu_div_corr_m: division-correction (unsmoothed) matrix]
	[- ifu_interp_corr_m: interpolation-correction (smoothed) matrix]
	'''

	# DEPENDENCIES
	import numpy as np
	import matplotlib.pyplot as plt
	import pyfits as fits

	# im contains the IFU data in layer [0].data
	ifu = im[0].data

	central_gal_fiber = 52
	central_det_fiber = int(np.floor(np.median(range(len(fiberflat[0].data)))))
	goodxllim = 200
	goodxhlim = np.shape(fiberflat[0].data)[1] - 100
	pc = np.polyfit(range(goodxllim, goodxhlim), fiberflat[0].data[central_det_fiber][goodxllim:goodxhlim], 4)
	rowcenter = len(fiberflat[0].data[central_det_fiber])/2.

	if verbose == True:
		print np.shape(fiberflat[0].data)

		# display the un-corrected fiberflat
		#fig = plt.figure(figsize=(10,20))
		from matplotlib.colors import LogNorm
		plt.imshow(fiberflat[0].data, aspect = 20, origin = 'lower', interpolation = 'None', cmap = 'binary')
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
		#plt.tight_layout()
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
	    p = np.polyfit(range(goodxllim, goodxhlim), row[goodxllim:goodxhlim], 4)
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

def z_test(im, fiberflat, sdss, z = 0.):
	'''
	display a plot of the interp-corrected and SDSS spectra, for the purposes of lining up lines and finding a redshift

	Arguments:
	 - im: .fits HDU for the full science IFU frame
	 - fiberflat: .fits HDU of the full fiberflat IFU frame
	 - sdss: SDSS spectrum object, as generated by astroML.datasets.fetch_sdss_spectrum
	 - z (0.): redshift (I recommend testing redshifts by hand before generating a correction)

	Returns:
	 - scaled SDSS spectrum

	'''

	ifu = im[0].data

	fig = plt.figure(figsize=(12,6), dpi=200)
	NAXIS1 = im[0].header['NAXIS1']
	CRVAL1 = im[0].header['CRVAL1']
	CDELT1 = im[0].header['CDELT1']
	print 'NAXIS1:', NAXIS1; print 'CRVAL1:', CRVAL1; print 'CDELT1:', CDELT1
	# native sparsepak wavelength range (no de-redshifting)
	sparsepak_wavelength = 1. / (1. + z) * np.linspace(CRVAL1, CRVAL1 + CDELT1 * NAXIS1, NAXIS1, endpoint = True )

	h = 6.62606957e-27 #erg sec
	c = 2.99792458e18 #Angstrom/sec
	t = 1200. #sec
	a = 9.6e4 #cm^2
	dl = CDELT1

	fluxcal = 1e17 * h*c/( (sparsepak_wavelength) * dl * t * a )
	
	#calculate and plot the theoretical flux-calibrated SparsePak spectrum for the central fiber
	sparsepak_spectrum = fluxcal * interp_corr(im, fiberflat)[47]

	plt.plot(sparsepak_wavelength, sparsepak_spectrum, linewidth = 0.5, label = 'Interpolation-corrected')

	# plot the SDSS spectrum
	plt.plot(sdss.wavelength(), sdss.spectrum/np.max(sdss.spectrum) * np.percentile(sparsepak_spectrum, 98.), linewidth = 0.5, color = 'grey', label = 'SDSS Spectrum', zorder = 0)

	plt.legend(loc = 'best', title = 'Central Fiber', fontsize = 10)
	plt.ylabel('$F_{\lambda}$ [$10^{-17} erg/s/\AA/cm^2$]', size = 16)
	plt.xlabel('$\lambda_{rest}$ [$\AA$]', size = 16)
	plt.xlim([3800., 6600.])
	plt.title('Test SED', size = 18)
	plt.show()

	return sdss.spectrum/np.max(sdss.spectrum) * np.percentile(sparsepak_spectrum, 98.)

def sdss_cal(im, fiberflat, sdss, z, verbose = False, fiber = 47):
	'''
	using a known z, calculate a correction to a SparsePak spectrum based on an SDSS spectrum of the same image.

	Arguments:
	 - im: .fits HDU for the full science IFU frame
	 - fiberflat: .fits HDU of the full fiberflat IFU frame
	 - sdss: SDSS spectrum object, as generated by astroML.datasets.fetch_sdss_spectrum
	 - z: redshift (I recommend testing redshifts by hand before generating a correction)
	 - verbose (False): display plots?
	 - fiber (47): plot which fiber?

	Returns:
	 - ifu_corr: corrected IFU science frame (in np array format)
	'''

	# load properties of un-corrected science frame
	NAXIS1 = im[0].header['NAXIS1']
	CRVAL1 = im[0].header['CRVAL1']
	CDELT1 = im[0].header['CDELT1']
	if verbose == True: print 'Un-corrected science frame'; print 'NAXIS1:', NAXIS1; print 'CRVAL1:', CRVAL1; print 'CDELT1:', CDELT1

	sparsepak_wavelength = 1. / (1. + z) * np.linspace(CRVAL1, CRVAL1 + CDELT1 * NAXIS1, NAXIS1, endpoint = True )

	h = 6.62606957e-27 #erg sec
	c = 2.99792458e18 #Angstrom/sec
	t = 1200. #sec
	a = 9.6e4 #cm^2
	dl = CDELT1

	fluxcal = 1e17 * h*c/( (sparsepak_wavelength) * dl * t * a )
	sparsepak_spectrum = fluxcal * interp_corr(im, fiberflat)[47]

	sdss_spectrum = sdss.spectrum/np.max(sdss.spectrum) * np.percentile(sparsepak_spectrum, 98.)
	sdss_wavelength = sdss.wavelength()

	sdss_spectrum_interp = np.interp(sparsepak_wavelength[sparsepak_wavelength < 6500.], sdss_wavelength, sdss_spectrum)
	corr_arr = sdss_spectrum_interp/sparsepak_spectrum[sparsepak_wavelength < 6500.]
	corr_poly = np.polyfit(sparsepak_wavelength[sparsepak_wavelength < 6500.], corr_arr, 40)

	if verbose == True:
		plt.plot(sparsepak_wavelength[sparsepak_wavelength < 6500.], corr_arr, color = 'r', linewidth = 0.25, label = 'Correction array')
		plt.plot(sparsepak_wavelength, np.polyval(corr_poly, sparsepak_wavelength), label = 'Correction polynomial')
		plt.axhline(1., linestyle = '--', color = 'k')
		plt.tick_params(axis='both', which='major', labelsize=16)
		plt.xlabel('wavelength [$\AA$]', size = 18)
		plt.ylabel('correction factor', size = 18)
		plt.title('SparsePak-SDSS correction', size = 20)
		plt.legend(loc = 'best')
		plt.show()

	ifu_corr = corr_arr * interp_corr(im, fiberflat)

	if verbose == True:
		plt.plot(sparsepak_wavelength, ifu_corr[fiber], color = 'r', linewidth = 0.5, label = 'Corrected SparsePak spectrum')
		plt.plot(sparsepak_wavelength, sparsepak_spectrum[fiber], color = 'grey', linewidth = 0.25, label = 'Original SparsePak spectrum')
		plt.plot(sparsepak_wavelength, interp_corr(im, fiberflat)[fiber], color = 'b', linewidth = 0.25, label = 'Flat-corrected SparsePak spectrum')
		plt.axhline(1., linestyle = '--', color = 'k')
		plt.tick_params(axis='both', which='major', labelsize=16)
		plt.xlabel('wavelength [$\AA$]', size = 18)
		plt.ylabel('$F_{\lambda} [10^{-17}erg/s/cm^{-2}/\AA]$', size = 18)
		plt.title('Corrected spectrum', size = 20)
		plt.legend(loc = 'best')
		plt.show()

def write_corr_frame(ifu_corr, im, z, objname):
	'''
	write out a new .fits file for a flux-calibrated SparsePak spectrum (corrected relative to SDSS spectrum)
	
	Arguments:
	 - ifu_corr: final science frame (np array: one row is one fiber)
	 - im: original science frame (.fits HDU), used for wavelength calibration
	 - z: correct z found by trial-and-error (must be the same as the one plugged into sdss_cal)
	 - objname: string name for object
	'''

im = fits.open('NGC2558.msobj.fits')
fiberflat = fits.open('NGC2558.fiberflat.fits')

plate = 1615
mjd = 53166
fiber = 513
sdss = fetch_sdss_spectrum(plate, mjd, fiber)

sdss_cal(im, fiberflat, sdss, .0095, verbose = True)