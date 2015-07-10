'''
SparsePak Correction Tools (SPaCT)
Author: Zach Pace (U Wisc-Madison)
License: GNU GPLv2

Some tools for correcting un-flux-calibrated SparsePak data when there's
	no standard star.

Performance not guaranteed. Success dependent on correct redshift,
	correct centering of SparsePak & SDSS.

Changelog:
Version 0.1 (Nov 2014)
	- correct using interpolated fiberflat
Version 0.2 (Dec 2014)
	- move from NGC2558 ipynb into standalone file
	- add correction with SDSS spectrum
	- add output to <OBJNAME>_fluxcal.fits

Version 0.2.1 (Dec 2014)
	- add z_test function

Version 0.2.2 (Jan 2015)
	- pre-blur SDSS and SparsePak spectra before making correction array,
		to ensure just the overall shape of the curve is fit
	- pre-load reduced science and fiberflat frame

Version 0.3.0 (Jan 2015)
	- corrected to include both the survey-found z and a dz
	- added a blur parameter, along with a blur visualizer
	- provided SED in erg/s/cm^2/A, eliminating 10^-17 factor to make
		it more machine-readable

Version 0.4.0 (Mar 2015)
	- transcribed pPXF fit routine
	- added SNR calculations

Version 0.5.0 (Apr 2015)
	- corrected pPXF wrapper per M. Cappellari's suggestions
'''

# SYNTAX EXAMPLE

'''
im, fiberflat = SPaCT.load_ifus_precorrection('NGC2558')

plate = 1615
mjd = 53166
fiber = 513
sdss = SPaCT.fetch_sdss_spectrum(plate, mjd, fiber)

z, ifu_corr = SPaCT.sdss_cal(im, fiberflat, sdss, .0095, verbose = True)
SPaCT.write_corr_frame(ifu_corr, im, z, 'NGC2558')
'''


def plusminus(s):
    import numpy as np

    if s == 'p':
        return 1.
    if s == 'm':
        return -1.
    else:
        print 'bad input: plusminus takes either \'p\' or \'m\'!'


def flux_to_counts(flux, l, h=6.62606957e-27, c=2.99792458e18, t=1200.,
                   a=9.6e4, dl=1.4):
    import numpy as np
    '''
	Convert a flux measurement (in erg/s/cm^2/AA) into a photon counts
	measurement. Acts on a single data row at a time.
	'''

    fluxcal = h*c/((l) * dl * t * a)
    return flux / fluxcal


def counts_to_flux(counts, l, h=6.62606957e-27, c=2.99792458e18, t=1200.,
                   a=9.6e4, dl=1.4):
    import numpy as np
    '''
	Convert a counts measurement into a flux measurement (in erg/s/cm^2/AA).
	Acts on a single data row at a time.
	'''

    fluxcal = h*c/((l) * dl * t * a)
    return counts * fluxcal


def interp_corr(im, fiberflat, verbose=False, full=False):
    '''
    Use fiberflat image to provide a small correction to the detector counts.
    Applied to all fibers at once.

    Arguments:
     - im: .fits HDU for the full science IFU frame
     - fiberflat: .fits HDU for the dome flat frame
     - verbose (False): decides whether plots and supplementary information
            are displayed, defaults to False (no plots)
     - full (False): decides whether to output intermediate steps like
        polynomial fits

    Returns: [] denotes only for full = True mode
     - ifu_interp_corr: interpolation-corrected science frame
    [- ifu_div_cor: division-corrected science frame]
    [- ifu_div_corr_m: division-correction (unsmoothed) matrix]
    [- ifu_interp_corr_m: interpolation-correction (smoothed) matrix]
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits

    # im contains the IFU data in layer [0].data
    ifu = im[0].data

    central_gal_fiber = 52
    central_det_fiber = int(np.floor(np.median(range(len(fiberflat[0].data)))))
    goodxllim = 200
    goodxhlim = np.shape(fiberflat[0].data)[1] - 100
    pc = np.polyfit(range(goodxllim, goodxhlim), fiberflat[0].data[
                    central_det_fiber][goodxllim:goodxhlim], 4)
    rowcenter = len(fiberflat[0].data[central_det_fiber])/2.

    if verbose == True:
        print np.shape(fiberflat[0].data)

        # display the un-corrected fiberflat
        # fig = plt.figure(figsize=(10,20))
        from matplotlib.colors import LogNorm
        plt.imshow(fiberflat[0].data, aspect=20,
                   origin='lower', interpolation='None', cmap='binary')
        plt.axhline(central_gal_fiber+0.5, c='r', linestyle='--')
        plt.axhline(central_gal_fiber-0.5, c='r', linestyle='--')
        plt.text(250, central_gal_fiber+1, 'Bundle Center', size=14, color='r')
        plt.axhline(central_det_fiber+0.5, c='k', linestyle='--')
        plt.axhline(central_det_fiber-0.5, c='k', linestyle='--')
        plt.text(
            250, central_det_fiber+1, 'Detector Center', size=14, color='k')
        plt.axvline(goodxllim, c='b')
        plt.axvline(goodxhlim, c='b')
        plt.xlabel('detector position', size=16)
        plt.ylabel('fiber num', size=16)
        plt.colorbar(shrink=0.8)
        # plt.tight_layout()
        plt.title('Un-corrected fiberflat frame', size=18)
        plt.show()

        # now overplot all the fiberflat spectra
        fig = plt.figure(figsize=(8, 6))
        fibernums = range(0, len(fiberflat[0].data))
        for fiber in fibernums:
            plt.plot(
                fiberflat[0].data[fiber], label='Fiber ' + str(fiber),
                linewidth=0.5)
        plt.axvline(goodxllim, c='gray')
        plt.axvline(goodxhlim, c='gray')
        plt.ylim([-0.1, np.percentile(fiberflat[0].data, 99.9)])
        plt.title('fiberflat response, by fiber', size=18)
        plt.xlabel('detector position', size=14)
        plt.show()

        # finally, plot the response of the central fiber
        fig = plt.figure(figsize=(8, 6))
        plt.plot(fiberflat[0].data[central_det_fiber],
                 label='Detector counts')
        plt.plot(
            range(0, np.shape(fiberflat[0].data)[1]), np.polyval(
                pc, range(0, np.shape(fiberflat[0].data)[1])), label='Fit')
        plt.plot(
            np.abs(
                np.polyval(pc, range(
                    0, np.shape(fiberflat[0].data)[1])) -
                fiberflat[0].data[central_det_fiber]), label='Residual')
        plt.ylim([-.1, 1.5])
        plt.axvline(goodxllim, c='gray', label='lower lim')
        plt.axvline(goodxhlim, c='gray', label='upper lim')
        plt.axvline(rowcenter, c='orange', label='central pixel')
        plt.legend(loc='best')
        plt.xlabel('detector position', size=14)
        plt.title('central fiber flat-field spectrum')
        plt.show()

    skyfibers = [16, 80, 2, 22, 70, 54, 37]

    # now  make a list of interp coefficients for each row
    row_interp_coeffs = []
    row_interp = []
    for row in fiberflat[0].data:
        p = np.polyfit(
            range(goodxllim, goodxhlim), row[goodxllim:goodxhlim], 4)
        row_interp_coeffs.append(p)
        row_interp.append(
            np.polyval(p, np.arange(0., np.shape(fiberflat[0].data)[1], 1.)))

    row_mult = row_interp[central_det_fiber][rowcenter] / row_interp

    # neglect the sky fibers, and return
    ifu_interp_corr = ifu * np.delete(row_mult, skyfibers, axis=0)

    if full == True:
        ifu_interp_corr_m = np.delete(row_mult, skyfibers, axis=0)
        ifu_div_corr_m = np.delete(np.asarray(fiberflat[0].data))
        ifu_div_cor = ifu * ifu_div_corr_m
        return ifu_interp_corr, ifu_div_corr, ifu_interp_corr_m, ifu_div_corr_m
    else:
        return ifu_interp_corr


def z_test(im, fiberflat, sdss, fiber=47, dz=0.):
    '''
    display a plot of the interp-corrected and SDSS spectra,
        for the purposes of lining up lines and finding a redshift

    Arguments:
     - im: .fits HDU for the full science IFU frame
     - fiberflat: .fits HDU of the full fiberflat IFU frame
     - sdss: SDSS spectrum object, as generated by
        astroML.datasets.fetch_sdss_spectrum
     - dz (0.): redshift correction wrt SDSS (I recommend testing redshifts
        by hand before generating a correction)
     - fiber (47): which fiber to use

    Returns:
     - scaled SDSS spectrum
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits

    ifu = im[0].data

    fig = plt.figure(figsize=(12, 6), dpi=200)
    NAXIS1 = im[0].header['NAXIS1']
    CRVAL1 = im[0].header['CRVAL1']
    CDELT1 = im[0].header['CDELT1']
    print 'NAXIS1:', NAXIS1
    print 'CRVAL1:', CRVAL1
    print 'CDELT1:', CDELT1
    # native sparsepak wavelength range (no de-redshifting)
    sparsepak_wavelength = 1. / \
        (1. + dz) * \
        np.linspace(CRVAL1, CRVAL1 + CDELT1 * NAXIS1, NAXIS1, endpoint=True)

    h = 6.62606957e-27  # erg sec
    c = 2.99792458e18  # Angstrom/sec
    t = 1200.  # sec
    a = 9.6e4  # cm^2
    dl = CDELT1

    fluxcal = 1e17 * h*c/((sparsepak_wavelength) * dl * t * a)

    # calculate and plot the theoretical flux-calibrated SparsePak spectrum
    # for the central fiber
    sparsepak_spectrum = fluxcal * interp_corr(im, fiberflat)[fiber]
    plt.plot(sparsepak_wavelength, sparsepak_spectrum,
             linewidth=0.5, label='Interpolation-corrected')

    # plot the SDSS spectrum
    plt.plot(
        sdss.wavelength(), sdss.spectrum/np.max(sdss.spectrum) *
        np.percentile(sparsepak_spectrum, 98.), linewidth=0.5,
        color='grey', label='SDSS Spectrum', zorder=0)

    plt.legend(loc='best', title='Central Fiber', fontsize=10)
    plt.ylabel('$F_{\lambda}$ [$10^{-17} erg/s/\AA/cm^2$]', size=16)
    plt.xlabel('$\lambda_{rest}$ [$\AA$]', size=16)
    plt.xlim([3800., 6600.])
    title_text = 'Test SED, fiber ' + str(fiber)
    plt.title(title_text, size=18)
    plt.show()

    return sdss.spectrum/np.max(sdss.spectrum) * \
        np.percentile(sparsepak_spectrum, 98.)


def sdss_cal(im, fiberflat, sdss, dz, z=0, verbose=False, fiber=47,
             blur=100, full_output=False):
    '''
    using a known dz, calculate a correction to a SparsePak spectrum
        based on an SDSS spectrum of the same object.

    Arguments:
     - im: .fits HDU for the full science IFU frame
     - fiberflat: .fits HDU of the full fiberflat IFU frame
     - sdss: SDSS spectrum object, as generated by
        astroML.datasets.fetch_sdss_spectrum
     - dz: redshift offset from SDSS (I recommend testing redshifts
        by hand before generating a correction)
     - z: optional (but recommended) parameter that lets you check
        wavelength solution with known Balmer lines
     - verbose (False): display plots?
     - fiber (47): plot which fiber?
     - blur (50): size of gaussian blur
     - full_output (False): return the correction array?

    Returns:
     - dz: redshift offset (same as argument)
     - ifu_corr: corrected IFU science frame (in np array format)
     - corr_poly: polynomial-interpolated correction array
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import scipy as sp
    import colorpy.ciexyz as ciexyz
    import colorpy.colormodels as cmodels

    # load properties of un-corrected science frame
    NAXIS1 = im[0].header['NAXIS1']
    CRVAL1 = im[0].header['CRVAL1']
    CDELT1 = im[0].header['CDELT1']
    if verbose == True:
        print 'Un-corrected science frame'
        print 'NAXIS1:', NAXIS1
        print 'CRVAL1:', CRVAL1
        print 'CDELT1:', CDELT1

    sparsepak_wavelength = 1. / \
        (1. + dz) * \
        np.linspace(CRVAL1, CRVAL1 + CDELT1 * NAXIS1, NAXIS1, endpoint=True)

    h = 6.62606957e-27  # erg sec
    c = 2.99792458e18  # Angstrom/sec
    t = 1200.  # sec
    a = 9.6e4  # cm^2
    dl = CDELT1

    # theoretical flux-calibration for the observed counts
    fluxcal = h*c/((sparsepak_wavelength) * dl * t * a)
    sparsepak_spectrum = counts_to_flux(
        interp_corr(im, fiberflat), sparsepak_wavelength)
    sparsepak_spectrum_ctr = sparsepak_spectrum[fiber]

    sdss_spectrum = sdss.spectrum / \
        np.max(sdss.spectrum) * np.percentile(sparsepak_spectrum, 98.)
    sdss_wavelength = sdss.wavelength()

    sdss_spectrum_interp = np.interp(
        sparsepak_wavelength[sparsepak_wavelength < 6500.],
        sdss_wavelength, sdss_spectrum)
    # smooth both the resampled sdss spectrum and the sparsepak spectrum, and
    # divide to get a correction array
    corr_arr = sp.ndimage.filters.gaussian_filter1d(
        sdss_spectrum_interp, blur)/sp.ndimage.filters.gaussian_filter1d(
        sparsepak_spectrum_ctr[sparsepak_wavelength < 6500.], blur)
    corr_poly = np.polyfit(
        sparsepak_wavelength[sparsepak_wavelength < 6500.], corr_arr, 40)

    if verbose:
        # plot the smoothed SDSS spectrum and the smoothed SparsePak spectrum
        plt.figure(figsize=(8, 6))
        ax1 = plt.subplot(111)
        ax1.plot(sdss.wavelength(), sdss.spectrum,
                 label='SDSS raw', linewidth=0.25, c='b')
        ax1.plot(sdss.wavelength(), sp.ndimage.filters.gaussian_filter1d(
            sdss.spectrum, blur), c='b', linewidth=2, label='SDSS shape')
        ax1.set_xlabel('$\lambda [\AA]$', size=14)
        ax1.set_ylabel('SDSS $F_{\lambda}$', color='b', size=14)
        [t.set_color('blue') for t in ax1.yaxis.get_ticklines()]
        [t.set_color('blue') for t in ax1.yaxis.get_ticklabels()]

        ax2 = plt.twinx(ax1)
        ax2.plot(
            sparsepak_wavelength[sparsepak_wavelength < 6500.],
            sparsepak_spectrum_ctr[sparsepak_wavelength < 6500.], c='k',
            linewidth=0.25, label='SparsePak raw')
        ax2.plot(sparsepak_wavelength[sparsepak_wavelength < 6500.],
                 sp.ndimage.filters.gaussian_filter1d(
            sparsepak_spectrum_ctr[sparsepak_wavelength < 6500.], blur),
            c='k', linewidth=2, label='SparsePak shape')
        ax2.set_ylabel('SparsePak $F_{\lambda}$', size=14, color='k')
        [t.set_color('k') for t in ax1.yaxis.get_ticklines()]
        [t.set_color('k') for t in ax1.yaxis.get_ticklabels()]

        ax1.legend(loc='upper left')
        ax2.legend(loc='upper right')
        plt.show()

    if verbose:
        # plot fit and error
        ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)  # axis for fit plot
        ax2 = plt.subplot2grid((2, 2), (1, 0))  # axis for fit error plot
        ax3 = plt.subplot2grid((2, 2), (1, 1))  # axis for fit error histogram

        ax1.plot(
            sparsepak_wavelength,
            np.polyval(corr_poly, sparsepak_wavelength), zorder=0,
            label='Correction polynomial')
        ax1.plot(
            sparsepak_wavelength[sparsepak_wavelength < 6500.], corr_arr,
            color='r', linewidth=0.25, zorder=1, label='Correction array')
        ax1.axhline(1., linestyle='--', color='k')
        ax1.tick_params(axis='both', which='major', labelsize=16)
        ax1.set_xlabel('wavelength [$\AA$]', size=18)
        ax1.set_ylabel('correction factor', size=18)
        ax1.set_title('SparsePak-SDSS correction', size=20)
        ax1.legend(loc='best')

        ax2.plot(
            sparsepak_wavelength[sparsepak_wavelength < 6500.],
            np.abs(np.polyval(
                corr_poly,
                sparsepak_wavelength[sparsepak_wavelength < 6500.]) -
                corr_arr) / corr_arr, linewidth=0.5)
        ax2.set_xlabel('wavelength [$\AA$]', size=14)
        ax2.set_ylabel(
            'Correction function relative error [$\\frac{F_c - A_c}{A_c}$]',
            size=14)

        ax3.hist(
            np.log10(
                np.abs(
                    np.polyval(
                        corr_poly,
                        sparsepak_wavelength[sparsepak_wavelength < 6500.]) -
                    corr_arr) / corr_arr), bins=20)
        ax3.set_xlabel('$\log\\frac{F_c - A_c}{A_c}$', size=14)
        print 'Median log-error:', np.median(
            np.log10(
                np.abs(
                    np.polyval(
                        corr_poly,
                        sparsepak_wavelength[sparsepak_wavelength < 6500.]) -
                    corr_arr) / corr_arr))

        # for diagnostic purposes, plot Balmer series
        plt.tight_layout()
        plt.show()

    if verbose:
        # plot SDSS spectrum and new corrected spectrum
        ax4 = plt.subplot(111)
        ax4.plot(
            sparsepak_wavelength,
            np.polyval(corr_poly, sparsepak_wavelength) * fluxcal *
            interp_corr(im, fiberflat)[fiber], c='b', linewidth=1.,
            label='SDSS-flux-calibrated SparsePak spectrum')
        ax4.plot(
            sdss_wavelength[sdss_wavelength < 6500.],
            sdss_spectrum[sdss_wavelength < 6500.], linewidth=0.25,
            c='r', label='SDSS spectrum')
        ax4.tick_params(axis='both', which='major', labelsize=16)
        ax4.set_xlabel('wavelength [$\AA$]', size=18)
        ax4.set_ylabel('$F_{\lambda} [erg/s/cm^{-2}/\AA]$', size=18)
        ax4.set_title('Corrected spectrum', size=20)
        ax4.legend(loc='best')
        ax4.set_ylim([0., np.max(sdss_spectrum)*1.1])

        if z:
            balmer = [['$\\alpha$', 6563.], ['$\\beta$', 4861.],
                      ['$\\gamma$', 4341.], ['$\\delta$', 4102.]]
            for line in balmer:
                line_c = cmodels.irgb_string_from_xyz(
                    ciexyz.xyz_from_wavelength(line[1]/10.))
                ax4.axvline(line[1] * (1 + dz + z), color=line_c)
                ax4.annotate(
                    'H-' + line[0],
                    xy=(line[1], 0.9*ax4.get_ylim()[1]),
                    xytext=(line[1]+10., 0.95*ax4.get_ylim()[1]), size=14)
            metals = [['Mg[I]', 5175.], ['Ca[I]', 4307.],
                      ['Ca H', 3968.], ['Ca K', 3933.]]
            for line in metals:
                line_c = cmodels.irgb_string_from_xyz(
                    ciexyz.xyz_from_wavelength(line[1]/10.))
                ax4.axvline(line[1] * (1 + dz + z), color=line_c)
                ax4.annotate(
                    line[0], xy=(line[1], 0.9*ax4.get_ylim()[1]),
                    xytext=(line[1]+10., 0.95*ax4.get_ylim()[1]), size=14)
        plt.show()

    ifu_corr = sparsepak_spectrum * np.polyval(corr_poly, sparsepak_wavelength)

    if full_output == True:
        return dz, ifu_corr, np.polyval(corr_poly, sparsepak_wavelength)
    else:
        return dz, ifu_corr


def write_corr_frame(ifu_corr, im, z, dz, objname, verbose=False):
    '''
    write out a new .fits file for a flux-calibrated SparsePak spectrum
        (corrected relative to SDSS spectrum)

    Arguments:
     - ifu_corr: final science frame (np array: one row is one fiber)
     - im: original science frame (.fits HDU), used for wavelength calibration
     - z: accepted (SDSS) redshift
     - dz: correct z offset found by trial-and-error
        (must be the same as the one plugged into sdss_cal)
     - objname: string name for object
     - verbose: display basic diagnostic information about file to be written

    Returns: NONE
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits

    NAXIS1 = im[0].header['NAXIS1']
    CRVAL1 = im[0].header['CRVAL1']
    CDELT1 = im[0].header['CDELT1']
    if verbose == True:
        print 'Un-corrected science frame'
        print 'NAXIS1:', NAXIS1
        print 'CRVAL1:', CRVAL1
        print 'CDELT1:', CDELT1

    new_header = fits.Header()
    new_header.set('SIMPLE', 'T', 'conforms to FITS standard')
    new_header.set('BITPIX', -64, 'array data type')
    new_header.set('NAXIS', 2, 'number of array dimensions')
    new_header.set('EXTEND', 'T')
    new_header.set('NAXIS1', NAXIS1)
    new_header.set('NAXIS2', 75)
    new_header.set('CRVAL1', CRVAL1)
    new_header.set('CDELT1', CDELT1)
    new_header.set('BUNIT', 'Data Value')
    new_header.set('Z', z + dz)

    sparsepak_wavelength = np.linspace(
        CRVAL1, CRVAL1 + CDELT1 * NAXIS1, NAXIS1, endpoint=True)

    CRVAL1 = np.min(sparsepak_wavelength)
    CDELT1 = 1.4
    NAXIS1 = len(ifu_corr[47])

    hdu_flux_cal = fits.PrimaryHDU(data=ifu_corr, header=new_header)
    hdulist = fits.HDUList([hdu_flux_cal])
    hdulist.writeto(objname + '_fluxcal.fits', clobber=True)

    if verbose == True:
        print fits.open(objname + '_fluxcal.fits')[0].header


def load_ifus_precorrection(objname):
    '''
    Load in the science frame (un-corrected) and the fiberflat

    Arguments:
     - objname: string name for object (i.e., the root filename,
        onto which .fiberflat.fits and .fits are appended)

    Returns:
     - im: .fits HDU of original, un-corrected (but reduced) science frame
     - fiberflat: .fits HDU of fiberflat frame
        (obtained through reduction pipeline)
    '''

    from astropy.io import fits

    return fits.open(
        objname + '.msobj.fits'), fits.open(objname + '.fiberflat.fits')


def noise_edgefibers(ifu, width, fiberlist, verbose=False):
    '''
        Find the noise in a sparsepak spectrum by median-combining edge
                fibers, then taking their statistics

        Make sure to feed it a FULL IFU ARRAY in flux units.
                It will output noise in flux units, as well.

        This will work if fibers are similar in response and there's
                not any appreciable non-uniformity in sky subtraction

        NOW DEPRECATED
    '''

    import pandas as pd
    import numpy as np

    if fiberlist == None:
        fiberlist = [51, 4, 1, 8, 10, 65, 71, 53]
    # print fiberlist

    if verbose == True:
        print 'width:', width
    # first select th e correct fibers from ifu and median-combine them
    ifu = np.median(ifu[fiberlist], axis=0)
    # now slide an aperture of some width over the spectrum
    # and take the stdev of the data in that window
    noise_array_1 = pd.rolling_std(
        ifu, window=width, center=True, min_periods=1)
    noise_array_2 = pd.rolling_std(
        ifu[::-1], window=width, center=True, min_periods=1)
    # requires a little reshaping
    noise = np.nanmean(
        np.concatenate((noise_array_1, noise_array_2[::-1])).reshape(
            (2, len(noise_array_1))), axis=0)
    return noise


def lin_remap(old, oldrange, newrange):
    '''
    linearly remaps an array to a different range
    '''
    import numpy as np

    newmin, newmax = newrange
    oldmin, oldmax = oldrange
    new = ((old - oldmin) / (oldmax - oldmin)) * (newmax - newmin) + newmin
    return new


def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
    Input data
    s : scalar or array_like, shape (n, )
            Radius of circle in data scale (ie. in data unit)
    c : color or sequence of color, optional, default : 'b'
            `c` can be a single color format string, or a sequence of color
            specifications of length `N`, or a sequence of `N` numbers to be
            mapped to colors using the `cmap` and `norm` specified via kwargs.
            Note that `c` should not be a single numeric RGB or
            RGBA sequence because that is indistinguishable from an array of
            values to be colormapped.  `c` can be a 2-D array in which the
            rows are RGB or RGBA, however.
    ax : Axes object, optional, default: None
            Parent axes of the plot. It uses gca() if not specified.
            vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with `norm` to normalize
            luminance data.  If either are `None`, the min and max of the
            color array is used.  (Note if you pass a `norm` instance, your
            settings for `vmin` and `vmax` will be ignored.)

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Other parameters
    ----------------
    kwargs : `~matplotlib.collections.Collection` properties
            eg. alpha, edgecolors, facecolors, linewidths,
            linestyles, norm, cmap

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    """
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection
    import pylab as plt
    import numpy as np
    # import matplotlib.colors as colors

    if ax is None:
        ax = plt.gca()

    if isinstance(c, basestring):
        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
    else:
        color = None  # use cmap, norm after collection is created
    kwargs.update(color=color, edgecolor=None)

    if isinstance(x, (int, long, float)):
        patches = [Circle((x, y), s), ]
    elif isinstance(s, (int, long, float)):
        patches = [Circle((x_, y_), s) for x_, y_ in zip(x, y)]
    else:
        patches = [Circle((x_, y_), s_) for x_, y_, s_ in zip(x, y, s)]
    collection = PatchCollection(patches, **kwargs)

    if color is None:
        collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    return collection


def gal_im_fiber_plot(objname, fibers, quantity=None, qty_dets='',
                      fibersize=4.687, text=False, save=False, offset=0.,
                      oe='even'):
    '''
    fibers is an astropy table mostly in format of `fiberdata.dat`,
    with an additional row based on what is being plotted
    '''
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib as mpl
    import plotting_tools

    plt.close('all')

    fig = plt.figure(figsize=(4, 4), dpi=300)

    # you can get custom-sized images at http://skyserver.sdss3.org/dr10/en/tools/chart/image.aspx
    # explore at http://skyserver.sdss3.org/dr10/en/tools/explore/default.aspx

    im = mpimg.imread(objname + '.png')
    x = y = np.linspace(-40., 40., np.shape(im)[0])
    X, Y = np.meshgrid(x, y)
    galim = plt.imshow(
        im, extent=[-40, 40, -40, 40], origin='lower',
        interpolation='nearest', zorder=0, aspect='equal')

    fibers = fibers[fibers['sky'] != 1]
    fibers = fibers[np.isnan(fibers[quantity]) == False]
    fibers[quantity] += offset
    # print fibers

    # set up cmap based on whether the quantity is even or odd
    # if odd, then it enforces symmetry around `offset` values
    if oe == 'odd':
        cmap = mpl.cm.get_cmap('RdBu_r')
    else:
        cmap = 'cubehelix_r'

    if (quantity != None) and (quantity not in ['Z', 't']):
        if quantity == 'V':
            cir = circles(fibers['ra'], fibers['dec'], s=fibersize/2.,
                          c=fibers[quantity], alpha=0.8, zorder=2,
                          cmap=cmap, vmin=-300., vmax=300.)
        elif quantity == 'sigma':
            cir = circles(fibers['ra'], fibers['dec'], s=fibersize/2.,
                          c=fibers[quantity], alpha=0.8, zorder=2,
                          cmap=cmap, vmin=0.0, vmax=300.)
        else:
            cir = circles(fibers['ra'], fibers['dec'], s=fibersize/2.,
                          c=fibers[quantity], alpha=0.8, zorder=2,
                          cmap=cmap)
        cbar = fig.colorbar(cir, shrink=0.8)
        cbar.set_label(quantity + qty_dets, size=12)
        if quantity == 'V':
            cbar.cmap.set_over('r')
            cbar.cmap.set_under('b')
        elif oe == 'even':
            cbar.cmap.set_over('k')
            cbar.cmap.set_under('w')

    elif quantity in ['Z', 't']:  # since these are percentile bounds
        q = fibers[quantity]
        cir = circles(fibers['ra'], fibers['dec'], s=fibersize/2.,
                      c=q, alpha=0.8, zorder=2, cmap=cmap)
        cbar = fig.colorbar(cir, shrink=0.8)
        cbar.set_label(quantity + qty_dets, size=12)
    else:
        cir = circles(fibers['ra'], fibers['dec'], s=fibersize/2.,
                      edgecolor='None', alpha=0.8, zorder=2)

    if text == True:
        for row in fibers:
            if row['sky'] != 1:
                plt.text(
                    row['ra']-1., row['dec']-1.75,
                    s=str(row['fiber']) + '\n' + str(row['row']),
                    color='g', size=14, zorder=3)

    plt.xlim([-40, 40])
    plt.xticks([-40., -20., 0., 20., 40.])
    plt.ylim([-40, 40])
    plt.yticks([-40., -20., 0., 20., 40.])
    plt.xlabel('RA offset [arcsec]', fontsize=12)
    plt.ylabel('Dec offset [arcsec]', fontsize=12)
    plt.title(objname, fontsize=14)
    plt.tight_layout()

    if (save == True) and (quantity != None):
        try:
            os.makedirs(objname)
        except OSError:
            pass
        plt.savefig(objname + '/' + quantity + 'map.png')
    else:
        plt.show()


def gal_rad_dep_plot(objname, fibers, quantity=None, qty_dets='',
                     save=False, offset=0., fit=False):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import fitting_tools

    plt.close('all')
    plt.figure(figsize=(4, 4), dpi=200)

    fibers = fibers[fibers['sky'] != 1]
    fibers = fibers[np.isnan(fibers[quantity]) == False]
    fibers[quantity] += offset

    if (quantity == 'V') or (quantity == 'sigma'):
        plt.errorbar(fibers['r'], fibers[quantity],
                     yerr=fibers['d' + quantity], fmt='x', alpha=0.5)
        yllim = np.min(np.insert(np.array(fibers[quantity]), 0, -1000.))
        yulim = np.max(np.insert(np.array(fibers[quantity]), 0, 1000.))
        plt.ylim([0.8*yllim, 1.2*yulim])
    elif (quantity == 'Z') or (quantity == 't'):
        plt.scatter(fibers['r'], fibers[quantity],
                    marker='x', alpha=0.5)
    else:
        print 'ERROR: UNKNOWN QUANTITY... PLOT WILL BE BLANK'

    plt.xlim([-5., 50.])

    if fit == True:
        if (quantity == 'V') or (quantity == 'sigma'):
            dq = fibers['d' + quantity]
            dq[dq == 0.] = np.min(dq[dq != 0.])
        else:
            dq = None

        '''
        popt, pcov = fitting_tools.linear_fit(
            fibers['r'], fibers[quantity], sigma=dq)

        m, b = popt[0], popt[1]
        if type(pcov) == np.ndarray:
            dm, db = np.sqrt(np.diagonal(pcov)[0]), np.sqrt(
                np.diagonal(pcov)[1])

            x_display = np.linspace(fibers['r'].min(), fibers['r'].max(), 10)
            y_display = fitting_tools.linear(x_display, m, b)
            line_label = r'$m={0:.3f} \pm {1:.3f}$; $b={2:.2f} \pm {3:.2f}$'\
                .format(m, dm, b, db)
            plt.plot(x_display, y_display, linestyle='--', c='k', alpha=0.75,
                     label=line_label)
            plt.legend(loc='best', prop={'size': 8})
        '''

    plt.ylabel(quantity + qty_dets, size=16)
    plt.xlabel('radius [arcsec]', size=16)
    plt.title(objname, fontsize=18)
    plt.tight_layout()

    if (save == True) and (quantity != None):
        try:
            os.makedirs(objname)
        except OSError:
            pass

        plt.savefig(objname + '/' + quantity + 'grad.png')
    else:
        plt.show()


def rebin(x, factor):
    """
    Rebin a one-dimensional vector by averaging
    in groups of "factor" adjacent values

    From Cappellari's ppxf_simulation_example.py

    """
    return np.mean(x.reshape(-1, factor), axis=1)


def pPXF_MC(pp, lam):
    '''
    INCOMPLETE

    A semi-optional step, post-fitting, to simulate pPXF result for a
    given fit, to better estimate kinematics errors

    Fashioned after Cappellari et al's ppxf_simulation_example.py

    Arguments:
     - pp: ppxf object
     - lam: array of wavelengths (same as fed into pPXF earlier)

    '''

    import numpy as np
    import ppxf_util as util
    from scipy import ndimage, signal
    import astropy.io.fits as fits

    bestfit = pp.bestfit


def SP_pPXF(ifu, fiber, l_summ, z, template_set='MILES', verbose=False,
            noise_plots=False, fit_plots=False, reddening=None, age_lim=13.6,
            n_moments=4, bias=None,	objname='', clean=True, quiet=False,
            oversample=False, save_fits=False, gas_comps=None, regul=100.):
    '''
    Run Cappellari et al.'s pPXF on a SparsePak fiber

    Arguments:
     - ifu: full IFU science array, put through all the reduction &
        correction ordeals
     - fiber: which fiber data row to use
     - l_summ: tuple of (CRVAL1, CDELT1, NAXIS1)
     - z: approximate starting redshift (using the quoted
        SDSS value is pretty safe)
     - template_set ('MILES'): which set of templates to use
        (either MILES or Vazdekis--MILES recommended)
     - noise_plots (False): plot noise characteristics of the fiber?
     - verbose (False): provide updates & diagnostics on fit
     - fit_plots (False): plot total fit & population diagram?
     - reddening (None): E(B-V) estimate
     - gas_comps (None): how many gas components to use (if it's
        `None`, then pPXF just fits the stars; the first gas component is
        **always** bound to the stellar velocity field)
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import signaltonoise
    from astropy.io import fits

    import sys
    import ppxf_util as util
    from ppxf import ppxf
    import glob as glob
    from scipy import ndimage
    from time import clock

    from astroML.datasets import fetch_sdss_spectrum
    import scipy as sp
    import colorpy.ciexyz as ciexyz
    import colorpy.colormodels as cmodels
    import warnings

    import os

    CRVAL1, CDELT1, NAXIS1 = l_summ
    lamRange1 = np.array([CRVAL1, CRVAL1 + CDELT1 * (NAXIS1 - 1)])
    FWHM_gal = 4.877  # FWHM from SparsePak simulator

    lamRange1 = lamRange1
    FWHM_gal = FWHM_gal / (1 + z)

    # print z

    galaxy, logLam1, velscale = util.log_rebin(lamRange1, ifu[fiber])
    lam_obs = np.exp(logLam1)
    logLam1 -= z
    lamRange1 = np.array([np.exp(logLam1[0]), np.exp(logLam1[-1])])

    # keep track of which pixels are negative, to mask them out later
    negative_pixels = np.where(galaxy < 0.)[0]
    galaxy = np.abs(galaxy)

    if verbose == True:
        print velscale[0], 'km/s/pix'

    # which template set is being used? make some decisions based on that
    if template_set == 'MILES':
        template_files = glob.glob('miles_models/Mun*.fits')
        template_files = [
            i for i in template_files if float(i[-12:-5]) < age_lim]
        FWHM_tem = 1.5  # Miles spectra have FWHM of 1.5A
    elif template_set == 'vazdekis':
        template_files = glob.glob(directory + 'Rbi1.30z*.fits')
        FWHM_tem = 1.8  # Vazdekis spectra have a resolution FWHM of 1.8A.
        galaxy, logLam1 = galaxy[:10], logLam1[:-10]
    else:
        raise NameError(
            'Invalid template set. Only \'vazdekis\' and \'MILES\' \
            are valid at this time.')

    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SparsePak galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.

    hdu = fits.open(template_files[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0., h2['CDELT1']*(h2['NAXIS1']-1)])
    sspNew, logLam2, velscale = util.log_rebin(
        lamRange2, ssp, velscale=velscale)
    templates = np.empty((sspNew.size, len(template_files)))

    if verbose == True:
        print 'SparsePak range:', lamRange1
        print 'Template range:', lamRange2

    # Pare down the SparsePak spectrum so that its wavelength range is
    # contained within the template spectrum's wavelength range

    galaxy = galaxy[(logLam1 >= logLam2[0]) * (logLam1 <= logLam2[-1])]
    logLam1 = logLam1[(logLam1 >= logLam2[0]) * (logLam1 <= logLam2[-1])]
    # the MILES templates run into issues at the red end of the spectra b/c of
    # correction in interpolation

    l = np.exp(logLam1)

    # munge the templates to get characteristics
    ssp_rows = []
    from astropy.table import Table, Column

    for template in template_files:
        template = template.rstrip('.fits').split('/')[1]

        spectral_range = template[0]
        IMF_type = template[1:3]
        IMF_slope = float(template[3:7])
        Z = plusminus(template[8])*float(template[9:13])
        T = float(template[14:])

        # print template + ':', spectral_range, IMF_type, IMF_slope, Z, T
        ssp_i = [template, spectral_range, IMF_type, IMF_slope, Z, T]
        ssp_rows.append(ssp_i)

    ssps = Table(
        map(
            list, zip(*ssp_rows)),
        names=['template name', 'spectral range', 'IMF type',
               'IMF slope', 'Z', 't'])
    ssps.sort(['Z', 't'])

    template_prefix = template_files[0].split('/')[0]
    # print template_prefix
    template_files = ssps['template name']

    if verbose == True:
        print ssps
        plt.figure(figsize=(4, 6))
        plt.scatter(ssps['Z'], ssps['t'])
        plt.xlabel('Z', size=16)
        plt.ylabel('age', size=16)
        plt.title('SSPs used', size=16)
        plt.show()

    # Convolve the whole SSP library of spectral templates
    # with the quadratic difference between the SAURON and the
    # Vazdekis instrumental resolution. Logarithmically rebin
    # and store each template as a column in the array TEMPLATES.

    # Quadratic sigma difference in pixels SSP --> instrument FWHM
    # The formula below is rigorously valid if the shapes of the
    # instrumental spectral profiles are well approximated by Gaussians.

    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
    sigma = FWHM_dif/2.355/h2['CDELT1']  # Sigma difference in pixels

    for j in range(len(template_files)):
        hdu = fits.open(template_prefix + '/' + template_files[j] + '.fits')
        ssp = hdu[0].data
        ssp = ndimage.gaussian_filter1d(ssp, sigma)
        sspNew, logLam2, velscale = util.log_rebin(
            lamRange2, ssp, velscale=velscale)
        templates[:, j] = sspNew

    templates /= np.median(templates)  # Normalizes templates

    # The galaxy and the template spectra do not have the same starting wavelength.
    # For this reason an extra velocity shift DV has to be applied to the template
    # to fit the galaxy spectrum. We remove this artificial shift by using the
    # keyword VSYST in the call to PPXF below, so that all velocities are
    # measured with respect to DV. This assume the redshift is negligible.
    # In the case of a high-redshift galaxy one should de-redshift its
    # wavelength to the rest frame before using the line below (see above).

    c = 299792.458
    dl = np.exp(logLam2[0])-np.exp(logLam1[0])
    dv = c*(logLam2[0] - logLam1[0])  # km/s
    # print dl, dv
    # print 'dv:', dv, 'km/s'

    vel = 0.
    # print 'SDSS velocity:', c*z, 'km/s'
    goodPixels = util.determine_goodpixels(logLam1, lamRange2, vel)
    # mask out the OI sky line and whatever weirdness goes along with it
    # goodPixels = np.array(
    #    [pixel for pixel in goodPixels if not 1183 <= pixel <= 1200])
    goodPixels = np.array(
        [pixel for pixel in goodPixels
         if not 5563.2 <= lam_obs[pixel] <= 5587.])
    # mask out last 100 pixels with weird end behavior
    goodPixels = np.array(
        [pixel for pixel in goodPixels if not len(galaxy) - 100
         <= pixel <= len(galaxy)])
    # mask out pixels that were originally negative
    goodPixels = np.array(
        [pixel for pixel in goodPixels if pixel not in negative_pixels])

    # Here the actual fit starts. The best fit is plotted on the screen.
    # Gas emission lines are excluded from the pPXF fit using the GOODPIXELS
    # keyword.

    start = [vel, 2.5*velscale]  # (km/s), starting guess for [V,sigma]

    t = clock()

    log_ages = np.log10(ssps['t'])
    d_log_ages = log_ages[1] - log_ages[0]
    dz = np.abs(np.mean(ssps['Z'][1:] - ssps['Z'][:-1]))

    # now reshape the templates according to the number of Z and t values
    # available
    nZ, nt, nalpha = len(np.unique(ssps['Z'])), len(
        np.unique(ssps['t'])), len(np.unique(ssps['IMF slope']))
    templates = np.reshape(templates, (len(templates[:, 0]), nZ, nt, nalpha))
    # print np.shape(templates)

    reg_dim = templates.shape[1:]
    if reg_dim[-1] == 1.:
        reg_dim = reg_dim[:-1]
    # print reg_dim

    templates = templates.reshape(templates.shape[0], -1)
    templates /= np.median(templates)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)

        lam = np.exp(logLam1)

        print 'ICs:', start

        # first pass: initial fit of parameters, no noise, no regul, no gas

        pp = ppxf(templates=templates, galaxy=galaxy/np.median(galaxy),
                  noise=galaxy*0.+1., velScale=velscale, start=start,
                  goodpixels=goodPixels, plot=False, moments=n_moments,
                  vsyst=dv, reddening=reddening,
                  lam=lam, quiet=True)
        # print 'First pass successful'
        print 'First pass results'
        print '\t', pp.sol
        # this will be anomalous, since uniform noise spectrum is assumed
        print '\tChi2/DOF =', pp.chi2

        # new noise spectrum will be residual from 1st fit
        noise = np.abs(pp.bestfit - galaxy/np.median(galaxy))

        pp = ppxf(templates=templates, galaxy=galaxy/np.median(galaxy),
                  noise=noise, velScale=velscale,
                  start=pp.sol[:2], goodpixels=goodPixels, clean=False,
                  moments=n_moments, vsyst=dv, quiet=True,
                  reddening=reddening, lam=lam, bias=bias,
                  oversample=oversample, regul=regul, reg_dim=reg_dim)
        print 'Second pass results'
        print '\t', pp.sol
        print '\tChi2/DOF =', pp.chi2

        # print np.mean(np.abs(pp.bestfit - galaxy/np.median(galaxy)))

        print 'rescaling by', pp.chi2

        noise *= np.sqrt(pp.chi2)

        # rescale noise by factor sqrt(pp.chi2)

        # add a gas component, if called for
        # this requires making `moments` into lists with length equal to
        # the number of components, as seen in `ppxf_population_gas_example.py`

        nStar_templates = templates.shape[-1]
        nGas_templates = 0
        # re-initialize the start values, in preparation for another run, if
        # necessary
        start = pp.sol[:2]

        if gas_comps not in (0, None):
            # BE CAREFUL!!!
            #`nGas_templates` is the number of gas emission lines that
            # are given for each kinematic component
            #`gas_comps` is the number of gas components
            # if you confuse these two, you're gonna have a bad time

            # also, gas_comps > 1 untested!

            gas_templates, line_names, line_wave = util.emission_lines(
                logLam2, lamRange1, FWHM_gal)
            nGas_templates = gas_templates.shape[-1]
            gas_templates = np.tile(gas_templates, (1, gas_comps))
            # print gas_templates.shape
            # print templates.shape
            all_templates = np.hstack([templates, gas_templates])

            # now modify the starting kinematics to permit multiple components
            moments = [n_moments] + [-2] + [2] * (gas_comps - 1)
            start = [start] * (gas_comps + 1)
        else:
            gas_comps = 0
            all_templates = templates

        # print all_templates.shape

        print 'gas components:', gas_comps
        print 'emission lines per gas component:', nGas_templates

        component = [0] * nStar_templates

        for gascomp_number in range(1, gas_comps + 1):
            component.extend([gascomp_number] * nGas_templates)

        # print len(component)

        pp = ppxf(templates=all_templates, galaxy=galaxy/np.median(galaxy),
                  noise=noise, velScale=velscale, start=start,
                  goodpixels=goodPixels, clean=clean, plot=False,
                  moments=n_moments, vsyst=dv, reddening=reddening,
                  lam=lam, bias=bias, quiet=quiet, oversample=oversample,
                  regul=regul, reg_dim=reg_dim, component=component)

        print z, lam[0], lam[-1]

        print 'Final pass results'
        if gas_comps > 0:
            print '\t', pp.sol[0]
            # shift lam slightly, based on sol
            #(similar statement for gas-less case is in else statement)
            lam = np.exp(logLam1 + pp.sol[0][0]/c)
        else:
            print '\t', pp.sol
            lam = np.exp(logLam1 + pp.sol[0]/c)
        print '\tChi2/DOF =', pp.chi2

        if gas_comps not in (0, None):
            for i in range(1, gas_comps + 1):
                gas = pp.matrix[:, -nGas_templates:].dot(
                    pp.weights[-nGas_templates:])
                # Extract weights of gas emissions
                w = np.where(np.array(component) == gas_comps)[0]
                print 'Gas component', i
                print '\t', pp.sol[i][:2]
                print 'Emission lines peak intensity:'
                for name, weight, line in zip(
                        line_names, pp.weights[w], pp.matrix[:, w].T):
                    print('\t %12s: %.3g' % (name, weight*np.max(line)))

        # print 'Third pass successful'

        if fit_plots == True:

                # create a figure
            plt.close('all')
            fig = plt.figure(figsize=(8, 6))

            import matplotlib.gridspec as gridspec
            gs = gridspec.GridSpec(4, 1, height_ratios=[3, 1, 1.25, 2])
            ax1 = plt.subplot(gs[0])

            # first plot the input spectrum
            ax1.plot(pp.lam, pp.galaxy, c='k', label='galaxy')
            ax1.fill_between(
                pp.lam,
                (galaxy - noise)/np.median(galaxy),
                (galaxy + noise)/np.median(galaxy), edgecolor='#ff5f00',
                facecolor='coral', alpha=0.5)
            # best fit
            ax1.plot(pp.lam, pp.bestfit, c='r', linewidth=2, label='pPXF fit')
            # residuals
            mn = np.min(pp.bestfit[pp.goodpixels])
            mx = np.max(pp.bestfit[pp.goodpixels])
            resid = mn + pp.galaxy - pp.bestfit
            mn1 = np.min(resid[pp.goodpixels])
            ax1.plot(pp.lam[pp.goodpixels], resid[pp.goodpixels],
                     marker='.', markersize=2, c='cyan',
                     markeredgecolor='cyan', linestyle='None', zorder=1)
            ax1.plot(pp.lam[pp.goodpixels], pp.goodpixels*0 + mn,
                     marker=',', c='k', zorder=0)

            w = np.where(np.diff(pp.goodpixels) > 1)[0]
            if w.size > 0:
                for wj in w:
                    x = np.arange(pp.goodpixels[wj], pp.goodpixels[wj+1])
                    ax1.plot(pp.lam[x], resid[x], 'indigo')
                w = np.hstack([0, w, w+1, -1])  # Add first and last point
            else:
                w = [0, -1]
            for gj in pp.goodpixels[w]:
                ax1.plot([pp.lam[gj], pp.lam[gj]], [mn, pp.bestfit[gj]],
                         color='orange', linewidth=0.5)

            # turn off tick labels for x axis
            ax1.spines['bottom'].set_visible(False)
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_ylabel("Counts", fontsize=16)
            ax1.legend(loc='best')

            # set up a twin axis to display pixel positions
            # DO NOT CHANGE XLIMS OF ANYTHING!!!!

            ax1_pix = ax1.twiny()
            ax1_pix.plot(np.arange(NAXIS1), np.zeros(NAXIS1))
            ax1_pix.set_xlim([0, NAXIS1 - 1])
            ax1_pix.set_xlabel('Pixel')

            ax1.set_ylim([mn1, mx] + np.array([-0.05, 0.05])*(mx-mn1))

            # set up an axis to display residuals vs noise

            ax1_res = plt.subplot(gs[1], sharex=ax1)
            ax1_res.plot(pp.lam[pp.goodpixels],
                         (noise/galaxy)[pp.goodpixels], marker='.',
                         c='coral', linestyle='None', markersize=2,
                         label='noise')
            ax1_res.plot(pp.lam[pp.goodpixels],
                         (resid/pp.bestfit)[pp.goodpixels], marker='.',
                         c='cyan', linestyle='None', markersize=2,
                         label='resid', alpha=0.5)
            ax1_res.set_xlabel(r'$\lambda_r ~ [\AA]$')

            ax1_res.legend(loc='best', prop={'size': 8})

            _ = [tick.label.set_fontsize(8) for tick in
                 ax1_res.yaxis.get_major_ticks()]

            ax1.set_xlim([np.min(pp.lam), np.max(pp.lam)])

            # this is just a dummy axis for spacing purposes
            pad_ax = plt.subplot(gs[2])
            pad_ax.axis('off')

            ax2 = plt.subplot(gs[3])
            # print s

            # make a weights array for display, s.t. # rows = # Z vals
            # and # cols = # age vals

            # extract the kinematics of the stars first
            weights = np.reshape(
                pp.weights[:nStar_templates],
                (len(np.unique(ssps['Z'])), len(np.unique(ssps['t']))))
            if weights.sum() > 0.:
                weights /= weights.sum()
            # print weights

            plt.imshow(
                weights, origin='lower', interpolation='nearest',
                cmap='cubehelix_r',
                vmin=0.0, extent=(
                    np.log10(np.min(ssps['t'])) - d_log_ages/2.,
                    np.log10(np.max(ssps['t'])) + d_log_ages/2.,
                    np.min(ssps['Z']) - dz/2.,
                    np.max(ssps['Z']) + dz/2.))

            plt.colorbar()
            plt.title("Mass Fraction", size=16)
            plt.xlabel(r'$\log_{10} \tau ~ [\mathrm{Gyr}]$', size=16)
            plt.ylabel(r'$[M/H]$', size=16)

            plt.suptitle(
                objname + ' pPXF fit: fiber ' + str(fiber),
                size=18)
            plt.subplots_adjust(hspace=0.01, top=0.85)
            # plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        '''
		print("Formal errors:")
		print("     dV    dsigma   dh3      dh4")
		print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

		print('Elapsed time in PPXF: %.2f s' % (clock() - t))
		print 'Best-fitting redshift z:', (z + 1)*(1 + pp.sol[0]/c) - 1, \
		'+/-', np.abs((pp.error*np.sqrt(pp.chi2)/c))[0]
		'''

    ssps.add_column(Column(name='fits', data=pp.weights[:nStar_templates]))

    # ssps.pprint(max_lines = -1)
    # print ssps#[ssps['fits'] != 0.]

    if fit_plots == True:
        if save_fits == False:
            plt.show()
        else:
            try:
                os.makedirs(objname)
            except OSError:
                pass

                print 'writing to ' + objname + '/' + str(fiber) + '.png'
                plt.savefig(objname + '/' + str(fiber) + '.png')

    return pp, ssps


def pPXF_make_derived_plots(objname, v_offset=0.):
    import astropy.io.ascii as ascii

    fiberdata = ascii.read(objname + '/fiberfits.dat')

    gal_im_fiber_plot(
        objname=objname, fibers=fiberdata, quantity='t',
        qty_dets='[Gyr]', save=True)
    gal_im_fiber_plot(
        objname=objname, fibers=fiberdata, quantity='Z',
        qty_dets='[M/H]', save=True)
    gal_im_fiber_plot(
        objname=objname, fibers=fiberdata, quantity='V',
        qty_dets='[km/s]', save=True, offset=v_offset, oe='odd')
    gal_im_fiber_plot(
        objname=objname, fibers=fiberdata,
        quantity='sigma', qty_dets='[km/s]', save=True)

    gal_rad_dep_plot(
        objname=objname, fibers=fiberdata,
        quantity='t', qty_dets='[Gyr]', save=True, fit=True)
    gal_rad_dep_plot(
        objname=objname, fibers=fiberdata,
        quantity='Z', qty_dets='[M/H]', save=True, fit=True)
    gal_rad_dep_plot(
        objname=objname, fibers=fiberdata,
        quantity='V', qty_dets='[km/s]', save=True, offset=v_offset)
    gal_rad_dep_plot(
        objname=objname, fibers=fiberdata,
        quantity='sigma', qty_dets='[km/s]', save=True, fit=True)


def pPXF_run_galaxy(objname, first_few=None, gas_comps=None, regul=100.):
    import ppxf
    from astroML.datasets import fetch_sdss_spectrum
    import warnings
    import matplotlib.pyplot as plt
    import astropy.io.fits as fits
    import astropy.io.ascii as ascii
    import astropy.table as table
    import numpy as np
    import os
    from plotting_tools import rejection_sample_2d

    warnings.filterwarnings(
        'ignore', message="Polyfit may be poorly conditioned")
    warnings.filterwarnings(
        "ignore", message="using a non-integer number instead of an integer \
        will result in an error in the future")

    plt.rcParams['figure.figsize'] = 16, 12

    objects = {
        'NGC2558': {'z': 0.0166076, 'EBV': .0417,
                    'sdss_spec': (1927, 53321, 290)},
        'NGC2557': {'z': 0.0161257, 'EBV': .0366,
                    'sdss_spec': (1927, 53321, 387)},
        'NGC2562': {'z': 0.0168859, 'EBV': .0376,
                    'sdss_spec': (1927, 53321, 186)},
        'NGC2563': {'z': 0.0150407, 'EBV': .0381,
                    'sdss_spec': (1927, 53321, 194)},
        'NGC4061': {'z': 0.0244529, 'EBV': .0290,
                    'sdss_spec': (2608, 54474, 514)},
        'NGC6109': {'z': 0.0298147, 'EBV': .0196,
                    'sdss_spec': (1057, 52522, 391)},
        'NGC5846': {'z': None, 'EBV': .0477,
                    'sdss_spec': (None, None, None)},
        'UGC04344': {'z': 0.0167887, 'EBV': .0439,
                     'sdss_spec': (1927, 53321, 230)},
        'UGC04329': {'z': 0.0136599, 'EBV': .0440,
                     'sdss_spec': (4479, 55592, 204)}
    }

    # reddenings from IRSA (Schlafly & Finkbeiner 2011), using 5' radius from
    # center

    z = objects[objname]['z']
    EBV = objects[objname]['EBV']
    dz = 0.

    fiberdata = ascii.read('fiberdata.dat')
    # print fiberdata

    im, fiberflat = load_ifus_precorrection(objname)

    plate = objects[objname]['sdss_spec'][0]
    mjd = objects[objname]['sdss_spec'][1]
    sdss_fiber = objects[objname]['sdss_spec'][2]
    sdss = fetch_sdss_spectrum(plate, mjd, sdss_fiber)

    dz, ifu_corr, corr_poly = sdss_cal(
        im, fiberflat, sdss, dz=0., z=z, verbose=False,
        blur=75, full_output=True)
    write_corr_frame(ifu_corr, im, z, dz, objname, verbose=False)

    fiberdata.add_column(
        table.Column(name='Z', data=np.nan*np.ones(len(fiberdata['row']))))
    fiberdata.add_column(
        table.Column(name='t', data=np.nan*np.ones(len(fiberdata['row']))))
    fiberdata.add_column(
        table.Column(name='V', data=np.nan*np.ones(len(fiberdata['row']))))
    fiberdata.add_column(
        table.Column(name='sigma', data=np.nan*np.ones(len(fiberdata['row']))))
    fiberdata.add_column(
        table.Column(name='dV', data=np.nan*np.ones(len(fiberdata['row']))))
    fiberdata.add_column(
        table.Column(
            name='dsigma', data=np.nan*np.ones(len(fiberdata['row']))))

    for i, fiber_entry in enumerate(fiberdata[fiberdata['row'] != -9999]):
        # print fiber_entry
        fiber = fiber_entry['row']
        # if fiber != -9999:
        print 'Running fiber', fiber

        n_moments = 2

        ifu = fits.open(objname + '_fluxcal.fits')[0].data
        h = fits.open(objname + '_fluxcal.fits')[0].header
        NAXIS1 = h['NAXIS1']
        CRVAL1 = h['CRVAL1']
        CDELT1 = h['CDELT1']

        save_fits = True

        try:
            pp, ssps = SP_pPXF(
                (ifu.T/np.median(ifu, axis=1)).T, fiber=fiber,
                l_summ=(CRVAL1, CDELT1, NAXIS1), z=z + dz, verbose=False,
                noise_plots=False, fit_plots=True, save_fits=save_fits,
                clean=True, quiet=True, age_lim=13.5, n_moments=n_moments,
                bias=None, objname=objname, oversample=False, reddening=EBV,
                gas_comps=gas_comps, regul=regul)

            Zs, ts = np.unique(ssps['Z']), np.unique(
                ssps['t'])  # all possible Z and t values
            # Zs are in axis 0, ts in axis 1

            # print pp.weights.shape, len(Zs), len(ts)

            '''pop_samples = rejection_sample_2d(
                    Zs, ts, pp.weights.reshape(len(Zs), len(ts)))
                # no need to normalize here
                Z_marg = pp.weights.reshape(len(Zs), len(ts)).sum(axis=1)
                t_marg = pp.weights.reshape(len(Zs), len(ts)).sum(axis=0)
                Z_samples = pop_samples[:, 0]
                t_samples = pop_samples[:, 1]

                fiber_Z = np.array(np.percentile(Z_samples, [16., 50., 84.]))
                fiber_t = np.array(np.percentile(t_samples, [16., 50., 84.]))'''

            fiber_Z_avg = np.log10(
                np.average(10.**ssps['Z'], weights=ssps['fits']))
            fiber_age_avg = np.average(ssps['t'], weights=ssps['fits'])
            fiber_V = pp.sol[0]
            fiber_sigma = pp.sol[1]
            fiber_V_err = pp.error[0]
            fiber_sigma_err = pp.error[1]

            # fiberdata['Z'] and ['t'] have 16th, 50th, and 84th percentiles

            fiberdata['Z'][np.where(fiberdata['row'] == fiber)] = fiber_Z_avg
            fiberdata['t'][np.where(fiberdata['row'] == fiber)] = fiber_age_avg

            fiberdata['V'][np.where(fiberdata['row'] == fiber)] = fiber_V
            fiberdata['sigma'][
                np.where(fiberdata['row'] == fiber)] = fiber_sigma
            fiberdata['dV'][np.where(fiberdata['row'] == fiber)] = fiber_V_err
            fiberdata['dsigma'][
                np.where(fiberdata['row'] == fiber)] = fiber_sigma_err

            print 'Fiber', fiber, 'done!'

        except ZeroDivisionError:
            print 'No fit in fiber', fiber

        except TypeError:  # occasionally you get a bad mpfit
            print 'probably velScale/mpfit error in fiber', fiber

        if first_few:
            if i >= (first_few - 1):
                break

    if not first_few:
        try:
            os.makedirs(objname)
        except OSError:
            pass
        fiberdata.write(objname + '/fiberfits.dat', format='ascii')
        print 'Written to', objname + '/fiberfits.dat'


def find_voffset(objname):
    import astropy.io.ascii as ascii
    import lts_planefit  # cappellari's lts_planefit
    import matplotlib.pyplot as plt
    import numpy as np

    fiberfits = ascii.read(objname + '/fiberfits.dat')
    plt.close('all')

    fiberfits = fiberfits[fiberfits['sky'] != 1]
    fiberfits = fiberfits[np.isnan(fiberfits['V']) != True]

    # calculate the 1-sigma error bars for each fiber position
    # 68.3% of light is contained within dr**2. = dx**2 + dy**2.
    dx = np.sqrt(0.683 * (4.687/2)**2.) * np.ones(len(fiberfits))
    dy = dx

    pf = lts_planefit.lts_planefit(
        fiberfits['ra'], fiberfits['dec'], fiberfits['V'],
        sigx=dx, sigy=dy, sigz=fiberfits['dV'],
        pivotx=0., pivoty=0., plot=True, text=False, frac=0.5)

    plt.savefig(objname + '/V_plane.png')

    # just return `a`
    return pf.abc[0]
