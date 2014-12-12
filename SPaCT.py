'''
SparsePak Correction Tools (SPaCT)
Author: Zach Pace (U Wisc-Madison)

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