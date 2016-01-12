import numpy as np
from astroquery.simbad import Simbad
import montage_wrapper as m
import astropy.units as u
import os
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import table
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import aplpy
import pywcsgrid2
#from mpl_toolkits.axes_grid1.axes_rbg import imshow_rgb

objnames = ['UGC04329', 'UGC04344', 'NGC2558', 'NGC2557', 'NGC2562',
            'NGC2563']

objnames2 = ['NGC2560', 'NGC2569', 'IC2293', 'UGC04324', 'UGC04386',
             'IC2339', 'IC2341']

objnames3 = ['CGCG119-040', 'CGCG119-047', 'CGCG119-051', 'CGCG119-069']
objs3 = SkyCoord(
    ra=[124.62096, 124.82825, 124.76342, 125.36387,125.36387],
    dec=[20.89526, 21.76629, 20.82625, 21.03625, 21.03625],
    unit=(u.deg, u.deg), frame='icrs')

def get_obj_coords(names):
    for i, objname in enumerate(names):
        #print objname
        if i == 0:
            res_table = Simbad.query_object(objname, wildcard=True)
        else:
            table_to_add = Simbad.query_object(objname, wildcard=True)
            res_table = table.vstack([res_table, table_to_add])
    return SkyCoord(ra=res_table['RA'], dec=res_table['DEC'],
                    unit=(u.hourangle, u.deg), frame='icrs')


if __name__ == '__main__':
    #os.system('mHdr -p 2 NGC2563 2 NGC2563.hdr')
    #os.system('mExec -c -o NGC2563B.fits -f NGC2563.hdr DSS DSS2B tempdir')
    #os.system('mExec -c -o NGC2563R.fits -f NGC2563.hdr DSS DSS2R tempdir')
    #os.system('mExec -c -o NGC2563IR.fits -f NGC2563.hdr DSS DSS2IR tempdir')
    #os.system('mExec -c -o NGC2563K.fits -f NGC2563.hdr 2MASS Ks tempdir')

    fig = plt.figure(figsize=(10, 7), dpi=600)
    fo = fits.open('NGC2563R.fits')
    ax = pywcsgrid2.axes([0.105, 0.075, 0.6, 0.9], header=fo[0].header)
    ax.imshow(fo[0].data, origin='lower', cmap='binary',
              aspect='equal', zorder=0,
              norm=LogNorm(vmin=np.percentile(fo[0].data, 15.)))
    ax.grid(zorder=1)

    objs = get_obj_coords(objnames)
    objs2 = get_obj_coords(objnames2)

    ax['fk5'].scatter(objs.ra.deg, objs.dec.deg, facecolor='None',
                      edgecolor='g', s=50, zorder=2, marker='H')
    ax['fk5'].scatter(objs2.ra.deg, objs2.dec.deg, marker='s',
                      edgecolor='r', facecolor='None', s=40, zorder=2)
    ax['fk5'].scatter(objs3.ra.deg, objs3.dec.deg, marker='s',
                      edgecolor='r', facecolor='None', s=40, zorder=2)

    print objs.to_string('hmsdms')

    txt = ['{}: {}'.format(n, str(o.to_string('hmsdms')))
           for n, o in zip(objnames, objs)]
    txt2 = ['{}: {}'.format(n, str(o.to_string('hmsdms')))
            for n, o in zip(objnames2, objs2)]
    txt3 = ['{}: {}'.format(n, str(o.to_string('hmsdms')))
            for n, o in zip(objnames3, objs3)]

    fig.text(.725, .55, '\n'.join(txt), color = 'g', size=8)
    fig.text(.725, .30, '\n'.join(txt2 + txt3), color = 'r', size=8)
    d = 61. #dist [Mpc] to NGC2563
    ds = np.array([0.1, 0.2, 0.5, 1., 2.]) #projected radius
    print ds/d*206265./3600.

    NGC2563 = objs[-1]
    #circles = [plt.Circle((NGC2563.ra.deg, NGC2563.dec.deg), s, edgecolor='k',
    #                      facecolor='None') for s in ds/d/3600.]

    #[ax['fk5'].add_artist(c) for c in circles];
    for r, l in zip(ds/d*206265./3600., ds):
        ax['fk5'].add_artist(
            plt.Circle(
                (NGC2563.ra.deg, NGC2563.dec.deg), r, facecolor='None',
                edgecolor='k', linestyle='dashed'))

        ax['fk5'].annotate(xy=(NGC2563.ra.deg + .1, NGC2563.dec.deg + r + .01),
                           s='R$_{{proj}}$ = {:.2} Mpc'.format(l), size=8)

    plt.suptitle('NGC2563 Group Diagram')

    plt.savefig('field_diagram.png')
