"""
Module to overplot gaia stars on an image and then use them to
determine a WCS.
"""

from typing import Tuple
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord
import astropy.units as u
import pylab as plt
from astroquery.gaia import Gaia


def search_gaia_archives(ra: float, dec: float, height:float, width:float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Does a box search centered on ra and dec with a given height and width.
    """
    coord = SkyCoord(ra = ra*u.deg, dec=dec*u.deg)
    width = u.Quantity(width*u.deg)
    height = u.Quantity(height*u.deg)
    results = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    return np.array(list(results['ra'])), np.array(list(results['dec']))


if __name__ == '__main__':
    INFILE = '/home/tlambert/Downloads/g_band/scamp_test/iff0158c1.wcs.fits'
    #INFILE = '/home/tlambert/Downloads/g_band/SCIENCE/RAWDATA/iff0158c1.wcs.fits'

    hdul = fits.open(INFILE)
    data = hdul[0].data
    header = hdul[0].header
    wcs = WCS(header)

    center_y_pix = data.shape[0]/2
    center_x_pix = data.shape[1]/2

    # Assuming square pixels.
    search_width = data.shape[1] * np.abs(header['PC1_1'])
    search_height = data.shape[0] * np.abs(header['PC1_1'])

    center_ra, center_dec = wcs.pixel_to_world_values(center_x_pix, center_y_pix)

    gaia_ra, gaia_dec = search_gaia_archives(center_ra, center_dec, width = search_width, height=search_height)
    gaia_x, gaia_y = wcs.world_to_pixel_values(gaia_ra, gaia_dec)

    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    fig = plt.figure()
    ax = fig.add_subplot(projection=wcs)
    ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray')
    ax.scatter(gaia_x, gaia_y, facecolor='None', edgecolor='r')
    plt.show()
