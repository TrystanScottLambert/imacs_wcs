"""
Does a final refinement based on cross matching. 
"""

import numpy as np
import matplotlib.pylab as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import utils
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils.detection import DAOStarFinder

from manual_astrometry import ChipImage

class RefinedAlignment(ChipImage):
    """
    Performs cross matching to refine the wcs alignment with gaia.
    """

    def get_star_positions(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Determines the xy positions of the 2d data array.
        """
        _, median, std = sigma_clipped_stats(self.data, sigma=3.0)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
        sources = daofind(self.data - median)
        x_pos, y_pos = np.array(sources['xcentroid']), np.array(sources['ycentroid'])
        return x_pos, y_pos

    def match_to_gaia(self) -> tuple[np.array, SkyCoord]:
        """
        Performs a cross match to gaia returning the pixel values
        of the image and the sky values of gaia.
        """
        x_image, y_image = self.get_star_positions()
        image_ra, image_dec = WCS(self.header).pixel_to_world_values(x_image, y_image)
        c_image = SkyCoord(ra = image_ra * u.deg, dec = image_dec * u.deg)
        c_gaia = SkyCoord(ra = self.gaia_ra * u.deg, dec = self.gaia_dec * u.deg)
        idx, d2d, _ = c_gaia.match_to_catalog_sky(c_image)

        # Remove outliers, in the cases where gaia sources weren't detected.
        cut = np.where(d2d < 10 * u.arcsec)[0] # hard limit
        idx = idx[cut]

        pix_coords = np.array([x_image[idx], y_image[idx]])
        sky_coords = SkyCoord([
            (self.gaia_ra[cut][i], self.gaia_dec[cut][i]) for i,_ in enumerate(self.gaia_ra[cut])
            ], unit = u.deg)
        return pix_coords, sky_coords

    def calcualte_aligned_wcs(self) -> WCS:
        """
        works out the correct wcs based on the alignment to gaia.
        """
        pix_coords, sky_coords = self.match_to_gaia()
        wcs = utils.fit_wcs_from_points(pix_coords, sky_coords)
        return wcs

    def align(self) -> None:
        """
        Aligns the current image with the correct wcs based on alignment with gaia.
        """
        wcs = self.calcualte_aligned_wcs()
        self.hdul[0].header.update(wcs.to_header())
        self.hdul.writeto(self.file_name.split('.fits')[0] + '.refined.fits', overwrite=True)

    def _show_matches(self) -> None:
        """
        Shows a plot comparing the positions of the gaia sources 
        and the cross-matched stars.
        """

        pix_coords, sky_coords = self.match_to_gaia()
        ra, dec = sky_coords.ra.value, sky_coords.dec.value
        gaia_x, gaia_y = WCS(self.header).world_to_pixel_values(ra, dec)

        image_x, image_y = pix_coords[0], pix_coords[1]
        zscale = ZScaleInterval()
        v_min, v_max = zscale.get_limits(self.data)
        fig = plt.figure()
        ax_diagnostic = fig.add_subplot(projection=WCS(self.header))
        ax_diagnostic.imshow(self.data, vmin=v_min, vmax=v_max)
        ax_diagnostic.scatter(
            image_x, image_y, marker='s', color='b', facecolor=None, label='Detected Sources')
        ax_diagnostic.scatter(gaia_x, gaia_y, facecolor=None, color='r', label='GAIA')
        ax_diagnostic.legend()
        plt.show()


def align_image_to_gaia(file_name: str) -> None:
    """ 
    Performs a cross match to gaia and uses that to determine
    the correct wcs.
    """
    alignment = RefinedAlignment(file_name)
    alignment.align()
