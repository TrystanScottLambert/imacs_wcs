"""
Does a final refinement based on cross matching. 
"""

import glob
from rich.progress import track
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import utils
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astroquery.astrometry_net import AstrometryNet

from manual_astrometry import ChipImage

class RefinedAlignment(ChipImage):
    """
    Performs cross matching to refine the wcs alignment with gaia.
    """

    def find_stars(self, fwhm: float = 3.0, sigma: float = 5.0) -> Table:
        """
        Performs DOA star finding on the image and returns
        an astropy table with the results.

        fwhm should be given as the number of pixels that a
        star in the image has. Sigma is the number of standard
        deviations above the noise the sources should be.
        """
        _, median, std = sigma_clipped_stats(self.data, sigma=3.0)
        daofind = DAOStarFinder(fwhm=fwhm, threshold=sigma*std)
        sources = daofind(self.data - median)
        return sources

    def get_star_positions(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Determines the xy positions of the 2d data array.
        """
        sources = self.find_stars()
        x_pos, y_pos = np.array(sources['xcentroid']), np.array(sources['ycentroid'])
        return x_pos, y_pos

    def match_to_astrometry(
            self, api_key: str, fwhm: float = 3.0, sigma: float = 5.0) -> fits.Header:
        """
        Uploads the catalog of detected stars to astrometry.net to determine wcs
        """
        ast = AstrometryNet()
        ast.api_key = api_key
        sources = self.find_stars(fwhm, sigma)
        sources.sort('flux')
        sources.reverse()
        image_width = self.data.shape[1]
        image_height = self.data.shape[0]
        wcs_header = ast.solve_from_source_list(sources['xcentroid'], sources['ycentroid'],
                                        image_width, image_height,
                                        solve_timeout=300)
        return wcs_header


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

    def calculate_aligned_wcs(self) -> fits.Header:
        """
        works out the correct wcs based on the alignment to gaia.
        """
        pix_coords, sky_coords = self.match_to_gaia()
        wcs = utils.fit_wcs_from_points(pix_coords, sky_coords)
        return wcs.to_header()

    def write_wcs_header(self, wcs_header: fits.Header) -> None:
        """
        Writes a wcs fits header object to the header of the 
        image.
        """
        self.hdul[0].header.update(wcs_header)
        self.hdul.writeto(self.file_name.split('.fits')[0] + '.refined.test.fits', overwrite=True)

if __name__ == '__main__':
    api_key = 'juseocalskdsoavj'
    DIR = '/home/tlambert/Desktop/IMACS_analysis/IMACS_RAWDATA/ut211024_25/SCIENCE/'
    files = np.sort(glob.glob(f'{DIR}*.wcs.fits'))
    done_files = np.sort(glob.glob(f'{DIR}*wcs.refined.test.fits'))
    done_files = [file.replace('.refined.test.','') for file in done_files]
    to_do_files = np.setdiff1d(files, done_files)

    for file in track(to_do_files, description='Solving with astrometry.net'):
        align = RefinedAlignment(file)
        wcs_thing = align.match_to_astrometry(api_key=api_key)
        align.write_wcs_header(wcs_thing)
