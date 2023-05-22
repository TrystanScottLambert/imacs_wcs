"""
Does a final refinement based on cross matching. 
"""

import glob
import keyring
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
from sip_tpv import sip_to_pv

from manual_astrometry import ChipImage

class RefinedAlignment(ChipImage):
    """
    Performs cross matching to refine the wcs alignment with gaia.
    """
    def __init__(self, file_name):
        super().__init__(file_name)
        self.wcs = WCS(self.header)
        self.ra, self.dec = self.wcs.pixel_to_world_values(self.data.shape[1]/2, self.data.shape[0]/2)
        self.star_catalog = self.find_stars()

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
        sources = self.star_catalog
        x_pos, y_pos = np.array(sources['xcentroid']), np.array(sources['ycentroid'])
        return x_pos, y_pos

    def match_to_astrometry(self) -> fits.Header:
        """
        Uploads the catalog of detected stars to astrometry.net to determine wcs
        """
        ast = AstrometryNet()
        ast.api_key = keyring.get_password('astroquery:astrometry_net', 'IMACS')
        sources = self.star_catalog
        sources.sort('flux')
        # Reverse to get descending order
        sources.reverse()

        image_width = self.data.shape[1]
        image_height = self.data.shape[0]
        wcs_header = ast.solve_from_source_list(
            sources['xcentroid'], sources['ycentroid'],
            image_width, image_height, solve_timeout=120,
            center_ra = float(self.ra), center_dec = float(self.dec), radius = 0.2,
            allow_commercial_use='n')
        return wcs_header

    '''def match_to_astrometry(self) -> fits.Header:
        """
        Uploads the catalog of detected stars to astrometry.net to determine wcs
        """
        ast = AstrometryNet()
        ast.api_key = keyring.get_password('astroquery:astrometry_net', 'IMACS')
        wcs_header = ast.solve_from_image(
            self.file_name, force_image_upload=True, solve_timeout=900,
            center_ra = float(self.ra), center_dec = float(self.dec), radius = 0.2, allow_commercial_use='n')
        return wcs_header'''

    def match_to_gaia(self) -> tuple[np.array, SkyCoord]:
        """
        Performs a cross match to gaia returning the pixel values
        of the image and the sky values of gaia.
        """
        x_image, y_image = self.get_star_positions()
        image_ra, image_dec = self.wcs.pixel_to_world_values(x_image, y_image)
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
        wcs = utils.fit_wcs_from_points(pix_coords, sky_coords, sip_degree=7)
        header = wcs.to_header(relax=True)
        header['CD1_1'] = header['PC1_1']
        header['CD1_2'] = header['PC1_2']
        header['CD2_1'] = header['PC2_1']
        header['CD2_2'] = header['PC2_2']

        del header['PC1_1']
        del header['PC1_2']
        del header['PC2_1']
        del header['PC2_2']
        del header['CDELT1']
        del header['CDELT2']
        sip_to_pv(header)
        return header

    def write_wcs_header(self, wcs_header: fits.Header) -> None:
        """
        Writes a wcs fits header object to the header of the 
        image.
        """
        self.hdul[0].header.update(wcs_header)
        del self.hdul[0].header['PC1_1']
        del self.hdul[0].header['PC1_2']
        del self.hdul[0].header['PC2_1']
        del self.hdul[0].header['PC2_2']
        del self.hdul[0].header['CDELT1']
        del self.hdul[0].header['CDELT2']

        self.hdul.writeto(self.file_name.split('.fits')[0] + '.refined.fits', overwrite=True)

def refine_alignment(file_name: str, method: str) -> None:
    """
    Aligns an image using a given method. 

    Methods available:
        - Doing a cross match to gaia sources. 
        - Uploading image to astrometry.net (To use this method you need to get an API 
            the website https://nova.astrometry.net/ and run the command
            import keyring 
            keyring.set_password('astroquery:astrometry_net', 'IMACS', 'apikeyhere'))
    """
    align = RefinedAlignment(file_name)
    if method == 'Gaia':
        wcs_header = align.calculate_aligned_wcs()
    elif method == 'Astrometry':
        wcs_header = align.match_to_astrometry()
    else:
        raise ValueError('Method name is incorrect.')
    align.write_wcs_header(wcs_header)

if __name__ == '__main__':
    DIR = '/home/tlambert/Desktop/IMACS_analysis/IMACS_RAWDATA/ut211024_25/SCIENCE/'
    files = np.sort(glob.glob(f'{DIR}*.wcs.fits'))
    done_files = np.sort(glob.glob(f'{DIR}*wcs.refined.test.fits'))
    done_files = [file.replace('refined.test.','') for file in done_files]
    to_do_files = np.setdiff1d(files, done_files)

    for file in track(to_do_files, description='Solving with astrometry.net'):
        print(file)
        align = RefinedAlignment(file)
        wcs_thing = align.match_to_astrometry()
        align.write_wcs_header(wcs_thing)
