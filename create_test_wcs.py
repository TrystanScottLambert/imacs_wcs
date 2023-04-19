"""
Manually creating a world coordinate system that we can test against. 
This wcs was created by visually identifying the stars on the chips and 
manually finding them in aladin. It was not fun.
"""

from typing import Tuple
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
import astropy.units as u

def read_wcs_catalog(file_name: str) -> Tuple:
    """
    Reads the wcs catalogs which have hour angles and pixels.
    Returns the world_coordinates and pix coordinates in the format
    that wcs.utils need to create a WCS.
    """
    with open(file_name, encoding='utf8') as file:
        lines = file.readlines()

    separated_lines = [line.split('\n')[0].split(' ') for line in lines]
    coordinates_string = [(line[0], line[1]) for line in separated_lines]
    world_coordinates = SkyCoord(coordinates_string, frame='icrs', unit=(u.hourangle, u.deg))

    x_values = [int(line[2]) for line in separated_lines]
    y_values = [int(line[3]) for line in separated_lines]

    pix_coordinates = np.array([x_values, y_values])
    return world_coordinates, pix_coordinates

def get_chip_number_from_file_name(file_name:str) -> int:
    """
    Gets the chip number which need to match to the correct fits file.
    """
    return int(file_name[-5])

FITS_FILES = {
    1:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c1.fits',
    2:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c2.fits',
    3:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c3.fits',
    4:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c4.fits',
    5:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c5.fits',
    6:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c6.fits',
    7:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c7.fits',
    8:'/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c8.fits',
}

def apply_wcs(star_wcs_directory: str) -> None:
    """
    Works out the wcs object based on the give star_wcs file then applies
    to the correct files.
    """
    world_coords, pix_coords = read_wcs_catalog(star_wcs_directory)
    wcs = utils.fit_wcs_from_points(pix_coords, world_coords)
    hdul = fits.open(FITS_FILES[get_chip_number_from_file_name(star_wcs_directory)])
    hdul[0].header.update(wcs.to_header())
    hdul.writeto(FITS_FILES[get_chip_number_from_file_name(star_wcs_directory)].split('.fits')[0] + '.wcs.fits', overwrite=True)

if __name__ == '__main__':
    INFILE = '/home/tlambert/Downloads/g_band/RAW_SCIENCE/iff0158c1.fits'
    STAR_WCS_DIRECTORIES = glob.glob('/home/tlambert/Downloads/g_band/*.txt')
    for file in STAR_WCS_DIRECTORIES:
        apply_wcs(file)
    

