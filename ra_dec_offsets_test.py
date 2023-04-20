"""
Script to investigate the possible offsets in RA and DEC using 
an image which whose WCS was manually determined as a test.
"""

import glob
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import position_angle, angular_separation
import header_reader

if __name__ == '__main__':
    files = np.sort(glob.glob('/home/tlambert/Downloads/g_band/scamp_test/*158*wcs.fits*'))
    reference_files = np.sort(glob.glob('/home/tlambert/Downloads/g_band/RAW_SCIENCE/*158*.fits'))
    angular_offsets_of_chips = []
    pa_chip_offsets = []
    ra_bad  = 62.6110246
    dec_bad  = -9.3008723

    ra_good = 62.6110216
    dec_good = -9.3060079
    ra_diff = ra_good - ra_bad
    dec_diff = dec_good - dec_bad
    print('ra_diff: ', ra_diff)
    print('dec_diff: ', dec_diff)
    for i, file in enumerate(files):
        hdul = fits.open(file)
        wcs = WCS(hdul[0].header)
        header_information = header_reader.HeaderInformation(reference_files[i])
        telescope_pa = header_information.pa *u.rad

        center_pix = hdul[0].data.shape
        center_x = hdul[0].header['CRPIX1']
        center_y = hdul[0].header['CRPIX2']
        center_ra, center_dec = wcs.pixel_to_world_values(center_x, center_y)
        telescope_ra = hdul[0].header['RA-D']
        telescope_dec = hdul[0].header['DEC-D']
        sep = angular_separation(center_ra *u.deg, center_dec *u.deg, telescope_ra *u.deg, telescope_dec *u.deg)
        total_pa = position_angle(center_ra, center_dec, telescope_ra, telescope_dec)
        chip_pa = total_pa - telescope_pa
        pa_chip_offsets.append(chip_pa.value)

        angular_offsets_of_chips.append(sep.to(u.deg).value)
        print('Difference RA: ', center_ra-telescope_ra)
        print('Difference Dec: ', center_dec-telescope_dec)
        print('Position Angle: ', position_angle(center_ra, center_dec, telescope_ra, telescope_dec))
        print('Angular Separation: ', angular_separation(center_ra, center_dec, telescope_ra, telescope_dec))
        print()
    print(angular_offsets_of_chips)
    print(pa_chip_offsets)
