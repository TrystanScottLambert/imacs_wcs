"""
Semi-automated method for determining accurate WCS.
This method will rely on manually correcting one chip and then 
assume the correction is the same for the remaining chips. The remaining 
chips will then be automatically matched with gaia stars. If at some point there isn't 
a good match, manual correction will be required and the process continues.
"""

import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from manual_astrometry import ChipImage

def set_up_files(directory: str, chip_number: int) -> tuple[str,list[str]]:
    """
    reads in the files, sorts them, and determines the primary file.
    """
    # .wcs.fits is the extension for rough wcs alignment.
    files = np.sort(glob.glob(f'{directory}*{chip_number}*.wcs.fits'))
    return files[0], files[1:]


if __name__ == '__main__':
    DIRECTORY = '/home/tlambert/Downloads/g_band/SCIENCE/'
    chip1 = ChipImage(f'{DIRECTORY}iff0158c1.wcs.fits')
    wcs = chip1.determine_wcs_with_gaia()
    chip1.update_wcs(wcs)

    chip2 = ChipImage(f'{DIRECTORY}iff0159c1.wcs.fits')
    chip2.determine_wcs_with_gaia(manual=False, x_pix_offset=chip1.gaia_offset_x, y_pix_offset=chip1.gaia_offset_y, show_alignment=True)