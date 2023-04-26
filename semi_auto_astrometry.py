"""
Semi-automated method for determining accurate WCS.
This method will rely on manually correcting one chip and then 
assume the correction is the same for the remaining chips. The remaining 
chips will then be automatically matched with gaia stars. If at some point there isn't 
a good match, manual correction will be required and the process continues.
"""

import glob
import numpy as np
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
    files = np.sort(glob.glob(f'{DIRECTORY}iff*c2*.wcs.fits'))
    done_files = np.sort(glob.glob(f'{DIRECTORY}iff*c2*.wcs.wcs_aligned.fits'))
    done_files = [done_file.replace('wcs_aligned.','') for done_file in done_files]

    to_do_files = np.setdiff1d(files, done_files)

    chip1 = ChipImage(to_do_files[0])
    wcs = chip1.determine_wcs_manually()
    offset_x, offset_y = chip1.gaia_offset_x, chip1.gaia_offset_y
    chip1.update_wcs(wcs)

    for file in to_do_files[1:]:
        chip2 = ChipImage(file)
        wcs = chip2.determine_wcs_with_offsets(x_pix_offset = offset_x, y_pix_offset = offset_y)
        chip2.update_wcs(wcs)
        offset_x = chip2.gaia_offset_x
        offset_y = chip2.gaia_offset_y
