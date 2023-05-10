"""
Semi-automated method for determining accurate WCS.
This method will rely on manually correcting one chip and then 
assume the correction is the same for the remaining chips. The remaining 
chips will then be automatically matched with gaia stars. If at some point there isn't 
a good match, manual correction will be required and the process continues.
"""

import glob
import numpy as np

from .manual_astrometry import ChipImage


def do_semi_automation_for_chip(directory: str, chip_number: int) -> None:
    """
    Runs the semi automation on a select chip.
    """
    if chip_number<1 or chip_number>8:
        raise ValueError('IMACS only has 8 chips. Chip number must be between 1 and 8.')

    files = np.sort(glob.glob(f'{directory}*c{chip_number}*.wcs.fits'))
    if len(files) == 0:
        raise FileExistsError('No rough wcs files found. Please run imacs_wcs first.')

    done_files = np.sort(glob.glob(f'{directory}*c{chip_number}*.wcs.wcs_aligned.fits'))
    done_files = [done_file.replace('wcs_aligned.','') for done_file in done_files]

    # Ignoring any files that have already been aligned.
    to_do_files = np.setdiff1d(files, done_files)

    if len(to_do_files) ==0:
        print(f'{chip_number} already done.')
        return None
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

def do_semi_automation_of_every_chip(directory: str) -> None:
    """
    Aligns every rough wcs image using the semi-automated approach.
    """
    for chip_number in range(1, 9):
        print(f'Aligning chip: {chip_number}')
        do_semi_automation_for_chip(directory, chip_number)

if __name__ == '__main__':
    DIRECTORY = '../SCIENCE/'
    do_semi_automation_for_chip(DIRECTORY, 4)
    do_semi_automation_of_every_chip(DIRECTORY)
