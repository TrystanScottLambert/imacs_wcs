"""
WCS pipeline for IMACS data. Module for performing all steps for 
adding wcs information to a directory of science IMACS data.
"""

import os
import glob
from rich.progress import track
import warnings
import logging
import numpy as np

from imacs_wcs import add_wcs_to_fits
from semi_auto_astrometry import do_semi_automation_of_every_chip
from refine_wcs import align_image_to_gaia

#Turn off logging and warnings used by other packages
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore")

def get_unprocessed_files(raw_files: list[str], reduced_files: list[str]) -> np.ndarray:
    """
    Takes the raw files and determines which of them still need to be reduced.
    The raw files which don't have counter parts in the reduced files are the 
    ones that still need to be reduced.
    """
    raw_file_names = [raw_file.split('.fits')[0] for raw_file in raw_files]
    missing_indicies = []
    for i, raw_file_name in enumerate(raw_file_names):
        done = False
        for reduced_file in reduced_files:
            if raw_file_name in reduced_file:
                done = True
        if not done:
            missing_indicies.append(i)

    return raw_files[np.array(missing_indicies).astype(int)]


def get_raw_from_all(all_files: list[str]) -> list[str]:
    """Returns the raw files from all the .fits files"""
    raw = [file for file in all_files if file.split('.fits')[0][-1].isdigit()]
    return np.sort(raw)


def correct_rough_wcs(files: list[str]) -> None:
    """
    Uses the iraf script routine to add rough initial 
    wcs information.
    """
    for file in track(files, description='Adding rough wcs information to files...'):
        add_wcs_to_fits(file)

def refine_wcs(files: list[str]) -> None:
    """
    Refines the wcs based on matching to gaia.
    """
    for file in track(files, description='Refining wcs information...'):
        align_image_to_gaia(file)

def clean_up_directory(directory: str) -> None:
    """
    Removes file extensions and intermittent files.
    """

    refined_files = np.sort(glob.glob(directory + '*refined.fits'))
    for file in refined_files:
        os.rename(file, file.replace('.wcs.wcs_aligned.refined', '.wcs_corrected.'))

    os.system('mkdir temp')
    os.system('mv *.wcs.* temp/')
    os.system('mv temp/*wcs_corrected* ./')
    os.system('rm -r temp/')


def solve_wcs_for_directory(directory: str, clean: bool = True) -> None:
    """
    Performs the wcs pipeline on the given directory.
    Directory must include the last '/'
    """

    all_files = np.sort(glob.glob(directory + '*.fits'))
    raw_files = get_raw_from_all(all_files)
    rough_files = np.sort(glob.glob(directory + '*.wcs.fits'))
    to_do = get_unprocessed_files(raw_files, rough_files)

    correct_rough_wcs(to_do)
    do_semi_automation_of_every_chip(directory)

    semi_auto_files = np.sort(glob.glob(directory + '*.wcs_aligned.fits'))
    refined_files = np.sort(glob.glob(directory + '*refined.fits'))
    to_do = get_unprocessed_files(semi_auto_files, refined_files)
    refine_wcs(to_do)

    if clean:
        clean_up_directory(directory)

if __name__ == '__main__':
    DIR = '/home/tlambert/Desktop/IMACS_analysis/IMACS_RAWDATA/ut211024_25/SCIENCE/'
    solve_wcs_for_directory(DIR)
