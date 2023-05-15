"""
WCS pipeline for IMACS data. Module for performing all steps for 
adding wcs information to a directory of science IMACS data.
"""

import glob
from rich.progress import track

from .imacs_wcs import add_wcs_to_fits
from .semi_auto_astrometry import do_semi_automation_of_every_chip
from .refine_wcs import align_image_to_gaia

def correct_rough_wcs(files: list[str]) -> None:
    """
    Uses the iraf script routine to add rough initial 
    wcs information.
    """
    for file in track(files, description='Adding rough wcs information to files'):
        add_wcs_to_fits(file)

def refine_wcs(files: list[str]) -> None:
    """
    Refines the wcs based on matching to gaia.
    """
    for file in files:
        align_image_to_gaia(file)

def clean_up_directory(directory: str) -> None:
    """
    Removes file extensions and intermittent files.
    """
    pass

def solve_wcs(directory: str, clean: bool = True) -> None:
    """
    Performs the wcs pipeline on the given directory.
    Directory must include the last '/'
    """
    files = glob.glob(directory + '*.fits')
    correct_rough_wcs(files)
    do_semi_automation_of_every_chip(directory)
    refine_wcs(files)

    if clean:
        clean_up_directory(directory)
