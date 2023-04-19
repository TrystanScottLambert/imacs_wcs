"""
Module to read necessary information from the fits header and calculate correct values.
"""

import numpy as np
from astropy.io import fits
from astropy import wcs

from instrument_constants import Constants as c
from chip import chips

def read_chip_number(fits_file_name: str) -> str:
    """Finds the chip number from the raw imacs file name."""
    chip_number = fits_file_name.split('.fits')[0].split('c')[-1]
    return f'c{chip_number}'

def calculate_wcs_matrix(
        sxz: int, syz: int, xb: int, yb: int, pixscale:float, pa:float, camera_type:str):
    """Calculates the wcs matrix values depending on camera."""
    if camera_type == 'Long':
        cd1_1 = sxz * xb * pixscale * np.cos(pa)/3600.
        cd2_2 = syz * yb * pixscale * np.cos(pa)/3600.
        cd1_2 =-sxz * xb * pixscale * np.sin(pa)/3600.
        cd2_1 = syz * yb * pixscale * np.sin(pa)/3600.

    elif camera_type == 'Short':
        cd1_1 = sxz * xb * pixscale * np.cos(pa)/3600.
        cd2_2 =-syz * yb * pixscale * np.cos(pa)/3600.
        cd1_2 =-sxz * xb * pixscale * np.sin(pa)/3600.
        cd2_1 =-syz * yb * pixscale * np.sin(pa)/3600.

    return [cd1_1, cd1_2, cd2_1, cd2_2]

class HeaderInformation:
    """Reads in all appropriate header information."""
    def __init__(self, fits_file_name: str) -> None:
        """Reading in the imacs file without wcs."""
        self.hdu_list = fits.open(fits_file_name)
        self.chip_number = read_chip_number(fits_file_name)
        self.chip = chips[self.chip_number]
        self.header = self.hdu_list[0].header

        self.camera = self.header['INSTRUME']
        binning = self.header['BINNING']
        self.rotan_d = self.header['ROTANGLE']
        self.ra1 = self.header['RA-D']
        self.dec1 = self.header['DEC-D']

        self.xb = int(binning[0])
        self.yb = int(binning[2])
        self.pa = self.north * np.pi/180

    @property
    def camera_type(self):
        """Determines if the camera is Long or Short."""
        if 'Short' in self.camera:
            camera_type = 'Short'
        elif 'Long' in self.camera:
            camera_type = 'Long'
        return camera_type

    @property
    def dor(self):
        """Determines the dewar orientation."""
        if 'Normal' in self.header['DEWARORI']:
            dor_value = 0
        elif 'Nod' in self.header['DEWARORI']:
            dor_value = -90
        return dor_value

    @property
    def north(self):
        """Determines value for north depending on camera type."""
        if self.camera_type == 'Short':
            north_value = 0.0 - (self.rotan_d - (90.0+c.iroa_d)) + self.dor
        elif self.camera_type == 'Long':
            north_value = 180.0 + (self.rotan_d - (90.0+c.iroa_d)) + self.dor
        return north_value

    @property
    def x0(self):
        """Calculates the x0 value depending on camera type and chip number."""
        if self.camera_type == 'Long':
            x0_value = chips[self.chip_number].x0_SITe
        elif self.camera_type == 'Short':
            x0_value = chips[self.chip_number].x0_E2V
        return x0_value

    @property
    def y0(self):
        """Calculates the y0 value depending on camera type and chip number."""
        if self.camera_type == 'Long':
            y0_value = chips[self.chip_number].y0_SITe
        elif self.camera_type == 'Short':
            y0_value = chips[self.chip_number].y0_SITe
        return y0_value

    @property
    def mscale(self):
        if self.camera_type == 'Long':
            value =  66.66667
        elif self.camera_type == 'Short':
            value = 67.18
        return value

    @property
    def pixscale(self):
        if self.camera_type == 'Long':
            value = 0.111
        elif self.camera_type == 'Short':
            value = 0.202*66.66667/self.mscale
        return value

    @property
    def corner_pixels(self):
        """ Works out the corner pixels and returns (crp1, crp2)"""
        nxm = c.nx_num/self.xb
        nym = c.ny_num/self.yb
        crp1_big = nxm/2
        crp2_big = nym/2
        crp1_value = (crp1_big*self.xb-c.offx-self.x0*self.mscale-self.chip.xz)/(self.xb*self.chip.sxz)
        crp2_value = (crp2_big*self.yb-c.offy-self.y0*self.mscale-self.chip.yz)/(self.yb*self.chip.syz)

        return [crp1_value, crp2_value]

    @property
    def wcs_matrix(self):
        """Calculates the indicies of the WCS matrix depending on camera."""
        return calculate_wcs_matrix(
            self.chip.sxz, self.chip.syz, self.xb, self.yb, self.pixscale, self.pa, self.camera_type)

    def generate_wcs_object(self):
        """writes the WCS object."""
        wcs_dict = {
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRVAL1': self.ra1,
            'CRVAL2': self.dec1, 
            'CRPIX1': self.corner_pixels[0], 
            'CRPIX2': self.corner_pixels[1], 
            'CD1_1': self.wcs_matrix[0],
            'CD1_2': self.wcs_matrix[1],
            'CD2_1': self.wcs_matrix[2],
            'CD2_2': self.wcs_matrix[3],
            'NAXIS1': self.header['NAXIS1'],
            'NAXIS2': self.header['NAXIS2'],
            'CUNIT1': 'deg',
            'CUNIT2': 'deg',
        }

        return wcs.WCS(wcs_dict)
