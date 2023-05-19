"""
Script for working out the 
"""

import os
import glob
from dataclasses import dataclass
from rich.progress import track
import numpy as np
from astropy.io import fits

DEFAULT_SEX_CONFIG = 'default.sex'
DEFAULT_SEX_PARAMS = 'default.param'
DEFAULT_SCAMP_CONFIG = 'default.scamp'

@dataclass
class ConfigurationFiles:
    """Location of all the configuration files"""
    sex_config_file: str
    sex_param_file: str
    scamp_config_file: str

def sed(infile: str, old_string: str, new_string: str, outfile: str) -> None:
    """Pythonic version of sed with 's/string/newstring/g' flags."""
    with open(infile, encoding='utf8') as file:
        lines = file.readlines()
    lines = [line.replace(f'{old_string}', f'{new_string}') for line in lines]
    with open(outfile, 'w', encoding='utf8') as new_file:
        for line in lines:
            new_file.write(line)

def create_sex_config(image_name: str, sex_config_file: str) -> None:
    """Changes the name of default.cat file to the image_name.cat."""
    sed(sex_config_file,'test.cat', image_name.replace('.fits','.cat'), f'{image_name.replace(".fits", ".sex")}')

def run_sextractor(image_name: str, config_file: str) -> None:
    """Runs sextractor on the given image using the given config files."""
    os.system(f'sex {image_name} -c {config_file}')

def run_scamp(sextractor_catalog: str, config_file: str) -> None:
    """Runs scamp on the given image  with the given config file"""
    os.system(f'scamp {sextractor_catalog} -c {config_file}')

def create_head_file(image_name: str, config_files: ConfigurationFiles) -> None:
    """Performs the steps to make a .head file which can be used by swarp."""
    create_sex_config(image_name, config_files.sex_config_file)
    run_sextractor(image_name, image_name.replace('.fits','.sex'))
    run_scamp(image_name.replace('.fits','.cat'), image_name.replace('.fits','.scamp'))

    os.remove(image_name.replace('.fits','.sex'))
    os.remove(image_name.replace('.fits','.cat'))

if __name__ == '__main__':
    DIR = '/home/tlambert/Desktop/IMACS_analysis/IMACS_RAWDATA/ut211024_25/SCIENCE/'
    files = np.sort(glob.glob(f'{DIR}*.wcs_aligned.fits'))
    settings = ConfigurationFiles('default.sex','default.param', 'default.scamp')
    for file in track(files):
        create_head_file(file, settings)
