# imacs_wcs

We have written a simple procedure based on the iraf module "iwcs.cl" to add WCS information into the headers of raw IMACS photometric images. This then needs to be further refined using GAIA.

## Warning
This is a first attempt at doing this and is soley based on the iwcs.cl. No tests have been written and the python code has been written in such a way as to keep with the naming conventions used in the iwcs.cl. There may be fringe cases were this method does not work that we are not aware of. If so please let us know. This is also a first attempt at generating a wcs image.


## Installation

```python
pip install imacs_wcs
```

## Making rough wcs correction
First correct the imacs images using the routine defined in iwcs.cl. Note that the module has been written for IMACS naming conventions (e.g. ift1111c1.fits).
imacs_wcs has a method for apply rough wcs information for individual images.
```python
from imacs_wcs import imacs_wcs
infile = 'ift1067c4.fits'
imacs_wcs.add_wcs_to_fits(infile)
```
You can do this for the entire directory using glob if you wish.

```python
import glob
files = glob.glob('science_directory/*.fits')
for file in files:
  imacs_wcs.add_wcs_to_fits(file)
```

## Making a full wcs correction
From testing we find that the rough wcs is about 20" offset. To correct for this the user must align the images with gaia. To do this the manual_astrometry and semi_auto_astrometry modules have been included. These modules will align all images in a directory. For example, if all the imacs data was in a directory called "science" (and have been roughly corrected) then we could align them like: 

```python
import imacs_wcs.semi_auto_astrometry as semi
directory = 'science/'

# for a single chip number (if you want to)
semi.do_semi_automation_for_chip(directory)

# for every chip in a directory (reccomended)
semi.do_semi_automation_of_every_chip(directory)
```

The user will be promted if the program cannot automatically fit a solution. If, for some reason, the process if interrupted it can be run again because the alignment will only run on images which haven't already been aligned.

 ## Future
 
 If we can confirm that this method is working correctly then there are some things we would like to implement 
  - write tests
  - write further auto alignment
  
 It is also possible that this may form part of a purely python-based pipeline for reducing IMACS data.
