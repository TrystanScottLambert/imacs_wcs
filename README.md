# imacs_wcs

We have written a simple procedure based on the iraf module "iwcs.cl" (included here for reference) to add WCS information into the headers of raw IMACS photometric images. 

## Warning
This is a first attempt at doing this and is soley based on the iwcs.cl. No tests have been written and the python code has been written in such a way as to keep with the naming conventions used in the iwcs.cl. There may be fringe cases were this method does not work that we are not aware of. If so please let us know. This is also a first attempt at generating a wcs image. It is strongly recommended (i.e., you need to do this) to run scamp or some other kind of astrometric correction software. Users have been fairly warned.  


## How to use
Simply download the three python scripts and use the "add_wcs_to_fits" function in the imacs_wcs.py script. 
eg: 

```python
infile = 'ift1067c4.fits'
add_wcs_to_fits(infile)
```

You can add an output name if you wish but by default the file with the outputted name will be appended with a .wcs. i.e., ift1067c4.fits becomes ift1057c4.wcs.fits. 

If you want to update an entire folder you can: 

```python
import glob
files = glob.glob('*.fits')
for file in files:
  add_wcs_to_fits(file)
 ```
 
 ## Future
 
 If we can confirm that this method is working correctly then there are some things we would like to implement 
  - write tests
  - rewrite in a more pythonic way
  - publish as a python package

  
 It is also possible that this may form part of a purely python-based pipeline for reducing IMACS data.
