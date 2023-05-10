from distutils.core import setup

setup(
  name = 'imacs_wcs',
  packages = ['imacs_wcs'],
  version = '0.0.2',
  license='MIT',
  description = 'Library for determining and writing wcs information for IMACS images.',
  author = 'Trystan Scott Lambert',
  author_email = 'trystanscottlambert@gmail.com',
  url = 'https://github.com/TrystanScottLambert',
  download_url = 'https://github.com/TrystanScottLambert/imacs_wcs/archive/refs/tags/v0.0.2.tar.gz',
  keywords = ['astronomy', 'wcs', 'imacs'],
  install_requires=[
    'numpy',
    'photutils',
    'astropy',
    'scipy',
    'matplotlib',
    'astroquery',
    ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
  ],
)
