"""
GAIA Alignment
"""

from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from astroquery.gaia import Gaia


def search_gaia_archives(ra: float, dec: float, height:float, width:float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Does a box search centered on ra and dec with a given height and width.
    """
    coord = SkyCoord(ra = ra*u.deg, dec=dec*u.deg)
    width = u.Quantity(width*u.deg)
    height = u.Quantity(height*u.deg)
    results = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    return np.array(list(results['ra'])), np.array(list(results['dec']))

def get_fits_center(hdu) -> Tuple[float, float]:
    """
    Determines the centeral pixels of the image.
    """
    local_wcs = WCS(hdu.header)
    y,x = hdu.shape
    ra_center, dec_center = local_wcs.pixel_to_world_values(x/2,y/2)
    return float(ra_center), float(dec_center)

def get_fits_dimensions(hdu):
    """
    Calculates the width and height of a fits image
    """
    local_wcs = WCS(hdu.header)
    pixHeight, pixWidth = hdu.data.shape
    x_pixels = np.array([0,0,pixWidth,pixWidth])
    y_pixels = np.array([0,pixHeight,0,pixHeight])
    ra_pixels,dec_pixels = local_wcs.pixel_to_world_values(x_pixels,y_pixels)
    return np.abs(ra_pixels[0] - ra_pixels[2]), np.abs(dec_pixels[0]- dec_pixels[1])

def locateImageHDUValue(hdulist):
    for i in range(len(hdulist)):
        try:
            hdulist[i].header['CRVAL1'] = hdulist[i].header['CRVAL1']
            value = i
        except:
            pass
    return value

def splitTupples(listOfTupples):
	firstElements = [tupple[0] for tupple in listOfTupples]
	secondElements = [tupple[1] for tupple in listOfTupples]
	return firstElements, secondElements

def makeTupples(firstElementsArray, secondElementsArray):
	tupples = np.array([(firstElementsArray[i], secondElementsArray[i]) for i in range(len(firstElementsArray))])
	return tupples

def sumColumnsAndRows(array2D):
	columnSum = np.sum(array2D,axis=0) # axis 0 gives columns  
	rowSum = np.sum(array2D,axis=1)
	return columnSum, rowSum


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gaussFit(x, y):
	mean = sum(x * y) / sum(y)
	sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
	popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
	return popt

def cutStar(xPosition,yPosition,imageData):
	xPositionRounded, yPositionRounded = int(round(xPosition)), int(round(yPosition))
	padding = 50
	imageDataPostageStamp = imageData[yPositionRounded-padding:yPositionRounded+padding,xPositionRounded-padding:xPositionRounded+padding]


	localYPosition, localXPosition = getStarPosition(imageDataPostageStamp)
	accurateXPosition, accurateYPosition = xPositionRounded-padding + localYPosition,yPositionRounded-padding + localXPosition

	return imageDataPostageStamp, accurateYPosition, accurateXPosition

def getStarPosition(star):
	columnSum, rowSum = sumColumnsAndRows(star)
	xValues = np.arange(len(columnSum))
	columnPopt = gaussFit(xValues,columnSum)
	rowPopt = gaussFit(xValues,rowSum)
	return columnPopt[2], rowPopt[2]


def prettifyPlot(xlabel,ylabel):
	plt.xlabel(xlabel,fontsize=30)
	plt.ylabel(ylabel,fontsize=30)
	plt.tick_params(axis='both',labelsize=20)
	plt.minorticks_on()
	plt.tick_params(which='both', width=2,direction='in') #this has to be a separate line because width can't be set when using which='minor'
	plt.tick_params(which='major', length=8, direction='in') #If you want tick marks on the opposite side also, add right=True
	plt.tick_params(which='minor', length=4)



class GaiaAlignment:
	def __init__(self, inputFile):
		self.hdulist = fits.open(inputFile)
		hduValue = locateImageHDUValue(self.hdulist)
		self.hdu = self.hdulist[hduValue]
		self.hduHeader = self.hdu.header 
		self.hduData = self.hdu.data 
		self.wcs = WCS(self.hduHeader)


		centerRa,centerDec = get_fits_center(self.hdu)
		width,height = get_fits_dimensions(self.hdu)
		self.gaiaRA, self.gaiaDec = search_gaia_archives(centerRa,centerDec,height,width)
		self.gaiaX, self.gaiaY = self.wcs.world_to_pixel_values(self.gaiaRA,self.gaiaDec)

		self.plotImage()
		self.generateWorldCoordinates()
		self.calculateXYOffset()


	def plotImage(self):
		self.coords = []
		fig = plt.figure()
		interval = ZScaleInterval()

		ax = fig.add_subplot()
		limits = interval.get_limits(self.hduData)	
		ax.imshow(self.hduData, cmap='Greys', origin='lower',interpolation='nearest',vmin=limits[0],vmax=limits[1])
		ax.scatter(self.gaiaX,self.gaiaY,marker='*',s=100,facecolors='none',edgecolors='cyan')

		def onclick(event):
			if event.dblclick:
				star,columnPix, rowPix = cutStar(event.xdata,event.ydata,self.hduData)
				if np.abs(rowPix-event.xdata) < 15 and np.abs(columnPix-event.ydata) < 15:
					print(f'x = {event.xdata}, y = {event.ydata}')
					plt.scatter(rowPix, columnPix,facecolors='none',edgecolors='r',s=50)
					fig.canvas.draw()
					self.coords.append((rowPix,columnPix))
		cid = fig.canvas.mpl_connect('button_press_event', onclick)
		plt.show()

	def generateWorldCoordinates(self):
		x ,y = splitTupples(self.coords)
		self.starWorldCoordinates = self.wcs.pixel_to_world_values(x,y)

	def calculateXYOffset(self):
		counter = 0
		raOffsets = []
		decOffsets = []
		print(f'\t Gaia_Star_RA, Gaia_Star_Dec, HST_Star_RA, HST_Star_Dec, RA_Offset, Dec_Offset')
		for i in range(len(self.coords)):
			distances = np.sqrt((self.coords[i][0] - self.gaiaX)**2 + (self.coords[i][1] - self.gaiaY)**2)
			idxMatch = np.where(distances == np.min(distances))[0]

			matchingGaiaRA, matchingGaiaDec = self.gaiaRA[idxMatch], self.gaiaDec[idxMatch]
			currentRAOffset, currentDecOffset = matchingGaiaRA - self.starWorldCoordinates[0][i], matchingGaiaDec - self.starWorldCoordinates[1][i]
			raOffsets.append(currentRAOffset)
			decOffsets.append(currentDecOffset)

			print(f'Gaia Star {i+1}: {matchingGaiaRA}, {matchingGaiaDec}, {self.starWorldCoordinates[0][i]}, {self.starWorldCoordinates[1][i]}, {currentRAOffset}, {currentDecOffset}')

		self.raOffset = np.mean(raOffsets)
		self.decOffset = np.mean(decOffsets)	
		self.raSTD = np.std(np.abs(raOffsets))
		self.decSTD = np.std(np.abs(decOffsets))

		print(f' \n\n TOTAL OFFSET: \n RA: {self.raOffset} +- {self.raSTD} \n Dec: {self.decOffset} +- {self.decSTD}')

	def applySimpleCorrection(self,outputFileName):
		self.hduHeader['CRVAL1'] = self.hduHeader['CRVAL1'] + self.raOffset
		self.hduHeader['CRVAL2'] = self.hduHeader['CRVAL2'] + self.decOffset
		self.hdulist.writeto(outputFileName,overwrite=True)


if __name__ == '__main__':
	infile1 = 'data/ASCImages/raw/hst_13641_07_wfc3_ir_f105w_sci.fits'
	infile2 = 'data/ASCImages/raw/hst_13641_07_wfc3_ir_f125w_sci.fits'
	infile3 = 'data/ASCImages/raw/hst_13641_07_wfc3_ir_f160w_sci.fits'
	infile4 = 'data/ASCImages/raw/hst_13641_07_wfc3_ir_total_sci.fits'
	infiles = [infile1,infile2,infile3,infile4]

	for infile in infiles:
		ga = GaiaAlignment(infile)
		#ga.applySimpleCorrection('data/ASCImages/gaiacorrected/' + infile.split('/')[-1].split('.fit')[0] + '_gaia_corrected.fits')
Footer