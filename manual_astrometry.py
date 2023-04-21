"""
Module to overplot gaia stars on an image and then use them to
determine a WCS.
"""

from typing import Tuple
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.wcs import utils
import astropy.units as u
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia


def search_gaia_archives(
        ra: float, dec: float, height:float, width:float)\
              -> Tuple[np.ndarray, np.ndarray]:
    """
    Does a box search centered on ra and dec with a given height and width.
    """
    coord = SkyCoord(ra = ra*u.deg, dec=dec*u.deg)
    width = u.Quantity(width*u.deg)
    height = u.Quantity(height*u.deg)
    results = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    return np.array(list(results['ra'])), np.array(list(results['dec']))

def sum_columns_and_rows(array_2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """ Sums up the rows and the columns so that we can fit gaussians to them."""
    column_sum = np.sum(array_2d,axis=0) # axis 0 gives columnsaccurate_y_position
    row_sum = np.sum(array_2d,axis=1)
    return column_sum, row_sum

def gauss(x_array: np.ndarray, constant:float, amplitude:float, mean:float, sigma:float):
    """Definiton of a gaussian."""
    return constant + amplitude * np.exp(-(x_array - mean) ** 2 / (2 * sigma ** 2))

def gauss_fit(x_array: np.ndarray, y_array:np.ndarray) -> np.ndarray:
    """Fit a 1-D gaussian to x and y data."""
    mean = np.sum(x_array * y_array) / np.sum(y_array)
    sigma = np.sqrt(np.sum(y_array * (x_array - mean) ** 2) / np.sum(y_array))
    popt, _ = curve_fit(gauss, x_array, y_array, p0=[np.min(y_array), np.max(y_array), mean, sigma])
    return popt

def get_star_position(star):
    """Determines the two 1-d gaussian fit positoin of the cut out star."""
    column_sum, row_sum = sum_columns_and_rows(star)
    x_values = np.arange(len(row_sum))
    row_popt = gauss_fit(x_values,row_sum)
    column_popt = gauss_fit(x_values,column_sum)
    return column_popt[2], row_popt[2]

def cut_star(x_position,y_position,image_data):
    """Cuts out star and determines the center of said star."""
    x_position_rounded, y_position_rounded = int(round(x_position)), int(round(y_position))
    padding = 20
    image_postage_stamp = image_data[
         y_position_rounded-padding:y_position_rounded+padding,
         x_position_rounded-padding:x_position_rounded+padding]

    if 0 in image_postage_stamp.shape:
        return None

    try:
        local_x_position, local_y_position = get_star_position(image_postage_stamp)
        accurate_x_position = x_position_rounded - padding + local_x_position
        accurate_y_position = y_position_rounded - padding + local_y_position
        return image_postage_stamp, accurate_x_position, accurate_y_position
    except (RuntimeError, ValueError):
        return None


class DraggableScatter():

    epsilon = 5

    def __init__(self, scatter):
        "Initiliazing."
        self.scatter = scatter
        self._ind = None
        self.ax = scatter.axes
        self.canvas = self.ax.figure.canvas
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        plt.show()


    def get_ind_under_point(self, event):
        "Idenifying the scatter point we are clickin."
        xy = np.asarray(self.scatter.get_offsets())
        xyt = self.ax.transData.transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]

        d = np.sqrt((xt - event.x)**2 + (yt - event.y)**2)
        ind = d.argmin()

        if d[ind] >= self.epsilon:
            ind = None
        return ind, xy[ind]

    def button_press_callback(self, event):
        """Clicking the button"""
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind, pos = self.get_ind_under_point(event)
        self.x_ref = pos[0]
        self.y_ref = pos[1]


    def button_release_callback(self, event):
        """Letting go of the button"""
        if event.button != 1:
            return
        self._ind = None
        
        xy = np.asarray(self.scatter.get_offsets())
        x_offset = event.xdata - self.x_ref
        y_offset = event.ydata - self.y_ref

        #xy[self._ind] = np.array([x, y])
        for idx, _ in enumerate(xy):
            xy[idx] = np.array([xy[idx][0]+x_offset, xy[idx][1]+y_offset])

        self.scatter.set_offsets(xy)
        self.canvas.draw_idle()
        self.positions = xy


    def motion_notify_callback(self, event):
        """Moving the mouse"""
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return


def get_usable_gaia(gaia_x_array: np.ndarray, gaia_y_array: np.ndarray, image: np.ndarray):
    """
    Searches the positions of the x and y coordinates for
    a bright star that can be fit.
    """

    accurate_x = []
    accurate_y = []
    indicies = []
    for i, _ in enumerate(gaia_x_array):
        star_cutout = cut_star(gaia_x_array[i], gaia_y_array[i], image)
        if star_cutout is not None:    
            _, x_point, y_point = star_cutout
            accurate_x.append(x_point)
            accurate_y.append(y_point)
            indicies.append(i)

    return np.array(accurate_x), np.array(accurate_y), np.array(indicies)


class ChipImage:
    """Main class for chip data."""

    def __init__(self, chip_name:str):
        """Initilizing."""
        self.file_name = chip_name
        self.hdul = fits.open(chip_name)
        self.data = self.hdul[0].data
        self.header = self.hdul[0].header
        self.current_wcs = WCS(self.header)
        current_gaia_info = self.query_gaia()
        self.gaia_ra = current_gaia_info[0]
        self.gaia_dec = current_gaia_info[1]
        self.current_gaia_x = current_gaia_info[2]
        self.current_gaia_y = current_gaia_info[3]
    
    def query_gaia(self):
        """Gets the gaia positions in the chip frame."""
        center_y_pix = self.data.shape[0]/2
        center_x_pix = self.data.shape[1]/2

        # Assuming square pixels.
        search_width = self.data.shape[1] * np.abs(self.header['PC1_1'])
        search_height = self.data.shape[0] * np.abs(self.header['PC1_1'])

        center_ra, center_dec = self.current_wcs.pixel_to_world_values(center_x_pix, center_y_pix)
        gaia_ra, gaia_dec = search_gaia_archives(center_ra, center_dec, width = search_width, height=search_height)
        gaia_x, gaia_y = self.current_wcs.world_to_pixel_values(gaia_ra, gaia_dec)

        return gaia_ra, gaia_dec, gaia_x, gaia_y
    
    def manually_determine_wcs_with_gaia(self, show_alignment=False):
        """Starts the draggable plot for alignment."""

        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(self.data)

        fig = plt.figure()
        ax = fig.add_subplot(projection=self.current_wcs)
        ax.imshow(self.data, vmin=vmin, vmax=vmax, cmap='gray')
        scatter = ax.scatter(
            self.current_gaia_x, self.current_gaia_y, facecolor='None', edgecolor='r', marker='s', s=80, picker=True)
        DraggableScatter(scatter)
        plt.show()

        updated_gaia_x, updated_gaia_y = scatter.get_offsets().T
        accurate_gaia_x, accurate_gaia_y, msk = get_usable_gaia(updated_gaia_x, updated_gaia_y, self.data)
        diff = np.hypot(updated_gaia_x[msk] - accurate_gaia_x, updated_gaia_y[msk] - accurate_gaia_y)
        _, diff_median, diff_std = sigma_clipped_stats(diff)
        cut = np.where(diff < diff_median + 1*diff_std)[0]

        #print('Difference in x: ', np.mean(accurate_gaia_x[cut] - self.current_gaia_x[msk][cut]), np.std(accurate_gaia_x[cut] - self.current_gaia_x[msk][cut]))
        #print('Difference in y: ', np.mean(accurate_gaia_y[cut] - self.current_gaia_y[msk][cut]), np.std(accurate_gaia_y[cut] - self.current_gaia_y[msk][cut]))

        gaia_coords = [(self.gaia_ra[msk][cut][i], self.gaia_dec[msk][cut][i]) for i, _ in enumerate(self.gaia_ra[msk][cut])]
        gaia_skycoords = SkyCoord(gaia_coords, unit=(u.deg, u.deg))
        gaia_pixcoords = np.array([accurate_gaia_x[cut], accurate_gaia_y[cut]])

        wcs = utils.fit_wcs_from_points(gaia_pixcoords, gaia_skycoords)

        if show_alignment:

            fig = plt.figure()
            ax = fig.add_subplot(projection=self.current_wcs)
            ax.imshow(self.data, vmin=vmin, vmax=vmax, cmap='gray')
            ax.scatter(accurate_gaia_x[cut], accurate_gaia_y[cut], facecolor='None', edgecolor='r', marker='s', s=100, picker=True)
            plt.show()

        return wcs

    def manually_update_wcs_with_gaia(self, show_alignment=False):
        """Aligns the image with gaia and determines accurate WCS information."""
        wcs = self.manually_determine_wcs_with_gaia(show_alignment)
        self.hdul[0].header.update(wcs.to_header())
        self.hdul.writeto(self.file_name.split('.fits')[0] + '.wcs_aligned.fits', overwrite=True)

if __name__ == '__main__':
    files = glob.glob('/home/tlambert/Downloads/g_band/RAW_SCIENCE/*c1*.wcs.fits')
    done_files = glob.glob('/home/tlambert/Downloads/g_band/RAW_SCIENCE/*c1*.wcs_aligned*')
    done_file_originals = [file.replace('.wcs_aligned','') for file in done_files]
    to_be_done_files = np.setdiff1d(files, done_file_originals)
    for file in to_be_done_files:
        chip = ChipImage(file)
        chip.manually_update_wcs_with_gaia()
