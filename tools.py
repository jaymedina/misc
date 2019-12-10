import numpy as np
import os

from astropy.io import fits


def bpix_image(bpixtab, path='/grp/hst/wfc3j/jmedina/bpixtab_test/'):
  """ This function takes a UVIS bad pixel table and turns it into an image.

  Parameter(s)
  ------------
  bpixtab : str
    The bad pixel table you want to turn into an image.
    Format: `[filename]_bpx.fits`

  Returns
  -------
  A FITS file with the chip images of the unstable pixels on each UVIS CCD
  """
  # Getting data
  data = fits.getdata(bpixtab)

  # Getting the dimensions of one chip
  cols = 4096
  rows = 2051

  # Initializing an empty array that will be our image
  chip1 = np.zeros((rows, cols))
  chip2 = np.zeros((rows, cols))

  # Filling this blank image with the unstable pixel locations
  for chip, row, col, pixval in zip(data['CCDCHIP'], data['PIX2'], \
                                    data['PIX1'], data['VALUE']):
      if chip == 1:
          chip1[row-1, col-1] = pixval
      elif chip == 2:
          chip2[row-1, col-1] = pixval

  # Now let's throw it into a FITS file we can display on DS9
  hdu0 = fits.PrimaryHDU()
  hdu1 = fits.ImageHDU([chip2])
  hdu2 = fits.ImageHDU([chip1])

  hdulist = fits.HDUList([hdu0, hdu1, hdu2])

  rootname = bpixtab.split('_')[0]
  hdulist.writeto(os.path.join(path, rootname+'_img.fits'), overwrite=False)

def make_fits(array, filename, path=''):
    """ This function will put your array in a FITS file that you can open in
    DS9 for visual inspection, or any other purpose.

    Parameter(s)
    -----------
    array : array
        The 2-D array you'd like to put in a FITS file.
    filename : str
        The name of your file. Cannot include spaces.
    path : str
        The path where you want to save your array.
    """

    hdu0 = fits.PrimaryHDU([])
    hdu1 = fits.ImageHDU([array])

    hdulist = fits.HDUList([hdu0, hdu1])

    if path=='':
        path = os.getcwd()
        hdulist.writeto(path+filename+'.fits', overwrite=False)
    else:
        hdulist.writeto(path+filename+'.fits', overwrite=False)

def group_visits(wdir):
    """ Organizes the files in your working directory based on visit number.
    Generates a dictionary that sorts the files based on visit number.

    Parameters
    ---------
    wdir : str
        Your working directory

    Returns
    -------
    group : dict
        A dictionary whose keywords are the visit numbers, and the keyword values
        are the lists of files generated from that visit.
    """
    all_files = glob(os.path.join(wdir, '*flc.fits'))
    group = dict()
    for file in all_files:
        visit = fits.getheader(file)['LINENUM'].split('.')[0]
        if visit not in group:
            group[str(visit)] = [str(file)]
        elif visit in group:
            group[str(visit)].append(str(file))

    return group
