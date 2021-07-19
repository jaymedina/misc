"""
Functions in this module:
*bpix_image()
*make_fits()
*group_visits()
*get_slope
"""
import numpy as np
import os

from astropy.io import fits

def bpix_kw(bpixtab):
    """ This function prints out header keywords as part of BPIXTAB verification
    procedure.

    Parameter(s)
    ------------
    bpixtab : str
      The bad pixel table you want to turn into an image.
      Format: `[filename]_bpx.fits` or `path/to/[filename]_bpx.fits`

    """
    print('Verifying the header keywords of UVIS bad pixel table {}...'.format(bpixtab))
    print('USEAFTER:')
    print(fits.getheader(bpixtab)['USEAFTER'])
    print(' ')
    print('PEDIGREE:')
    print(fits.getheader(bpixtab)['PEDIGREE'])
    print(' ')
    print('DESCRIP:')
    print(fits.getheader(bpixtab)['DESCRIP'])
    print(' ')
    print('COMMENT:')
    print(fits.getheader(bpixtab)['COMMENT'])
    print(' ')
    print('HISTORY:')
    print(fits.getheader(bpixtab)['HISTORY'])

def bpix_image(bpixtab, path='/grp/hst/wfc3j/jmedina/bpixtab_test/'):
  """ This function takes a UVIS bad pixel table and turns it into an image.

  Parameter(s)
  ------------
  bpixtab : str
    The bad pixel table you want to turn into an image.
    Format: `[filename]_bpx.fits` or `path/to/[filename]_bpx.fits`

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

  if '/' in bpixtab:
      bpixtab = bpixtab.split('/')[-1]

  rootname = bpixtab.split('_')[0]
  hdulist.writeto(os.path.join(path, rootname+'_img.fits'), overwrite=False)
  print('New file created!')
  print('View it by entering `ds9 {}{}` in a fresh terminal.'.format(\
                                                            path, \
                                                            rootname+'_img.fits'))

def bpixtab_test(bpixtab, path='/grp/hst/wfc3j/jmedina/bpixtab_test/'):
    """ The main function for the UVIS bad pixel table verification procedure.
    It wraps the following functions:
    *bpix_kw
    *bpix_image

    Parameter(s)
    ------------
    bpixtab : str
      The bad pixel table you want to turn into an image.
      Format: `[filename]_bpx.fits` or `path/to/[filename]_bpx.fits`
    """
    # Verifying the header keywords
    bpix_kw(bpixtab)

    # Generating an image of the bad pixels using the bad pixel table
    # which can be inspected using DS9
    bpix_image(bpixtab, path) # uses default path

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

def get_slope(x, y, deg=1, err=[]):
    """ Calculate the slope of the polynomial that best fits your data.

    Parameters
    ----------
    x : arr, list
    y : arr, list
    deg : int
        The degree of your polynomial fit (default is deg=1, i.e. a line)

    Returns
    -------
    z and p

    """
    inverse_error = []
    for i in err:
        inv = 1/i
        inverse_error.append(i)

    if len(err)>0:
        z = np.polyfit(x, y, deg, w=inverse_error)
    else:
        z = np.polyfit(x, y, deg)

    m, b = z
    p = np.poly1d(z)

    return m, b, p
