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

"""
## -- IMPORTS
import glob
from multiprocessing import Pool
import os
import shutil
import socket
import sys
import warnings

from astropy.io import fits
#from astropy.utils.exceptions import AstropyUserWarning
#from bs4 import BeautifulSoup
#import ftplib
#from ftplib import FTP
#from getpass import getpass
import itertools
from urllib.request import urlopen

from pyql.database.ql_database_interface import Master
from pyql.database.ql_database_interface import session
from pyql.download.download_tools import get_partition
from pyql.download.download_tools import build_unlooked_directory
from pyql.file_system.set_permissions import set_permissions

## -- FUNCTIONS

def main(cores=30, mp_sort=True):
    """
    Main function to retrieve the datasets from staging
    and move them to the unlooked directory structure.
    Parameters
    ----------
    cores : int, optional
        The number of corse to parallel procces with.
    mp_sort : bool, optional
        Whether or not to multiprocess the data shuffling.
    """

    user = input('Please input your MAST user -- WITHOUT @stsci.edu : \n')
    password = getpass()
    collect_data_from_staging(user, password)

    p = Pool(cores)
    # Collect the list of files newly downloaded
    # And send them to their final destination
    unlooked_files = glob.glob('/grp/hst/wfc3a/unlooked_temp/*.fits')
    print('{} files will be moved.'.format(len(unlooked_files)))

    if mp_sort:
        p.map(sort_data, unlooked_files)
    else:
        for unlooked in unlooked_files:
            sort_data(unlooked)


def collect_data_from_staging(user, password, dest='/grp/hst/wfc3j/jmedina/', cores=30):
    """
    Writes the data from the MAST ftp stagin area to a
    temporary unlooked directory (or one of your choosing.)
    Parameters
    ----------
    user : str
        the user that made the MAST request.
    password : str
        The AD password.
    dest : str, optional
        The path to write the data.
    cores : int
        The number of cores to use when parallel.
    """

    try:
        # Initialize ftp connection with the archive
        ftp = FTP('stdatu.stsci.edu')
        ftp.login(user='{}@stsci.edu'.format(user), passwd=password)
        ftp.cwd('/stage/{}/'.format(user))

        # Traverse the structure of every ftp request
        dir_info = []
        ftp.dir(dir_info.append)
        dirs = [user + dir_dat.split(user)[-1] for dir_dat in dir_info]

        # Collect every file in staging
        files = []
        for data_dir in dirs:
            ftp.cwd(data_dir)
            temp_files = ftp.nlst()
            files += [data_dir + '/' + temp_file for temp_file in temp_files]
            ftp.cwd('/stage/{}/'.format(user))
        ftp.quit()

        # Collect a list of every file in the ql_db or unlooked_temp
        ql_results = session.query(Master.ql_root).all()
        ql_db_files = [result[0] for result in ql_results]
        present_files = [present_file.split('/')[-1][:8] for present_file in glob.glob(dest + '*')]

        # Remove anything already in Quicklook or unlooked_temp
        new_files = []
        for f in files:
            filename = f.split('/')[-1][:8]
            if filename not in ql_db_files + present_files and '.fits' in f and filename[0] == 'i':
                new_files.append(f)
        print('{} files will be downloaded.'.format(len(new_files)))

        # Copy them to unlooked_temp
        # Parallel processing engaged
        p = Pool(cores)
        if len(new_files) > 0:
            mapping_args = zip(new_files, itertools.repeat(user), itertools.repeat(password), itertools.repeat(dest))
            p.map(download_data, mapping_args)

    # If the ftp lib throws a password error restart the func.
    except ftplib.error_perm as e:
        if 'No such file or directory' in str(e):
            print('There are no files for you in staging.')

        else:
            print('Your user/password are wrong. Starting over...')
            main()


def download_data(mapped_args):
    """
    Opens up a NEW AND EXCITING connection to the MAST ftp
    stagin' area to download one file. Optimized for
    parallel momes.
    Parameters
    ----------
    mapped_args : tup
        Tuple containing (staged file, MAST user, MAST
        password, and destination to write the file).
    """

    # Unpack zipped args
    staged_file, user, password, dest = mapped_args

    try:

        # Initialize ftp connection
        ftp = FTP('stdatu.stsci.edu')
        ftp.login(user='{}@stsci.edu'.format(user), passwd=password)
        ftp.cwd('/stage/{0}/{1}'.format(user, staged_file.split('/')[-2]))

        # Write new file to dest
        outfile = open(dest + staged_file.split('/')[-1], 'wb')
        ftp.retrbinary('RETR ' + staged_file.split('/')[-1], outfile.write)
        outfile.close()
        ftp.quit()

        new_file_path = dest + staged_file.split('/')[-1]

        if test_buffer(new_file_path):
            print('Downloaded file from staging area to {}.'.format(dest + staged_file.split('/')[-1]))

        else:
            print("The file was corrupted in the process. We'll try again.")
            print('Restart download was triggered.')
            download_data(mapped_args)

    # These errors indicate the file might have been
    # corrupted while downloading and should be done again
    # OR MAST got spooped by so many open connections and
    # shut it down.
    except (ftplib.error_temp, OSError) as e:
        print('Restart download was triggered.')
        download_data(mapped_args)


def find_proposal_id(infile):
    """ Finds the proposal ID of a given proposal.
    Matches to it or another file's FITS header.
    Parameters
    ----------
    infile : str
        The file path.
    Returns
    -------
    prop_id : str of int
        The proposal ID.
    """

    prop_id = True
    try:
        # Find a matching file if it's a trl or hlet
        if ('trl.fits' in infile) or ('hlet.fits' in infile):
            root = infile.split('/')[-1].split('_')[0]
            files = glob.glob('/'.join(infile.split('/')[:-1]) + '/{}*raw.fits'.format(root))
            ql_results = session.query(Master.dir).filter(Master.ql_root == root[:-1]).all()

            if len(files) > 0:
                infile = files[0]

            elif len(ql_results) > 0:
                file_path = ql_results[0][0]
                prop_id = file_path.split('/')[5]
                return prop_id

            else:
                prop_id = False

        if prop_id:
            with fits.open(infile) as hdu:
                hdr = hdu[0].header.copy()
                prop_id = hdr['proposid']

    except (KeyError, TypeError, OSError) as e:
        print('Something is still wrong with this file : {}.'.format(infile))
        os.remove(infile)
        print('{} has been removed. Try downloading it again later...'.format(infile))
        prop_id = False

    return prop_id


def find_proposal_partition(prop):
    """ Finds the proposal type of a given proposal.
    Scrapes the proposal info page.
    Parameters
    ----------
    prop : str of int
        The proposal ID.
    Returns
    -------
    prop_partition : str
        'wfc3a' or 'wfc3f' depending on which type.
    """

    # Open up the proposal status page
    url = 'http://www.stsci.edu/cgi-bin/get-proposal-info?id={}&submit=Go&observatory=HST'.format(prop)
    html = urlopen(url).read()
    soup = BeautifulSoup(html, 'lxml').findAll('a')

    prop_type = 'web error'

    # Scoop up the proposal type
    for turtle in soup:
        if 'prop_type' in str(turtle):
            prop_type = str(turtle).split('prop_type">')[-1].split('</a>')[0]

    if 'CAL' in prop_type:
        prop_partition = 'wfc3a'
    else:
        prop_parition = 'wfc3f'

    return prop_partition


def sort_data(infile):
    """
    Moves a file to its place in the unlooked structure.
    Parameters
    ----------
    file : str
        Full file path.
    """

    print('Working on file {}'.format(infile))
    prop = find_proposal_id(infile)

    if prop == False:
        print('This awkward file-type -- {} -- has no companion just yet.'.format(infile))

    else:
        prop_partition = get_partition(prop)
        visit = infile.split('/')[-1][4:6].upper()

        visit_path = '/grp/hst/{0}/unlooked/{1}/Visit{2}/'.format(prop_partition, prop, visit)
        new_path = '/grp/hst/{0}/unlooked/{1}/Visit{2}/{3}'.format(prop_partition, prop, visit, infile.split('/')[-1])
        build_unlooked_directory(visit_path)
        shutil.move(infile, new_path)
        set_permissions(new_path)
        print('File moved from {0} to {1}.'.format(infile, new_path))


def test_buffer(infile):
    """
    Tests for a whole and complete fits file buffer.
    An incomplete buffer usually signals and incomplete
    file and will break things down the line. This just
    excepts the AstropyUserWarning.
    Parameters
    ----------
    infile : str
        The path of the file to test.
    Returns
    -------
    complete_file : bool
        Wether or not the file is okay.
    """

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        try:
            with fits.open(infile) as hdu:
                hdu.verify()
                complete_file = True
        except AstropyUserWarning:
            complete_file = False

    return complete_file


## -- RUN

if __name__ == "__main__":
    main()
"""
