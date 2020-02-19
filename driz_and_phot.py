""" Drizzle and photometry tools (wrappers) used for the WFC3/UVIS
time-dependent zeropoint calculations.

Functions
---------
- drizzle
- get_fluxes

Author(s)
---------
Jennifer V. Medina (2019)
"""

import os
import numpy as np

import matplotlib.pyplot as plt
import pysynphot as s
from astropy.io import fits
from drizzlepac import astrodrizzle
from glob import glob
from photutils.detection import DAOStarFinder
from phot_wrapper import uvis_photometry as photometry
from stsci.skypac import pamutils

def drizzle(input, work_dir):
    """ Drizzles the input using new drizzle parameters provided by Jennifer Mack.

    Parameters
    ----------
    input : FLC list or single ASN
        The input may be a list of FLCs, or an ASN FITS file.
    work_dir : str
        The directory in which you want to pull the ``input`` from, and where
        you want the drizzled product to reside.
    """
    os.cwd(work_dir)

    # If the input is an []_asn.fits file...
    if isinstance(input, str) == True:
        # getting the rootname for the ``output`` parameters
        name = input.split('_')[0]
        astrodrizzle.AstroDrizzle(input, output='new_'+name, skymethod='match',
                              driz_sep_bits='65535', combine_type='mindmed',
                              driz_cr_snr='4.5 4.0', final_bits='80',
                              driz_cr_scale='3.0 2.5', build=True)

    # otherwise, it must be a list of []_flc.fits files...
    else:
        for i in input:
            # getting the rootname for the ``output`` parameters
            name = input.split('_')[0]
            astrodrizzle.AstroDrizzle(i, output='new_'+name, skymethod='match',
                                  driz_sep_bits='65535', combine_type='mindmed',
                                  driz_cr_snr='4.5 4.0', final_bits='80',
                                  driz_cr_scale='3.0 2.5', build=True)

def get_fluxes(file_list, filt, targ, chip, fwhm=1.675, threshold=1000):
    """ Calculates the photometry of the correct sources and does a ratio.

    Parameters
    --------
    file_list : list
        List of files (1 drc file + n flc files in visit)
    filt : str
        Filter
    targ : str
        Target
    chip : int
        1 or 2
    fwhm : float
    threshold : float

    Returns
    -------
    filenames : list
        the list of files
    fluxes : list
        the list of photometric fluxes
    """

    # get bandpass
    bandpass = s.ObsBandpass('wfc3,uvis{},{},aper#0.3962'.format(1, filt))
    # get spectrum
    spectrum = s.FileSpectrum('/grp/hst/cdbs/calspec/{}_stisnic_006.fits'.format(targ))
    # get observation
    observations = s.Observation(spectrum, bandpass, force='extrap')
    # find the source(s)
    countrate = observations.countrate()

    dsf = DAOStarFinder(threshold, fwhm)
    # doing the calculations
    delta_dict = {}
    fwhms, xs, ys = [], [], []
    filenames, fluxes = [], []

    for file in file_list:

        # getting the right SCI ext
        if chip == 1:
            data = fits.getdata(file, 4)
        elif chip == 2:
            data = fits.getdata(file, 1)

        # correcting for distortion using PAM
        if 'flc' in file:

            if chip == 1:
                pamutils.pam_from_file(file, ext=4, output_pam='pam_chip_'+str(chip)+str(file))
                PAM = fits.getdata('pam'+str(file))
                data = data*PAM
            elif chip == 2:
                pamutils.pam_from_file(file, ext=1, output_pam='pam_chip_'+str(chip)+str(file))
                PAM = fits.getdata('pam_chip_'+str(chip)+str(file))
                data = data*PAM

        # continue the calculations
        src_tbl = dsf.find_stars(data)
        crd = np.array([src_tbl['xcentroid'], src_tbl['ycentroid']])
        phot_tbl = photometry(data=data, coords=crd)
        exptime = fits.getheader(file)['EXPTIME']

        if 'flc.fits' in file:
            flux = phot_tbl['flux']/exptime

        else:
            flux = phot_tbl['flux']

        phot_tbl['synratio'] = np.abs((flux/countrate) - 1.) # < - this should be zero
        phot_tbl.sort('synratio') # < - this basically puts the TRUE source (our white dwarf)
                                  # at the top of the list, so that we aren't
                                  # comparing: Cosmic Ray vs. White Dwarf
                                  # we're comapring: WD vs. WD

        print('Almost done. Some quality checks...')
        print('Making sure the order is what I expect:')
        print(phot_tbl['synratio'][0:10])
        print('This is the synratio Obs/Syn (should be close to 0.)')
        print(phot_tbl['synratio'][0])
        print('This is the flux of the EE that was chosen, and the file')
        print(phot_tbl['flux'][0], file)
        print('These are the coordinates of the EE used to calculate the flux')
        print('(X, Y)')
        print(phot_tbl['X'][0], phot_tbl['Y'][0])

        filenames.append(file)
        if 'flc.fits' in file:
            fluxes.append(phot_tbl['flux'][0]/exptime)

        else:
            fluxes.append(phot_tbl['flux'][0])

    return filenames, fluxes
