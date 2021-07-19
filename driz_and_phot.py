""" Drizzle and photometry tools (wrappers) used for the WFC3/UVIS
time-dependent zeropoint calculations.

Functions
---------
- drizzle
- get_fluxes
- wcs_hack
- obs_syn_ratio

Author(s)
---------
Jennifer V. Medina (2019)
"""

import os
import numpy as np

import matplotlib.pyplot as plt
import pysynphot as s
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
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

def wcs_hack(input_flcs, input_target, input_filter):
    """ This modified the WCS for all the input FLCs so that their rotations
    are matched. Initially this was performed so we would be able to drizzle
    GRW +70 FLCs that were taken in different visits. Now this is used to
    before drzzling all GRW +70 FLCs to create a `mega-drizzle` which is
    then used to calculate the Encircled Energy for a given filter at
    `infinity` (or 6.0 arcseconds, for WFC3 UVIS). This function keeps the
    rotations matched so that all the PSF wings are aligned before drizzling.

    What this function does:
    1. Arbitrarily selects an FLC whose WCS will be used as the `reference WCS`.
       This is the WCS whose CD matrix will be copied. Each amp should have its
       own `reference FLC`.

    Next, for each FLC...
    > 2. Sets the reference pixels CRPIX(1/2) to the star's centroid (X,Y) which
         is pulled from Clare's photometry catalogs.
    > 3. Sets CRVAL(1/2) for every FLC to the RA,DEC of the target.
    > 4. Sets the CD matrix terms CD(1_1, 1_2, 2_1, 2_2) for all the FLCs, to
         the same ones from the `reference FLC`. This is so all the images have
         the same rotation.


    Parameter(s)
    -----------
    input_flcs : list
        The list of FLCs you want to drizzle together. *IMPORTANT*: The WCS
        of these FLCs will be permanently modified.
    input_target : str
        The target in the FLCs.
    input_filter : str
        The filter observed with.
    """
    # Use one set of values for all the GRW images, you can use these, or the simbad ones or whatever
    # these are just some value near the star.
    if input_target=='GRW70':
        ra = 204.7063652201
        dec = 70.285363975590002
    elif input_target=='GD153':
        ra = 194.2583
        dec = 22.0314
    elif input_target=='P330E':
        ra, dec = 247.8891622, 30.14621295
    elif input_target=='G191B2B':
        ra, dec = 76.37752048, 52.83078003

    # path to Clare's tables
    path_to_catalog = '/grp/hst/wfc3p/cshanahan/phot_group_work/'+\
                      'analysis_pipelines/standard_star_photometry_pipeline/output/photcats/'+\
                      '{}_{}.dat'.format(input_target, input_filter)

    tbl = Table.read(path_to_catalog, format='ascii')


    # **Each amp needs it's own reference WCS**, I just chose the first of the images for that amp
    # Separate images by target/amp, in this case I only have amp A

    hdu = fits.open(input_flcs[0])

    # Make reference WCS object
    w1 = WCS(hdu[1].header, hdu) # get the WCS from which we'll copy the CD matrix.  Do this for each amp independently
    hdu.close()

    bad_flcs = []
    for im in input_flcs:
        row = [row for row in tbl if row['rootname'] in im]

        if (len(row) == 0):
            print('this FLC is not in Clare`s catalog. we won`t be drizzling for now')
            print(os.path.basename(im))
            bad_flcs.append(im)
            continue

        row = row[0]

        x, y = row['xcenter'], row['ycenter']

        hdu = fits.open(im, mode='update')
        hdu[1].header['CRPIX1'] = x # Set the reference pixel to the star's centroid
        hdu[1].header['CRPIX2'] = y
        hdu[1].header['CRVAL1'] = ra # set the RA/DEC to the same thing for all the images
        hdu[1].header['CRVAL2'] = dec
        hdu[1].header['CD1_1'] = w1.wcs.cd[0,0] # Copy the CD matrix terms to make all the images have the same rotation
        hdu[1].header['CD1_2'] = w1.wcs.cd[0,1]
        hdu[1].header['CD2_1'] = w1.wcs.cd[1,0]
        hdu[1].header['CD2_2'] = w1.wcs.cd[1,1]
        hdu[1].header['WCSNAME'] = 'JANK' # lol
        hdu.close()

    return bad_flcs

def obs_syn_ratio(flux, target, input_filter):

    target = 'grw_70d5824' if target=='GRW70' else target
    # get bandpass
    bandpass = s.ObsBandpass('wfc3,uvis1,{},aper#0.3962'.format(input_filter.lower()))
    # get spectrum
    spectrum = s.FileSpectrum('/grp/hst/cdbs/calspec/{}_stisnic_006.fits'.format(target.lower()))
    # get observation
    observations = s.Observation(spectrum, bandpass, force='extrap')
    # get synthetic flux
    countrate = observations.countrate()
    # make the calculations
    ratio = flux/countrate

    return ratio
