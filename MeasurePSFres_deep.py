#!/usr/bin/env python
# coding: utf-8

# In[17]:


#Make a catalogue with real and model PSFs + magnitude etc, for PSF testing script

#! /usr/bin/env python

from __future__ import print_function
import os
import numpy as np
#from read_psf_cats import read_data, band_combinations
import fitsio
import treecorr
import matplotlib
import matplotlib
matplotlib.use('Agg') # needs to be done before import pyplot
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table,join
import h5py as h

from __future__ import print_function
import os
import sys
import shutil
import logging
import datetime
import traceback
import numpy as np
import copy
import glob
import time
import fitsio
#import pixmappy
import pandas
#import galsim
#import galsim.des
#import piff
import ngmix
import wget


# In[14]:


#read in list of stars made from Sextractor and PSFEx
star_file= "/global/homes/a/aamon/DES/DEStests/DEEP/deeppsfs/UltraVista/UVISTA_J_21_01_16_psfex-starlist.fits"

dat = fits.open(star_file)
cols = dat[2].columns
print(cols)

  # This has the following columns:
    # id: The original id from the SExtractor catalog
    # x: The x position
    # y: The y position
    # sky: The local sky value
    # noise: The estimated noise.  But these are all 0, so I think this isn't being calculated.
    # size_flags: Error flags that occurred when estimating the size
    # mag: The magnitude from SExtractor
    # sg: SExtractor's star/galaxy estimate.  Currently SPREAD_MODEL.  (Actually, currently none)
    # sigma0: The shapelet sigma that results in a b_11 = 0 shapelet parameter.
    # star_flag: 1 if findstars thought this was a star, 0 otherwise.


# In[12]:


# Download the files we need:

image_file = wget('ftp://ftp.star.ucl.ac.uk/whartley/ultraVISTA/UVISTA_J_21_01_16_allpaw_skysub_015_dr3_rc_v5.fits.gz')
row['root'] = root
row['image_file'] = image_file

#usually weight is in image file but in this case, it's a separate file
weight_file = wget('ftp://ftp.star.ucl.ac.uk/whartley/ultraVISTA/UVISTA_J_21_01_16_allpaw_skysub_015_dr3_rc_v5.weight.fits.gz')


# In[ ]:


#Not sure this is necessary, but having this information might be useful for further tests

def read_image_header(row, img_file):
    """Read some information from the image header and write into the df row.
    """
    hdu = 0

    # Note: The next line usually works, but fitsio doesn't support CONTINUE lines, which DES
    #       image headers sometimes include.
    #h = fitsio.read_header(img_file, hdu)
    # I don't  care about any of the lines the sometimes use CONITNUE (e.g. OBSERVER), so I
    # just remove them and make the header with the rest of the entries.
    f = fitsio.FITS(img_file)
    header_list = f[hdu].read_header_list()
    header_list = [ d for d in header_list if 'CONTINUE' not in d['name'] ]
    h = fitsio.FITSHDR(header_list)
    try:
        date = h['DATE-OBS']
        date, time = date.strip().split('T',1)

        filter = h['FILTER']
        filter = filter.split()[0]

        sat = h['SATURATE']
        fwhm = h['FWHM']

        ccdnum = int(h['CCDNUM'])
        detpos = h['DETPOS'].strip()

        telra = h['TELRA']
        teldec = h['TELDEC']
        telha = h['HA']
        if galsim.__version__ >= '1.5.1':
            telra = galsim.Angle.from_hms(telra) / galsim.degrees
            teldec = galsim.Angle.from_dms(teldec) / galsim.degrees
            telha = galsim.Angle.from_hms(telha) / galsim.degrees
        else:
            telra = galsim.HMS_Angle(telra) / galsim.degrees
            teldec = galsim.DMS_Angle(teldec) / galsim.degrees
            telha = galsim.HMS_Angle(telha) / galsim.degrees

        airmass = float(h.get('AIRMASS',-999))
        sky = float(h.get('SKYBRITE',-999))
        sigsky = float(h.get('SKYSIGMA',-999))

        tiling = int(h.get('TILING',0))
        hex = int(h.get('HEX',0))

    except Exception as e:
        logger.info("Caught %s",e)
        logger.info("Cannot read header information from %s", img_file)
        raise

    row['date'] = date
    row['time'] = time
    row['sat'] = sat
    row['fits_filter'] = filter
    row['fits_fwhm'] = fwhm
    row['fits_ccdnum'] = ccdnum
    row['telra'] = telra
    row['teldec'] = teldec
    row['telha'] = telha
    row['airmass'] = airmass
    row['sky'] = sky
    row['sigsky'] = sigsky
    row['tiling'] = tiling
    row['hex'] = hex


# In[ ]:


read_image_header(row, image_file)


# In[9]:


#put the stars data into a dataframe 

def read_findstars(star_file, img_file):
    """Read the findstars output file
    """
    if not os.path.exists(star_file):
        return None

    # Read the output and make a DataFrome with the contents *********something buggy here
    data = fitsio.read(star_file)
    data = data.astype(data.dtype.newbyteorder('='))
    print(data) 
    df = pandas.DataFrame(data)
    print(df)
    ntot = len(df)
    ######nstars = df['star_flag'].sum()
    #logger.info('   found %d stars',ntot)
    print('   found %d stars',ntot)

    #print('mag range = ',np.min(df['mag']), np.max(df['mag']))
    #####is_star = df['star_flag'] == 1
    #print('star mag range = ',np.min(df['mag'][is_star]), np.max(df['mag'][is_star]))
    #print('zero point = ',magzp)
    #####df['mag'] += magzp - 25.
    #print('star mag range => ',np.min(df['mag'][is_star]), np.max(df['mag'][is_star]))

    #Add on some extra information from the sextractor catalog
    #INSTEAD I'LL USE THE WCS AND THE X,Y TO GET RA AND DEC.
    image = galsim.fits.read(img_file)
    wcs = image.wcs
    world = w.wcs_pix2world((x,y))
    print(world)
    df['ra'] = world[:,0]
    df['dec'] = world[:,1]
    print(df)
    return df


# In[10]:


df= read_findstars(star_file,image_file)


# In[31]:


#read in psf model file

psfex_file= "/global/homes/a/aamon/DES/DEStests/deeppsfs/UltraVista/UVISTA_J_21_01_16_psfcat.psf"
dat = fits.open(psf_file)
print(dat.info()) 
print(dat[1].header)
data= dat[1].data

#if args.get_psfex:
#  if not (args.use_existing and os.path.exists(psfex_File)):
#   psfex_file = wget(url_base, base_path + '/psf/', wdir, root + '_psfexcat.psf', logger)
#    logger.info('psfex_file = %s',psfex_file)
#    row['psfex_file'] = psfex_file
#    keep_files.append(psfex_file)


# In[ ]:


#neither my starlist nor Mike's has an obs_flux in the starlist?

def measure_psfex_shapes(df, psfex_file, image_file, noweight, wcs, fwhm): #, logger):
    """Measure shapes of the PSFEx solution at each location.
    """
    #logger.info('Read in PSFEx file: %s',psfex_file)

    #ignore fact that I have no star_file for now
    ind = df.index[df] 
    #ind = df.index[df['star_flag'] == 1]
    #logger.info('ind = %s',ind)
    #n_psf = len(ind)
    #logger.info('n_psf = %s',n_psf)

    df['psfex_dx'] = [ -999. ] * len(df)
    df['psfex_dy'] = [ -999. ] * len(df)
    df['psfex_e1'] = [ -999. ] * len(df)
    df['psfex_e2'] = [ -999. ] * len(df)
    df['psfex_T'] = [ -999. ] * len(df)
    df['psfex_flux'] = [ -999. ] * len(df)
    df['psfex_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'psfex_flag'] = 0

    full_image = galsim.fits.read(image_file, hdu=0)

    if wcs is not None:
        full_image.wcs = wcs

    if not noweight:
        print("I'm using a weight)")
        full_weight = galsim.fits.read(image_file, hdu=0)
        full_weight.array[full_weight.array < 0] = 0.

    stamp_size = 48

    for i in ind:
        x = df['X_IMAGE'].iloc[i]
        y = df['Y_IMAGE'].iloc[i]
        
        #print('Measure PSFEx model shape at ',x,y)
        image_pos = galsim.PositionD(x,y)
        psf_i = psf.getPSF(image_pos)

        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2,
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        im = psf_i.drawImage(image=im, method='no_pixel')
        im *= df['obs_flux'].iloc[i]    

        if noweight:
            wt = None
        else:
            wt = full_weight[b]
            var = wt.copy()
            var.invertSelf()
            im.addNoise(galsim.VariableGaussianNoise(rng, var))
        
        dx, dy, e1, e2, T, flux, flag = ngmix_fit(im, wt, fwhm, x, y, logger)
        
        if np.any(np.isnan([dx,dy,e1,e2,T,flux])):
            logger.info(' *** NaN detected (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
            flag |= BAD_MEASUREMENT
        else:
            df.loc[i, 'psfex_dx'] = dx
            df.loc[i, 'psfex_dy'] = dy
            df.loc[i, 'psfex_e1'] = e1
            df.loc[i, 'psfex_e2'] = e2
            df.loc[i, 'psfex_T'] = T
            df.loc[i, 'psfex_flux'] = flux
        df.loc[i, 'psfex_flag'] |= flag
    logger.info('final psfex_flag = %s',df['psfex_flag'][ind].values)
    #print('df[ind] = ',df.loc[ind].describe())
    flag_outliers(df, ind, 'psfex', 4., logger)


# In[ ]:


measure_psfex_shapes(df, psfex_file, image_file, noweight=False, wcs, fwhm) #, logger)   


# In[30]:


def measure_star_shapes(df, image_file, noweight, wcs, fwhm): #, logger):
    """Measure shapes of the raw stellar images at each location.
    """
    #logger.info('Read in stars in file: %s',image_file)

    ind = df.index[df] #['star_flag'] == 1]
    #logger.info('ind = %s',ind)
    n_psf = len(ind)
    #logger.info('n_psf = %s',n_psf) #ignore logger for now
    print('n_psf = %s',n_psf)

    df['obs_dx'] = [ -999. ] * len(df)
    df['obs_dy'] = [ -999. ] * len(df)
    df['obs_e1'] = [ -999. ] * len(df)
    df['obs_e2'] = [ -999. ] * len(df)
    df['obs_T'] = [ -999. ] * len(df)
    df['obs_flux'] = [ -999. ] * len(df)
    df['obs_flag'] = [ NOT_STAR ] * len(df)
    df.loc[ind, 'obs_flag'] = 0

    full_image = galsim.fits.read(image_file, hdu=0)

    if wcs is not None:
        full_image.wcs = wcs

    if not noweight:
        full_weight = galsim.fits.read(image_file, hdu=2)
        full_weight.array[full_weight.array < 0] = 0.

    stamp_size = 48

    for i in ind:
        x = df['x'].iloc[i]
        y = df['y'].iloc[i]

        #print('Measure shape for star at ',x,y)
        b = galsim.BoundsI(int(x)-stamp_size/2, int(x)+stamp_size/2,
                           int(y)-stamp_size/2, int(y)+stamp_size/2)
        b = b & full_image.bounds
        im = full_image[b]

        if noweight:
            wt = None
        else:
            wt = full_weight[b]

        
        dx, dy, e1, e2, T, flux, flag = ngmix_fit(im, wt, fwhm, x, y, logger)

        #logger.info('ngmix measurement: (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
        if np.any(np.isnan([dx,dy,e1,e2,T,flux])):
            logger.info(' *** NaN detected (%f,%f,%f,%f,%f,%f).',dx,dy,e1,e2,T,flux)
            flag |= BAD_MEASUREMENT
        else:
            df.loc[i, 'obs_dx'] = dx
            df.loc[i, 'obs_dy'] = dy
            df.loc[i, 'obs_e1'] = e1
            df.loc[i, 'obs_e2'] = e2
            df.loc[i, 'obs_T'] = T
            df.loc[i, 'obs_flux'] = flux
        df.loc[i, 'obs_flag'] |= flag
    logger.info('final obs_flag = %s',df['obs_flag'][ind].values)
    #print('df[ind] = ',df.loc[ind].describe())
    flag_outliers(df, ind, 'obs', 4., logger)

    # Any stars that weren't measurable here, don't use for PSF fitting.
    df.loc[df['obs_flag']!=0, 'use'] = False
    
measure_star_shapes(df, image_file, noweight, wcs, fwhm, logger)


# In[ ]:


exp_cat_file = os.path.join(wdir, 'exp_psf_cat_%d.fits'%exp)
        with fitsio.FITS(exp_cat_file,'rw',clobber=True) as f:
            f.write_table(exp_stars_df.to_records(index=False), extname='stars')
            f.write_table(exp_info_df.to_records(index=False), extname='info')


# In[ ]:





# In[ ]:





# In[ ]:




