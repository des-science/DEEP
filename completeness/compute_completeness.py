import numpy as np
import healsparse as hs
import desmasks
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord
import sys


def is_masked(mask_map, ra, dec):
        """
        check if the input positions are masked
        """

        mask_values = mask_map._mask_map.get_values_pos(ra, dec, lonlat=True)
        bounds_values = mask_map._bounds_map.get_values_pos(ra, dec, lonlat=True)

        return (
            (mask_values > 0) | (bounds_values == 0)
        )


def match_cats(ra1, dec1, ra2, dec2):
    # ra1, dec1 are the deep cat.
    c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
    catalogue = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
    idx, d2d, _ = c.match_to_catalog_sky(catalogue)
    return np.where(d2d.value<1./3600.)
    

fields = {'X3':
         {'ban_list':[1,2,3,4,5,9,25,30,31,37,38,41,43,44,48,50,61],
          'coadd_pattern':"SN-X3_C00_r3688p01_i_cat.fits",
          'deep_pattern':"SN-X3_C00_r3901p01_i_cat.fits"},
         'C3':
         {'ban_list':[2,3,8,15,18,21,22,25,28,31,32,45,46,61],
          'coadd_pattern':"SN-C3_C00_r3688p01_i_cat.fits",
          'deep_pattern':"SN-C3_C00_r3901p01_i_cat.fits"},
         'E2':
         {'ban_list':[2,3,5,6,9,10,16,25,31,35,38,47,53,54,58,59,61,62],
          'coadd_pattern':"SN-E2_C00_r3688p01_i_cat.fits",
          'deep_pattern':"SN-E2_C00_r3899p01_i_cat.fits"}}


field = sys.argv[1]
mag_bins = np.arange(20.,28.6,0.2)

deep_counts = np.zeros(len(mag_bins)-1)
coadd_counts = np.zeros(len(mag_bins)-1)

for chip in range(17):
    if chip+1 not in fields[field]['ban_list']:
        hs_bound = "SN-{}_C{:02d}-griz-bounds-healsparse.fits".format(field,chip+1)
        hs_mask = "SN-{}_C{:02d}-griz-healsparse.fits".format(field,chip+1)
        tile_mask = desmasks.TileMask(mask_fname=hs_mask, bounds_fname=hs_bound)

        fdeep = fields[field]['deep_pattern'].split('00')
        fcoadd = fields[field]['coadd_pattern'].split('00')

        deep = pyfits.open(fdeep[0]+"{:02d}".format(chip+1)+fdeep[1])[1].data
        coadd = pyfits.open(fcoadd[0]+"{:02d}".format(chip+1)+fcoadd[1])[1].data

        mdeep = is_masked(tile_mask,deep['ALPHAWIN_J2000'], deep['DELTAWIN_J2000'])
        mcoadd = is_masked(tile_mask,coadd['ALPHAWIN_J2000'], coadd['DELTAWIN_J2000'])

        hdeep = np.histogram(deep['MAG_DETMODEL'][mdeep==False], bins=mag_bins)[0]
        deep_counts += hdeep

        # now find which ones have a match in the coadd set
        matched = match_cats(deep['ALPHAWIN_J2000'][mdeep==False], deep['DELTAWIN_J2000'][mdeep==False],
                        coadd['ALPHAWIN_J2000'][mcoadd==False], coadd['DELTAWIN_J2000'][mcoadd==False])
        hcoadd = np.histogram(deep['MAG_DETMODEL'][mdeep==False][matched], bins=mag_bins)[0]
        coadd_counts += hcoadd


compl = coadd_counts/deep_counts
[print((mag_bins[i]+mag_bins[i+1])/2., compl[i]) for i in range(len(compl))]
        
