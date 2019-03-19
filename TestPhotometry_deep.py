#!/usr/bin/env python
# coding: utf-8

# In[1]:


#! /usr/bin/env python
# Test photometry
# Simple plots eg. colour-colour
# Match and compare deep data to wide, per galaxy

from __future__ import print_function
import os
import numpy as np
#from read_psf_cats import read_data, band_combinations
import fitsio
import treecorr
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py as h
from astropy.table import Table,join


# In[2]:


#read in Erin's catalogue
cosmosfile='/global/homes/a/aamon/DES/DEStests/DEEP/MOFcats/SN-C3/run-vd03-SN-C3_C01_r3688p01.fits'

data = fits.open(cosmosfile) 
data.info() 
print(data.info)
deep=Table(data[1].data)
print(min(deep['ra']),max(deep['ra']))

cols = data[1].columns
print(cols)

print(deep['bdf_mag'][:,2]) #ugriz YJHKs
#print(len(deep['mag_auto']))

deepra=deep['ra']
deepdec=deep['dec']
deepu=deep['bdf_mag'][:,0]
deepr=deep['bdf_mag'][:,2]
deepi=deep['bdf_mag'][:,3]

deepra=deepra[deepr>0]
deepdec=deepdec[deepr>0]
deepu=deepu[deepr>0]
deepi=deepi[deepr>0]
deepr=deepr[deepr>0]
print(len(deepra))
deepra=deepra[deepi>0]
deepdec=deepdec[deepi>0]
deepu=deepu[deepi>0]
deepr=deepr[deepi>0]
deepi=deepi[deepi>0]
print(len(deepra))

#print(deep['id'])
#print(deep['ra'])

#ASSUMING MAGS ARE UGRI YJHKS


# In[3]:


#colour-colour
plt.scatter(deep['bdf_mag'][2]-deep['bdf_mag'][3],deep['bdf_mag'][1]-deep['bdf_mag'][2], marker='.') #,markersize=10 )
plt.xlabel('r-i')
plt.ylabel('g-r')


# In[4]:


#read in gold 
catname = '/global/cscratch1/sd/troxel/cats_des_y3/Y3_mastercat_v2_6_20_18.h5'

f = h.File(catname,'r')
print(f['catalog'].keys())
print(f['catalog/gold'].keys())
print(f['catalog/gold/coadd_object_id'])

#FLAGS_GOLD  If you are using SExtractor quantities, you should add (FLAGS_GOLD & 1111000) = 0; and (FLAGS_BADREGIONS & 01) = 0

#gflag = np.array(f['catalog/gold/flags_gold'])
#print(gflag)

ra = np.array(f['catalog/gold/ra'])#[star_mask] 
dec = np.array(f['catalog/gold/dec'])#[star_mask]
print(len(ra))

"""#cosmos only
ra=ra[(ra<53)] # & (ra>49)]
#dec=dec[cosmosonly]
print(len(ra))
ra=ra[(ra>51)]
print(len(ra))"""

zeropt=30
r = zeropt- 2.5*np.log10(np.array(f['catalog/gold/sof_cm_flux_corrected_r']))#[star_mask]
i = zeropt- 2.5*np.log10(np.array(f['catalog/gold/sof_cm_flux_corrected_i']))#[star_mask]
z = zeropt- 2.5*np.log10(np.array(f['catalog/gold/sof_cm_flux_corrected_z']))#[star_mask]

print(len(r))

"""cosmosonly=np.where(  (ra<max(deep['ra'])) & (ra>min(deep['ra']))
               &  (dec<max(deep['dec'])) & (dec>min(deep['dec']))   )
i=i[cosmosonly]
z=z[cosmosonly]
r=r[cosmosonly]
print(len(r))"""
print(min(ra)) 
print(max(ra))
ra=ra[np.where((i<30) & (i>0))]
dec=dec[np.where((i<30)& (i>0))]
z=z[np.where((i<30)& (i>0))]
r=r[np.where((i<30)& (i>0))]
i=i[np.where((i<30)& (i>0))]
print(len(ra))

ra[ra > 180] -= 360

gold=np.column_stack((ra,dec,r,i,z))
print(gold)
#gold = gold[gold[:,0].argsort()][:1000000]

gold=gold[np.where(ra>min(deep['ra']))]
#,max(deep['ra'])

print(gold)
print(len(gold))
goldra=gold[:,0]
golddec=gold[:,1]
goldr=gold[:,2]
goldi=gold[:,3]
goldz=gold[:,4]


# In[5]:


#match galaxies by ra and dec

from astropy.coordinates import SkyCoord
from astropy import units as u

goldcat = SkyCoord(ra=goldra*u.degree, dec=golddec*u.degree)  
catalog = SkyCoord(ra=deepra*u.degree, dec=deepdec*u.degree)  
idx, d2d, d3d = catalog.match_to_catalog_sky(goldcat, nthneighbor=1) 

print(goldra[idx])


# In[6]:


print(len(goldra)) 
print(len(deep['bdf_mag'])) 
print(len(d2d))  
print(d2d)
print(d2d.arcsecond)
plt.hist(d2d.arcsecond, 50, range=(0, 20)) #20 is the max matching range in arcmin
plt.xlabel('d2d (arcsec)')
print(deepra[d2d.arcsecond < 10])


# In[34]:


matchlim=1
plt.scatter(deepra[np.where(d2d.arcsecond < matchlim)]-goldra[idx][np.where(d2d.arcsecond < matchlim)],deepra[np.where(d2d.arcsecond < matchlim)], marker='.')
plt.xlabel('DEEP RA-GOLD RA')
plt.ylabel('DEEP RA')
plt.ticklabel_format(useOffset=False)
#plt.xlim(min(goldra[idx][np.where(d2d < matchlim)]),max(goldra[idx][np.where(d2d < matchlim)]) )
plt.ylim(min(deepra[np.where(d2d.arcsecond < matchlim)]),max(deepra[np.where(d2d.arcsecond < matchlim)]) )
plt.xlim(min(deepra[np.where(d2d.arcsecond < matchlim)]-goldra[idx][np.where(d2d.arcsecond < matchlim)]),max(deepra[np.where(d2d.arcsecond < matchlim)]-goldra[idx][np.where(d2d.arcsecond < matchlim)]) )
print(min(goldra[idx]),max(goldra[idx]) )
print(min(deepra),max(deepra))


# In[35]:


#plot magnitudes
#print(len(goldr[idx][np.where(d2d.arcsecond < matchlim)]))
#print(len(deep['bdf_mag'][np.where(d2d.arcsecond < matchlim)]))

#3372/14324  #gold matches/all deep ~quarter
print("percentage matched: ", float(len(deepr[np.where(d2d.arcsecond < matchlim)]))/float(len(deepr))*100.)
#fit = np.polyfit(goldr[idx][np.where(d2d.arcsecond < matchlim)], deepr[np.where(d2d.arcsecond < matchlim)], 1)
#fit_fn = np.poly1d(fit) 
# fit_fn is now a function which takes in x and returns an estimate for y
print(len(deepr))
print(len(deepr[deepr>0]))
#plt.plot(goldr[idx][np.where(d2d.arcsecond < matchlim)], fit_fn(goldr[idx][np.where(d2d.arcsecond < matchlim)]), '--k')
#x = np.linspace(14, 40, 1000)
#plt.plot(x,x,color='red')
plt.scatter(goldr[idx][np.where(d2d.arcsecond < matchlim)], deepr[np.where(d2d.arcsecond < matchlim)]-goldr[idx][np.where(d2d.arcsecond < matchlim)], marker='.', facecolors='lightblue', color='blue',alpha=0.5)
plt.xlim(14,27)
plt.axhline(y=0, color='red')
#plt.ylim(16,38)
plt.xlabel('GOLD r')
plt.ylabel('DEEP r - GOLD r')


# In[28]:


#plot magnitudes

#fit = np.polyfit(goldr[idx][np.where(d2d.arcsecond < matchlim)], deepr[np.where(d2d.arcsecond < matchlim)], 1)
#fit_fn = np.poly1d(fit) 
# fit_fn is now a function which takes in x and returns an estimate for y
#plt.plot(goldr[idx][np.where(d2d.arcsecond < matchlim)], fit_fn(goldr[idx][np.where(d2d.arcsecond < matchlim)]), '--k')
#x = np.linspace(14, 40, 1000)
#plt.plot(x,x,color='red')
plt.scatter(goldi[idx][np.where(d2d.arcsecond < matchlim)], deepi[np.where(d2d.arcsecond < matchlim)]-goldi[idx][np.where(d2d.arcsecond < matchlim)], marker='.', facecolors='lightblue', color='blue',alpha=0.5)
plt.xlim(15,31)
#plt.ylim(15,30)
plt.axhline(y=0, color='red')
plt.xlabel('GOLD i')
plt.ylabel('DEEP i - GOLD i')


# In[29]:


plt.hist(deepi, 50, range=(15, 30))
plt.hist(deepi[np.where(d2d.arcsecond < matchlim)], 50, range=(15, 30))
plt.xlabel('i')


# In[30]:


plt.hist(deepi[np.where(d2d.arcsecond < matchlim)]-goldi[idx][np.where(d2d.arcsecond < matchlim)], 50, range=(-1, 1))
plt.axvline(x=0, color='red')
#plt.hist(deepi[np.where(d2d.arcsecond < matchlim)], 50, range=(15, 30))
plt.xlabel('delta i')


# In[31]:


plt.hist(deepr[np.where(d2d.arcsecond < matchlim)]-goldr[idx][np.where(d2d.arcsecond < matchlim)], 50, range=(-1, 1))
plt.axvline(x=0, color='red')
#plt.hist(deepi[np.where(d2d.arcsecond < matchlim)], 50, range=(15, 30))
plt.xlabel('delta r')


# In[ ]:




