#######################################
#
# snapdragons2fits.py
#
# This snippet converts snapdragons 
#  output to a fits file readable with 
#  standard astro tools
#
# Input: needs to be ascii (a bit slow)
# Output: FITS
#
#######################################

# input file name
inputfile='../output/GeneratedStars.dat'
# output file
outfile='../output/snapdragons_gaia.fits'


from astropy.io import fits as pyfits
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import constants as const
import struct
from galpy.util import bovy_coords

# reading the data
rdata=np.loadtxt(inputfile)

# radians -> degrees
rdata[:,0]=180.0*rdata[:,0]/np.pi
rdata[:,1]=180.0*rdata[:,1]/np.pi
rdata[:,6]=180.0*rdata[:,6]/np.pi
rdata[:,7]=180.0*rdata[:,7]/np.pi
rdata[:,12]=180.0*rdata[:,12]/np.pi
rdata[:,13]=180.0*rdata[:,13]/np.pi
# arcsec -> mas
rdata[:,2]=rdata[:,2]*1000.0
rdata[:,8]=rdata[:,8]*1000.0
rdata[:,14]=rdata[:,14]*1000.0
# arcsec/yr -> mas/yr
rdata[:,3]=rdata[:,3]*1000.0
rdata[:,4]=rdata[:,4]*1000.0
rdata[:,9]=rdata[:,9]*1000.0
rdata[:,10]=rdata[:,10]*1000.0
rdata[:,15]=rdata[:,15]*1000.0
rdata[:,16]=rdata[:,16]*1000.0

# Galactic coordinates
# True l and b
RA_true=rdata[:,0]
DEC_true=rdata[:,1]
Tllbb=bovy_coords.radec_to_lb(RA_true,DEC_true,degree=True,epoch=2000.0)
GLON_true=Tllbb[:,0]
GLAT_true=Tllbb[:,1]
# True pmGLON, pmGLAT
pmRA_true=rdata[:,3]
pmDEC_true=rdata[:,4]
Tpmllbb=bovy_coords.pmrapmdec_to_pmllpmbb(pmRA_true,pmDEC_true \
   ,RA_true,DEC_true,degree=True,epoch=2000.0)
pmGLON_true=Tpmllbb[:,0]
pmGLAT_true=Tpmllbb[:,1]
# observed 
RA_obs=rdata[:,6]
DEC_obs=rdata[:,7]
pmRA_obs=rdata[:,9]
pmDEC_obs=rdata[:,10]
Tpmllbb=bovy_coords.pmrapmdec_to_pmllpmbb(pmRA_obs,pmDEC_obs \
   ,RA_true,DEC_true,degree=True,epoch=2000.0)
pmGLON_obs=Tpmllbb[:,0]
pmGLAT_obs=Tpmllbb[:,1]
# error
e_pmRA=rdata[:,12]
e_pmDEC=rdata[:,13]
Tpmllbb=bovy_coords.pmrapmdec_to_pmllpmbb(e_pmRA,e_pmDEC \
   ,RA_true,DEC_true,degree=True,epoch=2000.0)
e_pmGLON=Tpmllbb[:,0]
e_pmGLAT=Tpmllbb[:,1]

tbhdu = pyfits.BinTableHDU.from_columns([\
  pyfits.Column(name='RA_true',unit='(degree)',format='D',array=rdata[:,0]),\
  pyfits.Column(name='DEC_true',unit='(degree)',format='D',array=rdata[:,1]),\
  pyfits.Column(name='Plx_true',unit='(mas)',format='D',array=rdata[:,2]),\
  pyfits.Column(name='pmRA_true',unit='(mas/yr)',format='D',array=rdata[:,3]),\
  pyfits.Column(name='pmDEC_true',unit='(mas/yr)',format='D',array=rdata[:,4]),\
  pyfits.Column(name='HRV_true',unit='(km/s)',format='D',array=rdata[:,5]),\
# observed
  pyfits.Column(name='RA_obs',unit='(degree)',format='D',array=rdata[:,6]),\
  pyfits.Column(name='DEC_obs',unit='(degree)',format='D',array=rdata[:,7]),\
  pyfits.Column(name='Plx_obs',unit='(mas)',format='D',array=rdata[:,8]),\
  pyfits.Column(name='pmRA_obs',unit='(mas/yr)',format='D',array=rdata[:,9]),\
  pyfits.Column(name='pmDEC_obs',unit='(mas/yr)',format='D',array=rdata[:,10]),\
  pyfits.Column(name='HRV_obs',unit='(km/s)',format='D',array=rdata[:,11]),\
# error
  pyfits.Column(name='e_RA',unit='(degree)',format='D',array=rdata[:,12]),\
  pyfits.Column(name='e_DEC',unit='(degree)',format='D',array=rdata[:,13]),\
  pyfits.Column(name='e_Plx',unit='(mas)',format='D',array=rdata[:,14]),\
  pyfits.Column(name='e_pmRA',unit='(mas/yr)',format='D',array=rdata[:,15]),\
  pyfits.Column(name='e_pmDEC',unit='(mas/yr)',format='D',array=rdata[:,16]),\
  pyfits.Column(name='e_HRV',unit='(km/s)',format='D',array=rdata[:,17]),\
# True
  pyfits.Column(name='G_true',unit='(mag)',format='D',array=rdata[:,18]),\
  pyfits.Column(name='BP_RP_true',unit='(mag)',format='D',array=rdata[:,19]),\
# Observed
  pyfits.Column(name='G_obs',unit='(mag)',format='D',array=rdata[:,20]),\
  pyfits.Column(name='BP_RP_obs',unit='(mag)',format='D',array=rdata[:,21]),\
# Error
  pyfits.Column(name='e_G',unit='(mag)',format='D',array=rdata[:,22]),\
  pyfits.Column(name='e_BP_RP',unit='(mag)',format='D',array=rdata[:,23]),\
# True
  pyfits.Column(name='V_true',unit='(mag)',format='D',array=rdata[:,24]),\
  pyfits.Column(name='dist_true',unit='(kpc)',format='D',array=rdata[:,25]),\
  pyfits.Column(name='x_true',unit='(kpc)',format='D',array=rdata[:,26]),\
  pyfits.Column(name='y_true',unit='(kpc)',format='D',array=rdata[:,27]),\
  pyfits.Column(name='z_true',unit='(kpc)',format='D',array=rdata[:,28]),\
  pyfits.Column(name='vx_true',unit='(km/s)',format='D',array=rdata[:,29]),\
  pyfits.Column(name='vy_true',unit='(km/s)',format='D',array=rdata[:,30]),\
  pyfits.Column(name='vz_true',unit='(km/s)',format='D',array=rdata[:,31]),\
  pyfits.Column(name='vr_true',unit='(km/s)',format='D',array=rdata[:,50]),\
  pyfits.Column(name='vphi_true',unit='(km/s)',format='D',array=rdata[:,51]),\
  pyfits.Column(name='mass',unit='(M_sol)',format='D',array=rdata[:,32]),\
  pyfits.Column(name='G_RVS',unit='(mag)',format='D',array=rdata[:,33]),\
  pyfits.Column(name='V-I',unit='(mag)',format='D',array=rdata[:,34]),\
  pyfits.Column(name='Av_true',unit='(mag)',format='D',array=rdata[:,35]),\
  pyfits.Column(name='[Fe/H]_true',unit='(dex)',format='D',array=np.log10((rdata[:,38]/1000000.)/0.0152)),\
  pyfits.Column(name='age_true',unit='(Gyr)',format='D',array=10.**(rdata[:,39])),\
# Observed
  pyfits.Column(name='dist_obs',unit='(kpc)',format='D',array=rdata[:,40]),\
  pyfits.Column(name='x_obs',unit='(kpc)',format='D',array=rdata[:,41]),\
  pyfits.Column(name='y_obs',unit='(kpc)',format='D',array=rdata[:,42]),\
  pyfits.Column(name='z_obs',unit='(kpc)',format='D',array=rdata[:,43]),\
  pyfits.Column(name='vx_obs',unit='(km/s)',format='D',array=rdata[:,44]),\
  pyfits.Column(name='vy_obs',unit='(km/s)',format='D',array=rdata[:,45]),\
  pyfits.Column(name='vz_obs',unit='(km/s)',format='D',array=rdata[:,46]),\
  pyfits.Column(name='vr_obs',unit='(km/s)',format='D',array=rdata[:,47]),\
  pyfits.Column(name='vphi_obs',unit='(km/s)',format='D',array=rdata[:,48]),\
# GLON and GLAT 
  pyfits.Column(name='GLON_true',unit='(degree)',format='D',array=GLON_true),\
  pyfits.Column(name='GLAT_true',unit='(degree)',format='D',array=GLAT_true),\
# pmGLON and pmGLAT true
  pyfits.Column(name='pmGLON_true',unit='(mas/yr)',format='D',\
                array=pmGLON_true),\
  pyfits.Column(name='pmGLAT_true',unit='(mas/yr)',format='D',\
                array=pmGLAT_true),\
# observed
  pyfits.Column(name='pmGLON_obs',unit='(mas/yr)',format='D',\
                array=pmGLON_obs),\
  pyfits.Column(name='pmGLAT_obs',unit='(mas/yr)',format='D',\
                array=pmGLAT_obs),\
# Errors
  pyfits.Column(name='e_pmGLON',unit='(mas/yr)',format='D',\
    array=e_pmGLON),\
  pyfits.Column(name='e_pmGLAT',unit='(mas/yr)',format='D', array=e_pmGLAT)
])
tbhdu.writeto(outfile,clobber=True)

