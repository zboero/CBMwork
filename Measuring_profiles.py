#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:40:33 2023

@author: ezequiel
"""

# ------------------------------------------------------------------------------
# Program:
# -------
#
# This program produces a Healpix maps with the superposition of a radial profile
# associated to each galaxy of a given sample.
# The profile is an angular profile (projected on the sky) and it is assumed
# symmetric.
#
# ------------------------------------------------------------------------------

host = 'local'

if host == 'local':
    #root_home = '/home/ezequiel/'                                  # Local
    root_home = '/Users/ezequielboero/'                            # Local
    fileProj  = root_home+'Projects/CMB/'                          # Folder of work for the article
    data      = fileProj+'data/'
    graficos  = fileProj+'graficos/'                               # Folder with the plots
    Plnk      = root_home+'Planck_dataproducts/PR3/'
    fileMask  = Plnk+'Masks/'

elif host == 'IATE':
    root_home = '/home/zboero/'                                    # Clemente y IATE
    fileProj  = root_home+'Projects/CMB/'                          # Folder of work for the article
    data      = fileProj+'data/'                                   # Folder with the data
    graficos  = fileProj+'graficos/'                               # Folder with the plots
    #
    Plnk      = root_home+'Documentos/Planck+Data+Product+Release_2019/'
    fileMask  = Plnk+'Masks/'                                      # Aquí están los archivos del release de Planck 2018

SMICA       = root_home+'Planck_dataproducts/PR3/CMB_Maps/COM_CMB_IQU-smica_2048_R3.00_full.fits')
common_mask = fileMask+'COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
####################################################################################
####################################################################################

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# Read map
map, header = hp.read_map( map_name, nest=False, h=True, field=0 )    # field=0 --> Tempearture array

# Set the size of the map
nside = 128
#map   = hp.ud_grade( map, nside_out=nside, order_in='RING', order_out='RING')

# Angles of each pixel of the map
n_pix = hp.nside2npix(nside)
th_map, phi_map = hp.pix2ang( nside, np.arange(n_pix) )

# Define the radial profile function
def T_dist_disc( th, phi, th_gxs, phi_gxs, th_max, n_gxs):
    ''' Radial profile around the position of galaxies
    Parameters
    ----------
    th : numpy array
        co-latitude (0,pi) in radians for each pixel of a healpix map.
    phi : numpy array
        longitude (0,2*pi) in radians for each pixel of a healpix map.
    th_gxs : numpy array
        co-latitude (0,pi) in radians for each galaxy in the sample.
    phi_gxs : numpy array
        longitude (0,2*pi) in radians for each galaxy in the sample.
    th_max : float
        maximum angular distance for the radial profie

    Returns
    -------
    A numpy array with the 2d Histogram radius vs Temperature of the pixels in the disc
    centered around a galaxy
    '''
    import numpy as np

    # cos_gamma note that the formula below corresponds to -1 <= cos(th) <= 1 convention
    # which implies that 0 <= th <= 180 in deg or 0 <= th <= pi in rad...This is th must be a co-latitude
    cos_sep = np.cos(th_gxs) * np.cos(th) + np.sin(th_gxs) * np.sin(th) * np.cos(phi - phi_gxs)
    omega   = np.arccos( cos_sep )
    disc    = np.where( omega > th_max, 0.0, omega )
    id_disc = np.nonzero(disc)[0]
    bins_r  = 10
    bins_T  = 1000
    r_min   = 0.0
    r_max   = th_max
    T_min   = -500
    T_max   = 500
    dr      = np.rad2deg(th_max) / bins_r
    r, T    = np.rad2deg(omega[id_disc]), 0.5 * map[id_disc]/(np.rad2deg(omega[id_disc]) * np.pi * dr )
    Hist2d  = np.histogram2d( r, T/n_gxs,  bins=[bins_r, bins_T], range=[[r_min, r_max], [T_min, T_max]], density=False )[0]
            
    return Hist2d

# Coordinates of the galaxies in the sample
th_gxs    = np.deg2rad(90.0 - df_gxs["b"].to_numpy()) #np.deg2rad(90.0 - df_gxs["DECdeg"].to_numpy())
phi_gxs   = np.deg2rad(df_gxs["l"].to_numpy())#np.deg2rad(df_gxs["RAdeg"].to_numpy())

n_gxs     = len(df_gxs)
th_max    = np.deg2rad(10)
prof      = np.zeros( (10,1000) )
for i in range(n_gxs):
    prof_i = T_dist_disc( th_map, phi_map , th_gxs[i], phi_gxs[i], th_max, n_gxs )
    prof += prof_i


# Plot of the Histogram
r_edges = np.linspace(0, th_max, 10)
T_edges = np.linspace(-500, 500, 1000)
prof = prof.T

fig = plt.figure( figsize=(10, 5) )
plt.pcolormesh(r_edges, T_edges, prof, cmap='rainbow')
plt.colorbar()
plt.title( ' Distribution of temperature around galaxies ')
plt.xlabel( 'Angle [rad]', fontsize='x-large')
plt.ylabel( 'Temperature [muK]', fontsize='x-large')
plt.axes
#plt.legend(loc='upper right', fontsize='large', markerscale=3.0)
#plt.savefig(output_file)
plt.tight_layout()
plt.grid()
plt.show()


########################################################################################################

