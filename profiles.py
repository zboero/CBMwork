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

host = 'IATE'

if host == 'local':
    #root_home = '/home/ezequiel/'                                  # Local
    root_home = '/Users/ezequielboero/'                            # Local
    fileProj  = root_home+'Projects/CMB/'                          # Folder of work for the article
    data      = fileProj+'data/'
    graficos  = fileProj+'graficos/'                               # Folder with the plots

elif host == 'IATE':
    root_home = '/home/zboero/'                                    # Clemente y IATE
    fileProj  = root_home+'Projects/CMB/'                          # Folder of work for the article
    data      = fileProj+'data/'                                   # Folder with the data
    graficos  = fileProj+'graficos/'                               # Folder with the plots
    #
    Plnk      = root_home+'Documentos/Planck+Data+Product+Release_2019/'
    fileMask  = Plnk+'Masks/'                                      # Aquí están los archivos del release de Planck 2018


common_mask = fileMask+'COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
####################################################################################
####################################################################################

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# Define the number of galaxies and the size of the map
n_gxs = 1000
nside = 512

# Generate random positions for the galaxies
lon_array   = np.random.uniform( low= 0.0 , high= 2*np.pi, size= n_gxs )
costheta    = np.random.uniform( low= -1.0, high= 1.0    , size= n_gxs )
colat_array = np.arccos( costheta )
lat_array = colat_array - 0.5*np.pi

th_gxs, phi_gxs = colat_array, lon_array

# Create an empty map
n_pix = hp.nside2npix(nside)
map = np.zeros(n_pix)
th_map, phi_map = hp.pix2ang( nside, np.arange(n_pix) )

# Define the radial profile function
def radial_profile( th, phi, th_gxs, phi_gxs ):
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
    Returns
    -------
    A numpy array with the value of the temperature profile coorresponding to each pixel of the map.
    In other words the output is a heaply map as a numpy array type.
    '''
    import numpy as np

    # cos_gamma note that the formula below corresponds to -1 <= cos(th) <= 1 convention
    # which implies that 0 <= th <= 180 in deg or 0 <= th <= pi in rad...This is th must be a co-latitude
    cos_sep = np.cos(th_gxs) * np.cos(th) + np.sin(th_gxs) * np.sin(th) * np.cos(phi - phi_gxs)
    omega = np.arccos( cos_sep )
    
    a1 = np.deg2rad(0.3)
    a2 = np.deg2rad(1.0)
    a3 = np.deg2rad(10.0)
    
    Amp_1 = -8.66
    Amp_2 = 7.0 * np.log10( np.rad2deg(omega) ) - 5.0
    Amp_3 = 5.0 * np.log10( np.rad2deg(omega) ) - 5.0
    Amp_4 = np.zeros( len(cos_sep) )
    
    prof = np.where( omega  > a3, Amp_4, Amp_3 )
    prof = np.where( omega <= a2, Amp_2, prof )
    prof = np.where( omega <= a1, Amp_1, prof )
    
    return prof

# Add the radial profiles to the map
for i in range(n_gxs):
    profile_i = radial_profile( th_map, phi_map , th_gxs[i], phi_gxs[i] )
    map += profile_i

# Normalize the map
map -= np.mean(map)


# Visualize the map
Tmin = np.min(map)
Tmax = np.max(map)
title = 'profile'
output_file = 'profile_.png'
def sky_plot(map, title, Tmin, Tmax, output_file):
    ''' Produces a Mollweide map
        Requires healpy
    Parameters
    ----------
    map : numpy array
        The temperatue values (in units of \muK).
    Tmin : float
        The lower bound for T (values lesser than Tmin will be painted with the
    minimun value of the colorscale)
    Tmax : float
        The upper bound for T (values grater than Tmax will be painted with the
    maximun value of the colorscale)
    output_file : str
        The output file (.png in general)
    Returns
    -------
    The figure with the map masked with the common mask used by Planck.
    '''
    import numpy as np
    import healpy as hp
    import matplotlib.pyplot as plt


    colormap = 'inferno'          #'viridis', 'plasma', 'magma', 'cividis'
    #colormap = 'viridis'
    mask = hp.read_map( common_mask, nest=False )
    mask = hp.ud_grade( mask, nside)
    map = np.where(mask, map, hp.UNSEEN)
    plt.figure()
    hp.mollview(map, coord='G', unit='Temperature [$\mu$K]', xsize=800,
                     title=title, cbar=True, cmap=colormap,
                     min=Tmin, max=Tmax, badcolor='black')
    plt.savefig(output_file)
    plt.close()


sky_plot(map, title, Tmin, Tmax, output_file)






########################################################################################################

# Appendix:
# --------
#

# Radial profile used to fit large and isolated spiral in reference 
# Luparello et. al, MNRAS, Volume 518, Issue 4, February 2023, Pages 5643–5652.
#
#   ------------------------------------------------------------------------------------------------
#     cos_sep = np.cos(th_gxs) * np.cos(th) + np.sin(th_gxs) * np.sin(th) * np.cos(phi - phi_gxs)
#     omega = np.arccos( cos_sep )

#     a1 = np.deg2rad(0.3)
#     a2 = np.deg2rad(1.0)
#     a3 = np.deg2rad(10.0)
    
#     Amp_1 = -8.66
#     Amp_2 = 7.0 * np.log10( np.rad2deg(omega) ) - 5.0
#     Amp_3 = 5.0 * np.log10( np.rad2deg(omega) ) - 5.0
#     Amp_4 = np.zeros( len(cos_sep) )
    
#     prof = np.where( omega  > a3, Amp_4, Amp_3 )
#     prof = np.where( omega <= a2, Amp_2, prof )
#     prof = np.where( omega <= a1, Amp_1, prof )
#   ------------------------------------------------------------------------------------------------
