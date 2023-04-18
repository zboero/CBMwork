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


common_mask = fileMask+'COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
####################################################################################
####################################################################################

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# Define the number of galaxies and the size of the map
nside = 64#512

# n_gxs = 1000
# # Generate random positions for the galaxies
# def rdm_sample(n_gxs):
#     ''' Produces a uniform random sample of points in the unit sphere
        
#     Parameters
#     ----------
#     n_events : int
#         The total number of events of the map
        
#     Returns
#     -------
#     A pandas dataframe with the events
#     '''

#     lon_array   = np.random.uniform( low= 0.0 , high= 2*np.pi, size= n_gxs )
#     costheta    = np.random.uniform( low= -1.0, high= 1.0    , size= n_gxs )
#     colat_array = np.arccos( costheta )         # the range is [0, pi]
#     lat_array   = 0.5*np.pi - colat_array

#     cols_df_rdm = { 'l (rad)': lon_array,
#                     'colat (rad)': colat_array,
#                     'b (rad)': lat_array
#                   }
#     df_rdm = pd.DataFrame( data= cols_df_rdm )

#     return df_rdm

# df_rdm = rdm_sample(n_gxs)
# th_gxs, phi_gxs = df_rdm['colat (rad)'], df_rdm['l (rad)']


# Define the radial profile function
def radial_profile( th, phi, th_gxs, phi_gxs, flag, size, z, dist_5gxy ):
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
    flag : string
        Type of profile: 'original_paper_I' or 'Frode_model' for the moment...
    size : numpy array
        galaxy optical size in kpc
    z : numpy array
        redshift of the galaxy
    dist_5gxy : numpy array
        Distance to the 5th galaxy (This usually comes from a bigger catalogue than the sample)

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


    if (flag == 'original_paper_I'):
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
        
    elif (flag == 'Frode_model'):
        a_deep   = - 30.0
        size_ref = 8.5
        dist_ref = 3.0#np.deg2rad(3.0)
        a_scale  = ( size / size_ref )**2 * ( dist_5gxy / dist_ref )**4

        z_ref    = 0.01
        ext_prof = np.deg2rad(2.0)
        b_dist   = ext_prof *  ( z_ref / z ) * ( dist_ref / dist_5gxy )**(2.8)
        
        prof     = np.zeros( len(cos_sep) )
        lin      = ( a_deep * a_scale ) * ( 1.0 - (1.0 / b_dist) * omega )
        prof     = np.where( omega  > b_dist, prof, lin )
        
    return prof


# Create an empty map
n_pix = hp.nside2npix(nside)
map = np.zeros(n_pix)
th_map, phi_map = hp.pix2ang( nside, np.arange(n_pix) )

# Add the radial profiles to the map
#flag      = 'original_paper_I'
flag      = 'Frode_model'
H0        = 100.0
c_luz     = 299792.458                         # in km/s
#asec2rad  = np.pi/( 180.0*3600.0 )          # Convertion factor between [arcseconds] and [radians].
size      = df_gxs["rad"].to_numpy()   #(df_gxs['v']/H0 * 10**( df_gxs['r_ext'] ) * 1000 * asec2rad).to_numpy()
th_gxs    = np.deg2rad(90.0 - df_gxs["b"].to_numpy()) #np.deg2rad(90.0 - df_gxs["DECdeg"].to_numpy())
phi_gxs   = np.deg2rad(df_gxs["l"].to_numpy())#np.deg2rad(df_gxs["RAdeg"].to_numpy())
z         = df_gxs["z"].to_numpy() #df_gxs["v"].to_numpy() / c_luz
dist_5gxy = df_gxs["d5"].to_numpy() #dist5gxs
n_gxs     = len(df_gxs)

for i in range(n_gxs): #n_gxs
    profile_i = radial_profile( th_map, phi_map , th_gxs[i], phi_gxs[i], flag, size[i], z[i], dist_5gxy[i] )
    map += profile_i

# Normalize the map
map -= np.mean(map)


# Visualize the map
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

    colormap = 'jet'                         # 'inferno', 'viridis', 'plasma', 'magma', 'cividis', etc...
    mask = hp.read_map( common_mask, nest=False )
    mask = hp.ud_grade( mask, nside)
    map = np.where(mask, map, hp.UNSEEN)
    plt.figure()
    hp.mollview(map, coord='G', unit='Temperature [$\mu$K]', xsize=800,
                     title=title, cbar=True, cmap=colormap,
                     min=Tmin, max=Tmax, badcolor='black')
    #hp.visufunc.graticule( dpar=15, dmer=30, coord='G')
    plt.savefig(output_file)
    plt.close()

Tmin = np.min(map)
Tmax = np.max(map)
title = 'profile'
output_file = 'profile_FrodeModel.png'
sky_plot(map, title, Tmin, Tmax, output_file)



# Save the map as a file .fits
filename = 'profile_FrodeModel.fits'
hp.write_map( filename, map, nest=False, coord='G')



def sky_plot( th, phi, th_str, phi_str, coord_sys, title, output_file ):
    ''' Produces a scatter map in Mollweide projection with the position
    of events in the sky.
        
    Parameters
    ----------
    th, phi : numpy arrays
        The angular coordinates on the sphere of the events: co-latitude and longitude
        Units must be in radians.
    th_str, phi_str : strings
        The labels of the coordinates (for example, 'RA', 'Dec' or 'lon', 'lat')
    coord_sys : str
        The coordinate system under use: 'Galactic', 'Equatorial', etc
    title : str
        The title of the plot
    output_file : str
        The output file (.png in general)
        
    Returns
    -------
    The figure with the scatter plot in Mollweide projection.
    '''

    org =0                                                 # Origin of the map
    projection ='mollweide'                                # 2D projection
    x = np.remainder( phi + 2.0*np.pi - org, 2.0*np.pi)    # shift phi (RA) values
    ind = x > np.pi
    x[ind] -= 2.0*np.pi                                    # scale conversion to [-180, 180]
    x = -x                                                 # reverse the scale: East to the left
    y = 0.5*np.pi - th                                     # From colatitude to latitude
    
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder( tick_labels + 360 + org, 360 )
    
    fig = plt.figure( figsize=(10, 5) )
    ax = fig.add_subplot( 111, projection=projection )
    ax.scatter( x, y, s=1.5 )
    ax.set_xticklabels( tick_labels )                      # we add the scale on the x axis
    ax.set_title( title )
    ax.title.set_fontsize(15)
    ax.set_xlabel( phi_str )
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel( th_str )
    ax.yaxis.label.set_fontsize(12)
    ax.grid('True', linestyle='--')
    ax.text( 5.5, -1.3 , coord_sys , fontsize=12)
    plt.savefig( output_file )
    #plt.savefig(graficos+'vMF_dist.png')


# sky map with the events
th_str, phi_str = 'b (deg)', 'l (deg)'
th_gxs  = np.deg2rad(90.0 - df_gxs["b"].to_numpy())#np.deg2rad( 90.0 - th_gxs )
phi_gxs = np.deg2rad(df_gxs["l"].to_numpy())#np.deg2rad( phi_gxs )
coord_sys = 'Galactic'
title   = 'Only spirals'
output_file = graficos+'skymap_FrodeModel.png'
sky_plot( th_gxs, phi_gxs, th_str, phi_str, coord_sys, title, output_file )

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


# Radial profile used to fit enlarged samples of galaxies
# Hansen et. al
#
#   ------------------------------------------------------------------------------------------------
#     a_deep   = - 30.0
#     size_ref = 8.5
#     dist_ref = 3.0#np.deg2rad(3.0)
#     a_scale  = ( size / size_ref )**2 * ( dist_5gxy / dist_ref )**4
#
#     z_ref    = 0.01
#     ext_prof = np.deg2rad(2.0)
#     b_dist   = ext_prof *  ( z_ref / z ) * ( dist_ref / dist_5gxy )**(2.8)
#
#     prof     = np.zeros( len(cos_sep) )
#     lin      = ( a_deep * a_scale ) * ( 1.0 - (1.0 / b_dist) * omega )
#     prof     = np.where( omega  > b_dist, prof, lin )
#   ------------------------------------------------------------------------------------------------

