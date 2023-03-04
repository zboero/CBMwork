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
# This program computes the power spectrum of a set syntehtic maps that are generated
# with a given spectrum as input and afterward modified with a temperature radial
# profile around the position of galaxies.
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
    PS_folder = Plnk+'Cosmology_by_Planck/Power_Spectrae/'         # Ruta donde están los Cl's medidos por (Planck 2018)


calPlanck = 0.1000442*(10**1)                                      # [y_{\rm cal}] Esto tiene que multiplicar

PS      = PS_folder+'COM_PowerSpect_CMB-TT-full_R3.01.txt'         # String del espectro de Planck
PS_bfit = PS_folder+'COM_PowerSpect_CMB-base-plikHM-TTTEEE-'\
          +'lowl-lowE-lensing-minimum-theory_R3.01.txt'            # String del best fit teorico dado por Planck

####################################################################################
####################################################################################

import numpy as np
import healpy as hp


# Load the .fits files
load_maps       = data+'maps_for_ezequiel/'
foregr          = load_maps+'foregroundmap_'
path_Mod1       = foregr+'model1_nside512_muK_ringordering.fits'
path_Mod1extz   = foregr+'model1_extendedz_nside512_muK_ringordering.fits'
path_Mod1extz01 = foregr+'model1_extendedz0.1_nside512_muK_ringordering.fits'
path_Mod1extz10 = foregr+'model1_extendedz1.0_nside512_muK_ringordering.fits'
path_Mod2       = foregr+'model2_nside512_muK_ringordering.fits'
path_Mod2extz   = foregr+'model2_extendedz_nside512_muK_ringordering.fits'
path_Mod2extz01 = foregr+'model2_extendedz0.1_nside512_muK_ringordering.fits'
path_Mod3       = foregr+'model3_nside512_muK_ringordering.fits'
path_Mod3extz   = foregr+'model3_extendedz_nside512_muK_ringordering.fits'
path_Mod4       = foregr+'model4_nside512_muK_ringordering.fits'
path_Mod4extz   = foregr+'model4_extendedz_nside512_muK_ringordering.fits'
path_Mod4extz01 = foregr+'model4_extendedz0.1_nside512_muK_ringordering.fits'


# Read the power spectrum released by Planck2018
Dls_data2018_table = np.genfromtxt(PS     , skip_header=1, unpack=True)        # measured
Dls_bfit2018_table = np.genfromtxt(PS_bfit, skip_header=1, unpack=True)        # best-fit based model LCMD

# Get the TT spectrum (Dls)
Dls_data2018 = Dls_data2018_table[1]                                           # in units of muK**2
Dls_bfit2018 = Dls_bfit2018_table[1]                                           # in units of muK**2

# Get the multipole index
ell_idx = Dls_data2018_table[0]                                                # ell's from 2 to l_max
ell     = np.append( np.array([0.0, 1.0]), ell_idx )                           # ell's from 0 to l_max

# Get the TT spectrum (Cls)
mon_dip = np.array([0.0, 0.0])                                                 # Vanishing monopole and dipole
clfactor = 2.0*np.pi/( ell_idx * (ell_idx + 1) )
Cls_data  = np.append( mon_dip, Dls_data2018 * clfactor )
Cls_bfit  = np.append( mon_dip, Dls_bfit2018 * clfactor )

# ---------------------------------------------------------------------------------
# 10) Bloque de procesamientos
# ---------------------------------------------------------------------------------

# Choose the model
path_Mod = path_Mod1extz10
Mod, header = hp.read_map( path_Mod, nest=False, h=True, field=(0) )

# Set the simulation of maps
n_maps = 10
nside  = 512                                                                   # nisde resolution of the healpix map

lMax = 3*nside - 1                                                             # max multipole of the spectrae

# Generate maps with synfast, modify it and recover the modified spectrum with anafast
def Modified_spectrum(Cl_input, n_maps, nside, hp_prof):
    ''' Given an initial spectrum "Cl_input" it produces a number "n_maps"
        of maps with resolution "nside".
        Each map is then modified with the map "hp_prof" and their spectrae are measured.
        The output consist of the averaged spectrum over the n_maps and their statistical
        uncertainty for each multipole (ell).
    
    Parameters
    ----------
    Cl_input : numpy array
        The spectrum of input.
    n_maps : float
        The numbers of maps to generate and the total number of spectrum to average
    nside : float
        The nside resolution of the healpix map
    hp_prof: numpy array
        The map that modifies the synthetic ones genereated with Cl_input.
        
    Returns: list with two arrays
    -------
        The final averaged Cls, namely "Cl_output" and the error estimated trough the
        standard deviation for each ell, namely "std_Cl".
    '''

    import numpy as np
    import healpy as hp

    it_e = 3                                                           # Number of iter of the algorithm (standard is 3)
    l_max = 3*nside-1                                                  # highest multipole allowed to be recovered
    
    Cls_list     = []
    Cls_Plk_list = []
    # loop over each map
    for i in range(0, n_maps):
        # Generate the map and modify it
        Map_syn = hp.synfast( Cl_input, nside, pol=False, alm=False, pixwin=False, fwhm=0.0 )
        Map_mod = Map_syn + hp_prof
        
        # For the modified map...
        # Recover the spectrum and append the result in a list (we will use it later to get the mean and std)
        Cls_recov = hp.anafast( Map_mod, nspec=None, lmax=l_max, iter=it_e, alm=False)        
        Cls_list.append(Cls_recov)

        # For the synthetic map previously generated...
        # Recover the spectrum and append the result in a list (we will use it later to get the mean and std)
        Cls_recov_Plk = hp.anafast( Map_syn, nspec=None, lmax=l_max, iter=it_e, alm=False)
        Cls_Plk_list.append(Cls_recov_Plk)

        
    Cl_output         = np.mean(Cls_list, axis=0)
    Std_output        = np.std( Cls_list, axis=0)
    Cl_output_Plk     = np.mean(Cls_Plk_list, axis=0)
    Std_output_Plk   = np.std( Cls_Plk_list, axis=0)

    return(Cl_output, Std_output, Cl_output_Plk, Std_output_Plk)


# Compute the modified spectrum (Cls, err_Cls, etc...)
Cls, err_Cls, Cls_Plk, err_Cls_Plk = Modified_spectrum( Cls_bfit, n_maps, nside, Mod )         # We use the bfit spectrum
ell = ell[ 2 : lMax ]

Dls         = Cls[ 2 : lMax ] * (2.*np.pi)**(-1.) * ( ell * ( ell + 1 ) )
err_Dls     = err_Cls[ 2 : lMax ] * (2.*np.pi)**(-1.) * ( ell * ( ell + 1 ) )
Dls_Plk     = Cls_Plk[ 2 : lMax ] * (2.*np.pi)**(-1.) * ( ell * ( ell + 1 ) )
err_Dls_Plk = err_Cls_Plk[ 2 : lMax ] * (2.*np.pi)**(-1.) * ( ell * ( ell + 1 ) )

# Save the data in .dat format
stack = np.column_stack( ( ell, Dls, err_Dls, Dls_Plk, err_Dls_Plk, Dls_bfit2018[2:lMax] ) )
header_dat = "       ell        ,       Dls        ,      err_Dls        ,     Dls_recov      ,   err_Dls_recov     ,   Dls_bfitPlanck "

output_file = data+'foregroundmap_'
output_file = output_file+'model1extendedz10_nside512_muK_ringordering.dat'
np.savetxt( output_file, stack, header=header_dat, \
            fmt=['%18.f','%18.12E','%18.12E','%18.12E','%18.12E','%18.12E']
            )


# Plot the spectrum and its modification
def spectrum_plot(ell, Dls, err_Dls, Dls_Plk, lMax, title, output_png):

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig1  = plt.figure(constrained_layout=True)
    spec1 = gridspec.GridSpec(ncols=1, nrows=7, figure=fig1)

    f1_ax1 = fig1.add_subplot(spec1[0:4,:])
    plt.errorbar(ell, Dls, yerr=err_Dls, fmt='o', marker='s', markersize=2, capsize=0.8, color='black', label='Modified')
    plt.semilogx( ell, Dls_Plk, '.', markersize=2, color='orange', label='Bfit Planck 2018')
    plt.setp( f1_ax1.get_xticklabels(), visible=False )
    plt.title( title, loc='center')
    plt.axes
    plt.xscale('log')
    plt.xlim(0,lMax+1)
    plt.ylim(1,7000)
    plt.ylabel('$ ( \ell(\ell+1)/2\pi)C_\ell^\mathrm{TT} [\mu K^2]$',fontsize='small')
    plt.legend()

    f1_ax2 = fig1.add_subplot(spec1[4:7,:])
    plt.semilogx( ell, (Dls_Plk - Dls), '.',markersize=2, label='Recov - Modif')
    plt.axhline(y=0.00, color='grey', linewidth=1.0)
    plt.xlim(0,lMax+1)
    plt.xlabel('$\ell$ (Multipole)', fontsize='small')
    plt.ylabel('$ ( \ell(\ell+1)/2\pi) \Delta C_\ell^\mathrm{TT} [\mu K^2]$',fontsize='small')
    plt.legend()

    plt.savefig(output_png)

    plt.close()


# Make the plot
title = 'Power Spectrum Model 1 extendedz1.0'
output_png = 'spectrum_.png'
spectrum_plot(ell, Dls, err_Dls, Dls_Plk, lMax, title, output_png)
