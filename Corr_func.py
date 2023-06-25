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
# This program computes the correlation function from the 
# power spectrum of a set syntehtic maps that are generated 
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

path_testMod1       = foregr+'testmodel_extendedz0.05_nside512_muK_ringordering.fits'
path_testMod2       = foregr+'testmodel2_extendedz0.05_nside512_muK_ringordering.fits'

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


# Define the two-point angular correlation function in terms of the spectrum
def TT_Corr(Dl, gamma, lMax):
    import numpy as np
    from scipy.special import eval_legendre

    sum=0
    for n in range(2,lMax):
        sum= sum + (2*n+1)/(4*np.pi)*(2*np.pi/(n*(n+1))*Dl[n-2])*eval_legendre(n,np.cos(gamma))
    return sum

ang_corr = TT_Corr(Dl, gamma, lMax)
# Plot of the TT_corr function
plt.figure()
plt.plot( gamma*np.pi/180.0, ang_corr, label='data' )
plt.title( "Angular two-point correlation function", loc='center' )
plt.axes
plt.xlim(0,180)
plt.ylim(-400,400)
plt.xlabel('$\gamma \; [degree]$', fontsize='x-large')
plt.ylabel('$Cor(\gamma) \; [\mu K^2]$',fontsize='x-large')    
plt.axhline(y=0, color='grey', linewidth=0.01)
plt.legend()
plt.savefig(graficos+'angCorr_.png')
plt.show()
plt.close()
