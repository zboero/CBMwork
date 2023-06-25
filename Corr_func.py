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


