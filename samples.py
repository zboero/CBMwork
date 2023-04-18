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
# This program is used to genereate samples of galaxies from the 2MRS catalogue
#
# ------------------------------------------------------------------------------

host = 'local'

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


# Samples of galaxies....
#data_2MRS   = data+'2MASS/2mrs_1175_done.dat'                      # 2MRS catalog
data_2MRS   = data+'2mrsgalaxies_z0.001_z0.02_correctedxWISE.dat'
#path_out_Sb = data+'Galaxy_Samples/Sb_sample.pickle'   # Path to the Filename of Sb sample
#path_out_Sc = data+'Galaxy_Samples/Sc_sample.pickle'   # Path to the Filename of Sc+Sd sample
#data_test   = data+'tstcat.txt'

#common_mask = fileMask+'COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'

####################################################################################
####################################################################################

import numpy  as np
import pandas as pd
import healpy as hp
import pickle as pkl

head_2MRS = ['RAdeg', 'DECdeg', 'l', 'b', 'k_c', 'h_c', 'j_c', 'k_tc', 'h_tc', 'j_tc',\
                               'e_k', 'e_h', 'e_j', 'e_kt', 'e_ht', 'e_jt', 'e_bv', 'r_iso', 'r_ext',\
                               'b/a', 'flgs', 'type', 'ts', 'v', 'e_v', 'c']

head_smp1 = ['l', 'b', 'rad', 'd5', 'z', 'type', 'RAdeg', 'DECdeg', 'W1mag',
                               'W2mag', 'W3mag', 'e_W1mag'   'e_W2mag'   'e_W3mag']

df_2MRS = pd.read_table(data_2MRS, skiprows=10,\
                        names= head_smp1, sep="\s+",\
                        index_col=False)

# Redshift cut ....
z_min = 0.004
z_max = 0.02
c_luz = 299792.458                         # in km/s

#df_gxs      = df_2MRS[ df_2MRS["v"] > c_luz * z_min ]
#df_gxs      = df_gxs[ df_gxs["v"] < c_luz * z_max ]
df_gxs      = df_2MRS[ df_2MRS["z"] > z_min ]
df_gxs      = df_gxs[ df_gxs["z"] < z_max ]

# Size cut ....
sz_min = 0.0                               # in kpc
sz_max = 20.0                              # in kpc
asec2rad = np.pi/( 180.0*3600.0 )          # Convertion factor between [arcseconds] and [radians].
H0 = 100.0                                 # in km/s/Mpc

#df_gxs = df_gxs[ ( df_gxs['v']/H0 * 10**( df_gxs['r_ext'] ) * asec2rad > (sz_min / 1000) ) ]
#df_gxs = df_gxs[ ( df_gxs['v']/H0 * 10**( df_gxs['r_ext'] ) * asec2rad < (sz_max / 1000) ) ]
df_gxs = df_gxs[ df_gxs['rad'] < 20.0  ]
df_gxs = df_gxs[ df_gxs['rad'] > 0.0   ]

# Filtering gxs type ....
# 1: Sa, 2: Sab, 3: Sb, 4: Sbc, 5: Sc, 6: Scd, 7: Sd, 8: Sdm, 9: Sm, 20: (S..., Sc-Irr, Unclassified Spiral)
# -1: (L+,SO+), -2: (L,SO), -3: (L-,SO-), -4: (E,SO), -5: (E,E dwarf), -6: (Compact Ellipt), -7: (Unclasiffied Ell)
Sp_list = np.array([3.0,4.0,5.0,6.0,7.0,8.0])#["1","2","3","4","5","6","7","8","9"]

kfail = 0
df_gxs = df_gxs.reset_index(drop=True)
for i in range( len(df_gxs) ):
    try:
        if df_gxs['type'][i] not in Sp_list:
            df_gxs = df_gxs.drop([i])
        elif df_gxs['type'][i] in Sp_list:
            continue
            
    except:
        print('Problem in index', i)
        kfail += 1                # kfail +1

print('Remaining spirals, ????: ', len(df_gxs), kfail )

d5_ref = 3.0
df_gxs = df_gxs[ df_gxs["d5"] < d5_ref ]


# # We compute the distances to the 5th galaxy...
# sigma5 = np.zeros( len(df_gxs) )                                               # array with len() equal to the sample

# for i in range( len(df_gxs) ):
#     df_dist1  = df_2MRS[ abs( df_2MRS['v'] - df_gxs['v'].iloc[i] ) < 500.0 ]
#     dist      = hp.rotator.angdist( [ np.deg2rad( 90.0 - df_gxs['b'].iloc[i] ), np.deg2rad( df_gxs['l'].iloc[i] ) ],
#                                     [ np.deg2rad( 90.0 - df_dist1['b'] ),  np.deg2rad( df_dist1['l'] )  ], lonlat=False )
#     dist_sort = np.sort( dist.tolist() )
#     sigma5[i] = dist_sort[4]


# d5_ref = np.deg2rad(3.0)
# #df_gxs = df_gxs[ sigma5 < st.median(sigma5)   ]
# df_gxs = df_gxs[ sigma5 < d5_ref ]


# # Save the final sample ...
# RA        = df_gxs["RAdeg"]
# DEC       = df_gxs["DECdeg"]
# l         = df_gxs["l"]
# b         = df_gxs["b"]
# r_ext     = df_gxs["r_ext"]
# flgs      = df_gxs["flgs"]
# v         = df_gxs["v"]
# dist5gxs  = np.array(sigma5[ sigma5 < d5_ref])

# stack = np.column_stack( ( RA, DEC, l, b, r_ext, v, dist5gxs ) )
# header_dat = " 'RAdeg', 'DECdeg', 'l', 'b', 'r_ext', 'v', 'dist5gxs'"

# output_file = data+'gxs_sampleFrode_allSp.dat'
# #output_file = output_file+'testmodel2_extendedz0.05_nside512_muK_ringordering.dat'
# np.savetxt( output_file, stack, header=header_dat, \
#             fmt=['%18.12E', '%18.12E', '%18.12E', '%18.12E', '%18.12E', '%18.12E', '%18.12E',]
#             )


