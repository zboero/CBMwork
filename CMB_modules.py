#c_luz = 2.99792458e+5		# Speed of light in Km/s.

# -------------------------------------------------------------------------------------
# It makes the conversion Velocity in Km/s to redshift.
# -------------------------------------------------------------------------------------

def v2z(x):		 					        
    ''' Converts the velocity of a gxs to redshift

    Parameters
    ----------
    x : float
        Velocity of the gxs from the catalog in Km/s. 
 
    Returns
    -------
    z : float
	The redshift via the relation z = v/c. 
    '''

    c_luz = 2.99792458e+5		# Speed of light in Km/s.
    z = x/c_luz

    return z


# -------------------------------------------------------------------------------------
# It makes the conversion redshift to Velocity in Km/s.
# -------------------------------------------------------------------------------------

def z2v(x):						      		
    ''' Converts the redshift of a gxs to velocity

    Parameters
    ----------
    x : float
        Redshift of the gxs from the catalog. 
 
    Returns
    -------
    v : float 
	The velocity via the relation v = c*z (in Km/s). 
    '''

    c_luz = 2.99792458e+5		# Speed of light in Km/s.
    v = x*c_luz
    return v


# -------------------------------------------------------------------------------------
# It computes the f_sky parameter.
# -------------------------------------------------------------------------------------

def f_sky(mask):     		
    ''' Computes the percentage of non-masked pixels with 
        respect to the total number of pixel in a map trought the
	"fsky" parameter.
	(Requires numpy library)

    Parameters
    ----------
    mask : np.array
        The mask array associated to a given map. 
 
    Returns
    -------
    fsky : float 
	The quotient (# of masked pixels)/(total # of pixels). 
    '''
    import numpy as np

    fsky = np.sum(mask) / ( len(mask) )

    return( fsky )


# -------------------------------------------------------------------------------------
# It reads the CMB maps and its associated mask.
# -------------------------------------------------------------------------------------

def Rd_map( filename, map_id, mask_id ):
    ''' Read the CMB maps from the .fits files.
	The function needs the "healpy" lybrary.

    Parameters
    ----------
    filename : string
        Redshift of the gxs from the catalog. 
    map_id : integer
        Numerical label of the maps within the .fits file.
        For instance: 0 ---> I, 
                      1 ---> Q,
                      2 ---> U.
    mask_id : integer
        Numerical label of the mask within the .fits file.
        For instance: 3 ---> Temperature mask,
                      4 ---> Polarization mask.

    Returns
    -------
    T_data : array
	Array of pixel (in RING ordering) with CMB temperature (in units of \muK) 
	or polarization data 
    T_mask : array
	Array of pixel (in RING ordering) with the CMB temperature or 
        polarization mask
    '''

    import numpy as np    			   
    import healpy as hp

    K_to_muK = 10**6		# Convertion factor between [Kelvin] and [microKelvin].
    
    CMB_T, header = hp.read_map( filename, nest=False, h=True, field=(map_id, mask_id) )
    T_data = CMB_T[0] * K_to_muK           	# Temperature (numpy) array 
    T_mask = CMB_T[1]         	      		# Mask (numpy) array with 0's and 1's.

    return( T_data, T_mask )



def Rd_map0(filename):			# It reads the CMB maps.
    ''' Rd_map0 is useful when we want to load a single map;
        for example, just the Intensity map or the common mask
        for all the CMB map cleaned of foregrounds.
	The function needs the "healpy" lybrary.

    Parameters
    ----------
    filename : string
        The path to the .fits file with the data. 

    Returns
    -------
    T_data : array
	Array of pixel (in RING ordering) with the map (in units of K if the 
        map is temperature)   
    '''

    import numpy as np    			   
    import healpy as hp

    T_data, header = hp.read_map( filename, nest=False, h=True, field=(0) )
    return(T_data) 

# -------------------------------------------------------------------------------------
# It choose the galaxy sample to be employed.
# -------------------------------------------------------------------------------------

def sample(gxs_type):
    ''' Choose the sample of galaxy according its morphologica type.

    Parameters
    ----------
    gxs_type : string
        The type of galaxy according to the Hubble classification.
	Only admit, 'Sa', 'Sb', 'Sc' and 'Sb+Sc'

    Returns
    -------
    Sample_in_use : Data frame
	Data frame with the info of gxs cominf from a catalog.
    '''
#    gxs_type_dict = {'Sa': 'Sa_sample.pickle', 'Sb': 'Sb_sample.pickle', 
#		     'Sc': 'Sc_sample.pickle', 'Sb+Sc': 'concat' } 
    import pickle as pkl

    root_home ='/home/zboero/'                                                  # Clemente y IATE
    fileProj  = root_home+'Projects/CMB/'					# Folder of work for the article
			
    if gxs_type == 'Sa': # ---- Sa Sample
        path_out_Sa = fileProj+'data/Galaxy_Samples/Sa_sample.pickle'   # Path to the Filename of Sa sample
        with open(path_out_Sa, 'rb') as f_Sa:
            Sample_in_use = pkl.load(f_Sa)

    elif gxs_type == 'Sb': # ---- Sb Sample
        path_out_Sb = fileProj+'data/Galaxy_Samples/Sb_sample.pickle'   # Path to the Filename of Sb sample
        with open(path_out_Sb, 'rb') as f_Sb:
            Sample_in_use = pkl.load(f_Sb)

    elif gxs_type == 'Sc': # ---- Sc Sample
        path_out_Sc = fileProj+'data/Galaxy_Samples/Sc_sample.pickle'   # Path to the Filename of Sc+Sd sample
        with open(path_out_Sc, 'rb') as f_Sc:
            Sample_in_use = pkl.load(f_Sc)

    elif gxs_type == 'Sb+Sc': # ---- Sb+Sc Sample
        path_out_Sb = fileProj+'data/Galaxy_Samples/Sb_sample.pickle'   # Path to the Filename of Sb sample
        with open(path_out_Sb, 'rb') as f_Sb:
            Sample_Sb = pkl.load(f_Sb)

        path_out_Sc = fileProj+'data/Galaxy_Samples/Sc_sample.pickle'    # Path to the Filename of Sc+Sd sample
        with open(path_out_Sc, 'rb') as f_Sc:
            Sample_Sc = pkl.load(f_Sc)

        Sample_in_use = pd.concat([Sample_Sb, Sample_Sc])

    return(Sample_in_use)


# -------------------------------------------------------------------------------------
# Makes a filter of the galaxy sample.
# -------------------------------------------------------------------------------------

def filt_sample(Sample, z_min, z_max, Bright_cut):
    ''' Makes a filtering of galaxies according two values in redshit.
        It also includes a filtering according if:
	Rext_cut   = False       # If True it cut galaxies with 10**(r_ext) > 120.0 arcsec.
	Bright_cut = True        # If True we take those with size > 8.5Kpc
        (revisar en otra oportunidad)
 
    Parameters
    ----------
    Sample : Data Frame
        The type of galaxy according to the Hubble classification.
	Only admit, 'Sa', 'Sb', 'Sc' and 'Sb+Sc'
    z_min : float
        The lower bound in redshift 
    z_max : float
        The upper bound in redshift

    Returns
    -------
    Sample_filt : Data frame
	Data frame with the final sample of galaxy after redshift filtering.
    '''
    import numpy as np

    asec2rad = np.pi/( 180.0*3600.0 )	      # Convertion factor between [arcseconds] and [radians].
    H0 = 100.0				      # Km/s/Mpc

    Sample_filt = Sample[ ( Sample['v'] < z2v(z_max) ) & ( Sample['v'] > z2v(z_min) ) ]

#    if Rext_cut == True:
#        Sample_filt = Sample_filt[ ( 10**(Sample_filt['r_ext']) < 120  ) ]

    if Bright_cut == True:
        Sample_filt = Sample_filt[ ( (Sample_filt['v']/H0) * (10**(Sample_filt['r_ext'])) * asec2rad > 0.0085  ) ]

    return(Sample_filt)


# -------------------------------------------------------------------------------------
# Makes a filter of the galaxy sample.
# -------------------------------------------------------------------------------------

def CMB_MwPlot(T_data, Tmin, Tmax, idMap, id_Sample):
    ''' Produces a Mollweide maps of the CMB maps
        Requires healpy
    Parameters
    ----------
    T_data : numpy array 
        The temperatue values (in units of \muK).
    Tmin : float
        The lower bound for T (values lesser than Tmin will be painted with the 
	minimun value of the colorscale)
    Tmax : float
        The upper bound for T (values grater than Tmax will be painted with the 
	maximun value of the colorscale)
    idMap : str 
        The label corresponding to the map (SMICA, SEVEM, NILC, Commander)
    id_Sample : str
        The label corresponding to the sample in use
        (Sa, Sb, Sc, Sb+Sc)   
    Returns
    -------
    The figure with the CMB map masked with the galaxy sample.
    '''

    import numpy as np    			   
    import healpy as hp
    import matplotlib.pyplot as plt

    root_home ='/home/zboero/'                                    # Clemente y IATE
    fileCMB_Maps = root_home+'Projects/CMB/graficos/'             # All the graphics and figures.
    colormap = 'inferno'#'viridis', 'plasma', 'magma', 'cividis'

    plt.close()
    plt.figure()
    
    hp.mollview( T_data, coord='G',unit='Temperature [$\mu$K]', xsize=800,\
                 title='CMB '+str(idMap)+' masked map '+'with '+str(id_Sample)+' gxs',\
                 cbar=True, cmap=colormap, min=Tmin, max=Tmax, badcolor='black')	#bgcolor='white',

#    plt.show()
    plt.savefig(fileCMB_Maps+str(idMap)+'_masked_with_'+str(id_Sample)+'_gxs.png')
    plt.close()


# -------------------------------------------------------------------------------------
# Make a measurement of the power spectrum from synthetic maps built with a given
#  input set of Cl's. This is to test the accuracy in which one recovers the
#  Cl's of input. 
# -------------------------------------------------------------------------------------

def Cl_recovery(Cl_input, n_maps, nside, lMax):
    ''' Given an initial spectrum "Cl_input" it produces a number "n_maps" of maps with
        resolution "nside" and from them it measure the spectrum in each map and
        average them to get final reconstrution of the original spectrum.
        The function uses synfast to prodce synthetic maps and anafast to measure the
        spectrum on the maps.
	
    Parameters
    ----------
    Cl_input : numpy array 
        The spectrum of input.
    n_maps : float
        The numbers of maps to generate and the total number of spectrum to average
    nside : float
        The nside resolution of the healpix map
    lMax : float
        The maximum multipole ell for the recovered spectrum.
    Returns
    -------
    The final averaged Cls, namely "Cl_output" and the error estimated trough the 
    standard deviation for each ell, namely "std_Cl".
    '''

    import numpy as np    			   
    import healpy as hp

    it_e = 3		                        # Number of iterations in the algorithm (standard is 3)
    Sum_Cl = np.zeros( lMax + 1 )
    Sum_Cl_sq = np.zeros( lMax + 1 )
    for i in range(0, n_maps):
        Map_syn = hp.synfast( Cl_input, nside, lmax=lMax, pol=False, alm=False, pixwin=False, fwhm=0.0 )

        Cls    = hp.anafast( Map_syn, nspec=None, lmax=3*nside-1, iter=it_e, alm=False)
        Cls_sq = Cls**2
	
        Sum_Cl    = Sum_Cl + Cls[:(lMax + 1)]
        Sum_Cl_sq = Sum_Cl_sq + Cls_sq[:(lMax + 1)]

    Cl_output = Sum_Cl/n_maps 
    Std_output = Sum_Cl_sq/n_maps - Cl_output**2

    return(Cl_output, Std_output)


# -------------------------------------------------------------------------------------
# Make a measurement of the power spectrum from synthetic maps built with a given
#  input set of Cl's. This is to test the accuracy in which one recovers the
#  Cl's of input. 
# -------------------------------------------------------------------------------------

def Cl_recovery_with_profile(Cl_input, n_maps, nside, lMax, hp_prof):
    ''' Given an initial spectrum "Cl_input" it produces a number "n_maps" of maps with
        resolution "nside"; the map is modified with the map "hp_prof" and from them, 
        the spectrum is measured in each map and average over the n_maps to get final 
        reconstrution of the original spectrum.
        The function uses synfast to prodce synthetic maps and anafast to measure the
        spectrum on the maps.
	
    Parameters
    ----------
    Cl_input : numpy array 
        The spectrum of input.
    n_maps : float
        The numbers of maps to generate and the total number of spectrum to average
    nside : float
        The nside resolution of the healpix map
    lMax : float
        The maximum multipole ell for the recovered spectrum.
    hp_prof: numpy array
        The map that modifies the synthetic ones genereated with Cl_input.
    Returns: list with two arrays
    -------
        The final averaged Cls, namely "Cl_output" and the error estimated trough the 
        standard deviation for each ell, namely "std_Cl".
    '''

    import numpy as np
    import healpy as hp

    it_e = 3		                        # Number of iterations in the algorithm (standard is 3)
    Sum_Cl = np.zeros( lMax + 1 )
    Sum_Cl_sq = np.zeros( lMax + 1 )
    Sum_Cl_Planck = np.zeros( lMax + 1 )
    Sum_Cl_Planck_sq = np.zeros( lMax + 1 )
    for i in range(0, n_maps):
        print('Mapa nro = ', i ) 
        Map_syn = hp.synfast( Cl_input, nside, lmax=lMax, pol=False, alm=False, pixwin=False, fwhm=0.0 )
        Map_mod = Map_syn + hp_prof

        Cls    = hp.anafast( Map_mod, nspec=None, lmax=3*nside-1, iter=it_e, alm=False)
        Cls_sq = Cls**2
        Sum_Cl    = Sum_Cl + Cls[:(lMax + 1)]
        Sum_Cl_sq = Sum_Cl_sq + Cls_sq[:(lMax + 1)]

        Cls_Planck = hp.anafast( Map_syn, nspec=None, lmax=3*nside-1, iter=it_e, alm=False)
        Cls_Planck_sq = Cls_Planck**2
        Sum_Cl_Planck = Sum_Cl_Planck + Cls_Planck[:(lMax+1)]
        Sum_Cl_Planck_sq = Sum_Cl_Planck_sq + Cls_Planck_sq[:(lMax+1)]
        
    Cl_output = Sum_Cl/n_maps 
    Std_output = Sum_Cl_sq/n_maps - Cl_output**2
    Cl_output_Planck = Sum_Cl_Planck/n_maps 
    Std_output_Planck = Sum_Cl_Planck_sq/n_maps - Cl_output_Planck**2


    return(Cl_output, Std_output, Cl_output_Planck, Std_output_Planck)


# -------------------------------------------------------------------------------------
# Making the profile that modifies the CMB maps
# -------------------------------------------------------------------------------------

def hp_profile(gxs_sampl, nside, lMax):
    ''' Given an initial spectrum "Cl_input" it produces a number "n_maps" of maps with
        resolution "nside"; the map is modified with the map "hp_prof" and from them, 
        the spectrum is measured in each map and average over the n_maps to get final 
        reconstrution of the original spectrum.
        The function uses synfast to prodce synthetic maps and anafast to measure the
        spectrum on the maps.
	
    Parameters
    ----------
    gxs_sampl : pandas dataframe
        The dataframe with the galaxy sample.
    nside : float
        The nside resolution of the healpix map
    lMax : float
        The maximum multipole ell for the recovered spectrum.
    Returns
    -------
    hp_prof : array
        The map that modifies the synthetic ones genereated with Cl_input.
    '''

    import numpy as np
    import healpy as hp

    min2deg = 1./60.

    nside_high = 2048
    nside_low  = 512
    nside_nube = 256
    hp_prof_high = np.zeros( hp.nside2npix( nside_high ) )  #Array con longitud igual al nro de pixeles en resolución 2048
                                                            # Este array con el valor del perfil impuesto a los pixeles 
                                                            # alrededdor de las gxs es que vamos a sumar a los mapas 
    hp_prof_low  = np.zeros( hp.nside2npix( nside_low ) )   #Array con longitud igual al nro de pixeles en resolución 512
                                                            # cumple la misma función que el anterior
    hp_prof_nube = np.zeros( hp.nside2npix( nside_nube ) )  #Array con longitud igual al nro de pixeles en resolución 128
                                                            # cumple la misma función que el anterior

    vec = hp.ang2vec( gxs_sampl['b'] * np.pi/180.0, gxs_sampl['l'] * np.pi/180.0, 
                      lonlat=False )                        # Array con las posiciones de las galaxias, es independiente de los
                                                            #  pixeles y de la resolución de los mapas...

    rpr = np.logspace( start= np.log10(50), stop= np.log10(70000), num= 30 )/60.0     # Array con los valores del bineado radial 
                                                                                      #  que se usa para hacer los perfiles. 
                                                                                      # La escala está en minutos de arco.

    rpfit    = rpr[0:27]                                    # Esto define donde va a ser no nulo el perfil impuesto...
                                                            #  los bines del 27 al 30 de rpr se mantienen a zero, los demás no.  
#    rpfit_gr = rpr[0:27]*0.0166667                          # Se pasa la escala de los bines a grados...
    rpfit_gr = rpr[0:27] * min2deg                          # Se pasa la escala de los bines a grados... 
##### Debajo es para 7 log tehta - 5 ......
#    rpfit = rpr[0:24]                                       #corto en cero para el perfil modelado
#    rpfit_gr = rpr * 0.0166667                              #esto está en grados!!! 
                                                            #lo paso a grados porque el ajuste esta en grados

    # defino el perfil a aplicar:
    #alpha = 1.0 # Ponemos alfa 3 para amplificar el efecto
    #Mfit = np.zeros( len(rpfit_gr) )
    # Debajo hacemos la siguiente asignación: cuando estamos a 1 grado o menos le asignamos 
    #   alfa * ( 7 * log10( theta ) - 5.0 );  cuando estamos a más de 1 grado  le asignamos
    #   alfa * ( 5 * log10( theta ) - 5.0 );
    Mfit = np.zeros( len(rpfit_gr) )
    Mfit[ np.log10(rpfit_gr) <= 0 ] = 7 * np.log10( rpfit_gr[ np.log10( rpfit_gr ) <= 0.0 ] ) - 5.0
    Mfit[ np.log10(rpfit_gr)  > 0 ] = 5 * np.log10( rpfit_gr[ np.log10( rpfit_gr )  > 0.0 ] ) - 5.0
    Mfit[ rpfit_gr < 0.3 ]          = -8.12898859
    
    
    for k in range( len(rpfit) ):
        print('Radio =', rpfit[k], k )
        rp = rpfit[k]
        rp_ant = rpfit[k-1]

        Kd=[]

        if ( rp > 50 ):
            nside_sel = nside_low
            nside_nube = nside_nube  

        if ( rp < 50 ):
            nside_sel = nside_high

        for i in range( len(vec) ):
            print('Galaxy ------------------ ', i )
        #for i in range( 0,1 ): 
            listpix = hp.query_disc( nside_sel, vec[i], rp*60.0/206264.99992965, inclusive=True,
                                     fact=4, nest=False )                            #radianes (pasando a seg de arco)
            listpix_prev = hp.query_disc( nside_sel, vec[i], rp_ant*60.0/206264.99992965, inclusive=True,
                                     fact=4, nest=False ) 
            listpix_nube = hp.query_disc( nside_nube, vec[i], rp*60.0/206264.99992965, inclusive=True,
                                     fact=4, nest=False )                            #radianes (pasando a seg de arco)
            listpix_prev_nube = hp.query_disc( nside_nube, vec[i], rp_ant*60.0/206264.99992965, inclusive=True,
                                     fact=4, nest=False ) 
            if ( k == 0 ):
                listpix_prev = []
                listpix_prev_nube = []

            listpix_mask = []
            listpix_mask_nube = []

            for j in range( len(listpix) ):
                if ( listpix[j] not in listpix_prev ):
                    listpix_mask.append( listpix[j] )

            for j in range( len(listpix_nube) ):
                if ( listpix_nube[j] not in listpix_prev_nube ):
                    listpix_mask_nube.append( listpix_nube[j] )

            if ( rp > 50 ):
                Rdm_noise = np.random.uniform( -0.5, 0.5  , size= len(listpix_mask_nube) )
                #hp_prof_low[listpix_mask] = Mfit[k] + hp_prof_low[listpix_mask]
                hp_prof_low[listpix_mask] = Mfit[k] + hp_prof_low[listpix_mask] 
                hp_prof_nube[listpix_mask_nube] = hp_prof_nube[listpix_mask_nube] +\
                    ( 20.0 * np.sqrt(abs(Mfit[k])) ) * ( Rdm_noise )
#                hp_prof_nube[listpix_mask_nube] = Mfit[k] + hp_prof_nube[listpix_mask_nube] +\
#                                                            ( 30.0 * np.sqrt(abs(Mfit[k])) ) * ( Rdm_noise - 0.5  )
            if ( rp < 50 ):
                #Rdm_noise = np.random.uniform( 0.0, 1.0  , size= len(listpix_mask) )
                hp_prof_high[listpix_mask] = Mfit[k] + hp_prof_high[listpix_mask]   
                #hp_prof_high[listpix_mask] = Mfit[k] + hp_prof_high[listpix_mask] + ( 10.0 * np.sqrt(abs(Mfit[k])) ) * ( Rdm_noise - 0.5  )  

    hp_prof_up = hp.pixelfunc.ud_grade( hp_prof_low, nside_high, pess=True, order_out='RING', power=0.0 )     #le subo la resolución al mapa de 512
    hp_prof_nube_up = hp.pixelfunc.ud_grade( hp_prof_nube, nside_high, pess=True, order_out='RING', power=0.0 ) #le subo la resolución al mapa de 128
    hp_prof    = hp_prof_high + hp_prof_up + hp_prof_nube_up  # mapa con los perfiles sobre todas las galaxias en todo el rango de escala
    mean_map   = np.mean( hp_prof )
    hp_prof    = hp_prof - mean_map/ hp.nside2npix( nside_high )

    return hp_prof #hp_prof_low, hp_prof_high




# -------------------------------------------------------------------------------------
# Making the profile that modifies the CMB maps with a random smaple of gxs
# -------------------------------------------------------------------------------------

def hp_profile_rdm(lon_array, lat_array, nside, lMax):
    ''' Given an initial spectrum "Cl_input" it produces a number "n_maps" of maps with
        resolution "nside"; the map is modified with the map "hp_prof" that is built 
	with the random sample of galaxies; and from them, 
        the spectrum is measured in each map and average over the n_maps to get final 
        reconstrution of the original spectrum.
        The function uses synfast to produce synthetic maps and anafast to measure the
        spectrum on the maps.
	
    Parameters
    ----------
    lon_array : integer
        The number of galaxies in the random sample.
    lat_array : float
        The nside resolution of the healpix map
    nside : float
        The maximum multipole ell for the recovered spectrum.
    lMax : float
        The maximum multipole ell for the recovered spectrum.
    Returns
    -------
    hp_prof : array
        The map that modifies the synthetic ones genereated with Cl_input.
    '''

    import numpy as np
    import pandas as pd
    import healpy as hp

    nside_high = 2048
    nside_low = 512
    hp_prof_high = np.zeros( hp.nside2npix( nside_high ) )  #perfil de resolución 2048 hasta aprox 1 grado
    hp_prof_low  = np.zeros( hp.nside2npix( nside_low ) )   #perfil de resolución 512 para el resto de la escala

#    lon_array = np.random.uniform(  0.0, 360.0, size= n_gxs )
#    costheta  = np.random.uniform( -1.0, 1.0  , size= n_gxs )
#    lat_array = np.arccos( costheta )

    dict_gxs_sampl = {'l': lon_array, 'b': lat_array }
    gxs_sampl = pd.DataFrame( data= dict_gxs_sampl )

    vec = hp.ang2vec( gxs_sampl['b'], gxs_sampl['l'] * np.pi/180.0, lonlat=False ) #vector con la posicion de las galaxias

    rpr = np.logspace( start= np.log10(50), stop= np.log10(70000), num= 30 )/60.0      #esto está en minutos!!! 
                                                                                       #bines en eje x del perfil real

    rpfit = rpr[0:27]
    rpfit_gr = rpr[0:27]*0.0166667
##### Debajo es para 7 log tehta - 5 ......
#    rpfit = rpr[0:24]                                       #corto en cero para el perfil modelado
#    rpfit_gr = rpr * 0.0166667                              #esto está en grados!!! 
                                                            #lo paso a grados porque el ajuste esta en grados
    # defino el perfil a aplicar:
    alpha = 1.0 # Ponemos alfa 3 para amplificar el efecto
    Mfit = np.zeros( len(rpfit_gr) )
#    Mfit = (7.0 * np.log10( rpfit_gr ) - 5.0) * alpha	    
    Mfit[ np.log10(rpfit_gr) <= 0 ] = alpha * ( 7 * np.log10( rpfit_gr[ np.log10(rpfit_gr) <= 0.0 ] ) - 5.0 )
    Mfit[ np.log10(rpfit_gr)  > 0 ] = alpha * ( 5 * np.log10( rpfit_gr[ np.log10(rpfit_gr)  > 0.0 ] ) - 5.0 )

    for k in range( len(rpfit) ):
        #print(k)
        rp = rpfit[k]
        rp_ant = rpfit[k-1]

        Kd=[]

        if ( rp > 50 ):
            nside_sel = nside_low

        if ( rp < 50 ):
            nside_sel = nside_high

        for i in range( len(vec) ):
            listpix = hp.query_disc( nside_sel, vec[i], rp*60.0/206264.99992965, inclusive=True,
                                     fact=4, nest=False )                            #radianes (pasando a seg de arco)
            listpix_prev = hp.query_disc( nside_sel, vec[i], rp_ant*60.0/206264.99992965, inclusive=True,
                                     fact=4, nest=False ) 
            if ( k == 0 ):
                listpix_prev = []

            listpix_mask = []

            for j in range( len(listpix) ):
                if ( listpix[j] not in listpix_prev ):
                    listpix_mask.append( listpix[j] )

            if ( rp > 50 ) :
                hp_prof_low[listpix_mask] = Mfit[k] + hp_prof_low[listpix_mask]    #el ajuste está en grados
            if ( rp < 50 ):
                hp_prof_high[listpix_mask] = Mfit[k] + hp_prof_high[listpix_mask]  #el ajuste está en grados

    hp_prof_up = hp.pixelfunc.ud_grade( hp_prof_low, nside_high, pess=True )     #le subo la resolución al mapa de 512
    hp_prof    = hp_prof_high + hp_prof_up  # mapa con los perfiles sobre todas las galaxias en todo el rango de escala
    mean_map   = np.mean( hp_prof )
    hp_prof    = hp_prof - mean_map

    return hp_prof

# -------------------------------------------------------------------------------------
# Averages in bins of lenght N.
# -------------------------------------------------------------------------------------

def Bin_average( Cls, bin_n ):
    ''' Performs an average in bins of width "bin_n" of the array "Cls"

    Parameters
    ----------
    Cls : array
        The array with the powerspectrum Cl's.
    bin_n : integer
        The width of the bin chosen to average.
	Let us note that bin_n mut be an odd integer.

    Returns
    -------
    Cls_averaged : array
        The averaged power spectrum in bins.
    '''

    import numpy as np 

    Cls_averaged = Cls.copy()
    j_ini = (bin_n - 1)/2
    Cls_ite = np.zeros( len(Cls_bin) )

    for k in range( 0, j_ini + 1 ) :
        Cls_ite[ j_ini : - j_ini ] = Cls_ite[ j_ini : - j_ini ] + Cls[ (j_ini - k) : - (j_ini + k) ]
           
    Cls_averaged[ j_ini : j_end ] = Cls_ite[ j_ini : j_end ] / bin_n 

    return(Cls_averaged)
    


# -------------------------------------------------------------------------------------
# Measuring the power spectrum on a masked mask
# -------------------------------------------------------------------------------------

def Pseudo_Cls( CMB_map, CMB_mask, ):
    ''' Performs a measurement of the power spectrum of a masked map
            employing the so-called Pseudo Cl's technique.
            (It needs healpy camb libraries)

    Parameters
    ----------
    CMB_map : array
        The whole sky temperature map without mask.
    CMB_mask : array
        The array defining the mask (only containing 1's and 0's)

    Returns
    -------
    Cls : array
        The power spectrum measured with the method implemented here.
    '''
    
    import numpy as np
    import healpy as hp
    import camb 

    K_to_muK = 10**6		# Convertion factor between [Kelvin] and [microKelvin].
    Nside = hp.pixelfunc.get_nside(CMB_map) 

    M_l_lp = camb.mathutils.scalar_coupling_matrix(CMB_mask, 2499)     # Coupling matrix computed with the mask	
    M_inv = np.linalg.inv(M_l_lp)                                      # Inverse of the coupling matrix

    bare_Cl = hp.anafast( CMB_map, lmax=3*Nside-1, iter=3, 
                          alm=False, pol=False)                        # Cl's of the non-masked map
    CMB_mapMasked = np.where( CMB_map, T_data, hp.UNSEEN )             # Masked map
    Pseudo_Cls = hp.anafast( CMB_mapMasked, lmax=3*Nside-1, iter=3, 
                             alm=False, pol=False )                    # Cl's of the masked map; i.e. Pseudo Cl's
    Cls = np.dot( M_inv, Pseudo_Cls[0:2500] )                          # We revert the system Pseudo = M * Cl's;
                                                                       # so that this is our estimate of the Cls  

    return(Cls)

# -------------------------------------------------------------------------------------
# Rotation of a given map
# -------------------------------------------------------------------------------------

def Rot_phi( CMB_map, phi ):
    ''' Performs a rotation of a given map "CMB_map" along an axis
            oriented in the direction z.
            (It needs healpy library)

    Parameters
    ----------
    CMB_map : array
        The whole sky temperature map without mask.
    phi : float
        The rotation angle in degrees

    Returns
    -------
    Cls : array
        The power spectrum measured with the method implemented here.
    '''

    import healpy as hp

    Rot = hp.Rotator( rot=(phi, 0.0, 0.0), coord=None, inv=None, deg=True, eulertype='X')
    rotated_map = Rot.rotate_map_alms( CMB_map )

    return( rotated_map )


# -------------------------------------------------------------------------------------
# The power spectrum of the galaxy sample D_l
# -------------------------------------------------------------------------------------

def PS_gxs_sampl( gxs_sampl ):
    ''' Computes the power spectrum of the galaxy distribution.

    Parameters
    ----------
    gxs_sampl : pandas dataframe
        The sample of gxs.

    Returns
    -------
    Dl : array
        The power spectrum measured.
    '''

    import numpy as np
    import pandas as pd
    from scipy.special import eval_legendre

    th  = gxs_sampl['b'] * np.pi/180.0          # Latitude in radians 
    ph  = gxs_sampl['l'] * np.pi/180.0          # Longitude in radians

    cos_th = np.cos(th)  
    cos_ph = np.cos(ph)
    sin_th = np.sin(th)
    sin_ph = np.sin(ph)
    cos_ph_ph = np.outer( cos_ph, cos_ph ) + np.outer( sin_ph, sin_ph )

    cos_gamma = np.outer( cos_th, cos_th ) + np.outer( sin_th, sin_th ) * cos_ph_ph
    
    Dl = []
    for l in range( 0, 250 ):
        P_l_cos_gamma = eval_legendre( l, cos_gamma)
        Dl_idx = np.sum( P_l_cos_gamma )
        Dl.append( Dl_idx)


    return( Dl )


# -------------------------------------------------------------------------------------
# Distribution of gxs on a Mollweide projection
# -------------------------------------------------------------------------------------

def Gxs_MwPlot(l, b, title):
    ''' Produces a Mollweide map with the gxs distribution inprinted on it.
        It also allows for a mask placed on the sky.
    Parameters
    ----------
    l : array
        The longitude of galaxies in deg (0,360).
    b : pandas dataframe
        The latitude of galaxies in deg (-90,90).
    Mask : array
        Array with the mask of the full sky map
    title : str
        Title of the graphic
    id_Sample : str
        The label corresponding to the sample in use
        (Sa, Sb, Sc, Sb+Sc)   
    Returns
    -------
    The figure with the CMB map masked with the galaxy sample.
    '''

    import numpy as np    			   
    import matplotlib.pyplot as plt

    root_home ='/home/zboero/'                                    # Clemente y IATE
    fileout = root_home+'Projects/CMB/graficos/'             # All the graphics and figures.
#    colormap = 'inferno'#'viridis', 'plasma', 'magma', 'cividis'

#    ''' RA, Dec are arrays of the same length.
#    RA takes values in [0,360), Dec in [-90,90], which represent angles in degrees.
#    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
#    title is the title of the figure.
#    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
#    '''
    RA  = l #gxs_sampl['l']
    Dec = b #gxs_sampl['b']
    org =0
    projection ='mollweide' 
    x = np.remainder(RA+360-org,360) # shift RA values
    ind = x > 180
    x[ind] -= 360    # scale conversion to [-180, 180]
    x = -x    # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder( tick_labels+360+org, 360 )
    fig = plt.figure( figsize=(10, 5) )
    ax = fig.add_subplot( 111, projection=projection ) #, axisbg ='LightCyan'
    ax.scatter( np.radians(x), np.radians(Dec), s=6 )  # convert degrees to radians
    ax.set_xticklabels( tick_labels )                  # we add the scale on the x axis
    ax.set_title( title )
    ax.title.set_fontsize(15)
    ax.set_xlabel("l")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("b")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)


##################
###### Correlation function Planck and class
##################
###### GRAFICO 1

def Planck_Corr(Dl, gamma, lMax):
    import numpy as np
    from scipy.special import eval_legendre

    sum=0
    for n in range(2,lMax):
        sum= sum + (2*n+1)/(4*np.pi)*(2*np.pi/(n*(n+1))*Dl[n-2])*eval_legendre(n,np.cos(gamma))
    return sum

def Planck_CorrBF(Dlbf, gamma, lMax):
    import numpy as np
    from scipy.special import eval_legendre

    sum=0
    for n in range(2,lMax):
         sum= sum + (2*n+1)/(4*np.pi)*(2*np.pi/(n*(n+1))*Dlbf[n-2])*eval_legendre(n,np.cos(gamma))
    return sum


def class_Corr(Dlclass, gamma, lMax):
    import numpy as np
    from scipy.special import eval_legendre

    sum=0
    for n in range(2,lMax):
         sum= sum + (2*n+1)/(4*np.pi)*(2*np.pi/(n*(n+1))*Dlclass[n-2]*1e+13)*eval_legendre(n,np.cos(gamma))
    return sum




