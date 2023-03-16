"""
    Retrieve the spectra for the Aringer et al. (2009) model grid of carbon stars
    from their VizieR table:
    https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/503/913/spectra
"""
import numpy as np
from astropy.table import Table
from astropy import units, constants
import pyvo as vo
import collections

def get_VizieR_table():
    """
    Download the VizieR table containing the parameters of the grid
    as well as the names of the files containing the spectra
    """
    queryurl = 'http://tapvizier.u-strasbg.fr/TAPVizieR/tap'
    Vizier_service = vo.dal.TAPService(queryurl)
    table = 'J/A+A/503/913/spectra'
    query = """ select * from "{}" """.format(table)
    r = Vizier_service.search(query).to_table()
    return r

def get_model_spec(t = None, ModelNum = None, Teff = 3000.0, logg = 0.0, Mass = 1.0, Z = 1.0, \
                   C2O = 2.0, xi = 2.5, SpFileName = None, Dkpc = None):
    """
    Retrieve the model spectrum for the specified set of parameters. If
        no model found for the exact combination, find one with parameters
        closest to the ones specified.
        The proximity search is done in the following order: TEFF, LOGG, MASS, Z,
            C2O, XI.
    Arguments:
    t -- astropy Table containing at least one of the stellar parameters required
        to retrieve spectra. If relevant parameters are missing, the default values
        are used instead.
        One spectrum is returned for each row in the table.
    #
    The following arguments can be single values as in the default, or they can be
        lists or arrays (in which case their lengths should be identical and equal
        to the number of rows in the input table if provided).
    MODELNUM -- the number(s) of the model(s) from J/A+A/503/913/spectra. If provided,
        all other inputs are ignored and the specific model(s) is (are) retrieved.
    TEFF -- the effective temperature(s) of the model(s) in K.
    LOGG -- log10(surface gravity in cm/s2).
    MASS -- stellar mass(es) in Msun.
    Z -- metallicity (metallicities) in Zsun.
    C2O -- the C/O ratio(s).
    XI -- the microturbulent velocity (velocities).
    SPFILENAME -- the name(s) of the spectrum file. If provided, all other inputs
        are ignored and the specific model(s) is (are) retrieved.
    #
    Dkpc -- distances to the sources in kpc. If provided, the output spectrum is returned
        in Jy. If not set, the spectrum is returned in W/Hz.
    Output:
    spec -- astropy Table containing the model spectrum.
    """

    r = get_VizieR_table()

    if not(isinstance(ModelNum, type(None))):
        print('ModelNum specified!')
        if not(isinstance(ModelNum, collections.Sequence)):
            ModelNum = [ModelNum]
        # return spectra for the ModelNum values specified
        k = np.nonzero(r['Mod'] == ModelNum)[0]
        if len(k) == 0:
            raise ValueError("ModelNum not found in grid!")
        else:
            pass
    if not(isinstance(SpFileName, type(None))):
        # return spectra for the SpFileName values specified
        print('SpFileName specified!')
        pass

    pars = {'Teff': Teff, 'logg': logg, 'Mass': Mass, 'Z': Z, 'C/O': C2O, 'xi': xi}
    if isinstance(t, type(None)):
        # Generate a table with a single row and fill it with
        #   the default parameter values
        print('Input table not provided, generating one.')
        t = Table(list(np.array([list(pars.values())]).T), names = tuple(pars.keys()))
    else:
        print('Input table provided!')
        # Fill in the default values of any parameters that are missing in the
        #   input table
        for k, v in pars.items():
            if k not in t.colnames:
                t[k] = v

    # Retrieve the model spectra file names
    t['model'] = ' ' * 100
    for tt in t:
        rr = r.copy()
        # For each parameter, find the closest value in the grid
        #   to the value specified, then find all the models
        #   with that value.
        for p in pars:
            kk = np.abs(rr[p] - tt[p]).argmin()
            k = np.nonzero(rr[p] == rr[p][kk])[0]
            if len(k) != 0:
                rr = rr[k].copy()
            else:
                print("No matches found for parameter '{}'!".format(p))
        # If more than one model spectrum is found for the parameter
        #   combination, select the first one
        tt['model'] = rr['SpFileName'][0]
        # if len(rr) > 1:
        #     tt['model'] = rr['SpFileName'][0]
        # else:
        #     tt['model'] = rr['SpFileName']

    # Now, retrieve the spectra
    spectra_dir = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/503/913/spec/'
    spec = []
    for tt in t:
        spec.append(np.loadtxt(spectra_dir + tt['model'] + '.gz', unpack = False).T)

    spec = np.array(spec)
    wavelength = spec[:, 0, :] * units.Angstrom
    nulnu = spec[:, 2, :] * 1e-7 * units.W
    nu = wavelength.to('Hz', equivalencies = units.spectral())
    lnu = nulnu / nu

    if not(isinstance(Dkpc, type(None))):
        if len(Dkpc) != len(t):
            raise ValueError("Dkpc must have same length as input table!")
        print("Distances provided. Spectra will be in Jy.")
        if isinstance(Dkpc, list):
            Dkpc = np.array(Dkpc)
        fnu = (lnu / (4 * np.pi * (np.repeat(Dkpc[:, np.newaxis], lnu.shape[1], axis = 1) * units.kpc)**2)).to('Jy')
    else:
        print("Distances not provided. Spectra will be in W/Hz.")
        fnu = lnu

    t['wspec'] = wavelength.to('micron')
    t['trace'] = spec[:, 1, :]
    t['fspec'] = fnu

    return t

def get_logg(t):
    """
    Given the table of stellar parameters for the sample, return the log(g[cm/s^2])
        using log(g_sun[cm/s^2]) = 4.44.
    Arguments:
    t -- astropy Table object containing the stellar parameters in columns named
        'Lbol', 'Teff', and 'Mass'
    Output:
    logg -- log(g[cm/s^2])
    """
    Teff_sun = np.sqrt(np.sqrt((constants.L_sun / (4 * np.pi * constants.sigma_sb * constants.R_sun**2)).value))
    logR2 = np.log10(t['Lbol'].data) - 4 * np.log10(t['Teff'].data / Teff_sun)
    logg = np.log10(t['Mass'].data) - logR2 + np.log10((constants.G * constants.M_sun / constants.R_sun**2).to('cm s-2').value)

    return logg

def get_model_params():
    """
    Create an astropy Table object containing the stellar parameters from
        Table 2.1 in Itzel's thesis
    Two entries for Y CVn, because two mass estimates
    """
    ID = ['U Cam', 'R Scl', 'TX Psc', 'U Ant', 'U Hya', 'RT Cap', 'Y CVn']#, 'Y CVn']
    IRASPSC = ['03374+6229' , '01246-3248', '23438+0312', '10329-3918', \
               '10350-1307', '20141-2128', '12427+4542']#, '12427+4542']
    Lbol_Kastner = [11800.0, 6750.0, 5300.0, 7100.0, 2800.0, 10900.0, 10600.0]#, 10600.0]
    Lbol_bib = [13900.0, 7700.0, 6000.0, 6060.0, 2960.0, 8500.0, 9100.0]#, 9100.0]
    Teff = [2695.0, 2625.0, 3000.0, 2600.0, 2800.0, 2200.0, 2200.0]#, 2200.0]
    Mass = [3.6, 3.6, 2.0, 2.0, 2.4, 3.2, 2.8]#, 3.6]
    C2O = [1.30,1.34,1.07,1.44,1.043,1.10,1.087]
    Dkpc = [0.532, 0.266, 0.275, 0.268, 0.208, 0.291, 0.321]
    theta = [7.3, 19.5, 22, 42.5, 110, 94, 190]
    Rs = (np.array(Dkpc) * units.kpc * np.array(theta) * units.arcsec.to('radian')).to('pc')
    dR = np.array([0.0036, 0.0282, 0.0013, 0.015, 0.0123, 0.016, 0.0648]) * units.pc

    t = Table([ID, IRASPSC, Lbol_Kastner, Lbol_bib, Teff, Mass, C2O], \
              names = ('ID', 'IRASPSC', 'Lbol_Kastner', 'Lbol_bib', 'Teff', 'Mass', 'C/O'), \
              units = ('', '', 'L_sun', 'L_sun', 'K', 'M_sun', ''))

    t['Lbol'] = t['Lbol_Kastner']
    t['logg_Kastner'] = get_logg(t)
    t['logg_Kastner'].unit = ''
    t['Lbol'] = t['Lbol_Kastner']
    t['logg_Kastner'] = get_logg(t)
    t['logg_Kastner'].unit = ''
    t['logg'] = t['logg_Kastner'] # Just pick one
    t['logg'].unit = ''
    t['Rmin'] = Rs - dR / 2.0
    t['Rmax'] = Rs + dR / 2.0

    return t
