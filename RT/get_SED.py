from astropy.table import Table, join
from astropy import units
import numpy as np
from hyperion.util.constants import rsun, au, msun, sigma, pc, lsun


def get_stellar_params(model = None, dist = 268 * units.pc):
    """
    Retrieve stellar parameters for U Ant from Dharmawardena et al. (2019)
    Arguments:
    model -- input Hyperion model object
    dist -- distance in pc to star (defaults to U Ant distance)
    Returns:
    m -- input Hyperion model object
    lam -- wavelength grid of stellar photosphere in micron
    fnu -- flux grid of stellar photosphere in Jy
    """
    # Photosphere spectrum
    specfile = 'Aringer+2009_StellarSpectrum.ascii'
    wav, fcont, nulnu = np.loadtxt(specfile, unpack = True)
    nu = (wav * units.Angstrom).to('Hz', equivalencies = units.spectral())
    # the model has nu * L_nu, so the following is technically L_nu
    # We'll scale it with the distance if set, if not use L_nu.
    #   Hyperion scales it knowing the total luminosity.
    #conversion should be 4 pi*distance**2
    fnu = ((nulnu * units.erg / units.s) / nu / (4 * np.pi * dist.to('m')**2)).to('Jy')

    lam = (wav * units.Angstrom).to('micron').value

    if isinstance(model, type(None)):
        m = None
        """
        # From Section 3.2.2 of Dharmawardena et al. 2019
        m.star.luminosity = 7000 * lsun
        mass = 2.0 # in Msun
        m.star.mass = mass * msun
        logg = -0.5 # surface gravity in cm/s^2
        logg_sun = 4.44 # solar surface gravity in the same units
        m.star.radius = np.sqrt(mass / 10**(logg - logg_sun)) * units.Rsun.to('cm')
        m.star.spectrum = (nu.value, fnu.value)
        # m.star.dstar = 260 * pc
        # Only set the temperature if the spectrum isn't input
        # m.star.temperature = 2600.
        """
    else:
        m = model

    return m, lam, fnu.value

def get_lams_and_fluxes():
    """
    Retrieve SEDs for the entire NESS sample
    Function for internal use only, use get_SED if you want
    to extract SEDs for a subsample of NESS sources!!!
    """
    ness = Table.read('NESS_SEDs_good.csv', format = 'csv')

    # connect the pivot wavelengths to the flux columns
    f = Table.read('filters.vot', format = 'votable')
    colnames = np.array(ness.colnames)
    k = np.nonzero([((c[-5:] == '_flux') or ('UnWISE_FW' in c)) and ('error' not in c) for c in colnames])[0]
    t_1 = Table([[fc.split('_')[-2] for fc in colnames[k]], colnames[k]], names = ('fluxcols', 'fullfluxcols'))
    f['fluxcols'] = [fn.split('/')[1] for fn in f['filterName']]
    t_2 = join(f['filterName', 'fluxcols', 'lpivot'], t_1, keys = 'fluxcols', join_type = 'inner')
    # fix the UnWISE column names, which are mapped to WISE names
    for c in ['WISE.W1_flux', 'WISE.W2_flux']:
        k = np.nonzero(t_2['fullfluxcols'] == c)[0][0]
        t_2['fullfluxcols'][k] = 'UnWISE_FW' + c[6]
    # ignore the WISE.W1 and WISE.W2 columns
    k = np.nonzero([('WISE.W1' not in c) & ('WISE.W2' not in c) for c in t_2['fullfluxcols']])[0]
    t_2 = t_2[k].copy()
    # sort according to wavelength
    ks = np.argsort(t_2['lpivot'])
    t_2 = t_2[ks].copy()

    # add a description to each flux/error column with the correct SVO filter name
    for tt in t_2:
        ness[tt['fullfluxcols']].description = tt['filterName']
        ness[tt['fullfluxcols'] + '_error'].description = tt['filterName']

    # pivot wavelength information goes into the meta attribute
    ness.meta['filterName'] = t_2['filterName'].data.data
    ness.meta['lpivot'] = t_2['lpivot'].data.data * t_2['lpivot'].unit

    # collect the fluxes into a matrix, ignoring the WISE.W1 and WISE.W2 points
    fluxes = []
    flux_errors = []
    for c in t_2['fullfluxcols']:
        fluxes.append(ness[c].data.data)
        flux_errors.append(ness[c + '_error'].data.data)
    fluxes = np.array(fluxes).T * units.Jy # makes an array of shape Nsources x Nwavelengths
    flux_errors = np.array(flux_errors).T * units.Jy

    return ness, fluxes, flux_errors

def get_SED(IRASPSC = ['03374+6229' , '01246-3248', '23438+0312', '10329-3918', \
               '10350-1307', '20141-2128', '12427+4542']):
    """
    Retrieve the SED of a list of sources from the NESS SEDs table.

    Arguments:
        IRASPSC -- list containing one or more IRASPSC identifiers (without the 'IRAS ' prefix)

    Returns:
        sed -- astropy Table object containing the following columns for each valid entry in IRASPSC:
            flux -- array containing the photometry in Jy
            flux_err -- array containing the photometric uncertainties in Jy
            The sed.meta attribute also contains the following columns:
                sed.meta['filterName'] -- the names of the broadband filters
                sed.meta['lpivot'] -- the pivot wavelengths of these filters in um
    """

    if IRASPSC == []:
        raise ValueError("")
    else:
        phot_table, f, df = get_lams_and_fluxes()
        nan = np.repeat(np.nan, f.shape[1]) # same length as number of filters
        nan = np.tile(np.nan, (len(IRASPSC), f.shape[1]))
        sed = Table([IRASPSC, nan, nan], names = ('IRASPSC', 'flux', 'flux_err'))
        k1 = np.nonzero(np.isin(IRASPSC, phot_table['IRASPSC']))[0]
        k2 = np.nonzero(np.isin(phot_table['IRASPSC'], IRASPSC))[0]
        sed['flux'][k1] = f[k2, :]
        sed['flux_err'][k1] = df[k2, :]
        sed.meta['filterName'] = phot_table.meta['filterName']
        sed.meta['lpivot'] = phot_table.meta['lpivot']

    return sed
