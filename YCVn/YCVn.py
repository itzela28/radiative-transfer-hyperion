"""
Model for U Cam with stellar parameters from get_SED.py
and photosphere spectrum from Aringer+2009_StellarSpectrum.ascii
"""
import numpy as np
from hyperion.util.constants import rsun, au, msun, sigma, pc, lsun
from astropy import units
from get_SED import *
from get_Aringeretal2009_model import *

#env_mass=[1.70422005e-05, 2.08293561e-05, 2.65100896e-05, 3.21908231e-05, 3.78715566e-05]
env_mass=[1.17506284e-05, 1.95843807e-05, 2.7418133e-05, 3.52518853e-05, 4.30856376e-05]

data = get_SED(IRASPSC = ['03374+6229',
'01246-3248',
'23438+0312',
'10329-3918',
'10350-1307',
'20141-2128',
'12427+4542'])
data['flux'].unit = units.Jy
data['flux_err'].unit = units.Jy

t = table2_1()
spec = get_spec(t,Dkpc=[0.532,0.266,0.275,0.268,0.208,0.291,0.321],C2O=[1.30,1.34,1.07,1.44,1.043,1.10,1.087])

all_data = join(data, spec, keys = 'IRASPSC', join_type = 'inner')

class PowerLawEnvelope():
    from hyperion.model import AnalyticalYSOModel
    from hyperion.densities.power_law_envelope import PowerLawEnvelope

    # Initalize the model
    m = AnalyticalYSOModel()
    
    # the lam and fnu can be used to plot the stellar photosphere
    # lam is in micron, fnu is in Jy
    m = get_stellar_params(m, dist=321*(np.sqrt(2))*units.pc)
    wave = all_data['wspec'][6]
    frec = (wave*units.um).to('Hz', equivalencies=units.spectral()).value
    fnu = all_data['fspec'][6]
    m.star.spectrum = (frec, fnu)
    import pdb; pdb.set_trace()

    # Add an envelope
    env = m.add_power_law_envelope()
    env.mass = (env_mass[0]) * msun          # Envelope mass
    env.rmin = (0.2529)*pc                  # Inner radius, 34au agb detached shell
    env.rmax = (0.3385)*pc      # Outer radius
    env.power = -2                 # Radial power
    #env.rho_0 = 2 only needed if env.mass is not given
    env.r_0 = env.rmin 
    env.dust = 'Cstar_DustProperties.hdf5'
    #import pdb; pdb.set_trace()
    # Use raytracing to improve s/n of thermal/source emission
    m.set_raytracing(True)

    # Use the modified random walk
    m.set_mrw(True, gamma=2.)

    # Set up grid
    m.set_spherical_polar_grid_auto(100, 1, 1)

    # Set up SED for 10 viewing angles
    sed = m.add_peeled_images(sed=True, image=False)
    sed.set_uncertainties(uncertainties=True)
    sed.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
    sed.set_wavelength_range(100, 0.3, 1500.)
    sed.set_track_origin('basic')

    # Set number of photons
    m.set_n_photons(initial=1e5, imaging=1e6,
                raytracing_sources=1e4, raytracing_dust=1e6)

    # Set number of temperature iterations
    m.set_n_initial_iterations(10)
    m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

    # Write out file
    m.write('YCVn_m2.rtin', overwrite=True)
    m.run('YCVn_m2.rtout', mpi=False, overwrite=True)
