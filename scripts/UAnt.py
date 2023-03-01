"""
Model for U Ant with stellar parameters from get_SED.py
and photosphere spectrum from Aringer+2009_StellarSpectrum.ascii
"""
import numpy as np
from hyperion.util.constants import rsun, au, msun, sigma, pc, lsun
from astropy import units
from get_SED import get_stellar_params

class PowerLawEnvelope():
    from hyperion.model import AnalyticalYSOModel
    from hyperion.densities.power_law_envelope import PowerLawEnvelope

    # Initalize the model
    m = AnalyticalYSOModel()

    # Set the stellar parameters
    # m.star.mass = 2. * msun
    # m.star.radius = 309.7 * rsun
    # m.star.temperature = 2800.
    # m.star.luminosity = 8000. * lsun
    # # m.star.dstar=260*pc
    
    # the lam and fnu can be used to plot the stellar photosphere
    # lam is in micron, fnu is in Jy
    m, lam, fnu = get_stellar_params(m)

    # Add an envelope
    env = m.add_power_law_envelope()
    env.mass = (2.e-5) * msun          # Envelope mass
    env.rmin = (0.0504205)*pc                  # Inner radius, 34au agb detached shell
    env.rmax =  (env.rmin)*(1.3)          # Outer radius
    env.power = -2                 # Radial power
    #env.rho_0 = 2 only needed if env.mass is not given
    env.r_0 = (0.0504205)*pc #34 au for detached shell agb
    #env.dust = 'kmh_lite.hdf5'
    env.dust = 'Cstar_DustProperties.hdf5' # 'kmh_lite.hdf5'

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
    sed.set_wavelength_range(10, 0.3, 1500.)
    sed.set_track_origin('basic')

    # Set number of photons
    m.set_n_photons(initial=1e5, imaging=1e6,
                raytracing_sources=1e4, raytracing_dust=1e6)

    # Set number of temperature iterations
    m.set_n_initial_iterations(5)
    m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

    # Write out file
    m.write('sed_UAnt_Crich.rtin')
    m.run('sed_UAnt_Crich.rtout', mpi=False)
