#blackbody with different radial power dependence
#rho~r^(-b), b=[0,2]

import numpy as np
from astropy import log as logger
import six
#from .core import Envelope
#class hyperion.densities.PowerLawEnvelope
class PowerLawEnvelope():
    from hyperion.util.constants import rsun, au, msun, sigma, pc, lsun
    from hyperion.model import AnalyticalYSOModel
    from hyperion.densities.power_law_envelope import PowerLawEnvelope

    # Initalize the model
    m = AnalyticalYSOModel()

    # Set the stellar parameters
    m.star.mass = 2. * msun
    m.star.radius = 309.7 * rsun
    m.star.temperature = 20000. 
    m.star.luminosity = 7000. * lsun
    # m.star.dstar=1000*pc  distance to star is given when the SED is ploted 

    # Add an envelope
    env = m.add_power_law_envelope()
    env.mass = (2.0e-5) * msun          # Envelope mass
    env.rmin = 34 * au                  # Inner radius, 4au regular abg, 34au agb detached shell
    env.rmax = 80 * au          # Outer radius 46 before
    env.power = 0                 # Radial power, values: 0,-1,-2
    #env.rho_0 = 2 only needed if env.mass is not given
    env.r_0 = 34 * au #34 au for detached shell agb
    env.dust = 'kmh_lite.hdf5'

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
    m.set_n_initial_iterations(5)
    m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

    # Write out file
    m.write('sed_rho_0_T20.rtin')
    m.run('sed_rho_0_T20.rtout', mpi=True)