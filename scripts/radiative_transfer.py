"""
"""
import numpy as np
from hyperion.util.constants import rsun, au, msun, sigma, pc, lsun
from astropy import units 
# from get_SED import get_stellar_params
from get_SED import *
from get_Aringeretal2009_model import *

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
    
    data = get_SED(IRASPSC = ['03374+6229',
    '01246-3248',
    '23438+0312',
    '10329-3918',
    '10350-1307',
    '20141-2128',
    '12427+4542'])
    
    spec = get_spec(table2_1(),
        Dkpc=[0.532,0.266,0.275,0.268,0.208,0.291,0.321],C2O=[1.30,1.34,1.07,1.44,1.043,1.10,1.087])

    #3 es para U Ant. Necesita ser automatizado para escoger una fuente
    # the lam and fnu can be used to plot the stellar photosphere
    # lam is in micron, fnu is in Jy
    i = 3
    m.star.luminosity = spec['Lbol_Kastner'][i] * spec['Lbol_Kastner'].unit.to('erg/s')
    # m.star.luminosity = (spec['Lbol_Kastner'].quantity)[i] # 7000 * lsun
    m.star.mass = spec['Mass'][i] * spec['Mass'].unit.to('g')
    logg = spec['logg'][i] # surface gravity in cm/s^2
    logg_sun = 4.44 # solar surface gravity in the same units
    m.star.radius = np.sqrt(m.star.mass / 10**(logg - logg_sun))
    wave = (spec['wspec'].quantity)[i]
    m.star.spectrum = (wave.to('Hz', equivalencies=units.spectral()).value, 
    spec['fspec'][i])
    m, lam, fnu = get_stellar_params(m)


    env_mass=[1.11757156e-05, 1.3659208e-05, 1.73844465e-05,
		2.11096851e-05, 2.48349236e-05]

    # Add an envelope    
    #env = m.add_power_law_envelope()
             # Envelope mass
    #env.rmin =  (spec['Rs'][i]-spec['dR'][i]/2) * spec['Rs'].unit.to('cm')                  # Inner radius, 34au agb detached shell
    #env.rmax =  (spec['Rs'][i]+spec['dR'][i]/2) * spec['Rs'].unit.to('cm')        # Outer radius
    #env.mass = env_mass[4] * units.solMass.to('g') 
    #env.mass = (tau/kappa)*(rmax*rmin) * units.solMass.to('g') 
    #write env.mass in terms of opt. depth, rmax,rmin, opacity
    #env.power = -2                 # Radial power
    #env.rho_0 = 2 only needed if env.mass is not given
    #env.r_0 = env.rmin #34 au for detached shell agb
    #env.dust = 'kmh_lite.hdf5'
    #env.dust = 'Cstar_DustProperties.hdf5' # 'kmh_lite.hdf5'
    #print(env.rmin, env.rmax, env.mass)
    env = m.add_power_law_envelope()
    env.mass = (1.12e-5) * msun          # Envelope mass
    env.rmin = (0.0504205)*pc                  # Inner radius, 34au agb detached shell
    env.rmax =  (env.rmin)*(1.3)          # Outer radius
    env.power = -2                 # Radial power
    #env.rho_0 = 2 only needed if env.mass is not given
    env.r_0 = (0.0504205)*pc #34 au for detached shell agb
    #env.dust = 'kmh_lite.hdf5'
    env.dust = 'Cstar_DustProperties.hdf5'


    # Use raytracing to improve s/n of thermal/source emission
    m.set_raytracing(True)

    # Use the modified random walk
    m.set_mrw(True, gamma=2.)

    # Set up grid
    m.set_spherical_polar_grid_auto(100, 1, 1)

    # Set up SED for 10 viewing angles
    waverange = np.array([0.3, 1500.]) * units.micron
    sed = m.add_peeled_images(sed=True, image=False)
    sed.set_uncertainties(uncertainties=True)
    sed.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
    sed.set_wavelength_range(100, waverange[0].value, waverange[1].value)
    sed.set_track_origin('basic')

    # Set number of photons
    m.set_n_photons(initial=1e5, imaging=1e6,
                raytracing_sources=1e4, raytracing_dust=1e6)

    # Set number of temperature iterations
    m.set_n_initial_iterations(5)
    m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

    # Write out file
    #import pdb; pdb.set_trace()
    m.write('sed4.rtin')
    m.run('sed4.rtout', mpi=False)
