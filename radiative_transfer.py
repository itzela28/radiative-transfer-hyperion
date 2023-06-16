"""
"""
import numpy as np
from hyperion.util.constants import rsun, au, msun, sigma, pc, lsun
from astropy import units 
import matplotlib.pyplot as plt
# from get_SED import get_stellar_params
from get_SED import *
from get_Aringeretal2009_model import *
import h5py, h5glance
from h5glance import H5Glance
from hyperion.dust import SphericalDust

ref_wave=11.3
opt_depth=[0.9e-4, 1.1e-4, 1.4e-4, 1.7e-4, 2.0e-4]
d = [532, 266, 275, 268, 208, 291, 321]
i = 0 #number of the source
j = 0 #position of the value in opt_depth
dust_file_name="Cstar_DustProperties.hdf5"

#dust_file2=h5py.File("AmC_Rouleau1991_plus_SiC_Pegourie1988.hdf5", "r+")


def opacity(ref_wave, dust_file_name):
    f=h5py.File(dust_file_name, "r")
    albedo=f['optical_properties']['albedo']
    chi=f['optical_properties']['chi']
    k_sca=chi*albedo
    k_abs=chi*(1-albedo)
    frec=f['optical_properties']['nu']*units.Hz
    wave=frec.to('um',equivalencies=units.spectral())
    r_wave = ref_wave * units.um
    cls_wave=wave[np.argmin(np.abs(wave - r_wave))]
    chi_cw=chi[np.argmin(np.abs(wave - r_wave))]
    return cls_wave, chi_cw

def extrap_opt_c(in_dust_file, min_wave, max_wave, out_dust_file):
    #min_wave, max_wave in micron
    d=SphericalDust()
    d.read(in_dust_file)
    d.optical_properties.extrapolate_wav(min_wave, max_wave)
    d.write(out_dust_file)
    return




def extrap(wave, fspec, trace, output='photosphere_extrap.vot'):
    """
    Extrapolate the stellar spectrum as follows:
    1. Compute the stellar continuum using fspec and trace
    2. Choose only the fluxes for wavelengths beyond 1.7 µm for mid-IR and before 0.5 µm for UV
    3. Fit this spectrum with a polynomial. IMPORTANT: fit log(flux) to log(wavelength)
    4. Use the fit coefficients to predict the fluxes for wavelengths in the range (24, 1000) and (0.2, 0.4447)
    5. This is the new input spectrum (add these fluxes to the existing grid)

    wave=wavelenght
    fspec=flux (stellar spectrum)
    """
    fcont = fspec/trace #flux continumm
    wave_log=np.log10(wave)
    fcont_log=np.log10(fcont)

    #fit & prediction for UV 
    w_filter1=wave_log[wave<=0.5] #filter for UV range
    fcont_filter1=fcont_log[wave<=0.5]
    coeff1=np.polyfit(w_filter1, fcont_filter1, 2) 
    polynom1=np.poly1d(coeff1)
    fcont_fit1=polynom1(w_filter1)
    f=h5py.File("Cstar_DustProperties.hdf5", "r+")
    frec_dust=f['optical_properties']['nu']*units.Hz
    wave_dust=frec_dust.to('um',equivalencies=units.spectral())
    w_min=wave_dust.min()
    w_max=wave_dust.max()
    w_pred1 = np.geomspace(w_min.value*1.05, 0.4447, 20)
    fcont_pred1= 10**(polynom1(np.log10(w_pred1)))
    trace1=np.ones(20)*trace[0]

    #fit & prediction for mid-IR
    w_filter2=wave_log[wave>=10] #filter for mid-IR range
    fcont_filter2=fcont_log[wave>=10]
    coeff2=np.polyfit(w_filter2, fcont_filter2, 1) 
    polynom2=np.poly1d(coeff2)
    fcont_fit2=polynom2(w_filter2)
    w_pred2 = np.geomspace(24, w_max.value*0.95, 20)
    fcont_pred2 = 10**(polynom2(np.log10(w_pred2)))
    trace2=np.ones(20)

    #add the prediction to the data
    w_ex=np.append(w_pred1, np.append(wave,w_pred2))
    fcont_ex=np.append(fcont_pred1, np.append(fcont, fcont_pred2))
    tr_ex=np.append(trace1, np.append(trace,trace2))
    fspec_ex=fcont_ex*tr_ex

    t = Table([w_ex, fspec_ex, tr_ex], names=('wave','spec','trace'))
    t.write(output, format='votable', overwrite=True)
    return w_ex, fspec_ex, tr_ex



class PowerLawEnvelope():
    from hyperion.model import AnalyticalYSOModel
    from hyperion.densities.power_law_envelope import PowerLawEnvelope

    # Initalize the model
    m = AnalyticalYSOModel()
    
    data = get_SED(IRASPSC = ['03374+6229',
    '01246-3248',
    '23438+0312',
    '10329-3918',
    '10350-1307',
    '20141-2128',
    '12427+4542'])
    data['flux'].unit = units.Jy
    data['flux_err'].unit = units.Jy

    spec = get_spec(table2_1(),
        Dkpc=[0.532,0.266,0.275,0.268,0.208,0.291,0.321],C2O=[1.30,1.34,1.07,1.44,1.043,1.10,1.087])

    fspec1 = spec['fspec'][i]
    trace1 = spec['trace'][i]
    wave1 = spec['wspec'][i]

    wave_ex, fspec_ex, trace_ex = extrap(wave1, fspec1, trace1)

    """
    wave2, chi1 = opacity(ref_wave, dust_file_name)
    plt.plot(wave2, chi1, 'r-')
    plt.plot(wave_ex, fspec_ex, 'b-')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    import pdb; pdb.set_trace()
    """
    # the lam and fnu can be used to plot the stellar photosphere
    # lam is in micron, fnu is in J
    m.star.luminosity = spec['Lbol_Kastner'][i] * spec['Lbol_Kastner'].unit.to('erg/s')
    m.star.mass = spec['Mass'][i] * spec['Mass'].unit.to('g')
    logg = spec['logg'][i] # surface gravity in cm/s^2
    logg_sun = 4.44 # solar surface gravity in the same units
    m.star.radius = np.sqrt(m.star.mass / 10**(logg - logg_sun))
    wave = wave_ex
    #m.star.spectrum = (wave.to('Hz', equivalencies=units.spectral()).value, 
    #spec['fspec'][i])
    m.star.spectrum = (wave, fspec_ex)
    m, lam, fnu = get_stellar_params(m, dist=d[i]*units.pc)

    #extrapolate opt. constants onto the wave grid of photospheric spectra
    dust_extrap_file_name=dust_file_name[:-5] + '_extrap.hdf5'
    extrap_opt_c(dust_file_name, wave[0]*0.9, wave[-1]*1.05, dust_extrap_file_name)
    # Add an envelope
    close_wave, chi_cl = opacity(ref_wave, dust_extrap_file_name)    
    
    env = m.add_power_law_envelope()
             # Envelope mass
    env.rmin =  (spec['Rs'][i]-(spec['dR'][i]/2)) * spec['Rs'].unit.to('cm') # Inner radius, 34au agb detached shell
    env.rmax =  (spec['Rs'][i]+(spec['dR'][i]/2)) * spec['Rs'].unit.to('cm') # Outer radius 
    #env.mass in terms of opt. depth, rmax,rmin, opacity
    env.mass = (4*np.pi)*(opt_depth[j]/chi_cl)*(env.rmax*env.rmin) 
    env.power = -2 # Radial power
    #env.rho_0 = 2 #only needed if env.mass is not given
    env.r_0 = env.rmin #34 au for detached shell agb
    env.dust = dust_extrap_file_name
    #import pdb; pdb.set_trace()

    # Use raytracing to improve s/n of thermal/source emission
    m.set_raytracing(True)

    # Use the modified random walk
    m.set_mrw(True, gamma=2.)

    # Set up grid
    m.set_spherical_polar_grid_auto(100, 1, 1)

    # Set up SED for 10 viewing angles
    waverange = np.array([wave[0], wave[-1]]) * units.micron
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
    m.write('sed_UCam_tau0.9_test.rtin')
    m.run('sed_UCam_tau0.9_test.rtout', mpi=False)