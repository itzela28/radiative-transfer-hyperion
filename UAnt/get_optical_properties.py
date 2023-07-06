import numpy as np
from hyperion.dust import OpticalProperties
import h5py
from astropy import units

def get_opacity(file = 'Cstar_DustProperties.hdf5', wave_ref = 1 * units.micron):
	"""
	Return the opacity for a reference wavelength
		by reading information from an input HDF5
		dust properties file
	Arguments:
	file -- name (with extension) of HDF5 dust properties file
	wave_ref -- reference wavelength WITH UNITS
	"""
	f = h5py.File(file, 'r')

	d = OpticalProperties()
	d.from_hdf5_group(f)

	wave = (d.nu * units.Hz).to('micron', equivalencies = units.spectral())

	k = np.argmin(np.abs(wave - wave_ref))
	
	return d.chi[k] * units.cm**2 / units.g

#funcion para calcular masa, usando valor de k (kappa) y agregando array con diferentes opt depth (tau)
def get_env_mass(kap, tau, rmin, rmax):
	"""
	kap -- opacity in cm2/g, from get_opacity 
	tau -- optical depth
	rmin -- Rint of the shell convert to cm
	rmax -- Rext of the shell convert to cm
	"""
	m=4*np.pi*(tau/kap)*(rmin*rmax)
	return m

env_mass=[1.11757156e-05, 1.3659208e-05, 1.73844465e-05,
 		2.11096851e-05, 2.48349236e-05]*units.solMass
