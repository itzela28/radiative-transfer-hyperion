#extract and plot SEDs from output files .rout

import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Open the model
mO = ModelOutput('sed_rho-2_Odust.rtout')
mC = ModelOutput('sed_rho-2_Cdust.rtout')

fig = plt.figure(1, figsize=(12,8))
plt.subplots_adjust(wspace=0.05, bottom=0.105, 
                    right=0.985, left=0.1, 
                    top=0.94)

sed_O = mO.get_sed(inclination=0, aperture=-1, units='Jy', distance=260 * pc)
sed_C = mC.get_sed(inclination=0, aperture=-1, units='Jy', distance=260 * pc)
plt.plot(sed_O.wav, sed_O.val, color='purple',linestyle='--', linewidth=2.2,label='O-rich')
plt.plot(sed_C.wav, sed_C.val, color='teal',linestyle='-.', linewidth=2.2,label='C-rich')

plt.xscale('log')
plt.yscale('log')
plt.legend(loc=4, fontsize = 20)
plt.title(r'$\rho\propto r^{-2}$, $M_{\rm{env}}=2\times10^{-4} \rm{M}_{\odot}$', fontsize = 25)

plt.xlabel(r'$\lambda$ [$\mu$m]', fontsize = 25, fontproperties='serif')
plt.ylabel(r'$F_\lambda$ [Jy]', fontsize = 25, fontproperties='serif')
plt.tick_params(axis='both', which='major', labelsize=20)

# Set view limits
plt.xlim(0.4, 25.)
plt.ylim(1e-4,5e4)

# Write out the plot
plt.show()
fig.savefig('sed_dif_dust.pdf')