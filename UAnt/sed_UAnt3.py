""""
SED for U Ant with different M_env mass
photometry obtained from get_SED.py
"""

from tkinter import font
import matplotlib.pyplot as plt
from get_SED import *
from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Distance to the source
Dpc = 268 * pc
# Open the model
m1 = ModelOutput('UAnt_tau0.2.rtout')
m2 = ModelOutput('UAnt_tau0.3.rtout')
m3 = ModelOutput('UAnt_tau0.4.rtout')
m4 = ModelOutput('UAnt_tau0.5.rtout')
m5 = ModelOutput('UAnt_tau0.6.rtout')
#m1_5 = ModelOutput('sed_UAnt_1_5.rtout')
#m2 = ModelOutput('sed_UAnt_2.rtout')
#m2_5 = ModelOutput('sed_UAnt_2_5.rtout')
#m3 = ModelOutput('sed_UAnt_3.rtout')
# Create the plot
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_subplot(1, 1, 1)

# Extract the SED for the smallest inclination and largest aperture, and
# scale to 300pc. In Python, negative indices can be used for lists and
# arrays, and indicate the position from the end. So to get the SED in the
# largest aperture, we set aperture=-1.
#total sed 0 deg
sed = get_SED(IRASPSC = ['10329-3918'])
sed_1 = m1.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#import pdb; pdb.set_trace()
sed_2 = m2.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
sed_3 = m3.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
sed_4 = m4.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
sed_5 = m5.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#total sed 90 deg
#sed_full90 = m.get_sed(inclination=9, aperture=-1, units='Jy', distance=Dpc)

# Plot the SED. The loglog command is similar to plot, but automatically
# sets the x and y axes to be on a log scale.
_, lam, fnu = get_stellar_params()
ax.loglog(lam,fnu, linestyle='--', label='Fotósfera')
ax.loglog(sed_1.wav, sed_1.val, color='green',linestyle='-', label=r'$\tau=0.2\times10^{-3}, M_{\rm{env}}=9.62 \times10^{-6}$ M$_{\odot}$')
#ax.loglog(sed_1_12.wav, sed_1_12.val, 'ro', label=r'$M_{env}=1 \times10^{-5}$ M$_{\odot}$')
#print("Number of wavelengths: {}".format(len(sed_1_12.wav)))
ax.loglog(sed_2.wav, sed_2.val, color='teal',linestyle='-', label=r'$\tau=0.3\times10^{-3}, M_{\rm{env}}=1.44 \times10^{-5}$ M$_{\odot}$')
ax.loglog(sed_3.wav, sed_3.val, color='purple',linestyle='-', label=r'$\tau=0.4\times10^{-3}, M_{\rm{env}}=1.92 \times10^{-5}$ M$_{\odot}$')
ax.loglog(sed_4.wav, sed_4.val, color='deeppink',linestyle='-', label=r'$\tau=0.5\times10^{-3}, M_{\rm{env}}=2.40 \times10^{-5}$ M$_{\odot}$')
ax.loglog(sed_5.wav, sed_5.val, color='tab:cyan',linestyle='-', label=r'$\tau=0.6\times10^{-3}, M_{\rm{env}}=2.89 \times10^{-5}$ M$_{\odot}$')
#ax.plot(sed.meta['lpivot'].value, sed['flux'][0], 'ro')
ax.errorbar(sed.meta['lpivot'].value, sed['flux'][0], yerr = sed['flux_err'][0], fmt ='bo', ecolor='gray',capsize=2)


plt.legend(loc=3, fontsize = 18)
plt.title(r'U Ant', fontsize = 25)
# Add some axis labels (we are using LaTeX here)
#ax.set_xticks(fontsize=20)
ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize = 25)
ax.set_ylabel(r'$F_\lambda$ [Jy]', fontsize = 25)


# Set view limits
ax.set_xlim(0.1, 3000.)
ax.set_ylim(1e-6,5e4)

ax.tick_params(axis='both', which='major', labelsize=20)
#ax.xticks(fontsize=20)
#ax.rc(axis='ytick', labelsize=20)

# Write out the plot
plt.show()
#fig.savefig('UAnt_model.pdf')

#.zshrc