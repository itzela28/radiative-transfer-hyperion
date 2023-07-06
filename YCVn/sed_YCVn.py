""""
SED for U RScl with different M_env mass
photometry obtained from get_SED.py
"""

from tkinter import font
import matplotlib.pyplot as plt
from astropy import units as u
from get_SED import *
from get_Aringeretal2009_model import *
from hyperion.model import ModelOutput
from hyperion.util.constants import pc


# Distance to the source
Dpc = 321*(np.sqrt(2)) * pc
# Open the model
m1 = ModelOutput('YCVn_m1.rtout')
m2 = ModelOutput('YCVn_m2.rtout')
#m3 = ModelOutput('UCam_m3.rtout')
#m4 = ModelOutput('UCam_m4.rtout')
#m5 = ModelOutput('UCam_m5.rtout')

data = get_SED(IRASPSC = ['03374+6229',
'01246-3248',
'23438+0312',
'10329-3918',
'10350-1307',
'20141-2128',
'12427+4542'])
data['flux'].unit = u.Jy
data['flux_err'].unit = u.Jy

t = table2_1()
spec = get_spec(t,Dkpc=[0.532,0.266,0.275,0.268,0.208,0.291,0.321],C2O=[1.30,1.34,1.07,1.44,1.043,1.10,1.087])

all_data = join(data, spec, keys = 'IRASPSC', join_type = 'inner')


# Create the plot
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_subplot(1, 1, 1)

# Extract the SED for the smallest inclination and largest aperture, and
# scale to 300pc. In Python, negative indices can be used for lists and
# arrays, and indicate the position from the end. So to get the SED in the
# largest aperture, we set aperture=-1.
#total sed 0 deg
sed = get_SED(IRASPSC = ['12427+4542'])
sed_1 = m1.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#import pdb; pdb.set_trace()
sed_2 = m2.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#sed_3 = m3.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#sed_4 = m4.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#sed_5 = m5.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#total sed 90 deg
#sed_full90 = m.get_sed(inclination=9, aperture=-1, units='Jy', distance=Dpc)

# Plot the SED. The loglog command is similar to plot, but automatically
# sets the x and y axes to be on a log scale.
#_, lam, fnu = get_stellar_params()
#ax.loglog(lam,fnu, linestyle='--', label='Fotósfera')

ax.loglog(all_data['wspec'][6], all_data['fspec'][6], 
        color = 'gray', linestyle='--', label='Modelo espectral')
ax.errorbar(data.meta['lpivot'].value, data['flux'][6], 
            yerr = data['flux_err'][6], fmt='o', 
            color='navy', markersize=5, ecolor='gray',
            elinewidth=1,capsize=2)

ax.loglog(sed_1.wav, sed_1.val, color='green',linestyle='-', label=r'$\tau=0.3\times10^{-5}, M_{\rm{env}}=1.17 \times10^{-5}$ M$_{\odot}$')
#print("Number of wavelengths: {}".format(len(sed_1_12.wav)))
ax.loglog(sed_2.wav, sed_2.val, color='teal',linestyle='-', label=r'$\tau=0.5\times10^{-3}, M_{\rm{env}}=9.47 \times10^{-6}$ M$_{\odot}$')
#ax.loglog(sed_3.wav, sed_3.val, color='purple',linestyle='-', label=r'$\tau=0.7\times10^{-3}, M_{\rm{env}}=1.32 \times10^{-5}$ M$_{\odot}$')
#ax.loglog(sed_4.wav, sed_4.val, color='deeppink',linestyle='-', label=r'$\tau=0.9\times10^{-3}, M_{\rm{env}}=1.70 \times10^{-5}$ M$_{\odot}$')
#ax.loglog(sed_5.wav, sed_5.val, color='tab:cyan',linestyle='-', label=r'$\tau=1.1\times10^{-3}, M_{\rm{env}}=2.08 \times10^{-5}$ M$_{\odot}$')
#ax.plot(sed.meta['lpivot'].value, sed['flux'][0], 'ro')
#ax.errorbar(sed.meta['lpivot'].value, sed['flux'][0], yerr = sed['flux_err'][0], fmt ='bo', ecolor='gray',capsize=2)


plt.legend(loc=3, fontsize = 18)
plt.title(r'Y CVn', fontsize = 25)
# Add some axis labels (we are using LaTeX here)
#ax.set_xticks(fontsize=20)
ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize = 25)
ax.set_ylabel(r'$F_\lambda$ [Jy]', fontsize = 25)


# Set view limits
ax.set_xlim(0.2, 2.5e3)
ax.set_ylim(5e-4, 2e4)

ax.tick_params(axis='both', which='major', labelsize=20)
#ax.xticks(fontsize=20)
#ax.rc(axis='ytick', labelsize=20)

# Write out the plot
plt.show()
#fig.savefig('YCVn_model.pdf', overwrite=True)

#.zshrc
