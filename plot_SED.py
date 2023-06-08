"""
plot SED for diferent sources using file.rtout
photometry obtained from get_SED.py
"""
from tkinter import font
import matplotlib.pyplot as plt
from get_SED import *
from hyperion.model import ModelOutput
from hyperion.util.constants import pc

i=0
id = ['U Cam', 'R Scl', 'TX Psc', 'U Ant', 'U Hya', 'RT Cap', 'Y CVn']#, 'Y CVn']
iras = ['03374+6229' , '01246-3248', '23438+0312', '10329-3918', \
               '10350-1307', '20141-2128', '12427+4542']
d = [532, 266, 275, 268, 208, 291, 321]
# Distance to the source
Dpc = d[i] * pc
# Open the model
m1 = ModelOutput('sed_UAnt_m1.rtout')
m2 = ModelOutput('sed_UAnt_m2.rtout')
m3 = ModelOutput('sed_UAnt_m3.rtout')
m4 = ModelOutput('sed_UAnt_m4.rtout')
m5 = ModelOutput('sed_UCam_tau0.9.rtout')
# Create the plot
fig = plt.figure(1, figsize=(12,8))
ax = fig.add_subplot(1, 1, 1)

# Extract the SED for the smallest inclination and largest aperture, and
# scale to 300pc. In Python, negative indices can be used for lists and
# arrays, and indicate the position from the end. So to get the SED in the
# largest aperture, we set aperture=-1.
#total sed 0 deg
sed = get_SED(IRASPSC = [iras[i]]) #03374+6229=UCam
sed1 = m1.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#import pdb; pdb.set_trace()
sed2 = m2.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
sed3 = m3.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
sed4 = m4.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
sed5 = m5.get_sed(inclination=0, aperture=-1, units='Jy', distance=Dpc)
#total sed 90 deg
#sed_full90 = m.get_sed(inclination=9, aperture=-1, units='Jy', distance=Dpc)

# Plot the SED. The loglog command is similar to plot, but automatically
# sets the x and y axes to be on a log scale.
_, lam, fnu = get_stellar_params()
ax.loglog(lam,fnu, linestyle='--', label='Fotósfera')
#ax.loglog(sed1.wav, sed1.val, color='green',linestyle='-', label=r'$M_{env}=1.12 \times10^{-5}$ M$_{\odot}$')
#ax.loglog(sed5.wav, sed5.val, 'ro', label=r'$M_{env}=1.12 \times10^{-5}$ M$_{\odot}$')
ax.loglog(sed5.wav, sed5.val, 'r-', label=r'$\tau=0.9\times10^{-4}$')
#print("Longest wavelength : {}".format(sed5.wav[0]))
#print("Flux at longest wavelength : {}".format(sed5.val[0]))
#ax.loglog(sed2.wav, sed2.val, color='teal',linestyle='-', label=r'$M_{env}=1.36 \times10^{-5}$ M$_{\odot}$')
#ax.loglog(sed2.wav, sed2.val, color='teal',linestyle='-', label=r'test 2')
#ax.loglog(sed3.wav, sed3.val, color='purple',linestyle='-', label=r'$M_{env}=1.74 \times10^{-5}$ M$_{\odot}$')
#ax.loglog(sed4.wav, sed4.val, color='deeppink',linestyle='-', label=r'$M_{env}=2.11 \times10^{-5}$ M$_{\odot}$')
#ax.plot(sed.meta['lpivot'].value, sed['flux'][0], 'ro')
ax.errorbar(sed.meta['lpivot'].value, sed['flux'][0], yerr = sed['flux_err'][0], fmt ='bo', ecolor='gray',capsize=2)
#ax.loglog(sed5.wav, sed5.val, color='tab:cyan',linestyle='-', label=r'$M_{env}=2.48 \times10^{-5}$ M$_{\odot}$')

plt.legend(loc=3, fontsize = 18)
plt.title(r'{}'.format(id[i]), fontsize = 25)
# Add some axis labels (we are using LaTeX here)
#ax.set_xticks(fontsize=20)
ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize = 25)
ax.set_ylabel(r'$F_\lambda$ [Jy]', fontsize = 25)


# Set view limits
ax.set_xlim(0.1, 5000.)
ax.set_ylim(1e-6,10e4)

ax.tick_params(axis='both', which='major', labelsize=20)
#ax.xticks(fontsize=20)
#ax.rc(axis='ytick', labelsize=20)

# Write out the plot
plt.show()
#fig.savefig('UAnt_range1.pdf')

#.zshrc
