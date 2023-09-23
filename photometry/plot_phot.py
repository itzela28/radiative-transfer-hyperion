"""""
SEDs with photometric data from NESS
Model spectrum obtained from GRAMS (get_Aringeretal2009_model)
"""""
import numpy as np
import astropy
from astropy import units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt
from astropy.table import Table, join
from get_SED import *
from get_Aringeretal2009_model import *

# data = Table.read('sed7.vot', format = 'votable')
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
spec = get_spec(t,Dkpc=[0.532,0.266,0.275,0.268,0.208,0.528,0.321],C2O=[1.30,1.34,1.07,1.44,1.043,1.10,1.087])

all_data = join(data, spec, keys = 'IRASPSC', join_type = 'inner')
#import pdb; pdb.set_trace()
#print(data['flux_err'])

i=1 #numero de la fuente
x=data.meta['lpivot'].value
y=data['flux'][i]
y_err=data['flux_err'][i]
umbral_error=200

filters=np.nonzero(y_err < umbral_error)[0]
x_fil=x[filters]
y_fil=y[filters]
err=y_err[filters]

fig = plt.figure(1, figsize=(10,8))
plt.subplots_adjust(wspace=0.05, bottom=0.1, 
                    right=0.985, left=0.12, 
                    top=0.95)
plt.errorbar(x_fil, y_fil, 
            yerr = err, fmt='o', 
            color='navy', markersize=5, ecolor='gray',
            elinewidth=1,capsize=2)
plt.plot(all_data['wspec'][i], all_data['fspec'][i], 
        color = 'gray', linestyle='--', label='Modelo espectral')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0, 5e3)
plt.legend(loc=1, fontsize = 20)
plt.xlabel(r'$\lambda$ [$\mu$m]', fontsize = 25, fontproperties='serif')
plt.ylabel(r'$F_\lambda$ [Jy]', fontsize = 25, fontproperties='serif')
plt.tick_params(axis='both', which='major', labelsize=20)
#plt.title(r'R Scl', fontproperties='serif', fontsize = 25)
plt.show()
#fig.savefig('RTCap.pdf')