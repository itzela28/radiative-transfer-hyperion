To generate the radiative transfer model:

1. Modify the code UAnt.py. Input parameters:
a) m, lam, fnu, uses get_stellar_params(m) function from get_SED.py (This file gets the SED of each source and its photosphere model). Inside this function, the stellar parameters for U Ant are specified. It takes the spectrum file Aringer+2009_StellarSpectrum.ascii. This file is for U Ant stellar parameters. This function returns the values for the stellar spectrum.
b) Envelope parameters:
env.rmin, env.rmax (both in parsec), env.mass (Takes a velue for the array env_mass, in solar mass. These values are computed in a separate code, using values for the optical depth, opacity, R_int and R_out of the shell). Code for the envelope mass: https://colab.research.google.com/drive/1t6oIEgHLu57XQxtfjKulCoQCb3sFzDVd?usp=sharing
c) dust_file_name, name of the file with the dust properties
d)m.write('file_name.rtin')
e)m.run('file_name.rtout'), generates the models with extension .rtin and .rtout.

2. Use the sed_UAnt3.py file to extract the SED from the models. This file uses the model with the .rtout extension to extract the SED and create the plot. Modify:
a) Dpc = distance to the star, in parsec.
b) sed= get_SED(IRASPSC=['iras_id_of_the_source']).
c) m1, m2, etc. Write the name of the .rtout file
d) change the plot limits in set_xlim, set_ylim if necessary.
e) rename the .pdf file to create the plot in fig.savefig()
