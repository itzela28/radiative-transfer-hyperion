# radiative-transfer-hyperion
To generate the radiative transfer model:
  1. Modify the code radiative_transfer.py.
     Input parameters:
      a) reference wavelength ref_wave
      b) array for optical depth opt_depth=[ ]
      c) number of the source i (in table2_1() of the file get_Aringeretal2009_model.py, the data of each source is found).
      d) j, position of the optical depth value of the opt_depth array
      e) dust_file_name, name of the file with the dust properties
      f)m.write('file_name.rtin')
      g)m.run('file_name.rtout')
     This file uses the functions to extrapolate the spectrum and obtain the optical constants, as well as the function to extrapolate the optical constants. 
     It gets the SED of each source and its photosphere model from the get_SED.py file, which uses the file Aringer+2009_StellarSpectrum.ascii. radiative_transfer.py generates the models with extension .rtin and .rtout.
  
  2. Use the plot_SED.py file to extract the SED from the models. This file uses the model with the .rtout extension to extract the SED and create the plot. Modify:
    a) i, source number. This is used to choose the position of the source in the id, iras and d arrays.
    b) m1, m2, etc. Write the name of the .rtout file
    c) change the chart limits in set_xlim, set_ylim if necessary.
    d) rename the .pdf file to create the plot in fig.savefig()
