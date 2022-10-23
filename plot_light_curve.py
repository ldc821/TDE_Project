import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10, cos, pi
import constCGS
import os
from TDE_parameters import *
from observation_setup import *

##### read data from file

Ldnuarr = np.loadtxt(os.path.join(folder,'specific_luminosity.txt'), skiprows=1)

##### discretization 

# observation frequency
lambobsarr = np.linspace(lambobsmin, lambobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobsarr/constCGS.cm2um) # in Hz

# observation time
tobsmin = robsmin*constCGS.pc2cm/constCGS.C_LIGHT
tobsmax = robsmax*constCGS.pc2cm/constCGS.C_LIGHT
tobsarr = np.logspace(log10(tobsmin), log10(tobsmax), Ntobs)
tobsarrd = tobsarr/constCGS.d2sec
tobsarryr = tobsarrd/constCGS.yr2d

##### plot
colors = cm.viridis(np.linspace(0, 1, Nnuobs))[::-1]
fig, axs = plt.subplots(figsize=(8, 6))

for i_nuobs in range(0, Nnuobs, 2):
    axs.loglog(tobsarrd, Ldnuarr[i_nuobs, :], linewidth=2, color = colors[i_nuobs],\
        label='$\lambda_{obs}$'+' = {:.2f} '.format(lambobsarr[i_nuobs])+'$\mu$m')
axs.set_xlabel('$t_{obs}$ (d)', size=15)
axs.set_ylabel('$L_\\nu\ (erg^{-1} Hz^{-1})$', size=15)
axs.legend(fontsize=10)
axs.grid()
axs.tick_params(axis='both', which='major', direction='in', length=12, width=1.5)
axs.tick_params(axis='both', which='minor', direction='in', length=8, width=1)

plt.show()