from turtle import width
import numpy as np
import matplotlib.pyplot as plt
import constCGS
from TDE_parameters import *

##### read data from file

Ldnuarr = np.loadtxt('specific_luminosity.txt')

##### discretization 

# observation frequency
lambobs = np.linspace(nuobsmin, nuobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobs/constCGS.cm2um) # in Hz

# observation time
tobsmin = robsmin*constCGS.pc2cm/constCGS.C_LIGHT
tobsmax = robsmax*constCGS.pc2cm/constCGS.C_LIGHT
tobsarr = np.linspace(tobsmin, tobsmax, Ntobs)
tobsarryr = tobsarr/constCGS.d2sec/365.25

##### plot
fig, axs = plt.subplots(figsize=(10, 8))

for i_nuobs in range(Nnuobs):
    axs.loglog(tobsarryr, Ldnuarr[i_nuobs, :], label='$\\nu_{obs}$'+' = {:.3e} Hz'.format(nuobsarr[i_nuobs]))
axs.set_xlabel('$t_{obs}$ (yr)', size=15)
axs.set_ylabel('$L_\\nu\ (erg^{-1} Hz^{-1})$', size=15)
axs.legend()
axs.grid()
axs.tick_params(axis='both', which='major', direction='in', length=10, width=2)
axs.tick_params(axis='both', which='minor', direction='in', length=8, width=1)

plt.show()