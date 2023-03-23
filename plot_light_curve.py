import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10, cos, pi
import constCGS
import os
from TDE_parameters import *
from TDE_observation_setup import *
# from PTF09ge_parameters import *
# from PTF09ge_observation_setup import *
# from ASASSN14li_parameters import *
# from ASASSN14li_observation_setup import *

##### read data from file
with open(os.path.join(folder, 'specific_luminosity_tobsmax{:.2e}.txt'.format(tobsmax)), 'r') as file:
    Ldnuarr_shape = file.readline().strip('\n').split()
#     Lbol = float(file.readline().strip('\n').split()[1])
    Ldnuarr = np.loadtxt(file)

##### discretization 

# observation frequency
lambobsarr = np.linspace(lambobsmin, lambobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobsarr/constCGS.cm2um) # in Hz

# observation time
# tobsmax = robsmax*constCGS.pc2cm/constCGS.C_LIGHT
# tobsmin = robsmin*constCGS.pc2cm/constCGS.C_LIGHT*(1-cos(max(0, 4*pi/180)))
# tobsmax = robsmax*constCGS.pc2cm/constCGS.C_LIGHT
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
axs.set_ylabel('$L_\\nu\ (erg s^{-1} Hz^{-1})$', size=15)
axs.set_xlim([1, 1e5])
axs.set_ylim(ymin=1e25)
axs.legend(fontsize=10)
axs.grid()
axs.tick_params(axis='both', which='major', direction='in', length=12, width=1.5)
axs.tick_params(axis='both', which='minor', direction='in', length=8, width=1)
axs.annotate('params: Ttde{:.2e}K_tdur{:.2e}s_nH0{:.2e}_nHpower{:.2f}_lamb0{:d}um'.format(T_tde, tdur, nH0, densprof, lamb0),\
        xy=(0.5, 0.01), xycoords='figure fraction',\
        xytext=(0.5, 0), textcoords= 'figure fraction', \
        horizontalalignment='center', verticalalignment='bottom', size=13)
# axs.annotate("$L_{IR}$"+ " = {:.3e} erg/s".format(Lbol),\
#         xy=(0.3, 0.2), xycoords='axes fraction',\
#         xytext=(0.3, 0.2), textcoords= 'axes fraction', \
#         horizontalalignment='center', verticalalignment='bottom', size=12)
axs.set_title('Dust Echo Light Curve', size=15)

plt.savefig(os.path.join(folder, 'light_curve.jpg'))
plt.show()