import numpy as np
import os
from math import pi, log10, cos, sin, acos
import constCGS
import matplotlib.pyplot as plt
from matplotlib import cm

##### import parameters

from TDE_parameters import *
from TDE_observation_setup import *

##### read data from file

# arrays
with open(os.path.join(folder, 'arrays_info.txt'), 'r') as file:
    file.readline()
    tarr_info = file.readline().strip('\n').split()
    rarr_info = file.readline().strip('\n').split()
    aarr_info = file.readline().strip('\n').split()

# emissvity 
with open(os.path.join(folder, 'emissivity.txt'), 'r') as file:
    jdnuarr_shape = file.readline().strip('\n').split()
    jdnuarr = np.loadtxt(file)
jdnuarr_shape = [int(shape) for shape in jdnuarr_shape]
jdnuarr = np.reshape(jdnuarr, jdnuarr_shape)

# tarr 
tst, tend, Nt = float(tarr_info[2]), float(tarr_info[3]), int(tarr_info[4])
if tarr_info[1] == 'linear':
    tarr = np.linspace(tst, tend, Nt)
elif tarr_info[1] == 'logarithmic':
    tarr = np.logspace(log10(tst), log10(tend), Nt)

# rarr
rst, rend, Nr = float(rarr_info[2]), float(rarr_info[3]), int(rarr_info[4])
if rarr_info[1] == 'linear':
    rarr = np.linspace(rst, rend, Nr)
elif rarr_info[1] == 'logarithmic':
    rarr = np.logspace(log10(rst), log10(rend), Nr)

# observation frequency
lambobsarr = np.linspace(lambobsmin, lambobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobsarr/constCGS.cm2um) # in Hz

##### plot
colors = cm.viridis(np.linspace(0, 1, Nnuobs))[::-1]
fig, axs = plt.subplots(2, 1, figsize=(7, 8))

i_r = 5
for i_nuobs in range(0, Nnuobs, 2):
#     axs.semilogy(tobsarrd, Ldnuarr[i_nuobs, :], linewidth=2, color = colors[i_nuobs],\
#         label='$\lambda_{obs}$'+' = {:.2f} '.format(lambobsarr[i_nuobs])+'$\mu$m')
    axs[0].semilogy(tarr, jdnuarr[:, i_r, i_nuobs], linewidth=2, color = colors[i_nuobs],\
        label='$\lambda_{obs}$'+' = {:.2f} '.format(lambobsarr[i_nuobs])+'$\mu$m')

axs[0].set_xlabel('t (s)', size=10)
axs[0].set_ylabel('$j_\\nu\ (erg s^{-1} cm^{-3} sr^{-1} Hz^{-1})$', size=10)
# axs.set_xlim([1, 1e5])
# axs.set_ylim(ymin=1e25)
axs[0].legend(fontsize=10)
axs[0].grid()
axs[0].tick_params(axis='both', which='major', direction='in', length=12, width=1.5)
axs[0].tick_params(axis='both', which='minor', direction='in', length=8, width=1)
# axs[0].annotate(' Ttde{:.2e}K_tdur{:.2e}s_nH0{:.2e}_nHpower{:.2f}_lamb0{:d}um'.format(T_tde, tdur, nH0, densprof, lamb0),\
#         xy=(0.5, 0.01), xycoords='figure fraction',\
#         xytext=(0.5, 0), textcoords= 'figure fraction', \
#         horizontalalignment='center', verticalalignment='bottom', size=13)
# axs[0].annotate("$L_{IR}$"+ " = {:.3e} erg/s".format(Lbol),\
#         xy=(0.3, 0.2), xycoords='axes fraction',\
#         xytext=(0.3, 0.2), textcoords= 'axes fraction', \
#         horizontalalignment='center', verticalalignment='bottom', size=12)
# axs[0].set_title('Dust Echo Light Curve', size=15)

i_t = 15
for i_nuobs in range(0, Nnuobs, 2):
#     axs.semilogy(tobsarrd, Ldnuarr[i_nuobs, :], linewidth=2, color = colors[i_nuobs],\
#         label='$\lambda_{obs}$'+' = {:.2f} '.format(lambobsarr[i_nuobs])+'$\mu$m')
    axs[1].semilogy(rarr, jdnuarr[i_t, :, i_nuobs], linewidth=2, color = colors[i_nuobs],\
        label='$\lambda_{obs}$'+' = {:.2f} '.format(lambobsarr[i_nuobs])+'$\mu$m')
axs[1].set_xlabel('r (pc)', size=10)
axs[1].set_ylabel('$j_\\nu\ (erg s^{-1} cm^{-3} sr^{-1} Hz^{-1})$', size=10)
# axs.set_xlim([1, 1e5])
# axs.set_ylim(ymin=1e25)
axs[1].legend(fontsize=10)
axs[1].grid()
axs[1].tick_params(axis='both', which='major', direction='in', length=12, width=1.5)
axs[1].tick_params(axis='both', which='minor', direction='in', length=8, width=1)
axs[1].annotate('Ttde{:.2e}K_tdur{:.2e}s_nH0{:.2e}_nHpower{:.2f}_lamb0{:d}um'.format(T_tde, tdur, nH0, densprof, lamb0),\
        xy=(0.5, 0.01), xycoords='figure fraction',\
        xytext=(0.5, 0), textcoords= 'figure fraction', \
        horizontalalignment='center', verticalalignment='bottom', size=10)
# axs[1].annotate("$L_{IR}$"+ " = {:.3e} erg/s".format(Lbol),\
#         xy=(0.3, 0.2), xycoords='axes fraction',\
#         xytext=(0.3, 0.2), textcoords= 'axes fraction', \
#         horizontalalignment='center', verticalalignment='bottom', size=12)
# axs[1].set_title('Dust Echo Light Curve', size=15)
fig.suptitle('Emissivity', fontsize=15)

plt.savefig(os.path.join(folder, 'emissivity.jpg'))
plt.show()