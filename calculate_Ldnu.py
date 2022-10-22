from matplotlib.pyplot import contour
import numpy as np
import os
from math import pi, log10, cos
import constCGS

##### import parameters

from TDE_parameters import *
from observation_setup import *

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

##### discretization

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

# aarr
ast, aend, Na = float(aarr_info[2]), float(aarr_info[3]), int(aarr_info[4])
if aarr_info[1] == 'linear':
    aarr = np.linspace(ast, aend, Na)
elif aarr_info[1] == 'logarithmic':
    aarr = np.logspace(log10(ast), log10(aend), Na)

# observation frequency
lambobs = np.linspace(lambobsmin, lambobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobs/constCGS.cm2um) # in Hz

# observation time
tobsmin = robsmin*constCGS.pc2cm/constCGS.C_LIGHT
tobsmax = robsmax*constCGS.pc2cm/constCGS.C_LIGHT
tobsarr = np.logspace(log10(tobsmin), log10(tobsmax), Ntobs)

# specific luminosity
Ldnuarr = np.zeros((Nnuobs, Ntobs), dtype=float)

##### functions

# linear interpolation of the emissivity at a given distance
def jdnu_intp(t, i_r, i_nuobs):
    i_floor = np.argmin(np.abs(t - tarr))
    if t > tarr[-1] + (tarr[-1] - tarr[-2]):
        return 0
    if tarr[i_floor] == t:
        return jdnuarr[i_floor, i_r, i_nuobs]
    if (tarr[i_floor] > t and i_floor > 0) or i_floor == Nt-1:
        i_floor -= 1
    slope = (jdnuarr[i_floor+1, i_r, i_nuobs] - jdnuarr[i_floor, i_r, i_nuobs])/(tarr[i_floor+1] - tarr[i_floor])
    return max(jdnuarr[i_floor, i_r, i_nuobs] + slope * (t - tarr[i_floor]), 0)

# calculate mu from tobs and t, Eq. 27
def mu(rcm, tobs, i_t):
    return 1 - constCGS.C_LIGHT*(tobs - tarr[i_t])/rcm

# solve tobs at a given r from mu
def musolve_t(mu, rcm, tobs):
    return tobs - rcm/constCGS.C_LIGHT*(1 - mu)

# specifc luminosity, Eq.28
def lum_dnu(tobs, i_nuobs):
    r_integral = 0
    for i_r in range(Nr-1):
        rcm = rarr[i_r]*constCGS.pc2cm
        drcm = (rarr[i_r+1] - rarr[i_r])*constCGS.pc2cm
        # print('min: {:.3f}pc, r: {:.3f}pc'.format(constCGS.C_LIGHT*(tobs-tmin)/constCGS.pc2cm, rarrlog[i_r]))
        # light echo has already passed
        if rcm < max(0, 2*constCGS.C_LIGHT*(tobs-tmax)):
            continue 
        # mumin= max(1 - constCGS.C_LIGHT*tobs/rcm, -1)
        # 1e-10 is added to account for calculation instability that may cause t < 0
        mumin = 1 - constCGS.C_LIGHT*tobs/rcm + 1e-10
        # mumax = min(1 - constCGS.C_LIGHT*(tobs - tdur)/rcm, 1)
        mumax = 1 - constCGS.C_LIGHT*(tobs - tdur)/rcm
        mu_integral = 0
        muarr = np.linspace(mumin, mumax, Nmu)
        for i_mu in range(Nmu - 1):
            dmu = muarr[i_mu+1] - muarr[i_mu]
            t = musolve_t(muarr[i_mu], rcm, tobs)
            mu_integral += dmu * jdnu_intp(t, i_r, i_nuobs)
        r_integral += mu_integral * rcm**2 * drcm
    return 8 * pi**2 * r_integral

##### calculate specific luminosity

print('Calculating...')
progress = 0

for i_tobs in range(Ntobs):
    tobs = tobsarr[i_tobs]
    for i_nuobs in range(Nnuobs):
        Ldnuarr[i_nuobs, i_tobs] = lum_dnu(tobs, i_nuobs)
    if i_tobs/Ntobs > progress:
        print('{:.2%}'.format(i_tobs/Ntobs), end='...')
        progress += 0.1

##### save data

print('\n\nSaving data...')

Ldnuarr_shape = '{}\t{}\n'.format(Nnuobs, Ntobs)
with open(os.path.join(folder, 'specific_luminosity.txt'), 'w') as file:
    file.write(Ldnuarr_shape)
    np.savetxt(file, Ldnuarr)

print('All done!')