import numpy as np
from math import pi, log10, exp
import constCGS
import os

##### import parameters

from TDE_parameters import *
from TDE_observation_setup import *
# from PTF09ge_parameters import *
# from PTF09ge_observation_setup import *
# from ASASSN14li_parameters import *
# from ASASSN14li_observation_setup import *

##### read data from file

# arrays
with open(os.path.join(folder, 'arrays_info.txt'), 'r') as file:
    file.readline()
    tarr_info = file.readline().strip('\n').split()
    rarr_info = file.readline().strip('\n').split()
    aarr_info = file.readline().strip('\n').split()

# sublimation radius
with open(os.path.join(folder, 'sublimation_radius.txt'), 'r') as file:
    asubarr = np.loadtxt(file, skiprows=1)

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

# dust temperature
with open(os.path.join(folder, 'dust_temperature.txt'), 'r') as file:
    Tarr_shape = file.readline().strip('\n').split()
    Tarr = np.loadtxt(file)
Tarr_shape = [int(shape) for shape in Tarr_shape]
Tarr = np.reshape(Tarr, Tarr_shape)

# observation frequency
lambobsarr = np.linspace(lambobsmin, lambobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobsarr/constCGS.cm2um) # in Hz

# emissivity
jdnuarr = np.zeros((Nt, Nr, Nnuobs), dtype=float)

##### functions

# H number density 
# power law, assuming spherical symmetry
def nH(r):
    return nH0*(r/rmin)**(densprof)

# normalization constant 
def n0(r):
    return nH(r)*n02nH 

# emissivity, Eq.25
def jdnu(nu, i_r, i_t):
    lamb = constCGS.C_LIGHT/nu*constCGS.cm2um   
    asub = asubarr[i_t, i_r]
    integral = 0
    for i in range(Na-1):
        # if aarr[i] <= asub or Tarr[i_t, i_r, i] < 100:
        if aarr[i] <= asub:
            continue
        aum = aarr[i]
        daum = aarr[i+1] - aarr[i]
        try:
            integral += daum*aum**(-0.5)/(aum + (lamb/lamb0)**2)/(exp(constCGS.H_PLANCK*nu/constCGS.K_B/Tarr[i_t, i_r, i]) - 1)
        except OverflowError: # the temperature is too low
            continue 
    jdnu = integral*2*pi*constCGS.H_PLANCK*nu*n0(rarr[i_r])/lamb**2
    return jdnu

##### calculate emissivity

print('Calculating...')
progress = 0

for i_t in range(Nt):
    for i_nu in range(Nnuobs):
        for i_r in range(Nr):
            jdnuarr[i_t, i_r, i_nu] = jdnu(nuobsarr[i_nu], i_r, i_t)
    if i_t/Nt > progress:
        print('{:.2%}'.format(i_t/Nt))
        progress += 0.1
print('{:.2%}'.format(1))


##### save data 

print('\nSaving data...')

jdnuarr_shape = '{}\t{}\t{}\n'.format(Nt, Nr, Nnuobs)
with open(os.path.join(folder, 'emissivity.txt'), 'w') as file:
    file.write(jdnuarr_shape)
    np.savetxt(file, jdnuarr.reshape(Nt, -1))

print('All done!')