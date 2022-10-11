import numpy as np
from math import pi, log10, sqrt, log, exp, sqrt, atan
import constCGS

##### import parameters

from TDE_parameters import *

print('Parameters setup\n')

##### read data from file

# arrays
with open('arrays_info.txt', 'r') as file:
    file.readline()
    tarr_info = file.readline().strip('\n').split()
    rarr_info = file.readline().strip('\n').split()
    aarr_info = file.readline().strip('\n').split()

# sublimation radius
asubarr = np.loadtxt('sublimation_radius.txt', skiprows=1)

##### discretization

# tarr 
tst, tend, Nt = float(tarr_info[2]), float(tarr_info[3]), int(tarr_info[4])
if tarr_info[1] == 'linear':
    tarr = np.linspace(tst, tend, Nt)
elif tarr_info[1] == 'logarithmic':
    tarr = np.logspace(log10(tst), log10(tend), Nt)
# print(tarr)

# rarr
rst, rend, Nr = float(rarr_info[2]), float(rarr_info[3]), int(rarr_info[4])
if rarr_info[1] == 'linear':
    rarr = np.linspace(rst, rend, Nr)
elif rarr_info[1] == 'logarithmic':
    rarr = np.logspace(log10(rst), log10(rend), Nr)
rarr = np.linspace(rst, rend, Nr)
# print(rarr)

# aarr
ast, aend, Na = float(aarr_info[2]), float(aarr_info[3]), int(aarr_info[4])
if aarr_info[1] == 'linear':
    aarr = np.linspace(ast, aend, Na)
elif aarr_info[1] == 'logarithmic':
    aarr = np.logspace(log10(ast), log10(aend), Na)

# dust temperature
with open('dust_temperature.txt', 'r') as file:
    Tarr_shape = file.readline().strip('\n').split()
Tarr_shape = [int(shape) for shape in Tarr_shape]
Tarr = np.loadtxt('dust_temperature.txt', skiprows=1)
Tarr = np.reshape(Tarr, Tarr_shape)

# observation frequency
lambobs = np.linspace(nuobsmin, nuobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobs/constCGS.cm2um) # in Hz

# emissivity
jdnuarr = np.zeros((Nt, Nr, Nnuobs), dtype=float)

print('Discretizations setup\n')

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

print('Functions ready\n')

##### calculate emissivity

print('Calculating...')
progress = 0

for i_t in range(Nt):
    for i_nu in range(Nnuobs):
        for i_r in range(Nr):
            jdnuarr[i_t, i_r, i_nu] = jdnu(nuobsarr[i_nu], i_r, i_t)
    if i_t/Nt > progress:
        print('{:.2%}'.format(i_t/Nt), end='...')
        progress += 0.1


##### save data 

print('\n\nSaving data...')

jdnuarr_shape = '{}\t{}\t{}\n'.format(Nt, Nr, Nnuobs)
with open('emissivity.txt', 'w') as file:
    file.write(jdnuarr_shape)
    np.savetxt('emissivity.txt', jdnuarr.reshape(Nt, -1))

print('All done!')