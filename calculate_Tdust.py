import numpy as np
from math import pi, log10, sqrt, log, exp, sqrt, atan
import constCGS

##### import parameters

from TDE_parameters import *

print('Paramters setup')

##### discretization 

# grain size
aarrlin = np.linspace(amin, amax, Na)   

# distance from the center  
# rarrlog = np.logspace(log10(rmin), log10(rmax), Nr)
rarrlin = np.linspace(rmin, rmax, Nr)

# time 
tarr = np.linspace(tmin, tmax, Nt, endpoint=False)
tarr += (tarr[1] - tarr[0])/2.
tarrday = tarr/constCGS.d2sec

# source frequency
hnumin, hnumax = 0.1, 50     # [eV]
numin, numax = hnumin*constCGS.eV2erg/constCGS.H_PLANCK, hnumax*constCGS.eV2erg/constCGS.H_PLANCK
Nnu = 100
nuarrlog = np.logspace(log10(numin), log10(numax), Nnu)

Tarr = np.zeros((Nt, Nr, Na), dtype=float)   # dust temperature
asubarr = np.full((Nt, Nr), amin)    # sublimation radii
taudarr = np.zeros((Nnu, Nr), dtype=float)   # dust extinction optical depth at each nu

print('Discretizations setup\n')

##### functions

# H number density 
# power law, assuming spherical symmetry
def nH(r):
    return nH0*(r/rmin)**(densprof)

# normalization constant 
def n0(r):
    return nH(r)*n02nH 

# absorption efficiency factor, Eq.11
# dimensionless 
def Qabs(nu, aum):
    lamb = constCGS.C_LIGHT/nu*constCGS.cm2um   
    return 1./(1 + (lamb/lamb0)**2/aum) 

# luminosity function, cgs units
# blackbody radiation
def lumV(t, nu):
    if t > tdur:
        return 0.
    try:
        B_nuT = 2*constCGS.H_PLANCK*nu**3/constCGS.C_LIGHT**2/(exp(constCGS.H_PLANCK*nu/constCGS.K_B/T_tde) - 1)
    except OverflowError:
        B_nuT = float(0)
    surfa = 4*pi*(lumr*constCGS.au2cm)**2
    return B_nuT*surfa*4*pi

# dust extinction optical depth, Eq.24
# quantities are converted cgs units
# note that amax, asub are dimensionless here
def taud(nu, i_t, i_r, amax, asubarr):
    lamb = constCGS.C_LIGHT/nu*constCGS.cm2um   
    xmax = amax*(lamb0/lamb)**2
    integral = 0
    for i in range(i_r-1):
        r = rarrlin[i]
        xsub = asubarr[i_t, i]*(lamb0/lamb)**2
        dr = rarrlin[i+1] - rarrlin[i]
        integral += dr*constCGS.pc2cm*n0(r)*(atan(sqrt(0.5*xmax)) - atan(sqrt(0.5*xsub)))
    tau = 2*sqrt(2)*pi*(lamb0/lamb)*integral/(constCGS.cm2um**2)
    return tau

# heating rate, Eq.14
def qdot_h(t, i_r, aum, taudarr):
    qdot_h = 0
    rcm = rarrlin[i_r]*constCGS.pc2cm
    acm = aum/constCGS.cm2um
    for i_nu in range(Nnu-1):
        dnu = nuarrlog[i_nu+1] - nuarrlog[i_nu]
        Cabs = pi*acm**2*Qabs(nuarrlog[i_nu],aum)
        qdot_h += dnu*lumV(t, nuarrlog[i_nu])*exp(-taudarr[i_nu, i_r])/(4*pi*rcm**2)*Cabs
    return qdot_h

# dust temperature, Eq.19
def T_d(qdot_h_cgs, aum):
    y = qdot_h_cgs/(7.12576*aum**2)
    if y >= (31.5/aum)**2/12:
        return 3240/sqrt(aum) 
    xi = sqrt((31.5/aum)**2 - 12*y) + 31.5/aum
    T_d = sqrt((2*y**2/(3*xi))**(1/3) + (xi*y/18)**(1/3))
    return T_d*1e3

# sublimation temperature, Eq.21
def T_sub(t, aum):
    return 1.76e3*(1 - 0.025*log(t/1e6/aum))

# sublimation radius, Eq.22
# the first half of the iteration
def r_sub(i_t, Tarr, taudarr, asubarr):
    asub = amin
    t = tarr[i_t]
    for i_r in range(Nr):
        for i_a in range(Na):
            aum = aarrlin[i_a]
            T = T_d(qdot_h(t, i_r, aum, taudarr), aum)
            if T > T_sub(t, aum):
                asub = aum
                Tarr[i_t, i_r, i_a] = 0.
            else:
                Tarr[i_t, i_r, i_a] = T
        if i_t == 0:
            asubarr[i_t, i_r] = asub
        else:
            asubarr[i_t, i_r] = max(asub, asubarr[i_t-1, i_r])
    return
    
# update dust extinction
# the second half of the iteration
def taud_update(i_t, asubarr, taudarr):
    for i_nu in range(Nnu):
        nu = nuarrlog[i_nu]
        for i_r in range(Nr):
            taudarr[i_nu, i_r] = taud(nu, i_t, i_r, amax, asubarr)
    return

print('Functions ready\n')

##### calculate dust temprature

print('Calculating...')
progress = 0

for i_t in range(Nt):
    frac_diff = 1
    n_iter = 0
    while frac_diff > tol:
        n_iter += 1
        asubold = np.copy(asubarr[i_t])
        r_sub(i_t, Tarr, taudarr, asubarr)
        taud_update(i_t, asubarr, taudarr)
        frac_diff = 0
        for j in range(Nr):
            frac_diff = max(frac_diff, abs(asubold[j] - asubarr[i_t, j])/asubarr[i_t, j])
        if i_t/Nt > progress:
            print('{:.2%}'.format(i_t/Nt), end='...')
            progress += 0.1

##### save data to a file

print('\n\nSaving data...')

# arrays information
with open('arrays_info.txt', 'w') as file:
    
    # information on the range of arrays
    file.write('{:>16}{:>16}{:>16}{:>16}{:>16}\n'.format('Array', 'Type', 'Start', 'End', 'Length'))
    file.write('{:>16}{:>16}{:>16f}{:>16f}{:>16d}\n'.format('tarr', 'linear', tmin, tmax, Nt))
    file.write('{:>16}{:>16}{:>16f}{:>16f}{:>16d}\n'.format('rarr', 'linear', rmin, rmax, Nr))
    file.write('{:>16}{:>16}{:>16f}{:>16f}{:>16d}\n'.format('aarr', 'linear', amin, amax, Na))

# Td_shape = '{:>5d}{:>5d}{:>5d}'.format(Nt, Nr, Na)
Td_shape = '{}\t{}\t{}\n'.format(Nt, Nr, Na)
with open('dust_temperature.txt', 'w') as file:
    file.write(Td_shape)
    np.savetxt(file, Tarr.reshape(Nt, -1))

# asub_shape = '{:>5d}{:>5d}'.format(Nr, Na)
asub_shape = '{}\t{}\n'.format(Nr, Na)
with open('sublimation_radius.txt', 'w') as file:
    file.write(asub_shape)
    np.savetxt(file, asubarr)

print('All done!')