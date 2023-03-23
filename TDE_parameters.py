##### model parameters #####

# tidal disruption events
T_tde = 3e4             # [K], temperature of TDE
lumr = 100              # [au], luminous radius of the blackhole
tdur = 1e6              # [s], duration of the burst

###### ASASSN-14li
# T_tde = 3.5e4
# lumr = 100
# tdur = 178*24*60
# t0 = 60*24*60

###### PTF-09ge
# T_tde = 2.2e4
# lumr = 100
# tdur = 1e6

# dust distribution
nH0 = 1                # [cm^{-3}], H number density
densprof = -0.5         # the exponent of the density profile
n02nH = 1.45e-15        # n0 over nH(r)
lamb0 = 2               # [um], critical wavelength

##### simulation setup #####

amin, amax = 0.01, 0.3      # [um] min/max grain radius
Na = 30

rmin, rmax = 0.4, 100.      # [pc] radial layers
Nr = 100

tmin, tmax = 0, tdur        # [s]  time since start of TDE events
Nt = 50

hnumin, hnumax = 0.1, 50    # [eV] source frequency
Nnu = 100

tol = 0.1                   # fractional tolerance for sublimation radius

##### folder to store the outputs

name = "data"
# name = "ASASSN14li"
# name = "PTF09ge"
folder = '{}_Ttde{:.1e}_tdur{:.1e}_nH0{:.2e}^{:.2f}_lamb0{:d}um'.format(name, T_tde, tdur, nH0, densprof, lamb0)