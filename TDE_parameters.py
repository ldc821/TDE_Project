##### model parameters #####

# tidal disruption events
T_tde = 30000           # [K], temperature of TDE
lumr = 1500             # [au], luminous radius of the blackhole
tdur = 1e6              # [s], duration of the burst

# dust distribution
nH0 = 1                 # [cm^{-3}], H number density
n02nH = 1.45e-15        # n0 over nH(r)
lamb0 = 2               # [um], critical wavelength
densprof = -0.5         # the exponent of the density profile

##### simulation setup #####

amin, amax = 0.01, 0.3      # [um] min/max grain radius
Na = 50

rmin, rmax = 0.4, 100.      # [pc] radial layers
Nr = 50

tmin, tmax = 0, tdur        # [s]  time since start of TDE events
Nt = 50

hnumin, hnumax = 0.1, 50    # [eV] source frequency
Nnu = 100

tol = 0.1                   # fractional tolerance for sublimation radius

##### folder to store the outputs
folder = 'data_Ttde{:.2e}K_nH0{:.2e}_lamb0{:d}um'.format(T_tde, nH0, lamb0)