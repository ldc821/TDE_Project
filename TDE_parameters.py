#  parameters used by calculate_Tdust.py and calculate_jdnu.py

T_tde = 30000           # [K], temperature of TDE
nH0 = 1                 # [cm^{-3}], H number density
n02nH = 1.45e-15        # n0 over nH(r)
lamb0 = 2               # [um], critical wavelength
lumr = 1500             # [au], luminous radius of the blackhole
tdur = 1e6              # [s], duration of the burst
densprof = -0.5         # the exponent of the density profile
tol = 0.0001            # tolerance for asubarr

##### discretization 

# dust grain size
amin, amax = 0.01, 0.3      # [um] min/max grain radius 
Na = 500

# distance from the center  
rmin, rmax = 0.4, 100.      # [pc] radial layers
Nr = 50  

# time 
tmin, tmax = 0, tdur
Nt = 50

# observation

Nnuobs = 10
nuobsmin, nuobsmax = 2, 3   # in um

Nmu = 100

Ntobs = 100
robsmin, robsmax = 10, rmax