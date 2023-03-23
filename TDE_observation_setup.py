import constCGS

##### observation setup #####

lambobsmin, lambobsmax = 2, 4  # [um] observation wavelength
Nnuobs = 10

Nmu = 100                   

Ntobs = 100
robsmin, robsmax = 0, 200
tobsmin = 100
tobsmax = 200*constCGS.pc2cm/constCGS.C_LIGHT