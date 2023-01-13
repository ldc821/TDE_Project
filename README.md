## Dust Echo IR Emission from TDE

### Goal 
This code is intended to simulate dust echo emission in the infrared from a tidal disruption event

### Model assumptions
- Dust distribution: a spherically symmetric distribution of dust grains around the black hole
- Dust density profile: exponential decay with respect to distance to the black hole


### Scheme 
- First, we calculate the sublimation radius for dust grains of different sizes
- Then we integrate over all the grains that have not sublimated to get the volumetric emissivity at the observing wavelength
- Finally, we account for time-travel delays and generate the expected light curve given observer’s viewing angle and observing wavelength

### Guide to use the code
1. Run `calculate_Tdust_asub.py` to generate two files: 
    - `sublimation_radius.txt`, which contains the sublimation radius for dust of different sizes
    -  `dust_temperature.txt`, which contains the temperature for all grains at each time step, if the grain is sublimated, the temperature will be set to 0

2. Run `calculate_jdnu.py` to get the volumetric luminosity at oberserving wavelength
3. Run `calculate_Ldnu.py` to calculate the light curve at designated viewing angle and wavelength

### Visualization 
We also provide two files for visualization of the data, our plots require the matplotlib package
- `plot_dust_sublimation.py`: plot temprature of dust of different sizes, a vertical line indicates that the dust has sublimated
- `plot_light_curve.py`: plot the light curve of dust echo

This project is based on the dust echo model for binary neutron star mergers (arXiv:2108.04243)
