import numpy as np
import matplotlib.pyplot as plt
import os

# Put path to the folder with the dll here (the whole thing)
dll_path = r"C:\Users\... ...\Comp_SED\SED\.libs"
os.add_dll_directory(dll_path)

import SED

lday = 2.998e8 * 3600*24    # lday to m conversion
Ang = 1e-10                 # Ang to m conversion
deg = np.pi/180             # deg to rad conversion

r = np.array([1, 2, 3]) * lday
dr = np.array([1, 1, 1]) * lday
T = np.array([1e5, 5e4, 1e4])

nr = len(r)
nw = 100

wl_min = 5e2 * Ang
wl_max = 1e5 * Ang

SED_wl_array = np.logspace(np.log10(wl_min), np.log10(wl_max), nw)  # Note - log spacing is not necessary
                                                                    # Can also work with linear spacing

incl = 45 * deg

# Note - Please call the function as follows, otherwise things tend to break
Lnu = SED.sed(
    nr      = nr,           # Number of radial elements
    r_r     = r,            # Radius of the element (in m)
    dr_r    = dr,           # Radial width of the element (in m)
    t_r     = T,            # Temperature (in K)
    nw      = nw,           # Number of wavelengths to compute SED at
    wl      = SED_wl_array, # List of wl values to compute SED at (in m)
    incl    = incl)         # inclination of the object (0 being disc face-on)


fig, ax = plt.subplots(1,1, figsize=(8,6))

ax.loglog(SED_wl_array/Ang, Lnu, linewidth=2, c='C0')
ax.set_xlabel('Wavelength (Å)', fontsize=14)
ax.set_ylabel(r'$L_\nu$ (erg/s/Hz)', fontsize=14)
ax.set_xlim(400, 130000)
ax.grid(alpha=0.5)
ax.tick_params(which='both', top=True, right=True, direction='in', labelsize=12)
ax.set_title('Disc Spectrum (Continuum)', fontsize=18) 

plt.show()