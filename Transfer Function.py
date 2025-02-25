import numpy as np
import matplotlib.pyplot as plt
import psi
from scipy.integrate import simpson

import importlib
importlib.reload(psi)

bhmass = 10**8   # Msol
mdot = 1         # Msol/yr
inc = 0           # inclination, degrees

taus = np.arange(.01,150,.05)
#taus = np.logspace(-2,2,20)
combi_psi = psi.pytfb_sub(taus,1e8,mdot,4000,inc,norm=1)
# normalise 
normi = simpson(combi_psi,x=taus)
plt.plot(taus,combi_psi/normi,'green',ls='--',label='Disc',marker='.')
plt.xscale('log')
