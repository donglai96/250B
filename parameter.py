'''
This File is the parameter of wave growth fild
'''
import numpy as np
import mpmath as mp
import scipy.special
import numpy as np
from astropy import units as u
import plasmapy
import matplotlib.pyplot as plt
from plasmapy.constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, e)

from plasmapy.physics.distribution import Maxwellian_1D

emass = 9.10938291e-31
pmass = 1.67262178e-27
echarge = 1.60217657e-19
permittivity = 8.854187817e-12
permeability = 4 * np.pi * 1e-7
cspeed = 299792458
boltzmann = 1.3806488e-23

p_dens = Maxwellian_1D(v=1 * u.m / u.s,
                       T=30 * u.K,
                       particle='e',
                       v_drift=0 * u.m / u.s)
print(p_dens)

T = 3e4 * u.K
dv = 10 * u.m / u.s
v = np.arange(-5e6, 5e6, 10) * u.m / u.s


for particle in ['p', 'e']:
    pdf = Maxwellian_1D(v, T=T, particle=particle)
    integral = (pdf).sum() * dv
    print(f"Integral value for {particle}: {integral}")
    plt.plot(v, pdf, label=particle)
plt.legend()
plt.show()

