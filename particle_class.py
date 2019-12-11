import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

from plasmapy.physics import parameters
from sympy import *
import math
from astropy import constants as const
from scipy import optimize
import scipy.special as spl

from astropy import units as u
from plasmapy.physics import parameters
import mpmath

from scipy import integrate
from scipy.special import gamma
from sympy import besselj,jn
from sympy import DiracDelta

class particle():
    def __init__(self,name,distribution,density,magnetic_field, T_perp, T_para):
        self.name = name
        self.distribution = distribution
        self.density = density
        self.magnetic_field = magnetic_field
        self.gyrof = parameters.gyrofrequency(B=magnetic_field, particle=name, signed=True)
        self.plasmaf = parameters.plasma_frequency(density,particle=name)
        self.vth_perp = parameters.thermal_speed(T_perp, name)
        self.vth_para = parameters.thermal_speed(T_para, name)
    def F(self):
        """
        The distribution of velocity, now is only bi-maxwellian distribution
        @return: A function with vx and vy
        """
        vx = symbols('vx')
        vy = symbols('vy')
        coeff = ((self.vth_perp.value**2 * np.pi) ** (-1 / 2) * (self.vth_para.value**2 * np.pi) ** (-1 / 2))
        expterm = exp(-vx**2 / self.vth_perp.value**2) * exp(-vy**2 / self.vth_para.value**2)

        return coeff*expterm

    def g1_term_1(self):
        """
        From Kennel 1966

        G_{1}^{\pm}=\frac{\partial F^{\pm}}{\partial v_{\perp}}-\frac{k_{\|}}{\omega}\left(v_{\|}
        \frac{\partial F^{\pm}}{\partial v_{\perp}}-v_{\perp} \frac{\partial F^{\pm}}{\partial v_{\|}}\right)
        @return: A function with vx and vy
        """
        # vx is v_perpendicular and vy is v_parallel
        vx = symbols('vx')
        vy = symbols('vy')
        return diff(self.F(),vx)
    def g1_term_2(self):
        vx = symbols('vx')
        vy = symbols('vy')

        return vy * diff(self.F(),vx) - vx * diff(self.F(),vy)
    def Jm(self,m):
        """
        The argument of Jm is k_perp * v_perp(vx)/gyro
        @param m: the order of Bessel function
        @return: Jm(k_perp * vx/Omega)

        """
        k_perp = symbols('k_perp')
        vx = symbols('vx')
        argument = k_perp * vx /self.gyrof.value
        return besselj(m,argument)
    def delta_function(self,m):
        """
        The argument of Dirac is v_para - (w -m*Omega/k_parallel)
        @param m: The order
        @return: Dirac function of argument
        """
        vy = symbols('vy')
        k_para = symbols('k_para')
        w = symbols('w')
        argument = vy - (w - m * self.gyrof.value)/k_para
        return DiracDelta(argument)








species = ['e', 'p']
B = 1e-5 * u.T
n = [1e10 * u.m ** -3, 1e10 * u.m ** -3]
T = 300*u.K



ion = particle('p','bi',n[0], B,T,T)
electron = particle('e', 'bi', n[0], B,T,2*T)
print (ion.name)
print(ion.gyrof)
print(ion.plasmaf)
print(electron.name)
print(electron.gyrof)
print(electron.plasmaf)
print(electron.F())
vx = symbols('vx')
print(diff(vx *electron.F(),vx))
print('test',electron.g1_term_2())

from sympy.abc import y,n

print("test", electron.Jm(0))
print('gyro',electron.gyrof)
print('Delta', electron.delta_function(1))