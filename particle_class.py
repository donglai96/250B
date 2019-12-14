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
from sympy import cos,pi
import dispersion_solve

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
        vth_x = symbols('vth_x')
        vth_y = symbols('vth_y')

        # coeff = ((self.vth_perp.value**2 * np.pi) ** (-1 / 2) * (self.vth_para.value**2 * np.pi) ** (-1 / 2))
        # The vx means v perpendicular!!
        # coeff = pi**(-3/2)*(self.vth_perp.value**-2)*self.vth_para.value**-1
        # expterm = exp(-vx**2 / self.vth_perp.value**2) * exp(-vy**2 / self.vth_para.value**2)

        # coeff = pi ** (-3 / 2) * (vth_x ** -2) * vth_y ** -1
        # expterm = exp(-vx**2 / vth_x**2) * exp(-vy**2 / vth_y**2)

        coeff = pi ** (-2 / 2) * (vth_x ** -1) * vth_y ** -1
        expterm = exp(-vx ** 2 / vth_x ** 2) * exp(-vy ** 2 / vth_y ** 2)

        return coeff*expterm

    def g1_term(self):
        """
        From Kennel 1966

        G_{1}^{\pm}=\frac{\partial F^{\pm}}{\partial v_{\perp}}-\frac{k_{\|}}{\omega}\left(v_{\|}
        \frac{\partial F^{\pm}}{\partial v_{\perp}}-v_{\perp} \frac{\partial F^{\pm}}{\partial v_{\|}}\right)
        @return: A function with vx and vy
        """
        # vx is v_perpendicular and vy is v_parallel
        vx = symbols('vx')
        vy = symbols('vy')
        k_para = symbols('k_para')
        omega = symbols('w')
        return diff(self.F(),vx) - k_para * (vy * diff(self.F(),vx) - vx * diff(self.F(),vy))/omega
    def g1_term_2(self):
        vx = symbols('vx')
        vy = symbols('vy')
        k_para = symbols('k_para')
        omega = symbols('w')
        return diff(self.F(),vx)

    def g1_term_p(self):
        """
        change vx to r
        @return:
        """
        k_perp = symbols('k_perp')
        vx = symbols('vx')
        r = symbols('r')
        vv_perp = r * self.gyrof.value/ k_perp
        c = self.g1_term().subs([(vx,vv_perp)])
        return c
    def g1_term_parallel(self):
        vx = symbols('vx')
        vy = symbols('vy')
        omega = symbols('w')
        k = symbols('k')
        v_para = (omega + self.gyrof.value)/k

        g1 = diff(self.F(), vx)
        g1_para = g1.subs([(vy,v_para )])
        return g1_para
    def f_term_not(self):
        vx = symbols('vx')
        vy = symbols('vy')
        omega = symbols('w')
        k = symbols('k')
        v_para = (omega + self.gyrof.value)/k

        g1 = self.F()
        g1_para = g1.subs([(vy,v_para )])
        return g1_para




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
    def Jm_p(self,m):
        """
        The argument of Jm change to r
        r = k_perp*v_perp/Omega
        @param m: the order of Bessel function
        @return:
        """
        r = symbols('r')
        argument = r
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
    def weight_function(self,m):
        """
        from Kennel 1966

        @param m:
        @return: function with theta
        """
        theta = symbols('theta')
        gamma_m = (((1 + cos(theta))*self.Jm(m+1) + (1-cos(theta))*self.Jm(m-1))/(2*cos(theta)))**2
        return gamma_m
    def weight_function_p(self,m):
        theta = symbols('theta')
        gamma_m = (((1 + cos(theta))*self.Jm_p(m+1) + (1-cos(theta))*self.Jm_p(m-1))/(2*cos(theta)))**2

        return gamma_m




#
# species = ['e', 'p']
# B = 1e-5 * u.T
# n = [1e10 * u.m ** -3, 1e10 * u.m ** -3]
# T = 300*u.K
#
#
#
# ion = particle('p','bi',n[0], B,T,T)
# electron = particle('e', 'bi', n[0], B,T,2*T)
# print (ion.name)
# print(ion.gyrof)
# print(ion.plasmaf)
# print(electron.name)
# print(electron.gyrof)
# print(electron.plasmaf)
# print(electron.F())
# vx = symbols('vx')
# print(diff(vx *electron.F(),vx))
# print('test',electron.g1_term_2())
#
# from sympy.abc import y,n
#
# print("test", electron.Jm(0))
# print('gyro',electron.gyrof)
# print('Delta', electron.delta_function(1))
# print('Weight',electron.weight_function(1))