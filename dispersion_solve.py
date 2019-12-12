"""
This file is the function of solving the Dispersion relationship
The function is based on Kennel 1966 Low-Frequency Whistler Mode

"""

"""
Author : Donglai Ma  UCLA AOS
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

from plasmapy.physics import parameters
from sympy import *
import math
from astropy import constants as const
from scipy import optimize
import scipy.special as spl

def cold_plasma_LRP(B, species, n, omega ,flag= True ):
    """

    @param B: Magnetic field in z direction
    @param species: ['e','p'] means electrons and protons
    @param n: The density of different species
    @param omega: wave_frequency
    @param flag: if flag == 1 then omega means to be a number, if flag is not True the omega means a symbol
    @return: (l, R, P)

    """

    if flag:
        wave_frequency = omega
        w = wave_frequency.value
    else:
        wave_frequency = symbols('w')
        w = wave_frequency
    L, R, P = 1, 1, 1

    for s, n_s in zip(species, n):
        omega_c = parameters.gyrofrequency(B=B, particle=s, signed=True)
        omega_p = parameters.plasma_frequency(n=n_s, particle=s)

        L += - omega_p.value ** 2 / (w * (w - omega_c.value))
        R += - omega_p.value ** 2 / (w * (w + omega_c.value))
        P += - omega_p.value ** 2 / w ** 2
    return L, R, P



def dispersion_matrix(B, species, n, theta,  omega = 0 , k_wave = 0, mode = 'w'):
    """

    @param B: Magnetic field in z direction
    @param species: ['e','p'] means electrons and protons
    @param n: The density of different species
    @param omega: wave_frequency
    @param k: wave_number
    @param mode: if mode = w it means the dispersion relationship is in wave_frequency space
                If mode = k it means the dispersion relationship is in wave_number space

    @return: Matrix of D_0
    """
    ## In the wave frequency space, the result is represent in k

    if mode =='k':
        k = k_wave

        L, R, P = cold_plasma_LRP(B,species,n,omega,flag = False)
    elif mode == 'w':
        k = symbols('k')

        L, R, P = cold_plasma_LRP(B,species,n,omega,flag = True)

    print(k)
    k_unit = (u.rad/u.m)*const.c/omega
    nn = k * k_unit

    print(nn)

    theta_rad = theta * np.pi / 180

    m_11 = 2 * (R - nn ** 2 + 0.5 * nn ** 2 * (math.sin(theta_rad)) ** 2)
    m_12 = nn ** 2 * math.sin(theta_rad) ** 2
    m_13 = nn ** 2 * math.cos(theta_rad) * math.sin(theta_rad)
    m_21 = m_12
    m_22 = 2 * (L - nn ** 2 + 0.5 * nn ** 2 * (math.sin(theta_rad)) ** 2)
    m_23 = m_13
    m_31 = m_13
    m_32 = m_13
    m_33 = P - nn ** 2 * (math.sin(theta_rad) ** 2)

    D_0 = Matrix([[m_11, m_12, m_13], [m_21, m_22, m_23], [m_31, m_32, m_33]])
    return D_0


# species = ['e','p']
# B = 1e-5 * u.T
# n = [1e10*u.m**-3, 1e10*u.m**-3]
# gyro_frequency = parameters.gyrofrequency(B,'e')
# wave_frequency = 0.1 * gyro_frequency
# print(cold_plasma_LRP(B,species,n,wave_frequency,flag=True))
# print(dispersion_matrix(B,species,n, theta = 30, k_wave= 1,mode = 'k'))






