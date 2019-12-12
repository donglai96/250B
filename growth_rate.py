import dispersion_solve
import particle_class
import math
from astropy import units as u
from sympy import symbols
from sympy import sympify
from sympy import cos
from sympy import Integral
from sympy import integrate
from sympy import oo



def growth_rate_electron(B, n, T_perp, T_para):
    """
    Calculate the groth rate based on Kennel 1966
    For electron only which means the effect due to resonant ions
    @param B: The magnetic field
    @param species: The species of plasma, Here use ['e', 'p']
    @param n: The density of each species

    @param T_perp: T_perpendicular
    @param T_para: T_parallel


    @return: cold wave growth rate
    """

    electron = particle_class.particle('e','bi_maxwellian',n[0],B,T_perp, T_para)
    omega = symbols('w')

    k_wave = symbols('k')
    vx = symbols('vx')
    vy = symbols('vy')
    sum_term = 0
    for m in range(-10,10):
        sum_term  = sum_term + electron.weight_function(m)*electron.g1_term()*electron.delta_function(m)
    gama_term = (math.pi ** 2 * abs(electron.gyrof.value) * omega / k_wave)
    inter_term = vy**2 * sum_term
    # inter = Integral(f,(vx,0,oo),(vy,-oo,oo))


    return gama_term, inter_term

#
#
# species = ['e', 'p']
# B = 1e-5 * u.T
# n = [1e10 * u.m ** -3, 1e10 * u.m ** -3]
# T = 300*u.K
# T_perp = T
# T_para = 2*T
# electron = particle_class.particle('e','bi_maxwellian',n[0],B,T_perp, T_para)
# theta_sub = 30
# theta = symbols('theta')
#
# m = 1
# f = electron.weight_function(m)*electron.g1_term_2()
# print(f)
# c= cos(90)
# print(c.evalf())
# vx = symbols('vx')
# vy = symbols('vy')
# omega = symbols('w')
# k_wave = symbols('k')
# sum_term = 0
# for m in range(-10, 10):
#     sum_term = sum_term + electron.weight_function(m) * electron.g1_term() * electron.delta_function(m)
# print(sympify(sum_term))
# gama_term = (math.pi ** 2 * abs(electron.gyrof.value) * omega / k_wave)
# f = vy ** 2 * sum_term
# inter = Integral(f, (vx, 0, oo), (vy, -oo, oo))
# print(sympify(inter*gama_term))
# gama = growth_rate_electron(B,n,T_perp, T_para)
# print(gama)