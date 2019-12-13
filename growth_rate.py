import dispersion_solve
import particle_class
import math
from astropy import units as u
from sympy import symbols
from sympy import sympify,pi
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
    for m in range(-1,2):
        sum_term  = sum_term + electron.weight_function(m)*electron.g1_term()*electron.delta_function(m)
    gama_term = (pi ** 2 * (-electron.gyrof.value) * omega / k_wave)
    inter_term = vx**2 * sum_term
    test_term = electron.g1_term()
    # inter = Integral(f,(vx,0,oo),(vy,-oo,oo))


    return gama_term, inter_term, test_term


def growth_rate_part(B, n , T_perp, T_para, m):
    """
    This is for calculation of different wave growth part which means Landau resonance and cyclotron resonance

    In this part will not do the integration.
    The argument in Jm was k_perp * v_perp /Omega and it must be convert into r because it's hard for python to
    do the integration like Jm(0,a*x)^2 but only Jm(0,x)
    @param B:
    @param n:
    @param T_perp:
    @param T_para:
    @param m: m = -1,0,1...
    @return: The growth rate befor integration
    """
    electron = particle_class.particle('e','bi_maxwellian',n[0],B,T_perp, T_para)
    omega = symbols('w')
    r = symbols('r')
    k_wave = symbols('k')
    k_perp = symbols('k_perp')
    vx = r * electron.gyrof.value/k_perp
    vy = symbols('vy')


    coeff_term =  (pi ** 2 * (-electron.gyrof.value) * omega / k_wave)
    integral_term = vx**2 * electron.weight_function_p(m)*electron.g1_term_p()*electron.delta_function(m)

    return coeff_term, integral_term