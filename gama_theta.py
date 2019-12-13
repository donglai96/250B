import dispersion_solve
import growth_rate
import particle_class
from astropy import units as u
from plasmapy.physics import parameters
from sympy.solvers import solve
from sympy import symbols
from sympy import Matrix
from sympy import cos,sin,pi
from sympy import simplify, integrate, oo,exp
import math
from sympy import solveset, S
from astropy import constants as const
import numpy as np
from sympy import DiracDelta,besselj

# First set the background value of this plasma question
def dispersion_matrix():
    L = symbols('L')
    R = symbols('R')
    P = symbols('P')
    theta_rad = symbols('theta')
    nn = symbols('nn',positive = True)
    m_11 = 2 * (R - nn ** 2 + 0.5 * nn ** 2 * (sin(theta_rad)) ** 2)
    m_12 = nn ** 2 * sin(theta_rad) ** 2
    m_13 = nn ** 2 * cos(theta_rad) * sin(theta_rad)
    m_21 = m_12
    m_22 = 2 * (L - nn ** 2 + 0.5 * nn ** 2 * (sin(theta_rad)) ** 2)
    m_23 = m_13
    m_31 = m_13
    m_32 = m_13
    m_33 = P - nn ** 2 * (sin(theta_rad) ** 2)
    D_0 = Matrix([[m_11, m_12, m_13], [m_21, m_22, m_23], [m_31, m_32, m_33]])
    return D_0

def solve_dispersion(B, species, n, wave_frequency, theta_input):
    L = symbols('L')
    R = symbols('R')
    P = symbols('P')
    theta = symbols('theta')



    # Get the function
    theta_deg = theta_input *pi/180
    D_0 = dispersion_matrix()
    L_wave, R_wave, P_wave = dispersion_solve.cold_plasma_LRP(B, species, n, wave_frequency)
    f = simplify(D_0.det())
    print("D_0.det is ")
    print(f.subs([(L, L_wave), (P, P_wave), (R, R_wave), (theta, theta_deg)]))
    D_0_sim = f.subs([(L, L_wave), (P, P_wave), (R, R_wave), (theta, theta_deg)])
    refractive_index = solve(D_0_sim)
    return refractive_index


# Solve the nn
#
species = ['e','p']
T_perp = 100 * u.K
T_para = 50 * u.K
B = 1e-6 * u.T
n = [300e6 * u.m ** -3, 300e6 * u.m ** -3]
theta_input = 30
theta_deg = theta_input *pi/180
gyro_frequency = parameters.gyrofrequency(B, 'e')
wave_frequency = 0.05 * gyro_frequency
print(wave_frequency)

# solve L,R,P

re_index = solve_dispersion(B,species,n,wave_frequency,theta_input)
print(re_index[0])
nn = re_index[0]
#
# # solve the gamma
gamma_term, inter_term, test_term = growth_rate.growth_rate_electron(B,n,T_perp, T_para)
# print(inter_term.free_symbols)
# print('gamma_term',gamma_term)
# print('inter_term',inter_term)
theta = symbols('theta')
k_para = symbols('k_para')
k_perp = symbols('k_perp')
k = symbols('k')
w = symbols('w')
vx = symbols('vx')
vy = symbols('vy')
#
kk = (nn * wave_frequency/const.c).value

kk_para = kk * cos(theta_deg)
kk_perp = kk * sin(theta_deg)
#
inter_term_new = inter_term.subs([(theta,theta_deg),(k_para,kk_para),(k,kk),(k_perp,kk_perp),(w,wave_frequency.value)])
inter_term_new = simplify(inter_term_new)
print(inter_term_new)
print("test_term ",test_term)
print(simplify(test_term))
integrate_gamma_vy = integrate(inter_term_new,(vy,-oo,oo))
print(integrate_gamma_vy)
print("success")
# intergrate_gamma_vx = integrate(inter_term_new,(vx,0,oo))
# print(intergrate_gamma_vx)
# print("success!!")
a = symbols('a')
b = symbols('b')
print("start a new test!!!")



testvx = vy*test_term *DiracDelta(vy-1)
test_inte = integrate(testvx,(vy,-oo,oo))
print(test_inte)
# test for calculation of bessel integration
testvxx = besselj(1, vx)*besselj(0,vx) *exp(-vx**2*a)*vx**3*b
print(testvxx)
test_vxx_result = integrate(testvxx,(vx,0,oo))
print(test_vxx_result)
print("success!!!")
print('test_result',test_vxx_result.subs([(a,1)]).evalf())

####
####
# After test, for sympy when you calculate integral like j_m(0,ax)^2 it not work but if you remove a it will work
# so the normalization is important. So back to before I need to change vx to kxvx/Omega
# def growth_inter(m):

###########################
# start a new test
m = 0
r = symbols('r')
coeff_term, int_term = growth_rate.growth_rate_part(B,n,T_perp, T_para,m)
print("The new test is ")
print("***************")
print(int_term)
int_term_vy = integrate(int_term,(vy,-oo,oo))
print("int vy sucess, the result is ", int_term_vy)
int_term_vx = integrate(int_term_vy,(r,-oo,0))
print("int vx sucess, the result is ", int_term_vx)