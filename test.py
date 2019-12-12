
from sympy import symbols
import math
from sympy import symbols
from sympy import Matrix
from sympy import cos,sin
from sympy import simplify

L =symbols('L')
R = symbols('R')
P = symbols('P')
theta_rad = symbols('theta')
nn = symbols('nn')



m_11 = 2*(R - nn**2 +0.5*nn**2*(sin(theta_rad))**2)
m_12 = nn**2 *sin(theta_rad)**2
m_13 = nn**2 *cos(theta_rad)*sin(theta_rad)
m_21 = m_12
m_22 = 2*(L - nn**2 +0.5*nn**2*(sin(theta_rad))**2)
m_23 = m_13
m_31 = m_13
m_32 = m_13
m_33 = P - nn**2 *(sin(theta_rad)**2)
D_0 = Matrix([[m_11, m_12, m_13], [m_21, m_22, m_23], [m_31, m_32, m_33]])
print(simplify(D_0.det()))