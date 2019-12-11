"""
Distribution of plasma (bi_maxwellian distribution)

"""
import astropy as astropy
from astropy import units as u
from plasmapy.physics import parameters
import numpy as np
from scipy import integrate
from scipy.special import gamma
def _v_drift_units(v_drift):
    # Helper method to assign units to  v_drift if it takes a default value
    if (v_drift == 0 and
        not isinstance(v_drift, astropy.units.quantity.Quantity)):
        v_drift = v_drift * u.m / u.s
    else:
        v_drift = v_drift.to(u.m / u.s)
    return v_drift

def bi_maxwellian(vx,
                  vy,
                  T_per,
                  T_par,
                  particle="e",
                  vx_drift=0,
                  vy_drift=0,
                  vTh_per=np.nan,
                  vTh_par=np.nan,
                  units="units"):
    if units == "units":
        # unit checks and conversions
        # checking velocity units
        vx = vx.to(u.m / u.s)
        vy = vy.to(u.m / u.s)
        # catching case where drift velocities have default values, they
        # need to be assigned units
        vx_drift = _v_drift_units(vx_drift)
        vy_drift = _v_drift_units(vy_drift)
        # convert perpendicular temperature to Kelvins
        T_per = T_per.to(u.K, equivalencies=u.temperature_energy())
        if np.isnan(vTh_per):
            # get thermal velocity and thermal velocity squared
            vTh_per = (parameters.thermal_speed(T_per,
                                            particle=particle,
                                            method="most_probable"))
        elif not np.isnan(vTh_per):
            # check units of thermal velocity
            vTh_per = vTh_per.to(u.m / u.s)
        if np.isnan(vTh_par):
            # get thermal velocity and thermal velocity squared
            vTh_par = (parameters.thermal_speed(T_par,
                                                    particle=particle,
                                                    method="most_probable"))
        elif not np.isnan(vTh_par):
                # check units of thermal velocity
            vTh_par = vTh_par.to(u.m / u.s)
    elif np.isnan(vTh_per) and units == "unitless":
        # assuming unitless temperature is in Kelvins
        vTh_per = parameters.thermal_speed(T_per * u.K,
                                       particle=particle,
                                       method="most_probable").si.value
        # convert parallel temperature to Kelvins
        T_par = T_par.to(u.K, equivalencies=u.temperature_energy())





    # accounting for thermal velocity in 2D
    vThSq_per = vTh_per ** 2
    vThSq_par = vTh_par ** 2

    # Get square of relative particle velocity
    vSq_per = (vx - vx_drift) ** 2
    vSq_par = (vy - vy_drift) ** 2
    # calculating distribution function
    coeff = (vThSq_per * np.pi) ** (-1 / 2) * (vThSq_par * np.pi) ** (-1 / 2)
    expTerm = np.exp(-vSq_par / vThSq_par) * np.exp(-vSq_per / vThSq_per)
    distFunc = coeff * expTerm
    if units == "units":
        return distFunc.to((u.s / u.m)**2)
    elif units == "unitless":
        return distFunc

p_dens = bi_maxwellian(vx= 0*u.m/u.s,
                       vy= 0.5*u.m/u.s,
                       T_per = 30 *u.K,
                       T_par = 30 *u.K,
                       particle = 'e'
                       )
print(p_dens)

from plasmapy.physics.distribution import Maxwellian_velocity_2D

p_dens_2  = Maxwellian_velocity_2D(vx= 0*u.m/u.s,
                       vy= 0.5*u.m/u.s,
                       T = 30*u.K,
                       particle = 'e'
                       )
print(p_dens_2)

T = 3e4 * u.K
dv_par = 10 * u.m / u.s
v_parallel = np.arange(-5e6, 5e6, 10) * u.m / u.s
v_parallel = np.arange(0, 5e3, 10) * u.m / u.s
# Intergrate
sum_v = 0
# for v_perp in v_parallel:
#     pdf = bi_maxwellian(vx= v_perp, vy =v_parallel,T_per = T, T_par=T, particle='e')
#     print(pdf.sum())
# print(sum_v)
T_per = T
vth_per = parameters.thermal_speed(T = T_per ,particle ='e')
print(vth_per)
def bi_max(vx,vy,vth_per,vth_par):
    return ((vth_per**2 * np.pi) ** (-1 / 2) * (vth_par**2 * np.pi) ** (-1 / 2)) *\
           np.exp(-vx**2 / vth_per**2) * np.exp(-vy**2 / vth_par**2)
v, err = integrate.dblquad(bi_max, 0, np.inf, - np.inf, np.inf, args=(23,2))
print(v)

class particles:
    particle_name = 'ion'

c = particles
d = particles
d.particle_name = 'electron'
print(c.particle_name)
print(d.particle_name)