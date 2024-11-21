import numpy as np
import matplotlib.pyplot as plt
import astropy
import scipy.integrate as scint

#importing constants
from astropy import constants as const
from astropy import units as unit
from astropy.units.quantity_helper.function_helpers import solve

#need planck constant (h), speed of light (c), mass of an electron (m_e), mass of a proton (m_p)
#mu_e is the number of nucleons (of mass m_p) per electron. this is not a constant

mu_e = 2 #from question 1

R_0 = 7.72e8 / mu_e
M_0 = 5.67e33 / ((mu_e)**2)
rho_0 = (8 * np.pi / 3) * ((const.m_e * const.c / const.h)**3) * (const.m_p * mu_e)

rho_min = 10**-33
rad_wd_avg = 7e8 #7000km in cm since WD stars tend to be about the radius of the Earth and was are using cgs units
zero_approx = 10**-100
step = 1

def system_eq_wd(r, ystate):
    rho, m = ystate

    x = (rho) ** (1 / 3)
    gamma_x = (x ** 2) / (3 * np.sqrt(1 + x ** 2))
    # equation 8
    drho_dr = - (m * rho) / (gamma_x * r ** 2)
    # equation 9
    dm_dr = (r ** 2) * rho
    return drho_dr, dm_dr

#need to find the radius (range) of integration
def density_range(r, ystate, rho_min=rho_min):
    if ystate[0] < rho_min:
        return 0
    return ystate[0]

def solution_system_eq_wd(rho_max, rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx):
    range_density = [rho_max, 0]
    range_radius = [rad_min_approx, rad_wd]
    range_evaluation = np.linspace(rad_min_approx, rad_wd, step_size)

    return scint.solve_ivp(system_eq_wd, range_radius, range_density, events=density_range)

rho_max = 2.6e6
wd_sol = solution_system_eq_wd(rho_max, rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx)
print(wd_sol.y)
print(wd_sol.t)