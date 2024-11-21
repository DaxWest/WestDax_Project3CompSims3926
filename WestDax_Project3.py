import numpy as np
import matplotlib.pyplot as plt
import astropy

#importing constants
from astropy import constants as const
from astropy import units as unit
#need planck constant (h), speed of light (c), mass of an electron (m_e), mass of a proton (m_p)
#mu_e is the number of nucleons (of mass m_p) per electron. this is not a constant

mu_e = 2 #from question 1

R_0 = 7.72e8 / mu_e
M_0 = 5.67e33 / ((mu_e)**2)
rho_0 = (8 * np.pi / 3) * ((const.m_e * const.c / const.h)**3) * (const.m_p * mu_e)

rho_min = 10**-33

def system_eq_wd(r, ystate, R_0=1, M_0=1, rho_0=1):
    rho, m = ystate

    RHO = rho / rho_0
    x = (RHO) ** (1 / 3)
    M = m / M_0
    R = r / R_0
    gamma_x = (x ** 2) / (3 * np.sqrt(1 + x ** 2))
    # equation 8
    drho_dr = - (M * RHO) / (gamma_x * R ** 2)
    # equation 9
    dm_dr = (R ** 2) * RHO
    return drho_dr, dm_dr

#need to find the radius (range) of integration
def density_range(r, ystate, rho_min=0):
    if ystate[0] < rho_min:
        return 0
    return ystate[0]

