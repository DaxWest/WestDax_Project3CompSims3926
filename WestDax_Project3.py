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

x = (rho / rho_0)**(1/3)
M = m / M_0
R = r / R_0

gamma_x = (x**2) / (3 * np.sqrt(1 + x**2))

#equation 8
drho_dr = - (M * rho) / (gamma_x * r**2)

#equation 9
dm_dr = (R**2) * rho