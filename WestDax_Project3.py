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

R_sun = 695700000
M_sun = 1.989e30
M_ch = 5.836/ ((mu_e)**2)
R_0 = 7.72e6 / mu_e
M_0 = 5.67e30 / ((mu_e)**2)

m_e = 9.109e-31
m_p = 1.672e-27
c = 2.998e8
h = 6.626e-34

rho_0 = (8 * np.pi / 3) * ((m_e * c / h)**3) * (m_p * mu_e)

#Question 1
rho_min = 10**-1
rad_wd_avg = 7e8 #7000km in cm since WD stars tend to be about the radius of the Earth and was are using cgs units
zero_approx = 10**-10
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
    if ystate[0] - rho_min < 0:
        ystate[0] = 0
    return ystate[0]

def solution_system_eq_wd(rho_max, rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx, method='RK45'):
    range_density = [rho_max, 0]
    range_radius = [rad_min_approx, rad_wd]
    # range_evaluation = np.linspace(rad_min_approx, rad_wd, step_size) #not sure why i made this and i don't think i need it
    density_range.terminal = True
    return scint.solve_ivp(system_eq_wd, range_radius, range_density, events=density_range, method=method)

rho_c = np.logspace(0, 6.39, 10) #used trial and error to get upper limit. this was very close to upper bound given in the question

radius_sol = []
density_sol = []
mass_sol = []

for i in range(len(rho_c)):
    wd_sol = solution_system_eq_wd(rho_c[i], rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx)

    rad_sol = wd_sol.t * R_0 / R_sun
    #y[0] is density and y[1] is mass since system_eq_wd returns drho_dr first and dm_dr second
    den_sol = wd_sol.y[0] * rho_0
    m_sol = wd_sol.y[1] * M_0 / M_sun

    radius_sol.append(rad_sol)
    density_sol.append(den_sol)
    mass_sol.append(m_sol)

#Question 2
fig = plt.figure()
for i in range(len(density_sol)):
    plt.plot(mass_sol[i], radius_sol[i])
plt.show()

M_ch_estimation = []
for i in range(len(mass_sol)):
    M_ch_estimation.append(np.max(mass_sol[i]))

print(f'The estimation of the Chandrasekhar mass limit is: {np.max(M_ch_estimation)} Solar Masses')
print(f'The true value of the Chandrasekhar mass limit is: {M_ch} Solar Masses')

#Question 3
rho_c_3 = [rho_c[1], rho_c[4], rho_c[8]]
radius_sol_3 = []
density_sol_3 = []
mass_sol_3 = []

for i in range(len(rho_c_3)):
    wd_sol_3 = solution_system_eq_wd(rho_c_3[i], rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx, method='RK23')

    rad_sol_3 = wd_sol_3.t * R_0
    den_sol_3 = wd_sol_3.y[0] * rho_0
    m_sol_3 = wd_sol_3.y[1] * M_0

    radius_sol_3.append(rad_sol_3)
    density_sol_3.append(den_sol_3)
    mass_sol_3.append(m_sol_3)

fig2 = plt.figure()
for i in range(len(density_sol_3)):
    plt.plot(mass_sol_3[i], radius_sol_3[i])
plt.show()

#Question 4
