import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.integrate as scint
from decimal import Decimal

mu_e = 2 #from question 1

R_sun = 6.957e8
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
rad_wd_avg = 7e6 #7000km in m since this is the upper limit of the radii of WD stars
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

radius_sol_limit = []
density_sol_limit = []
mass_sol_limit = []

for i in range(len(rho_c)):
    wd_sol = solution_system_eq_wd(rho_c[i], rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx)

    rad_sol = wd_sol.t * R_0 / R_sun
    #y[0] is density and y[1] is mass since system_eq_wd returns drho_dr first and dm_dr second
    den_sol = wd_sol.y[0] * rho_0
    m_sol = wd_sol.y[1] * M_0 / M_sun

    radius_sol.append(rad_sol)
    density_sol.append(den_sol)
    mass_sol.append(m_sol)

    rad_sol_limit = wd_sol.t_events[0][0] * R_0 / R_sun
    den_sol_limit = wd_sol.y_events[0][0][0] * rho_0
    m_sol_limit = wd_sol.y_events[0][0][1] * M_0 / M_sun

    radius_sol_limit.append(rad_sol_limit)
    density_sol_limit.append(den_sol_limit)
    mass_sol_limit.append(m_sol_limit)

#Question 2
rho_legend = []
for i in range(len(rho_c)):
    legend = f'$\\rho$ = {Decimal(rho_c[i]):.2e}'
    rho_legend.append(legend)

fig = plt.figure(figsize=(8,6))
for i in range(len(density_sol)):
    plt.plot(mass_sol[i], radius_sol[i], ls=':')
plt.title('Mass vs Radius for Various Maximum Internal Densities')
plt.xlabel('Mass (solar masses)')
plt.ylabel('Radius (solar radii)')
plt.legend(rho_legend)
for i in range(len(density_sol_limit)):
    plt.scatter(mass_sol_limit[i], radius_sol_limit[i])
plt.show()

M_ch_estimation = []
for i in range(len(mass_sol)):
    M_ch_estimation.append(np.max(mass_sol[i]))

print(f'The estimation of the Chandrasekhar mass limit is: {Decimal(np.max(M_ch_estimation)):.4} Solar Masses')
print(f'The true value of the Chandrasekhar mass limit is: {M_ch} Solar Masses')

#Question 3
rho_c_3 = [rho_c[1], rho_c[4], rho_c[8]]
radius_sol_3 = []
density_sol_3 = []
mass_sol_3 = []

radius_sol_3_limit = []
density_sol_3_limit = []
mass_sol_3_limit = []

for i in range(len(rho_c_3)):
    wd_sol_3 = solution_system_eq_wd(rho_c_3[i], rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx, method='RK23')

    rad_sol_3 = wd_sol_3.t * R_0 / R_sun
    den_sol_3 = wd_sol_3.y[0] * rho_0
    m_sol_3 = wd_sol_3.y[1] * M_0 / M_sun

    radius_sol_3.append(rad_sol_3)
    density_sol_3.append(den_sol_3)
    mass_sol_3.append(m_sol_3)

    rad_sol_3_limit = wd_sol_3.t_events[0][0] * R_0 / R_sun
    den_sol_3_limit = wd_sol_3.y_events[0][0][0] * rho_0
    m_sol_3_limit = wd_sol_3.y_events[0][0][1] * M_0 / M_sun

    radius_sol_3_limit.append(rad_sol_3_limit)
    density_sol_3_limit.append(den_sol_3_limit)
    mass_sol_3_limit.append(m_sol_3_limit)

legend_method_comp = [f'$\\rho$ = {Decimal(rho_c[1]):.2e}', f'$\\rho$ = {Decimal(rho_c[4]):.2e}', f'$\\rho$ = {Decimal(rho_c[8]):.2e}']
fig2, (ax1, ax2) = plt.subplots(1,2, figsize=(15,5))
for i in range(len(density_sol_3)):
    ax1.plot(mass_sol_3[i], radius_sol_3[i], ls=':')
ax1.legend(legend_method_comp)
ax2.legend(legend_method_comp)
for i in range(len(density_sol_3_limit)):
    ax1.scatter(mass_sol_3_limit[i], radius_sol_3_limit[i])
ax2.plot(mass_sol[1], radius_sol[1], ls=':')
ax2.plot(mass_sol[4], radius_sol[4], ls=':')
ax2.plot(mass_sol[8], radius_sol[8], ls=':')
ax2.scatter(mass_sol_limit[1], radius_sol_limit[1])
ax2.scatter(mass_sol_limit[4], radius_sol_limit[4])
ax2.scatter(mass_sol_limit[8], radius_sol_limit[8])
ax1.set_title('Method: RK23')
ax2.set_title('Method: RK45')
ax1.legend(legend_method_comp)
ax2.legend(legend_method_comp)
ax1.set_xlabel('Mass (solar masses)')
ax2.set_xlabel('Mass (solar masses)')
ax1.set_ylabel('Radius (solar radii)')
ax2.set_ylabel('Radius (solar radii)')
plt.show()

#Question 4
wd_mass = []
wd_mass_unc = []
wd_radius = []
wd_radius_unc = []
with open('wd_mass_radius.csv', 'r') as wd_mass_data:
    for line in wd_mass_data.readlines()[1:]:
        cols = line.strip().split(',')
        wd_mass.append(float(cols[0]))
        wd_mass_unc.append(float(cols[1]))
        wd_radius.append(float(cols[2]))
        wd_radius_unc.append(float(cols[3]))

fig3 = plt.figure(figsize=(8,6))
for i in range(len(density_sol)):
    plt.plot(mass_sol[i], radius_sol[i], ls=':')
plt.legend(rho_legend)
plt.xlabel('Mass (solar masses)')
plt.ylabel('Radius (solar radii)')
plt.title('Observation vs Simulation Comparison')
for i in range(len(density_sol_limit)):
    plt.scatter(mass_sol_limit[i], radius_sol_limit[i])
plt.scatter(wd_mass, wd_radius)
plt.errorbar(wd_mass, wd_radius, xerr=wd_mass_unc, yerr=wd_radius_unc, color='k', ls='')
plt.show()