import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
from decimal import Decimal

mu_e = 2 #from question 1, number of nucleons per electron

R_sun = 6.957e8 #radius of the sun in m
M_sun = 1.989e30 #mass of the sun in kg
M_ch = 5.836/ ((mu_e)**2) #Chandrsekhar mass limit
R_0 = 7.72e6 / mu_e #radius in m/(number of nucleons per electron)
M_0 = 5.67e30 / ((mu_e)**2) #mass in kg/(number of nucleons per electron)^2

m_e = 9.109e-31 #mass of an electron in kg
m_p = 1.672e-27 #mass of a proton in kg
c = 2.998e8 #speed of light in m/s
h = 6.626e-34 #planck's constant

rho_0 = (8 * np.pi / 3) * ((m_e * c / h)**3) * (m_p * mu_e) #density given by known constants

#Question 1
rho_min = 10**-1 #minimum density
rad_wd_avg = 7e6 #7000km in m since this is the upper limit of the radii of WD stars
zero_approx = 10**-10 #zero approximation for termination function
step = 1 #timestep size

def system_eq_wd(r, ystate):
    '''
    Returns the results of equations 8 and 9 from the project description for use with scint.solve_ivp
    :param r: radius
    :param ystate: state vector for density and mass
    :return: equations 8 and 9
    '''
    rho, m = ystate

    #all equations used in this function were given in the project description
    x = (rho) ** (1 / 3)
    gamma_x = (x ** 2) / (3 * np.sqrt(1 + x ** 2))
    # equation 8
    drho_dr = - (m * rho) / (gamma_x * r ** 2)
    # equation 9
    dm_dr = (r ** 2) * rho
    return drho_dr, dm_dr

#need to find the radius (range) of integration
def density_range(r, ystate, rho_min=rho_min):
    '''
    Termination function for use in scint.solve_ivp since density will never become negative
    :param r: radius
    :param ystate: state vector for density and mass
    :param rho_min: zero approximation of density
    :return: returns checked ystate for rho and will return 0 when the value of rho gets very small
    '''
    if ystate[0] - rho_min < 0:
        ystate[0] = 0
    return ystate[0]

def solution_system_eq_wd(rho_max, rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx, method='RK45'):
    '''

    :param rho_max: maximum density of the white dwarf star
    :param rad_wd: maximum radius of a white dwarf star
    :param step_size: timestep
    :param rad_min_approx: zero approximation for radius
    :param method: solution method used by scint.solve_ivp (must be string)
    :return: use of scint.solve_ivp with the given values passed.
    '''
    range_density = [rho_max, 0] #range of densities that will be integrated over
    range_radius = [rad_min_approx, rad_wd] #range of radii
    density_range.terminal = True #terminates the solver when 0 is returned by "density_range"
    return scint.solve_ivp(system_eq_wd, range_radius, range_density, events=density_range, method=method)

rho_c = np.logspace(0, 6.39, 10) #used trial and error to get upper limit. this was very close to upper bound given in the question

#initialized storage for results for .t and .y
radius_sol = []
density_sol = []
mass_sol = []
#initialized storage for results for .t_events and .y_events
radius_sol_limit = []
density_sol_limit = []
mass_sol_limit = []

#finding solutions based on the parameters in question 2
for i in range(len(rho_c)):
    wd_sol = solution_system_eq_wd(rho_c[i], rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx)

    rad_sol = wd_sol.t * R_0 / R_sun #changing to solar units
    #y[0] is density and y[1] is mass since system_eq_wd returns drho_dr first and dm_dr second
    den_sol = wd_sol.y[0] * rho_0
    m_sol = wd_sol.y[1] * M_0 / M_sun #changing to solar units

    radius_sol.append(rad_sol)
    density_sol.append(den_sol)
    mass_sol.append(m_sol)

    #this finds the endpoint where radius and mass are maximum and density is minimum
    rad_sol_limit = wd_sol.t_events[0][0] * R_0 / R_sun #changing to solar units
    den_sol_limit = wd_sol.y_events[0][0][0] * rho_0
    m_sol_limit = wd_sol.y_events[0][0][1] * M_0 / M_sun #changing to solar units

    radius_sol_limit.append(rad_sol_limit)
    density_sol_limit.append(den_sol_limit)
    mass_sol_limit.append(m_sol_limit)

#Question 2
rho_legend = []
for i in range(len(rho_c)):
    legend = f'$\\rho$ = {Decimal(rho_c[i]):.2e}'
    rho_legend.append(legend)

fig = plt.figure(figsize=(8,6))

#plotting curve of how mass and radius are related
for i in range(len(density_sol)):
    plt.plot(mass_sol[i], radius_sol[i], ls=':')

plt.title('Mass vs Radius for Various Maximum Internal Densities')
plt.xlabel('Mass (solar masses)')
plt.ylabel('Radius (solar radii)')
plt.legend(rho_legend)

#plotting final mass vs radius value to show the ultimate structure of the simulated star
for i in range(len(density_sol_limit)):
    plt.scatter(mass_sol_limit[i], radius_sol_limit[i])
plt.savefig('westdax_project3_question2_fig')

#checking if the results approach the Chandrasekhar mass limit
M_ch_estimation = []
for i in range(len(mass_sol)):
    M_ch_estimation.append(np.max(mass_sol[i]))
print(f'The estimation of the Chandrasekhar mass limit is: {Decimal(np.max(M_ch_estimation)):.4} Solar Masses')
print(f'The true value of the Chandrasekhar mass limit is: {M_ch} Solar Masses')

#Question 3
rho_c_3 = [rho_c[1], rho_c[4], rho_c[8]] #these three value for initial density were randomly chosen

#initialized storage for results for .t and .y
radius_sol_3 = []
density_sol_3 = []
mass_sol_3 = []
#initialized storage for results for .t_events and .y_events
radius_sol_3_limit = []
density_sol_3_limit = []
mass_sol_3_limit = []

for i in range(len(rho_c_3)):
    #finding solution with alternate solution method
    wd_sol_3 = solution_system_eq_wd(rho_c_3[i], rad_wd=rad_wd_avg, step_size=step, rad_min_approx=zero_approx, method='RK23')

    rad_sol_3 = wd_sol_3.t * R_0 / R_sun #changing to solar units
    den_sol_3 = wd_sol_3.y[0] * rho_0
    m_sol_3 = wd_sol_3.y[1] * M_0 / M_sun #changing to solar units

    radius_sol_3.append(rad_sol_3)
    density_sol_3.append(den_sol_3)
    mass_sol_3.append(m_sol_3)

    # this finds the endpoint where radius and mass are maximum and density is minimum
    rad_sol_3_limit = wd_sol_3.t_events[0][0] * R_0 / R_sun #changing to solar units
    den_sol_3_limit = wd_sol_3.y_events[0][0][0] * rho_0
    m_sol_3_limit = wd_sol_3.y_events[0][0][1] * M_0 / M_sun #changing to solar units

    radius_sol_3_limit.append(rad_sol_3_limit)
    density_sol_3_limit.append(den_sol_3_limit)
    mass_sol_3_limit.append(m_sol_3_limit)

legend_method_comp = [f'$\\rho$ = {Decimal(rho_c[1]):.2e}', f'$\\rho$ = {Decimal(rho_c[4]):.2e}', f'$\\rho$ = {Decimal(rho_c[8]):.2e}']
fig2, (ax1, ax2) = plt.subplots(1,2, figsize=(15,5))

#alternate solution plotting of mass vs radius curve
for i in range(len(density_sol_3)):
    ax1.plot(mass_sol_3[i], radius_sol_3[i], ls=':')

ax1.legend(legend_method_comp)
ax2.legend(legend_method_comp)

#plotting final mass vs radius value to show the ultimate structure of the simulated star
for i in range(len(density_sol_3_limit)):
    ax1.scatter(mass_sol_3_limit[i], radius_sol_3_limit[i])

#plotting same values from original solution method for the sake of comparison
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
plt.savefig('westdax_project3_question3_fig')

#Question 4
#initializing storage of csv file data
wd_mass = []
wd_mass_unc = []
wd_radius = []
wd_radius_unc = []

#opening file and sorting/stripping it to make it usable. This was done by opening it manually in an alternate program and looking at the structure
with open('wd_mass_radius.csv', 'r') as wd_mass_data:
    for line in wd_mass_data.readlines()[1:]:
        cols = line.strip().split(',')
        wd_mass.append(float(cols[0]))
        wd_mass_unc.append(float(cols[1]))
        wd_radius.append(float(cols[2]))
        wd_radius_unc.append(float(cols[3]))

fig3 = plt.figure(figsize=(8,6))

#this is the same data/plotting as question 2
for i in range(len(density_sol)):
    plt.plot(mass_sol[i], radius_sol[i], ls=':')

plt.legend(rho_legend)
plt.xlabel('Mass (solar masses)')
plt.ylabel('Radius (solar radii)')
plt.title('Observation vs Simulation Comparison')

#csv data plotted for comparison with simulation results
for i in range(len(density_sol_limit)):
    plt.scatter(mass_sol_limit[i], radius_sol_limit[i])
plt.scatter(wd_mass, wd_radius)
plt.errorbar(wd_mass, wd_radius, xerr=wd_mass_unc, yerr=wd_radius_unc, color='k', ls='')

plt.savefig('westdax_project3_question4_fig')