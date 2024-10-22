import openmc
import matplotlib.pyplot as plt
import numpy as np
from FYS4580_project import erange


""" PART 1 exercise E """

sp = openmc.StatePoint('statepoint.100.h5')
# available_tallies = sp.tallies
# print("Available tallys:", list(available_tallies.keys()))
tally = sp.tallies[1]
flux = tally.get_slice(scores=['flux'])
prompt = tally.get_slice(scores=['prompt-nu-fission'])
flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
prompt.std_dev.shape = (100, 100)
prompt.mean.shape = (100, 100)
fig = plt.subplot(121, title = "Neutron flux distribution")
fig.imshow(flux.mean)
fig2 = plt.subplot(122, title = "Fission site heatmap")
fig2.imshow(prompt.mean)
plt.show()

''' PART 1 EXERCISE F '''

tally2 = sp.tallies[2]
flx = tally2.mean.ravel()
plt.loglog(erange[:-1], flx)
plt.grid()
plt.xlabel("Energy eV")
plt.ylabel("Flux [n/cm-src]")
plt.title("Neutron energy spectrum")
plt.show()


''' PART 1 EXERCISE I '''

# fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
# moderator_therm_abs_rate = sp.get_tally(name='moderator therm. abs. rate')
# cladding_therm_abs_rate = sp.get_tally(name='cladding therm. abs. rate')

# # absorption rates
# fuel_abs_rate = fuel_therm_abs_rate.mean.sum()
# moderator_abs_rate = moderator_therm_abs_rate.mean.sum()
# cladding_abs_rate = cladding_therm_abs_rate.mean.sum()

# # Calculating total non-fuel thermal absorption
# non_fuel_abs_rate = moderator_abs_rate + cladding_abs_rate

# # Calculating thermal utilization factor
# thermal_utilization_factor = fuel_abs_rate / (fuel_abs_rate + non_fuel_abs_rate)

# print(f"Thermal Utilization Factor f:", thermal_utilization_factor)