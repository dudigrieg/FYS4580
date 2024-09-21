import openmc
import matplotlib.pyplot as plt
import numpy

""" SECTION 2 exercise E """


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