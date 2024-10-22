import openmc
import numpy as np
import matplotlib.pyplot as plt
import os

# create urox (UO2)
urox = openmc.Material(name = "Uranium dioxide")
urox.add_nuclide("U235", 0.02)
urox.add_nuclide("U238", 0.98)
urox.add_nuclide("O16", 2.0)
urox.set_density("g/cm3",10.97)
urox.add_s_alpha_beta("c_O_in_UO2")

# create light water (H20)
light_water = openmc.Material(name = "Light Water")
light_water.add_nuclide('H1', 2)
light_water.add_nuclide('O16', 1)
# light_water.add_element("H", 2.0)
# light_water.add_element("O", 1.0)
light_water.set_density("g/cm3", 0.997)
light_water.add_s_alpha_beta('c_H_in_H2O')

# create Helium (He)
helium = openmc.Material(name = "Helium")
helium.add_element("He", 1.0)
helium.set_density("g/cm3",0.1786)


# create zircaloy
zircaloy = openmc.Material(name= "Zircaloy")
zircaloy.add_element('Zr', 0.9825)
zircaloy.add_element('Sn', 0.0145)
zircaloy.add_element('Fe', 0.0015)
zircaloy.add_element('Cr', 0.0015)
zircaloy.set_density('g/cm3', 6.56)

''' create heavy water '''

# create heavy water
heavy_water = openmc.Material(name = "Heavy Water")
heavy_water.add_nuclide('H2', 2)
heavy_water.add_nuclide('O16', 1)
heavy_water.set_density("g/cm3", 1.1)
heavy_water.add_s_alpha_beta('c_D_in_D2O')

'''.'''


''' create Beryllium '''

beryllium = openmc.Material(name="Beryllium")
beryllium.add_nuclide('Be9', 1.0)
beryllium.set_density("g/cm3", 1.85)
beryllium.add_s_alpha_beta('c_Be_in_BeO')

'''.'''



# Create box with cylinder inside
box = openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type = "vacuum")
#box = openmc.model.RectangularPrism(width=4.2, height=4.2, boundary_type = "vacuum")

z_max = openmc.ZPlane(1, boundary_type = "vacuum")
z_min = openmc.ZPlane(-1, boundary_type = "vacuum")

# create the three cylinders inside the box
fuel_cylinder = -openmc.ZCylinder(r=0.4, boundary_type = "transmission")
helium_gap = -openmc.ZCylinder(r = 0.42, boundary_type = "transmission")
cladding_region = -openmc.ZCylinder(r = 0.45, boundary_type = "transmission")

# fuel_cylinder = -openmc.ZCylinder(r=1.2, boundary_type = "transmission")
# helium_gap = -openmc.ZCylinder(r = 1.26, boundary_type = "transmission")
# cladding_region = -openmc.ZCylinder(r = 1.35, boundary_type = "transmission")

fuel_cylinder = fuel_cylinder &+z_min &-z_max

# create box that doesn't include area of the cylinders
box_outside_cylinder = -box &+ z_min &- z_max &~ cladding_region

# make sure the cylinders don't overlap
helium_gap = helium_gap &+z_min &-z_max &~fuel_cylinder
cladding_region = cladding_region &+z_min &-z_max &~helium_gap &~fuel_cylinder


# create cells
pin_universe = openmc.Universe(name = "Pincell")
fuel_cell = openmc.Cell(fill = urox, region = fuel_cylinder)
helium_cell = openmc.Cell(fill = helium, region = helium_gap)
cladding = openmc.Cell(fill = zircaloy, region = cladding_region)
moderator = openmc.Cell(fill = light_water, region = box_outside_cylinder)



pin_universe.add_cells([fuel_cell, helium_cell, cladding, moderator])


geometry = openmc.Geometry(pin_universe)
geometry.root_universe = pin_universe #What universe do we run?
geometry.export_to_xml()

pin_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
    colors = {urox: "green", helium: "brown",
            zircaloy: "gray", light_water : "blue"})
#plt.show()
#plt.savefig(fname = "My_Pincell")




# Run the simulation

cross_sections_file = "/Users/claudiagrieg/Documents/MASTER/FYS4580/endfb-vii.1-hdf5/cross_sections.xml"
os.environ["OPENMC_CROSS_SECTIONS"] = cross_sections_file


settings = openmc.Settings()
settings.batches = 100 #How many iterations of neutron generations
settings.inactive = 10 #How many generations used to determine sample spots
settings.particles = 10000 #Number of neutrons in each generation

materials = openmc.Materials([urox, helium, zircaloy, light_water, heavy_water])
materials.export_to_xml()

#The following section tells OpenMC where to sample particle sites.
bounds = [-0.7,-0.7,-1,0.7,0.7,1]
#bounds = [-1.4, -1.4, -1, 1.4, 1.4, 1]
# bounds = [-2.1, -2.1, -1, 2.1, 2.1, 1]

uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()


mesh = openmc.RegularMesh() # Making a mesh instance
mesh.dimension = [100,100] # Setting the bins of the mesh to 100 x 100 bins
mesh.lower_left = [-0.7,-0.7] # Setting lower left coordinates
mesh.upper_right = [0.7,0.7] # Setting upper right coordinates
# mesh.lower_left = [-1.4, -1.4]
# mesh.upper_right = [1.4, 1.4]
# mesh.lower_left = [-2.1, -2.1]
# mesh.upper_right = [2.1, 2.1]
mesh_filter = openmc.MeshFilter(mesh) # Making a meshfilter
tallies_file = openmc.Tallies() # For our tallies to be included we have to put them into an inp
tally = openmc.Tally(name = "Tally") # Making a tally instance
tally.filters = [mesh_filter] # Designatingen the mesh for the tally
tally.scores = ["flux", 'prompt-nu-fission'] # Designating the objects to be tallied
tallies_file.append(tally) # Appending the tally to our tallies file
tallies_file.export_to_xml()

if __name__ == "__main__":
    openmc.run(threads = 4) #Threads set for heplab primarily