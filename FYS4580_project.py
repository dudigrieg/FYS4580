import openmc
import numpy as np
import matplotlib.pyplot as plt
import os

''' PART 1 '''

## EXERCISE A

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
# cross section?

# create zircaloy
zircaloy = openmc.Material(name= "Zircaloy")
zircaloy.add_element('Zr', 0.9825)
zircaloy.add_element('Sn', 0.0145)
zircaloy.add_element('Fe', 0.0015)
zircaloy.add_element('Cr', 0.0015)
zircaloy.set_density('g/cm3', 6.56)
# cross section ?

# create heavy water
heavy_water = openmc.Material(name = "Heavy Water")
heavy_water.add_nuclide('H2', 2)
heavy_water.add_nuclide('O16', 1)
heavy_water.set_density("g/cm3", 1.1)
heavy_water.add_s_alpha_beta('c_D_in_D2O')

# create graphite
graphite = openmc.Material(name='Graphite')
graphite.add_nuclide('C0', 1.0)
graphite.set_density('g/cm3', 1.7)

# # PART 2: Define control rod material
# boron_carbide = openmc.Material(name="Boron Carbide (B4C)")
# boron_carbide.add_element('B', 4.0)
# boron_carbide.add_nuclide('C0', 1.0)
# boron_carbide.set_density('g/cm3', 2.52)



## EXERCISE B

# Create box with cylinder inside
box = openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type = "reflective")

height = 300 #cm
z_max = openmc.ZPlane(1, boundary_type = "reflective")
z_min = openmc.ZPlane(-1, boundary_type = "reflective")

# create the three cylinders inside the box
fuel_cylinder = -openmc.ZCylinder(r=0.4, boundary_type = "transmission")
helium_gap = -openmc.ZCylinder(r = 0.42, boundary_type = "transmission")
cladding_region = -openmc.ZCylinder(r = 0.45, boundary_type = "transmission")

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
moderator = openmc.Cell(fill = heavy_water, region = box_outside_cylinder)


pin_universe.add_cells([fuel_cell, helium_cell, cladding, moderator])

# # PART 2: Define control rod geometry (replace fuel with boron carbide)
controlpin_cell = openmc.Cell(fill=boron_carbide, region=fuel_cylinder)
controlpin_gap_cell = openmc.Cell(fill=helium, region=helium_gap)
controlpin_clad_cell = openmc.Cell(fill=zircaloy, region=cladding_region)
controlpin_moderator_cell = openmc.Cell(fill=heavy_water, region=box_outside_cylinder)

controlpin_universe = openmc.Universe(name = "Control Pin")
controlpin_universe.add_cells([controlpin_cell, controlpin_gap_cell, controlpin_clad_cell, controlpin_moderator_cell])

pinx = 1.4
piny = 1.4
pinz = 2
assembly_x = 17
assembly_y = 17
assembly_z = 150

universii = np.zeros((assembly_z,assembly_x,assembly_y), dtype = openmc.Universe)

#Making array for universes
for z in range(0,assembly_z):
    for i in range (0,assembly_x):
        for j in range(0,assembly_y):
            if i in [x,x,x] and j in [x,x,x]:
                universii[z][i][j] = c
            else:
                universii[z][i][j] = f
#Decide placement of controlpins
#Place controlpins
#Place fuel pins

#making the lattice
f_assembly = openmc.RectLattice(name = "fuel_assembly")
f_assembly.lower_left = (-assembly_x*abs(pinx)/2,-assembly_y*abs(piny)/2,-assembly_z)
f_assembly.universes = universii
f_assembly.pitch = (1.4,1.4,2)

#making outer universe
water_cell = openmc.Cell(fill = water)
outer_uni = openmc.Universe(cells = [water_cell])
f_assembly.outer = outer_uni

#Putting the lattice into a simulation
outedge = "reflective" #Infinite reactor or singular assembly
outcap= openmc.ZPlane(assembly_z, boundary_type = outedge)
outbot = openmc.ZPlane(-assembly_z, boundary_type = outedge)
outcell = openmc.model.rectangular_prism(width = 23.8, height = 23.8, boundary_type = outedge)
outcell = outcell &+outbot &-outcap

f_cell = openmc.Cell(fill = f_assembly, region = outcell) #making cell for assembly
cell_list = [f_cell]
assembly_uni = openmc.Universe(cells = cell_list) #Putting the assembly into a universe





geometry = openmc.Geometry()
geometry.root_universe = assembly_uni
geometry.root_universe = pin_universe #What universe do we run?
geometry.export_to_xml()

pin_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
    colors = {urox: "green", helium: "brown",
            zircaloy: "gray", heavy_water : "blue"})
plt.show()
#plt.savefig(fname = "My_Pincell")



## EXERCISE C + D

# # Run the simulation

''' changed from vacuum to reflective. what happens to k?'''

''' ============================>     RESULTS  VACUUM   <============================

    k-effective (Collision)     = 0.00956 +/- 0.00005
 k-effective (Track-length)  = 0.00951 +/- 0.00001
 k-effective (Absorption)    = 0.00988 +/- 0.00016
 Combined k-effective        = 0.00950 +/- 0.00001
 Leakage Fraction            = 0.99549 +/- 0.00007 '''

'''  ============================>     RESULTS  REFLECTIVE   <============================

 k-effective (Collision)     = 1.29520 +/- 0.00136
 k-effective (Track-length)  = 1.29373 +/- 0.00148
 k-effective (Absorption)    = 1.29335 +/- 0.00102
 Combined k-effective        = 1.29349 +/- 0.00089
 Leakage Fraction            = 0.00000 +/- 0.00000 '''

cross_sections_file = "/Users/claudiagrieg/Documents/MASTER/FYS4580/endfb-vii.1-hdf5/cross_sections.xml"
os.environ["OPENMC_CROSS_SECTIONS"] = cross_sections_file


settings = openmc.Settings()
settings.batches = 100 #How many iterations of neutron generations
settings.inactive = 10 #How many generations used to determine sample spots
settings.particles = 10000 #Number of neutrons in each generation

materials = openmc.Materials([light_water, urox, helium, zircaloy, heavy_water, graphite])
materials.export_to_xml()

#The following section tells OpenMC where to sample particle sites.
bounds = [-0.7,-0.7,-1,0.7,0.7,1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()


## EXERCISE E

''' in file called post.py '''

mesh = openmc.RegularMesh() # Making a mesh instance
mesh.dimension = [100,100] # Setting the bins of the mesh to 100 x 100 bins
mesh.lower_left = [-0.7,-0.7] # Setting lower left coordinates
mesh.upper_right = [0.7,0.7] # Setting upper right coordinates
mesh_filter = openmc.MeshFilter(mesh) # Making a meshfilter
tallies_file = openmc.Tallies() # For our tallies to be included we have to put them into an inp
tally = openmc.Tally(name = "Tally") # Making a tally instance
tally.filters = [mesh_filter] # Designatingen the mesh for the tally
tally.scores = ["flux", 'prompt-nu-fission'] # Designating the objects to be tallied
tallies_file.append(tally) # Appending the tally to our tallies file



## EXERCISE F

''' in file called post.py'''

erange = np.logspace(-5, 7.5, 501)
energy_filter = openmc.EnergyFilter(erange) # create energy filter

tally2 = openmc.Tally(name = "Energy") #Making a new tally only looking at flux
tally2.filters = [energy_filter] #Designating the energy filter
tally2.scores = ["flux"]
tallies_file.append(tally2)
#tallies_file.export_to_xml()


## EXERCISE G 

''' exchanged light water with heavy water and graphite. What happens to k? settings: reflective boundary, default size. 
light water:
 k-effective (Collision)     = 1.29298 +/- 0.00133
 k-effective (Track-length)  = 1.29138 +/- 0.00143
 k-effective (Absorption)    = 1.29339 +/- 0.00101
 Combined k-effective        = 1.29295 +/- 0.00087.

heavy_water: 
 k-effective (Collision)     = 0.86686 +/- 0.00111
 k-effective (Track-length)  = 0.86796 +/- 0.00129
 k-effective (Absorption)    = 0.86889 +/- 0.00097
 Combined k-effective        = 0.86787 +/- 0.00082

graphite:
 k-effective (Collision)     = 0.52479 +/- 0.00072
 k-effective (Track-length)  = 0.52483 +/- 0.00071
 k-effective (Absorption)    = 0.52423 +/- 0.00083
 Combined k-effective        = 0.52459 +/- 0.00059
 
changed size to 2x bigger:
light water:
 k-effective (Collision)     = 1.28187 +/- 0.00140
 k-effective (Track-length)  = 1.28111 +/- 0.00178
 k-effective (Absorption)    = 1.28293 +/- 0.00106
 Combined k-effective        = 1.28251 +/- 0.00095

heavy water:
 k-effective (Collision)     = 0.92415 +/- 0.00122
 k-effective (Track-length)  = 0.92544 +/- 0.00132
 k-effective (Absorption)    = 0.92392 +/- 0.00099
 Combined k-effective        = 0.92426 +/- 0.00088

graphite:
 k-effective (Collision)     = 0.56519 +/- 0.00078
 k-effective (Track-length)  = 0.56505 +/- 0.00081
 k-effective (Absorption)    = 0.56569 +/- 0.00112
 Combined k-effective        = 0.56529 +/- 0.00074

changed size to 4x bigger:
light water:
 k-effective (Collision)     = 1.19951 +/- 0.00128
 k-effective (Track-length)  = 1.20008 +/- 0.00155
 k-effective (Absorption)    = 1.20117 +/- 0.00101
 Combined k-effective        = 1.20079 +/- 0.00092

heavy water:
 k-effective (Collision)     = 1.00246 +/- 0.00128
 k-effective (Track-length)  = 1.00210 +/- 0.00141
 k-effective (Absorption)    = 1.00175 +/- 0.00100
 Combined k-effective        = 1.00197 +/- 0.00093

graphite:
 k-effective (Collision)     = 0.62165 +/- 0.00095
 k-effective (Track-length)  = 0.62163 +/- 0.00103
 k-effective (Absorption)    = 0.62027 +/- 0.00102
 Combined k-effective        = 0.62102 +/- 0.00077
'''

## EXERCISE H

''' based on what i saw in exercise G, i used heavy water as moderator material with cylinder size 3.5x the original size, and increased the fuel enrichment to 3.5%, which lead to a Combined k-effective = 1.00322 +/- 0.00097. I finetuned this by setting the enrichment to 3.45% and got k = 1.00043 +/- 0.00091. Variables in four factor formula changed:'''

## EXERCISE I

## calculate thermal utilization factor f in post processing file

# fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
# fuel_therm_abs_rate.scores = ['absorption']
# fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625e-6]),
#                                openmc.CellFilter([fuel_cell.id])]
# tallies_file.append(fuel_therm_abs_rate)

# moderator_therm_abs_rate = openmc.Tally(name='moderator therm. abs. rate')
# moderator_therm_abs_rate.scores = ['absorption']
# moderator_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625e-6]),
#                                     openmc.CellFilter([moderator.id])]
# tallies_file.append(moderator_therm_abs_rate)

# cladding_therm_abs_rate = openmc.Tally(name='cladding therm. abs. rate')
# cladding_therm_abs_rate.scores = ['absorption']
# cladding_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625e-6]),
#                                    openmc.CellFilter([cladding.id])]
# tallies_file.append(cladding_therm_abs_rate)

tallies_file.export_to_xml()

if __name__ == "__main__":
    openmc.run(threads=4)



''' PART TWO '''

## EXERCISE A

# define control rod material, define control rod geometry