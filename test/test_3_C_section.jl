using PurlinLine
using StructuresKit

design_code = "AISI S100-16 nominal"

# Define the properties of each purlin segment along the line.
                #length, section_properties, material_properties
segments = [(25.0*12, 1, 1)]

# Define the purlin spacing.
spacing = 60;  #in.

# Define the roof slope.
roof_slope = 0.0;   #degrees

#Define the purlin cross-section type.
cross_section_dimensions = [("C", 0.059, 0.91, 2.5, 8.0, 2.5, 0.91, -90.0, 0.0, 90.0, 180.0, -90.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059)]

material_properties = [(29500.0, 0.30, 55.0, 70.0)]

#type="screw-fastened", thickness, fastener spacing, fastener diameter, fastener_shear_strength
deck_details = ("screw-fastened", 0.018, 6.0, 0.212, 2.50)

deck_material_properties = (29500.0, 0.30, 55.0, 70.0)

frame_flange_width = 24.0

support_locations = [0.0, 25.0*12]

bridging_locations =[ ]

#Calculate purlin line design variables from user inputs and store them in the data structure.
purlin_line = PurlinLine.build(design_code, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, bridging_locations)
# purlin_line.applied_pressure = 10 / 1000 / 144 #psf to kips/in^2

purlin_line.cross_section_data[1].node_geometry

using Plots

plot(purlin_line.cross_section_data[1].node_geometry[:,1], purlin_line.cross_section_data[1].node_geometry[:,2], seriestype = :scatter)

plot(purlin_line.free_flange_cross_section_data[1].node_geometry[:,1], purlin_line.free_flange_cross_section_data[1].node_geometry[:,2], seriestype = :scatter)


purlin_line.loading_direction = "gravity"

#Perform a test to collapse.
purlin_line = PurlinLine.test(purlin_line)

purlin_line.applied_pressure * 1000 * 144

purlin_line.cross_section_data[1].plastic_section_properties

using Plots
plot(purlin_line.model.z, purlin_line.free_flange_model.u)

plot(purlin_line.model.z, purlin_line.free_flange_model.u)

plot(purlin_line.model.z, purlin_line.flexure_torsion_demand_to_capacity.action_Mxx)
plot!(purlin_line.model.z, purlin_line.flexure_torsion_demand_to_capacity.action_Myy)
plot!(purlin_line.model.z, purlin_line.flexure_torsion_demand_to_capacity.action_B)
plot!(purlin_line.model.z, purlin_line.flexure_torsion_demand_to_capacity.action_Myy_freeflange)

purlin_line.failure_limit_state

purlin_line.failure_location

plot(purlin_line.model.z, purlin_line.free_flange_internal_forces.P)

plot(purlin_line.model.z, purlin_line.free_flange_model.kx)

calculate_free_flange_shear_flow_properties(purlin_line)


#' Plot the purlin local buckling mode shape.
using StructuresKit
mode_index = 2
scale_x = 1.0
scale_y = 1.0
node = purlin_line.local_buckling_xx_pos[1].CUFSM_data.node
elem = purlin_line.local_buckling_xx_pos[1].CUFSM_data.elem
shapes = purlin_line.local_buckling_xx_pos[1].CUFSM_data.shapes

StructuresKit.CUFSM.view_multi_branch_section_mode_shape(node, elem, shapes, mode_index, scale_x, scale_y)


#' Plot the purlin local buckling signature curve.
num_lengths = length(purlin_line.local_buckling_xx_pos[1].CUFSM_data.lengths)
half_wavelengths = [purlin_line.local_buckling_xx_pos[1].CUFSM_data.curve[i,1][1] for i=1:num_lengths]
load_factors = [purlin_line.local_buckling_xx_pos[1].CUFSM_data.curve[i,1][2] for i=1:num_lengths]
plot(half_wavelengths, load_factors, legend = false)


#' Plot the purlin local buckling mode shape.
using StructuresKit
mode_index = 3
scale_x = 1.0
scale_y = 1.0
node = purlin_line.distortional_buckling_xx_pos[1].CUFSM_data.node
elem = purlin_line.distortional_buckling_xx_pos[1].CUFSM_data.elem
shapes = purlin_line.distortional_buckling_xx_pos[1].CUFSM_data.shapes

StructuresKit.CUFSM.view_multi_branch_section_mode_shape(node, elem, shapes, mode_index, scale_x, scale_y)



#' Plot the purlin local buckling signature curve.
num_lengths = length(purlin_line.distortional_buckling_xx_pos[1].CUFSM_data.lengths)
half_wavelengths = [purlin_line.distortional_buckling_xx_pos[1].CUFSM_data.curve[i,1][1] for i=1:num_lengths]
load_factors = [purlin_line.distortional_buckling_xx_pos[1].CUFSM_data.curve[i,1][2] for i=1:num_lengths]
plot(half_wavelengths, load_factors, legend = false)



