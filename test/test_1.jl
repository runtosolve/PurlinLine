using PurlinLine
using StructuresKit

design_code = "AISI S100-16 ASD"

# Define the properties of each purlin segment along the line.
                #length, section_properties, material_properties
segments = [(23.0*12, 1, 1),
                (2.0*12, 2, 1),
                (2.0*12, 2, 1),
                (23.0*12, 1, 1)]

# Define the purlin spacing.
spacing = 60;  #in.

# Define the roof slope.
roof_slope = 0.0;   #degrees

#Define the purlin cross-section type.
cross_section_dimensions = [("Z", 0.059, 0.91, 2.5, 8.0, 2.5, 0.91, -50.0, 0.0, 90.0, 0.0, -50.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059),
                                        ("Z", 0.059*2, 0.91, 2.5, 8.0, 2.5, 0.91, -50.0, 0.0, 90.0, 0.0, -50.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059)]

material_properties = [(29500.0, 0.30, 55.0, 70.0),
                              (29500.0, 0.30, 55.0, 70.0)]

#type="screw-fastened", thickness, fastener spacing, fastener diameter, fastener_shear_strength

#type="standing seam", thickness, clip spacing, clip stiffness
deck_details = ("screw-fastened", 0.0179, 12.0, 0.212, 2.50)

deck_material_properties = (29500.0, 0.30, 55.0, 70.0)

frame_flange_width = 24.0

support_locations = [0.0, 25.0*12, 50.0*12]

purlin_frame_connections = "bottom flange connection"
# purlin_frame_connections = "anti-roll clip"

bridging_locations =[0.0, 10.0*12, 50.0*12]

#Calculate purlin line design variables from user inputs and store them in the data structure.
purlin_line = PurlinLine.build(design_code, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, purlin_frame_connections, bridging_locations);
# purlin_line.applied_pressure = 10 / 1000 / 144 #psf to kips/in^2

purlin_line.loading_direction = "gravity"

#Perform a test to collapse.
purlin_line = PurlinLine.test(purlin_line);

purlin_line.applied_pressure*1000*144

using Plots
plot(purlin_line.model.z, purlin_line.internal_forces.T)








