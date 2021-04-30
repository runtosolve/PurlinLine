using PurlinLine
using StructuresKit

design_code = "AISI S100-16 ASD"

loading_direction = "gravity"   #or "uplift"

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
                                        ("Z", 0.059, 0.91, 2.5, 8.0, 2.5, 0.91, -50.0, 0.0, 90.0, 0.0, -50.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059)]

material_properties = [(29500.0, 0.30, 55.0, 70.0),
                              (29500.0, 0.30, 55.0, 70.0)]

#type="screw-fastened", thickness, fastener spacing, fastener diameter, fastener_shear_strength fastener_top_flange_location
#type="standing seam", thickness, clip spacing, clip stifness
deck_details = ("screw-fastened", 0.0179, 12.0, 0.212, 2.50, 1.25)

deck_material_properties = (29500.0, 0.30, 55.0, 70.0)

frame_flange_width = 24.0

support_locations = [0.0, 25.0*12, 50.0*12]

bridging_locations =[0.0, 50.0*12]


purlin_line = PurlinLine.define(design_code, loading_direction, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, bridging_locations)


#Write PurlinLine to ThinWalledBeam interface.

#inputs

#Ix Iy Ixy J Cw

num_purlin_sections = size(purlin_line.inputs.cross_section_data)[1]

section_properties = Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}(undef, num_purlin_sections)

for i = 1:num_purlin_sections

    #Note -Ixy here since +y in ThinWalledBeam formulation is pointing down.
    section_properties[i] = (purlin_line.cross_section_data[i].section_properties.Ixx, purlin_line.cross_section_data[i].section_properties.Iyy, -purlin_line.cross_section_data[i].section_properties.Ixy, purlin_line.cross_section_data[i].section_properties.J, purlin_line.cross_section_data[i].section_properties.Cw)

end

num_purlin_materials = size(purlin_line.inputs.material_properties)[1]

material_properties = Vector{Tuple{Float64, Float64}}(undef, num_purlin_materials)

for i = 1:num_purlin_materials

    material_properties[i] = (purlin_line.inputs.material_properties[i][1], purlin_line.inputs.material_properties[i][2])

end

#Calculate the distance from the purlin shear center to the applied load point which is assumed to be the center of the top flange.

num_purlin_sections = size(purlin_line.cross_section_data)[1]

load_location= Vector{Tuple{Float64, Float64}}(undef, num_purlin_sections)

for i = 1:num_purlin_sections

    center_top_flange_node_index = sum(purlin_line.cross_section_data[i].n[1:3]) + sum(purlin_line.cross_section_data[i].n_radius[1:3]) + floor(Int,purlin_line.cross_section_data[i].n[4]/2) + 1

    ax = purlin_line.cross_section_data[i].node_geometry[center_top_flange_node_index, 1] - purlin_line.cross_section_data[1].section_properties.xs

    t = purlin_line.inputs.cross_section_dimensions[i][2]

    ay = (purlin_line.cross_section_data[i].node_geometry[center_top_flange_node_index, 2] + t/2) - purlin_line.cross_section_data[1].section_properties.ys

    load_location[i] = (ax, ay)

end

#Define the lateral and rotational stiffness magnitudes in each purlin segment.
num_purlin_segments = size(purlin_line.bracing_data)[1]

spring_stiffness = Vector{Tuple{Float64, Float64}}(undef, num_purlin_segments)

for i=1:num_purlin_segments

    spring_stiffness[i] = (purlin_line.bracing_data[i].kx, purlin_line.bracing_data[i].kϕ)

end

#Define the y-distance from the purlin shear center to the lateral translational spring.

#Assume ay_kx = ay calculate for the load location.

num_purlin_sections = size(purlin_line.cross_section_data)[1]

spring_location = Vector{Tuple{Float64}}(undef, num_purlin_sections)

for i = 1:num_purlin_sections

    spring_location[i] = (load_location[i][2],)  #check if tuple is needed here

end

#Define the purlin segment properties.
num_purlin_segments = size(purlin_line.bracing_data)[1]

member_definitions = Vector{Tuple{Float64, Float64, Int64, Int64, Int64, Int64, Int64}}(undef, num_purlin_segments)

for i=1:num_purlin_segments

    L = purlin_line.inputs.segments[i][1]
    dL = L / 13
    section_id = purlin_line.inputs.segments[i][2]
    material_id = purlin_line.inputs.segments[i][3]

    #L dL section_properties material_properties load_location spring_stiffness spring_location
    #This is a little tricky. The load_location goes with the purlin cross-section.   The spring_stiffness goes with the purlin_segment.  The spring_location goes with the purlin cross-section.
    member_definitions[i] = (L, dL, section_id, material_id, section_id,i, section_id)

end

#Define purlin line support locations.
#location where u=v=ϕ=0
supports = purlin_line.inputs.support_locations

#Define purlin line end boundary conditions.

end_boundary_conditions = Array{Int64}(undef, 2)

purlin_line_length = sum([purlin_line.inputs.segments[i][1] for i=1:size(purlin_line.inputs.segments)[1]])

#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)

#z=0 (left) end
if supports[1] == 0.0
    end_boundary_conditions[1] = 1 #pin
else
    end_boundary_conditions[1] = 3  #cantilever
end

#z=purlin_line_length (right) end
if supports[end] == purlin_line_length
    end_boundary_conditions[2] = 1
else
    end_boundary_conditions[2] = 3  #cantilever
end


load = (0.0, 0.001) #kips/in.   


z, u, v, ϕ, beam_properties = ThinWalledBeam.solve(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, load, load_location, end_boundary_conditions)

using Plots
plot(z, ϕ)

Mxx = InternalForces.moment(z, beam_properties.dm, -v, beam_properties.E, beamProperties.Ix)

Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)

Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)

B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)





#MAPPING LAYER


#ANALYSIS LAYER

#run a structural 



#RESULTS LAYER

#stiffness provided by deck to purlin
#This needs section properties and material properties.
#Calculate for every unique segment along the span.



#Go through the purlin cross-sections to calculate lateral and rotational stiffness.

#Define the spacing between top flange connections that restrain distortional buckling.
