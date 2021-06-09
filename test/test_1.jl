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
                                        ("Z", 0.049, 0.91, 2.5, 8.0, 2.5, 0.91, -50.0, 0.0, 90.0, 0.0, -50.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059)]

material_properties = [(29500.0, 0.30, 55.0, 70.0),
                              (29500.0, 0.30, 55.0, 70.0)]

#type="screw-fastened", thickness, fastener spacing, fastener diameter, fastener_shear_strength

#type="standing seam", thickness, clip spacing, clip stiffness
deck_details = ("screw-fastened", 0.0179, 12.0, 0.212, 2.50)

deck_material_properties = (29500.0, 0.30, 55.0, 70.0)

frame_flange_width = 24.0

support_locations = [0.0, 25.0*12, 50.0*12]

bridging_locations =[0.0, 10.0*12, 50.0*12]

#Add inputs to data structure.
purlin_line = PurlinLine.define(design_code, loading_direction, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, bridging_locations)

#Prepare ThinWalledBeam inputs.
section_properties, material_properties, spring_stiffness, spring_location, supports, end_boundary_conditions = PurlinLine.thin_walled_beam_interface(purlin_line)

#Work on loads tomorrow!


#Define loads.
dz, z, dm = Mesh.define_line_element(member_definitions)
num_nodes = length(z)
#qx qy
loads = [(0.001 * ones(num_nodes)), (0.001 * ones(num_nodes))] #kips/in.  
#ax ay
load_location = []

#Calculate the distance from the purlin shear center to the applied load point which is assumed to be the center of the top flange.

num_purlin_sections = size(purlin_line.inputs.cross_section_dimensions)[1]

ax_purlin_section = Vector{Float64}(undef, num_purlin_sections)
ay_purlin_section = Vector{Float64}(undef, num_purlin_sections)

for i = 1:num_purlin_sections

    center_top_flange_node_index = sum(purlin_line.cross_section_data[i].n[1:3]) + sum(purlin_line.cross_section_data[i].n_radius[1:3]) + floor(Int,purlin_line.cross_section_data[i].n[4]/2) + 1

    ax_purlin_section[i] = purlin_line.cross_section_data[i].node_geometry[center_top_flange_node_index, 1] - purlin_line.cross_section_data[1].section_properties.xs

    t = purlin_line.inputs.cross_section_dimensions[i][2]

    ay_purlin_section[i] = (purlin_line.cross_section_data[i].node_geometry[center_top_flange_node_index, 2] + t/2) - purlin_line.cross_section_data[1].section_properties.ys

    
end

dz, z, dm = Mesh.define_line_element(purlin_line.inputs.member_definitions)
ax = Mesh.create_line_element_property_array(purlin_line.model.member_definitions, purlin_line.model.dm, dz, ax_purlin_section, 3, 1)
ay = Mesh.create_line_element_property_array(purlin_line.model.member_definitions, purlin_line.model.dm, dz, ay_purlin_section, 3, 1)

load_location = [ax, ay]



















#Define the ThinWalledBeam model.
model = ThinWalledBeam.define(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, loads, load_location, end_boundary_conditions)

#Solve the ThinWalledBeam model.
model = ThinWalledBeam.solve(model)

#Add second order solution to PurlinLine data structure.
purlin_line.model = model




using Plots
plot(purlin_line.model.z, purlin_line.model.ϕ)

Mxx = InternalForces.moment(purlin_line.model.z, purlin_line.model.dm, -purlin_line.model.v, purlin_line.model.E, purlin_line.model.Ix)
Myy = InternalForces.moment(purlin_line.model.z, purlin_line.model.dm, -purlin_line.model.u, purlin_line.model.E, purlin_line.model.Iy)
Vyy = InternalForces.shear(purlin_line.model.z, purlin_line.model.dm, -purlin_line.model.v, purlin_line.model.E, purlin_line.model.Ix)
B = InternalForces.bimoment(purlin_line.model.z, purlin_line.model.dm, purlin_line.model.ϕ, purlin_line.model.E, purlin_line.model.Cw)

P_free_flange = calculate_free_flange_axial_force(Mxx, purlin_line)

function calculate_free_flange_axial_force(Mxx, purlin_line)

    dz = diff(purlin_line.model.z)

    #Define base metal thickness along the purlin line.
    t = Mesh.create_line_element_property_array(purlin_line.model.member_definitions, purlin_line.model.dm, dz, purlin_line.inputs.cross_section_dimensions, 3, 2)

    #Define out-to-out purlin web depth.
    H = Mesh.create_line_element_property_array(purlin_line.model.member_definitions, purlin_line.model.dm, dz, purlin_line.inputs.cross_section_dimensions, 3, 5)

    #Define the centroid location of the free flange from the bottom of the purlin.
    num_purlin_sections = size(purlin_line.inputs.cross_section_dimensions)[1]
    bottom_flange_centroid_location = [purlin_line.free_flange_cross_section_data[i].section_properties.yc for i=1:num_purlin_sections]
    ycf = Mesh.create_line_element_property_array(purlin_line.model.member_definitions, purlin_line.model.dm, dz, bottom_flange_centroid_location, 3, 1)

    #Approximate the axial force in free flange by assuming a force couple between top and bottom flanges.   Assume for now that the purlin has the same top flange and bottom flange centroids.
    P = -Mxx ./ (H .- 2 * abs.(ycf))
   
    return P

end

function beam_column_interface(purlin_line)

    dz = diff(purlin_line.model.z)

    #Define the number of purlin cross-sections.
    num_purlin_sections = size(purlin_line.inputs.cross_section_dimensions)[1]

    #Initialize an array of tuples to hold the free flange section properties.
    section_properties = Vector{Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,}}(undef, num_purlin_sections)

    for i = 1:num_purlin_sections

        Af = purlin_line.free_flange_cross_section_data[i].section_properties.A
        Ixf = purlin_line.free_flange_cross_section_data[i].section_properties.Ixx
        Iyf = purlin_line.free_flange_cross_section_data[i].section_properties.Iyy
        Jf = purlin_line.free_flange_cross_section_data[i].section_properties.J
        Cwf = purlin_line.free_flange_cross_section_data[i].section_properties.Cw
        xcf = purlin_line.free_flange_cross_section_data[i].section_properties.xc
        ycf = purlin_line.free_flange_cross_section_data[i].section_properties.yc
        xsf = purlin_line.free_flange_cross_section_data[i].section_properties.xs
        ysf = purlin_line.free_flange_cross_section_data[i].section_properties.ys

        section_properties[i] = (Af, Ixf, Iyf, Jf, Cwf, xcf, ycf, xsf, ysf)

    end


    #Define the purlin segment properties.
    num_purlin_segments = size(purlin_line.bracing_data)[1]

    member_definitions = Vector{Tuple{Float64, Float64, Int64, Int64}}(undef, num_purlin_segments)

    for i=1:num_purlin_segments

        L = purlin_line.inputs.segments[i][1]
        dL = L / 13
        section_id = purlin_line.inputs.segments[i][2]
        material_id = purlin_line.inputs.segments[i][3]

        #L dL section_properties material_properties
        member_definitions[i] = (L, dL, section_id, material_id)

    end

    #Define purlin line end boundary conditions.

    end_boundary_conditions = Array{Int64}(undef, 2)

    purlin_line_length = sum([purlin_line.inputs.segments[i][1] for i=1:size(purlin_line.inputs.segments)[1]])

    #type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)

    #z=0 (left) end
    if purlin_line.inputs.support_locations[1] == 0.0
        end_boundary_conditions[1] = 1 #pin
    else
        end_boundary_conditions[1] = 3  #cantilever
    end

    #z=purlin_line_length (right) end
    if purlin_line.inputs.support_locations[end] == purlin_line_length
        end_boundary_conditions[2] = 1
    else
        end_boundary_conditions[2] = 3  #cantilever
    end

    #Define supports.   Combine frame supports and intermediate bridging here.
    supports = sort(unique([purlin_line.inputs.support_locations; purlin_line.inputs.bridging_locations]))


    #Define springs.  These are assigned differently than ThinWalledBeam.  I guess kind of an experiment to see which style works better in the long run.   The spring info is assigned along the beam discretization, node by node.

    #Define kxf along the purlin line.
    free_flange_lateral_spring = [purlin_line.free_flange_data[i].kxf for i=1:num_purlin_segments]
    kxf = Mesh.create_line_element_property_array(member_definitions, purlin_line.model.dm, dz, free_flange_lateral_spring, 3, 1)



    #Define kϕf along the purlin line.
    free_flange_rotational_spring = [purlin_line.free_flange_data[i].kϕf for i=1:num_purlin_segments]
    kϕf = Mesh.create_line_element_property_array(member_definitions, purlin_line.model.dm, dz, free_flange_rotational_spring, 3, 1)    

    #Assume the lateral spring acts at the free flange centroid.  This means hy = 0.
    num_nodes = length(purlin_line.model.z)
    #kx ky kϕ hx hy
    springs = [(kxf), (0.0 .* ones(num_nodes)), (kϕf), (0.0 .* ones(num_nodes)), (0.0 .* ones(num_nodes))]


    #Define the purlin material properties.  Just E and ν needed here.
    num_purlin_materials = size(purlin_line.inputs.material_properties)[1]

    material_properties = Vector{Tuple{Float64, Float64}}(undef, num_purlin_materials)

    for i = 1:num_purlin_materials

        material_properties[i] = (purlin_line.inputs.material_properties[i][1], purlin_line.inputs.material_properties[i][2])

    end


    return member_definitions, section_properties, material_properties, springs, end_boundary_conditions, supports

end


member_definitions, section_properties, material_properties, springs, end_boundary_conditions, supports = beam_column_interface(purlin_line)


#Define shear flow force in free flange.

    #Define the purlin segment properties.
    num_purlin_segments = size(purlin_line.bracing_data)[1]

    num_nodes = length(purlin_line.model.z)

    dz = diff(purlin_line.model.z)

    free_flange_kH = [purlin_line.free_flange_data[i].kH for i=1:num_purlin_segments]
    kH = Mesh.create_line_element_property_array(member_definitions, purlin_line.model.dm, dz, free_flange_kH, 3, 1)

        #Define the number of purlin cross-sections.
        num_purlin_sections = size(purlin_line.inputs.cross_section_dimensions)[1]

    free_flange_yc = [purlin_line.free_flange_cross_section_data[i].section_properties.yc for i=1:num_purlin_sections]
    ycf = Mesh.create_line_element_property_array(member_definitions, purlin_line.model.dm, dz, free_flange_yc, 3, 1)

    q = 0.001

    qx = q .* kH

loads = [P_free_flange, qx, (0.0*ones(num_nodes)),(0.0*ones(num_nodes)), ycf]


model = BeamColumn.define(member_definitions, section_properties, material_properties, springs, supports, loads, end_boundary_conditions)


    u, v, ϕ, properties = BeamColumn.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports)







#MAPPING LAYER


#ANALYSIS LAYER

#run a structural 



#RESULTS LAYER

#stiffness provided by deck to purlin
#This needs section properties and material properties.
#Calculate for every unique segment along the span.



#Go through the purlin cross-sections to calculate lateral and rotational stiffness.

#Define the spacing between top flange connections that restrain distortional buckling.
