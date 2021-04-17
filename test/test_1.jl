using PurlinLine
using StructuresKit

#INPUT LAYER


mutable struct purlin_line_object

    design_code::String
    loading_direction::String
    segments::Vector{Tuple{Float64, Int64, Int64}}
    spacing::Float64
    roof_slope::Float64
    cross_section_dimensions::Vector{Tuple{String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}
    material_properties::Vector{NTuple{4, Float64}}
    deck_details::Tuple{String, Float64, Float64, Float64, Float64, Float64}
    deck_material_properties::NTuple{4, Float64}
    frame_flange_width::Float64
    support_locations::Vector{Float64}
    bridging_locations::Vector{Float64}

    section_properties::Vector{StructuresKit.CrossSection.SectionProperties}
    free_flange_section_properties::Vector{StructuresKit.CrossSection.SectionProperties}

    cross_section_node_geometry::Vector{Array{Float64}}
    cross_section_element_connectivity::Vector{Array{Float64}}

    kp::Vector{Float64}
    kϕ::Vector{Float64}
    kϕ_dist::Vector{Float64}
    L_crd::Vector{Float64}
    Lm::Float64
    kx::Vector{Float64}

    purlin_line_object() = new()

end


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

#MAPPING LAYER

#Create the data structure.
purlin_line = purlin_line_object()

#Map user inputs to the data structure.
purlin_line.design_code = design_code
purlin_line.loading_direction = loading_direction
purlin_line.segments = segments
purlin_line.spacing = spacing
purlin_line.roof_slope = roof_slope
purlin_line.cross_section_dimensions = cross_section_dimensions
purlin_line.material_properties = material_properties
purlin_line.deck_details = deck_details
purlin_line.deck_material_properties = deck_material_properties
purlin_line.frame_flange_width = frame_flange_width
purlin_line.support_locations = support_locations
purlin_line.bridging_locations = bridging_locations

#DEFINITIONS LAYER

##Calculate properties associated with each purlin line segment.

#Calculate the purlin section properties.
purlin_line.section_properties, purlin_line.free_flange_section_properties, purlin_line.cross_section_node_geometry, purlin_line.cross_section_element_connectivity = PurlinLine.calculate_purlin_section_properties(purlin_line.cross_section_dimensions)

purlin_line.kp, purlin_line.kϕ, purlin_line.kϕ_dist, purlin_line.kx, purlin_line.Lm = PurlinLine.define_deck_bracing_properties(purlin_line)


P = 1.0
Mxx = 0.0
Mzz = 0.0
M11 = 0.0
M22 = 0.0
A = purlin_line.section_properties[1].A
xcg = purlin_line.section_properties[1].xc
zcg = purlin_line.section_properties[1].yc
Ixx = purlin_line.section_properties[1].Ixx
Izz = purlin_line.section_properties[1].Iyy
Ixz = purlin_line.section_properties[1].Ixy
thetap = rad2deg(purlin_line.section_properties[1].θ)
I11 = purlin_line.section_properties[1].I1
I22 = purlin_line.section_properties[1].I2
unsymm = 1

num_cross_section_nodes = size(purlin_line.cross_section_node_geometry[1])[1]
node = zeros(Float64, (num_cross_section_nodes, 8))
node[:, 2:3] .= purlin_line.cross_section_node_geometry[1]

node = CUFSM.stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

#Calculate
#local buckling
#calculate along with cross-section

#bracing stiffnesses
#distortional buckling
#flexural strength
#shear strength
#web crippling strength

#ANALYSIS LAYER

#RESULTS LAYER

#stiffness provided by deck to purlin
#This needs section properties and material properties.
#Calculate for every unique segment along the span.



#Go through the purlin cross-sections to calculate lateral and rotational stiffness.

#Define the spacing between top flange connections that restrain distortional buckling.
