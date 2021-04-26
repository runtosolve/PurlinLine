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

    free_flange_cross_section_node_geometry::Vector{Array{Float64}}
    free_flange_cross_section_element_connectivity::Vector{Array{Float64}}

    kp::Vector{Float64}
    kϕ::Vector{Float64}
    kϕ_dist::Vector{Float64}
    Lcrd_AISI::Vector{Float64}
    Lm::Float64
    kx::Vector{Float64}

    Lcrd_pos_CUFSM::Vector{Float64}
    Lcrd_neg_CUFSM::Vector{Float64}
    Mcrd_pos::Vector{Float64}
    Mcrd_neg::Vector{Float64}

    Lcrℓ_xx_pos::Vector{Float64}
    Lcrℓ_xx_neg::Vector{Float64}
    Mcrℓ_xx_pos::Vector{Float64}
    Mcrℓ_xx_neg::Vector{Float64}
    Lcrℓ_yy_pos::Vector{Float64}
    Lcrℓ_yy_neg::Vector{Float64}
    Mcrℓ_yy_pos::Vector{Float64}
    Mcrℓ_yy_neg::Vector{Float64}
    
    CUFSM_local_xx_pos_data::Array{CUFSM.data}
    CUFSM_local_xx_neg_data::Array{CUFSM.data}
    CUFSM_local_yy_pos_data::Array{CUFSM.data}
    CUFSM_local_yy_neg_data::Array{CUFSM.data}
    CUFSM_dist_pos_data::Array{CUFSM.data}
    CUFSM_dist_neg_data::Array{CUFSM.data}

    Sxx_pos::Vector{Float64}
    Sxx_neg::Vector{Float64}
    My_xx_pos::Vector{Float64}
    My_xx_neg::Vector{Float64}
    My_xx::Vector{Float64}
    Mne_xx::Vector{Float64}
    Mnℓ_xx_pos::Vector{Float64}
    Mnℓ_xx_neg::Vector{Float64}
    eMnℓ_xx_pos::Vector{Float64}
    eMnℓ_xx_neg::Vector{Float64}

    Syy_pos::Vector{Float64}
    Syy_neg::Vector{Float64}
    My_yy_pos::Vector{Float64}
    My_yy_neg::Vector{Float64}
    My_yy::Vector{Float64}
    Mne_yy::Vector{Float64}
    Mnℓ_yy_pos::Vector{Float64}
    Mnℓ_yy_neg::Vector{Float64}
    eMnℓ_yy_pos::Vector{Float64}
    eMnℓ_yy_neg::Vector{Float64}

    Syy_pos_free_flange::Vector{Float64}
    Syy_neg_free_flange::Vector{Float64}
    My_yy_pos_free_flange::Vector{Float64}
    My_yy_neg_free_flange::Vector{Float64}
    My_yy_free_flange::Vector{Float64}
    Mne_yy_free_flange::Vector{Float64}
    Mnℓ_yy_pos_free_flange::Vector{Float64}
    Mnℓ_yy_neg_free_flange::Vector{Float64}
    eMnℓ_yy_pos_free_flange::Vector{Float64}
    eMnℓ_yy_neg_free_flange::Vector{Float64}

    Bn::Vector{Float64}
    eBn::Vector{Float64}

    h_flat::Vector{Float64}
    Fcrv::Vector{Float64}
    kv::Vector{Float64}
    Vcr::Vector{Float64}
    Vn::Vector{Float64}
    eVn::Vector{Float64}

    web_crippling::Array{PurlinLine.WebCripplingData}

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

#CALCULATIONS LAYER

##Calculate properties associated with each purlin line segment.

#Calculate purlin and purlin free flange section properties.
purlin_line.section_properties, purlin_line.free_flange_section_properties, purlin_line.cross_section_node_geometry, purlin_line.cross_section_element_connectivity, purlin_line.free_flange_cross_section_node_geometry, purlin_line.free_flange_cross_section_element_connectivity = PurlinLine.calculate_purlin_section_properties(purlin_line.cross_section_dimensions)


#Calculate deck bracing properties. 
purlin_line.kp, purlin_line.kϕ, purlin_line.kϕ_dist, purlin_line.kx, purlin_line.Lm, purlin_line.Lcrd_AISI = PurlinLine.define_deck_bracing_properties(purlin_line)

#Calculate the critical elastic local buckling and distortional buckling properties for each purlin line segment.

purlin_line.CUFSM_local_xx_pos_data, purlin_line.CUFSM_local_xx_neg_data, purlin_line.CUFSM_local_yy_pos_data, purlin_line.CUFSM_local_yy_neg_data, purlin_line.CUFSM_dist_pos_data, purlin_line.CUFSM_dist_neg_data, purlin_line.Mcrℓ_xx_pos, purlin_line.Mcrℓ_xx_neg, purlin_line.Mcrℓ_yy_pos, purlin_line.Mcrℓ_yy_neg, purlin_line.Mcrd_pos, purlin_line.Mcrd_neg, purlin_line.Lcrℓ_xx_pos, purlin_line.Lcrℓ_xx_neg, purlin_line.Lcrℓ_yy_pos, purlin_line.Lcrℓ_yy_neg, purlin_line.Lcrd_pos_CUFSM, purlin_line.Lcrd_neg_CUFSM = calculate_elastic_buckling_properties(purlin_line)


#Calculate the flexural strengths for each purlin line segment.   
purlin_line.Sxx_pos, purlin_line.Sxx_neg, purlin_line.My_xx_pos, purlin_line.My_xx_neg, purlin_line.My_xx, purlin_line.Mne_xx, purlin_line.Mnℓ_xx_pos, purlin_line.Mnℓ_xx_neg, purlin_line.eMnℓ_xx_pos, purlin_line.eMnℓ_xx_neg, purlin_line.Syy_pos, purlin_line.Syy_neg, purlin_line.My_yy_pos, purlin_line.My_yy_neg, purlin_line.My_yy, purlin_line.Mne_yy, purlin_line.Mnℓ_yy_pos, purlin_line.Mnℓ_yy_neg, purlin_line.eMnℓ_yy_pos, purlin_line.eMnℓ_yy_neg, purlin_line.Syy_pos_free_flange, purlin_line.Syy_neg_free_flange, purlin_line.My_yy_pos_free_flange, purlin_line.My_yy_neg_free_flange, purlin_line.My_yy_free_flange, purlin_line.Mne_yy_free_flange, purlin_line.Mnℓ_yy_neg_free_flange, purlin_line.eMnℓ_yy_pos_free_flange, purlin_line.eMnℓ_yy_neg_free_flange = calculate_flexural_strength(purlin_line)

#Calculate torsion strength for each purlin line segment.
purlin_line.Bn, purlin_line.eBn = calculate_torsion_strength(purlin_line)

#Calculate shear strength for each purlin line segment.
purlin_line.h_flat, purlin_line.Fcrv, purlin_line.kv, purlin_line.Vcr, purlin_line.Vn, purlin_line.eVn = calculate_shear_strength(purlin_line)


#Calculate web crippling strength at each support.
purlin_line.web_crippling = calculate_web_crippling_strength(purlin_line)






#ANALYSIS LAYER

#RESULTS LAYER

#stiffness provided by deck to purlin
#This needs section properties and material properties.
#Calculate for every unique segment along the span.



#Go through the purlin cross-sections to calculate lateral and rotational stiffness.

#Define the spacing between top flange connections that restrain distortional buckling.
