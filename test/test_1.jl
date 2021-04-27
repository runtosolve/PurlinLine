using PurlinLine
using StructuresKit


mutable struct purlin_line_object

    inputs::PurlinLine.Inputs

    cross_section_data::Array{PurlinLine.CrossSectionData}
    free_flange_cross_section_data::Array{PurlinLine.CrossSectionData}

    bracing_data::Array{PurlinLine.BracingData}

    local_buckling_xx_pos::Array{PurlinLine.ElasticBucklingData}
    local_buckling_xx_neg::Array{PurlinLine.ElasticBucklingData}
    local_buckling_yy_pos::Array{PurlinLine.ElasticBucklingData}
    local_buckling_yy_neg::Array{PurlinLine.ElasticBucklingData}
    distortional_buckling_xx_pos::Array{PurlinLine.ElasticBucklingData}
    distortional_buckling_xx_neg::Array{PurlinLine.ElasticBucklingData}

    local_global_flexural_strength_xx::Array{PurlinLine.LocalGlobalFlexuralStrengthData}
    local_global_flexural_strength_yy::Array{PurlinLine.LocalGlobalFlexuralStrengthData}
    local_global_flexural_strength_free_flange_yy::Array{PurlinLine.LocalGlobalFlexuralStrengthData}

    torsion_strength::Array{PurlinLine.TorsionStrengthData}

    shear_strength::Array{PurlinLine.ShearStrengthData}

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

#Capture inputs.
purlin_line.inputs = PurlinLine.Inputs(design_code, loading_direction, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, bridging_locations)

#CALCULATIONS LAYER

#Calculate purlin and purlin free flange section properties.
purlin_line.cross_section_data, purlin_line.free_flange_cross_section_data = PurlinLine.calculate_purlin_section_properties(purlin_line)

#Calculate deck bracing properties. 
purlin_line.bracing_data = PurlinLine.define_deck_bracing_properties(purlin_line)

#Calculate the critical elastic local buckling and distortional buckling properties for each purlin line segment.
purlin_line.local_buckling_xx_pos, purlin_line.local_buckling_xx_neg, purlin_line.local_buckling_yy_pos, purlin_line.local_buckling_yy_neg, purlin_line.distortional_buckling_xx_pos, purlin_line.distortional_buckling_xx_neg  = calculate_elastic_buckling_properties(purlin_line)

#Calculate the local-global flexural strengths for each purlin line segment.   
purlin_line.local_global_flexural_strength_xx, purlin_line.local_global_flexural_strength_yy, purlin_line.local_global_flexural_strength_free_flange_yy = calculate_flexural_strength(purlin_line)

#Need distortional buckling strength!

#Calculate torsion strength for each purlin line segment.
purlin_line.torsion_strength = calculate_torsion_strength(purlin_line)

#Calculate shear strength for each purlin line segment.
purlin_line.shear_strength = calculate_shear_strength(purlin_line)

#Calculate web crippling strength at each support.
purlin_line.web_crippling = calculate_web_crippling_strength(purlin_line)






#ANALYSIS LAYER

#RESULTS LAYER

#stiffness provided by deck to purlin
#This needs section properties and material properties.
#Calculate for every unique segment along the span.



#Go through the purlin cross-sections to calculate lateral and rotational stiffness.

#Define the spacing between top flange connections that restrain distortional buckling.
