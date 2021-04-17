using PurlinLineDesign
using StructuresKit

#INPUT LAYER


mutable struct purlin_line_object

    design_code::String
    loading_direction::String
    span_geometry::Vector{Tuple{Float64, Int64, Int64}}
    spacing::Float64
    roof_slope::Float64
    cross_section_dimensions::Vector{Tuple{String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}
    material_properties::Vector{NTuple{4, Float64}}
    deck_details::Tuple{String, Float64, Float64, Float64, Float64, Float64}
    deck_material_properties::NTuple{4, Float64}
    frame_flange_width::Float64
    support_locations::Vector{Float64}
    bridging_locations::Vector{Float64}

    purlin_line_object() = new()

end


design_code = "AISI S100-16 ASD"

loading_direction = "gravity"   #or "uplift"

# Define the properties of each purlin segment along the line.
                #length, section_properties, material_properties
span_geometry = [(23.0*12, 1, 1),
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
purlin_line.span_geometry = span_geometry
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
purlin_section_properties, purlin_free_flange_section_properties = PurlinLineDesign.calculate_purlin_section_properties(purlin_cross_section_dimensions)

plot(purlin_free_flange_section_properties[1].node_geometry[:, 1],purlin_free_flange_section_properties[1].node_geometry[:,2], seriestype = scatter)

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

num_purlin_sections = size(purlin_cross_section_dimensions)[1]


#Go through the purlin cross-sections to calculate lateral and rotational stiffness.

#Define the spacing between top flange connections that restrain distortional buckling.
Lm = deck_purlin_fastener_spacing

function define_purlin_line_stiffness(purlin_line)

if deck_details[1] == "screw-fastened"

    #Define the deck to purlin screw-fastened connection spacing.
    deck_purlin_fastener_spacing = deck_details[3]

    #Define the deck to purlin screw diameter.
    deck_purlin_fastener_diameter = deck_details[4]

    #Define the nominal shear strength of the typical screw.
    Fss = deck_details[5]

    #Define the roof deck base metal thickness.
    t_roof_deck = deck_details[2]

    #Define roof deck steel elastic modulus.
    E_roof_deck = deck_material_properties[1]

    #Define roof deck steel ultimate yield stress.
    Fu_roof_deck = deck_material_properties[4]

    #Define purlin steel elastic modulus.
    E_purlin = purlin_material_properties[1][1]

    #Define purlin steel Poisson's ratio.
    μ_purlin = purlin_material_properties[1][2]

    #Define purlin steel ultimate yield stress.
    Fu_purlin =



    for i = 1:num_purlin_sections

        #Define the purlin top flange width.
        b_top = purlin_cross_section_dimensions[i][6]

        #Define purlin base metal thickness.
        t_purlin = purlin_cross_section_dimensions[i][2]

        #Define out-to-out purlin web depth.
        ho = purlin_cross_section_dimensions[i][5]

        #Define purlin top flange lip length.
        d_top = purlin_cross_section_dimensions[i][7]

        #Define purlin top flange lip angle from the horizon, in degrees.
        θ_top = purlin_cross_section_dimensions[i][11] - purlin_cross_section_dimensions[i][12]

        #Define the location from the purlin top flange pivot point to the fastener.  Assume the fastener is centered in the flange.
        c = b_top/2

        #Define the distance from the purlin web to the screw location in the top flange.  Assume it is centered in the top flange.
        deck_purlin_fastener_top_flange_location = b_top/2

        #Define the deck fastener pull-through plate stiffness.  Assume the fastener is centered between two panel ribs.
        kp = deck_pull_through_fastener_stiffness(deck_material_properties, b_top, t_roof_deck)

        #Calculate the rotational stiffness provided to the purlin by the screw-fastened connection between the deck and the purlin.  It is assumed that the deck flexural stiffness is much higher than the connection stiffness.
        kϕ = Connections.cfs_rot_screwfastened_k(b_top, c, deck_purlin_fastener_spacing, t_purlin, kp, E_purlin, CorZ)

        #Calculate the purlin distortional buckling half-wavelength.

        if deck_details[1] == "Z"
            CorZ = 1
        elseif deck_details[1] == "C"
            CorZ = 0
        end

        #Calculate top flange + lip section properties.
        Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof,  hxf, hyf, yof = AISIS10016.table23131(CorZ, t_purlin, b_top, d_top, θ_top)

        #Calculate the purlin distortional buckling half-wavelength.
        Lcrd, L = AISIS10016.app23334(ho, μ_purlin, t_purlin, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

        #If Lcrd is longer than the fastener spacing, then the distortional buckling will be restrained by the deck.
        if Lcrd >= Lm
            kϕ_dist = kϕ
        else
            kϕ_dist = 0.0
        end

        #Approximate the lateral stiffness provided to the top of the purlin by the screw-fastened connection between the deck and the purlin.

        #Calculate the stiffness of a single screw-fastened connection.
        Ka, ψ, α, β, Ke = Connections.cfs_trans_screwfastened_k(t_roof_deck, t_purlin, E_roof_deck, E_purlin, Fss, Fu_roof_deck, Fu_purlin, deck_purlin_fastener_diameter)

        #Convert the discrete stiffness to a distributed stiffness, divide by the fastener spacing.
        kx = Ke / deck_purlin_fastener_spacing
