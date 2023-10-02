# using Pkg

# Pkg.activate(".")

using PurlinLine

design_code = "AISI S100-16 ASD"

# Define the properties of each purlin segment along the line.
                #length, section_properties, material_properties
segments = [(25.0*12, 25.0, 1, 1)]

# Define the purlin spacing.
spacing = 60;  #in.

# Define the roof slope.
roof_slope = 0.0;   #degrees

#Define the purlin cross-section type.
# cross_section_dimensions = [("C", 0.059, 0.91, 2.5, 8.0, 2.5, 0.91, -90.0, 180.0, 90.0, 0.0, -90.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059)]

cross_section_dimensions = [("Z", 0.059, 0.91, 2.5, 8.0, 2.5, 0.91, -55.0, 0.0, 90.0, 0.0, -55.0, 3*0.059, 3*0.059, 3*0.059, 3*0.059)]

material_properties = [(29500.0, 0.30, 55.0, 70.0)]

#type="vertical leg standing seam", clip spacing"
deck_details = ("vertical leg standing seam", 18.0)

deck_material_properties = (29500.0, 0.30, 55.0, 70.0)

frame_flange_width = 10.0

support_locations = [0.0, 25.0*12]

purlin_frame_connections = "bottom flange connection"
# purlin_frame_connections = "anti-roll clip"

bridging_locations = Float64[ ]

#Calculate purlin line design variables from user inputs and store them in the data structure.

inputs = PurlinLine.Inputs(design_code, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, purlin_frame_connections, bridging_locations)

purlin_line = PurlinLine.build(inputs)

purlin_line.loading_direction = "gravity"

#Perform a test to collapse.
purlin_line = PurlinLine.test(purlin_line)

purlin_line.applied_pressure * 1000 * 144


using Plots 
plot(purlin_line.model.inputs.z, purlin_line.expected_strengths.eMnd_xx)

plot(purlin_line.model.inputs.z, purlin_line.Β_distortional_gradient_factor)



v = [1, 2, 3, 5] #example list
x = 5 #value to insert
index = searchsortedfirst(v, x)



using Dierckx, S100AISI


function calculate_distortional_buckling_gradient_factor(Mxx, z, Lcrd)

    
    spl = Spline1D(z, Mxx)

    num_halfwavelengths = floor(maximum(z)/Lcrd)

    z_Lcrd = 0.0:maximum(z)/num_halfwavelengths:maximum(z)

    z_Lcrd_start = z_Lcrd[1:end-1]

    z_Lcrd_end = z_Lcrd[2:end]

    M_start = [spl(z_Lcrd_start[i]) for i in eachindex(z_Lcrd_start)]

    M_end = [spl(z_Lcrd_end[i]) for i in eachindex(z_Lcrd_end)]

    M1 = [minimum([abs(M_start[i]), abs(M_end[i])]) for i in eachindex(M_start)]
    M2 = -[maximum([abs(M_start[i]), abs(M_end[i])]) for i in eachindex(M_start)]

    Lm = Lcrd  .* ones(Float64, length(M1))
    L = Lcrd  .* ones(Float64, length(M1))

    Β_range = S100AISI.v16.app23333.(L, Lm, M1, M2)

    Β = Array{Float64}(undef, length(z))

    for i in eachindex(z)

        index = searchsortedfirst(z_Lcrd_start, z[i])

        if index > length(z)
            index = index - 1
        end

        Β[i] = Β_range[index]

    end

    index = findall(x->isnan(x), Β)
    Β[index] .= 1.3  # fix NaN problem when M2 is zero

    return Β

end



β = calculate_distortional_buckling_gradient_factor(purlin_line.internal_forces.Mxx, purlin_line.model.inputs.z, purlin_line.distortional_buckling_xx_pos[1].Lcr)


plot(purlin_line.model.inputs.z, β)
