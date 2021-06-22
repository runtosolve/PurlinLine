using PurlinLine
using StructuresKit

design_code = "AISI S100-16 ASD"

# loading_direction = "gravity"   #or "uplift"

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

pressure = -10 / 1000 / 144 #psf to kips/in^2

#Calculate purlin line design variables from user inputs and store them in the data structure.
purlin_line = PurlinLine.define(design_code, pressure, segments, spacing, roof_slope, cross_section_dimensions, material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, bridging_locations)

#Translate purlin_line design variables to ThinWalledBeam design variables.
z, m, member_definitions, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports = PurlinLine.thin_walled_beam_interface(purlin_line)

#Set up ThinWalledBeam model.
model = ThinWalledBeam.define(z, m, member_definitions, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)

#Solve ThinWalledBeam model.
model = ThinWalledBeam.solve(model)

#Add second order solution to PurlinLine data structure.
purlin_line.model = model


Mxx = InternalForces.moment(purlin_line.model.z, purlin_line.model.m, -purlin_line.model.v, purlin_line.model.E, purlin_line.model.Ix)
Myy = InternalForces.moment(purlin_line.model.z, purlin_line.model.m, -purlin_line.model.u, purlin_line.model.E, purlin_line.model.Iy)
Vyy = InternalForces.shear(purlin_line.model.z, purlin_line.model.m, -purlin_line.model.v, purlin_line.model.E, purlin_line.model.Ix)
B = InternalForces.bimoment(purlin_line.model.z, purlin_line.model.m, purlin_line.model.ϕ, purlin_line.model.E, purlin_line.model.Cw)


z, m, member_definitions, section_properties, material_properties, kxf, kyf, kϕf, kH, hxf, hyf, axf, ayf, end_boundary_conditions, supports = beam_column_interface(purlin_line)

#Calculate axial force in free flange
Pf = calculate_free_flange_axial_force(Mxx, member_definitions, purlin_line)

#Apply the shear flow based on the y-direction load along the purlin line.
qxf = qy .* kH

#The y-direction load is assumed to be zero.
num_nodes = length(z)
qyf = zeros(Float64, num_nodes)

free_flange_model = BeamColumn.define(z, m, member_definitions, section_properties, material_properties, kxf, kyf, kϕf, hxf, hyf, qxf, qyf, Pf, axf, ayf, end_boundary_conditions, supports)

free_flange_model = BeamColumn.solve(free_flange_model)

purlin_line.free_flange_model = free_flange_model

 #distance from neutral axis to (-x or left) outer fiber
#Positive moment is applied when this outer fiber is compressed.
Myy_free_flange = InternalForces.moment(purlin_line.model.z, purlin_line.model.m, -purlin_line.free_flange_model.u, purlin_line.model.E, purlin_line.free_flange_model.Iy)

#find flexure+torsion D/C



function bending_torsion_DC(Mxx, Myy, B, Myy_free_flange, eMnℓ_xx, eMnℓ_yy, eBn, eMnℓ_yy_free_flange)

    #check bending + torsion interaction
    #ActionM1, ActionM2, ActionB, Interaction
    action_Mxx, action_Myy, action_B, action_Myy_free_flange, interaction = AISIS10024.h42(Mxx, Myy, B, Myy_free_flange, eMnℓ_xx, eMnℓ_yy, eBn, eMnℓ_yy_free_flange)

    DC = interaction / 1.15   #Consider updating this 1.15 in the future.

    return action_Mxx, action_Myy, action_B, action_Myy_free_flange, interaction, DC

end


function calculate_flexural_capacity_envelope(m, eMn_pos, eMn_neg, M)

    num_nodes = length(m)

    eMn_pos_all = [eMn_pos[m[i]] for i=1:num_nodes]
    eMn_neg_all = [eMn_neg[m[i]] for i=1:num_nodes]

    eMn_all = zeros(Float64, num_nodes)

    for i in eachindex(M)

        if M[i] >= 0.0

            eMn_all[i] = eMn_pos_all[i]
            
        elseif M[i] < 0.0

            eMn_all[i] = eMn_neg_all[i]

        end

    end

    return eMn_all

end

function calculate_bending_torsion_DC(purlin_line, Mxx, Myy, B, Myy_free_flange)

    num_purlin_segments = size(purlin_line.inputs.segments)[1]
    eMnℓ_xx_pos_range = [purlin_line.local_global_flexural_strength_xx[i].eMnℓ_pos for i=1:num_purlin_segments]
    eMnℓ_xx_neg_range = [purlin_line.local_global_flexural_strength_xx[i].eMnℓ_neg for i=1:num_purlin_segments]
    eMnℓ_xx_all = calculate_flexural_capacity_envelope(m, eMnℓ_xx_pos_range, eMnℓ_xx_neg_range, Mxx)

    eMnℓ_yy_pos_range = [purlin_line.local_global_flexural_strength_yy[i].eMnℓ_pos for i=1:num_purlin_segments]
    eMnℓ_yy_neg_range = [purlin_line.local_global_flexural_strength_yy[i].eMnℓ_neg for i=1:num_purlin_segments]
    eMnℓ_yy_all = calculate_flexural_capacity_envelope(m, eMnℓ_yy_pos_range, eMnℓ_yy_neg_range, Myy)

    eBn_range = [purlin_line.torsion_strength[i].eBn for i=1:num_purlin_segments]
    eBn_all = zeros(Float64, num_nodes)
    eBn_all .= eBn_range[m]

    #There is no positive or negative capacity here because a first yield criteria is used to determine strength.  Local buckling is not considered.
    eMnℓ_yy_free_flange_range = [purlin_line.yielding_flexural_strength_free_flange_yy[i].eMy for i=1:num_purlin_segments]
    eMnℓ_yy_free_flange_all = zeros(Float64, num_nodes)
    eMnℓ_yy_free_flange_all .= eMnℓ_yy_free_flange_range[m]

    results = AISIS10024.h42.(Mxx, Myy, B, Myy_free_flange, eMnℓ_xx_all, eMnℓ_yy_all, eBn_all, eMnℓ_yy_free_flange_all)

    interaction = [results[i][5] for i=1:num_nodes]

    DC = interaction ./ 1.15   #Consider updating this 1.15 in the future based on AISI COS discussions.

    return DC

end


DC_flexure_torsion = calculate_bending_torsion_DC(purlin_line, Mxx, Myy, B, Myy_free_flange)
DC_distortional = calculate_distortional_buckling_DC(purlin_line, Mxx)
DC_flexure_shear = calculate_flexure_shear_DC(purlin_line, Vyy, Mxx)
DC_biaxial_bending = calculate_biaxial_bending_DC(purlin_line, Mxx, Myy)

#check web crippling
#check web crippling + bending (need to differentiate between two Zs and one web)



#find distortional buckling D/C
function calculate_distortional_buckling_DC(purlin_line, Mxx)

    num_purlin_segments = size(purlin_line.inputs.segments)[1]
    eMnd_xx_pos_range = [purlin_line.distortional_flexural_strength_xx[i].eMnd_pos for i=1:num_purlin_segments]
    eMnd_xx_neg_range = [purlin_line.distortional_flexural_strength_xx[i].eMnd_neg for i=1:num_purlin_segments]
    eMnd_xx_all = calculate_flexural_capacity_envelope(m, eMnd_xx_pos_range, eMnd_xx_neg_range, Mxx)

    #check distortional buckling
    DC = abs.(Mxx./eMnd_xx_all)

    return DC

end

#find flexure+shear D/C
function calculate_flexure_shear_DC(purlin_line, Vyy, Mxx)

    num_purlin_segments = size(purlin_line.inputs.segments)[1]
    eMnℓ_xx_pos_range = [purlin_line.local_global_flexural_strength_xx[i].eMnℓ_pos for i=1:num_purlin_segments]
    eMnℓ_xx_neg_range = [purlin_line.local_global_flexural_strength_xx[i].eMnℓ_neg for i=1:num_purlin_segments]
    eMnℓ_xx_all = calculate_flexural_capacity_envelope(m, eMnℓ_xx_pos_range, eMnℓ_xx_neg_range, Mxx)

    eVn_range = [purlin_line.shear_strength[i].eVn for i=1:num_purlin_segments]
    eVn_all = zeros(Float64, num_nodes)
    eVn_all .= eVn_range[m]

    interaction = AISIS10016.h21.(Mxx, Vyy, eMnℓ_xx_all, eVn_all)

    DC = interaction

    return DC

end


#find biaxial bending D/C
function calculate_biaxial_bending_DC(purlin_line, Mxx, Myy)

    #no axial force for now
    Pbar=zeros(Float64, length(Mxx))
    Pa=ones(Float64, length(Mxx))

    num_purlin_segments = size(purlin_line.inputs.segments)[1]
    eMnℓ_xx_pos_range = [purlin_line.local_global_flexural_strength_xx[i].eMnℓ_pos for i=1:num_purlin_segments]
    eMnℓ_xx_neg_range = [purlin_line.local_global_flexural_strength_xx[i].eMnℓ_neg for i=1:num_purlin_segments]
    eMnℓ_xx_all = calculate_flexural_capacity_envelope(m, eMnℓ_xx_pos_range, eMnℓ_xx_neg_range, Mxx)

    eMnℓ_yy_pos_range = [purlin_line.local_global_flexural_strength_yy[i].eMnℓ_pos for i=1:num_purlin_segments]
    eMnℓ_yy_neg_range = [purlin_line.local_global_flexural_strength_yy[i].eMnℓ_neg for i=1:num_purlin_segments]
    eMnℓ_yy_all = calculate_flexural_capacity_envelope(m, eMnℓ_yy_pos_range, eMnℓ_yy_neg_range, Myy)

    results = AISIS10016.h121.(Pbar, Mxx, Myy, Pa, eMnℓ_xx_all, eMnℓ_yy_all)

    actionP = [x[1] for x in results]
    actionMxx = [x[2] for x in results]
    actionMyy = [x[3] for x in results]
    interaction = [x[4] for x in results]

    DC = interaction

    return DC

end