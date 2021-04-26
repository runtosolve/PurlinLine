module PurlinLine


using StructuresKit


export calculate_purlin_section_properties, deck_pull_through_fastener_stiffness, define_deck_bracing_properties, calculate_elastic_buckling_properties, calculate_flexural_strength, calculate_torsion_strength, calculate_shear_strength, calculate_web_crippling_strength, WebCripplingData


struct WebCripplingData

    support_condition::String
    flange_condition::String
    load_case::String
    load_location::String
    C::Float64 
    C_R::Float64
    R::Float64 
    C_N::Float64 
    N::Float64
    C_h::Float64 
    ϕ_w::Float64 
    Ω_w::Float64 
    Pn::Float64
    ePn::Float64

end


function define_purlin_cross_section(cross_section_type, t, d1, b1, h, b2, d2, α1, α2, α3, α4, α5, r1, r2, r3, r4, n, n_radius)

    if cross_section_type == "Z"

        #Define the top Z purlin surface.   For the top flange, this means the out-to-out dimensions.  For the bottom flange, the interior outside dimensions should be used.

        #First calculate the correction on out-to-out length to go from the outside surface to the inside bottom flange surface.
        delta_lip_bottom = t / tan((π - deg2rad(abs(α2 - α1))) / 2)
        delta_web_bottom = t / tan((π - deg2rad(abs(α3 - α2))) / 2)

        #Note here that the bottom flange and lip dimensions are smaller here.
        ΔL = [d1 - delta_lip_bottom, b1 - delta_lip_bottom - delta_web_bottom, h - delta_web_bottom, b2, d2]
        θ = [α1, α2, α3, α4, α5]

        #Note that the outside radius is used at the top flange, and the inside radius is used for the bottom flange.
        radius = [r1 - t, r2 - t, r3, r4]

        closed_or_open = 1

        purlin = CrossSection.Feature(ΔL, θ, n, radius, n_radius, closed_or_open)

        #Calculate the out-to-out purlin surface coordinates.
        xcoords_out, ycoords_out = CrossSection.get_xy_coordinates(purlin)

        #Calculate centerline purlin coordinates.
        unitnormals = CrossSection.surface_normals(xcoords_out, ycoords_out, closed_or_open)
        nodenormals = CrossSection.avg_node_normals(unitnormals, closed_or_open)
        xcoords_center, ycoords_center = CrossSection.xycoords_along_normal(xcoords_out, ycoords_out, nodenormals, -t/2)

        #Shift y coordinates so that the bottom purlin face is at y = 0.
        ycoords_center = ycoords_center .- minimum(ycoords_center) .+ t/2

        #Shift x coordinates so that the purlin web centerline is at x = 0.
        index = floor(Int, length(xcoords_center)/2)
        xcoords_center = xcoords_center .- xcoords_center[index]

    end

    #Add C section here at some point.

    #Package nodal geometry.
    node_geometry = [xcoords_center ycoords_center]

    #Define cross-section element connectivity and thicknesses.
    num_cross_section_nodes = length(xcoords_center)
    element_info = [1:(num_cross_section_nodes - 1) 2:num_cross_section_nodes ones(num_cross_section_nodes - 1) * t]

    return node_geometry, element_info

end



function define_purlin_free_flange_cross_section(cross_section_type, t, d1, b1, h, α1, α2, α3, r1, r2, n, n_radius)

    if cross_section_type == "Z"

        #Define the top Z purlin surface.   For the top flange, this means the out-to-out dimensions.  For the bottom flange, the interior outside dimensions should be used.

        #First calculate the correction on out-to-out length to go from the outside surface to the inside bottom flange surface.
        delta_lip_bottom = t / tan((π - deg2rad(abs(α2 - α1))) / 2)
        delta_web_bottom = t / tan((π - deg2rad(abs(α3 - α2))) / 2)

        #Note here that the bottom flange and lip dimensions are smaller here.
        #Use 1/5 of the web height.
        ΔL = [d1 - delta_lip_bottom, b1 - delta_lip_bottom - delta_web_bottom, h/5 - delta_web_bottom]
        θ = [α1, α2, α3]

        #The inside radius is used for the bottom flange.
        radius = [r1 - t, r2 - t]

        closed_or_open = 1

        purlin = CrossSection.Feature(ΔL, θ, n, radius, n_radius, closed_or_open)

        #Calculate the out-to-out purlin surface coordinates.
        xcoords_out, ycoords_out = CrossSection.get_xy_coordinates(purlin)

        #Calculate centerline purlin coordinates.
        unitnormals = CrossSection.surface_normals(xcoords_out, ycoords_out, closed_or_open)
        nodenormals = CrossSection.avg_node_normals(unitnormals, closed_or_open)
        xcoords_center, ycoords_center = CrossSection.xycoords_along_normal(xcoords_out, ycoords_out, nodenormals, -t/2)

        #Shift y coordinates so that the bottom purlin face is at y = 0.
        ycoords_center = ycoords_center .- minimum(ycoords_center) .+ t/2

        #Shift x coordinates so that the purlin web centerline is at x = 0.
        index = length(xcoords_center)
        xcoords_center = xcoords_center .- xcoords_center[index]

    end

    #Add C section here at some point.

    #Package nodal geometry.
    node_geometry = [xcoords_center ycoords_center]

    #Define cross-section element connectivity and thicknesses.
    num_cross_section_nodes = length(xcoords_center)
    element_info = [1:(num_cross_section_nodes - 1) 2:num_cross_section_nodes ones(num_cross_section_nodes - 1) * t]

    return node_geometry, element_info

end


function calculate_purlin_section_properties(purlin_cross_section_dimensions)

    num_purlin_sections = size(purlin_cross_section_dimensions)[1]

    purlin_section_properties = Array{CrossSection.SectionProperties, 1}(undef, num_purlin_sections)
    purlin_free_flange_section_properties = Array{CrossSection.SectionProperties, 1}(undef, num_purlin_sections)

    purlin_node_geometry = Vector{Array{Float64}}(undef, num_purlin_sections)
    purlin_element_info = Vector{Array{Float64}}(undef, num_purlin_sections)

    purlin_free_flange_node_geometry = Vector{Array{Float64}}(undef, num_purlin_sections)
    purlin_free_flange_element_info = Vector{Array{Float64}}(undef, num_purlin_sections)

    for i=1:num_purlin_sections

        cross_section_type = purlin_cross_section_dimensions[i][1]
        t = purlin_cross_section_dimensions[i][2]
        d1 = purlin_cross_section_dimensions[i][3]
        b1 = purlin_cross_section_dimensions[i][4]
        h = purlin_cross_section_dimensions[i][5]
        b2 = purlin_cross_section_dimensions[i][6]
        d2 = purlin_cross_section_dimensions[i][7]
        α1 = purlin_cross_section_dimensions[i][8]
        α2 = purlin_cross_section_dimensions[i][9]
        α3 = purlin_cross_section_dimensions[i][10]
        α4 = purlin_cross_section_dimensions[i][11]
        α5 = purlin_cross_section_dimensions[i][12]
        r1 = purlin_cross_section_dimensions[i][13]
        r2 = purlin_cross_section_dimensions[i][14]
        r3 = purlin_cross_section_dimensions[i][15]
        r4 = purlin_cross_section_dimensions[i][16]


        #Define the purlin cross-section discretization.
        n = [4, 4, 4, 4, 4]
        n_radius = [4, 4, 4, 4]

        #Define the purlin cross-section nodes and elements.
        purlin_node_geometry[i], purlin_element_info[i] = define_purlin_cross_section(cross_section_type, t, d1, b1, h, b2, d2, α1, α2, α3, α4, α5, r1, r2, r3, r4, n, n_radius)

        #Define the purlin free flange discretization.
        n = [4, 4, 4]
        n_radius = [4, 4]

        #Define the purlin free flange cross-section nodes and elements.
        purlin_free_flange_node_geometry[i], purlin_free_flange_element_info[i] = define_purlin_free_flange_cross_section(cross_section_type, t, d1, b1, h, α1, α2, α3, r1, r2, n, n_radius)

        #Calculate the purlin cross-section properties.
        purlin_section_properties[i] = CrossSection.CUFSMsection_properties(purlin_node_geometry[i], purlin_element_info[i])

        #Calculate the purlin free flange cross-section properties.
        purlin_free_flange_section_properties[i] = CrossSection.CUFSMsection_properties(purlin_free_flange_node_geometry[i], purlin_free_flange_element_info[i])

    end

    return purlin_section_properties, purlin_free_flange_section_properties, purlin_node_geometry, purlin_element_info, purlin_free_flange_node_geometry, purlin_free_flange_element_info

end

function deck_pull_through_fastener_stiffness(deck_material_properties, b_top, t_roof_deck)

    #Define the fastener distance from a major R-panel rib.  Hard coded for now.
    x = 25.4 #mm

    #Need this if statement to make sure units are treated propertly here.
    if deck_material_properties[1] == 203255.0

        #Define the fastener distance from the purlin pivot point.  The purlin pivot point for a Z section is the top flange - web intersection. Assume the fastener location is centered in the purlin top flange.
        c = b_top / 2 #mm

        #Define the roof deck thickness in metric.
        tw = t_roof_deck  #mm

        #Approximate the roof panel base metal thickness.
        kp = Connections.cfs_pull_through_plate_stiffness(x, c, tw)

    elseif deck_material_properties[1] == 29500.0

        #Define the fastener distance from the purlin pivot point.  The purlin pivot point for a Z section is the top flange - web intersection. Assume the fastener location is centered in the purlin top flange.
        c = b_top * 25.4 / 2 #mm

        #Define the roof deck thickness in metric.
        tw = t_roof_deck * 25.4  #mm

        #Approximate the roof panel base metal thickness.
        kp = Connections.cfs_pull_through_plate_stiffness(x, c, tw)

        #Convert kp from N/mm to kips/in.
        kp = kp / 1000 / 4.448 * 25.4

    else
        error("Set the deck elastic modulus to 29500.0 ksi or 203225.0 MPa.")

    end

    return kp

end


function define_deck_bracing_properties(purlin_line)

    num_purlin_segments = size(purlin_line.segments)[1]

    if purlin_line.deck_details[1] == "screw-fastened"

        #Define the deck to purlin screw-fastened connection spacing.
        deck_purlin_fastener_spacing = purlin_line.deck_details[3]

        #Define the deck to purlin screw diameter.
        deck_purlin_fastener_diameter = purlin_line.deck_details[4]

        #Define the nominal shear strength of the typical screw.
        Fss = purlin_line.deck_details[5]

        #Define the roof deck base metal thickness.
        t_roof_deck = purlin_line.deck_details[2]

        #Define roof deck steel elastic modulus.
        E_roof_deck = purlin_line.deck_material_properties[1]

        #Define roof deck steel ultimate yield stress.
        Fu_roof_deck = purlin_line.deck_material_properties[4]

        #Define the distance between fasteners as the distortional discrete bracing length.
        Lm = deck_purlin_fastener_spacing

        #Initialize storage vectors.
        kp = zeros(Float64, num_purlin_segments)
        kϕ = zeros(Float64, num_purlin_segments)
        kϕ_dist = zeros(Float64, num_purlin_segments)
        kx = zeros(Float64, num_purlin_segments)
        Lcrd = zeros(Float64, num_purlin_segments)

        #Loop over all the purlin segments in the line.
        for i = 1:num_purlin_segments

            #Define the section property index associated with purlin segment i.
            section_index = purlin_line.segments[i][2]

            #Define the material property index associated with purlin segment i.
            material_index = purlin_line.segments[i][3]

            #Define purlin steel elastic modulus.
            E_purlin = purlin_line.material_properties[material_index][1]

            #Define purlin steel Poisson's ratio.
            μ_purlin = purlin_line.material_properties[material_index][2]

            #Define purlin steel ultimate stress.
            Fu_purlin = purlin_line.material_properties[material_index][4]

            #Define the purlin top flange width.
            b_top = purlin_line.cross_section_dimensions[section_index][6]

            #Define purlin base metal thickness.
            t_purlin = purlin_line.cross_section_dimensions[section_index][2]

            #Define out-to-out purlin web depth.
            ho = purlin_line.cross_section_dimensions[section_index][5]

            #Define purlin top flange lip length.
            d_top = purlin_line.cross_section_dimensions[section_index][7]

            #Define purlin top flange lip angle from the horizon, in degrees.
            θ_top = purlin_line.cross_section_dimensions[section_index][11] - purlin_line.cross_section_dimensions[section_index][12]

            #Define the location from the purlin top flange pivot point to the fastener.  Assume the fastener is centered in the flange.
            c = b_top/2

            #Define the distance from the purlin web to the screw location in the top flange.  Assume it is centered in the top flange.
            deck_purlin_fastener_top_flange_location = b_top/2

            #Define the deck fastener pull-through plate stiffness.  Assume the fastener is centered between two panel ribs.
            kp[i] = deck_pull_through_fastener_stiffness(purlin_line.deck_material_properties, b_top, t_roof_deck)

            #Apply Cee or Zee binary.
            if purlin_line.cross_section_dimensions[section_index][1] == "Z"
                CorZ = 1
            elseif purlin_line.cross_section_dimensions[section_index][1] == "C"
                CorZ = 0
            end

            #Calculate the rotational stiffness provided to the purlin by the screw-fastened connection between the deck and the purlin.  It is assumed that the deck flexural stiffness is much higher than the connection stiffness.
            kϕ[i] = Connections.cfs_rot_screwfastened_k(b_top, c, deck_purlin_fastener_spacing, t_purlin, kp[i], E_purlin, CorZ)

            #Calculate the purlin distortional buckling half-wavelength.

            #Calculate top flange + lip section properties.
            Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof,  hxf, hyf, yof = AISIS10016.table23131(CorZ, t_purlin, b_top, d_top, θ_top)

            #Calculate the purlin distortional buckling half-wavelength.
            Lcrd[i], L = AISIS10016.app23334(ho, μ_purlin, t_purlin, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

            #If Lcrd is longer than the fastener spacing, then the distortional buckling will be restrained by the deck.
            if Lcrd[i] >= Lm
                kϕ_dist[i] = kϕ[i]
            else
                kϕ_dist[i] = 0.0
            end

            #Approximate the lateral stiffness provided to the top of the purlin by the screw-fastened connection between the deck and the purlin.

            #Calculate the stiffness of a single screw-fastened connection.
            Ka, ψ, α, β, Ke = Connections.cfs_trans_screwfastened_k(t_roof_deck, t_purlin, E_roof_deck, E_purlin, Fss, Fu_roof_deck, Fu_purlin, deck_purlin_fastener_diameter)

            #Convert the discrete stiffness to a distributed stiffness, divide by the fastener spacing.
            kx[i] = Ke / deck_purlin_fastener_spacing

        end

    end

    return kp, kϕ, kϕ_dist, kx, Lm, Lcrd

end

function get_elastic_buckling(prop, node, elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    node_with_stress = CUFSM.stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    curve, shapes = CUFSM.strip(prop, node_with_stress, elem, lengths, springs, constraints, neigs)

    data = CUFSM.data(prop, node_with_stress, elem, lengths, springs, constraints, neigs, curve, shapes)

    half_wavelength = [curve[i,1][1] for i=1:length(lengths)]
    load_factor = [curve[i,1][2] for i=1:length(lengths)]

    Mcr = minimum(load_factor)

    min_index = findfirst(x->x==minimum(load_factor), load_factor)    

    Lcr = half_wavelength[min_index]

    return data, Mcr, Lcr

end

function calculate_elastic_buckling_properties(purlin_line)

    num_purlin_segments = size(purlin_line.segments)[1]

    CUFSM_local_xx_pos_data = Array{CUFSM.data, 1}(undef, num_purlin_segments)
    CUFSM_local_xx_neg_data = Array{CUFSM.data, 1}(undef, num_purlin_segments)
    CUFSM_local_yy_pos_data = Array{CUFSM.data, 1}(undef, num_purlin_segments)
    CUFSM_local_yy_neg_data = Array{CUFSM.data, 1}(undef, num_purlin_segments)
    CUFSM_dist_pos_data = Array{CUFSM.data, 1}(undef, num_purlin_segments)
    CUFSM_dist_neg_data = Array{CUFSM.data, 1}(undef, num_purlin_segments)
    Mcrd_pos = Array{Float64, 1}(undef, num_purlin_segments)
    Mcrd_neg = Array{Float64, 1}(undef, num_purlin_segments)
    Mcrℓ_xx_pos = Array{Float64, 1}(undef, num_purlin_segments)
    Mcrℓ_xx_neg = Array{Float64, 1}(undef, num_purlin_segments)
    Mcrℓ_yy_pos = Array{Float64, 1}(undef, num_purlin_segments)
    Mcrℓ_yy_neg = Array{Float64, 1}(undef, num_purlin_segments)
    Lcrd_pos_CUFSM = Array{Float64, 1}(undef, num_purlin_segments)
    Lcrd_neg_CUFSM = Array{Float64, 1}(undef, num_purlin_segments)
    Lcrℓ_xx_pos = Array{Float64, 1}(undef, num_purlin_segments)
    Lcrℓ_xx_neg = Array{Float64, 1}(undef, num_purlin_segments)
    Lcrℓ_yy_pos = Array{Float64, 1}(undef, num_purlin_segments)
    Lcrℓ_yy_neg = Array{Float64, 1}(undef, num_purlin_segments)

    #Loop over all the purlin segments in the line.
    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = purlin_line.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = purlin_line.segments[i][3]
        
        #Map section properties to CUFSM.
        A = purlin_line.section_properties[section_index].A
        xcg = purlin_line.section_properties[section_index].xc
        zcg = purlin_line.section_properties[section_index].yc
        Ixx = purlin_line.section_properties[section_index].Ixx
        Izz = purlin_line.section_properties[section_index].Iyy
        Ixz = purlin_line.section_properties[section_index].Ixy
        thetap = rad2deg(purlin_line.section_properties[section_index].θ)
        I11 = purlin_line.section_properties[section_index].I1
        I22 = purlin_line.section_properties[section_index].I2
        unsymm = 0  #Sets Ixz=0 if unsymm = 0

        #Define the number of cross-section nodes.
        num_cross_section_nodes = size(purlin_line.cross_section_node_geometry[section_index])[1]

        #Initialize CUFSM node matrix.
        node = zeros(Float64, (num_cross_section_nodes, 8))

        #Add node numbers to node matrix.
        node[:, 1] .= 1:num_cross_section_nodes

        #Add nodal coordinates to node matrix.
        node[:, 2:3] .= purlin_line.cross_section_node_geometry[section_index]

        #Add nodal restraints to node matrix.
        node[:, 4:7] .= ones(num_cross_section_nodes,4)

        #Define number of cross-section elements.
        num_cross_section_elements = size(purlin_line.cross_section_element_connectivity[section_index])[1]

        #Initialize CUFSM elem matrix.
        elem = zeros(Float64, (num_cross_section_elements, 5))

        #Add element numbers to elem matrix.
        elem[:, 1] = 1:num_cross_section_elements

        #Add element connectivity and thickness to elem matrix.
        elem[:, 2:4] .= purlin_line.cross_section_element_connectivity[section_index]

        #Add element material reference to elem matrix.
        elem[:, 5] .= ones(num_cross_section_elements) * 100
                                
                        #lip curve bottom_flange curve web curve top_flange
        center_top_flange_node = 4 + 4 + 4 + 4 + 4 + 4 + 2 + 1
        springs = [1 center_top_flange_node 0 purlin_line.kx[i] 0 0 purlin_line.kϕ_dist[i] 0 0 0]
        constraints = 0

        E = purlin_line.material_properties[material_index][1]
        ν = purlin_line.material_properties[material_index][2]
        G = E / (2 *(1 + ν))
        prop = [100 E E ν ν G]

        neigs = 1  #just need the first mode 

        ###Local buckling - xx axis, positive 

        #Add reference stress to node matrix.

        #Define reference loads.  
        P = 0.0
        Mxx = 1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        h = purlin_line.cross_section_dimensions[section_index][5]  #this is a little dangerous
        length_inc = 5
        lengths = collect(0.25*h:0.75*h/length_inc:1.0*h)   #define to catch the local minimum

        CUFSM_local_xx_pos_data[i], Mcrℓ_xx_pos[i], Lcrℓ_xx_pos[i] = get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)   
        
        #Needed this deepcopy here to make struct work correctly.  Otherwise 'node' just kept changing.


        ###Local buckling - xx axis, negative 

        #Add reference stress to node matrix.

        #Define reference loads.  
        P = 0.0
        Mxx = -1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        h = purlin_line.cross_section_dimensions[section_index][5]  #this is a little dangerous
        length_inc = 5
        lengths = collect(0.25*h:0.75*h/length_inc:1.0*h)   #define to catch the local minimum

        CUFSM_local_xx_neg_data[i], Mcrℓ_xx_neg[i], Lcrℓ_xx_neg[i] = get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)


        ###local buckling - yy axis, positive
        
        #Define reference loads.  
        P = 0.0
        Mxx = 0.0  
        Mzz = 1.0  #assume centroidal moment always for now
        M11 = 0.0
        M22 = 0.0

        #Try Lcrd as a guide for finding the half-wavelength of the flange and lip (unstiffened element).
        length_inc = 5
        lengths = collect(0.25 * purlin_line.Lcrd_AISI[i]:1.0/length_inc:1.25 * purlin_line.Lcrd_AISI[i])   #define to catch the local minimum

        CUFSM_local_yy_pos_data[i], Mcrℓ_yy_pos[i], Lcrℓ_yy_pos[i] = get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

  
        ###local buckling - yy axis, negative
        
        #Define reference loads.  
        P = 0.0
        Mxx = 0.0  
        Mzz = -1.0  #assume centroidal moment always for now
        M11 = 0.0
        M22 = 0.0
    
        length_inc = 5
        #Try Lcrd as a guide for finding the half-wavelength of the flange and lip (unstiffened element).
        lengths = collect(0.25 * purlin_line.Lcrd_AISI[i]:1.0/length_inc:1.25 * purlin_line.Lcrd_AISI[i])   #define to catch the local minimum

        CUFSM_local_yy_neg_data[i], Mcrℓ_yy_neg[i], Lcrℓ_yy_neg[i] = get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

        ###Distortional buckling - xx axis, positive

        #Define reference loads.  
        P = 0.0
        Mxx = 1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        length_inc = 5
        lengths = collect(0.75 * purlin_line.Lcrd_AISI[i]:0.50/length_inc:1.25 * purlin_line.Lcrd_AISI[i])  #define to catch distortional minimum

        CUFSM_dist_pos_data[i], Mcrd_pos[i], Lcrd_pos_CUFSM[i] = get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

         ###Distortional buckling - xx axis, negative

        #Define reference loads.  
        P = 0.0
        Mxx = -1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        length_inc = 5
        lengths = collect(0.75 * purlin_line.Lcrd_AISI[i]:0.50/length_inc:1.25 * purlin_line.Lcrd_AISI[i])  #define to catch distortional minimum

        CUFSM_dist_neg_data[i], Mcrd_neg[i], Lcrd_neg_CUFSM[i] = get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    end

    return CUFSM_local_xx_pos_data, CUFSM_local_xx_neg_data, CUFSM_local_yy_pos_data, CUFSM_local_yy_neg_data, CUFSM_dist_pos_data, CUFSM_dist_neg_data, Mcrℓ_xx_pos, Mcrℓ_xx_neg, Mcrℓ_yy_pos, Mcrℓ_yy_neg, Mcrd_pos, Mcrd_neg, Lcrℓ_xx_pos, Lcrℓ_xx_neg, Lcrℓ_yy_pos, Lcrℓ_yy_neg, Lcrd_pos_CUFSM, Lcrd_neg_CUFSM


end



function calculate_flexural_strength(purlin_line)

    num_purlin_segments = size(purlin_line.segments)[1]

    Sxx_pos = zeros(Float64, num_purlin_segments)
    Sxx_neg = zeros(Float64, num_purlin_segments)
    My_xx_pos = zeros(Float64, num_purlin_segments)
    My_xx_neg = zeros(Float64, num_purlin_segments)
    My_xx = zeros(Float64, num_purlin_segments)
    Mne_xx = zeros(Float64, num_purlin_segments)
    Mnℓ_xx_pos = zeros(Float64, num_purlin_segments)
    Mnℓ_xx_neg = zeros(Float64, num_purlin_segments)
    eMnℓ_xx_pos = zeros(Float64, num_purlin_segments)
    eMnℓ_xx_neg = zeros(Float64, num_purlin_segments)

    Syy_pos = zeros(Float64, num_purlin_segments)
    Syy_neg= zeros(Float64, num_purlin_segments)
    My_yy_pos = zeros(Float64, num_purlin_segments)
    My_yy_neg = zeros(Float64, num_purlin_segments)
    My_yy = zeros(Float64, num_purlin_segments)
    Mne_yy = zeros(Float64, num_purlin_segments)
    Mnℓ_yy_pos = zeros(Float64, num_purlin_segments)
    Mnℓ_yy_neg = zeros(Float64, num_purlin_segments)
    eMnℓ_yy_pos = zeros(Float64, num_purlin_segments)
    eMnℓ_yy_neg = zeros(Float64, num_purlin_segments)

    Syy_pos_free_flange = zeros(Float64, num_purlin_segments)
    Syy_neg_free_flange = zeros(Float64, num_purlin_segments)
    My_yy_pos_free_flange = zeros(Float64, num_purlin_segments)
    My_yy_neg_free_flange = zeros(Float64, num_purlin_segments)
    My_yy_free_flange = zeros(Float64, num_purlin_segments)
    Mne_yy_free_flange = zeros(Float64, num_purlin_segments)
    Mnℓ_yy_pos_free_flange = zeros(Float64, num_purlin_segments)
    Mnℓ_yy_neg_free_flange = zeros(Float64, num_purlin_segments)
    eMnℓ_yy_pos_free_flange = zeros(Float64, num_purlin_segments)
    eMnℓ_yy_neg_free_flange = zeros(Float64, num_purlin_segments)

    if purlin_line.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif purlin_line.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    end

    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = purlin_line.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = purlin_line.segments[i][3]

        ###strong axis flexure, local-global interaction
        Fy = purlin_line.material_properties[material_index][3]
        Ixx = purlin_line.section_properties[section_index].Ixx
        ho = purlin_line.cross_section_dimensions[section_index][5]
        cy_bottom = purlin_line.section_properties[section_index].yc  #distance from neutral axis to bottom outer fiber
        cy_top = ho - cy_bottom #distance from neutral axis to top outer fiber
        Sxx_pos[i] = Ixx/cy_top
        Sxx_neg[i] = Ixx/cy_bottom
        My_xx_pos[i] = Fy*Sxx_pos[i]
        My_xx_neg[i] = Fy*Sxx_neg[i]
        My_xx[i] = minimum([My_xx_pos[i] My_xx_neg[i]])  #first yield criterion for AISI 
        Mne_xx[i] = My_xx[i]  #handle global buckling in the ThinWalledBeam second order analysis

        Mnℓ_xx_pos[i], eMnℓ_xx_pos[i] =  AISIS10016.f321(Mne_xx[i], purlin_line.Mcrℓ_xx_pos[i], ASDorLRFD)

        Mnℓ_xx_neg[i], eMnℓ_xx_neg[i] =  AISIS10016.f321(Mne_xx[i], purlin_line.Mcrℓ_xx_neg[i], ASDorLRFD)

        ###weak axis flexure, local-global interaction
        Iyy = purlin_line.section_properties[section_index].Ixx

        #distance from neutral axis to (-x or left) outer fiber
        #Positive moment is applied when this outer fiber is compressed.
        cx_minusx = purlin_line.section_properties[section_index].xc - minimum(purlin_line.cross_section_node_geometry[section_index][:,1])
        #distance from neutral axis to (+x or right) outer fiber
        #Negative moment is applied when this outer fiber is compressed.
        cx_plusx = maximum(purlin_line.cross_section_node_geometry[section_index][:,1]) - purlin_line.section_properties[section_index].xc 

        Syy_pos[i] = Iyy / cx_minusx
        Syy_neg[i] = Iyy / cx_plusx

        My_yy_pos[i] = Fy*Syy_pos[i]
        My_yy_neg[i] = Fy*Syy_neg[i]
        My_yy[i] = minimum([My_yy_pos[i] My_yy_neg[i]])  #first yield criterion for AISI 
        Mne_yy[i] = My_yy[i]

        Mnℓ_yy_pos[i], eMnℓ_yy_pos[i] = AISIS10016.f321(Mne_yy[i], purlin_line.Mcrℓ_yy_pos[i], ASDorLRFD)

        Mnℓ_yy_neg[i], eMnℓ_yy_neg[i] = AISIS10016.f321(Mne_yy[i], purlin_line.Mcrℓ_yy_neg[i], ASDorLRFD)


        ###free flange yy-axis, local-global interaction

        #define free flange properties
        Iyyf = purlin_line.free_flange_section_properties[section_index].Iyy

        #distance from neutral axis to (-x or left) outer fiber
        #Positive moment is applied when this outer fiber is compressed.
        cxf_minusx = purlin_line.free_flange_section_properties[section_index].xc - minimum(purlin_line.free_flange_cross_section_node_geometry[section_index][:,1])
        #distance from neutral axis to (+x or right) outer fiber
        #Negative moment is applied when this outer fiber is compressed.
        cxf_plusx = maximum(purlin_line.free_flange_cross_section_node_geometry[section_index][:,1]) - purlin_line.free_flange_section_properties[section_index].xc 

        Syy_pos_free_flange[i] = Iyyf / cxf_minusx
        Syy_neg_free_flange[i] = Iyyf / cxf_plusx

        My_yy_pos_free_flange[i] = Fy*Syy_pos_free_flange[i]
        My_yy_neg_free_flange[i] = Fy*Syy_neg_free_flange[i]
        My_yy_free_flange[i] = minimum([My_yy_pos_free_flange[i] My_yy_neg_free_flange[i]])  #first yield criterion for AISI 
        Mne_yy_free_flange[i] = My_yy_free_flange[i]

        #Assume no local buckling for now in the free flange strength calculation.  Set Mcrℓ to Mne times a big number. 

        Mnℓ_yy_pos_free_flange[i], eMnℓ_yy_pos_free_flange[i] = AISIS10016.f321(Mne_yy_free_flange[i], Mne_yy_free_flange[i] * 1000, ASDorLRFD)

        Mnℓ_yy_neg_free_flange[i], eMnℓ_yy_neg_free_flange[i] = AISIS10016.f321(Mne_yy_free_flange[i], Mne_yy_free_flange[i] * 1000, ASDorLRFD)


    end

    return Sxx_pos, Sxx_neg, My_xx_pos, My_xx_neg, My_xx, Mne_xx, Mnℓ_xx_pos, Mnℓ_xx_neg, eMnℓ_xx_pos, eMnℓ_xx_neg, Syy_pos, Syy_neg, My_yy_pos, My_yy_neg, My_yy, Mne_yy, Mnℓ_yy_pos, Mnℓ_yy_neg, eMnℓ_yy_pos, eMnℓ_yy_neg, Syy_pos_free_flange, Syy_neg_free_flange, My_yy_pos_free_flange, My_yy_neg_free_flange, My_yy_free_flange, Mne_yy_free_flange, Mnℓ_yy_neg_free_flange, eMnℓ_yy_pos_free_flange, eMnℓ_yy_neg_free_flange

end

function calculate_torsion_strength(purlin_line)

    num_purlin_segments = size(purlin_line.segments)[1]

    Bn = zeros(Float64, num_purlin_segments)
    eBn = zeros(Float64, num_purlin_segments)
    
    if purlin_line.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif purlin_line.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    end

    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = purlin_line.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = purlin_line.segments[i][3]
        
        Cw = purlin_line.section_properties[section_index].Cw
        Fy = purlin_line.material_properties[material_index][3]

        #This is the maximum magnitude of the warping stress function.  
        Wn = maximum(abs.(purlin_line.section_properties[section_index].wn))

        Bn[i], eBn[i] = AISIS10024.h411(Cw, Fy, Wn, ASDorLRFD)

    end

    return Bn, eBn

end


function calculate_shear_strength(purlin_line)

    num_purlin_segments = size(purlin_line.segments)[1]

    h_flat = zeros(Float64, num_purlin_segments)
    Fcrv = zeros(Float64, num_purlin_segments)
    kv = zeros(Float64, num_purlin_segments)
    Vcr = zeros(Float64, num_purlin_segments)
    Vn = zeros(Float64, num_purlin_segments)
    eVn = zeros(Float64, num_purlin_segments)

    if purlin_line.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif purlin_line.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    end

    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = purlin_line.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = purlin_line.segments[i][3]

        #Set a, the shear stiffener spacing, to the sum of the purlin segment lengths.  This assumes that shear stiffeners are not provided.
        sum_purlin_segments = sum([purlin_line.segments[i][1] for i=1:size(purlin_line.segments)[1]])
        a = sum_purlin_segments

        #Define base metal thickness.
        t = purlin_line.cross_section_dimensions[section_index][2]

        #Define material properties.
        E = purlin_line.material_properties[material_index][1]
        μ = purlin_line.material_properties[material_index][2]
        Fy = purlin_line.material_properties[material_index][3]

        #Depth of flat portion of web.
        full_web_depth = purlin_line.cross_section_dimensions[section_index][5]
        bottom_flange_web_outside_radius = purlin_line.cross_section_dimensions[section_index][14]
        top_flange_web_outside_radius = purlin_line.cross_section_dimensions[section_index][15]
        h_flat[i] = full_web_depth - bottom_flange_web_outside_radius - top_flange_web_outside_radius

        #Calculate plate buckling coefficient.
        kv[i]  = AISIS10016.g233(a, h_flat[i])

        #Calculate shear buckling stress.
        Fcrv[i] = AISIS10016.g232(E, μ, kv[i], h_flat[i], t)
        Vcr[i] = AISIS10016.g231(h_flat[i], t, Fcrv[i])

        #Calculate shear buckling strength.
        Vn[i], eVn[i] = AISIS10016.g21(E, h_flat[i], t, Fy, Vcr[i], ASDorLRFD)

    end

    return h_flat, Fcrv, kv, Vcr, Vn, eVn 

end 


#Calculate the web crippling strength at each support location.
function calculate_web_crippling_strength(purlin_line)

    ###Assumptions...
    #Purlin is always fastened to a support.
    #Purlin always has stiffened or partially stiffened flanges.
    #The loading is always a one-flange loading.

    if purlin_line.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif purlin_line.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    end

    #Define the number of supports along the purlin line.
    num_supports = length(purlin_line.support_locations)

    #Initialize a vector that will hold all the web crippling inputs and outputs.
    web_crippling = Array{WebCripplingData, 1}(undef, num_supports)

    #Define coordinates along purlin line where segment properties change.
    purlin_range = [0; cumsum([purlin_line.segments[i][1] for i=1:size(purlin_line.segments)[1]])]
          
    for i = 1:num_supports

        #Find purlin segment that coincides with a support.
        purlin_range_indices = findall(x->x <= purlin_line.support_locations[i], purlin_range)
        if purlin_range_indices == [1]
            segment_index = 1
        else
            segment_index = maximum(purlin_range_indices) - 1
        end 

        #Define if support is at the end or in the interior of a purlin line.
        if (purlin_line.support_locations[i] == purlin_range[1]) | (purlin_line.support_locations[i] == purlin_range[end])
            load_location = "End"
        else
            load_location = "Interior"
        end

        #Define section and material indices to use for web crippling definitions.
        section_index = purlin_line.segments[segment_index][2]
        material_index = purlin_line.segments[segment_index][3]

        t = purlin_line.cross_section_dimensions[section_index][2]
        Fy = purlin_line.material_properties[material_index][3]
       
        full_web_depth = purlin_line.cross_section_dimensions[section_index][5]
        bottom_flange_web_outside_radius = purlin_line.cross_section_dimensions[section_index][14]
        top_flange_web_outside_radius = purlin_line.cross_section_dimensions[section_index][15]
        h_flat = full_web_depth - bottom_flange_web_outside_radius - top_flange_web_outside_radius

        θ = purlin_line.cross_section_dimensions[section_index][10]  #angle between web plane and surface plane 

        if purlin_line.cross_section_dimensions[section_index][1] == "Z"

            #Use AISI S100-16 Table G5-3 for Z-sections.
            table_g53 = AISIS10016.table_g53()  
           

            web_crippling_coeff = filter(row -> row.support_condition == "Fastened to Support", table_g53)
            web_crippling_coeff = filter(row -> row.flange_condition == "Stiffened or Partially Stiffened Flanges", web_crippling_coeff)
            web_crippling_coeff = filter(row -> row.load_case == "One-Flange Loading or Reaction", web_crippling_coeff)
            web_crippling_coeff = filter(row -> row.load_location== load_location, web_crippling_coeff)

            C = web_crippling_coeff.C[1]
            C_R = web_crippling_coeff.C_R[1]
            R = purlin_line.cross_section_dimensions[section_index][14] - t  #inside radius
            C_N = web_crippling_coeff.C_N[1]
            N = purlin_line.frame_flange_width
            C_h = web_crippling_coeff.C_h[1]
            ϕ_w = web_crippling_coeff.LRFD[1]
            Ω_w = web_crippling_coeff.ASD[1]

            Pn, ePn = AISIS10016.g51(t, h_flat, Fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ASDorLRFD)

            web_crippling[i] = WebCripplingData(web_crippling_coeff.support_condition[1], web_crippling_coeff.flange_condition[1], web_crippling_coeff.load_case[1], web_crippling_coeff.load_location[1], C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, Pn, ePn)

        end

    end

    return web_crippling

end

end # module
