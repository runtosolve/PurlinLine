module PurlinLine


using StructuresKit


export calculate_purlin_section_properties, deck_pull_through_fastener_stiffness


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
        xcoords_center = xcoords_center .- xcoords_center[index] .- t/2

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
        xcoords_center = xcoords_center .- xcoords_center[index] .- t/2

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
        purlin_node_geometry, purlin_element_info = define_purlin_cross_section(cross_section_type, t, d1, b1, h, b2, d2, α1, α2, α3, α4, α5, r1, r2, r3, r4, n, n_radius)

        #Define the purlin free flange discretization.
        n = [4, 4, 4]
        n_radius = [4, 4]

        #Define the purlin free flange cross-section nodes and elements.
        purlin_free_flange_node_geometry, purlin_free_flange_element_info = define_purlin_free_flange_cross_section(cross_section_type, t, d1, b1, h, α1, α2, α3, r1, r2, n, n_radius)

        #Calculate the purlin cross-section properties.
        purlin_section_properties[i] = CrossSection.CUFSMsection_properties(purlin_node_geometry, purlin_element_info)

        #Calculate the purlin free flange cross-section properties.
        purlin_free_flange_section_properties[i] = CrossSection.CUFSMsection_properties(purlin_free_flange_node_geometry, purlin_free_flange_element_info)

    end

    return purlin_section_properties, purlin_free_flange_section_properties

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

end # module
