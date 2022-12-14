module UI

using PurlinLine, Plots, LinesCurvesNodes, CrossSection


function define_lap_section_types(purlin_size_span_assignment)

	#There are 3 possible combinations of purlin types at a lap: 1-1, 2-2, or 1-2.

	num_laps = length(purlin_size_span_assignment) - 1

	lap_section_types = Array{String}(undef, num_laps)

	for i=1:num_laps

		purlin_span_1 = purlin_size_span_assignment[i]
		purlin_span_2 = purlin_size_span_assignment[i+1]

		if (purlin_span_1 == purlin_span_2) & (purlin_span_1 == 1)

			lap_section_types[i] = "1-1"

		elseif (purlin_span_1 == purlin_span_2) & (purlin_span_1 == 2)

			lap_section_types[i] = "2-2"

		elseif purlin_span_1 != purlin_span_2

			lap_section_types[i] = "1-2"

		end

	end

	return lap_section_types

end

function define_lap_section_index(lap_section_types)

	num_lap_sections = length(lap_section_types)
	
	num_unique_lap_sections = length(unique(lap_section_types))

	unique_lap_sections = unique(lap_section_types)

	lap_section_index = Array{Int64}(undef, num_lap_sections)

	for i = 1:num_unique_lap_sections

		index = findall(x->x==unique_lap_sections[i], lap_section_types)

		lap_section_index[index] .= 2 + i

	end

	return lap_section_index

end

function define_lap_segments(purlin_laps, purlin_size_span_assignment)

	num_interior_supports = trunc(Int, length(purlin_laps)/2)

	lap_segments = Array{Tuple{Float64, Float64, Int64, Int64}, 1}(undef, num_interior_supports * 2)

	lap_section_types = define_lap_section_types(purlin_size_span_assignment)

	lap_section_index = define_lap_section_index(lap_section_types)

	for i = 1:num_interior_supports

		lap_segments[2*i - 1] = (purlin_laps[2*i - 1]*12.0, (purlin_laps[2*i - 1]*12.0)/6, lap_section_index[i], 1)
		lap_segments[2*i] = (purlin_laps[2*i]*12.0, (purlin_laps[2*i]*12.0)/6, lap_section_index[i], 1)

	end

	return lap_segments

end

function define_purlin_line_segments(span_segments, lap_segments)

	num_spans = length(span_segments)
	
	num_purlin_line_segments = length(span_segments) + length(lap_segments)

	purlin_line_segments = Array{Tuple{Float64, Float64, Int64, Int64}, 1}(undef, num_purlin_line_segments)

	segment_index = 1

	lap_segment_index = 1	

	for i = 1:num_spans

		if i == 1 #first span
	
			purlin_line_segments[i] = span_segments[1]

			if num_spans > 1
				purlin_line_segments[i+1] = lap_segments[1]
	
				lap_segment_index = lap_segment_index + 1
				segment_index = segment_index + 2
			end
	
		elseif (i > 1) & (i != num_spans) #interior span
	
			purlin_line_segments[segment_index] = lap_segments[lap_segment_index]
			purlin_line_segments[segment_index+1] = span_segments[i]
			purlin_line_segments[segment_index+2] = lap_segments[lap_segment_index+1]
	
			lap_segment_index = lap_segment_index + 2
			segment_index = segment_index + 3
	
		elseif i == num_spans  #end span
	
			purlin_line_segments[segment_index] = lap_segments[lap_segment_index]
			purlin_line_segments[segment_index+1] = span_segments[i]
	
		end

	end

	return purlin_line_segments

end
	
function define_purlin_line_cross_section_dimensions(purlin_line_segments, lap_section_types, purlin_data, purlin_type_1, purlin_type_2)

	purlin_line_cross_section_indices = [purlin_line_segments[i][3] for i=1:length(purlin_line_segments)]

	unique_purlin_line_cross_section_indices = sort(unique(purlin_line_cross_section_indices))

	if !isempty(lap_section_types) #for multiple spans only
	
		if isempty(findall(x->x==2, unique_purlin_line_cross_section_indices)) #if there is no second purlin type defined, then add it, it won't be used though
	
			unique_purlin_line_cross_section_indices = [unique_purlin_line_cross_section_indices[1]; 2; unique_purlin_line_cross_section_indices[2:end]]
	
		end

	end

	num_purlin_line_cross_sections = length(unique_purlin_line_cross_section_indices)

	purlin_line_cross_section_dimensions = Vector{Tuple{String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}(undef, num_purlin_line_cross_sections)

	purlin_index_1 = findfirst(==(purlin_type_1), purlin_data.section_name)
	purlin_index_2 = findfirst(==(purlin_type_2), purlin_data.section_name)
	
	purlin_type_indices = [purlin_index_1; purlin_index_2]

	for i = 1:num_purlin_line_cross_sections

		cross_section_index = unique_purlin_line_cross_section_indices[i]

		if cross_section_index <= 2

			if purlin_type_2 == "none" #not used, set cross-section equal to purlin_type_1
	
				cross_section_index = unique_purlin_line_cross_section_indices[1]
	
			end
			
			purlin_line_cross_section_dimensions[i] = tuple([purlin_data[purlin_type_indices[cross_section_index], :][i] for i=2:17]...)

		elseif cross_section_index > 2 #these are the laps

			lap_type = lap_section_types[cross_section_index - 2]

			if lap_type == "1-1"

				cross_section_index = unique_purlin_line_cross_section_indices[1]

				baseline_cross_section_dimensions = deepcopy(purlin_data[purlin_type_indices[cross_section_index], :])

				#multiply base metal thickness by 2
				baseline_cross_section_dimensions[3] = baseline_cross_section_dimensions[3] * 2.0
				
				purlin_line_cross_section_dimensions[i] = tuple([baseline_cross_section_dimensions[i] for i=2:17]...)

			elseif lap_type == "2-2"

				cross_section_index = unique_purlin_line_cross_section_indices[2]

				baseline_cross_section_dimensions = deepcopy(purlin_data[purlin_type_indices[cross_section_index], :])

				#multiply base metal thickness by 2
				baseline_cross_section_dimensions[3] = baseline_cross_section_dimensions[3] * 2.0

				purlin_line_cross_section_dimensions[i] = tuple([baseline_cross_section_dimensions[i] for i=2:17]...)

			elseif lap_type == "1-2"

				cross_section_index = unique_purlin_line_cross_section_indices[1]

				baseline_cross_section_dimensions = deepcopy(purlin_data[purlin_type_indices[cross_section_index], :])

				#add the purlin 1 and purlin 2 thicknesses together
				baseline_cross_section_dimensions[3] = baseline_cross_section_dimensions[3] + purlin_line_cross_section_dimensions[2][2]

				
				purlin_line_cross_section_dimensions[i] = tuple([baseline_cross_section_dimensions[i] for i=2:17]...)

			end

		end

	end

	return purlin_line_cross_section_dimensions

end

function define_span_segments(purlin_spans, purlin_laps, purlin_size_span_assignment)

	num_spans = length(purlin_spans)

	span_segments = Array{Tuple{Float64, Float64, Int64, Int64}, 1}(undef, num_spans)

	lap_index = 1

	for i = 1:num_spans

		if i == 1 #first spans

			if num_spans == 1 #for single span

				segment_length = purlin_spans[i]*12

			else #multiple spans
				
				segment_length = purlin_spans[i]*12 - purlin_laps[lap_index]*12

			end
			
			span_segments[i] = (segment_length, segment_length/18, purlin_size_span_assignment[i], 1)
	
			lap_index = lap_index + 1
	
		elseif (i > 1) & (i != num_spans) #interior span
	
			segment_length = purlin_spans[i]*12 - purlin_laps[lap_index]*12 - purlin_laps[lap_index+1]*12
			span_segments[i] = (segment_length, segment_length/18, purlin_size_span_assignment[i], 1)
	
			lap_index = lap_index + 2
	
		elseif i==num_spans #end span
	
			segment_length = purlin_spans[i]*12 - purlin_laps[lap_index]*12
			span_segments[i] = (segment_length, segment_length/18, purlin_size_span_assignment[i], 1)
	
		end

	end

	return span_segments

end

function roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width, purlin_type_1, purlin_type_2, purlin_size_span_assignment, design_code, Fy)

	# design_code = "AISI S100-16 ASD"

	span_segments = define_span_segments(purlin_spans, purlin_laps, purlin_size_span_assignment)

	lap_section_types = define_lap_section_types(purlin_size_span_assignment)

	lap_segments = define_lap_segments(purlin_laps, purlin_size_span_assignment)

	purlin_segments = define_purlin_line_segments(span_segments, lap_segments)

	purlin_spacing = purlin_spacing * 12.0

	roof_slope = rad2deg(atan(roof_slope))

	purlin_cross_section_dimensions = define_purlin_line_cross_section_dimensions(purlin_segments, lap_section_types, purlin_data, purlin_type_1, purlin_type_2)

	purlin_material_properties = [(29500.0, 0.30, Fy, 70.0)];  #E, ν, Fy, Fu

	deck_index = findfirst(==(existing_deck_type), existing_deck_data.deck_name)

	# existing_roof_panel_details = ("screw-fastened", existing_deck_data[deck_index, 2], existing_deck_data[deck_index, 3], existing_deck_data[deck_index, 4], existing_deck_data[deck_index, 5])

	if !ismissing(existing_deck_data.fastener_spacing[deck_index])
		existing_roof_panel_details = ("screw-fastened", existing_deck_data[deck_index, 2], existing_deck_data[deck_index, 3], existing_deck_data[deck_index, 4], existing_deck_data[deck_index, 5])
	elseif !ismissing(existing_deck_data.clip_spacing[deck_index])
		existing_roof_panel_details = ("vertical leg standing seam", existing_deck_data[deck_index, 7])
	end

	existing_roof_panel_material_properties = (29500.0, 0.30, 55.0, 70.0);  #E, ν, Fy, Fu

	support_locations = [0.0; collect(cumsum(purlin_spans .* 12.0))]

	# if purlin_frame_connection == "Clip-mounted"
		
		purlin_frame_connections = "anti-roll clip"
		
	# elseif purlin_frame_connection == "Direct"
		
		# purlin_frame_connections = "bottom flange connection"
		
	# end

	intermediate_bridging_locations = [ ]

	purlin_line = PurlinLine.build(design_code, purlin_segments, purlin_spacing, roof_slope, purlin_cross_section_dimensions, purlin_material_properties, existing_roof_panel_details, existing_roof_panel_material_properties, frame_flange_width, support_locations, purlin_frame_connections, intermediate_bridging_locations)

	#Run a gravity test.
	purlin_line_gravity = deepcopy(purlin_line)
	purlin_line_gravity.loading_direction = "gravity"
	purlin_line_gravity = PurlinLine.test(purlin_line_gravity)

    #Run an uplift test.
	purlin_line_uplift = deepcopy(purlin_line)
	purlin_line_uplift.loading_direction = "uplift"
	purlin_line_uplift = PurlinLine.test(purlin_line_uplift)

	return purlin_line_gravity, purlin_line_uplift

end