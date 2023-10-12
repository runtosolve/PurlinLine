#+echo=false
using CSV, DataFrames

#+echo=false
purlin_data = CSV.read("/Users/crismoen/.julia/dev/PurlinLine/database/Purlins.csv", DataFrame);

deck_data = CSV.read("/Users/crismoen/.julia/dev/PurlinLine/database/Existing_Deck.csv", DataFrame);



using PurlinLine

purlin_types = ["Z8x2.5 060"]

purlin_spans = ones(Float64, 4) .* 25.0  #ft

purlin_size_span_assignment = ones(Int, 8)

purlin_laps = ones(Float64, 3*2) .* 18.0/12

purlin_spacing = 5.0  #ft

frame_flange_width = 6.0  #in

roof_slope = 1.0/12

deck_type = "SSR Ultra-Dek 18 in. 24 ga";

purlin_line_gravity = UI.calculate_response(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, deck_type, deck_data, frame_flange_width, purlin_types, purlin_size_span_assignment);
