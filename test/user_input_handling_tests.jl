using BiologicalOscillations

# Test is_valid_connectivity
valid_connectivity = [0 0 1; 1 0 0; 0 -1 0]
result, errmsg = is_valid_connectivity(valid_connectivity)
@test result == true
@test errmsg == ""

invalid_connectivity = [0 0 1; 1 0 0; 0 -1 2]
result, errmsg = is_valid_connectivity(invalid_connectivity)
@test result == false
@test errmsg == "Only -1, 0, and 1 are allowed as connectivity values"

# Test connectivity_string_to_matrix
connectivity_string = "[0 0 1; 1 0 0; 0 -1 0]"
valid_connectivity = [0 0 1; 1 0 0; 0 -1 0]
connectivity_matrix = connectivity_string_to_matrix(connectivity_string)
@test connectivity_matrix == valid_connectivity