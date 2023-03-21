"""
        is_valid_connectivity(connectivity::AbstractMatrix)

Check if connectivity is a 2-dimensional square matrix filled with 0, 1, and -1 values.

# Arguments
- `connectivity::AbstractMatrix`: A matrix.

# Returns
- `Bool`: True if connectivity is 2-dimensional square matrix filled with 0, 1, and -1 values. False otherwise.
- `String`: Error message to display in case of connectivity failing to meet the requirements. Empty string if connectivity is appropriate.
"""
function is_valid_connectivity(connectivity::AbstractMatrix)
    result = false
    errmsg = ""
    if isempty(connectivity)
        errmsg = "Connectivity cannot be empty"
    elseif ndims(connectivity) != 2
        errmsg = "Connectivity has to be a 2-dimensional matrix"
    elseif size(connectivity, 1) != size(connectivity, 2)
        errmsg = "Connectivity has to be a square matrix"
    elseif any([v ∉ [-1 0 1] for v in connectivity])
        errmsg = "Only -1, 0, and 1 are allowed as connectivity values"
    else
        result = true
    end
    return [result, errmsg]
end