"""
    network_permutations(connectivity::AbstractMatrix)

    Returns a vector containing all possible node permutations of a network connectivity matrix.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network

# Returns
- `connectivity_vector::AbstractVector`: Vector containing all possible node permutations of a network connectivity matrix
"""
function network_permutations(connectivity::AbstractMatrix)
    connectivity_vector = []
    combinations = collect(permutations(1:size(connectivity,1)))
    for permutation in combinations
        permuted_connectivity = zeros(Int8, size(connectivity))
        permutation_matrix = zeros(Int8, size(connectivity))
        # Create permutation matrix
        for (i,p) in enumerate(permutation)
            permutation_matrix[i, p] = 1
        end
        # Calculate permutation
        permuted_connectivity = permutation_matrix' * connectivity * permutation_matrix
        push!(connectivity_vector, permuted_connectivity)
    end
    return connectivity_vector
end