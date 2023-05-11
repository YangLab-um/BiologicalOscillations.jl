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


"""
    is_same_network(connectivity1::AbstractMatrix, connectivity2::AbstractMatrix)

    Returns true if connectiviy1 is a permutation of connectivity2. False otherwise.

# Arguments (Required)
- `connectivity1::AbstractMatrix`: Connectivity matrix of a network
- `connectivity2::AbstractMatrix`: Connectivity matrix of a network

# Returns
- `result::Bool`: True if any of the permutations of connectivity 1 is equal to any of the permutations of connectivity2. False otherwise.
"""
function is_same_network(connectivity1::AbstractMatrix, connectivity2::AbstractMatrix)
    permutations1 = network_permutations(connectivity1)
    permutations2 = network_permutations(connectivity2)
    # Compare permutations
    for p1 in permutations1
        for p2 in permutations2
            if p1 == p2
                return true
            end
        end
    end
    return false
end


"""
    all_network_additions(connectivity::AbstractMatrix, number_of_edges::Int64)

    Returns a vector containing all possible network additions of a network connectivity matrix.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network
- `number_of_edges::Int64`: Number of edges to add to the network

# Returns
- `connectivity_vector::AbstractVector`: Vector containing all possible network additions of a network connectivity matrix
"""
function all_network_additions(connectivity::AbstractMatrix, number_of_edges::Int64)
    connectivity_vector = []
    # Additions will only occur at places where the connectivity is zero
    zero_indices = findall(connectivity .== 0)
    # Get all possible combinations of signs for the additions
    sign_combinations = reduce(vcat, collect(Iterators.product([[-1,1] for i in 1:number_of_edges]...)))
    # Find all combinations of zero indices where number_of_edges additions will happen
    zero_index_combinations = collect(combinations(zero_indices, number_of_edges))
    # Create all possible network additions
    for sign_combination in sign_combinations
        for zero_index_combination in zero_index_combinations
            # Create a copy of the connectivity matrix
            new_connectivity = copy(connectivity)
            # Add the edges
            for (i,sign) in enumerate(sign_combination)
                new_connectivity[zero_index_combination[i]] = sign
            end
            push!(connectivity_vector, new_connectivity)
        end
    end

    return connectivity_vector
end


"""
    unique_network_additions(connectivity::AbstractMatrix, number_of_edges::Int64)

    Returns a vector containing all unique network additions of a network connectivity matrix.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network
- `number_of_edges::Int64`: Number of edges to add to the network

# Returns
- `connectivity_vector::AbstractVector`: Vector containing all unique network additions of a network connectivity matrix
"""
function unique_network_additions(connectivity::AbstractMatrix, number_of_edges::Int64)
    connectivity_vector = []
    all_additions = all_network_additions(connectivity, number_of_edges)
    # Remove duplicates
    for addition in all_additions
        if !any(is_same_network(addition, c) for c in connectivity_vector)
            push!(connectivity_vector, addition)
        end
    end

    return connectivity_vector
end