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
   is_directed_cycle_graph(connectivity::AbstractMatrix)

   Returns true if the given connectivity matrix represents a directed cycle graph. Types of interactions (positive or negative) are ignored.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of an interaction network

# Returns
- `result::Bool`: True if the interaction network is directed cyclic
"""
function is_directed_cycle_graph(connectivity::AbstractMatrix)
    n = size(connectivity)[1]
    reference_circular_graph = zeros(Int64, n, n)    
    for i in 1:n
        reference_circular_graph[i, mod(i - 2, n) + 1] = 1
    end
    return is_same_network(abs.(connectivity), reference_circular_graph)
end


"""
    is_same_set_of_networks(connectivity_vector1::AbstractVector, connectivity_vector2::AbstractVector)

    Returns true if two sets have the same elements up to graph isomorphisms. False otherwise.

# Arguments (Required)
- `connectivity_vector1::AbstractVector`: Vector containing connectivity matrices
- `connectivity_vector2::AbstractVector`: Vector containing connectivity matrices. Must have the same length as the first vector

# Returns
- `result::Bool`: True if two sets are equal. False otherwise. This may include the case where a set is degenerate
"""
function is_same_set_of_networks(connectivity_vector1::AbstractVector, connectivity_vector2::AbstractVector)
    # Test 1: compare sizes of each set. If they are different, return false
    if length(connectivity_vector1) != length(connectivity_vector2)
        return false
    end

    # Test 2: exhaustively check if there is a 1-1 map between set elements. This automatically tests if there are redundant elements in each set
    n = length(connectivity_vector1)
    identity = zeros(Int64, n, n)
    for i in 1:n
        for j in 1:n
            if is_same_network(connectivity_vector1[i], connectivity_vector2[j])
                identity[i, j] = 1
            else
                identity[i, j] = 0
            end
        end
    end
    if (all(sum(identity, dims=1) .== 1) & all(sum(identity, dims=2) .== 1))
        # To have a 1-1 map between sets, every column and row of the identity matrix defined above must have exactly one 1 while the other entries being 0
        return true
    else
        return false
    end
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


"""
    unique_negative_feedback_networks(nodes::Int64)

    Returns a vector containing all unique negative feedback networks of a given size.

# Arguments (Required)
- `nodes::Int64`: Number of nodes in the network

# Returns
- `connectivity_vector::AbstractVector`: Vector containing all unique negative feedback networks of a given size
"""
function unique_negative_feedback_networks(nodes::Int64)
    # Initialize it as a vector of matrices with int64s
    connectivity_vector = Array{Matrix{Int64}}(undef, 0)
    # Calculate the combinations of positive and negative edges for a network of size nodes
    # If nodes is an odd number, then the possible number of positive edges is 2, 4, 6, ..., nodes-1
    # If nodes is an even number, then the possible number of positive edges is 1, 3, 5, ..., nodes-1
    if mod(nodes, 2) == 0
        possible_positive_edges = collect(1:2:nodes-1)
    else
        possible_positive_edges = collect(0:2:nodes-1)
    end
    # Iterate over number of positive edges
    for positive_edges in possible_positive_edges
        subset_connectivity_vector = Array{Matrix{Int64}}(undef, 0)
        signs = [[1 for i in 1:positive_edges]; [-1 for i in 1:nodes-positive_edges]]
        # Get all possible combinations of signs for the additions
        sign_combinations = unique(permutations(signs))
        # Create all possible networks
        for sign_combination in sign_combinations
            # Create a copy of the connectivity matrix
            new_connectivity = zeros(Int64, nodes, nodes)
            # Add the edges
            for (i,sign) in enumerate(sign_combination)
                new_connectivity[i, i] = sign
            end
            # Create negative feedback loop by shifting the rows one position up
            new_connectivity = circshift(new_connectivity, (0, -1))
            push!(subset_connectivity_vector, new_connectivity)
        end
        # Remove duplicates
        for addition in subset_connectivity_vector
            if !any(is_same_network(addition, c) for c in connectivity_vector)
                push!(connectivity_vector, addition)
            end
        end
    end

    return connectivity_vector
end


"""
   connectivity_to_binary(connectivity::AbstractMatrix)

   Returns a binary representation of a circular network.

A binary representation of a circular interaction network tracks the type of interaction (either positive or negative) around the cycle. The current implementation assigns 1 to a positive interaction and 0 to a negative one. For example, the Goodwin oscillator can be represented as `"011"`, `"101"`, or `"110"`. If the nodes are assumed to be indistinguishable, the representation is not unique. Among possible equivalent representations, the smallest number is returned; for the above example, `"011"`.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a circular interaction network

# Returns
- `String`: Binary presentation of a circular interaction network. `nothing` if the input is not a directed circular graph
"""
function connectivity_to_binary(connectivity::AbstractMatrix)
    if !is_directed_cycle_graph(connectivity)
        # Input must be a directed cycle graph to have a binary representation
        return
    end

    binary = ""
    i_next = 1
    while true
        i_curr = i_next
        i_next = findfirst(item -> item == 1, abs.(connectivity[:, i_next]))

        if connectivity[i_next, i_curr] == 1
            binary *= '1'
        elseif connectivity[i_next, i_curr] == -1
            binary *= '0'
        end

        if i_next == 1
            # Stop traversal upon return to the starting node
            break
        end
    end

    equivalent_binaries = find_all_binary_circular_permutations(binary)
    _, i = findmin(parse.(Int64, equivalent_binaries, base=2))
    # Find the smallest one from the binary numbers that represent the same network
    binary = equivalent_binaries[i]

    return binary
end


"""
   find_all_binary_circular_permutations(feedback_loop_in_binary::String)

   Returns a vector containing all unique circular permutations of a given binary representation.

# Arguments (Required)
- `String`: Binary representation of a circular interaction network.

# Returns
- `AbstractVector`: Vector containing all unique circular permutations of given representation
"""
function find_all_binary_circular_permutations(feedback_loop_in_binary::String)
    return unique([join(circshift(split(feedback_loop_in_binary, ""), i)) for i in 1:length(feedback_loop_in_binary)])
end


"""
   binary_to_connectivity(binary::String)

   Returns the connectivity matrix corresponding to the given binary representation of a circular interaction network.

# Arguments (Required)
- `String`: Binary representation of a circular interaction network.

# Returns
- `connectivity::AbstractMatrix`: Connectivity matrix corresponding to the binary representation of a circular interaction network
"""
function binary_to_connectivity(binary::String)
    n = length(binary)

    connectivity = zeros(Int64, n, n)

    for i in 1:n
        connectivity[i, mod(i, n) + 1] = (binary[i] == '1' ? 1 : -1)
    end

    return connectivity
end


"""
   unique_cycle_addition(connectivity::AbstractMatrix)

   Returns a vector containing all unique single-interaction additions to a circular network. Supports only the cases where the underlying network is circular and a single interaction is added to it.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a circular network

# Returns
- `connectivity_vector::AbstractVector` Vector containing all unique single-interaction additions to a circular network. An empty vector is returned if the given connectivity is not circular
"""
function unique_cycle_addition(connectivity::AbstractMatrix)
    connectivity_vector = []

    feedback_loop_in_binary = connectivity_to_binary(connectivity)
    binaries = find_all_binary_circular_permutations(feedback_loop_in_binary)
    # The 1st node can get an additional input from one of the 1st, 2nd, ..., n-1-th nodes

    for binary in binaries
        matrix = binary_to_connectivity(binary)

        for i in 1:length(binary) - 1
            for type in [-1, 1]
                addition = copy(matrix)
                addition[i, 1] = type

                push!(connectivity_vector, addition)
            end
        end
    end

    return connectivity_vector
end
