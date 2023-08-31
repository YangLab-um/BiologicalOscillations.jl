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
    if connectivity1 == connectivity2
        return true
    end
    # If the sum of edges is different, then the networks are different
    if sum(connectivity1) != sum(connectivity2)
        return false
    end
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
        comparison = [is_same_network(addition, c) for c in connectivity_vector]
        if !any(comparison)
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
    count_inputs_by_coherence(connectivity::AbstractMatrix)

    Returns the total number of coherent and incoherent inputs for a given network.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network

# Returns
- `input_counts::DataFrame`: DataFrame containing the number of coherent and incoherent inputs for a given network
"""
function count_inputs_by_coherence(connectivity::AbstractMatrix)
    # Initialize dataframe
    input_counts = DataFrame(coherent = 0, incoherent = 0)
    # Go through each row
    for row in eachrow(connectivity)
        node_coherence = calculate_node_coherence(row)

        input_counts.coherent .+= node_coherence.coherent
        input_counts.incoherent .+= node_coherence.incoherent
    end
    return input_counts
end


"""
    is_negative_feedback_network(connectivity::AbstractMatrix)

    Returns true if the network is a negative-feedback-only network. False otherwise.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network

# Returns
- `result::Bool`: True if the network is a negative-feedback-only network. False otherwise.
"""
function is_negative_feedback_network(connectivity::AbstractMatrix)
    non_zero_edges = connectivity[connectivity .!= 0]
    # If the number of edges is larger than the number of nodes, then it is not a negative-feedback-only network
    if sum(abs.(non_zero_edges)) > size(connectivity, 1)
        return false
    # If the product of the edges is non-negative, then it is not a negative-feedback-only network
    elseif prod(non_zero_edges) >= 0
        return false
    else
        return true
    end
end


"""
    calculate_node_coherence(node_inputs::AbstractVector)

    Returns the coherence of a node given its inputs as a DataFrame

# Arguments (Required)
- `node_inputs::AbstractVector`: Vector containing the inputs of a node

# Returns
- `coherence::DataFrame`: DataFrame containing the coherence of a node
"""
function calculate_node_coherence(node_inputs::AbstractArray)
    n_positive = sum(node_inputs .== 1)
    n_negative = sum(node_inputs .== -1)

    incoherent = min(n_positive, n_negative)
    coherent = (max(n_positive, n_negative) - incoherent) / 2
    coherent = floor(Int, coherent)

    coherence = DataFrame(coherent = coherent, incoherent = incoherent)
    return coherence
end


"""
    calculate_loop_length_and_type(connectivity::AbstractMatrix, loop_start::AbstractVector)

    Returns the length and type of a feedback loop given its start node.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network
- `loop_start::AbstractVector`: Vector containing the coordinates of the start node of a feedback loop

# Returns
- `loop_properties::DataFrame`: DataFrame containing the length and type of a feedback loop
"""
function calculate_loop_length_and_type(connectivity::AbstractMatrix, loop_start::AbstractVector)
    feedback_sign = connectivity[loop_start...]
    loop_length = 1
    loop_type = "unknown"
    initial_node = loop_start[1]
    current_node = loop_start[2]
    next_node = loop_start[2]
    # Backtrace the loop
    while next_node != initial_node
        current_node = next_node
        next_node = findall(connectivity[next_node, :] .!= 0)[1]
        feedback_sign *= connectivity[current_node, next_node]
        loop_length += 1
    end

    if feedback_sign == 1
        loop_type = "positive"
    elseif feedback_sign == -1
        loop_type = "negative"
    end
    loop_properties = DataFrame(length = loop_length, type = loop_type)
    return loop_properties
end


"""
    classify_single_addition(reference_connectivity::AbstractMatrix, one_added_connectivity::AbstractMatrix)

    Returns the coherence and feedback loop type of a single addition to a reference network.

# Arguments (Required)
- `reference_connectivity::AbstractMatrix`: Connectivity matrix used as a reference for comparison
- `one_added_connectivity::AbstractMatrix`: Connectivity matrix with one addition with respect to the reference connectivity

# Returns
- `addition_properties::DataFrame`: DataFrame containing the coherence and feedback loop type of a single addition to a reference network
"""
function classify_single_addition(reference_connectivity::AbstractMatrix, one_added_connectivity::AbstractMatrix)
    added_edge_indices = findall(reference_connectivity .!= one_added_connectivity)[1]
    loop_start = [added_edge_indices[1], added_edge_indices[2]]
    loop_properties = calculate_loop_length_and_type(one_added_connectivity, loop_start)

    node_inputs = one_added_connectivity[added_edge_indices[1], :]
    node_coherence = calculate_node_coherence(node_inputs)
    if node_coherence.coherent[1] == 1 && node_coherence.incoherent[1] == 0
        loop_coherence = "coherent"
    elseif node_coherence.incoherent[1] == 1 && node_coherence.coherent[1] == 0
        loop_coherence = "incoherent"
    else
        loop_coherence = "unknown"
    end

    addition_properties = DataFrame(loop_coherence = loop_coherence, 
                                    loop_length = loop_properties.length, 
                                    loop_type = loop_properties.type)
    return addition_properties
end