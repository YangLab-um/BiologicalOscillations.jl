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
   is_directed_cycle_graph(connectivity::AbstractMatrix)

   Returns true if the given connectivity matrix represents a directed cycle graph. Types of interactions (positive or negative) are ignored.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of an interaction network

# Returns
- `result::Bool`: True if the interaction network is directed cyclic
"""
function is_directed_cycle_graph(connectivity::AbstractMatrix)
    n = size(connectivity)[1]
    reference_circular_graph = binary_to_connectivity("1"^n)
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
    if (all(sum(identity, dims=1) .== 1) && all(sum(identity, dims=2) .== 1))
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
        if node_coherence == "coherent"
            input_counts.coherent .+= 1
        elseif node_coherence == "incoherent"
            input_counts.incoherent .+= 1
        end
    end
    return input_counts
end

  
 """
   connectivity_to_binary(connectivity::AbstractMatrix)

   Returns a binary representation of a circular network.

A binary representation of a circular interaction network tracks the type of interaction (either positive or negative) around the cycle. The current implementation assigns 1 to a positive interaction and 0 to a negative one. For example, the Goodwin oscillator can be represented as `"011"`, `"101"`, or `"110"`. If the nodes are assumed to be indistinguishable, the representation is not unique. Among possible equivalent representations, the smallest number is returned; for the above example, `"011"`.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a circular interaction network

# Returns
- `feedback_loop_in_binary::String`: Binary presentation of a circular interaction network. `nothing` if the input is not a directed circular graph
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
    feedback_loop_in_binary = equivalent_binaries[i]

    return feedback_loop_in_binary
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
   find_all_binary_circular_permutations(feedback_loop_in_binary::String)

   Returns a vector containing all unique circular permutations of a given binary representation.

# Arguments (Required)
- `feedback_loop_in_binary::String`: Binary representation of a circular interaction network.

# Returns
- `result::AbstractVector`: Vector containing all unique circular permutations of given representation
"""
function find_all_binary_circular_permutations(feedback_loop_in_binary::String)
    return unique([join(circshift(split(feedback_loop_in_binary, ""), i)) for i in 1:length(feedback_loop_in_binary)])
end


"""
   binary_to_connectivity(feedback_loop_in_binary::String)

   Returns the connectivity matrix corresponding to the given binary representation of a circular interaction network.

# Arguments (Required)
- `feedback_loop_in_binary::String`: Binary representation of a circular interaction network.

# Returns
- `connectivity::AbstractMatrix`: Connectivity matrix corresponding to the binary representation of a circular interaction network
"""
function binary_to_connectivity(feedback_loop_in_binary::String)
    n = length(feedback_loop_in_binary)

    connectivity = zeros(Int64, n, n)

    for i in 1:n
        connectivity[mod(i, n) + 1, i] = (feedback_loop_in_binary[i] == '1' ? 1 : -1)
    end

    return connectivity
end


"""
    calculate_node_coherence(node_inputs::AbstractVector)

    Returns the coherence of a node given its inputs. A coherent node has all inputs with the same sign. An incoherent node has at least one input with a different sign.

# Arguments (Required)
- `node_inputs::AbstractVector`: Vector containing the inputs of a node

# Returns
- `coherence::String`: Coherence of a node. Either "coherent", "incoherent" or "unknown"
"""
function calculate_node_coherence(node_inputs::AbstractArray)
    coherence = "unknown"

    n_positive = sum(node_inputs .== 1)
    n_negative = sum(node_inputs .== -1)

    if n_positive == 0 && n_negative > 1 || n_positive > 1 && n_negative == 0
        coherence = "coherent"
    elseif n_positive > 0 && n_negative > 0
        coherence = "incoherent"
    end

    return coherence
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
                addition[1, i] = type

                push!(connectivity_vector, addition)
            end
        end
    end

    return connectivity_vector
end


"""
   find_all_cycles_and_types(connectivity::AbstractMatrix)

   Returns a vector containing all unique cycles in a connectivity matrix. Each cycle is represented by a vector of nonduplicate node indices that define a cycle.a vector containing all unique cycles that can be found in a connectivity matrix. For example, if a Goodwin oscillator with a self-loop ([1 0 1; 1 0 0; 0 -1 0]) is given as an input, vectors [[1], [1, 2, 3]] and ["positive", "negative"] will be returned.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network

# Returns
- `all_cycles::AbstractVector`: Vector containing vectors of node indices
- `all_types::AbstractVector`: Vector containing types of respective cycles in `all_cycles`
"""
function find_all_cycles_and_types(connectivity::AbstractMatrix)
    all_cycles = []
    all_types = []
    n = size(connectivity)[1]

    # Following for statements iterate over all possible cycles in a complete directed graph (self-loops allowed) of the same size as the given graph
    for cycle_length in 1:n
        # Iterate over lengths of cycles
        for c in combinations(1:n, cycle_length)
            # With a given cycle length, iterate over possible combinations of nodes
            for p in permutations(c[2:cycle_length])
                # With a given combination of nodes, iterate over its all circular permutations

                path_candidate = vcat(c[1], p, c[1])
                # Note that the start (end) point of a cycle is repeated in the above list

                index_start_node = path_candidate[1:end - 1]
                index_end_node = path_candidate[2:end]

                path_in_given_connectivity = connectivity[CartesianIndex.(index_end_node, index_start_node)]

                if all(path_in_given_connectivity .!= 0)
                    # When possible path exists in given connectivity
                    push!(all_cycles, index_start_node)
                    push!(all_types, (prod(path_in_given_connectivity) == 1 ? "positive" : "negative"))
                end
            end
        end
    end
    return all_cycles, all_types
end
