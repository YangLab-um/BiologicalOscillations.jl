using BiologicalOscillations

# Test `network_permutations`
## 2 nodes
true_permutations = [
    [0 1; 2 0],
    [0 2; 1 0],
]
calculated_permutations = network_permutations([0 1; 2 0])
@test true_permutations == calculated_permutations
## 3 nodes
true_permutations = [
    [0 0 1; 2 0 0; 0 3 0],
    [0 1 0; 0 0 3; 2 0 0],
    [0 2 0; 0 0 1; 3 0 0],
    [0 0 3; 1 0 0; 0 2 0],
    [0 0 2; 3 0 0; 0 1 0],
    [0 3 0; 0 0 2; 1 0 0]
]
calculated_permutations = network_permutations([0 0 1; 2 0 0; 0 3 0])
@test true_permutations == calculated_permutations

# Test the is_same_network function
connectivity_1a = [1 0 -1; -1 0 0; 0 -1 0]
connectivity_1b = [0 0 -1; -1 1 0; 0 -1 0]
@test is_same_network(connectivity_1a, connectivity_1b) == true

connectivity_2a = [0 0 0 0 -1; 1 0 0 0 0; 0 1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0;]
connectivity_2b = [0 0 0 0 1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 1 0;]
@test is_same_network(connectivity_1a, connectivity_1b) == true

connectivity_3a = [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 -1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0;]
connectivity_3b = [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0;]
@test is_same_network(connectivity_3a, connectivity_3b) == false

# Test the all_network_additions function
base_connectivity_1 = [0 0 -1;-1 0 0; 0 -1 0]
true_additions_1 = [
    [1 0 -1;-1 0 0; 0 -1 0], [-1 0 -1;-1 0 0; 0 -1 0],
    [0 1 -1;-1 0 0; 0 -1 0], [0 -1 -1;-1 0 0; 0 -1 0],
    [0 0 -1;-1 1 0; 0 -1 0], [0 0 -1;-1 -1 0; 0 -1 0],
    [0 0 -1;-1 0 1; 0 -1 0], [0 0 -1;-1 0 -1; 0 -1 0],
    [0 0 -1;-1 0 0; 1 -1 0], [0 0 -1;-1 0 0; -1 -1 0],
    [0 0 -1;-1 0 0; 0 -1 1], [0 0 -1;-1 0 0; 0 -1 -1],
]

calculated_additions_1 = all_network_additions(base_connectivity_1, 1)
@test Set(true_additions_1) == Set(calculated_additions_1)

for n in 1:6
    true_number_of_additions = 2^n * binomial(6, n)
    calculated_number_of_additions = length(all_network_additions(base_connectivity_1, n))
    @test true_number_of_additions == calculated_number_of_additions
end

# Test the unique_network_additions function
base_connectivity = [0 0 -1;-1 0 0; 0 -1 0]
true_unique_additions = [
    [1 0 -1;-1 0 0; 0 -1 0], [-1 0 -1;-1 0 0; 0 -1 0],
    [0 1 -1;-1 0 0; 0 -1 0], [0 -1 -1;-1 0 0; 0 -1 0],
]
calculated_unique_additions = unique_network_additions(base_connectivity, 1)
# Test size is the is the same
@test length(true_unique_additions) == length(calculated_unique_additions)
for addition in true_unique_additions
    @test any([is_same_network(addition, calculated_addition) for calculated_addition in calculated_unique_additions])
end

# Test the unique_negative_feedback_networks function
true_unique_3_nodes = [
    [0 0 -1;-1 0 0; 0 -1 0], [0 0 1;1 0 0;0 -1 0]
]
calculated_unique_3_nodes = unique_negative_feedback_networks(3)

@test length(true_unique_3_nodes) == length(calculated_unique_3_nodes)
for network in true_unique_3_nodes
    @test any([is_same_network(network, calculated_network) for calculated_network in calculated_unique_3_nodes])
end

true_unique_4_nodes = [
    [0 0 0 -1; 1 0 0 0; 0 -1 0 0; 0 0 -1 0], [0 0 0 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0]
]
calculated_unique_4_nodes = unique_negative_feedback_networks(4)

@test length(true_unique_4_nodes) == length(calculated_unique_4_nodes)
for network in true_unique_4_nodes
    @test any([is_same_network(network, calculated_network) for calculated_network in calculated_unique_4_nodes])
end

true_unique_5_nodes = [
    [0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0], [0 0 0 0 -1; 1 0 0 0 0; 0 1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 -1 0 0 0; 0 0 1 0 0; 0 0 0 -1 0], [0 0 0 0 -1; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0]
]
calculated_unique_5_nodes = unique_negative_feedback_networks(5)

@test length(true_unique_5_nodes) == length(calculated_unique_5_nodes)
for network in true_unique_5_nodes
    @test any([is_same_network(network, calculated_network) for calculated_network in calculated_unique_5_nodes])
end

true_unique_6_nodes = [
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0], 
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0], 
]
calculated_unique_6_nodes = unique_negative_feedback_networks(6)

@test length(true_unique_6_nodes) == length(calculated_unique_6_nodes)
for network in true_unique_6_nodes
    @test any([is_same_network(network, calculated_network) for calculated_network in calculated_unique_6_nodes])
end

# Test the count_inputs_by_coherence function
connectivity = [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0]
input_counts = count_inputs_by_coherence(connectivity)
@test input_counts.coherent[1] == 0
@test input_counts.incoherent[1] == 0

connectivity = [0 1 1;1 0 0;0 -1 1]
input_counts = count_inputs_by_coherence(connectivity)
@test input_counts.coherent[1] == 1
@test input_counts.incoherent[1] == 1

connectivity = [-1 1 1;1 1 1;-1 1 1]
input_counts = count_inputs_by_coherence(connectivity)
@test input_counts.coherent[1] == 1
@test input_counts.incoherent[1] == 2

# Test is_negative_feedback_network
nf_networks = [
    [0 0 -1; -1 0 0; 0 -1 0],
    [0 0 1; 1 0 0; 0 -1 0],
    [0 0 0 -1; 1 0 0 0; 0 -1 0 0; 0 0 -1 0],
    [0 0 0 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0],
    [0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 -1 0 0 0; 0 0 1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0],
    [0 0 0 0 0 0 -1; -1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 -1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 -1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 -1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 -1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 -1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 -1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 -1 0 0; 0 0 0 0 0 1 0],
    [0 0 0 0 0 0 -1; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0],
]

non_nf_networks = [
    [0 0 1; -1 0 0; 0 -1 0],
    [0 0 1; 1 0 0; 0 1 0],
    [0 0 0 -1; 1 0 0 1; 0 -1 0 0; 0 0 -1 0],
    [0 0 0 -1; 1 0 1 0; 1 1 0 0; 0 0 1 0],
    [0 0 0 0 -1; -1 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 -1 0],
    [0 0 0 0 -1; 1 0 0 0 0; 0 -1 0 0 0; 0 0 1 0 0; 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 1 0 0; 0 0 -1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 -1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 1 0 0 1 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 -1 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0],
    [0 0 0 0 0 -1; 1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0],
]

for network in nf_networks
    @test is_negative_feedback_network(network) == true
end

for network in non_nf_networks
    @test is_negative_feedback_network(network) == false
end

# Test node_coherence
node_inputs = [0 0 0 0 0 0]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "unknown"

node_inputs = [0 0 1 -1 0 0]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "incoherent"

node_inputs = [0 0 1 -1 0 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "incoherent"

node_inputs = [1 0 1 -1 0 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "incoherent"

node_inputs = [1 0 1 1 0 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "coherent"

node_inputs = [1 -1 1 1 -1 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "incoherent"

node_inputs = [0 0 1 0 0 0]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence == "unknown"

# Test is_directed_cycle_graph function
@test is_directed_cycle_graph([0 1 0; 0 0 1; 1 0 0])  # forward cycle
@test is_directed_cycle_graph([0 0 1; 1 0 0; 0 1 0])  # backward cycle
@test !is_directed_cycle_graph([0 1 0; 0 0 1; 1 0 1])  # two cycles case 1
@test !is_directed_cycle_graph([0 0 1; 1 0 1; 0 1 0])  # two cycles case 2
@test !is_directed_cycle_graph([0 1 0; 0 0 1; 0 0 0])  # dangling end

# Test is_same_set_of_networks
@test is_same_set_of_networks([[1 0 1; 1 0 0; 0 1 0], [0 0 1; 1 -1 0; 0 1 0]], [[1 0 1; 1 0 0; 0 1 0], [-1 0 1; 1 0 0; 0 1 0]])  # same sets up to permutations
@test !is_same_set_of_networks([[1 0 1; 1 0 0; 0 1 0], [0 0 1; 1 -1 0; 0 1 0]], [[1 0 1; 1 0 0; 0 1 0]])  # inputs with different lengths
@test !is_same_set_of_networks([[1 0 1; 1 0 0; 0 1 0], [0 0 1; 1 1 0; 0 1 0]], [[1 0 1; 1 0 0; 0 1 0], [-1 0 1; 1 0 0; 0 1 0]])  # duplicate elements (up to permutations) in one set

# Test connectivity_to_binary and binary_to_connectivity functions
t0 = [0 0 -1; -1 0 0; 0 -1 0]
t2 = [0 0 -1; 1 0 0; 0 1 0]
s1 = [0 0 0 -1; -1 0 0 0; 0 -1 0 0; 0 0 1 0]
s3 = [0 0 0 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0]
p2 = [0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 1 0 0; 0 0 0 1 0]

@test connectivity_to_binary(t0) == "000"
@test connectivity_to_binary(t2) == "011"
@test connectivity_to_binary(s1) == "0001"
@test connectivity_to_binary(s3) == "0111"
@test connectivity_to_binary(p2) == "00011"

@test is_same_network(binary_to_connectivity("000"), t0)

@test is_same_network(binary_to_connectivity("011"), t2)
@test is_same_network(binary_to_connectivity("101"), t2)
@test is_same_network(binary_to_connectivity("110"), t2)

@test is_same_network(binary_to_connectivity("0001"), s1)
@test is_same_network(binary_to_connectivity("0010"), s1)
@test is_same_network(binary_to_connectivity("0100"), s1)
@test is_same_network(binary_to_connectivity("1000"), s1)

@test is_same_network(binary_to_connectivity("0111"), s3)
@test is_same_network(binary_to_connectivity("1011"), s3)
@test is_same_network(binary_to_connectivity("1101"), s3)
@test is_same_network(binary_to_connectivity("1110"), s3)

@test is_same_network(binary_to_connectivity("00011"), p2)
@test is_same_network(binary_to_connectivity("00110"), p2)
@test is_same_network(binary_to_connectivity("01100"), p2)
@test is_same_network(binary_to_connectivity("10001"), p2)
@test is_same_network(binary_to_connectivity("11000"), p2)

@test !is_same_network(binary_to_connectivity("00101"), p2)
@test !is_same_network(binary_to_connectivity("01001"), p2)
@test !is_same_network(binary_to_connectivity("01010"), p2)
@test !is_same_network(binary_to_connectivity("10010"), p2)
@test !is_same_network(binary_to_connectivity("10100"), p2)

# Test find_all_binary_circular_permutations function
@test Set(find_all_binary_circular_permutations("110")) == Set(["110", "101", "011"])
@test Set(find_all_binary_circular_permutations("1010")) == Set(["1010", "0101"])
@test Set(find_all_binary_circular_permutations("1010")) != Set(["1010", "0101", "0011", "0110", "1100", "1001"])

# Test the unique_network_additions and unique_cycle_addition functions
@test is_same_set_of_networks(unique_network_additions(t0, 1), unique_cycle_addition(t0))
@test is_same_set_of_networks(unique_network_additions(t2, 1), unique_cycle_addition(t2))
@test is_same_set_of_networks(unique_network_additions(s1, 1), unique_cycle_addition(s1))
@test is_same_set_of_networks(unique_network_additions(s3, 1), unique_cycle_addition(s3))
@test is_same_set_of_networks(unique_network_additions(p2, 1), unique_cycle_addition(p2))

# Test find_all_cycles_and_types. Note that the outputs of find_all_cycles_and_types are sorted by the cycle length (ascending order)
connectivity = [0 0 0 1 -1;-1 0 0 0 0;0 -1 0 0 0;0 0 -1 0 0;0 0 0 -1 0]
cycle, type = find_all_cycles_and_types(connectivity)
@test cycle[1] == [1, 2, 3, 4]
@test type[1] == "negative"

connectivity = [1 0 0 0 -1;-1 0 0 0 0;0 -1 0 0 0;0 0 -1 0 0;0 0 0 -1 0]
cycle, type = find_all_cycles_and_types(connectivity)
@test cycle[1] == [1]
@test type[1] == "positive"

connectivity = [0 0 0 0 -1;-1 0 0 1 0;0 -1 0 0 0;0 0 -1 0 0;0 0 0 -1 0]
cycle, type = find_all_cycles_and_types(connectivity)
@test cycle[1] == [2, 3, 4]
@test type[1] == "positive"

connectivity = [0 0 0 0 -1;-1 0 0 0 0;0 -1 0 1 0;0 0 -1 0 0;0 0 0 -1 0]
cycle, type = find_all_cycles_and_types(connectivity)
@test cycle[1] == [3, 4]
@test type[1] == "negative"

connectivity = [0 0 0 0 -1;-1 0 0 0 0;0 -1 0 -1 0;0 0 -1 0 0;0 0 0 -1 0]
cycle, type = find_all_cycles_and_types(connectivity)
@test cycle[1] == [3, 4]
@test type[1] == "positive"
