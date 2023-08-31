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
node_inputs = [0 0 1 -1 0 0]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence.coherent[1] == 0
@test node_coherence.incoherent[1] == 1

node_inputs = [0 0 1 -1 0 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence.coherent[1] == 0
@test node_coherence.incoherent[1] == 1

node_inputs = [1 0 1 -1 0 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence.coherent[1] == 1
@test node_coherence.incoherent[1] == 1

node_inputs = [1 0 1 1 0 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence.coherent[1] == 2
@test node_coherence.incoherent[1] == 0

node_inputs = [1 -1 1 1 -1 1]
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence.coherent[1] == 1
@test node_coherence.incoherent[1] == 2

node_inputs = [0 0 1 0 0 0]
node_coherence = calculate_node_coherence(node_inputs)
node_coherence = calculate_node_coherence(node_inputs)
@test node_coherence.coherent[1] == 0
@test node_coherence.incoherent[1] == 0

# Test calculate_loop_length_and_type
connectivity = [0 0 0 1 -1;-1 0 0 0 0;0 -1 0 0 0;0 0 -1 0 0;0 0 0 -1 0]
loop_start = [1, 4]
loop_properties = calculate_loop_length_and_type(connectivity, loop_start)
@test loop_properties.length[1] == 4
@test loop_properties.type[1] == "negative"

connectivity = [1 0 0 0 -1;-1 0 0 0 0;0 -1 0 0 0;0 0 -1 0 0;0 0 0 -1 0]
loop_start = [1, 1]
loop_properties = calculate_loop_length_and_type(connectivity, loop_start)
@test loop_properties.length[1] == 1
@test loop_properties.type[1] == "positive"

connectivity = [0 0 0 0 -1;-1 0 0 1 0;0 -1 0 0 0;0 0 -1 0 0;0 0 0 -1 0]
loop_start = [2, 4]
loop_properties = calculate_loop_length_and_type(connectivity, loop_start)
@test loop_properties.length[1] == 3
@test loop_properties.type[1] == "positive"

connectivity = [0 0 0 0 -1;-1 0 0 0 0;0 -1 0 1 0;0 0 -1 0 0;0 0 0 -1 0]
loop_start = [3, 4]
loop_properties = calculate_loop_length_and_type(connectivity, loop_start)
@test loop_properties.length[1] == 2
@test loop_properties.type[1] == "negative"


connectivity = [0 0 0 0 -1;-1 0 0 0 0;0 -1 0 -1 0;0 0 -1 0 0;0 0 0 -1 0]
loop_start = [3, 4]
loop_properties = calculate_loop_length_and_type(connectivity, loop_start)
@test loop_properties.length[1] == 2
@test loop_properties.type[1] == "positive"