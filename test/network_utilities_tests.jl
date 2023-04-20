using BiologicalOscillations

# Test that the correct network permutations are generated
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