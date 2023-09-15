using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations

# Test `protein_interaction_network`
## Error handling
@test_throws MethodError protein_interaction_network([])
@test_throws MethodError protein_interaction_network("ASD")
@test_throws DomainError protein_interaction_network(zeros((0,0)))
@test_throws DomainError protein_interaction_network([0 1 -1])
@test_throws DomainError protein_interaction_network(ones((1,2)))
@test_throws DomainError protein_interaction_network([0 1; 3 0])
## Correct reactions - Repressilator
@variables t
@species (X(t))[collect(1:3)]
@parameters α[1:3], β[1:3], γ[1:3], κ[1:3], η[1:3]

true_reactions = [Reaction(α[1]*(1.0 - X[1]), nothing, [X[1]]), 
                  Reaction(β[1], [X[1]], nothing), 
                  Reaction(α[2]*(1.0 - X[2]), nothing, [X[2]]), 
                  Reaction(β[2], [X[2]], nothing), 
                  Reaction(α[3]*(1.0 - X[3]), nothing, [X[3]]), 
                  Reaction(β[3], [X[3]], nothing), 
                  Reaction(hill(abs(X[3]),γ[1],κ[1],η[1]), [X[1]], nothing), 
                  Reaction(hill(abs(X[1]),γ[2],κ[2],η[2]), [X[2]], nothing), 
                  Reaction(hill(abs(X[2]),γ[3],κ[3],η[3]), [X[3]], nothing)]

@named true_repressilator = ReactionSystem(true_reactions, t)
generated_repressilator = protein_interaction_network([0 0 -1;-1 0 0;0 -1 0])
@test reactions(generated_repressilator) == reactions(true_repressilator)
## Correct reactions - Goodwin
true_reactions = [Reaction(α[1]*(1.0 - X[1]), nothing, [X[1]]), 
                  Reaction(β[1], [X[1]], nothing), 
                  Reaction(α[2]*(1.0 - X[2]), nothing, [X[2]]), 
                  Reaction(β[2], [X[2]], nothing), 
                  Reaction(α[3]*(1.0 - X[3]), nothing, [X[3]]), 
                  Reaction(β[3], [X[3]], nothing), 
                  Reaction((1.0 - X[1])*hill(abs(X[3]),γ[1],κ[1],η[1]), nothing, [X[1]]), 
                  Reaction((1.0 - X[2])*hill(abs(X[1]),γ[2],κ[2],η[2]), nothing, [X[2]]), 
                  Reaction(hill(abs(X[2]),γ[3],κ[3],η[3]), [X[3]], nothing)]

@named true_goodwin = ReactionSystem(true_reactions, t)
generated_goodwin = protein_interaction_network([0 0 1;1 0 0;0 -1 0])
@test reactions(generated_goodwin) == reactions(true_goodwin)

# Test `pin_nodes_edges`
connectivity_2_nodes_1_edge = [1 0; 0 0]
connectivity_2_nodes_2_edges = [1 0; -1 0]
connectivity_3_nodes_4_edges = [1 0 0; -1 1 0; 0 -1 0]
connectivity_4_nodes_16_edges = [1 1 -1 1;-1 1 1 -1;1 1 1 1;-1 -1 -1 -1]

models_to_test = [
    protein_interaction_network(connectivity_2_nodes_1_edge),
    protein_interaction_network(connectivity_2_nodes_2_edges),
    protein_interaction_network(connectivity_3_nodes_4_edges),
    protein_interaction_network(connectivity_4_nodes_16_edges),
]

true_nodes = [2, 2, 3, 4]
true_edges = [1, 2, 4, 16]

for i in eachindex(models_to_test)
    @test pin_nodes_edges(models_to_test[i])[1] == true_nodes[i]
    @test pin_nodes_edges(models_to_test[i])[2] == true_edges[i]
end

# Test `pin_parameters`
α = [1.0, 0.13, 2.3]
β = [44.3, 75.3, 0.10]
γ = [2472.4, 442.2, 4410.0]
κ = [0.27, 0.91, 0.29]
η = [3.0, 4.8, 3.3]
@test_throws MethodError pin_parameters(generated_goodwin, "ASD", β, γ, κ, η)
@test_throws MethodError pin_parameters("ASD", α, β, γ, κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, [1.0, 1.0], β, γ, κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, α, [0.1, 10.0], γ, κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, α, β, [0.1, 1.0], κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, α, β, γ, [5.0, 5.0], η)
@test_throws DomainError pin_parameters(generated_goodwin, α, β, γ, κ, [1.0, 1.0, 1.0, 1.0])

generated_parameters = pin_parameters(generated_goodwin, α, β, γ, κ, η)

@nonamespace generated_α = [generated_parameters[generated_goodwin.α[i]] for i=1:3]
@nonamespace generated_β = [generated_parameters[generated_goodwin.β[i]] for i=1:3]
@nonamespace generated_γ = [generated_parameters[generated_goodwin.γ[i]] for i=1:3]
@nonamespace generated_κ = [generated_parameters[generated_goodwin.κ[i]] for i=1:3]
@nonamespace generated_η = [generated_parameters[generated_goodwin.η[i]] for i=1:3]

@test generated_α == α
@test generated_β == β
@test generated_γ == γ
@test generated_κ == κ
@test generated_η == η

# Test `pin_timescale`
α = [1.0, 0.13, 2.3]
β = [44.3, 75.3, 0.10]
γ = [2472.4, 442.2, 4410.0]
timescale = pin_timescale(α, β, γ)
@test timescale == 1.0/0.1

# Test `pin_parameter_sets`
repressilator = protein_interaction_network([0 0 -1;-1 0 0;0 -1 0])
samples = 10
parameter_array = pin_parameter_sets(repressilator, samples)
n_parameters = length(parameters(repressilator))
@test size(parameter_array) == (samples, n_parameters)
@test all(parameter_array[:,1] .== 1.0)

# Test `pin_equilibration_times`
repressilator = protein_interaction_network([0 0 -1;-1 0 0;0 -1 0])
samples = 10
parameter_array = pin_parameter_sets(repressilator, samples)
equilibration_times = pin_equilibration_times(repressilator, parameter_array)
@test size(equilibration_times) == (samples,)
#TODO: Create a known solution and test that the equilibration times are correct

# Test `find_pin_oscillations`
samples = 5000
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
T0_hit_rate = 0.0154
pin_result_T0 = find_pin_oscillations(connectivity_T0, samples)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, pin_result_T0["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples*T0_hit_rate*0.8

samples = 1000
connectivity_T0_3 = [1 0 -1; -1 0 0; 0 -1 0]
T0_3_hit_rate = 0.075
pin_result_T0_3 = find_pin_oscillations(connectivity_T0_3, samples)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, pin_result_T0_3["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples*T0_3_hit_rate*0.8

samples = 1000
connectivity_P0_6 = [0 0 0 0 -1; -1 0 0 0 0; 1 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]
P0_6_hit_rate = 0.059
pin_result_P0_6 = find_pin_oscillations(connectivity_P0_6, samples)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, pin_result_P0_6["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples*P0_6_hit_rate*0.8

# Test `pin_hit_rate`
samples = 5000
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
T0_hit_rate = 0.0154
calculated_hit_rate = pin_hit_rate(connectivity_T0, samples; verbose=false)
@test calculated_hit_rate ≈ T0_hit_rate rtol=0.3

samples = 1800
connectivity_T0_3 = [1 0 -1; -1 0 0; 0 -1 0]
T0_3_hit_rate = 0.075
calculated_hit_rate = pin_hit_rate(connectivity_T0_3, samples; verbose=false)
@test calculated_hit_rate ≈ T0_3_hit_rate rtol=0.3

samples = 2000
connectivity_P0_6 = [0 0 0 0 -1; -1 0 0 0 0; 1 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]
P0_6_hit_rate = 0.059
calculated_hit_rate = pin_hit_rate(connectivity_P0_6, samples; verbose=false)
@test calculated_hit_rate ≈ P0_6_hit_rate rtol=0.3