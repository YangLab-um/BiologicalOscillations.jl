using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations
using DataFrames, Peaks

# Test `nodes_have_valid_split_enzyme_outputs`
invalid_connectivity_1 = [1 0 1; 0 1 0; 0 0 -1]
invalid_connectivity_2 = [1 -1 -1; 0 1 0; 0 0 -1]
invalid_connectivity_3 = [1 1 -1 0; 0 1 0 1; 0 0 -1 1;0 0 0 -1]
result, errmsg = BiologicalOscillations.nodes_have_valid_split_enzyme_outputs(invalid_connectivity_1)
@test result == false
@test errmsg == "Invalid node output found. Each column of the connectivity matrix must have the same sign"
result, errmsg = BiologicalOscillations.nodes_have_valid_split_enzyme_outputs(invalid_connectivity_2)
@test result == false
@test errmsg == "Invalid node output found. Each column of the connectivity matrix must have the same sign"
result, errmsg = BiologicalOscillations.nodes_have_valid_split_enzyme_outputs(invalid_connectivity_3)
@test result == false
@test errmsg == "Invalid node output found. Each column of the connectivity matrix must have the same sign"

valid_connectivity_1 = [1 0 1; 0 1 0; 0 0 1]
valid_connectivity_2 = [1 1 -1; 0 1 0; 0 0 -1]
valid_connectivity_3 = [1 1 -1 0; 0 1 0 1; 0 0 -1 1;0 0 0 1]
result, errmsg = BiologicalOscillations.nodes_have_valid_split_enzyme_outputs(valid_connectivity_1)
@test result == true
@test errmsg == ""
result, errmsg = BiologicalOscillations.nodes_have_valid_split_enzyme_outputs(valid_connectivity_2)
@test result == true
@test errmsg == ""
result, errmsg = BiologicalOscillations.nodes_have_valid_split_enzyme_outputs(valid_connectivity_3)
@test result == true
@test errmsg == ""

# Test `needs_constitutive_kinase`
needs_kinase_1 = [0 0 1; 1 0 0; 0 1 0]
needs_kinase_2 = [0 0 1; -1 0 0; 0 -1 0]
needs_kinase_3 = [0 0 1; 1 0 0; 0 -1 0]
result = BiologicalOscillations.needs_constitutive_kinase(needs_kinase_1)
@test result == true
result = BiologicalOscillations.needs_constitutive_kinase(needs_kinase_2)
@test result == true
result = BiologicalOscillations.needs_constitutive_kinase(needs_kinase_3)
@test result == true

doesnt_need_kinase_1 = [0 0 -1; -1 0 0; 0 -1 0]
doesnt_need_kinase_2 = [1 0 -1; 1 -1 0; 0 -1 0]
doesnt_need_kinase_3 = [1 0 -1; 1 -1 0; 1 -1 -1]
result = BiologicalOscillations.needs_constitutive_kinase(doesnt_need_kinase_1)
@test result == false
result = BiologicalOscillations.needs_constitutive_kinase(doesnt_need_kinase_2)
@test result == false
result = BiologicalOscillations.needs_constitutive_kinase(doesnt_need_kinase_3)
@test result == false

# Test `split_enzyme_network`
## Error handling
@test_throws MethodError split_enzyme_network([])
@test_throws MethodError split_enzyme_network("ASD")
@test_throws DomainError split_enzyme_network(zeros((0,0)))
@test_throws DomainError split_enzyme_network([0 1 -1])
@test_throws DomainError split_enzyme_network(ones((1,2)))
@test_throws DomainError split_enzyme_network([0 1;3 -1])
## Correct reactions - Repressilator
@variables t
@species (X(t))[1:3]
@parameters k_b[1:3], k_u[1:3], n[1:3], κ_tot[1:3], η, P̃

true_reactions = []
for i=1:3
    K = X[mod(i-2,3)+1] # Get the previous species in the loop
    monomer = BiologicalOscillations.split_enzyme_monomer_concentration(X[i], K, P̃, κ_tot[i], n[i], η)
    push!(true_reactions, Reaction(k_b[i] * monomer ^ n[i], nothing, [X[i]]))
    push!(true_reactions, Reaction(k_u[i], [X[i]], nothing))
end
@named true_repressilator = ReactionSystem(true_reactions, t)
generated_repressilator = split_enzyme_network([0 0 -1;-1 0 0;0 -1 0])
@test reactions(generated_repressilator) == reactions(true_repressilator)
## Correct reactions - Goodwin + Incoherent positive feedback loop
@parameters k_b[1:3], k_u[1:3], n[1:3], κ_tot[1:3], η, P̃, K̃
true_reactions = []
monomer_1 = BiologicalOscillations.split_enzyme_monomer_concentration(X[1], K̃, X[3] + P̃, κ_tot[1], n[1], η)
push!(true_reactions, Reaction(k_b[1] * monomer_1 ^ n[1], nothing, [X[1]]))
push!(true_reactions, Reaction(k_u[1], [X[1]], nothing))
monomer_2 = BiologicalOscillations.split_enzyme_monomer_concentration(X[2], K̃, X[1] + P̃, κ_tot[2], n[2], η)
push!(true_reactions, Reaction(k_b[2] * monomer_2 ^ n[2], nothing, [X[2]]))
push!(true_reactions, Reaction(k_u[2], [X[2]], nothing))
monomer_3 = BiologicalOscillations.split_enzyme_monomer_concentration(X[3], X[2], X[3] + P̃, κ_tot[3], n[3], η)
push!(true_reactions, Reaction(k_b[3] * monomer_3 ^ n[3], nothing, [X[3]]))
push!(true_reactions, Reaction(k_u[3], [X[3]], nothing))

@named true_goodwin = ReactionSystem(true_reactions, t)
generated_goodwin = split_enzyme_network([0 0 1;1 0 0;0 -1 1])
@test reactions(generated_goodwin) == reactions(true_goodwin)

## Test `sen_timescale`
k_u = [1, 2, 3]
k_b = [4, 5, 6]
n = [7, 8, 9]
κ_tot = [10, 11, 12]
timescale = BiologicalOscillations.sen_timescale(k_u, k_b, n, κ_tot)
@test timescale == 1.0

## Test `sen_parameter_sets`
repressilator = split_enzyme_network([0 0 -1;-1 0 0;0 -1 0])
samples = 10
parameter_array = sen_parameter_sets(repressilator, samples, 123)
n_parameters = length(parameters(repressilator))
@test size(parameter_array) == (samples, n_parameters)
n_values = parameter_array[:, 3:4:4*3]
@test all(n_values .== round.(n_values)) # Check that n is an integer

## Test `sen_equilibration_times`
repressilator = split_enzyme_network([0 0 -1;-1 0 0;0 -1 0])
samples = 10
parameter_array = sen_parameter_sets(repressilator, samples, 123)
equilibration_times = sen_equilibration_times(repressilator, parameter_array)
@test size(equilibration_times) == (samples,)
# TODO: Create a known solution and test that the equilibration times are correct

## Test `find_sen_oscillations`
# TODO: Postponing until runtime is reduced
"""
samples = 5000
connectivity_T2 = [0 0 1;1 0 0;0 -1 0]
T2_hit_rate = 0.0012
sen_result_T2 = find_sen_oscillations(connectivity_T2, samples)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, sen_result_T2["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples*T2_hit_rate*0.8

samples = 5000
connectivity_T2_12 = [0 0 1; 1 0 0; 0 -1 1] 
T2_12_hit_rate = 0.0024 
sen_result_T2_12 = find_sen_oscillations(connectivity_T2_12, samples)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, sen_result_T2_12["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples*T2_12_hit_rate*0.8

samples = 1000
connectivity_P0 = [0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]
P0_hit_rate = 0.04
sen_result_P0 = find_sen_oscillations(connectivity_P0, samples)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, sen_result_P0["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples*P0_hit_rate*0.8
"""

## Test `sen_hit_rate`
# TODO: Postponing until runtime is reduced and simulations validated
"""
samples = 5000
connectivity_T2 = [0 0 1;1 0 0;0 -1 0]
T2_hit_rate = 0.0012
calculated_hit_rate = sen_hit_rate(connectivity_T2, samples; verbose=false)
@test calculated_hit_rate ≈ T2_hit_rate rtol=0.3

samples = 5000
connectivity_T2_12 = [0 0 1; 1 0 0; 0 -1 1]
T2_12_hit_rate = 0.0024
calculated_hit_rate = sen_hit_rate(connectivity_T2_12, samples; verbose=false)
@test calculated_hit_rate ≈ T2_12_hit_rate rtol=0.3

samples = 1000
connectivity_P0 = [0 0 0 0 -1; -1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]
P0_hit_rate = 0.04 
calculated_hit_rate = sen_hit_rate(connectivity_P0, samples; verbose=false)
@test calculated_hit_rate ≈ P0_hit_rate rtol=0.3
"""

## Test that parameters from run and generated with a random seed are the same
samples = 200
random_seed = 5678
connectivity_T2 = [0 0 1;1 0 0;0 -1 0]
hyperparameters = deepcopy(DEFAULT_SEN_HYPERPARAMETERS)
hyperparameters["random_seed"] = random_seed
result = find_sen_oscillations(connectivity_T2, samples; 
                               hyperparameters=hyperparameters)

oscillatory_params = result["parameter_sets"]["oscillatory"]
fixed_point_params = result["parameter_sets"]["non_oscillatory"]
all_params = vcat(oscillatory_params, fixed_point_params)
all_params = sort(all_params, :parameter_index)
all_params = all_params[:, 2:end]
column_order = ["k_b[1]", "k_u[1]", "n[1]", "κ_tot[1]",
                "k_b[2]", "k_u[2]", "n[2]", "κ_tot[2]",
                "k_b[3]", "k_u[3]", "n[3]", "κ_tot[3]",
                "η", "P̃", "K̃"]
all_params = select(all_params, column_order)

model = result["model"]
parameter_limits = hyperparameters["parameter_limits"]
sampling_scales = hyperparameters["sampling_scales"]
sampling_style = hyperparameters["sampling_style"]

parameter_sets = sen_parameter_sets(model, samples, random_seed; 
                                    parameter_limits=parameter_limits,
                                    sampling_scales=sampling_scales,
                                    sampling_style=sampling_style)

for (idx, row) in enumerate(eachrow(all_params))
    params_1 = parameter_sets[idx, :]
    params_2 = collect(row)
    @test params_1 == params_2
end

# Test that parameters are truly oscillatory and non-oscillatory
# Values taken from run on 01/29/2024
samples = 10000
random_seed = 1832
connectivity_T2_12 = [0 0 1; 1 0 0; 0 -1 1]
hyperparameters = deepcopy(DEFAULT_SEN_HYPERPARAMETERS)
hyperparameters["random_seed"] = random_seed
solver = hyperparameters["solver"]
result = find_sen_oscillations(connectivity_T2_12, samples; 
                               hyperparameters=hyperparameters)
model = result["model"]
oscillatory_params = result["parameter_sets"]["oscillatory"]

peak_number_found = zeros(5)

for selected_idx in 1:5
    parameter_index = oscillatory_params[selected_idx, 1]
    final_states = select(result["simulation_result"], r"parameter_index|final_state_")
    final_states = filter(row -> row["parameter_index"] == parameter_index, final_states)
    initial_condition = collect(final_states[1, 1:end-1])
    final_time = 50000
    params = collect(oscillatory_params[selected_idx, 2:end])
    ode_problem = ODEProblem(model, initial_condition, 
                            (0.0, final_time), params)
    sol = solve(ode_problem, solver, saveat = final_time/5000)
    df = DataFrame(t = sol.t, x = sol[1, :], y = sol[2, :], z = sol[3, :])
    peak_num = size(findmaxima(df.x)[1], 1)
    peak_number_found[selected_idx] = peak_num
end

expected_peaks = [1770.0, 185.0, 285.0, 2335.0, 19.0]
@test peak_number_found == expected_peaks
# Test that parameters are truly non-oscillatory
non_oscillatory_params = result["parameter_sets"]["non_oscillatory"]

for selected_idx in 1:5
    parameter_index = non_oscillatory_params[selected_idx, 1]
    final_states = select(result["equilibration_result"], r"parameter_index|final_state_")
    final_states = filter(row -> row["parameter_index"] == parameter_index, final_states)
    final_condition = collect(final_states[1, 1:end-1])
    initial_condition = zeros(3)
    final_time = 50000
    parameters = collect(non_oscillatory_params[selected_idx, 2:end])
    ode_problem = ODEProblem(model, initial_condition, 
                            (0.0, final_time), parameters)
    sol = solve(ode_problem, solver, saveat = final_time/5000)
    df = DataFrame(t = sol.t, x = sol[1, :], y = sol[2, :], z = sol[3, :])
    final_state_simulation = collect(df[end, 2:end])
    @test final_state_simulation ≈ final_condition rtol=1e-4
end

# Test that the output can be customized
samples = 200
connectivity_T2 = [0 0 1;1 0 0;0 -1 0]
hyperparameters = deepcopy(DEFAULT_SEN_HYPERPARAMETERS)

sim_output_config = Dict(
    "model" => true,
)
hyperparameters["simulation_output"] = sim_output_config
result = find_sen_oscillations(connectivity_T2, samples; 
                               hyperparameters=hyperparameters)
@test result["model"] isa ReactionSystem

sim_output_config = Dict(
    "hyperparameters" => true,
)
hyperparameters["simulation_output"] = sim_output_config
result = find_sen_oscillations(connectivity_T2, samples; 
                               hyperparameters=hyperparameters)
@test result["hyperparameters"] isa Dict
@test keys(result["hyperparameters"]) == keys(hyperparameters)
for key in keys(hyperparameters)
    # don't compare initial_conditions because it's NaN
    if key == "initial_conditions"
        continue
    end
    @test result["hyperparameters"][key] == hyperparameters[key]
end

sim_output_config = Dict(
    "parameter_sets" => Dict(
        "oscillatory" => true,
        "non_oscillatory" => true
    )
)
hyperparameters["simulation_output"] = sim_output_config
result = find_sen_oscillations(connectivity_T2, samples; 
                               hyperparameters=hyperparameters)
@test result["parameter_sets"] isa Dict
@test result["parameter_sets"]["oscillatory"] isa DataFrame
@test result["parameter_sets"]["non_oscillatory"] isa DataFrame

sim_output_config = Dict(
    "equilibration_result" => Dict(
        "parameter_index" => true,
        "equilibration_time" => true,
        "final_velocity" => true,
        "final_state" => true,
        "frequency" => true,
        "is_steady_state" => true
    )
)
hyperparameters["simulation_output"] = sim_output_config
result = find_sen_oscillations(connectivity_T2, samples; 
                               hyperparameters=hyperparameters)
@test result["equilibration_result"] isa DataFrame
@test size(result["equilibration_result"], 2) == 8

sim_output_config = Dict(
    "simulation_result" => Dict(
        "parameter_index" => true,
        "simulation_time" => true,
        "final_state" => true,
        "frequency" => true,
        "fft_power" => true,
        "amplitude" => true,
        "peak_variation" => true,
        "trough_variation" => true,
        "is_oscillatory" => true
    )
)

hyperparameters["simulation_output"] = sim_output_config
result = find_sen_oscillations(connectivity_T2, samples; 
                               hyperparameters=hyperparameters)
@test result["simulation_result"] isa DataFrame
@test size(result["simulation_result"], 2) == 21

# Scientific test. Check that the model in [Kimchi et al. 2020](https://www.science.org/doi/10.1126/sciadv.abc1939) produces oscillatory solutions
# TODO: Postponing until runtime is reduced and simulations validated
"""
connectivity = [1 -1;1 -1]
samples = 5000
hyperparameters = deepcopy(DEFAULT_SEN_HYPERPARAMETERS)
hyperparameters["random_seed"] = 1234
hyperparameters["parameter_limits"] = Dict(
    "k_b" => [1e-2, 1e0],
    "k_u" => [1e-3, 1e3],
    "n" => [2.0, 2.0],
    "κ_tot" => [1e-3, 1e2],
    "η" => [1e-1, 1e1],
    "P̃" => [1e-8, 1e-2],
)
hyperparameters["sampling_scales"] = Dict(
    "k_b" => "log",
    "k_u" => "log",
    "n" => "linear",
    "κ_tot" => "log",
    "η" => "log",
    "P̃" => "log",
)
result = find_sen_oscillations(connectivity, samples; 
                               hyperparameters=hyperparameters)
# TODO: Add test
"""