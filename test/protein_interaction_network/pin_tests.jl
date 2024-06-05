using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations
using DataFrames, Peaks

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
@species (X(t))[1:3]
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
parameter_array = pin_parameter_sets(repressilator, samples, 123)
n_parameters = length(parameters(repressilator))
@test size(parameter_array) == (samples, n_parameters)
@test all(parameter_array[:,1] .== 1.0)

# Test `pin_equilibration_times`
repressilator = protein_interaction_network([0 0 -1;-1 0 0;0 -1 0])
samples = 10
parameter_array = pin_parameter_sets(repressilator, samples, 123)
equilibration_times = pin_equilibration_times(repressilator, parameter_array)
@test size(equilibration_times) == (samples,)
#TODO: Create a known solution and test that the equilibration times are correct

# Test `find_pin_oscillations`
samples_T0 = 5000
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
T0_hit_rate = 0.0154
pin_result_T0 = find_pin_oscillations(connectivity_T0, samples_T0)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, pin_result_T0["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples_T0 * T0_hit_rate * 0.8

samples_T0_3 = 1000
connectivity_T0_3 = [1 0 -1; -1 0 0; 0 -1 0]
T0_3_hit_rate = 0.075
pin_result_T0_3 = find_pin_oscillations(connectivity_T0_3, samples_T0_3)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, pin_result_T0_3["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples_T0_3 * T0_3_hit_rate * 0.8

samples_P0_6 = 1000
connectivity_P0_6 = [0 0 0 0 -1; -1 0 0 0 0; 1 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0]
P0_6_hit_rate = 0.059
pin_result_P0_6 = find_pin_oscillations(connectivity_P0_6, samples_P0_6)
oscillatory_df = filter(row -> row["is_oscillatory"] == true, pin_result_P0_6["simulation_result"])
oscillatory_solutions = size(oscillatory_df, 1)
@test oscillatory_solutions > samples_P0_6 * P0_6_hit_rate * 0.8

# Test `simulate_pin_parameter_perturbations`
samples = 1000
random_seed = 123
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
pin_result = find_pin_oscillations(connectivity_T0, samples)
model = protein_interaction_network(connectivity_T0)
parameter_set = pin_parameter_sets(model, samples, random_seed)
perturbation_percentage = 0.01
keep_constant = [1]

perturbed_parameter_set = create_random_parameter_set_perturbation(parameter_set, perturbation_percentage,
                                                                   random_seed, keep_constant=keep_constant)     

perturbation_result = simulate_pin_parameter_perturbations(pin_result, perturbed_parameter_set)
original_parameter_sets = pin_result_T0["parameter_sets"]["oscillatory"]
perturbed_parameter_sets = perturbation_result["parameter_sets"]

@test size(perturbed_parameter_sets) == size(original_parameter_sets)

for i in axes(original_parameter_sets, 1)
    original = Array(original_parameter_sets[i, :])
    perturbed = Array(perturbed_parameter_sets[i, :])
    # First parameter should be the same
    @test perturbed[keep_constant] == original[keep_constant]
    # Only one parameter should be different and the difference should be the perturbation percentage
    @test sum(perturbed .!= original) == 1
    difference = abs.(perturbed - original)
    change_idx = findfirst(difference .!= 0)
    @test sum(difference / original[change_idx]) ≈ perturbation_percentage
end

@test perturbed_parameter_sets[!, "parameter_index"] == original_parameter_sets[!, "parameter_index"]
#TODO: Pick a single simulation (one parameter set) and test that the perturbation result is correct

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

# Test that parameters from run and generated with a random seed are the same
samples = 200
random_seed = 4567
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
hyperparameters = deepcopy(DEFAULT_PIN_HYPERPARAMETERS)
hyperparameters["random_seed"] = random_seed
result = find_pin_oscillations(connectivity_T0, samples; 
                               hyperparameters=hyperparameters)

oscillatory_params = result["parameter_sets"]["oscillatory"]
fixed_point_params = result["parameter_sets"]["non_oscillatory"]
all_params = vcat(oscillatory_params, fixed_point_params)
all_params = sort!(all_params, :parameter_index)
all_params = all_params[:, 2:end]
column_order = [ "α[1]", "β[1]", "α[2]", "β[2]", "α[3]", "β[3]",
                 "γ[1]", "κ[1]", "η[1]", "γ[2]", "κ[2]", "η[2]",
                 "γ[3]", "κ[3]", "η[3]" ]
all_params = select(all_params, column_order)

model = result["model"]
parameter_limits = hyperparameters["parameter_limits"]
sampling_scales = hyperparameters["sampling_scales"]
sampling_style = hyperparameters["sampling_style"]
dimensionless_time = hyperparameters["dimensionless_time"]

parameter_sets = pin_parameter_sets(model, samples, random_seed; 
                                    dimensionless_time=dimensionless_time, 
                                    parameter_limits=parameter_limits, 
                                    sampling_scales=sampling_scales, 
                                    sampling_style=sampling_style)

for (idx, row) in enumerate(eachrow(all_params))
    params_1 = parameter_sets[idx, :]
    params_2 = collect(row)
    @test params_1 == params_2
end

# Test that parameters are truly oscillatory and non-oscillatory
# Values taken from run on 10/11/2023
samples = 200
random_seed = 4567
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
hyperparameters = deepcopy(DEFAULT_PIN_HYPERPARAMETERS)
hyperparameters["random_seed"] = random_seed
solver = hyperparameters["solver"]
result = find_pin_oscillations(connectivity_T0, samples; 
                               hyperparameters=hyperparameters)
model = result["model"]
oscillatory_params = result["parameter_sets"]["oscillatory"]

peak_number_found = zeros(2)
for selected_idx in 1:2
    parameter_index = oscillatory_params[selected_idx, 1]
    final_states = select(result["simulation_result"], r"parameter_index|final_state_")
    final_states = filter(row -> row["parameter_index"] == parameter_index, final_states)
    initial_condition = collect(final_states[1, 1:end-1])
    final_time = 10
    parameters = collect(oscillatory_params[selected_idx, 2:end])
    ode_problem = ODEProblem(model, initial_condition, 
                            (0.0, final_time), parameters)
    sol = solve(ode_problem, solver, saveat = final_time/5000)
    df = DataFrame(t = sol.t, x = sol[1, :], y = sol[2, :], z = sol[3, :])
    peak_num = size(findmaxima(df.x)[1], 1)
    peak_number_found[selected_idx] = peak_num
end

expected_peaks = [13, 44]
@test peak_number_found == expected_peaks

# Test that parameters are truly non-oscillatory
non_oscillatory_params = result["parameter_sets"]["non_oscillatory"]

for selected_idx in 1:2
    parameter_index = non_oscillatory_params[selected_idx, 1]
    final_states = select(result["equilibration_result"], r"parameter_index|final_state_")
    final_states = filter(row -> row["parameter_index"] == parameter_index, final_states)
    final_condition = collect(final_states[1, 1:end-1])
    initial_condition = 0.5 * ones(3)
    final_time = 10
    parameters = collect(non_oscillatory_params[selected_idx, 2:end])
    ode_problem = ODEProblem(model, initial_condition, 
                            (0.0, final_time), parameters)
    sol = solve(ode_problem, solver, saveat = final_time/5000)
    df = DataFrame(t = sol.t, x = sol[1, :], y = sol[2, :], z = sol[3, :])
    final_state_simulation = collect(df[end, 2:end])
    @test final_state_simulation ≈ final_condition rtol=1e-4
end


# Test that output can be customized
samples = 200
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
hyperparameters = deepcopy(DEFAULT_PIN_HYPERPARAMETERS)

sim_output_config = Dict(
    "model" => true,
)
hyperparameters["simulation_output"] = sim_output_config
result = find_pin_oscillations(connectivity_T0, samples; 
                               hyperparameters=hyperparameters)
@test result["model"] isa ReactionSystem

sim_output_config = Dict(
    "hyperparameters" => true,
)
hyperparameters["simulation_output"] = sim_output_config
result = find_pin_oscillations(connectivity_T0, samples; 
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
result = find_pin_oscillations(connectivity_T0, samples; 
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
result = find_pin_oscillations(connectivity_T0, samples; 
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
result = find_pin_oscillations(connectivity_T0, samples; 
                               hyperparameters=hyperparameters)
@test result["simulation_result"] isa DataFrame
@test size(result["simulation_result"], 2) == 21

# Test correct number of samples and nodes on simulation_result
samples = 10
connectivity_T0 = [0 0 -1;-1 0 0;0 -1 0]
result = find_pin_oscillations(connectivity_T0, samples)
simulated_samples = result["samples"]
simulated_nodes = result["nodes"]
@test simulated_samples == samples
@test simulated_nodes == 3