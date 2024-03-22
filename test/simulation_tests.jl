using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations, Statistics

connectivity = [0 0 1; 1 0 0; 0 -1 0]
model = protein_interaction_network(connectivity)

parameter_values = [1.0, 44.3, 0.13, 75.3, 2.3, 0.10, 2472.4, 0.27, 3.0, 442.2, 0.91, 4.8, 4410.0, 0.29, 3.3]
parameter_sets = transpose(repeat(parameter_values, 1, 3))
initial_conditions = [0.5*ones(3) for i=1:3]

α = parameter_values[1:2:5]
β = parameter_values[2:2:6]
γ = parameter_values[7:3:13]

equilibration_times = [pin_timescale(α, β, γ) for i=1:3]

equilibration_data = equilibrate_ODEs(model, parameter_sets, initial_conditions, equilibration_times)

@test abs(mean(equilibration_data["frequency"]) - 17.1) / 17.1 < 0.05
@test abs(mean(equilibration_data["final_velocity"]) - 19.3) / 19.3 < 0.05

@test abs(equilibration_data["final_state"][1][1] - 0.44) / 0.44 < 0.05
@test abs(equilibration_data["final_state"][1][2] - 0.14) / 0.14 < 0.05
@test abs(equilibration_data["final_state"][1][3] - 0.02) / 0.02 < 0.05

@test abs(equilibration_data["final_state"][2][1] - 0.44) / 0.44 < 0.05
@test abs(equilibration_data["final_state"][2][2] - 0.14) / 0.14 < 0.05
@test abs(equilibration_data["final_state"][2][3] - 0.02) / 0.02 < 0.05

@test abs(equilibration_data["final_state"][3][1] - 0.44) / 0.44 < 0.05
@test abs(equilibration_data["final_state"][3][2] - 0.14) / 0.14 < 0.05
@test abs(equilibration_data["final_state"][3][3] - 0.02) / 0.02 < 0.05

# Test `create_random_parameter_set_perturbation`
samples = 1000
random_seed = 123
parameter_set = pin_parameter_sets(model, samples, random_seed)

perturbation_percentage = 0.01
perturbed_parameter_set = create_random_parameter_set_perturbation(parameter_set, perturbation_percentage,
                                                                   random_seed)     
for i in 1:samples
    # Only one parameter should be different and the difference should be the perturbation percentage
    @test sum(perturbed_parameter_set[i, :] .!= parameter_set[i, :]) == 1
    difference = abs.(perturbed_parameter_set[i, :] - parameter_set[i, :])
    change_idx = findfirst(difference .!= 0)
    @test sum(difference / parameter_set[i, change_idx]) ≈ perturbation_percentage
end

keep_constant = [1]
perturbed_parameter_set = create_random_parameter_set_perturbation(parameter_set, perturbation_percentage,
                                                                   random_seed; keep_constant=keep_constant)
for i in 1:samples
    # First parameter should be the same
    @test perturbed_parameter_set[i, keep_constant] == parameter_set[i, keep_constant]
    # Only one parameter should be different and the difference should be the perturbation percentage
    difference = abs.(perturbed_parameter_set[i, :] - parameter_set[i, :])
    change_idx = findfirst(difference .!= 0)
    @test sum(difference / parameter_set[i, change_idx]) ≈ perturbation_percentage
end

# Test `feature_change_from_perturbation`
find_oscillations_result = Dict(
    "simulation_result" => DataFrame(Dict(
        "parameter_index" => [7, 8, 9],
        "is_oscillatory" => [true, false, true],
        "frequency_1" => [1, 0.0, 1],
        "frequency_2" => [2, 0.0, 3],
        "frequency_3" => [4, 0.0, 5],
        "amplitude_1" => [0.1, 0.0, 0.1],
        "amplitude_2" => [0.2, 0.0, 0.3],
        "amplitude_3" => [0.4, 0.0, 0.5],
    ))
)
perturbation_result = Dict(
    "simulation_result" => DataFrame(Dict(
        "parameter_index" => [1, 2],
        "is_oscillatory" => [true, false],
        "frequency_1" => [2, NaN],
        "frequency_2" => [3, 3],
        "frequency_3" => [5, NaN],
        "amplitude_1" => [0.2, 0.1],
        "amplitude_2" => [0.3, NaN],
        "amplitude_3" => [0.5, 0.5],
    ))
)
feature_change = feature_change_from_perturbation(find_oscillations_result, perturbation_result)
@test feature_change[!, "parameter_index"] ≈ [7, 9]
@test feature_change[1, "frequency_change"] ≈ 3/7
@test feature_change[1, "amplitude_change"] ≈ 3/7
@test isnan(feature_change[2, "frequency_change"])
@test isnan(feature_change[2, "amplitude_change"])

# Test `calculate_perturbed_parameter_index`

samples = 1000
random_seed = 123
parameter_set = pin_parameter_sets(model, samples, random_seed)

perturbation_percentage = 0.01
perturbed_parameter_set = create_random_parameter_set_perturbation(parameter_set, perturbation_percentage,
                                                                   random_seed)
perturbed_parameter_index = calculate_perturbed_parameter_index(parameter_set, perturbed_parameter_set)
for i in 1:samples
    # Only one parameter should be different and the difference should be the perturbation percentage
    @test sum(perturbed_parameter_set[i, :] .!= parameter_set[i, :]) == 1
    difference = abs.(perturbed_parameter_set[i, :] - parameter_set[i, :])
    change_idx = findfirst(difference .!= 0)
    @test perturbed_parameter_index[i] == change_idx
end