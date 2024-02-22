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
parameter_set = pin_parameter_sets(model, samples, random_seed; dimensionless_time=dimensionless_time)

perturbation_percentage = 0.01
perturbed_parameter_set = create_random_parameter_set_perturbation(parameter_set, perturbation_percentage,
                                                                   random_seed)     
for i in 1:samples
    # Only one parameter should be different and the difference should be the perturbation percentage
    @test sum(perturbed_parameter_set[i, :] .!= parameter_set[i, :]) == 1
    diff = abs(perturbed_parameter_set[i, :] - parameter_set[i, :])
    @test diff / parameter_set[i, :] == perturbation_percentage
end

# TODO: check that keep_constant parameters are respected