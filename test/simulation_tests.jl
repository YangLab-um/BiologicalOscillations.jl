using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations, Statistics

connectivity = [0 0 1; 1 0 0; 0 -1 0]
model = protein_interaction_network(connectivity)

parameter_values = [1.0, 44.3, 0.13, 75.3, 2.3, 0.10, 2472.4, 0.27, 3.0, 442.2, 0.91, 4.8, 4410.0, 0.29, 3.3]
parameter_sets = transpose(repeat(parameter_values, 1, 3))
initial_conditions = [0.5*ones(3) for i=1:3]
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