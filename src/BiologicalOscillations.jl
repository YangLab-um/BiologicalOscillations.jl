module BiologicalOscillations

using Catalyst, DifferentialEquations, ModelingToolkit, Latexify
using Statistics, DSP, Peaks, LatinHypercubeSampling, DataFrames
using Combinatorics, Luxor, Images, FileIO, Random

# Models
export elowitz_2000, guan_2008
# Simulation
export generate_parameter_sets, equilibrate_ODEs, simulate_ODEs, calculate_simulation_times
export calculate_oscillatory_status, generate_find_oscillations_output
# Feature calculation
export calculate_main_frequency, calculate_amplitude, is_ODE_oscillatory, calculate_sojourn_time_fractions_in_nullclines
# Protein interaction network
export protein_interaction_network, pin_parameters, pin_timescale, pin_parameter_sets
export pin_equilibration_times, find_pin_oscillations, pin_nodes_edges
export pin_hit_rate
# Gene regulatory network
export gene_regulatory_network, grn_parameters, grn_timescale, grn_parameter_sets
export grn_equilibration_times, find_grn_oscillations
# Network utilities
export network_permutations, is_same_network, all_network_additions, unique_network_additions
export is_directed_cycle_graph, is_same_set_of_networks, unique_cycle_addition
export unique_negative_feedback_networks, count_inputs_by_coherence, is_negative_feedback_network
export connectivity_to_binary, find_all_binary_circular_permutations, binary_to_connectivity
export calculate_node_coherence, find_all_cycles_and_types
# User input handling
export is_valid_connectivity, connectivity_string_to_matrix
# Default hyperparameters
export DEFAULT_PIN_HYPERPARAMETERS, DEFAULT_PIN_PARAMETER_LIMITS, DEFAULT_PIN_SAMPLING_SCALES
export DEFAULT_GRN_HYPERPARAMETERS, DEFAULT_GRN_PARAMETER_LIMITS, DEFAULT_GRN_SAMPLING_SCALES
export DEFAULT_SIMULATION_OUTPUT
# Visualization utilities
export draw_connectivity

include("default_hyperparameters.jl")
include("models.jl")
include("simulation.jl")
include("network_utilities.jl")
include("user_input_handling.jl")
include("feature_calculation.jl")
include("gene_regulatory_network.jl")
include("protein_interaction_network.jl")
include("draw_utilities.jl")

end
