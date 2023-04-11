module BiologicalOscillations

using Catalyst, DifferentialEquations, ModelingToolkit, Latexify
using Statistics, DSP, Peaks, LatinHypercubeSampling, DataFrames
using Combinatorics

export elowitz_2000
export generate_parameter_sets, equilibrate_ODEs, simulate_ODEs
export calculate_main_frequency, calculate_amplitude, is_ODE_oscillatory
export protein_interaction_network, pin_parameters, pin_timescale, pin_parameter_sets
export pin_equilibration_times, pin_simulation_times, pin_oscillatory_status, find_pin_oscillations
export pin_hit_rate
export gene_regulatory_network, grn_parameters, grn_timescale, grn_parameter_sets
export grn_equilibration_times, grn_simulation_times, grn_oscillatory_status, find_grn_oscillations

include("models.jl")
include("simulation.jl")
include("network_utilities.jl")
include("user_input_handling.jl")
include("feature_calculation.jl")
include("gene_regulatory_network.jl")
include("protein_interaction_network.jl")

end
