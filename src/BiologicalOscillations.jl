module BiologicalOscillations

using Catalyst, DifferentialEquations, ModelingToolkit, Latexify
using Statistics, DSP, Peaks

export elowitz_2000
export calculate_main_frequency, calculate_amplitude, is_ODE_oscillatory
export protein_interaction_network, pin_parameters, pin_timescale
export gene_regulatory_network, grn_parameters, grn_timescale

include("models.jl")
include("user_input_handling.jl")
include("feature_calculation.jl")
include("gene_regulatory_network.jl")
include("protein_interaction_network.jl")

end
