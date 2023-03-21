module BiologicalOscillations

using Catalyst, DifferentialEquations, ModelingToolkit, Latexify

export protein_interaction_network, pin_parameters
export gene_regulatory_network, grn_parameters
export elowitz_2000

include("models.jl")
include("user_input_handling.jl")
include("gene_regulatory_network.jl")
include("protein_interaction_network.jl")

end
