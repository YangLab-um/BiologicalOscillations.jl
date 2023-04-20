using Test

@testset "Protein Interaction Networks" begin include("protein_interaction_network/pin_tests.jl") end
@testset "Gene Regulatory Networks" begin include("gene_regulatory_network/grn_tests.jl") end
@testset "Feature Calculation" begin include("feature_calculation_tests.jl") end
@testset "Network Utilities" begin include("network_utilities_tests.jl") end
@testset "Model Simulation" begin include("simulation_tests.jl") end