using Test

@testset "Protein Interaction Networks" begin include("protein_interaction_network/pin_tests.jl") end
@testset "Gene Regulatory Networks" begin include("gene_regulatory_network/grn_tests.jl") end