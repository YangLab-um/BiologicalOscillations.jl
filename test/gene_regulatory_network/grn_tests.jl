using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations

# Test errors for incorrect connectivity input on gene_regulatory_network
@test_throws MethodError gene_regulatory_network([])
@test_throws MethodError gene_regulatory_network("ASD")
@test_throws DomainError gene_regulatory_network(zeros((0,0)))
@test_throws DomainError gene_regulatory_network([0 1 -1])
@test_throws DomainError gene_regulatory_network(ones((1,2)))
@test_throws DomainError gene_regulatory_network([0 1; 3 0])

# Test that gene_regulatory_network creates the correct ReactionSystems
# Repressilator
@variables t
@species (m(t))[collect(1:3)]
@species (X(t))[collect(1:3)]
@parameters α[1:3], β[1:3], δ[1:3], γ[1:3], κ[1:3], η[1:3]

true_reactions = [Reaction(1 - α[1]*m[1], nothing, [m[1]]), 
                  Reaction(β[1]*m[1], nothing, [X[1]]),
                  Reaction(δ[1], [X[1]], nothing),
                  Reaction(1 - α[2]*m[2], nothing, [m[2]]), 
                  Reaction(β[2]*m[2], nothing, [X[2]]),
                  Reaction(δ[2], [X[2]], nothing),
                  Reaction(1 - α[3]*m[3], nothing, [m[3]]), 
                  Reaction(β[3]*m[3], nothing, [X[3]]),
                  Reaction(δ[3], [X[3]], nothing),
                  Reaction(hillr(abs(X[3]),γ[1],κ[1],η[1]), nothing, [m[1]]),
                  Reaction(hillr(abs(X[1]),γ[2],κ[2],η[2]), nothing, [m[2]]),
                  Reaction(hillr(abs(X[2]),γ[3],κ[3],η[3]), nothing, [m[3]]),
]

@named true_repressilator = ReactionSystem(true_reactions, t)
generated_repressilator = gene_regulatory_network([0 0 -1;-1 0 0;0 -1 0])
@test reactions(generated_repressilator) == reactions(true_repressilator)

# Goodwin
true_reactions = [Reaction(1 - α[1]*m[1], nothing, [m[1]]), 
                  Reaction(β[1]*m[1], nothing, [X[1]]),
                  Reaction(δ[1], [X[1]], nothing),
                  Reaction(1 - α[2]*m[2], nothing, [m[2]]), 
                  Reaction(β[2]*m[2], nothing, [X[2]]),
                  Reaction(δ[2], [X[2]], nothing),
                  Reaction(1 - α[3]*m[3], nothing, [m[3]]), 
                  Reaction(β[3]*m[3], nothing, [X[3]]),
                  Reaction(δ[3], [X[3]], nothing),
                  Reaction(hill(abs(X[3]),γ[1],κ[1],η[1]), nothing, [m[1]]),
                  Reaction(hill(abs(X[1]),γ[2],κ[2],η[2]), nothing, [m[2]]),
                  Reaction(hillr(abs(X[2]),γ[3],κ[3],η[3]), nothing, [m[3]]),
]

@named true_goodwin = ReactionSystem(true_reactions, t)
generated_goodwin = gene_regulatory_network([0 0 1;1 0 0;0 -1 0])
@test reactions(generated_goodwin) == reactions(true_goodwin)

# Test that the correct parameter map is generated and inputs handled correctly
α = [1.0, 0.13, 2.3]
β = [44.3, 75.3, 0.10]
δ = [4.3, 5.3, 0.20]
γ = [2472.4, 442.2, 4410.0]
κ = [0.27, 0.91, 0.29]
η = [3.0, 4.8, 3.3]
@test_throws MethodError grn_parameters(generated_goodwin, "ASD", β, δ, γ, κ, η)
@test_throws MethodError grn_parameters("ASD", α, β, δ, γ, κ, η)
@test_throws DomainError grn_parameters(generated_goodwin, [1.0, 1.0], β, δ, γ, κ, η)
@test_throws DomainError grn_parameters(generated_goodwin, α, [0.1, 10.0], δ, γ, κ, η)
@test_throws DomainError grn_parameters(generated_goodwin, α, β, [0.1, 1.0], γ, κ, η)
@test_throws DomainError grn_parameters(generated_goodwin, α, β, δ, [0.1, 1.0], κ, η)
@test_throws DomainError grn_parameters(generated_goodwin, α, β, δ, γ, [5.0, 5.0], η)
@test_throws DomainError grn_parameters(generated_goodwin, α, β, δ, γ, κ, [1.0, 1.0, 1.0, 1.0])

generated_parameters = grn_parameters(generated_goodwin, α, β, δ, γ, κ, η)
@test generated_parameters[1] == α[1]
@test generated_parameters[2] == β[1]
@test generated_parameters[3] == δ[1]
@test generated_parameters[4] == α[2]
@test generated_parameters[5] == β[2]
@test generated_parameters[6] == δ[2]
@test generated_parameters[7] == α[3]
@test generated_parameters[8] == β[3]
@test generated_parameters[9] == δ[3]
@test generated_parameters[10] == γ[1]
@test generated_parameters[11] == κ[1]
@test generated_parameters[12] == η[1]
@test generated_parameters[13] == γ[2]
@test generated_parameters[14] == κ[2]
@test generated_parameters[15] == η[2]
@test generated_parameters[16] == γ[3]
@test generated_parameters[17] == κ[3]
@test generated_parameters[18] == η[3]