using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations

# Test errors for incorrect connectivity input on protein_interaction_network
@test_throws MethodError protein_interaction_network([])
@test_throws MethodError protein_interaction_network("ASD")
@test_throws DomainError protein_interaction_network(zeros((0,0)))
@test_throws DomainError protein_interaction_network([0 1 -1])
@test_throws DomainError protein_interaction_network(ones((1,2)))
@test_throws DomainError protein_interaction_network([0 1; 3 0])

# Test that protein_interaction_network creates the correct ReactionSystems
# Repressilator
@variables t
@species (X(t))[collect(1:3)]
@parameters α[1:3], β[1:3], γ[1:3], κ[1:3], η[1:3]

true_reactions = [Reaction(α[1]*(1.0 - X[1]), nothing, [X[1]]), 
                  Reaction(β[1], [X[1]], nothing), 
                  Reaction(α[2]*(1.0 - X[2]), nothing, [X[2]]), 
                  Reaction(β[2], [X[2]], nothing), 
                  Reaction(α[3]*(1.0 - X[3]), nothing, [X[3]]), 
                  Reaction(β[3], [X[3]], nothing), 
                  Reaction(hill(abs(X[3]),γ[1],κ[1],η[1]), [X[1]], nothing), 
                  Reaction(hill(abs(X[1]),γ[2],κ[2],η[2]), [X[2]], nothing), 
                  Reaction(hill(abs(X[2]),γ[3],κ[3],η[3]), [X[3]], nothing)]

@named true_repressilator = ReactionSystem(true_reactions, t)
generated_repressilator = protein_interaction_network([0 0 -1;-1 0 0;0 -1 0])
@test reactions(generated_repressilator) == reactions(true_repressilator)

# Goodwin

true_reactions = [Reaction(α[1]*(1.0 - X[1]), nothing, [X[1]]), 
                  Reaction(β[1], [X[1]], nothing), 
                  Reaction(α[2]*(1.0 - X[2]), nothing, [X[2]]), 
                  Reaction(β[2], [X[2]], nothing), 
                  Reaction(α[3]*(1.0 - X[3]), nothing, [X[3]]), 
                  Reaction(β[3], [X[3]], nothing), 
                  Reaction((1.0 - X[1])*hill(abs(X[3]),γ[1],κ[1],η[1]), nothing, [X[1]]), 
                  Reaction((1.0 - X[2])*hill(abs(X[1]),γ[2],κ[2],η[2]), nothing, [X[2]]), 
                  Reaction(hill(abs(X[2]),γ[3],κ[3],η[3]), [X[3]], nothing)]

@named true_goodwin = ReactionSystem(true_reactions, t)
generated_goodwin = protein_interaction_network([0 0 1;1 0 0;0 -1 0])
@test reactions(generated_goodwin) == reactions(true_goodwin)

# Test that the correct parameter map is generated and inputs handled correctly
α = [1.0, 0.13, 2.3]
β = [44.3, 75.3, 0.10]
γ = [2472.4, 442.2, 4410.0]
κ = [0.27, 0.91, 0.29]
η = [3.0, 4.8, 3.3]
@test_throws MethodError pin_parameters(generated_goodwin, "ASD", β, γ, κ, η)
@test_throws MethodError pin_parameters("ASD", α, β, γ, κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, [1.0, 1.0], β, γ, κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, α, [0.1, 10.0], γ, κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, α, β, [0.1, 1.0], κ, η)
@test_throws DomainError pin_parameters(generated_goodwin, α, β, γ, [5.0, 5.0], η)
@test_throws DomainError pin_parameters(generated_goodwin, α, β, γ, κ, [1.0, 1.0, 1.0, 1.0])

generated_parameters = pin_parameters(generated_goodwin, α, β, γ, κ, η)

@nonamespace generated_α = [generated_parameters[generated_goodwin.α[i]] for i=1:3]
@nonamespace generated_β = [generated_parameters[generated_goodwin.β[i]] for i=1:3]
@nonamespace generated_γ = [generated_parameters[generated_goodwin.γ[i]] for i=1:3]
@nonamespace generated_κ = [generated_parameters[generated_goodwin.κ[i]] for i=1:3]
@nonamespace generated_η = [generated_parameters[generated_goodwin.η[i]] for i=1:3]

@test generated_α == α
@test generated_β == β
@test generated_γ == γ
@test generated_κ == κ
@test generated_η == η

# Test that the correct timescale is obtained
timescale = pin_timescale(α, β, γ)

@test timescale == 1.0/0.1