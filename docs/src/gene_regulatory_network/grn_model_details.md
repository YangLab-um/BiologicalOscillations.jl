# [Model details](@id grn_model_details)

```@example grn_example
using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify, Plots

connectivity_interaction = [0 0 1; 1 0 0; 0 -1 0]
model_interaction = gene_regulatory_network(connectivity_interaction)

tspan = (0.0, 0.5)
initial_condition = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
α = [1.0, 0.13, 2.3]
β = [44.3, 75.3, 0.10]
δ = [44.3, 75.3, 0.10]
γ = [2472.4, 442.2, 4410.0]
κ = [0.27, 0.91, 0.29]
η = [3.0, 4.8, 3.3]
params = grn_parameters(model_interaction, α, β, δ, γ, κ, η)

ode_problem = ODEProblem(model_interaction, initial_condition, tspan, params)
solution = solve(ode_problem, RadauIIA5())
plot(solution)
```