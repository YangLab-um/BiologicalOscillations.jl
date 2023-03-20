# [Model details](@id pin_model_details)

The [`protein_interaction_network`](@ref) functionality of BiologicalOscillations.jl stems from the work of [Tsai et al. 2008](https://doi.org/10.1126/science.1156951) and [Li et al. 2017](https://doi.org/10.1016/j.cels.2017.06.013) where a general model of interacting proteins was used to study the robustness and tunability of biological oscillators. Proteins in this model have two possible states: active or inactive, and interactions with other proteins affect the rates at which proteins convert between these two states. 

Every node has an intrinsic activation rate determined by the parameter ``\alpha`` and an inactivation rate given by ``\beta``:
```@setup pin_example
using Latexify

struct LaTeXEquation
    content::String
end

function Base.show(io::IO, ::MIME"text/latex", x::LaTeXEquation)
    # Wrap in $$ for display math printing
    return print(io, "\$\$ " * x.content * " \$\$")
end

Latexify.set_default(; starred=true)
```
```@example pin_example
using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify

connectivity_no_interaction = [0 0 0; 0 0 0; 0 0 0]
model_no_interaction = protein_interaction_network(connectivity_no_interaction)
odesys_no_interaction = convert(ODESystem, model_no_interaction)
latexify(odesys_no_interaction)
eq = latexify(odesys_no_interaction) # hide
LaTeXEquation(eq) # hide
```
In the example above we are creating a network of 3 nodes (notice the 3x3 size of connectivity) without any interaction between nodes (all zeros in the connectivity matrix). Thus, each node only contains their intrinsic activation and deactivation rates. To display the created equations we convert our `ReactionSystem` into an `ODESystem` and print the equations with `latexify()`. The model only tracks the fraction of each protein on its active state. Therefore, each variable in the differential equations (``X(t)_i``) is bounded to be between 0 and 1.

When interactions between nodes are added, additional terms affect the rate at which the fraction of active protein evolves over time. These interactions between nodes can be positive or negative and modify the equations in the following way:

```@example pin_example
using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify

connectivity_interaction = [0 0 1; 1 0 0; 0 -1 0]
model_interaction = protein_interaction_network(connectivity_interaction)
odesys_interaction = convert(ODESystem, model_interaction)
latexify(odesys_interaction)
eq = latexify(odesys_interaction) # hide
LaTeXEquation(eq) # hide
```
where
```math
\text{hill}(X(t)_i, \gamma_j, \kappa_j, \eta_j) = \gamma_j \frac{X(t)_i^{\eta_j}}{X(t)_i^{\eta_j} + \kappa_j^{\eta_j}}
```
with `i` counting over each node and `j` each edge of the network. Every interaction has three parameters: ``\gamma_j`` encodes the strength of the interaction while ``\kappa_j`` and ``\eta_j`` determine the midpoint and sensitivity of the response respectively. Here we have created a Goodwin-like oscillator by simply feeding `protein_interaction_network` with the desired connectivity.

Thanks to [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) and [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/) the model can be simulated by simply converting it to an `ODEProblem` and specifying the network parameters and initial condition. Network parameters should be specified using the function [`pin_parameters`](@ref) to ensure that values are assigned in the correct order:

```@example pin_example
using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify, Plots

connectivity_interaction = [0 0 1; 1 0 0; 0 -1 0]
model_interaction = protein_interaction_network(connectivity_interaction)

tspan = (0.0, 0.5)
initial_condition = [0.6, 0.6, 0.6]
α = [1.0, 0.13, 2.3]
β = [44.3, 75.3, 0.10]
γ = [2472.4, 442.2, 4410.0]
κ = [0.27, 0.91, 0.29]
η = [3.0, 4.8, 3.3]
params = pin_parameters(model_interaction, α, β, γ, κ, η)

ode_problem = ODEProblem(model_interaction, initial_condition, tspan, params)
solution = solve(ode_problem, RadauIIA5())
plot(solution)
```

It is recommended to set ``\alpha_1 = 1.0`` to normalize the time units to work with a dimensionless model both in time and concentration.