# [Introduction](@id pin_introduction)

BiologicalOscillations.jl contains a set of functions implemented to study the oscillatory behavior of interacting proteins. Through the function [`protein_interaction_network`](@ref) it is possible to generate a system of differential equations that models the dynamics of a collection of interacting proteins. The function has a single input encoding how every node of the network is connected to each other. For example:
```julia
using BiologicalOscillations

connectivity = [0 0 -1; -1 0 0; 0 -1 0]
model = protein_interaction_network(connectivity)
```
creating a system of equations for a Repressilator-like network. See [Model details](@ref pin_model_details) for more information about the generated equations.