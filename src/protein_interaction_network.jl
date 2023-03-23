"""
    protein_interaction_network(connectivity::AbstractMatrix)

Creates a ReactionSystem of interacting proteins based on the provided `connectivity`.

# Arguments
- `connectivity::AbstractMatrix`: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network.
"""
function protein_interaction_network(connectivity::AbstractMatrix)
    # Check that input is correct
    is_valid, errmsg = is_valid_connectivity(connectivity)
    if !is_valid
        throw(DomainError(connectivity, errmsg))
    end
    # Number of nodes
    N = size(connectivity, 1) 
    # Number of edges
    E = count(!iszero, connectivity)
    @variables t
    @species (X(t))[collect(1:N)]
    @parameters α[1:N], β[1:N], γ[1:E], κ[1:E], η[1:E]
    rxs = Reaction[]
    # Node-specific reactions
    for i=1:size(connectivity, 1)
        # + Alpha * (1 - A)
        push!(rxs, Reaction(α[i]*(1 - X[i]), nothing, [X[i]]))
        # - Beta * A
        push!(rxs, Reaction(β[i], [X[i]], nothing))
    end
    # Reactions between nodes
    e = 1 
    for (i, row) in enumerate(eachrow(connectivity))
        for (j, val) in enumerate(row)
            if val == -1
            push!(rxs, Reaction(hill(abs(X[j]),γ[e],κ[e],η[e]), [X[i]], nothing))
            e += 1
            elseif val == 1
            push!(rxs, Reaction((1.0 - X[i]) * hill(abs(X[j]),γ[e],κ[e],η[e]), nothing, [X[i]]))
            e += 1
            end
        end
    end
    @named model = ReactionSystem(rxs, t)
    return model
end

"""
    pin_parameters(model:ReactionSystem, α::AbstractVector, β::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)

Creates an ordered parameter vector for a model created with the [`protein_interaction_network`](@ref) function to use in ODEProblem.

# Arguments
- `model::ReactionSystem`: Model generated with `protein_interaction_network`
- `α:AbstractVector`: Vector of intrinsic protein activation rates for each node
- `β:AbstractVector`: Vector of intrinsic protein deactivation rates for each node
- `γ:AbstractVector`: Vector of strengths for each network interaction
- `κ:AbstractVector`: Vector of midpoints for each network interaction
- `η:AbstractVector`: Vector of sensitivity values for each network interaction

# Note
It is assumed that parameters follow the same order as the connectivity matrix. Namely, the first node is encoded on the first row of the connectivty matrix and the first edge comes from the first nonzero element of the connectivity.
"""
function pin_parameters(model::ReactionSystem, α::AbstractVector, β::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)
    N = length(species(model))
    P = length(parameters(model))
    E = Int((P - 2N)/3)
    # Check that inputs are correct
    if size(α, 1) != N
        throw(DomainError(α, "α has to be a $(N)-element vector"))
    elseif  size(β, 1) != N
        throw(DomainError(β, "β has to be a $(N)-element vector"))
    elseif size(γ, 1) != E 
        throw(DomainError(γ, "γ has to be $(E)-element vector"))
    elseif size(κ, 1) != E 
        throw(DomainError(κ, "κ has to be $(E)-element vector"))
    elseif size(η, 1) != E
        throw(DomainError(η, "η has to be $(E)-element vector"))
    end
    
    new_map_vals = []
    for j in 1:N
        @nonamespace push!(new_map_vals, (model.α[j], α[j]))
        @nonamespace push!(new_map_vals, (model.β[j], β[j]))
    end
    for j in 1:E
        @nonamespace push!(new_map_vals, (model.γ[j], γ[j]))
        @nonamespace push!(new_map_vals, (model.κ[j], κ[j]))
        @nonamespace push!(new_map_vals, (model.η[j], η[j]))
    end
    return ModelingToolkit.varmap_to_vars(Dict(new_map_vals), parameters(model))
end


"""
    pin_timescale(α::AbstractVector, β::AbstractVector, γ::AbstractVector)

Calculate the slowest timescale of the protein interaction network

# Arguments
- `α:AbstractVector`: Vector of intrinsic protein activation rates for each node
- `β:AbstractVector`: Vector of intrinsic protein deactivation rates for each node
- `γ:AbstractVector`: Vector of strengths for each network interaction

# Returns
- `timescale::Real`
"""
function pin_timescale(α::AbstractVector, β::AbstractVector, γ::AbstractVector)
    rates = vcat(α, β, γ)
    return 1.0 / minimum(rates)
end