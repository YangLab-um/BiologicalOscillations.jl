"""
    gene_regulatory_network(connectivity::AbstractMatrix)

Creates a ReactionSystem of interacting mRNA and proteins based on the provided `connectivity`.

# Arguments
- `connectivity::AbstractMatrix`: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network.
"""
function gene_regulatory_network(connectivity::AbstractMatrix)
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
    @species (m(t))[collect(1:N)]
    @species (X(t))[collect(1:N)]
    @parameters α[1:N], β[1:N], δ[1:N], γ[1:E], κ[1:E], η[1:E]
    rxs = Reaction[]
    # Node-specific reactions for mRNA and Protein
    for i=1:size(connectivity, 1)
        # mRNA
        push!(rxs, Reaction(1 - α[i]*m[i], nothing, [m[i]]))
        # Protein
        push!(rxs, Reaction(β[i]*m[i], nothing, [X[i]]))
        push!(rxs, Reaction(δ[i], [X[i]], nothing))
    end
    # Reactions between nodes
    e = 1 
    for (i, row) in enumerate(eachrow(connectivity))
        for (j, val) in enumerate(row)
        if val == -1
            # Repression
            push!(rxs, Reaction(hillr(abs(X[j]),γ[e],κ[e],η[e]), nothing, [m[i]]))
            e += 1
        elseif val == 1
            # Activation
            push!(rxs, Reaction(hill(abs(X[j]),γ[e],κ[e],η[e]), nothing, [m[i]]))
            e += 1
        end
        end
    end
    @named model = ReactionSystem(rxs, t)
    return model
end


"""
    grn_nodes_edges(model::ReactionSystem)

Returns the number of nodes and edges of a model created with [`gene_regulatory_network`](@ref).

# Arguments
- `model::ReactionSystem`: Model generated with [`gene_regulatory_network`](@ref)

# Returns
- `AbstractVector`: Array with 2 elements representing the number of nodes and edges
"""
function grn_nodes_edges(model::ReactionSystem)
    N = Int(length(species(model)) / 2)
    P = length(parameters(model))
    E = Int((P - 3N)/3)

    return [N, E]
end


"""
    grn_parameters(model:ReactionSystem, α::AbstractVector, β::AbstractVector, δ::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)

Creates an ordered parameter vector for a model created with the [`gene_regulatory_network`](@ref) function to use in ODEProblem.

# Arguments
- `model::ReactionSystem`: Model generated with [`gene_regulatory_network`](@ref)
- `α:AbstractVector`: Vector of intrinsic mRNA degradation rates
- `β:AbstractVector`: Vector of intrinsic protein synthesis rates
- `δ:AbstractVector`: Vector of intrinsic protein degradation rates
- `γ:AbstractVector`: Vector of strengths for each network interaction
- `κ:AbstractVector`: Vector of midpoints for each network interaction
- `η:AbstractVector`: Vector of sensitivity values for each network interaction

# Note
It is assumed that parameters follow the same order as the connectivity matrix. Namely, the first node is encoded on the first row of the connectivty matrix and the first edge comes from the first nonzero element of the connectivity.
"""
function grn_parameters(model::ReactionSystem, α::AbstractVector, β::AbstractVector, δ::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)
    N, E = grn_nodes_edges(model)

    # Check that inputs are correct
    if size(α, 1) != N
        throw(DomainError(α, "α has to be a $(N)-element vector"))
    elseif  size(β, 1) != N
        throw(DomainError(β, "β has to be a $(N)-element vector"))
    elseif  size(δ, 1) != N
        throw(DomainError(δ, "δ has to be a $(N)-element vector"))
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
        @nonamespace push!(new_map_vals, (model.δ[j], δ[j]))
    end
    for j in 1:E
        @nonamespace push!(new_map_vals, (model.γ[j], γ[j]))
        @nonamespace push!(new_map_vals, (model.κ[j], κ[j]))
        @nonamespace push!(new_map_vals, (model.η[j], η[j]))
    end
    return ModelingToolkit.varmap_to_vars(Dict(new_map_vals), parameters(model))
end


"""
    grn_timescale(α::AbstractVector, δ::AbstractVector)

Calculate the slowest timescale of the gene regulatory network

# Arguments
- `α:AbstractVector`: Vector of intrinsic mRNA degradation rates
- `δ:AbstractVector`: Vector of intrinsic protein degradation rates

# Returns
- `timescale::Real`
"""
function grn_timescale(α::AbstractVector, δ::AbstractVector)
    rates = vcat(α, δ)
    return 1.0 / minimum(rates)
end


"""
    grn_parameter_sets(model::ReactionSystem, samples::Int; dimensionless_time=true, parameter_limits=Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), "δ" => (1e-2, 1e2), "γ" => (1e2, 1e4), "κ" => (0.2, 1.0), "η" => (1.0, 5.0)), sampling_scales=Dict("α" => "log", "β" => "log", "δ" => "log", "γ" => "log", "κ" => "linear", "η" => "linear"), sampling_style="lhc")

Generate a set of parameter values for a gene regulatory network model.

# Arguments (Required)
- `model::ReactionSystem`: Model generated with `gene_regulatory_network`
- `samples::Int`: Number of parameter sets to generate

# Arguments (Optional)
- `dimensionless_time::Bool=true`: If true, α₁ is set to 1.0 as well as the first N κ's with N=number of nodes. This is done to make the timescale of the system dimensionless.
- `parameter_limits::Dict`: Dict with keys "α", "β", "δ", "γ", "κ", "η" and values of tuples with lower and upper limits for each parameter. Default values are: 
```julia
parameter_limits = Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), "δ" => (1e-2, 1e2),
                        "γ" => (1e2, 1e4), "κ" => (0.2, 1.0), "η" => (1.0, 5.0))
```
- `sampling_scales::Dict`: Dict with keys "α", "β", "δ", "γ", "κ", "η" and values of strings with the sampling scale for each parameter. Default values are: 
```julia
sampling_scales = Dict("α" => "log", "β" => "log", "δ" => "log", "γ" => "log", "κ" => "linear", "η" => "linear")
```
- `sampling_style::String="lhc"`: Sampling style. Options are "lhc" for latin hypercube sampling and "random" for random sampling.
"""
function grn_parameter_sets(model::ReactionSystem, samples::Int; dimensionless_time=true, parameter_limits=Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), "δ" => (1e-2, 1e2), "γ" => (1e2, 1e4), "κ" => (0.2, 1.0), "η" => (1.0, 5.0)), sampling_scales=Dict("α" => "log", "β" => "log", "δ" => "log", "γ" => "log", "κ" => "linear", "η" => "linear"), sampling_style="lhc")
    N, E = grn_nodes_edges(model)

    α_limits = [parameter_limits["α"] for i=1:N]
    β_limits = [parameter_limits["β"] for i=1:N]
    δ_limits = [parameter_limits["δ"] for i=1:N]
    γ_limits = [parameter_limits["γ"] for i=1:E]
    κ_limits = [parameter_limits["κ"] for i=1:E]
    η_limits = [parameter_limits["η"] for i=1:E]
    parameter_limits = vcat(α_limits, β_limits, δ_limits, γ_limits, κ_limits, η_limits)

    α_scales = [sampling_scales["α"] for i=1:N]
    β_scales = [sampling_scales["β"] for i=1:N]
    δ_scales = [sampling_scales["δ"] for i=1:N]
    γ_scales = [sampling_scales["γ"] for i=1:E]
    κ_scales = [sampling_scales["κ"] for i=1:E]
    η_scales = [sampling_scales["η"] for i=1:E]
    sampling_scales = vcat(α_scales, β_scales, δ_scales, γ_scales, κ_scales, η_scales)

    parameter_array = generate_parameter_sets(samples, parameter_limits, sampling_scales)

    parameter_sets = []
    for i=1:samples
        α = parameter_array[i][1:N]
        β = parameter_array[i][N+1:2*N]
        δ = parameter_array[i][2*N+1:3*N]
        γ = parameter_array[i][3*N+1:3*N+E]
        κ = parameter_array[i][3*N+E+1:3*N+2*E]
        η = parameter_array[i][3*N+2*E+1:3*N+3*E]

        if dimensionless_time
            α[1] = 1.0
            for j in 1:N
                κ[j] = 1.0
            end
        end

        p = grn_parameter_set(model, α, β, δ, γ, κ, η)
        push!(parameter_sets, p)
    end

    return parameter_sets
end


"""
    grn_equilibration_times(model::ReactionSystem, parameter_sets::Dict; equilibration_time_multiplier=10)

Calculate the equilibration times for a group of parameter sets.

# Arguments (Required)
- `model::ReactionSystem`: Model generated with [`gene_regulatory_network`](@ref)
- `parameter_sets::Dict`: Parameter sets generated with [`grn_parameter_sets`](@ref)

# Arguments (Optional)
- `equilibration_time_multiplier::Real=10`: Multiplier for the equilibration time. The equilibration time is calculated as the timescale of the system multiplied by this value.

# Returns
- `equilibration_times::AbstractVector`: Equilibration times for each parameter set.
"""
function grn_equilibration_times(model::ReactionSystem, parameter_sets::Dict; equilibration_time_multiplier=10)
    N, ~ = grn_nodes_edges(model)

    equilibration_times = []
    for i in axes(parameter_sets, 1)
        p = parameter_sets[i]
        @nonamespace α = [p[model.α[i]] for i=1:N]
        @nonamespace δ = [p[model.δ[i]] for i=1:N]
        
        timescale = grn_timescale(α, δ)
        push!(equilibration_times, timescale * equilibration_time_multiplier)
    end
end