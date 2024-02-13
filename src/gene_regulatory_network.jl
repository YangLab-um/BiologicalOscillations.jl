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
    @species (m(t))[1:N]
    @species (X(t))[1:N]
    @parameters α[1:N], β[1:N], δ[1:N], γ[1:E], κ[1:E], η[1:E]
    rxs = Reaction[]
    # Node-specific reactions for mRNA and Protein
    for i=1:size(connectivity, 1)
        # mRNA
        push!(rxs, Reaction(1, nothing, [m[i]]))
        push!(rxs, Reaction(α[i], [m[i]], nothing))
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
- `random_seed::Int`: Random seed for the sampling algorithm

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
function grn_parameter_sets(model::ReactionSystem, samples::Int, random_seed::Int; dimensionless_time=true, parameter_limits=Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), "δ" => (1e-2, 1e2), "γ" => (1e-2, 1e2), "κ" => (0.2, 1.0), "η" => (1.0, 5.0)), sampling_scales=Dict("α" => "log", "β" => "log", "δ" => "log", "γ" => "log", "κ" => "linear", "η" => "linear"), sampling_style="lhc")
    N, E = grn_nodes_edges(model)

    limits = []
    for i=1:N
        push!(limits, parameter_limits["α"])
        push!(limits, parameter_limits["β"])
        push!(limits, parameter_limits["δ"])
    end
    for i=1:E
        push!(limits, parameter_limits["γ"])
        push!(limits, parameter_limits["κ"])
        push!(limits, parameter_limits["η"])
    end

    scales = []
    for i=1:N
        push!(scales, sampling_scales["α"])
        push!(scales, sampling_scales["β"])
        push!(scales, sampling_scales["δ"])
    end
    for i=1:E
        push!(scales, sampling_scales["γ"])
        push!(scales, sampling_scales["κ"])
        push!(scales, sampling_scales["η"])
    end

    parameter_array = generate_parameter_sets(samples, limits, scales, random_seed; sampling_style=sampling_style)

    if dimensionless_time
        # α₁ = 1.0
        parameter_array[:,1] .= 1.0
        # κ₁ = κ₂ = ... = κₙ = 1.0
        for i=1:N
            parameter_array[:,3*N+1+3*i-2] .= 1.0
        end
    end
    
    return parameter_array
end


"""
    grn_equilibration_times(model::ReactionSystem, parameter_sets::AbstractArray; equilibration_time_multiplier=10)

Calculate the equilibration times for a group of parameter sets.

# Arguments (Required)
- `model::ReactionSystem`: Model generated with [`gene_regulatory_network`](@ref)
- `parameter_sets::AbstractArray`: Parameter sets generated with [`grn_parameter_sets`](@ref)

# Arguments (Optional)
- `equilibration_time_multiplier::Real=10`: Multiplier for the equilibration time. The equilibration time is calculated as the timescale of the system multiplied by this value.

# Returns
- `equilibration_times::AbstractVector`: Equilibration times for each parameter set.
"""
function grn_equilibration_times(model::ReactionSystem, parameter_sets::AbstractArray; equilibration_time_multiplier=10)
    N, E = grn_nodes_edges(model)

    equilibration_times = Array{Float64}(undef, 0)
    for i in axes(parameter_sets, 1)
        p = parameter_sets[i, :]
        α = p[1:3:3*N]
        δ = p[3:3:3*N]
        timescale = grn_timescale(α, δ)
        push!(equilibration_times, timescale * equilibration_time_multiplier)
    end

    return equilibration_times
end


"""
    find_grn_oscillations(connectivity::AbstractMatrix, samples::Int; initial_conditions::AbstractVector=NaN)

Find oscillatory parameter sets in a gene regulatory network.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix for the gene regulatory network
- `samples::Int`: Number of parameter sets to test

# Arguments (Optional)
- `hyperparameters::Dict`: Dictionary of hyperparameters for the algorithm. Default values are defined in [`DEFAULT_GRN_HYPERPARAMETERS`](@ref).

# Returns
- `pin_result::Dict`: A ditionary containing the results of the oscillatory parameter set search. Output is encoded as:
```julia
pin_result = Dict("model" => "ReactionSystem of the gene regulatory network",
                  "parameter_sets" => "Dataframe of parameter sets",
                  "equilibration_result" => "Dataframe containing result metrics for the equilibration of parameter sets",
                  "simulation_result" => "Dataframe containing result metrics for the simulated parameter sets",)
```
"""
function find_grn_oscillations(connectivity::AbstractMatrix, samples::Int; hyperparameters=DEFAULT_GRN_HYPERPARAMETERS)
    # Unpack hyperparameters
    random_seed = hyperparameters["random_seed"]
    initial_conditions = hyperparameters["initial_conditions"]
    dimensionless_time = hyperparameters["dimensionless_time"]
    parameter_limits = hyperparameters["parameter_limits"]
    sampling_scales = hyperparameters["sampling_scales"]
    sampling_style = hyperparameters["sampling_style"]
    equilibration_time_multiplier = hyperparameters["equilibration_time_multiplier"]
    solver = hyperparameters["solver"]
    abstol = hyperparameters["abstol"]
    reltol = hyperparameters["reltol"]
    maxiters = hyperparameters["maxiters"]
    simulation_time_multiplier = hyperparameters["simulation_time_multiplier"]
    fft_multiplier = hyperparameters["fft_multiplier"]
    freq_variation_threshold = hyperparameters["freq_variation_threshold"]
    power_threshold = hyperparameters["power_threshold"]
    amp_variation_threshold = hyperparameters["amp_variation_threshold"]

    model = gene_regulatory_network(connectivity)
    N = length(species(model))
    parameter_sets = grn_parameter_sets(model, samples, random_seed;
                                        dimensionless_time=dimensionless_time,
                                        parameter_limits=parameter_limits,
                                        sampling_scales=sampling_scales,
                                        sampling_style=sampling_style)
    equilibration_times = grn_equilibration_times(model, parameter_sets;
                                                  equilibration_time_multiplier=equilibration_time_multiplier)
    if isnan(initial_conditions)
        initial_conditions = [10.0*ones(N) for i=1:samples]
    end
    # Equilibrate
    equilibration_data = equilibrate_ODEs(model, parameter_sets, initial_conditions, equilibration_times;
                                          solver=solver, abstol=abstol, reltol=reltol, maxiters=maxiters)
    # Filter out solutions with velocity vector smaller than the mean velocity obtained from the equilibration
    velocity = equilibration_data["final_velocity"]
    cutoff = 10 ^ mean(log10.(velocity))
    filter = velocity .> cutoff
    simulation_times = calculate_simulation_times(equilibration_data, equilibration_times;
                                                  simulation_time_multiplier=simulation_time_multiplier)
    # Simulate
    simulation_data = simulate_ODEs(model, parameter_sets[filter,:], equilibration_data["final_state"][filter], simulation_times[filter];
                                    solver=solver, abstol=abstol, reltol=reltol, maxiters=maxiters, fft_multiplier=fft_multiplier)
    # Check for oscillations
    oscillatory_status = calculate_oscillatory_status(simulation_data;
                                                      freq_variation_threshold=freq_variation_threshold,
                                                      power_threshold=power_threshold,
                                                      amp_variation_threshold=amp_variation_threshold)
    # Output
    result = generate_find_oscillations_output(model, parameter_sets, equilibration_data, equilibration_times,
                                               simulation_data, simulation_times, oscillatory_status, hyperparameters)
    return result
end


"""
    grn_hit_rate(connectivity::AbstractMatrix, initial_samples::Int; target_oscillations::Int=100, max_samples::Int=1000000, max_trials::Int=5, hyperparameters=DEFAULT_GRN_HYPERPARAMETERS)

Estimates the hit rate of the oscillatory parameter set search in a gene regulatory network.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of the gene regulatory network
- `initial_samples::Int`: Number of parameter sets to sample in the first trial

# Arguments (Optional)
- `target_oscillations::Int`: Target number of oscillatory parameter sets
- `max_samples::Int`: Maximum number of samples allowed
- `max_trials::Int`: Maximum number of trials allowed
- `hyperparameters::Dict`: Dictionary of hyperparameters for the algorithm. Default values are defined in [`DEFAULT_GRN_HYPERPARAMETERS`](@ref).

# Returns
- `hit_rate::Real`: Estimated hit rate
"""
function grn_hit_rate(connectivity::AbstractMatrix, initial_samples::Int; target_oscillations::Int=100, max_samples::Int=1000000, max_trials::Int=5, hyperparameters=DEFAULT_PIN_HYPERPARAMETERS)
    grn_result = find_grn_oscillations(connectivity, initial_samples, hyperparameters=hyperparameters)
    oscillatory = sum(grn_result["simulation_result"][!,"is_oscillatory"])
    println("Oscillatory: $oscillatory, Samples $initial_samples")
    hit_rate =  oscillatory / initial_samples

    if oscillatory > target_oscillations
        return hit_rate
    else
        samples = initial_samples
        trial = 1
        while oscillatory < target_oscillations && trial < max_trials && samples < max_samples
            if oscillatory == 0
                samples = 2*samples
            else
                samples = ceil(Int, samples * (target_oscillations + 10) / oscillatory)
            end
            grn_result = find_grn_oscillations(connectivity, samples, hyperparameters=hyperparameters)
            oscillatory = sum(grn_result["simulation_result"][!,"is_oscillatory"])
            println("Oscillatory: $oscillatory, Samples $samples")
            hit_rate = oscillatory / samples
            trial += 1
        end
        return hit_rate
    end
end
