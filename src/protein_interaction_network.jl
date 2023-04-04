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
    pin_nodes_edges(model::ReactionSystem)

Returns the node and edge count of a protein interaction network

# Arguments
- `model::ReactionSystem`: Model generated with [`protein_interaction_network`](@ref)

# Returns
- `AbstractVector`: Array with 2 elements representing the number of nodes and edges
"""
function pin_nodes_edges(model::ReactionSystem)
    N = length(species(model))
    P = length(parameters(model))
    E = Int((P - 2N)/3)

    return [N, E]
end


"""
    pin_parameters(model:ReactionSystem, α::AbstractVector, β::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)

Creates an ordered parameter vector for a model created with the [`protein_interaction_network`](@ref) function to use in ODEProblem.

# Arguments
- `model::ReactionSystem`: Model generated with [`protein_interaction_network`](@ref)
- `α:AbstractVector`: Vector of intrinsic protein activation rates for each node
- `β:AbstractVector`: Vector of intrinsic protein deactivation rates for each node
- `γ:AbstractVector`: Vector of strengths for each network interaction
- `κ:AbstractVector`: Vector of midpoints for each network interaction
- `η:AbstractVector`: Vector of sensitivity values for each network interaction

# Note
It is assumed that parameters follow the same order as the connectivity matrix. Namely, the first node is encoded on the first row of the connectivty matrix and the first edge comes from the first nonzero element of the connectivity.
"""
function pin_parameters(model::ReactionSystem, α::AbstractVector, β::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)
    N, E = pin_nodes_edges(model)
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
    # return ModelingToolkit.varmap_to_vars(Dict(new_map_vals), parameters(model))
    return Dict(new_map_vals)
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


"""
    pin_parameter_sets(model::ReactionSystem, samples::Int; dimensionless_time=true, parameter_limits=Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), "γ" => (1e2, 1e4), "κ" => (0.2, 1.0), "η" => (1.0, 5.0)), sampling_scales=Dict("α" => "log", "β" => "log", "γ" => "log", "κ" => "linear", "η" => "linear"), sampling_style="lhc")

Creates an array of parameter sets for a protein interaction network model.
    
# Arguments (Required)
- `model::ReactionSystem`: Model generated with [`protein_interaction_network`](@ref)
- `samples::Int`: Number of parameter sets to be generated

# Arguments (Optional)
- `dimensionless_time::Bool`: If `true`, α₁ is set to 1.0 for all parameter sets, making time dimensionless. Default value is `true`
- `parameter_limits::Dict`: Dictionary encoding the lower and upper limit of each parameter (in the form of a tuple). Default value is
```julia
parameter_limits = Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), 
                        "γ" => (1e2, 1e4), "κ" => (0.2, 1.0), "η" => (1.0, 5.0))
```
- `sampling_scales::Dict`: Dictionary with strings encoding the sampling scale for each individual parameter in a model. Accepted strings are "linear" and "log". Default value is
```julia
sampling_scales = Dict("α" => "log", "β" => "log", "γ" => "log", "κ" => "linear", "η" => "linear")
```
- `sampling_style::String`: Sampling style of the algorithm. Accepted strings are "lhc" and "random". Default value is "lhc".
"""
function pin_parameter_sets(model::ReactionSystem, samples::Int; dimensionless_time=true, parameter_limits=Dict("α" => (1e-2, 1e2), "β" => (1e-2, 1e2), "γ" => (1e2, 1e4), "κ" => (0.2, 1.0), "η" => (1.0, 5.0)), sampling_scales=Dict("α" => "log", "β" => "log", "γ" => "log", "κ" => "linear", "η" => "linear"), sampling_style="lhc")
    N, E = pin_nodes_edges(model)

    α_limits = [parameter_limits["α"] for i=1:N]
    β_limits = [parameter_limits["β"] for i=1:N]
    γ_limits = [parameter_limits["γ"] for i=1:E]
    κ_limits = [parameter_limits["κ"] for i=1:E]
    η_limits = [parameter_limits["η"] for i=1:E]
    parameter_limits = vcat(α_limits, β_limits, γ_limits, κ_limits, η_limits)

    α_scales = [sampling_scales["α"] for i=1:N]
    β_scales = [sampling_scales["β"] for i=1:N]
    γ_scales = [sampling_scales["γ"] for i=1:E]
    κ_scales = [sampling_scales["κ"] for i=1:E]
    η_scales = [sampling_scales["η"] for i=1:E]
    sampling_scales = vcat(α_scales, β_scales, γ_scales, κ_scales, η_scales)

    parameter_array = generate_parameter_sets(samples, parameter_limits, sampling_scales)

    parameter_sets = []
    for i=1:samples
        α = parameter_array[i][1:N]
        β = parameter_array[i][N+1:2*N]
        γ = parameter_array[i][2*N+1:2*N+E]
        κ = parameter_array[i][2*N+E+1:2*N+2*E]
        η = parameter_array[i][2*N+2*E+1:2*N+3*E]

        if dimensionless_time
            α[1] = 1.0
        end

        p = pin_parameters(model, α, β, γ, κ, η)
        push!(parameter_sets, p)
    end
    
    return parameter_sets
end


"""
    pin_equilibration_times(model::ReactionSystem, parameter_sets::AbstractVector; equilibration_time_multiplier=10)

Calculates the equilibration times for each dataset

# Arguments (Required)
- `model::ReactionSystem`: Model generated with [`protein_interaction_network`](@ref)
- `parameter_sets::AbstractVector`: Parameter sets generated with [`pin_parameter_sets`](@ref)

# Arguments (Optional)
- `equilibration_time_multiplier::Real`: Factor by which the slowest timescale in the parameter set is multiplied by in order to calculate the final equilibration time

# Returns
- `equilibration_times::AbstractVector`: Array of equilibration times
"""
function pin_equilibration_times(model::ReactionSystem, parameter_sets::AbstractVector; equilibration_time_multiplier=10)
    N, E = pin_nodes_edges(model)

    equilibration_times = []
    for i in axes(parameter_sets, 1)
        p = parameter_sets[i]
        @nonamespace α = [p[model.α[i]] for i=1:N]
        @nonamespace β = [p[model.β[i]] for i=1:N]
        @nonamespace γ = [p[model.γ[i]] for i=1:E]
        
        timescale = pin_timescale(α, β, γ)
        push!(equilibration_times, timescale * equilibration_time_multiplier)
    end

    return equilibration_times
end


"""
    pin_simulation_times(equilibration_data::Dict, equilibration_times::AbstractVector; simulation_time_multiplier=10)

Calculates the simulation times for each parameter set

# Arguments (Required)
- `equilibration_data::Dict`: Dictionary of equilibration data generated with [`equilibrate_ODEs`](@ref)
- `equilibration_times::AbstractVector`: Array of equilibration times generated with [`pin_equilibration_times`](@ref)

# Arguments (Optional)
- `simulation_time_multiplier::Real`: Factor by which the period (given by the estimated frequency from the equilibration) is multiplied by in order to calculate the final simulation time

# Returns
- `simulation_times::AbstractVector`: Array of simulation times
"""
function pin_simulation_times(equilibration_data::Dict, equilibration_times::AbstractVector; simulation_time_multiplier=10)
    simulation_times = []

    for i in axes(equilibration_data["frequency"], 1)
        freq = equilibration_data["frequency"][i]
        if isnan(freq)
            push!(simulation_times, equilibration_times[i])
        else
            push!(simulation_times, simulation_time_multiplier / freq )
        end
    end

    return simulation_times
end


"""
    pin_oscillatory_status(simulation_data::Dict)

Calculates the oscillatory status for each parameter set using the functionality of [`is_ODE_oscillatory`](@ref)

# Arguments (Required)
- `simulation_data::Dict`: Dictionary of simulation data generated with [`simulate_ODEs`](@ref)

# Returns
- `oscillatory_status::AbstractVector`: Boolean array of oscillatory status
"""
function pin_oscillatory_status(simulation_data::Dict)
    oscillatory_status = []

    for i in axes(simulation_data["frequency_data"], 1)
        freq = simulation_data["frequency_data"][i]
        amp = simulation_data["amplitude_data"][i]
        decision = is_ODE_oscillatory(freq, amp)
        push!(oscillatory_status, decision)
    end

    return oscillatory_status
end


"""
    find_pin_oscillations(connectivity::AbstractMatrix, samples::Int; initial_conditions=NaN)

Finds oscillatory parameter sets in a protein interaction network

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of the protein interaction network
- `samples::Int`: Number of parameter sets to sample

# Arguments (Optional)
- `initial_conditions::AbstractVector`: Initial conditions for the ODEs provided as a vector of vectors. If not specified, the initial conditions are set to 0.5 for all variables. 

# Returns
- `pin_result::Dict`: A ditionary containing the results of the oscillatory parameter set search. Output is encoded as:
```julia
pin_result = Dict('model' => "Model generated with protein_interaction_network()",
                  'parameter_sets' => "Parameter sets used for simulation",
                  'equilibration_times' => "Equilibration times used for simulation",
                  'equilibration_data => "Result from equilibrate_ODEs()", 
                  'velocity_cutoff' => "Cutoff used to filter equilibrated solutions",
                  'simulation_times' => "Simulation times used for simulation",
                  'simulation_data' => "Result from simulate_ODEs() on potential oscillatory solutions",
                  'oscilatory_status' => "Boolean array of oscilatory status for each parameter set that passed the equilibration filter")
```
"""
function find_pin_oscillations(connectivity::AbstractMatrix, samples::Int; initial_conditions=NaN)
    model = protein_interaction_network(connectivity)
    nodes = length(species(model))
    parameter_sets = pin_parameter_sets(model, samples)
    equilibration_times = pin_equilibration_times(model, parameter_sets)
    if isnan(initial_conditions)
        initial_conditions = [0.5*ones(nodes) for i=1:samples]
    end
    # Equilibrate
    equilibration_data = equilibrate_ODEs(model, parameter_sets, initial_conditions, equilibration_times)
    # Filter out solutions with velocity vector smaller than the mean velocity obtained from the equilibration
    velocity = equilibration_data["final_velocity"]
    cutoff = 10 ^ mean(log10.(velocity))
    filter = velocity .> cutoff
    simulation_times = pin_simulation_times(equilibration_data, equilibration_times)
    # Simulate
    simulation_data = simulate_ODEs(model, parameter_sets[filter], equilibration_data["final_state"][filter], simulation_times[filter])
    # Check for oscillations
    oscillatory_status = pin_oscillatory_status(simulation_data)
    # Return results
    pin_result = Dict("model" => model,
                      "parameter_sets" => parameter_sets,
                      "equilibration_times" => equilibration_times,
                      "equilibration_data" => equilibration_data,
                      "velocity_cutoff" => cutoff,
                      "simulation_times" => simulation_times,
                      "simulation_data" => simulation_data,
                      "oscillatory_status" => oscillatory_status)
    
    return pin_result
end


"""
    pin_hit_rate(connectivity::AbstractMatrix, initial_samples::Int; target_oscillations::Int=100, max_samples::Int=1000000, max_trials::Int=5)

Estimates the hit rate of the oscillatory parameter set search in a protein interaction network

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of the protein interaction network
- `initial_samples::Int`: Number of parameter sets to sample in the first trial

# Arguments (Optional)
- `target_oscillations::Int`: Target number of oscillatory parameter sets
- `max_samples::Int`: Maximum number of samples allowed
- `max_trials::Int`: Maximum number of trials allowed

# Returns
- `hit_rate::Real`: Estimated hit rate
"""
function pin_hit_rate(connectivity::AbstractMatrix, initial_samples::Int; target_oscillations::Int=100, max_samples::Int=1000000, max_trials::Int=5)
    pin_result = find_pin_oscillations(connectivity, initial_samples)
    oscillatory = sum(pin_result["oscillatory_status"])
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
            pin_result = find_pin_oscillations(connectivity, samples)
            oscillatory = sum(pin_result["oscillatory_status"])
            println("Oscillatory: $oscillatory, Samples $samples")
            hit_rate = oscillatory / samples
            trial += 1
        end
        return hit_rate
    end
end