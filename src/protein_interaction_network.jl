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

    limits = []
    for i=1:N
        push!(limits, parameter_limits["α"])
        push!(limits, parameter_limits["β"])
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
    end
    for i=1:E
        push!(scales, sampling_scales["γ"])
        push!(scales, sampling_scales["κ"])
        push!(scales, sampling_scales["η"])
    end

    parameter_array = generate_parameter_sets(samples, limits, scales)

    if dimensionless_time
        parameter_array[:,1] .= 1.0
    end
    
    return parameter_array
end


"""
    pin_equilibration_times(model::ReactionSystem, parameter_sets::AbstractArray; equilibration_time_multiplier=10)

Calculates the equilibration times for each dataset

# Arguments (Required)
- `model::ReactionSystem`: Model generated with [`protein_interaction_network`](@ref)
- `parameter_sets::AbstractArray`: Parameter sets generated with [`pin_parameter_sets`](@ref)

# Arguments (Optional)
- `equilibration_time_multiplier::Real`: Factor by which the slowest timescale in the parameter set is multiplied by in order to calculate the final equilibration time

# Returns
- `equilibration_times::Array{Float64}`: Array of equilibration times
"""
function pin_equilibration_times(model::ReactionSystem, parameter_sets::AbstractArray; equilibration_time_multiplier=10)
    N, E = pin_nodes_edges(model)

    equilibration_times = Array{Float64}(undef, 0)
    for i in axes(parameter_sets, 1)
        p = parameter_sets[i, :]
        α = p[1:2:2N]
        β = p[2:2:2N]
        γ = p[2N+1:3:2N+3E]
        
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
- `simulation_times::Array{Float64}`: Array of simulation times
"""
function pin_simulation_times(equilibration_data::Dict, equilibration_times::AbstractVector; simulation_time_multiplier=10)
    simulation_times = Array{Float64}(undef, 0)

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
- `oscillatory_status::Array{Bool}`: Boolean array of oscillatory status
"""
function pin_oscillatory_status(simulation_data::Dict)
    oscillatory_status = Array{Bool}(undef, 0)

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
pin_result = Dict("model" => "ReactionSystem of the protein interaction network",
                  "parameter_sets" => "Dataframe of parameter sets",
                  "equilibration_result" => "Dataframe containing result metrics for the equilibration of parameter sets",
                  "simulation_result" => "Dataframe containing result metrics for the simulated parameter sets",
                  "non_oscillatory_time_series" => "Dictionary with parameter indexes as keys and TimeSeries as values",
                  "oscillatory_ode_solutions" => "Dictionary with parameter indexes as keys and ODESolutions as values",)
```

# Notes
- Non oscillatory time series are encoded as a matrix where the first column corresponds to time and the rest correspond to the values of each variable
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
    simulation_data = simulate_ODEs(model, parameter_sets[filter,:], equilibration_data["final_state"][filter], simulation_times[filter])
    # Check for oscillations
    oscillatory_status = pin_oscillatory_status(simulation_data)

    # Create a dataframe with the parameter sets
    parameter_map = paramsmap(model)
    parameter_names = Array{String}(undef, length(parameter_map))
    for (k, v) in parameter_map
        parameter_names[v] = string(k)
    end

    # Create a dataframe with the equilibration result
    equilibration_result = Dict(
        "parameter_index" => collect(1:1:samples),
        "equilibration_times" => equilibration_times,
        "final_velocity" => equilibration_data["final_velocity"],
        "frequency" => equilibration_data["frequency"],
        "is_steady_state" => .!filter,)
    final_state = mapreduce(permutedims, vcat, equilibration_data["final_state"])
    for n=1:nodes
        equilibration_result["final_state_$(n)"] = final_state[:,n]
    end

    # Create a dataframe with the simulation result
    simulation_result = Dict(
        "parameter_index" => collect(1:1:samples)[filter],
        "simulation_times" => simulation_times[filter],
        "is_oscillatory" => oscillatory_status,)
    final_state = mapreduce(permutedims, vcat, simulation_data["final_state"])
    for n=1:nodes
        frequency = Array{Float64}(undef, 0)
        power = Array{Float64}(undef, 0)
        amplitude = Array{Float64}(undef, 0)
        peak_variation = Array{Float64}(undef, 0)
        trough_variation = Array{Float64}(undef, 0)
        for i=1:sum(filter)
            push!(frequency, simulation_data["frequency_data"][i]["frequency"][n])
            push!(power, simulation_data["frequency_data"][i]["power"][n])
            push!(amplitude, simulation_data["amplitude_data"][i]["amplitude"][n])
            push!(peak_variation, simulation_data["amplitude_data"][i]["peak_variation"][n])
            push!(trough_variation, simulation_data["amplitude_data"][i]["trough_variation"][n])
        end
        simulation_result["final_state_$(n)"] = final_state[:,n]
        simulation_result["frequency_$(n)"] = frequency
        simulation_result["fft_power_$(n)"] = power
        simulation_result["amplitude_$(n)"] = amplitude
        simulation_result["peak_variation_$(n)"] = peak_variation
        simulation_result["trough_variation_$(n)"] = trough_variation
    end

    pin_result = Dict("model" => model,
                      "parameter_sets" => DataFrame(parameter_sets, parameter_names),
                      "equilibration_result" => DataFrame(equilibration_result),
                      "simulation_result" => DataFrame(simulation_result),)

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
    oscillatory = sum(pin_result["simulation_result"][!,"is_oscillatory"])
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
            oscillatory = sum(pin_result["simulation_result"][!,"is_oscillatory"])
            println("Oscillatory: $oscillatory, Samples $samples")
            hit_rate = oscillatory / samples
            trial += 1
        end
        return hit_rate
    end
end