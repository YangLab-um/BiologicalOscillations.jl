function split_enzyme_monomer_concentration(E, K, P, monomer_tot, n_monomers, η)
    """
    Concentration of enzyme monomers under the assumption of Michaelis-Menten dynamics plus a separation of timescales with
    phosphorylation/dephosphorylation reactions equilibrating much faster than split enzyme self-assembly.
    Expression taken from [Kimchi et al. 2020](https://www.science.org/doi/10.1126/sciadv.abc1939)

    # Arguments
    - `E`: Concentration of full enzyme
    - `K`: Concentration of kinase sequestering monomers
    - `P`: Concentration of phosphatase freeing monomers
    - `monomer_tot`: Total concentration of enzyme monomers
    - `n_monomers`: Number of monomers in a single split enzyme
    - `η`: Ratio of specificities between kinase and phosphatase
    """
    return abs((monomer_tot - n_monomers * E) / (1 + η * (K / P)))
end


"""
    nodes_have_valid_split_enzyme_outputs(connectivity::AbstractMatrix)

Check if each node has only one type of output (kinase or phosphatase). This translates to checking that each column of the connectivity matrix has outputs of the same sign.

# Arguments
- `connectivity::AbstractMatrix`: A 2-dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network. 
    -1 indicates a kinase input, 1 indicates a phosphatase input, and 0 indicates no input. Please note that a node can only be
    a kinase or a phosphatase, not both.
"""
function nodes_have_valid_split_enzyme_outputs(connectivity::AbstractMatrix)
    result = true
    errmsg = ""
    # Check that each column has the same sign
    for i in axes(connectivity, 2)
        column = connectivity[:,i]
        if !(all(column .>= 0) || all(column .<= 0))
            result = false
            errmsg = "Invalid node output found. Each column of the connectivity matrix must have the same sign"
            return [result, errmsg]
        end
    end
    return [result, errmsg]
end


"""
    needs_constitutive_kinase(connectivity::AbstractMatrix)

Check if there is at least one node that doesn't have a kinase input (i.e. a node that doesn't have a -1 in its row).

# Arguments
- `connectivity::AbstractMatrix`: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network. 
    -1 indicates a kinase input, 1 indicates a phosphatase input, and 0 indicates no input. 
"""
function needs_constitutive_kinase(connectivity::AbstractMatrix)
    needs = false
    for i in axes(connectivity, 1)
        row = connectivity[i,:]
        if !any(row .== -1)
            needs = true
        end
    end
    return needs
end


"""
    split_enzyme_network(connectivity::AbstractMatrix)

Creates a ReactionSystem of interacting split enzymes based on the connectivity matrix.
The model is based on [Kimchi et al. 2020](https://www.science.org/doi/10.1126/sciadv.abc1939).
A constitutive phosphatase parameter (P̃) is added to all nodes. A constitutive kinase (K̃) is 
only added to those nodes that don't have a kinase input.

# Arguments
- `connectivity::AbstractMatrix`: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network. 
    -1 indicates a kinase input, 1 indicates a phosphatase input, and 0 indicates no input. Please note that a node can only be
    a kinase or a phosphatase, not both. This means that all the non-zero values of a column must have the same sign.
"""
function split_enzyme_network(connectivity::AbstractMatrix)
    # Check that the input is correct
    is_valid, errmsg = is_valid_connectivity(connectivity)
    has_correct_node_outputs, node_output_errmsg = nodes_have_valid_split_enzyme_outputs(connectivity)
    if !is_valid
        throw(DomainError(connectivity, errmsg))
    elseif !has_correct_node_outputs
        throw(DomainError(connectivity, node_output_errmsg))
    end
    # Number of nodes
    N = size(connectivity, 1)
    # Does the network need a constitutive kinase?
    constitutive_kinase_needed = needs_constitutive_kinase(connectivity)
    @variables t
    @species (X(t))[1:N]
    if constitutive_kinase_needed
        @parameters k_b[1:N], k_u[1:N], n[1:N], κ_tot[1:N], η, K̃, P̃
    else
        @parameters k_b[1:N], k_u[1:N], n[1:N], κ_tot[1:N], η, P̃
    end
    rxs = Reaction[]
    for i=1:N
        # Search for all inputs to node i
        K = 0
        P = P̃
        for j=1:N
            if connectivity[i,j] == -1
                K += X[j]
            elseif connectivity[i,j] == 1
                P += X[j]
            end
        end
        if iszero(K) && constitutive_kinase_needed
            K += K̃
        end
        # Define binding and unbinding reactions
        monomer = split_enzyme_monomer_concentration(X[i], K, P, κ_tot[i], n[i], η)
        push!(rxs, Reaction(k_b[i] * monomer ^ n[i], nothing, [X[i]]))
        push!(rxs, Reaction(k_u[i], [X[i]], nothing))
    end
    # Reorder parameters
    parameters = [[k_b[i], k_u[i], n[i], κ_tot[i]] for i=1:N]
    parameters = vcat(parameters...)
    parameters = vcat(parameters, η, P̃)
    if constitutive_kinase_needed
        parameters = vcat(parameters, K̃)
    end
    @named model = ReactionSystem(rxs, t, X, parameters)
    return model
end


"""
    sen_parameter_sets(model::ReactionSystem, samples::Int, random_seed::Int; parameter_limits=DEFAULT_SEN_PARAMETER_LIMITS, sampling_scales=DEFAULT_SEN_SAMPLING_SCALES, sampling_style="lhc")

Creates an array of paramter sets for a split enzyme network model.

# Arguments (Required)
- `model::ReactionSystem`: A split enzyme network model
- `samples::Int`: Number of parameter sets to generate
- `random_seed::Int`: Random seed for the sampling algorithm reproducibility

# Arguments (Optional)
- `parameter_limits::Dict`: Dictionary of parameter limits (lower and upper bound) for each type of parameter in the model.
- `sampling_scales::Dict`: Dictionary of sampling scales for each parameter in the model. Options are "log" or "linear".
- `sampling_style::String`: Sampling style. Options are "lhc" (Latin hypercube) or "random".

"""
function sen_parameter_sets(model::ReactionSystem, samples::Int, random_seed::Int; parameter_limits=DEFAULT_SEN_PARAMETER_LIMITS, sampling_scales=DEFAULT_SEN_SAMPLING_SCALES, sampling_style="lhc")
    N = length(species(model))
    P = length(parameters(model))

    limits = []
    for i in 1:N
        push!(limits, parameter_limits["k_b"])
        push!(limits, parameter_limits["k_u"])
        push!(limits, parameter_limits["n"])
        push!(limits, parameter_limits["κ_tot"])
    end

    push!(limits, parameter_limits["η"])

    for i in 1:(P-4*N-1)
        push!(limits, parameter_limits["P̃"])
    end

    scales = []
    for i in 1:N
        push!(scales, sampling_scales["k_b"])
        push!(scales, sampling_scales["k_u"])
        push!(scales, sampling_scales["n"])
        push!(scales, sampling_scales["κ_tot"])
    end

    push!(scales, sampling_scales["η"])

    for i in 1:(P-4*N-1)
        push!(scales, sampling_scales["P̃"])
    end

    parameter_array = generate_parameter_sets(samples, limits, scales, random_seed;sampling_style=sampling_style)

    # Force number of monomers to be integer
    parameter_array[:,3:4:4*N] = round.(parameter_array[:,3:4:4*N])

    return parameter_array
end


"""
    sen_timescale(k_u::AbstractVector, k_b::AbstractVector, n::AbstractVector, κ_tot::AbstractVector)

Calculate the slowest timescale of a split enzyme network model.

# Arguments
- `k_u::AbstractVector`: Vector of unbinding rate constants
- `k_b::AbstractVector`: Vector of binding rate constants
- `n::AbstractVector`: Vector of number of monomers in each split enzyme
- `κ_tot::AbstractVector`: Vector of total concentration of enzyme monomers

# Returns
- `timescale::Real`
"""
function sen_timescale(k_u::AbstractVector, k_b::AbstractVector, n::AbstractVector, κ_tot::AbstractVector)
    rates = vcat(k_u, k_b .* (κ_tot) .^ (n .- 1))
    return 1.0 / minimum(rates)
end


"""
    sen_equilibration_times(model::ReactionSystem, parameter_sets::AbstractArray; equilibration_time_multiplier=10)

Calculates the equilibration time for each dataset

# Arguments (Required)
- `model::ReactionSystem`: A split enzyme network model
- `parameter_sets::AbstractArray`: Array of parameter sets generated with [`sen_parameter_sets`](@ref)

# Arguments (Optional)
- `equilibration_time_multiplier::Real`: Factor by which to multiply the slowest timescale to get the equilibration time

# Returns
- `equilibration_times::Array{Float64}`: Vector of equilibration times for each parameter set
"""
function sen_equilibration_times(model::ReactionSystem, parameter_sets::AbstractArray; equilibration_time_multiplier=10)
    N = length(species(model))
    equilibration_times = []
    for i in axes(parameter_sets, 1)
        k_b = parameter_sets[i, 1:4:4*N]
        k_u = parameter_sets[i, 2:4:4*N]
        n = parameter_sets[i, 3:4:4*N]
        κ_tot = parameter_sets[i, 4:4:4*N]
        timescale = sen_timescale(k_u, k_b, n, κ_tot)
        push!(equilibration_times, timescale * equilibration_time_multiplier)
    end
    return equilibration_times
end


"""
    find_sen_oscillations(model::ReactionSystem, samples::Int; hyperparameters=DEFAULT_SEN_HYPERPARAMETERS)
"""
function find_sen_oscillations(connectivity::AbstractMatrix, samples::Int; hyperparameters=DEFAULT_SEN_HYPERPARAMETERS)
    # Unpack hyperparameters
    random_seed = hyperparameters["random_seed"]
    initial_conditions = hyperparameters["initial_conditions"]
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

    model = split_enzyme_network(connectivity)
    N = length(species(model))
    parameter_sets = sen_parameter_sets(model, samples, random_seed;
                                        parameter_limits=parameter_limits, 
                                        sampling_scales=sampling_scales, 
                                        sampling_style=sampling_style)
    equilibration_times = sen_equilibration_times(model, parameter_sets; 
                                                  equilibration_time_multiplier=equilibration_time_multiplier)

    if isnan(initial_conditions)
        initial_conditions = [zeros(N) for i in 1:samples]
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
    sen_hit_rate(connectivity::AbstractMAtrix, initial_samples::Int; target_oscillations::Int=100, max_samples::Int=1000000, max_trials::Int=5, hyperparameters=DEFAULT_SEN_HYPERPARAMETERS, verbose=true)

Estimates the hit rate of the oscillatory parameter set search in a split enzyme network model.

# Arguments (Required)
- `connectivity::AbstractMatrix`: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network. 
    -1 indicates a kinase input, 1 indicates a phosphatase input, and 0 indicates no input. Please note that a node can only be
    a kinase or a phosphatase, not both. This means that all the non-zero values of a column must have the same sign.
- `initial_samples::Int`: Number of initial parameter sets to generate

# Arguments (Optional)
- `target_oscillations::Int`: Number of oscillatory parameter sets to find
- `max_samples::Int`: Maximum number of samples allowed
- `max_trials::Int`: Maximum number of trials allowed
- `hyperparameters::Dict`: Dictionary of hyperparameters for the oscillatory parameter set search
- `verbose::Bool`: Whether to print the progress of the search
"""
function sen_hit_rate(connectivity::AbstractMatrix, initial_samples::Int; target_oscillations::Int=100, max_samples::Int=1000000, max_trials::Int=5, hyperparameters=DEFAULT_SEN_HYPERPARAMETERS, verbose=true)
    sen_result = find_sen_oscillations(connectivity, initial_samples; hyperparameters=hyperparameters)
    oscillatory = sum(sen_result["simulation_result"][!, "is_oscillatory"])
    if verbose
        println("Oscillatory: $oscillatory, Samples: $initial_samples")
    end
    hit_rate = oscillatory / initial_samples

    if oscillatory > target_oscillations
        return hit_rate
    else
        samples = initial_samples
        trial = 1
        while oscillatory < target_oscillations && trial < max_trials && samples < max_samples
            if oscillatory == 0
                samples *= 2
            else
                samples = ceil(Int, samples * (target_oscillations + 10) / oscillatory)
            end
            sen_result = find_sen_oscillations(connectivity, samples; hyperparameters=hyperparameters)
            oscillatory = sum(sen_result["simulation_result"][!, "is_oscillatory"])
            if verbose
                println("Oscillatory: $oscillatory, Samples: $samples")
            end
            hit_rate = oscillatory / samples
            trial += 1
        end
        return hit_rate
    end
end


"""
Goodwin-type negative feedback oscillator comprised of two split phosphatases, one split kinase,
one constitutive kinase, and one constitutive phosphatase. Model based on [Kimchi et al. 2020](https://www.science.org/doi/10.1126/sciadv.abc1939)
Connectivity: K phosphorylates P₁, P₁ dephosphorylates P₂, P₂ dephosphorylates K
"""
self_assembly_negative_feedback = @reaction_network Self_assembly_negative_feedback begin
    @species K(t)=0.5 P₁(t)=0.5 P₂(t)=0.5
    @parameters k_bκ=1e-1 k_uκ=1e0 n=2.0 κ_tot=3.0 k_bρ₁=1e0 k_uρ₁=1e-3 m=2.0 ρ₁_tot=7.0  k_bρ₂=1e-1 k_uρ₂=1e-1 l=2.0 ρ₂_tot=1.0 η=1.0 P̃=1e-5 K̃=1e-6
    # K reactions
    k_bκ * split_enzyme_monomer_concentration(K, K̃, P₂ + P̃, κ_tot, n, η) ^ n, ∅ --> K
    k_uκ, K --> ∅
    # P₁ reactions
    k_bρ₁ * split_enzyme_monomer_concentration(P₁, K, P̃, ρ₁_tot, m, η) ^ m, ∅ --> P₁
    k_uρ₁, P₁ --> ∅
    # P₂ reactions
    k_bρ₂ * split_enzyme_monomer_concentration(P₂, K̃, P₁ + P̃, ρ₂_tot, l, η) ^ l, ∅ --> P₂
    k_uρ₂, P₂ --> ∅
end
