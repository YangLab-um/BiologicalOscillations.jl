"""
    generate_parameter_sets(samples::Int, parameter_limits::AbstractVector, sampling_scales::AbstractVector, sampling_style::String="lhc")

Creates an array of parameter sets within the limits provided in `parameter_limits` using either Latin Hypercube Sampling or Random Sampling. Each parameter can be sampled linearly or logarithmically.

# Arguments (Required)
- `samples::Int``: Number of parameter sets to be generated
- `parameter_limits::AbstractVector`: Array of tuples defining the lower and upper limits of each individual parameter in a model. The length of this array is used to calculate the number of parameters in each parameter set.
- `sampling_scales::AbstractVector`: Array of strings containing the sampling scale for each individual parameter in a model. Accepted strings are "linear" and "log".
- `random_seed::Int`: Seed for the random number generator.

# Arguments (Optional)
- `sampling_style::String`: Sampling style of the algorithm. Accepted strings are "lhc" and "random".

# Returns
- `parameter_sets::AbstractArray`: Array of parameter sets of length `samples`.
"""
function generate_parameter_sets(samples::Int, parameter_limits::AbstractVector, sampling_scales::AbstractVector, random_seed::Int; sampling_style::String="lhc")
    Random.seed!(random_seed)
    number_of_parameters = length(parameter_limits)
    # Draw sample
    if lowercase(sampling_style) == "lhc"
        LHCplan = randomLHC(samples, number_of_parameters)
    elseif lowercase(sampling_style) == "random"
        LHCplan = rand([i for i in 1:samples], samples, number_of_parameters)
    end

    # Rescale parameters
    scaling_plan = Tuple{Float64,Float64}[]
    for (idx, limits) in enumerate(parameter_limits)
        if sampling_scales[idx] == "linear"
            scaling_plan = vcat(scaling_plan, [limits])
        elseif sampling_scales[idx] == "log"
            scaling_plan = vcat(scaling_plan, [(log10(limits[1]), log10(limits[2]))])
        end
    end

    parameter_sets = scaleLHC(LHCplan, scaling_plan)

    for (idx, scale) in enumerate(sampling_scales)
        if scale == "log"
            parameter_sets[:,idx] = 10 .^ parameter_sets[:,idx]
        end
    end
    
    return parameter_sets
end


"""
    equilibrate_ODEs(model::ReactionSystem, parameter_sets::AbstractArray, initial_conditions::AbstractVector, equilibration_times::AbstractVector; solver=RadauIIA5(), abstol::Real=1e-7, reltol::Real=1e-4, maxiters=1e7)

Simulates a `model` for the different `parameter_sets` until the predetermined `equilibration_times` are reached.

# Arguments (Required)
- `model::ReactionSystem`: A set of differential equations encoded as a reaction ReactionSystem
- `parameter_sets::AbstractArray`: Array where each row defines the parameter set for each simulation
- `initial_conditions::AbstractVector`: Vector with each element defining the initial condition for each parameter set
- `equilibration_times::AbstractVector`: Vector of simulation times

# Arguments (Optional)
- `solver`: Method for solving the differential equations. Default value `RadauIIA5()`
- `abstol::Real`: Absolute tolerance of the solver. Default value 1e-7
- `reltol::Real`: Relative tolerance of the solver. Default value 1e-4
- `maxiters::Int`: Number of maximum iterations for the solver. Default value 1e7

# Returns
- `equilibration_data::Dict`: Dictionary containing the final state and velocity of each simulation as well as the frequency of the time series. The velocit is calculated as the norm of the derivative for each equation evaluated at the final equilibration state. The output is encoded as:
```julia
equilibration_data = Dict('final_state' => [[final_state_pset1_var1, ..., final_state_pset1_varN], ..., [final_state_psetM_var1, ..., final_state_psetM_varN]], 
                          'final_velocity' => [velocity_pset1, ..., velocity_psetM],
                          'frequency' => [frequency_pset1, ..., frequency_psetM])
```
"""
function equilibrate_ODEs(model::ReactionSystem, parameter_sets::AbstractArray, initial_conditions::AbstractVector, equilibration_times::AbstractVector; solver=RadauIIA5(), abstol::Real=1e-7, reltol::Real=1e-4, maxiters=1e7)
    number_of_psets = size(parameter_sets, 1)
    # Define problem and output functions for EnsembleProblem
    function prob_function(prob::ODEProblem, i, ~)
        idx = Int(i)
        remake(prob, u0=initial_conditions[idx], tspan=(0.0, equilibration_times[idx]), p=parameter_sets[idx,:])
    end

    function output_func(sol::ODESolution, i)
        idx = Int(i)
        ode_problem = ODEProblem(model, initial_conditions[idx], (0.0, equilibration_times[idx]), parameter_sets[idx,:])
        frequency_data = calculate_main_frequency(sol, length(sol.t), length(sol.t))
        mean_frequency = mean(frequency_data["frequency"])
        final_state = sol.u[end]
        final_derivative = ode_problem.f(final_state, parameter_sets[idx,:], equilibration_times[idx])
        final_velocity = sqrt(sum(final_derivative.^2))
        return ([final_state, final_velocity, mean_frequency], false)
    end

    ode_problem = ODEProblem(model, initial_conditions[1], (0.0, equilibration_times[1]), parameter_sets[1,:])
    ensemble_problem = EnsembleProblem(ode_problem, prob_func=prob_function, safetycopy=false, output_func=output_func)
    simulation = solve(ensemble_problem, solver, EnsembleThreads(), trajectories=number_of_psets,
                       abstol=abstol, reltol=reltol, maxiters=maxiters, dense=true)
    
    equilibration_data = Dict("final_state" => [simulation.u[i][1] for i=axes(simulation.u, 1)],
                              "final_velocity" => [simulation.u[i][2] for i=axes(simulation.u, 1)],
                              "frequency" => [simulation.u[i][3] for i=axes(simulation.u, 1)])

    return equilibration_data
end


"""
    simulate_ODEs(model::ReactionSystem, parameter_sets::AbstractVector, initial_conditions::AbstractVector, simulation_times::AbstractVector, solver=RadauIIA5(), abstol=1e-7, reltol=1e-4, maxiters=1e7)

Simulates a `model` for the different `parameter_sets` until the predetermined `simulation_times` and calculates features relevant for the detection of oscillations such as frequency and amplitude data.

# Arguments (Required)
- `model::ReactionSystem`: A set of differential equations encoded as a reaction ReactionSystem
- `parameter_sets::AbstractArray`: Array where each row defines the parameter set for each simulation
- `initial_conditions::AbstractVector`: Vector with each element defining the initial condition for each parameter set
- `simulation_times::AbstractVector`: Vector of simulation times

# Arguments (Optional)
- `solver`: Method for solving the differential equations. Default value `RadauIIA5()`
- `abstol::Real`: Absolute tolerance of the solver. Default value 1e-7
- `reltol::Real`: Relative tolerance of the solver. Default value 1e-4
- `maxiters::Int`: Number of maximum iterations for the solver. Default value 1e7
- `fft_multiplier::Int`: Factor by which the number of time points of the solution is multiplied to be used as the number of points in the periodogram. Default value 100

# Returns
- simulation_data::Dict: Dictionary containing the ODESolution, frequency data, and amplitude data for each parameter set. Frequency data is the output of [`calculate_main_frequency`](@ref) and amplitude data the output of [`calculate_amplitude`](@ref). The final output is encoded as:
```julia
simulation_data = Dict('final_state' => [final_state_pset1, ..., final_state_psetM], 
                       'frequency_data' => [freq_data_pset1, ..., freq_data_psetM],
                       'amplitude_data' => [amp_data_pset1, ..., amp_data_psetM])
```

# Note
- It is important that the parameter sets and initial conditions given to this function to be pre-equilibrated. Use [`equilibrate_ODEs`](@ref) for this purpose.
"""
function simulate_ODEs(model::ReactionSystem, parameter_sets::AbstractArray, initial_conditions::AbstractVector, simulation_times::AbstractVector; solver=RadauIIA5(), abstol=1e-7, reltol=1e-4, maxiters=1e7, fft_multiplier=100)
    number_of_psets = size(parameter_sets, 1)
    # Define problem and output functions for EnsembleProblem
    function prob_function(prob::ODEProblem, i, ~)
        idx = Int(i)
        remake(prob, u0=initial_conditions[idx], tspan=(0.0, simulation_times[idx]), p=parameter_sets[idx,:])
    end

    function output_func(sol::ODESolution, i)
        frequency_data = calculate_main_frequency(sol, length(sol.t), fft_multiplier * length(sol.t))
        amplitude_data = calculate_amplitude(sol)
        return ([sol.u[end], frequency_data, amplitude_data], false)
    end

    ode_problem = ODEProblem(model, initial_conditions[1], (0.0, simulation_times[1]), parameter_sets[1,:])
    ensemble_problem = EnsembleProblem(ode_problem, prob_func=prob_function, safetycopy=false, output_func=output_func)
    simulation = solve(ensemble_problem, solver, EnsembleThreads(), trajectories=number_of_psets,
                       abstol=abstol, reltol=reltol, maxiters=maxiters, dense=true)
    
    simulation_data = Dict("final_state" => [simulation.u[i][1] for i=axes(simulation.u, 1)],
                           "frequency_data" => [simulation.u[i][2] for i=axes(simulation.u, 1)],
                           "amplitude_data" => [simulation.u[i][3] for i=axes(simulation.u, 1)])

    return simulation_data
end


"""
    calculate_simulation_times(equilibration_data::Dict, equilibration_times::AbstractVector; simulation_time_multiplier=10)

Calculates the simulation times for each parameter set in a model

# Arguments (Required)
- `equilibration_data::Dict`: Dictionary of equilibration data generated with [`equilibrate_ODEs`](@ref)
- `equilibration_times::AbstractVector`: Array of equilibration times generated with [`pin_equilibration_times`](@ref)

# Arguments (Optional)
- `simulation_time_multiplier::Real`: Factor by which the period (given by the estimated frequency from the equilibration) is multiplied by in order to calculate the final simulation time

# Returns
- `simulation_times::Array{Float64}`: Array of simulation times
"""
function calculate_simulation_times(equilibration_data::Dict, equilibration_times::AbstractVector; simulation_time_multiplier=10)
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
    calculate_oscillatory_status(simulation_data::Dict; freq_variation_threshold::Real=0.05, power_threshold::Real=1e-7, amp_variation_threshold::Real=0.05)

Calculates the oscillatory status for each parameter set using the functionality of [`is_ODE_oscillatory`](@ref)

# Arguments (Required)
- `simulation_data::Dict`: Dictionary of simulation data generated with [`simulate_ODEs`](@ref)
- `freq_variation_threshold::Real`: Maximum tolerated variation in frequency between species in the solution to be declared oscillatory.
- `power_threshold::Real`: Minimum spectral power that the main peak has to have to be declared oscillatory.
- `amp_variation_threshold::Real`: Maximum tolerated value for peak/trough variation to be declared oscillatory.

# Returns
- `oscillatory_status::Array{Bool}`: Boolean array indicating whether each parameter set is oscillatory or not.
"""
function calculate_oscillatory_status(simulation_data::Dict; freq_variation_threshold::Real=0.05, power_threshold::Real=1e-7, amp_variation_threshold::Real=0.05)
    oscillatory_status = Array{Bool}(undef, 0)

    for i in axes(simulation_data["frequency_data"], 1)
        freq = simulation_data["frequency_data"][i]
        amp = simulation_data["amplitude_data"][i]
        decision = is_ODE_oscillatory(freq, amp; 
                                      freq_variation_threshold=freq_variation_threshold, 
                                      power_threshold=power_threshold, 
                                      amp_variation_threshold=amp_variation_threshold)
        push!(oscillatory_status, decision)
    end

    return oscillatory_status
end


"""
    generate_find_oscillations_output(model::ReactionSystem, parameter_sets::AbstractArray, equilibration_data::Dict, equilibration_times::AbstractVector, simulation_data::Dict, simulation_times::AbstractVector, oscillatory_status::AbstractVector, hyperparameters::Dict)

Generates a dataframe with the results of the simulation and equilibration of a model. The final output is controlled by the `hyperparameters` dict.

# Arguments (Required)
- `model::ReactionSystem`: A set of differential equations encoded as a reaction ReactionSystem
- `parameter_sets::AbstractArray`: Array where each row defines the parameter set for each simulation
- `equilibration_data::Dict`: Dictionary of equilibration data generated with [`equilibrate_ODEs`](@ref)
- `equilibration_times::AbstractVector`: Array of equilibration times 
- `simulation_data::Dict`: Dictionary of simulation data generated with [`simulate_ODEs`](@ref)
- `simulation_times::AbstractVector`: Array of simulation times
- `oscillatory_status::AbstractVector`: Array of oscillatory status for each parameter set
- `hyperparameters::Dict`: Dictionary of hyperparameters for the simulation

# Arguments (Optional)
- `filter_results::Bool`: Whether to filter the results based on equilibration's final velocity. Default value `true`
"""
function generate_find_oscillations_output(model::ReactionSystem, parameter_sets::AbstractArray, 
                                           equilibration_data::Dict, equilibration_times::AbstractVector, 
                                           simulation_data::Dict, simulation_times::AbstractVector, 
                                           oscillatory_status::AbstractVector, hyperparameters::Dict;
                                           filter_results=true)
    result = Dict()
    samples = size(parameter_sets, 1)
    N = length(species(model))
    velocity = equilibration_data["final_velocity"]
    cutoff = 10 ^ mean(log10.(velocity))
    if filter_results
        filter = velocity .> cutoff
    else
        filter = ones(Bool, samples)
    end
    # Model
    if haskey(hyperparameters["simulation_output"], "model")
        if hyperparameters["simulation_output"]["model"] == true
            result["model"] = model
        end
    end
    # Hyperparameters
    if haskey(hyperparameters["simulation_output"], "hyperparameters")
        if hyperparameters["simulation_output"]["hyperparameters"] == true
            result["hyperparameters"] = hyperparameters
        end
    end
    # Parameter sets
    if haskey(hyperparameters["simulation_output"], "parameter_sets")
        parameter_map = paramsmap(model)
        parameter_names = Array{String}(undef, length(parameter_map))
        for (k, v) in parameter_map
            parameter_names[v] = string(k)
        end
        result["parameter_sets"] = Dict()
        if haskey(hyperparameters["simulation_output"]["parameter_sets"], "oscillatory") 
            if hyperparameters["simulation_output"]["parameter_sets"]["oscillatory"] == true
                oscillatory_params_summary = Dict()
                oscillatory_idxs = collect(1:1:samples)[filter][oscillatory_status]
                oscillatory_params_summary["parameter_index"] = oscillatory_idxs
                oscillatory_parameters = parameter_sets[oscillatory_idxs, :]
                for i=1:length(parameter_map)
                    oscillatory_params_summary["$(parameter_names[i])"] = oscillatory_parameters[:,i]
                end
                oscillatory_dataframe = DataFrame(oscillatory_params_summary)
                column_order = vcat(["parameter_index"], parameter_names)
                result["parameter_sets"]["oscillatory"] = select(oscillatory_dataframe, column_order)
            end
        end
        if haskey(hyperparameters["simulation_output"]["parameter_sets"], "non_oscillatory")
            if hyperparameters["simulation_output"]["parameter_sets"]["non_oscillatory"] == true
                fixed_point_params_summary = Dict()
                steady_state_idxs = collect(1:1:samples)[.!filter]
                non_oscillatory_idxs = collect(1:1:samples)[filter][.!oscillatory_status]
                fixed_point_idxs = vcat(steady_state_idxs, non_oscillatory_idxs)
                fixed_point_params_summary["parameter_index"] = fixed_point_idxs
                non_oscillatory_parameters = parameter_sets[fixed_point_idxs, :]
                for i=1:length(parameter_map)
                    fixed_point_params_summary["$(parameter_names[i])"] = non_oscillatory_parameters[:,i]
                end
                fixed_point_dataframe = DataFrame(fixed_point_params_summary)
                column_order = vcat(["parameter_index"], parameter_names)
                result["parameter_sets"]["non_oscillatory"] = select(fixed_point_dataframe, column_order)
            end
        end
        if haskey(hyperparameters["simulation_output"]["parameter_sets"], "all")
            if hyperparameters["simulation_output"]["parameter_sets"]["all"] == true
                all_params_summary = Dict()
                all_params_summary["parameter_index"] = collect(1:1:samples)
                for i=1:length(parameter_map)
                    all_params_summary["$(parameter_names[i])"] = parameter_sets[:,i]
                end
                all_dataframe = DataFrame(all_params_summary)
                column_order = vcat(["parameter_index"], parameter_names)
                result["parameter_sets"] = select(all_dataframe, column_order)
            end
        end
    end
    # Equilibration result
    if haskey(hyperparameters["simulation_output"], "equilibration_result")
        equilibration_summary = Dict()
        if haskey(hyperparameters["simulation_output"]["equilibration_result"], "parameter_index")
            if hyperparameters["simulation_output"]["equilibration_result"]["parameter_index"] == true
                equilibration_summary["parameter_index"] = collect(1:1:samples)
            end
        end
        if haskey(hyperparameters["simulation_output"]["equilibration_result"], "equilibration_time")
            if hyperparameters["simulation_output"]["equilibration_result"]["equilibration_time"] == true
                equilibration_summary["equilibration_time"] = equilibration_times
            end
        end
        if haskey(hyperparameters["simulation_output"]["equilibration_result"], "final_velocity")
            if hyperparameters["simulation_output"]["equilibration_result"]["final_velocity"] == true
                equilibration_summary["final_velocity"] = equilibration_data["final_velocity"]
            end
        end
        if haskey(hyperparameters["simulation_output"]["equilibration_result"], "final_state")
            if hyperparameters["simulation_output"]["equilibration_result"]["final_state"] == true
                final_state = mapreduce(permutedims, vcat, equilibration_data["final_state"])
                for i=1:N
                    equilibration_summary["final_state_$(i)"] = final_state[:,i]
                end
            end
        end
        if haskey(hyperparameters["simulation_output"]["equilibration_result"], "frequency")
            if hyperparameters["simulation_output"]["equilibration_result"]["frequency"] == true
                equilibration_summary["frequency"] = equilibration_data["frequency"]
            end
        end
        if haskey(hyperparameters["simulation_output"]["equilibration_result"], "is_steady_state")
            if hyperparameters["simulation_output"]["equilibration_result"]["is_steady_state"] == true
                equilibration_summary["is_steady_state"] = .!filter
            end
        end
        result["equilibration_result"] = DataFrame(equilibration_summary)
    end

    # Simulation result
    if haskey(hyperparameters["simulation_output"], "simulation_result")
        simulation_summary = Dict()
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "parameter_index")
            if hyperparameters["simulation_output"]["simulation_result"]["parameter_index"] == true
                simulation_summary["parameter_index"] = collect(1:1:samples)[filter]
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "simulation_time")
            if hyperparameters["simulation_output"]["simulation_result"]["simulation_time"] == true
                simulation_summary["simulation_time"] = simulation_times[filter]
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "final_state")
            if hyperparameters["simulation_output"]["simulation_result"]["final_state"] == true
                final_state = mapreduce(permutedims, vcat, simulation_data["final_state"])
                for i=1:N
                    simulation_summary["final_state_$(i)"] = final_state[:,i]
                end
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "is_oscillatory")
            if hyperparameters["simulation_output"]["simulation_result"]["is_oscillatory"] == true
                simulation_summary["is_oscillatory"] = oscillatory_status
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "frequency")
            if hyperparameters["simulation_output"]["simulation_result"]["frequency"] == true
                for i=1:N
                    frequency = Array{Float64}(undef, 0)
                    for j=1:sum(filter)
                        push!(frequency, simulation_data["frequency_data"][j]["frequency"][i])
                    end
                    simulation_summary["frequency_$(i)"] = frequency
                end
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "fft_power")
            if hyperparameters["simulation_output"]["simulation_result"]["fft_power"] == true
                for i=1:N
                    power = Array{Float64}(undef, 0)
                    for j=1:sum(filter)
                        push!(power, simulation_data["frequency_data"][j]["power"][i])
                    end
                    simulation_summary["fft_power_$(i)"] = power
                end
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "amplitude")
            if hyperparameters["simulation_output"]["simulation_result"]["amplitude"] == true
                for i=1:N
                    amplitude = Array{Float64}(undef, 0)
                    for j=1:sum(filter)
                        push!(amplitude, simulation_data["amplitude_data"][j]["amplitude"][i])
                    end
                    simulation_summary["amplitude_$(i)"] = amplitude
                end
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "peak_variation")
            if hyperparameters["simulation_output"]["simulation_result"]["peak_variation"] == true
                for i=1:N
                    peak_variation = Array{Float64}(undef, 0)
                    for j=1:sum(filter)
                        push!(peak_variation, simulation_data["amplitude_data"][j]["peak_variation"][i])
                    end
                    simulation_summary["peak_variation_$(i)"] = peak_variation
                end
            end
        end
        if haskey(hyperparameters["simulation_output"]["simulation_result"], "trough_variation")
            if hyperparameters["simulation_output"]["simulation_result"]["trough_variation"] == true
                for i=1:N
                    trough_variation = Array{Float64}(undef, 0)
                    for j=1:sum(filter)
                        push!(trough_variation, simulation_data["amplitude_data"][j]["trough_variation"][i])
                    end
                    simulation_summary["trough_variation_$(i)"] = trough_variation
                end
            end
        end
        result["simulation_result"] = DataFrame(simulation_summary)
    end
    return result
end


"""
    create_random_parameter_set_perturbation(parameter_sets::AbstractArray, perturbation_percentage::Real, random_seed::Int)

Creates a perturbed parameter by randomly choosing a single parameter from each parameter set and increasing it by a percentage of its value.

# Arguments (Required)
- `parameter_sets::AbstractArray`: Array where each row defines the parameter set for each simulation
- `perturbation_percentage::Real`: Percentage by which the parameter is perturbed
- `random_seed::Int`: Seed for the random number generator

# Arguments (Optional)
- `keep_constant::AbstractArray{Int}`: Array of indices of parameters that should not be perturbed

# Returns
- `perturbed_parameter_sets::AbstractArray`: Array of perturbed parameter sets
"""
function create_random_parameter_set_perturbation(parameter_sets::AbstractArray, perturbation_percentage::Real, random_seed::Int; keep_constant::AbstractArray{Int}=Int[])
    Random.seed!(random_seed)
    number_of_parameters = size(parameter_sets, 2)
    perturbed_parameter_sets = zeros(size(parameter_sets))
    for i in axes(parameter_sets, 1)
        perturbed_set = copy(parameter_sets[i, :])
        idx = rand(1:number_of_parameters)
        trials = 0
        while idx in keep_constant
            if trials > 100
                error("Could not find a parameter to perturb. Check the keep_constant array.")
            end
            idx = rand(1:number_of_parameters)
            trials += 1
        end
        perturbed_set[idx] = perturbed_set[idx] * (1 + perturbation_percentage)
        perturbed_parameter_sets[i, :] = perturbed_set
    end
    return perturbed_parameter_sets
end


"""
    feature_change_from_perturbation(find_oscillations_result::Dict, perturbation_result::Dict)

Calculate frequency and amplitude changes between the original and perturbed parameter sets

# Arguments (Required)
- `find_oscillations_result::Dict`: Results of the original parameter sets generated with [`find_pin_oscillations`](@ref)
- `perturbation_result::Dict`: Results of the perturbed parameter sets generated with [`simulate_pin_parameter_perturbations`](@ref)

# Returns
- `perturbation_analysis::DataFrame`: DataFrame containing the frequency and amplitude changes between the original and perturbed parameter sets
"""
function feature_change_from_perturbation(find_oscillations_result::Dict, perturbation_result::Dict)
    # Average frequency and amplitude across nodes
    ## Original result
    oscillatory_df = filter(row -> row["is_oscillatory"] == true, find_oscillations_result["simulation_result"])
    original_frequency_df = select(oscillatory_df, r"frequency_*")
    original_amplitude_df = select(oscillatory_df, r"amplitude_*")
    original_frequencies = mean.(eachrow(original_frequency_df))
    original_amplitudes = mean.(eachrow(original_amplitude_df))
    ## Perturbed result
    perturbed_frequency_df = select(perturbation_result["simulation_result"], r"frequency_*")
    perturbed_amplitude_df = select(perturbation_result["simulation_result"], r"amplitude_*")
    perturbed_frequencies = mean.(eachrow(perturbed_frequency_df))
    perturbed_amplitudes = mean.(eachrow(perturbed_amplitude_df))
    # Calculate changes
    freq_change = (perturbed_frequencies .- original_frequencies) ./ original_frequencies
    amp_change = (perturbed_amplitudes .- original_amplitudes) ./ original_amplitudes
    # Create dataframe
    parameter_index = oscillatory_df.parameter_index
    feature_change = DataFrame("parameter_index" => parameter_index,
                               "frequency_change" => freq_change,
                               "amplitude_change" => amp_change)
    return feature_change
end


"""
    calculate_perturbed_parameter_index(original_parameter_sets::AbstractArray, perturbed_parameter_sets::AbstractArray)

Calculates the perturbed parameter index for all parameter sets

# Arguments (Required)
- `original_parameter_set::AbstractArray`: Array of original parameter sets
- `perturbed_parameter_set::AbstractArray`: Array of perturbed parameter sets
"""
function calculate_perturbed_parameter_index(original_parameter_set::AbstractArray, perturbed_parameter_set::AbstractArray)
    perturbed_parameter_index = []
    for i in axes(original_parameter_set, 1)
        original_set = original_parameter_set[i, :]
        perturbed_set = perturbed_parameter_set[i, :]
        diff = abs.(original_set .- perturbed_set)
        idx = findfirst(diff .!= 0)
        push!(perturbed_parameter_index, idx)
    end
    return perturbed_parameter_index
end