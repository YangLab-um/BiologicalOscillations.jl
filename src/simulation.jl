"""
    equilibrate_ODEs(model::ReactionSystem, parameter_sets::AbstractVector, initial_conditions::AbstractVector, equilibration_times::AbstractVector, solver=RadauIIA5(), abstol=1e-7, reltol=1e-4, maxiters=1e7)

    Simulates a `model` for the different `parameter_sets` until the predetermined `equilibration_times` are reached.

# Arguments (Required)
- model::ReactionSystem: A set of differential equations encoded as a reaction ReactionSystem
- parameter_sets::AbstracVector: Vector with each element defining the parameter set for each simulation
- initial_conditions::AbstractVector: Vector with each element defining the initial condition for each parameter set
- equilibration_times::AbstractVector: Vector of simulation times

# Arguments (Optional)
- solver: Method for solving the differential equations. Default value `RadauIIA5()`
- abstol: Absolute tolerance of the solver. Default value 1e-7
- reltol: Relative tolerance of the solver. Default value 1e-4
- maxiters: Number of maximum iterations for the solver. Default value 1e7

# Returns
- equilibration_data::Dict: Dictionary containing the final state of each simulation and the frequency of the equilibration time series. The output is encoded as:
```julia
equilibration_data = Dict('final_state' => [[final_state_pset1_var1, ..., final_state_pset1_varN], ..., [final_state_psetM_var1, ..., final_state_psetM_varN]], 
                          'frequency' => [frequency_pset1, ..., frequency_psetM])
```
"""
function equilibrate_ODEs(model::ReactionSystem, parameter_sets::AbstractVector, initial_conditions::AbstractVector, equilibration_times::AbstractVector, solver=RadauIIA5(), abstol=1e-7, reltol=1e-4, maxiters=1e7)
    number_of_psets = size(parameter_sets, 1)
    # Define problem and output functions for EnsembleProblem
    function prob_function(prob::ODEProblem, i, ~)
        idx = Int(i)
        remake(prob, u0=initial_conditions[idx], tspan=(0.0, equilibration_times[idx]), p=parameter_sets[idx])
    end

    function output_func(sol::ODESolution, ~)
        frequency_data = calculate_main_frequency(sol, length(sol.t), length(sol.t))
        mean_frequency = mean(frequency_data["frequency"])
        endpoint = sol.u[end]
        return ([endpoint, mean_frequency], false)
    end

    ode_problem = ODEProblem(model, initial_conditions[1], (0.0, equilibration_times[1]), parameter_sets[1])
    ensemble_problem = EnsembleProblem(ode_problem, prob_func=prob_function, safetycopy=false, output_func=output_func)
    simulation = solve(ensemble_problem, solver, EnsembleThreads(), trajectories=number_of_psets,
                       abstol=abstol, reltol=reltol, maxiters=maxiters, dense=true)
    
    equilibration_data = Dict("final_state" => [simulation.u[i][1] for i=axes(simulation.u, 1)],
                            "frequency" => [simulation.u[i][2] for i=axes(simulation.u, 1)])

    return equilibration_data
end
