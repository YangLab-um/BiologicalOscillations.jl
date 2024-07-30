"""
    DEFAULT_SIMULATION_OUTPUT::Dict{String, Any}

Default selection of outputs for `find_pin_oscillations` and `find_grn_oscillations`.
Each key specifies what can be saved from the following options:
- "model": If value is set to true, saves the model that was simulated.
- "hyperparameters": If value is set to true, saves the hyperparameters that were used in the simulation.
- "parameter_sets": Saves the parameter sets that were simulated. A Dict is used to specify if the oscillatory and non-oscillatory parameter sets should be saved.
- "equilibration_result": Saves data regarding the equilibration of the parameter sets. What data is saved is specified by the bool value in a Dict. For all the options see below.
- "simulation_result": Saves data regarding the simulation of the parameter sets. What data is saved is specified by the bool value in a Dict. For all the options see below.

The "equilibration_result" output can have the following options:
- "parameter_index": Saves an identifier for each parameter set. Useful for indentifying parameter sets across result dataframes.
- "equilibration_times": Times used for equilibrating each parameter set.
- "final_velocity": Velocity of the final state after equilibration.
- "final_state": Final state after equilibration. Creates one column for each species.
- "frequency": Main frequency of the equilibration time-series. The value represents the average over all species.
- "is_steady_state": Final decision on whether the parameter set is in steady state or not.

The "simulation_result" output can have the following options:
- "parameter_index": Saves an identifier for each parameter set. Useful for indentifying parameter sets across result dataframes.
- "simulation_times": Times used for simulating each parameter set.
- "final_state": Final state after simulation. Creates one column for each species.
- "frequency": Main frequency of the simulation time-series. Creates one column for each species.
- "fft_power": Power of the main frequency of the simulation time-series. Creates one column for each species.
- "amplitude": Amplitude of the simulation time-series. Creates one column for each species.
- "peak_variation": Coefficient of variation of the peak values of the simulation time-series. Creates one column for each species.
- "trough_variation": Coefficient of variation of the trough values of the simulation time-series. Creates one column for each species.
- "is_oscillatory": Final decision on whether the parameter set is oscillatory or not.
"""
DEFAULT_SIMULATION_OUTPUT = Dict(
    "model" => true,
    "hyperparameters" => true,
    "parameter_sets" => Dict(
        "oscillatory" => true,
        "non_oscillatory" => true
    ),
    "equilibration_result" => Dict(
        "parameter_index" => true,
        "equilibration_time" => true,
        "final_velocity" => true,
        "final_state" => true,
        "frequency" => true,
        "is_steady_state" => true
    ),
    "simulation_result" => Dict(
        "parameter_index" => true,
        "simulation_time" => true,
        "final_state" => true,
        "frequency" => true,
        "fft_power" => true,
        "amplitude" => true,
        "peak_variation" => true,
        "trough_variation" => true,
        "is_oscillatory" => true
    )
)


"""
    DEFAULT_PIN_PARAMETER_LIMITS::Dict{String, Tuple{Float64, Float64}}

Default parameter limits for the generation of pin parameter sets. 
Keys are the names of the parameters, values are tuples of the form (lower_limit, upper_limit).
Expected keys are: α, β, γ, κ, η.
"""
DEFAULT_PIN_PARAMETER_LIMITS = Dict(
    "α" => (1e-2, 1e2), 
    "β" => (1e-2, 1e2), 
    "γ" => (1e2, 1e4), 
    "κ" => (1e-2, 1.0), 
    "η" => (1.0, 5.0)
)


"""
    DEFAULT_PIN_SAMPLING_SCALES::Dict{String, String}

Default sampling scales for the generation of pin parameter sets. 
Keys are the names of the parameters, values are strings indicating the sampling scale (log or linear).
Expected keys are: α, β, γ, κ, η.
"""
DEFAULT_PIN_SAMPLING_SCALES = Dict(
    "α" => "log", 
    "β" => "log", 
    "γ" => "log", 
    "κ" => "log", 
    "η" => "linear"
)


"""
    DEFAULT_PIN_HYPERPARAMETERS::Dict{String, Any}

Default hyperparameters for the various algorithms involved in the find_pin_oscillations pipeline.
The hyperparameters are the following:
- random_seed: Random seed used to generate parameter sets.
- initial_conditions: Set of initial conditions for `find_pin_oscillations`. Size should be equal to the number of samples indicated. If NaN, all species are initialized at 0.5.
- equilibration_time_multiplier: Multiplier applied to the slowest timescale to determine the equilibration time in `pin_equilibration_times`.
- dimensionless_time: Whether to use dimensionless time in `pin_parameter_sets`. If true, the α parameter of the first node is set to 1 for all parameter sets. If false, no changes are made to the parameter sets.
- parameter_limits: Default parameter limits for the generation of pin parameter sets.
- sampling_scales: Default sampling scales for the generation of pin parameter sets.
- sampling_style: Sampling style for `pin_parameter_sets`. Can be "lhc" (Latin Hypercube) or "random" (Uniform).
- solver: Solver used for the integration of the ODEs in `find_pin_oscillations`.
- abstol: Absolute tolerance of the ODE solver
- reltol: Relative tolerance of the ODE solver
- maxiters: Maximum number of iterations of the ODE solver
- fft_multiplier: Multiplier applied to the number of timepoints in an ODE solution to determine the main frequency. Used in `simulate_ODEs` when calling `calculate_main_frequency`.
- simulation_time_multiplier: Number of periods to simulate. Used in `calculate_simulation_times`
- freq_variation_threshold: Maximum tolerated variation in frequency between species in the solution to be declared oscillatory.
- power_threshold: Minimum spectral power that the main peak has to have to be declared oscillatory.
- amp_variation_threshold: Maximum tolerated value for peak/trough variation to be declared oscillatory.
- simulation_output: Default selection of outputs for `find_pin_oscillations`.
"""
DEFAULT_PIN_HYPERPARAMETERS = Dict(
    "random_seed" => 1234,
    "initial_conditions" => NaN,
    "equilibration_time_multiplier" => 10,
    "dimensionless_time" => false,
    "parameter_limits" => DEFAULT_PIN_PARAMETER_LIMITS,
    "sampling_scales" => DEFAULT_PIN_SAMPLING_SCALES,
    "sampling_style" => "lhc",
    "solver" => RadauIIA5(),
    "abstol" => 1e-7,
    "reltol" => 1e-4,
    "maxiters" => 1e7,
    "fft_multiplier" => 100,
    "simulation_time_multiplier" => 10,
    "freq_variation_threshold" => 0.05,
    "power_threshold" => 1e-7,
    "amp_variation_threshold" => 0.05,
    "simulation_output" => DEFAULT_SIMULATION_OUTPUT,
)


"""
    DEFAULT_GRN_PARAMETER_LIMITS::Dict{String, Tuple{Float64, Float64}}

Default parameter limits for the generation of GRN parameter sets. 
Keys are the names of the parameters, values are tuples of the form (lower_limit, upper_limit).
Expected keys are: α, β, δ, γ, κ, η.
"""
DEFAULT_GRN_PARAMETER_LIMITS = Dict(
    "α" => (1e-2, 1e2), 
    "β" => (1e-2, 1e2), 
    "δ" => (1e-2, 1e2), 
    "γ" => (1e-2, 1e2), 
    "κ" => (0.2, 1.0), 
    "η" => (1.0, 5.0)
)


"""
    DEFAULT_GRN_SAMPLING_SCALES::Dict{String, String}

Default sampling scales for the generation of GRN parameter sets. 
Keys are the names of the parameters, values are strings indicating the sampling scale (log or linear).
Expected keys are: α, β, δ, γ, κ, η
"""
DEFAULT_GRN_SAMPLING_SCALES = Dict(
    "α" => "log", 
    "β" => "log", 
    "δ" => "log", 
    "γ" => "log", 
    "κ" => "linear", 
    "η" => "linear"
)


"""
    DEFAULT_GRN_HYPERPARAMETERS::Dict{String, Any}

Default hyperparameters for the various algorithms involved in the find_grn_oscillations pipeline.
The hyperparameters are the following:
- random_seed: Random seed used to generate parameter sets.
- initial_conditions: Set of initial conditions for `find_grn_oscillations`. Size should be equal to the number of samples indicated. If NaN, all species are initialized at 10.
- equilibration_time_multiplier: Multiplier applied to the slowest timescale to determine the equilibration time in `grn_equilibration_times`.
- dimensionless_time: Whether to use dimensionless time in `grn_parameter_sets`. If true, the α parameter of the first node is set to 1 and the first N (N = number of nodes) κ are set to 1 for all parameter sets. If false, no changes are made to the parameter sets.
- parameter_limits: Default parameter limits for the generation of GRN parameter sets.
- sampling_scales: Default sampling scales for the generation of GRN parameter sets.
- sampling_style: Sampling style for `grn_parameter_sets`. Can be "lhc" (Latin Hypercube) or "random" (Uniform).
- solver: Solver used for the integration of the ODEs in `find_grn_oscillations`.
- abstol: Absolute tolerance of the ODE solver
- reltol: Relative tolerance of the ODE solver
- maxiters: Maximum number of iterations of the ODE solver
- fft_multiplier: Multiplier applied to the number of timepoints in an ODE solution to determine the main frequency. Used in `simulate_ODEs` when calling `calculate_main_frequency`.
- simulation_time_multiplier: Number of periods to simulate. Used in `calculate_simulation_times`
- freq_variation_threshold: Maximum tolerated variation in frequency between species in the solution to be declared oscillatory.
- power_threshold: Minimum spectral power that the main peak has to have to be declared oscillatory.
- amp_variation_threshold: Maximum tolerated value for peak/trough variation to be declared oscillatory.
- simulation_output: Default selection of outputs for `find_grn_oscillations`.
"""
DEFAULT_GRN_HYPERPARAMETERS = Dict(
    "random_seed" => 1234,
    "initial_conditions" => NaN,
    "equilibration_time_multiplier" => 10,
    "dimensionless_time" => false,
    "parameter_limits" => DEFAULT_GRN_PARAMETER_LIMITS,
    "sampling_scales" => DEFAULT_GRN_SAMPLING_SCALES,
    "sampling_style" => "lhc",
    "solver" => RadauIIA5(),
    "abstol" => 1e-7,
    "reltol" => 1e-4,
    "maxiters" => 1e7,
    "fft_multiplier" => 100,
    "simulation_time_multiplier" => 10,
    "freq_variation_threshold" => 0.05,
    "power_threshold" => 1e-7,
    "amp_variation_threshold" => 0.05,
    "simulation_output" => DEFAULT_SIMULATION_OUTPUT,
)