"""
    DEFAULT_PIN_PARAMETER_LIMITS::Dict{String, Tuple{Float64, Float64}}

Default parameter limits for the generation of pin parameter sets. Keys are the names of the parameters, values are tuples of the form (lower_limit, upper_limit).
"""
DEFAULT_PIN_PARAMETER_LIMITS = Dict(
    "α" => (1e-2, 1e2), 
    "β" => (1e-2, 1e2), 
    "γ" => (1e2, 1e4), 
    "κ" => (0.2, 1.0), 
    "η" => (1.0, 5.0)
)


"""
    DEFAULT_PIN_SAMPLING_SCALES::Dict{String, String}

Default sampling scales for the generation of pin parameter sets. Keys are the names of the parameters, values are strings indicating the sampling scale (log or linear).
"""
DEFAULT_PIN_SAMPLING_SCALES = Dict(
    "α" => "log", 
    "β" => "log", 
    "γ" => "log", 
    "κ" => "linear", 
    "η" => "linear"
)


"""
    DEFAULT_PIN_HYPERPARAMETERS::Dict{String, Any}

Default hyperparameters for the various algorithms involved in the find_pin_oscillations pipeline.
"""
DEFAULT_PIN_HYPERPARAMETERS = Dict(
    "initial_conditions" => NaN, # Set of initial conditions for `find_pin_oscillations`. Size should be equal to the number of samples indicated. If NaN, all species are initialized at 0.5.
    "equilibration_time_multiplier" => 10, # Multiplier applied to the slowest timescale to determine the equilibration time in `pin_equilibration_times`.
    "dimensionless_time" => true, # Whether to use dimensionless time in `pin_parameter_sets`. If true, the α parameter of the first node is set to 1 for all parameter sets. If false, no changes are made to the parameter sets.
    "parameter_limits" => DEFAULT_PIN_PARAMETER_LIMITS,
    "sampling_scales" => DEFAULT_PIN_SAMPLING_SCALES,
    "sampling_style" => "lhc", # Sampling style for `pin_parameter_sets`. Can be "lhc" (Latin Hypercube) or "random" (Uniform).
    "solver" => RadauIIA5(), # Solver used for the integration of the ODEs in `find_pin_oscillations`.
    "abstol" => 1e-7, # Absolute tolerance of the ODE solver
    "reltol" => 1e-4, # Relative tolerance of the ODE solver
    "maxiters" => 1e7, # Maximum number of iterations of the ODE solver
    "fft_multiplier" => 100, # Multiplier applied to the number of timepoints in an ODE solution to determine the main frequency. Used in `simulate_ODEs` when calling `calculate_main_frequency`.
    "simulation_time_multiplier" => 10, # Number of periods to simulate. Used in `calculate_simulation_times`
    "freq_variation_threshold" => 0.05, # Maximum tolerated variation in frequency between species in the solution to be declared oscillatory.
    "power_threshold" => 1e-7, # Minimum spectral power that the main peak has to have to be declared oscillatory.
    "amp_variation_threshold" => 0.05 # Maximum tolerated value for peak/trough variation to be declared oscillatory.
)


"""
    DEFAULT_GRN_PARAMETER_LIMITS::Dict{String, Tuple{Float64, Float64}}

Default parameter limits for the generation of GRN parameter sets. Keys are the names of the parameters, values are tuples of the form (lower_limit, upper_limit).
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

Default sampling scales for the generation of GRN parameter sets. Keys are the names of the parameters, values are strings indicating the sampling scale (log or linear).
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
"""
DEFAULT_GRN_HYPERPARAMETERS = Dict(
    "initial_conditions" => NaN,  # Set of initial conditions for `find_grn_oscillations`. Size should be equal to the number of samples indicated. If NaN, all species are initialized at 10.
    "equilibration_time_multiplier" => 10, # Multiplier applied to the slowest timescale to determine the equilibration time in `grn_equilibration_times`.
    "dimensionless_time" => true, # Whether to use dimensionless time in `grn_parameter_sets`. If true, the α parameter of the first node is set to 1 and the first N (N = number of nodes) κ are set to 1 for all parameter sets. If false, no changes are made to the parameter sets.
    "parameter_limits" => DEFAULT_GRN_PARAMETER_LIMITS,
    "sampling_scales" => DEFAULT_GRN_SAMPLING_SCALES,
    "sampling_style" => "lhc",  # Sampling style for `grn_parameter_sets`. Can be "lhc" (Latin Hypercube) or "random" (Uniform).
    "solver" => RadauIIA5(), # Solver used for the integration of the ODEs in `find_grn_oscillations`.
    "abstol" => 1e-7, # Absolute tolerance of the ODE solver
    "reltol" => 1e-4, # Relative tolerance of the ODE solver
    "maxiters" => 1e7, # Maximum number of iterations of the ODE solver
    "fft_multiplier" => 100, # Multiplier applied to the number of timepoints in an ODE solution to determine the main frequency. Used in `simulate_ODEs` when calling `calculate_main_frequency`.
    "simulation_time_multiplier" => 10, # Number of periods to simulate. Used in `calculate_simulation_times`
    "freq_variation_threshold" => 0.05, # Maximum tolerated variation in frequency between species in the solution to be declared oscillatory.
    "power_threshold" => 1e-7, # Minimum spectral power that the main peak has to have to be declared oscillatory.
    "amp_variation_threshold" => 0.05 # Maximum tolerated value for peak/trough variation to be declared oscillatory.
)