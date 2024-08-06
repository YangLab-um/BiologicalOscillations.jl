using BiologicalOscillations, Catalyst, ModelingToolkit, DifferentialEquations, Statistics

# Test that correct features are obtained for oscillatory solutions
# Brusselator
brusselator = @reaction_network begin
    k1, ∅ --> X
    k2, 2X + Y --> 3X
    k3, X --> Y
    k4, X --> ∅
end
p = [1.0, 1.0, 3.0, 1.0]
u0 = [1.0, 1.0]
tspan=(.0, 50.0)
oprob = ODEProblem(brusselator, u0, tspan, p)
brusselator_sol = solve(oprob, Tsit5())

signal_sampling = Int(length(brusselator_sol.t))
spectrum_sampling = 100 * signal_sampling
amplitude_data = calculate_amplitude(brusselator_sol)
frequency_data = calculate_main_frequency(brusselator_sol, signal_sampling, spectrum_sampling)

@test abs(frequency_data["frequency"][1] - 0.078)/0.078 < 0.05
@test abs(frequency_data["frequency"][2] - 0.078)/0.078 < 0.05

@test abs(frequency_data["power"][1] - 23.19)/23.19 < 0.05
@test abs(frequency_data["power"][2] - 189.5)/189.5 < 0.05

@test abs(amplitude_data["amplitude"][1] - 7.33)/7.33 < 0.05
@test abs(amplitude_data["amplitude"][2] - 8.18)/8.18 < 0.05

@test abs(amplitude_data["peak_variation"][1] - 0.00047)/0.00047 < 0.05
@test abs(amplitude_data["peak_variation"][2] - 0.0014)/0.0014 < 0.05

@test abs(amplitude_data["trough_variation"][1] - 4.5e-5)/4.5e-5 < 0.05
@test abs(amplitude_data["trough_variation"][2] - 5.6e-5)/5.6e-5 < 0.05

@test is_ODE_oscillatory(frequency_data, amplitude_data) == true

# Repressilator
repressilator = @reaction_network Repressilator begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ <--> ∅
    (δ,γ), m₂ <--> ∅
    (δ,γ), m₃ <--> ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end

pmap  = (:α => .5, :K => 40, :n => 2, :δ => log(2)/120,
         :γ => 5e-3, :β => log(2)/6, :μ => log(2)/60)
u₀map = [:m₁ => 43, :m₂ => 6.1, :m₃ => 6.5, :P₁ => 405, :P₂ => 97, :P₃ => 47]
tspan = (0., 10000.)
oprob = ODEProblem(repressilator, u₀map, tspan, pmap)
repressilator_sol = solve(oprob, Tsit5())

signal_sampling = Int(length(repressilator_sol.t))
spectrum_sampling = 2 * signal_sampling
amplitude_data = calculate_amplitude(repressilator_sol)
frequency_data = calculate_main_frequency(repressilator_sol, signal_sampling, spectrum_sampling)

@test abs(frequency_data["frequency"][1] - 0.00065)/0.00065 < 0.05
@test abs(frequency_data["frequency"][2] - 0.00065)/0.00065 < 0.05
@test abs(frequency_data["frequency"][3] - 0.00065)/0.00065 < 0.05
@test abs(frequency_data["frequency"][4] - 0.00065)/0.00065 < 0.05
@test abs(frequency_data["frequency"][5] - 0.00065)/0.00065 < 0.05
@test abs(frequency_data["frequency"][6] - 0.00065)/0.00065 < 0.05

@test abs(frequency_data["power"][1] - 1.8e6)/1.8e6 < 0.05
@test abs(frequency_data["power"][2] - 1.7e6)/1.7e6 < 0.05
@test abs(frequency_data["power"][3] - 1.8e6)/1.8e6 < 0.05
@test abs(frequency_data["power"][4] - 1.6e8)/1.6e8 < 0.05
@test abs(frequency_data["power"][5] - 1.5e8)/1.5e8 < 0.05
@test abs(frequency_data["power"][5] - 1.6e8)/1.6e8 < 0.05

@test abs(amplitude_data["amplitude"][1] - 40.2)/40.2 < 0.05
@test abs(amplitude_data["amplitude"][2] - 40.3)/40.3 < 0.05
@test abs(amplitude_data["amplitude"][3] - 40.3)/40.3 < 0.05
@test abs(amplitude_data["amplitude"][4] - 374.3)/374.3 < 0.05
@test abs(amplitude_data["amplitude"][5] - 374.2)/374.2 < 0.05
@test abs(amplitude_data["amplitude"][6] - 375.3)/375.3 < 0.05

@test abs(amplitude_data["peak_variation"][1] - 0.0099)/0.0099 < 0.05
@test abs(amplitude_data["peak_variation"][2] - 0.0056)/0.0056 < 0.05
@test abs(amplitude_data["peak_variation"][3] - 0.0078)/0.0078 < 0.05
@test abs(amplitude_data["peak_variation"][4] - 0.0081)/0.0081 < 0.05
@test abs(amplitude_data["peak_variation"][5] - 0.0098)/0.0098 < 0.05
@test abs(amplitude_data["peak_variation"][6] - 0.0065)/0.0065 < 0.05

@test abs(amplitude_data["trough_variation"][1] - 0.00077)/0.00077 < 0.05
@test abs(amplitude_data["trough_variation"][2] - 0.00093)/0.00093 < 0.05
@test abs(amplitude_data["trough_variation"][3] - 0.00058)/0.00058 < 0.05
@test abs(amplitude_data["trough_variation"][4] - 0.0011)/0.0011 < 0.05
@test abs(amplitude_data["trough_variation"][5] - 0.00069)/0.00069 < 0.05
@test abs(amplitude_data["trough_variation"][6] - 0.0011)/0.0011 < 0.05

@test is_ODE_oscillatory(frequency_data, amplitude_data) == true

# Steady state solutions should return have small power values and NaN amplitude
# Brusselator
p = [1.1, 1.2, 0.9, 1.3]
u0 = [1.0, 2.0]
tspan=(.0, 50.0)
oprob = ODEProblem(brusselator, u0, tspan, p)
equil_brusselator = solve(oprob, RadauIIA5())
oprob = remake(oprob, u0=equil_brusselator[end])
brusselator_sol = solve(oprob, RadauIIA5())

signal_sampling = Int(length(brusselator_sol.t))
spectrum_sampling = 100 * signal_sampling
amplitude_data = calculate_amplitude(brusselator_sol)
frequency_data = calculate_main_frequency(brusselator_sol, signal_sampling, spectrum_sampling)

@test mean(frequency_data["power"]) < 1e-5
@test all(isnan.(amplitude_data["amplitude"]))

@test is_ODE_oscillatory(frequency_data, amplitude_data) == false

# Repressilator
pmap  = (:α => 1, :K => 20, :n => 1, :δ => 0.5,
         :γ => 5e-3, :β => 0.5, :μ => 0.5)
u₀map = [:m₁ => 1.84, :m₂ => 1.84, :m₃ => 1.84, :P₁ => 1.84, :P₂ => 1.84, :P₃ => 1.84]
tspan = (0., 20.)
oprob = ODEProblem(repressilator, u₀map, tspan, pmap)
equil_repressilator = solve(oprob, Tsit5())
oprob = remake(oprob, u0=equil_repressilator[end])
repressilator_sol = solve(oprob, Tsit5(), abstol=1e-8, reltol=1e-6)

signal_sampling = Int(length(repressilator_sol.t))
spectrum_sampling = 2 * signal_sampling
amplitude_data = calculate_amplitude(repressilator_sol)
frequency_data = calculate_main_frequency(repressilator_sol, signal_sampling, spectrum_sampling)

@test mean(frequency_data["power"]) < 1e-5
@test all(isnan.(amplitude_data["amplitude"]))

@test is_ODE_oscillatory(frequency_data, amplitude_data) == false

# time scale separation

# sojourn time in nullclines, calculate_sojourn_time_fractions_in_nullclines
# parameter set #3 and #4 give spiky solutions

connectivity_interaction = [
    [1 0 -1; 1 0 0; 0 1 0],
    [1 0 -1; 1 0 0; 0 1 0],
    [1 0 -1; 1 0 0; 0 1 0],
    [1 0 -1; 1 0 0; 0 1 0],
    [1 0 -1; 1 0 0; 0 1 0],
    [0 0 -1; 1 0 0; 0 1 0],
    [0 0 -1; 1 0 0; 0 1 0],
    [0 0 -1; 1 0 0; 0 1 0],
    [0 0 -1; 1 0 0; 0 1 0],
    [0 0 -1; 1 0 0; 0 1 0],
]

α = [
    [1,0.215443,0.284804],
    [1,0.0278256,5.09414],
    [1,0.0159228,2.65609],
    [1,0.033516,0.453488],
    [1,0.599484,0.010975],
    [1,1.21958,2.1234],
    [1,0.0494845,0.0246626],
    [1,0.0214404,0.0779272],
    [1,0.135554,0.260941],
    [1,1.3571,0.0335469],
]

β = [
    [0.0305386,32.7455,3.51119],
    [0.070548,24.7708,17.0735],
    [14.1747,35.9381,4.64159],
    [29.8365,15.5568,57.2237],
    [1.14976,10.7227,8.90215],
    [0.0214009,11.8116,32.565],
    [0.06375,14.8155,30.3353],
    [0.614015,6.77147,24.8164],
    [6.09647,19.4057,37.1156],
    [0.0578197,30.6725,5.32933],
]

γ = [
    [159.228,291.505,7220.81,642.807],
    [559.081,422.924,1555.68,2364.49],
    [4977.02,2983.65,739.072,231.013],
    [1629.75,8302.18,4328.76,126.186],
    [486.26,9545.48,1291.55,774.264],
    [396.698,5055.57,234.119],
    [672.485,215.791,1255.19],
    [5447.2,8038.79,1207.56],
    [942.969,5462.28,1565.74],
    [336.863,4114.92,161.221],
]

κ = [
    [0.975758,0.450505,0.531313,0.765657],
    [0.862626,0.951515,0.692929,0.79798],
    [0.660606,0.733333,0.644444,0.749495],
    [0.668687,0.256566,0.507071,0.264646],
    [0.59596,0.216162,0.353535,0.757576],
    [0.824862,0.360336,0.816622],
    [0.526673,0.90015,0.586039],
    [0.827583,0.534513,0.589239],
    [0.631803,0.69861,0.20464],
    [0.947435,0.246405,0.524912],
]

η = [
    [4.79798,3.62626,3.90909,3.26263],
    [3.0202,4.75758,3.74747,3.54545],
    [1.88889,4.31313,1.28283,3.70707],
    [2.33333,4.07071,3.42424,3.74747],
    [4.11111,4.47475,2.29293,4.39394],
    [4.33993,3.81748,2.94179],
    [2.13491,4.39114,4.79118],
    [3.78468,3.38704,2.15292],
    [4.46355,2.20132,4.17912],
    [2.50375,3.12861,4.76478],
]

tspan = (0.0, 2.0)
initial_condition = [0.6, 0.6, 0.6]

max_sojourn_time = []

for (i, c) in enumerate(connectivity_interaction)
    model_interaction = protein_interaction_network(c)
    params = pin_parameters(model_interaction, α[i], β[i], γ[i], κ[i], η[i])
    equilibration = ODEProblem(model_interaction, initial_condition, tspan, params)
    sol_eq = solve(equilibration)
    ode_problem = ODEProblem(model_interaction, sol_eq.u[end], tspan, params)
    sol = solve(ode_problem)

    push!(max_sojourn_time, maximum(calculate_sojourn_time_fractions_in_nullclines(sol)))
end

@test max_sojourn_time[3] > max_sojourn_time[1]
@test max_sojourn_time[3] > max_sojourn_time[2]
@test max_sojourn_time[3] > max_sojourn_time[5]
@test max_sojourn_time[3] > max_sojourn_time[6]
@test max_sojourn_time[3] > max_sojourn_time[7]
@test max_sojourn_time[3] > max_sojourn_time[8]
@test max_sojourn_time[3] > max_sojourn_time[9]
@test max_sojourn_time[3] > max_sojourn_time[10]

@test max_sojourn_time[4] > max_sojourn_time[1]
@test max_sojourn_time[4] > max_sojourn_time[2]
@test max_sojourn_time[4] > max_sojourn_time[5]
@test max_sojourn_time[4] > max_sojourn_time[6]
@test max_sojourn_time[4] > max_sojourn_time[7]
@test max_sojourn_time[4] > max_sojourn_time[8]
@test max_sojourn_time[4] > max_sojourn_time[9]
@test max_sojourn_time[4] > max_sojourn_time[10]
