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