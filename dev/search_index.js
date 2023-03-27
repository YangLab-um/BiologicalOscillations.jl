var documenterSearchIndex = {"docs":
[{"location":"gene_regulatory_network/grn_introduction/#grn_introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"models/elowitz_2000/#elowitz_2000","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"","category":"section"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"elowitz_2000","category":"page"},{"location":"models/elowitz_2000/#BiologicalOscillations.elowitz_2000","page":"Elowitz & Leibler 2000 - Repressilator","title":"BiologicalOscillations.elowitz_2000","text":"elowitz_2000\n\nGeneralized Repressilator model based on Elowitz and Leibler 2000 Nature article. \n\n\n\n\n\n","category":"constant"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"using Catalyst, DifferentialEquations, Plots, Latexify\n\nelowitz_2000 = @reaction_network Elowitz_2000 begin\n    hillr(P₃,α₃,K₃,n₃), ∅ --> m₁\n    hillr(P₁,α₁,K₁,n₁), ∅ --> m₂\n    hillr(P₂,α₂,K₂,n₂), ∅ --> m₃\n    (δ₁,γ₁), m₁ <--> ∅\n    (δ₂,γ₂), m₂ <--> ∅\n    (δ₃,γ₃), m₃ <--> ∅\n    β₁, m₁ --> m₁ + P₁\n    β₂, m₂ --> m₂ + P₂\n    β₃, m₃ --> m₃ + P₃\n    μ₁, P₁ --> ∅\n    μ₂, P₂ --> ∅\n    μ₃, P₃ --> ∅\nend\n\nstruct LaTeXEquation\n    content::String\nend\n\nfunction Base.show(io::IO, ::MIME\"text/latex\", x::LaTeXEquation)\n    # Wrap in $$ for display math printing\n    return print(io, \"\\$\\$ \" * x.content * \" \\$\\$\")\nend\n\nLatexify.set_default(; starred=true)\n\npmap  = (:α₁ => 5e-1, :α₂ => 5e-1, :α₃ => 5e-1,\n         :K₁ => 40, :K₂ => 40, :K₃ => 40,  \n         :n₁ => 2, :n₂ => 2, :n₃ => 2, \n         :δ₁ => 2.5e-3, :δ₂ => 2.5e-3, :δ₃ => 2.5e-3,\n         :γ₁ => 5e-3, :γ₂ => 5e-3, :γ₃ => 5e-3, \n         :β₁ => 5e-2, :β₂ => 5e-2, :β₃ => 5e-2, \n         :μ₁ => 5e-3, :μ₂ => 5e-3, :μ₃ => 5e-3)\nu₀map = [:m₁ => 0., :m₂ => 0., :m₃ => 0., :P₁ => 20., :P₂ => 0., :P₃ => 0.]\ntspan = (0., 10000.)\noprob = ODEProblem(elowitz_2000, u₀map, tspan, pmap)\nsol = solve(oprob, Tsit5(), saveat=10.)","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"In contrast to the original model, this implementation has independent parameters for each mRNA and protein pair:","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"odesys = convert(ODESystem, elowitz_2000) # hide\neq = latexify(odesys) # hide\nLaTeXEquation(eq) # hide","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"where ","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"texthillr(P(t) alpha K n) = alpha fracK^nP(t)^n + K^n","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"We can obtain an oscillatory solution by using the following parameters ","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"       \nalpha_1 alpha_2 alpha_3 0.5 K_1 K_2 K_3 40 beta_1 beta_2 beta_3 0.05 n_1 n_2 n_3 2\ngamma_1 gamma_2 gamma_3 0.005 mu_1 mu_2 mu_3 0.005 delta_1 delta_2 delta_3 0.0025  ","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"and the initial condition","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"     \nm_1(0) m_2(0) m_3(0) 0.0 P_1(0) 20.0 P_2(0) P_3(0) 0.0","category":"page"},{"location":"models/elowitz_2000/","page":"Elowitz & Leibler 2000 - Repressilator","title":"Elowitz & Leibler 2000 - Repressilator","text":"plot(sol) # hide","category":"page"},{"location":"protein_interaction_network/pin_introduction/#pin_introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"protein_interaction_network/pin_introduction/","page":"Introduction","title":"Introduction","text":"BiologicalOscillations.jl contains a set of functions implemented to study the oscillatory behavior of interacting proteins. Through the function protein_interaction_network it is possible to generate a system of differential equations that models the dynamics of a collection of interacting proteins. The function has a single input encoding how every node of the network is connected to each other. For example:","category":"page"},{"location":"protein_interaction_network/pin_introduction/","page":"Introduction","title":"Introduction","text":"using BiologicalOscillations\n\nconnectivity = [0 0 -1; -1 0 0; 0 -1 0]\nmodel = protein_interaction_network(connectivity)","category":"page"},{"location":"protein_interaction_network/pin_introduction/","page":"Introduction","title":"Introduction","text":"creating a system of equations for a Repressilator-like network. See Model details for more information about the generated equations.","category":"page"},{"location":"gene_regulatory_network/grn_documentation/#grn_documentation","page":"Documentation","title":"Documentation","text":"","category":"section"},{"location":"gene_regulatory_network/grn_documentation/","page":"Documentation","title":"Documentation","text":"Modules = [BiologicalOscillations]\nPages = [\"gene_regulatory_network.jl\"]","category":"page"},{"location":"gene_regulatory_network/grn_documentation/#BiologicalOscillations.gene_regulatory_network-Tuple{AbstractMatrix{T} where T}","page":"Documentation","title":"BiologicalOscillations.gene_regulatory_network","text":"gene_regulatory_network(connectivity::AbstractMatrix)\n\nCreates a ReactionSystem of interacting mRNA and proteins based on the provided connectivity.\n\nArguments\n\nconnectivity::AbstractMatrix: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network.\n\n\n\n\n\n","category":"method"},{"location":"gene_regulatory_network/grn_documentation/#BiologicalOscillations.grn_parameters-Tuple{Catalyst.ReactionSystem, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T}","page":"Documentation","title":"BiologicalOscillations.grn_parameters","text":"grn_parameters(model:ReactionSystem, α::AbstractVector, β::AbstractVector, δ::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)\n\nCreates an ordered parameter vector for a model created with the gene_regulatory_network function to use in ODEProblem.\n\nArguments\n\nmodel::ReactionSystem: Model generated with protein_interaction_network\nα:AbstractVector: Vector of intrinsic mRNA degradation rates\nβ:AbstractVector: Vector of intrinsic protein synthesis rates\nδ:AbstractVector: Vector of intrinsic protein degradation rates\nγ:AbstractVector: Vector of strengths for each network interaction\nκ:AbstractVector: Vector of midpoints for each network interaction\nη:AbstractVector: Vector of sensitivity values for each network interaction\n\nNote\n\nIt is assumed that parameters follow the same order as the connectivity matrix. Namely, the first node is encoded on the first row of the connectivty matrix and the first edge comes from the first nonzero element of the connectivity.\n\n\n\n\n\n","category":"method"},{"location":"gene_regulatory_network/grn_documentation/#BiologicalOscillations.grn_timescale-Tuple{AbstractVector{T} where T, AbstractVector{T} where T}","page":"Documentation","title":"BiologicalOscillations.grn_timescale","text":"grn_timescale(α::AbstractVector, δ::AbstractVector)\n\nCalculate the slowest timescale of the gene regulatory network\n\nArguments\n\nα:AbstractVector: Vector of intrinsic mRNA degradation rates\nδ:AbstractVector: Vector of intrinsic protein degradation rates\n\nReturns\n\ntimescale::Real\n\n\n\n\n\n","category":"method"},{"location":"feature_calculation/fc_documentation/#fc_documentation","page":"Documentation","title":"Documentation","text":"","category":"section"},{"location":"feature_calculation/fc_documentation/","page":"Documentation","title":"Documentation","text":"Modules = [BiologicalOscillations]\nPages = [\"feature_calculation.jl\"]","category":"page"},{"location":"feature_calculation/fc_documentation/#BiologicalOscillations.calculate_amplitude-Tuple{SciMLBase.ODESolution}","page":"Documentation","title":"BiologicalOscillations.calculate_amplitude","text":"calculate_amplitude(ode_solution::ODESolution)\n\nCalculate the amplitude of each signal in ode_solution and the relative peak/trough variation.\n\nArguments\n\node_solution::ODESolution: The full output of the solve() function from DifferentialEquations.jl\n\nReturns\n\namplitude_data::Dict: Dictionary containing the amplitude for each signal and the peak/trough variation. The output is encoded as: \n\namplitude_data = Dict('amplitude' => [amp_signal_1, ..., amp_signal_N], \n                      'peak_variation' => [peak_var_signal_1, ..., peak_var_signal_N], \n                      'trough_variation' => [trough_var_signal_1, ..., trough_var_signal_N])\n\nNotes\n\nVariation of peaks and troughs is calculated as the standard deviation of all the peaks/troughs in the signal divided by the amplitude. This quantity allows to distinguish stable oscillatory solution from dampened ones.\nTo find the most relevant peaks/troughs, we discard any peak/trough that has a prominence less than 1/3 of the total amplitude of the signal (calculated as the maximum - minimum of the signal). Therefore, it is better to use this function with a signal that has little to none transient behavior.\n\n\n\n\n\n","category":"method"},{"location":"feature_calculation/fc_documentation/#BiologicalOscillations.calculate_main_frequency-Tuple{SciMLBase.ODESolution, Int64, Int64}","page":"Documentation","title":"BiologicalOscillations.calculate_main_frequency","text":"    calculate_main_frequency(ode_solution::ODESolution, sampling::Int, fft_points::Int)\n\nObtains the dominant frequencies and their powers for each signal on ode_solution via the fast fourier transform\n\nArguments\n\node_solution::ODESolution: The full output of the solve() function from DifferentialEquations.jl\nsampling::Int: Number of points to be used for uniform re-sampling of ode_solution. The same sampling is used for all signals\nfft_points::Int: Number of points in the frequency spectrum\n\nReturns\n\nfrequency_data::Dict: Dictionary containing the dominant frequencies and their powers for each signal. Here power means the amplitude of that frequency on the FFT spectrum. The output is encoded as:\n\nfrequency_data = Dict('frequency' => [freq_signal_1, ..., freq_signal_N], \n                      'power' => [power_signal_1, ..., power_signal_N])\n\nNote\n\nThe provided ode_solution has to be solved with the dense = true setting to be suitable for this function since uniform re-sampling is performed.\n\n\n\n\n\n","category":"method"},{"location":"feature_calculation/fc_documentation/#BiologicalOscillations.is_ODE_oscillatory","page":"Documentation","title":"BiologicalOscillations.is_ODE_oscillatory","text":"is_ODE_oscillatory(frequency_data::Dict, amplitude_data::Dict, freq_variation_threshold=0.05, power_threshold=1e-7, amp_variation_threshold=0.01)\n\nReturns true if an ODESolution is oscillatory based on the calculated frequency and amplitude. False otherwise.\n\nArguments\n\nfrequency_data::Dict: The output from calculate_main_frequency\namplitude_data::Dict: The output from calculate_amplitude\nfreq_variation_threhsold::Real: Maximum tolerated variation in frequency between species in the solution to be declared oscillatory. Default value 0.05.\npower_threshold::Real: Minimum spectral power that the main peak has to have to be declared oscillatory. Default value 1e-7.\namp_variation_threhsold::Real: Maximum tolerated value for peak/trough variation to be declared oscillatory. Default value 0.01.\n\n\n\n\n\n","category":"function"},{"location":"gene_regulatory_network/grn_model_details/#grn_model_details","page":"Model details","title":"Model details","text":"","category":"section"},{"location":"gene_regulatory_network/grn_model_details/","page":"Model details","title":"Model details","text":"using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify, Plots\n\nconnectivity_interaction = [0 0 1; 1 0 0; 0 -1 0]\nmodel_interaction = gene_regulatory_network(connectivity_interaction)\n\ntspan = (0.0, 0.5)\ninitial_condition = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6]\nα = [1.0, 0.13, 2.3]\nβ = [44.3, 75.3, 0.10]\nδ = [44.3, 75.3, 0.10]\nγ = [2472.4, 442.2, 4410.0]\nκ = [0.27, 0.91, 0.29]\nη = [3.0, 4.8, 3.3]\nparams = grn_parameters(model_interaction, α, β, δ, γ, κ, η)\n\node_problem = ODEProblem(model_interaction, initial_condition, tspan, params)\nsolution = solve(ode_problem, RadauIIA5())\nplot(solution)","category":"page"},{"location":"#BiologicalOscillations.jl-Documentation","page":"Home","title":"BiologicalOscillations.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"BiologicalOscillations is a computational package for researchers studying biological oscillations. The package implements published mathematical models of oscillatory phenomena in a wide variety of areas such as cell cycle dynamics, circadian rhythms, synthetic biology, and embryonic development. The implementation leverages several tools across the Julia SciML ecosystem such as Catalyst.jl, DSP.jl, and DifferentialEquations.jl. ","category":"page"},{"location":"protein_interaction_network/pin_documentation/#pin_documentation","page":"Documentation","title":"Documentation","text":"","category":"section"},{"location":"protein_interaction_network/pin_documentation/","page":"Documentation","title":"Documentation","text":"Modules = [BiologicalOscillations]\nPages = [\"protein_interaction_network.jl\"]","category":"page"},{"location":"protein_interaction_network/pin_documentation/#BiologicalOscillations.pin_parameters-Tuple{Catalyst.ReactionSystem, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T}","page":"Documentation","title":"BiologicalOscillations.pin_parameters","text":"pin_parameters(model:ReactionSystem, α::AbstractVector, β::AbstractVector, γ::AbstractVector, κ::AbstractVector, η::AbstractVector)\n\nCreates an ordered parameter vector for a model created with the protein_interaction_network function to use in ODEProblem.\n\nArguments\n\nmodel::ReactionSystem: Model generated with protein_interaction_network\nα:AbstractVector: Vector of intrinsic protein activation rates for each node\nβ:AbstractVector: Vector of intrinsic protein deactivation rates for each node\nγ:AbstractVector: Vector of strengths for each network interaction\nκ:AbstractVector: Vector of midpoints for each network interaction\nη:AbstractVector: Vector of sensitivity values for each network interaction\n\nNote\n\nIt is assumed that parameters follow the same order as the connectivity matrix. Namely, the first node is encoded on the first row of the connectivty matrix and the first edge comes from the first nonzero element of the connectivity.\n\n\n\n\n\n","category":"method"},{"location":"protein_interaction_network/pin_documentation/#BiologicalOscillations.pin_timescale-Tuple{AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{T} where T}","page":"Documentation","title":"BiologicalOscillations.pin_timescale","text":"pin_timescale(α::AbstractVector, β::AbstractVector, γ::AbstractVector)\n\nCalculate the slowest timescale of the protein interaction network\n\nArguments\n\nα:AbstractVector: Vector of intrinsic protein activation rates for each node\nβ:AbstractVector: Vector of intrinsic protein deactivation rates for each node\nγ:AbstractVector: Vector of strengths for each network interaction\n\nReturns\n\ntimescale::Real\n\n\n\n\n\n","category":"method"},{"location":"protein_interaction_network/pin_documentation/#BiologicalOscillations.protein_interaction_network-Tuple{AbstractMatrix{T} where T}","page":"Documentation","title":"BiologicalOscillations.protein_interaction_network","text":"protein_interaction_network(connectivity::AbstractMatrix)\n\nCreates a ReactionSystem of interacting proteins based on the provided connectivity.\n\nArguments\n\nconnectivity::AbstractMatrix: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network.\n\n\n\n\n\n","category":"method"},{"location":"protein_interaction_network/pin_model_details/#pin_model_details","page":"Model details","title":"Model details","text":"","category":"section"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"The protein_interaction_network functionality of BiologicalOscillations.jl stems from the work of Tsai et al. 2008 and Li et al. 2017 where a general model of interacting proteins was used to study the robustness and tunability of biological oscillators. Proteins in this model have two possible states: active or inactive, and interactions with other proteins affect the rates at which proteins convert between these two states. ","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"Every node has an intrinsic activation rate determined by the parameter alpha and an inactivation rate given by beta:","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"using Latexify\n\nstruct LaTeXEquation\n    content::String\nend\n\nfunction Base.show(io::IO, ::MIME\"text/latex\", x::LaTeXEquation)\n    # Wrap in $$ for display math printing\n    return print(io, \"\\$\\$ \" * x.content * \" \\$\\$\")\nend\n\nLatexify.set_default(; starred=true)","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify\n\nconnectivity_no_interaction = [0 0 0; 0 0 0; 0 0 0]\nmodel_no_interaction = protein_interaction_network(connectivity_no_interaction)\nodesys_no_interaction = convert(ODESystem, model_no_interaction)\nlatexify(odesys_no_interaction)\neq = latexify(odesys_no_interaction) # hide\nLaTeXEquation(eq) # hide","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"In the example above we are creating a network of 3 nodes (notice the 3x3 size of connectivity) without any interaction between nodes (all zeros in the connectivity matrix). Thus, each node only contains their intrinsic activation and deactivation rates. To display the created equations we convert our ReactionSystem into an ODESystem and print the equations with latexify(). The model only tracks the fraction of each protein on its active state. Therefore, each variable in the differential equations (X(t)_i) is bounded to be between 0 and 1.","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"When interactions between nodes are added, additional terms affect the rate at which the fraction of active protein evolves over time. These interactions between nodes can be positive or negative and modify the equations in the following way:","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify\n\nconnectivity_interaction = [0 0 1; 1 0 0; 0 -1 0]\nmodel_interaction = protein_interaction_network(connectivity_interaction)\nodesys_interaction = convert(ODESystem, model_interaction)\nlatexify(odesys_interaction)\neq = latexify(odesys_interaction) # hide\nLaTeXEquation(eq) # hide","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"where","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"texthill(X(t)_i gamma_j kappa_j eta_j) = gamma_j fracX(t)_i^eta_jX(t)_i^eta_j + kappa_j^eta_j","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"with i counting over each node and j each edge of the network. Every interaction has three parameters: gamma_j encodes the strength of the interaction while kappa_j and eta_j determine the midpoint and sensitivity of the response respectively. Here we have created a Goodwin-like oscillator by simply feeding protein_interaction_network with the desired connectivity.","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"Thanks to ModelingToolkit.jl and Catalyst.jl the model can be simulated by simply converting it to an ODEProblem and specifying the network parameters and initial condition. Network parameters should be specified using the function pin_parameters to ensure that values are assigned in the correct order:","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"using BiologicalOscillations, Catalyst, DifferentialEquations, Latexify, Plots\n\nconnectivity_interaction = [0 0 1; 1 0 0; 0 -1 0]\nmodel_interaction = protein_interaction_network(connectivity_interaction)\n\ntspan = (0.0, 0.5)\ninitial_condition = [0.6, 0.6, 0.6]\nα = [1.0, 0.13, 2.3]\nβ = [44.3, 75.3, 0.10]\nγ = [2472.4, 442.2, 4410.0]\nκ = [0.27, 0.91, 0.29]\nη = [3.0, 4.8, 3.3]\nparams = pin_parameters(model_interaction, α, β, γ, κ, η)\n\node_problem = ODEProblem(model_interaction, initial_condition, tspan, params)\nsolution = solve(ode_problem, RadauIIA5())\nplot(solution)","category":"page"},{"location":"protein_interaction_network/pin_model_details/","page":"Model details","title":"Model details","text":"It is recommended to set alpha_1 = 10 to normalize the time units to work with a dimensionless model both in time and concentration.","category":"page"}]
}
