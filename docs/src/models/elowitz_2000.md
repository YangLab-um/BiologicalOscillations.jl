# [Elowitz & Leibler 2000 - Repressilator](@id elowitz_2000)
```@docs
elowitz_2000
```
```@setup elowitz_2000
using Catalyst, DifferentialEquations, Plots, Latexify

elowitz_2000 = @reaction_network Elowitz_2000 begin
    hillr(P₃,α₃,K₃,n₃), ∅ --> m₁
    hillr(P₁,α₁,K₁,n₁), ∅ --> m₂
    hillr(P₂,α₂,K₂,n₂), ∅ --> m₃
    (δ₁,γ₁), m₁ <--> ∅
    (δ₂,γ₂), m₂ <--> ∅
    (δ₃,γ₃), m₃ <--> ∅
    β₁, m₁ --> m₁ + P₁
    β₂, m₂ --> m₂ + P₂
    β₃, m₃ --> m₃ + P₃
    μ₁, P₁ --> ∅
    μ₂, P₂ --> ∅
    μ₃, P₃ --> ∅
end

struct LaTeXEquation
    content::String
end

function Base.show(io::IO, ::MIME"text/latex", x::LaTeXEquation)
    # Wrap in $$ for display math printing
    return print(io, "\$\$ " * x.content * " \$\$")
end

Latexify.set_default(; starred=true)

pmap  = (:α₁ => 5e-1, :α₂ => 5e-1, :α₃ => 5e-1,
         :K₁ => 40, :K₂ => 40, :K₃ => 40,  
         :n₁ => 2, :n₂ => 2, :n₃ => 2, 
         :δ₁ => 2.5e-3, :δ₂ => 2.5e-3, :δ₃ => 2.5e-3,
         :γ₁ => 5e-3, :γ₂ => 5e-3, :γ₃ => 5e-3, 
         :β₁ => 5e-2, :β₂ => 5e-2, :β₃ => 5e-2, 
         :μ₁ => 5e-3, :μ₂ => 5e-3, :μ₃ => 5e-3)
u₀map = [:m₁ => 0., :m₂ => 0., :m₃ => 0., :P₁ => 20., :P₂ => 0., :P₃ => 0.]
tspan = (0., 10000.)
oprob = ODEProblem(elowitz_2000, u₀map, tspan, pmap)
sol = solve(oprob, Tsit5(), saveat=10.)
```

In contrast to the original model, this implementation has independent parameters for each mRNA and protein pair:

```@example elowitz_2000
odesys = convert(ODESystem, elowitz_2000) # hide
eq = latexify(odesys) # hide
LaTeXEquation(eq) # hide
```
where 
```math
\text{hillr}(P(t), \alpha, K, n) = \alpha \frac{K^n}{P(t)^n + K^n}
```
Using parameters 

```math
\begin{aligned}
\alpha_1 = \alpha_2 = \alpha_3 &= 0.5 \\
\beta_1 = \beta_2 = \beta_3 &= 0.05 \\
\gamma_1 = \gamma_2 = \gamma_3 &= 0.005 \\
\delta_1 = \delta_2 = \delta_3 &= 0.0025 \\
K_1 = K_2 = K_3 &= 40 \\
n_1 = n_2 = n_3 &= 2 \\
\mu_1 = \mu_2 = \mu_3 &= 0.005
\end{aligned}
```

Under the initial condition

```math
\begin{aligned}
m_1(0) = m_2(0) = m_3(0) &= 0.0 \\
P_1(0) &= 20.0 \\
P_2(0) = P_3(0) &= 0.0
\end{aligned}
```
we obtain:

```@example elowitz_2000
plot(sol) # hide
```
