"""
    elowitz_2000

Generalized Repressilator model based on [Elowitz and Leibler 2000 Nature article](https://www.nature.com/articles/35002125). 
"""
elowitz_2000 = @reaction_network Elowitz_2000 begin
    hillr(P₃,α₃,K₃,n₃), ∅ --> m₁
    hillr(P₁,α₁,K₁,n₁), ∅ --> m₂
    hillr(P₂,α₂,K₂,n₂), ∅ --> m₃
    (δ₁,γ₁), m₁ <--> ∅
    (δ₂,γ₂), m₂ <--> ∅
    (δ,γ₃), m₃ <--> ∅
    β₁, m₁ --> m₁ + P₁
    β₂, m₂ --> m₂ + P₂
    β₃, m₃ --> m₃ + P₃
    μ₁, P₁ --> ∅
    μ₂, P₂ --> ∅
    μ₃, P₃ --> ∅
end 


"""
    guan_2008

Embryonic cell cycle model based on [Guan et al. 2008 eLife article](https://elifesciences.org/articles/33549). 
"""
guan_2008 = @reaction_network Guan_2008 begin
    @species B(t)=70.0 C(t)=35.0
    @parameters k_S=1.0 a_D=0.01 b_D=0.04 n_D=17.0 K_D=32.0 a_T=0.16 b_T=0.8 n_T=11.0 K_T=30.0 a_W=0.08 b_W=0.4 n_W=3.5 K_W=35.0 r=1.5

    # Synthesis
    k_S, ∅ --> B
    k_S, ∅ --> C
    # Degradation 
    a_D + hill(abs(C), b_D, K_D, n_D), B --> ∅
    a_D + hill(abs(C), b_D, K_D, n_D), C --> ∅
    # Cdk1 activation by Cdc25
    (1/sqrt(r))*(a_T + hill(abs(C), b_T, K_T, n_T)), B --> B + C
    (1/sqrt(r))*(a_T + hill(abs(C), b_T, K_T, n_T)), C --> ∅
    # Cdk1 deactivation by Wee1
    sqrt(r)*(a_W + hillr(abs(C), b_W, K_W, n_W)), C --> ∅
end