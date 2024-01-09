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


"""
    kimchi_2020
Posttranslational oscillator based on [Kimchi et al. 2020](https://www.science.org/doi/10.1126/sciadv.abc1939)
"""
kimchi_2020 = @reaction_network Kimchi_2020 begin
    @species K(t)=0.5 P(t)=0.5
    # Parameter values estimated from Fig S2A
    @parameters p_tilde=1e-9 k_bk=1.0 k_tot=1.0 n=2.0 eta_k=1.0 k_uk=1e-3 k_bp=1e-1 p_tot=7.0 m=2.0 eta_p=1.0 k_up=1e-1
    # Kinase reactions
    k_bk * ((k_tot - n * K) / (1 + eta_k * (K / (P + p_tilde)))) ^ n, ∅ --> K
    k_uk, K --> ∅
    # Phosphatase reactions
    k_bp * ((p_tot - m * P) / (1 + eta_p * (K / (P + p_tilde)))) ^ m, ∅ --> P
    k_up, P --> ∅
end