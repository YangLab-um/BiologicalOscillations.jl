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