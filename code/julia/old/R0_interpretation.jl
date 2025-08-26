using Symbolics
using LinearAlgebra

@variables μ γᵥ η γᵣ 

V_VR = [
    (μ+γᵥ+η) 0      0        0
    -η       (μ+γᵥ) 0        0
    -γᵥ      0      (μ+γᵣ+η) 0
    0        -γᵥ    -η       (μ+γᵣ)
]

inv_V_VR = inv(V_VR)