# Set up analysis varying parameters one at a time
include("GCD_R0_num_calc.jl")

# Function to calculate GCD and R0 for exponential case
# (lQ, lL, lP, lG, sigma, pL, pP, pG) = B_vals_in
function GCD_func_exp(b)
    return(1/b)
end

function R0_func_exp(b)
    # (lQ, lL, lP, lG, sigma, pL, pP, pG) = B_vals_in
    # theta = 1/b
    J_eq = KJ * (1 - (varPhi * (rhoJ / (rhoJ + muJ) * (b / (b + mu)) * (1 / (gV + mu)) * (1 - (b / (b + mu)) * (gV / (gV + mu)) * (gR / (gR + mu)))^(-1)))^(-1))
    V_eq = (1 - (b / (b + mu)) * (gV / (gV + mu)) * (gR / (gR + mu)))^(-1) * (b / (b + mu)) * (1 / (gV + mu)) * rhoJ * J_eq 
    B_eq = ((gV + mu) / b) * V_eq

    F_mat = [
        0            0 bH*b 0 0 0 0 
        0            0 0    0 0 0 0
        0            0 0    0 0 0 0
        bB*b*B_eq/KH 0 0    0 0 0 0
        0            0 0    0 0 0 0
        0            0 0    0 0 0 0
        0            0 0    0 0 0 0        
    ]

    V_mat = [
        gH+muH 0        0    0         0     0         0 
        0      eta+b+mu 0    0         0     -gR       0
        0      -eta     b+mu 0         0     0         -gR
        0      -b       0    eta+gV+mu 0     0         0
        0      0        -b   -eta      gV+mu 0         0
        0      0        0    -gV       0     eta+gR+mu 0
        0      0        0    0         -gV   -eta      gR+mu
    ]

    K_mat = F_mat * inv(V_mat)

    R0 = maximum(real(eigen(K_mat).values))
    return(R0)    
end

# Set up ranges to vary each parameter
const base_params_flighty = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.9f0, 0.5f0, 0.5f0, 0.5f0]

const base_params_persistent = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.66f0, 0.7f0, 0.8f0, 0.9f0]

const theta_range = [0f0, 20f0 * 1440] # minutes
# Set up data frame varying each parameter one by one, keeping track of the varied parameter

theta_vec = LinRange(480f0, 2/mu, 101)
b_vec = [1/v for v in theta_vec]



# Calculate GCD and R0 across the varied parameters

# Save as CSV