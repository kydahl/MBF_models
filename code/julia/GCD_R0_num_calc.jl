# Numerical calculators for key outputs
# Load packages
using IntervalSets
using LinearAlgebra
using Symbolics
using Serialization
using DataFrames
using CSV
using ProgressBars

# !!! put these in a separate script called by this and GCD_R0_symbolic_calc.jl
const (KJ, rhoJ, muJ, varPhi) = [3E8, (1/(12 * 1440)), (1 / (20 * 1440)), (300 / 1440)]
const (gV, gR, mu) = [(1/(2 * 1440)), (1/(2 * 1440)), (1/(20 * 1440))]
const (bH, bB, eta, gH) = [1, 1, (1/(6 * 1440)), (1/(7 * 1440))]
const (muH, KH) = [1/(365.25 * 65 * 1440), 1E8]
const pQ = 1.0f0
const (bH, bB, eta, gH) = [1f0, 1f0, 1 / (6 * 1440), 1 / (7 * 1440)]

# Load in output functions
GCD_func = Serialization.deserialize("GCD_func.jls");
N_offspring_func = Serialization.deserialize("N_offspring_func.jls");
# R0_func = Serialization.deserialize("R0_func.jls");

function R0_func(B_vals_in, N_offspring_in)
    (lQ, lL, lP, lG, sigma, pL, pP, pG) = B_vals_in

    if N_offspring_in < 1
        return(0)
    else 
    one_vec_five = [1.0f0 1.0f0 1.0f0 1.0f0 1.0f0]
    alpha_vec_five = [1.0f0; 0.0f0; 0.0f0; 0.0f0; 0.0f0]
    one_vec_four = [1.0f0 1.0f0 1.0f0 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    f = 1 - sigma
    # Define subintensity matrices
    A_mat = [
        -lQ         lQ                  0f0   0f0   
        f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0f0 
        f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP   pP*lP
        f*(1-pG)*lG (1-f)*(1-pG)*lG     0f0   -lG
    ]

    A_tilde = [
        -lQ         lQ                  0f0   0f0   0f0
        f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0f0   0f0
        f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP   pP*lP 0f0
        f*(1-pG)*lG (1-f)*(1-pG)*lG     0f0   -lG   pG*lG
        0f0         0f0                 0f0   0f0   -gV
    ]

    r = (N_offspring_in - 1) * KJ * ((rhoJ + muJ)/ varPhi) * (mu + gV)
    # Distribution of mosquitoes across states at equilibrium
    B_prefactor = (N_offspring_in - 1) * KJ * ((rhoJ / N_offspring_in) + (gR / (mu + gR) * gV * (rhoJ + muJ) / varPhi))
    temp_inv = simplify(inv(mu * I - transpose(A_mat)))
    B_postfactor = temp_inv * alpha_vec_four

    B_star = B_prefactor * B_postfactor
    V_star = r / (mu + gV)

    B_vec = [B_star; V_star]
    B_tot = sum(B_star)

    betaH_mat = zeros(Float64, 5, 5); betaV_mat = zeros(Float64, 5, 5); LambdaH = zeros(Float64, 5, 5); LambdaV = zeros(Float64, 5, 5)
    # Rate of contact is rate of entrance into transmission compartments
    LambdaH[3,3] = lP * B_tot / KH
    LambdaV[4,4] = lG
    betaH_mat[3,3] = bH
    betaV_mat[4,4] = bB

    spec_mat = -alpha_vec_five * one_vec_five * transpose(A_tilde)
    term_1 = mu * I - transpose(A_tilde) - (gR / (gR + mu)) * spec_mat
    term_2 = eta * I - (eta / (eta + gR + mu)) * (gR / (gR + mu)) * spec_mat
    term_3 = (eta + mu) * I - transpose(A_tilde) - (gR / (eta + gR + mu) * spec_mat)

    R02 = (1 / (gH + muH)) * one_vec_five * betaH_mat * LambdaH * inv(term_1) * term_2 * inv(term_3) * betaV_mat * LambdaV * (B_vec / (B_tot))
    R0 = sqrt(R02[1])
    return(R0)
    end
end

# Define parameter ranges

# # Set resolution of variation
# resolution = 10

# # Probabilities
# generic_prob_vec = range(0f0, 1f0, resolution)
# sigma_vec = generic_prob_vec
# pQ_vec = 1.0f0
# pL_vec = generic_prob_vec
# pP_vec = generic_prob_vec
# pG_vec = generic_prob_vec

# Calculate key outputs over the parameter combinations

# Run each function once to pre-load it into Julia

GCD_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
N_offspring_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
temp = R0_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0], N_offspring_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0]))