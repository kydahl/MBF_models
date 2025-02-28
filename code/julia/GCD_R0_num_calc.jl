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
const (KJ, rhoJ, muJ, varPhi) = [0.75 * 1E3, (1/(12 * 1440)), (1 / (20 * 1440)), (3 / 1440)]
const (gV, gR, mu) = [(1/(5 * 1440)), (1/(2 * 1440)), (1/(21 * 1440))]
const (muH, KH) = [1/(365.25 * 65 * 1440), 1E3]
# const pQ = 1.0f0
const (bH, bB, eta, gH) = [0.5f0, 0.5f0, 1 / (7 * 1440), 1 / (5 * 1440)]

# Load in output functions
GCD_func = Serialization.deserialize("GCD_func.jls");
N_offspring_func = Serialization.deserialize("N_offspring_func.jls");
# R0_func = Serialization.deserialize("R0_func.jls");

function R0_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    # one_vec_five = [1.0f0 1.0f0 1.0f0 1.0f0 1.0f0]
    # alpha_vec_five = [1.0f0; 0.0f0; 0.0f0; 0.0f0; 0.0f0]
    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    f = 1 - sigma
    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ         pQ*lQ                  0f0   0f0   
        f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0f0 
        f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP   pP*lP
        f*(1-pG)*lG (1-f)*(1-pG)*lG     0f0   -lG
    ]

    # A_tilde = [
    #     -lQ         lQ                  0f0   0f0   0f0
    #     f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0f0   0f0
    #     f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP   pP*lP 0f0
    #     f*(1-pG)*lG (1-f)*(1-pG)*lG     0f0   -lG   pG*lG
    #     0f0         0f0                 0f0   0f0   -gV
    # ]

    tau = transpose(-A_mat * one_vec_four) * inv(mu * I - transpose(A_mat)) * alpha_vec_four
    tau = tau[1]
    rho = (gV / (mu + gV)) * (gR / (mu + gR)) * tau
    nG = 1.0f0 / (1.0f0 - rho)
    # Basic offspring number
    N_offspring = tau * (varPhi / (mu + gV)) * (rhoJ / (rhoJ + muJ)) * nG
    
    if N_offspring < 1
        return(0)
    else 
        # Distribution of mosquitoes across states at equilibrium
        B_prefactor = rhoJ * KJ *(N_offspring - 1) * nG / N_offspring
        temp_inv = inv(mu * I - transpose(A_mat))
        B_postfactor = temp_inv * alpha_vec_four

        B_star = B_prefactor * B_postfactor
        # V_star = r / (mu + gV)
        KB = sum(B_star)

        betaH_mat = zeros(Float64, 4,4); betaV_mat = zeros(Float64, 4,4); LambdaH = zeros(Float64, 4,4); LambdaV = zeros(Float64, 4,4)
        # Rate of contact is rate of entrance into transmission compartments
        LambdaH[3,3] = lP * KB / KH
        LambdaV[4,4] = lG
        betaH_mat[3,3] = bH
        betaV_mat[4,4] = bB

        spec_mat = alpha_vec_four * transpose(-A_mat * one_vec_four)

        GammaI = inv(mu * I - transpose(A_mat) + (gV / (mu + gV)) * (gR / (mu + gR)) * spec_mat)
        GammaE = inv((eta + mu) * I - transpose(A_mat) + (gR / (mu + gR + eta)) * (gV / (mu + gV + eta)) * spec_mat)

        complicated_probability = (gV/(mu+gV)) * ((eta / (mu+gV+eta)) * (gR/(mu+gR+eta)) + (eta/(mu+gR+eta) * (gR/(mu+gR))))
        tauE = (eta * I + complicated_probability * spec_mat) * GammaE

        sum_B_star = sum(B_star)

        R02 = (1 / (gH + muH)) * transpose(one_vec_four) * betaH_mat * LambdaH * GammaI * tauE * betaV_mat * LambdaV * (B_star / sum_B_star)
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
temp = R0_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
