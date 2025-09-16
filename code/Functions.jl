# Numerical calculators for key outputs
# Load packages
using IntervalSets
using LinearAlgebra
using Symbolics
using Serialization
using DataFrames
using CSV
using ProgressBars

# Define constants
# Immature mosquito parameters
const (KJ, rhoJ, muJ, varPhi) = [0.75 * 1E3, (1/(12 * 1440)), (1 / (20 * 1440)), (3 / 1440)]
# Adult mosquito parameters
const (gV, gR, mu) = [(1/(5 * 1440)), (1/(2 * 1440)), (1/(21 * 1440))]
# Host parameters
const (muH, KH) = [1/(365.25 * 65 * 1440), 1E3]
# Epidemiological parameters
const (bH, bB, eta, gH) = [0.5f0, 0.5f0, 1 / (7 * 1440), 1 / (5 * 1440)]

# Define functions to calculate key outputs
function GCD_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

    GCD = transpose(alpha_vec_four) * transpose(A_mat) * one_vec_four

    return(GCD)
end

function N_offspring_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

    tau = transpose(-A_mat * one_vec_four) * inv(mu * I - transpose(A_mat)) * alpha_vec_four
    tau = tau[1]
    rho = tau * (gV / (mu + gV)) * (gR / (mu + gR))
    nG = 1.0f0 / (1.0f0 - rho)
    # Basic offspring number
    N_offspring = tau * nG * (rhoJ / (rhoJ + muJ)) * (varPhi / (mu + gV))
    return(N_offspring)
end

function R0_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

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
        KB = sum(B_star)

        betaH_mat = zeros(Float64, 4,4)
        betaV_mat = zeros(Float64, 4,4)
        LambdaH = zeros(Float64, 4,4)
        LambdaV = zeros(Float64, 4,4)
        # Rate of contact is rate of entrance into transmission compartments
        LambdaH[3,3] = lP * KB / KH
        LambdaV[4,4] = lG
        betaH_mat[3,3] = bH
        betaV_mat[4,4] = bB

        spec_mat = alpha_vec_four * transpose(-A_mat * one_vec_four)

        GammaI = inv(mu * I - transpose(A_mat) + (gV / (mu + gV)) * (gR / (mu + gR)) * spec_mat)
        GammaE = inv((eta + mu) * I - transpose(A_mat) + (gR / (mu + gR + eta)) * (gV / (mu + gV + eta)) * spec_mat)

        complicated_probability = (gV/(mu+gV)) * ((1 - eta / (mu+gR+eta)) * (eta/(mu+gV+eta)) + (eta/(mu+gR+eta) * (gR/(mu+gR))))
        tauE = (eta * I + complicated_probability * spec_mat) * GammaE

        R02 = (KH / KB) * transpose(one_vec_four) * betaH_mat * LambdaH * GammaI * tauE * (1/KH) * (betaV_mat * LambdaV * B_star) * (1 / (gH + muH))
        R0 = sqrt(R02[1])
        return(R0)
    end
end

function repnums_func(B_vals_in)

    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in 

    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

    tau = transpose(-A_mat * one_vec_four) * inv(mu * I - transpose(A_mat)) * alpha_vec_four
    tau = tau[1]
    rho = (gV / (mu + gV)) * (gR / (mu + gR)) * tau
    nG = 1.0f0 / (1.0f0 - rho)

    # Basic offspring number
    N_offspring = tau * (varPhi / (mu + gV)) * (rhoJ / (rhoJ + muJ)) * nG

    if N_offspring < 1
            RHV_scalar = NaN
            RVH_scalar = NaN
            R0 = NaN
            R0_alt = NaN
    else
        # Distribution of mosquitoes across states at equilibrium
        B_prefactor = rhoJ * KJ *(N_offspring - 1) * nG / N_offspring
        temp_inv = inv(mu * I - transpose(A_mat))
        B_postfactor = temp_inv * alpha_vec_four

        B_star = B_prefactor * B_postfactor
        KB = sum(B_star)

        betaH_mat = zeros(Float64, 4,4)
        betaV_mat = zeros(Float64, 4,4)
        LambdaH = zeros(Float64, 4,4)
        LambdaV = zeros(Float64, 4,4)
        # Rate of contact is rate of entrance into transmission compartments
        LambdaH[3,3] = lP * KB / KH
        LambdaV[4,4] = lG
        betaH_mat[3,3] = bH
        betaV_mat[4,4] = bB

        spec_mat = alpha_vec_four * transpose(-A_mat * one_vec_four)

        GammaI = inv(mu * I - transpose(A_mat) - (gV / (mu + gV)) * (gR / (mu + gR)) * spec_mat)
        GammaE = inv((mu+eta)*I-transpose(A_mat) - (gR / (mu + gR + eta)) * (gV / (mu + gV + eta)) * spec_mat)

        complicated_probability = (gV/(mu+gV)) * ((1-(eta / (mu+gR+eta))) * (eta/(mu+gV+eta)) + (eta/(mu+gR+eta))) * (gR/(mu+gR))
        tauE = (eta * I + complicated_probability * spec_mat) * GammaE

        # Reproduction number matrices
        # RVH = FVH * big_inv
        # RHV = FHV * inv(VH)
        # # R02 = RVH * RHV

        # RVH_scalar = RVH[1,4]
        # RHV_scalar = RHV[4,1]
        # R0 =  sqrt(RVH_scalar * RHV_scalar)

        RHV_scalar = (1/KH) * (betaV_mat * LambdaV * B_star) * (1 / (gH + muH))
        R02_scalar = (KH / KB) * transpose(one_vec_four) * betaH_mat * LambdaH * GammaI * tauE * (1/KH) * (betaV_mat * LambdaV * B_star) * (1 / (gH + muH))
        R0 = sqrt(R02_scalar[1])
        RHV_scalar = RHV_scalar[4]
        RVH_scalar = R02_scalar[1] / RHV_scalar
    end

    return([N_offspring, RVH_scalar, RHV_scalar, R0])
end

# Run each function once to pre-load them into Julia
GCD_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
N_offspring_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
R0_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
repnums_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
