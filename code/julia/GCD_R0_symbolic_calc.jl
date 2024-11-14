# Define the functions to calculate key outputs: 
# - mean gonotrophic cycle duration (meanGCD)
# - basic offspring number (N_offspring)
# - basic reproduction number (R0)

using Symbolics
using LinearAlgebra
using Serialization

#### Create mean GCD functions ####

# Give variables appropriate names
@variables V_vars[1:3] B_vars[1:9] 
(gV, gR, mu) = V_vars
(sigma, lQ, lL, lP, lG, pQ, pL, pP, pG) = B_vars
f = 1 - sigma

# Define subintensity matrix
A_mat = [-lQ         lQ                  0     0
         f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0
         f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP pP*lP
         f*(1-pG)*lG (1-f)*(1-pG)*lG     0  -lG]

tA_inv = inv(transpose(-A_mat))

temp_inv = simplify(inv(mu * I - A_mat))

one_vec = [1 1 1 1]
alpha_vec = [1;0;0;0]

GCD_sym = simplify(one_vec * temp_inv * alpha_vec)
GCD_mu_zero_sym = simplify(one_vec * tA_inv * alpha_vec)

# Generate a Julia function from the symbolic expression
GCD_func = Symbolics.build_function(GCD_sym, V_vars, B_vars)
GCD_mu_zero_func = Symbolics.build_function(GCD_mu_zero_sym, V_vars, B_vars)

# Save the function to use later without recalculation

# Save the functions to a file
Serialization.serialize("GCD_func.jls", GCD_func)
Serialization.serialize("GCD_mu_zero_func.jls", GCD_mu_zero_func)

# Here's how to call it later:
# GCD_func = Serialization.deserialize("GCD_func.jls")

#### Create basic offspring and basic reproduction number functions ####

# Give variables appropriate names
@variables L_vars[1:4] V_vars[1:3] B_vars[1:9] Epi_vars[1:4] Host_vars[1:2]
(KJ, rhoJ, muJ, varPhi) = L_vars
(gV, gR, mu) = V_vars
(sigma, lQ, lL, lP, lG, pQ, pL, pP, pG) = B_vars
(bH, bB, eta, gH) = Epi_vars
(muH, KH) = Host_vars

f = 1 - sigma

# Define subintensity matrices
A_mat = [-lQ         lQ                  0f0   0f0   
         f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0f0 
         f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP   pP*lP
         f*(1-pG)*lG (1-f)*(1-pG)*lG     0f0   -lG
         ]

A_tilde = [-lQ       lQ                  0f0   0f0   0f0
         f*(1-pL)*lL -lL+(1-f)*(1-pL)*lL pL*lL 0f0   0f0
         f*(1-pP)*lP (1-f)*(1-pP)*lP     -lP   pP*lP 0f0
         f*(1-pG)*lG (1-f)*(1-pG)*lG     0f0   -lG   pG*lG
         0f0         0f0                 0f0   0f0   -gV
         ]

one_vec_five = [1f0 1f0 1f0 1f0 1f0]
alpha_vec_five = [1f0; 0f0; 0f0; 0f0; 0f0]
one_vec_four = [1f0 1f0 1f0 1f0]
alpha_vec_four = [1f0; 0f0; 0f0; 0f0]

# Pr(complete a gonotrophic cycle)
tau = -one_vec_four * transpose(A_mat) * inv(mu * I - A_mat) * alpha_vec_four
tau = tau[1]
rho = 1f0 - (gV / (mu + gV)) * (gR / (mu + gR)) * tau
nG = 1f0 / rho
# Basic offspring number
N_offspring = tau * (varPhi / (mu + gV)) * (rhoJ / ( rhoJ + muJ)) * nG

# Generate a Julia function from the symbolic expression
N_offspring_func = Symbolics.build_function(N_offspring, L_vars, V_vars, B_vars)

# Save the basic offspring number function to a file
Serialization.serialize("N_offspring_func.jls", N_offspring_func)

# Continue to calculate R0 --

# r
r = (N_offspring - 1) * KJ * ((rhoJ + muJ)/ varPhi) * (mu + gV)

# Distribution of mosquitoes across states at equilibrium
B_prefactor = ((N_offspring - 1f0) * KJ * rhoJ / N_offspring) * (1f0 + (gR / (mu + gR)) * (gV / (mu + gV)) * nG * tau) * (inv(mu * I - transpose(A_mat)))
B_prefactor = B_prefactor[1]
B_postfactor = B_prefactor * alpha_vec_four

B_star = (inv(mu * I - transpose(A_mat))) * B_postfactor
V_star = r / (mu + gV)

B_vec = [B_star; V_star]

zero_mat = zeros(Num, 5, 5)
betaH_mat = zero_mat; betaV_mat = zero_mat; LambdaH = zero_mat; LambdaV = zero_mat

# # Rate of contact is rate of entrance into transmission compartments
LambdaH[2,3] = A_mat[2,3]
LambdaV[3,4] = A_mat[3,4]
betaH_mat[3,3] = bH
betaV_mat[4,4] = bB

temp_mat = (-one_vec_five) * transpose(A_tilde)

spec_mat = alpha_vec_five * temp_mat
M1 = mu * I - transpose(A_tilde) - (gR / (mu + gR)) * spec_mat
M2 = eta * I + (eta / (mu + gR + eta)) * (1 / (mu + gR)) * spec_mat
M3 = (eta + mu) * I - transpose(A_tilde) - (gR / (mu + gR + eta)) * spec_mat
Q = inv(M1) * M2 * inv(M3)

sum_B_star = sum(B_star)

R02 = one_vec_five * betaH_mat * transpose(LambdaH) * Q * betaV_mat * transpose(LambdaV) * B_vec / (sum_B_star * (muH + gH))
R02 = R02[1]

R0 = sqrt(R02)

# Generate a Julia function from the symbolic expression
R0_func = Symbolics.build_function(R0, L_vars, V_vars, B_vars, Epi_vars)

# Save the basic offspring number function to a file
Serialization.serialize("R0_func.jls", R0_func)