using Symbolics
using LinearAlgebra

@variables sigma mu lQ lL lP lG pQ pL pP pG KH KB

# Function to substitute mu = 0 into expressions
function mu_zero_func(expr)
    new_expr = substitute(expr, Dict([mu => 0]))
    return(new_expr)
end

function simplify_matrix(mat)
    new_mat = mat
    length_mat = length(mat)
    for index in 1:length_mat
        new_mat[index] = simplify(mat[index]; expand = true)        
    end
    return(new_mat)
end

# Note to self: this level of detail is not necessary. 
# We can just use the mean blood-feeding stage duration neglecting mortality.
# i.e. mu = 0
# For the parameter values I've been using (mu = 1/(20 days)), including mu amounts to ~60 minute reduction in the duration when it is set to 1 day


A_mat = [-lQ+(1-pQ)*lQ       pQ*lQ               0     0
         (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0
         (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
         (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0     -lG]


# @variables b
# F_mat = [-b   0   0   0   0   0
#           0 -2b  2b   0   0   0
#           0   0 -2b   0   0   0
#           0   0   0 -3b  3b   0
#           0   0   0   0 -3b  3b
#           0   0   0   0   0 -3b]
# one_six = transpose([1 1 1 1 1 1])
# alpha_six = transpose([1/3 1/3 0 1/3 0 0])
# F_inv = simplify(inv(-F_mat))

# theta_F = transpose(alpha_six) * F_inv * one_six
# theta_F = theta_F[]
# simplify(theta_F)

# @variables lambda sigma1 sigma2
# D_mat = [-lambda (1-sigma1)*lambda 0
#           0      -lambda           (1-sigma2)*lambda
#           0      0                 -lambda]

# one_three = transpose([1 1 1])
# alpha_three = transpose([1 0 0])
# D_inv = simplify(inv(-D_mat))

# theta_D = transpose(alpha_three) * D_inv * one_three
# theta_D = theta_D[]
# simplify(theta_D)

# mu_theta_D = transpose(alpha_three) * inv(mu * I - D_mat) * one_three
# mu_theta_D = mu_theta_D[]
# simplify(mu_theta_D)
          
# rowsums_D = -transpose(one_three) * transpose(D_mat)

# tau_D = rowsums_D * inv(mu * I - D_mat) * alpha_three
# tau_D = tau_D[]
# simplify(tau_D)

# dist_D = transpose(alpha_three) * inv(mu * I - D_mat)


# @variables b1 b2 b3
# sigma_2 = 1 - ( (1 / (b2 * lambda)) * (1 - (mu + lambda)^2 * b1) * (lambda / (mu+lambda)) - (mu+lambda))
# sigma_1 = 1 - (( (mu + lambda) + (1 - sigma_2))^(-1)) * (1/lambda) * ((1/b1) - (mu+lambda)^2)

# eq_1 = b1 - (1 + (1 - sigma1) * (lambda / (mu+lambda)) + (1-sigma1) * (1-sigma2) * (lambda / (mu+lambda))^2)^(-1)
# test_1 = simplify(substitute(eq_1, Dict([sigma1 => sigma_1, sigma2 => sigma_2])), expand = true)

# # mat_dynamic_lQ = [mu+lQ*KH/KB -lQ*KH/KB 0 0
# #        -f*(1-pL)*lL mu+lL-(1-f)*(1-pL)*lL -pL*lL 0
# #        -f*(1-pP)*lP -(1-f)*(1-pP)*lP mu+lP -pP*lP
# #        -f*(1-pG)*lG -(1-f)*(1-pG)*lG 0 mu+lG]

neg_mat = -A_mat
# neg_mat = [-(mu+lQ) lQ 0 0
# f*(1-pL)*lL -(mu+lL)+(1-f)*(1-pL)*lL pL*lL 0
# f*(1-pP)*lP (1-f)*(1-pP)*lP -(mu+lP) pP*lP
# f*(1-pG)*lG (1-f)*(1-pG)*lG 0 -(mu+lG)]

temp_inv = simplify(inv(mu * I + neg_mat))

alpha = transpose([1 0 0 0])
one_vec = transpose([1 1 1 1])

row_sums = A_mat * one_vec
row_sums = simplify_matrix(row_sums)

@variables gV pR pE mu eta betaH[1:4] betaB[1:4] lambdaH[1:4] lambdaB[1:4]

tilde_A = hcat(A_mat[:,1], A_mat[:,2], A_mat[:,3], A_mat[:,4], -row_sums)
tilde_A = vcat(tilde_A, [0 0 0 0 -gV])

tilde_alpha = transpose([1 0 0 0 0])
tilde_ones = transpose([1 1 1 1 1])
tilde_row_sums = simplify_matrix(tilde_A * tilde_ones)

gammaI = simplify_matrix(mu * I - transpose(A_mat) - pR * alpha * transpose(row_sums))
GammaI = simplify_matrix(inv(gammaI))
tilde_gammaI = simplify_matrix(mu * I - transpose(tilde_A) - pR * tilde_alpha * transpose(tilde_row_sums))
tilde_GammaI = simplify_matrix(inv(tilde_gammaI))

gammaE = ((eta + mu) * I - transpose(A_mat) - pE * alpha * transpose(row_sums))
GammaE = (inv(gammaE))
tilde_gammaE = ((eta + mu) * I - transpose(tilde_A) - pE * tilde_alpha * transpose(tilde_row_sums))
tilde_GammaE = (inv(tilde_gammaE))

tauE_first = simplify_matrix( eta * I + pE * pR * alpha * transpose(row_sums))
tauE = tauE_first * GammaE
tilde_tauE_first = simplify_matrix( eta * I + pE * pR * tilde_alpha * transpose(tilde_row_sums))
tilde_tauE = tilde_tauE_first * tilde_GammaE


BetaH = diagm([betaH[1], betaH[2], betaH[3], betaH[4]])
BetaB = diagm([betaB[1], betaB[2], betaB[3], betaB[4]])
LambdaH = diagm([lambdaH[1], lambdaH[2], lambdaH[3], lambdaH[4]])
LambdaB = diagm([lambdaB[1], lambdaB[2], lambdaB[3], lambdaB[4]])
tilde_BetaH = diagm([betaH[1], betaH[2], betaH[3], betaH[4], 0])
tilde_BetaB = diagm([betaB[1], betaB[2], betaB[3], betaB[4], 0])
tilde_LambdaH = diagm([lambdaH[1], lambdaH[2], lambdaH[3], lambdaH[4], 0])
tilde_LambdaB = diagm([lambdaB[1], lambdaB[2], lambdaB[3], lambdaB[4], 0])

@variables C_scalar B[1:5]

B_prime = [B[1], B[2], B[3], B[4]]
tilde_B_prime = [B[1], B[2], B[3], B[4], B[5]]

R0 = (transpose(one_vec) * BetaH * LambdaH * GammaI * tauE * BetaB * LambdaB * B_prime)[]
tilde_R0 = (transpose(tilde_ones) * tilde_BetaH * tilde_LambdaH * tilde_GammaI * tilde_tauE * tilde_BetaB * tilde_LambdaB * tilde_B_prime)[]

R0_test = simplify(R0 - tilde_R0)

diff_test = simplify(R0 - R0_test)

# # Alternate calculation
# pre_vec = transpose(alpha) * temp_inv
# simple_pre_vec = pre_vec
# simple_pre_vec[1] = simplify(pre_vec[1])
# simple_pre_vec[2] = simplify(pre_vec[2])
# simple_pre_vec[3] = simplify(pre_vec[3])
# simple_pre_vec[4] = simplify(pre_vec[4])
# mean_dur_alt = simple_pre_vec * one_vec

# Total mean duration with mortality (mu > 0)
mean_dur = transpose(one_vec) * temp_inv * alpha
mean_dur = mean_dur[]
# Total mean duration with NO mortality (mu = 0 )
simple_mean_dur = simplify(substitute(mean_dur, Dict([mu => 0])))


theta = transpose(alpha) * inv(neg_mat) * one_vec
theta = theta[]
simple_theta = simplify(theta)

# # Try to calculate mean duration (theta) as a function of mosquito population density when lQ = lQ * KH / KB
# D = Differential(KB)
# KB_derive = simplify(expand_derivatives(D(mean_dur[])))

# linear_KB_check = simplify(substitute(mean_dur[], Dict([KB => 0])))

# Calculate aspects of R0

# New parameters related to transmission and oviposition
@variables eta gammaV gammaR Q_state L_state P_state G_state V_state

B_state_vec = [Q_state; L_state; P_state; G_state; V_state]

betaH_mat = [0 0 0 0 0
             0 0 0 0 0
             0 0 1 0 0
             0 0 0 0 0
             0 0 0 0 0]

betaV_mat = [0 0 0 0 0
             0 0 0 0 0
             0 0 0 0 0
             0 0 0 1 0
             0 0 0 0 0]
# Update the matrix as needed
out_rates = -transpose(one_vec) * transpose(A_mat)
out_rates[1] = expand(out_rates[1])
out_rates[2] = expand(out_rates[2])
out_rates[3] = expand(out_rates[3])
out_rates[4] = expand(out_rates[4])

A_tilde = [A_mat transpose(out_rates)]
A_tilde = [A_tilde; [0 0 0 0 -gammaV]]

ext_alpha = [1;0;0;0;0]
ext_one = [1;1;1;1;1]

temp_mat = simplify_matrix(ext_alpha * transpose(ext_one) * transpose(A_tilde))

full_dur_mu_zero = simplify_matrix(inv(- transpose(A_tilde) + temp_mat))
full_dur_mat = simplify_matrix(inv(mu * I - transpose(A_tilde) + (gammaR / (mu + gammaR)) * temp_mat))

dev_rate_mat = simplify_matrix(eta * I + (eta / (mu + gammaR + eta)) * (gammaR / (mu + gammaR)) * temp_mat)

dev_dur_mat = simplify_matrix(inv((mu + eta) * I - transpose(A_tilde) + (gammaR / (mu + gammaR)) * temp_mat))

dev_prob_mat = simplify_matrix(dev_rate_mat * dev_dur_mat)

mid_mat = simplify_matrix(full_dur_mat * dev_prob_mat)

part_R0_mat = transpose(ext_one) * betaH_mat * mid_mat * betaV_mat * B_state_vec

# Attempt to simplify this expression
# attempt_1 = expand(part_R0_mat)

# attempt_2 = simplify(part_R0_mat)


simple_attempt_1 = mu_zero_func(part_R0_mat)
simple_attempt_2 = simplify(substitute(part_R0_mat, Dict([pG => 0])))