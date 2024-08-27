using Printf
using Distributions
using EMpht
using StatsBase
using LinearAlgebra
using ProgressBars


# data = rand(Exponential(1/10), 1_000)  # Generate some data to fit 
# sample = EMpht.Sample(obs=data)        # Create an EMpht Sample object with this data
# ph = empht(sample, p=5)                # Fit the data using p=5 phases

# xGrid = range(0, 8, length=1_00)       # Create a grid to evaluate the density function over
# fitPDFs = pdf.(ph, xGrid)              # The probability density function of the fitted phase-type

# Helper function -----------------------------------------------------------------------
# Replace very small elements with zeros
function set_zeros(array_in)
    array_in[abs.(array_in) .< eps(eltype(array_in))] .= 0.0
    return(array_in)
end

# Set up parameters ---------------------------------------------------------------------
# Questing
pQ = 1
lQ = 1 / (8 * 60) # 8 hours = 480 minutes

# Landing
pL = 0.7
lL = 0.1 # 10 minutes

# Probing
pP = 0.8
lP = 0.2 # 5 minutes

# Ingesting
pG = 0.9 #0.75
lG = 1 # 1 minutes

# Fleeing
f = 0.66

π = [1.0 0.0 0.0 0.0]

T = [-lQ               lQ                           0.0     0.0
     f * (1- pL) * lL  -lL + (1 - f) * (1- pL) * lL pL * lL 0.0
     f * (1 - pP) * lP (1 - f) * (1 - pP) * lP      -lP     pP * lP
     f * (1 - pG) * lG (1 - f) * (1 - pG) * lG      0.0     -lG
     ]

one_vec = transpose([1 1 1 1])

t = -T * one_vec

function PH_dist_sampler(n_samps::Int64, subint_mat::Matrix, alpha_vec)
    defect = 1 - sum(alpha_vec)
    init_probs = [alpha_vec defect]
    p = size(subint_mat)[1]
    
    out_rates = -subint_mat * ones(p)
    out_rates[abs.(out_rates) .< (10 * eps(eltype(out_rates)))] .= 0.0
    int_mat = hcat(subint_mat, out_rates)
    # Set up discrete Distributions
    values = 1:(p+1)
    probabilities = vec(init_probs)
    d = Categorical(probabilities)

    n_vec = nothing
    n_vec = zeros(n_samps)
    for i in 1:n_samps
        j = values[rand(d)]
        while j != (p + 1)
            n_vec[i] = n_vec[i] .+ rand(Exponential(-1/int_mat[j, j]), 1)[1]
            new_probs = int_mat[j,1:end .!=j]
            new_probs = new_probs ./ sum(new_probs)
            new_values = setdiff(1:(p+1), j)
            j = new_values[rand(Categorical(new_probs))]
        end
    end
    return(n_vec)
end

# Generate samples ----------------------------------------------------------------------
num_samples = 100_000::Int64
test_data = PH_dist_sampler(num_samples, T, π)

using Plots
histogram(log.(test_data))


# Intermediate fitting functions -------------------------------------------------------

function create_J_matrix(data_point::Float64, subint_mat::Matrix, alpha_vec::Matrix, t_vec::Matrix, PH_dim::Int64)
    temp_subint = zeros(2 * PH_dim, 2 * PH_dim)
    # temp_subint[(PH_dim+1):(2 * PH_dim), 1:PH_dim] .= 0.0
    temp_subint[1:PH_dim, 1:PH_dim] = subint_mat
    temp_subint[(PH_dim+1):(2 * PH_dim), (PH_dim+1):(2 * PH_dim)] = subint_mat
    temp_subint[1:PH_dim, (PH_dim+1):(2 * PH_dim)] = t_vec * alpha_vec
    # temp_subint[abs.(temp_subint) .< eps(eltype(temp_subint))] .= 0.0

    temp_matexp = exp(temp_subint * data_point)
    J_k = temp_matexp[1:PH_dim,(PH_dim+1):(2 * PH_dim)]
    return(J_k)
end

# test = create_J_matrix(10000.0::Float64, T, π, t, 4::Int64)

# Expected numbers for state i and data point k
function get_EB(idx, data_point, subint_mat, alpha_vec, t_vec)
    PH_dim = size(subint_mat)[1]
    e_i_vec = (1:PH_dim) .== idx
    common_denom = (alpha_vec * exp(subint_mat * data_point) * t_vec)[1]
    if common_denom == 0.0
        EB_ik = 0.0
        EZ_ik = 0.0
        EN_ij_k = zeros(1, PH_dim)
        EN_ik = 0.0
    else
        J_k = create_J_matrix(data_point::Float64, subint_mat::Matrix, alpha_vec::Matrix, t_vec::Matrix, PH_dim::Int64)

        # E[B_i]: Expected number of processes starting in state i
        EB_ik = ((alpha_vec[idx] * e_i_vec' * exp(subint_mat * data_point) * t_vec) / common_denom)[1]

        # E[Z_i]: Expected total time spent in state i in all processes
        EZ_ik = J_k[idx, idx] / common_denom

        # E[N_ij]: Expected number of observed jumps from state i to state j
        EN_ij_k = zeros(1, PH_dim)
        for j in 1:PH_dim
            EN_ij_k[1, j] = subint_mat[idx, j] * J_k[j, idx] / common_denom
        end
        # EN_ij_k = set_zeros(EN_ij_k)

        # E[N_i]: Expected number of processes that exit to the absorbing state from state i
        EN_ik = (alpha_vec * exp(subint_mat * data_point) * e_i_vec * t_vec[idx] / common_denom)[1]
    end

    return(EB_ik, EZ_ik, EN_ij_k, EN_ik)
end


function inner_expect(row, y_vec, subint_mat, alpha_vec, t_vec)
    PH_dim = size(subint_mat)[1]
    tempEB = 0.0
    tempEZ = 0.0
    tempEN = 0.0
    tempENij = zeros(1, PH_dim)
    for data_point in y_vec
        EB_ik, EZ_ik, EN_ij_k, EN_ik = get_EB(row, data_point, subint_mat, alpha_vec, t_vec)
        tempEB += EB_ik
        tempEZ += EZ_ik
        tempEN += EN_ik
        tempENij += EN_ij_k        
    end
    # tempENij
    return(tempEB, tempEZ, tempEN, tempENij)
end

# kk = 1

# data_point = y_vec[kk]
# EB_ik, EZ_ik, EN_ij_k, EN_ik = get_EB(row, data_point, subint_mat, alpha_vec, t_vec)
# tempEB += EB_ik
# tempEZ += EZ_ik
# tempEN += EN_ik
# tempENij += EN_ij_k
# println(kk)
# println(data_point)
# println(tempENij)

# kk = kk + 1

# Sum up all the expected numbers over the data set
function PH_expectation_function(data_in, subint_mat, alpha_vec)
    # Intermediate quantities
    A_mat = subint_mat
    PH_dim = size(A_mat)[1]
    e_vec = ones(PH_dim, 1)
    t_vec = - A_mat * e_vec
    # t_vec = set_zeros(t_vec)
    y_vec = data_in

    EB = nothing
    EZ = nothing
    EN = nothing
    ENij = nothing
    # EB = Matrix{Float64}(undef, 1, PH_dim)
    EB = zeros(1, PH_dim)
    # EZ = Matrix{Float64}(undef, 1, PH_dim)
    EZ = zeros(1, PH_dim)
    # EN = Matrix{Float64}(undef, 1, PH_dim)
    EN = zeros(1, PH_dim)
    # ENij = Matrix{Float64}(undef, PH_dim, PH_dim)
    ENij = zeros(PH_dim, PH_dim)

    for row in 1:PH_dim
        tempEB, tempEZ, tempEN, tempENij = inner_expect(row, y_vec, subint_mat, alpha_vec, t_vec)
        EB[row] = (tempEB)
        EZ[row] = (tempEZ)
        EN[row] = (tempEN)
        ENij[row, :] = (tempENij)
    end
    # ENij = set_zeros(ENij)
    ENij = convert(Matrix, ENij)
    return(EB, EZ, EN, ENij)
end

function PH_update_function(data_in, EB, EZ, EN, ENij)
    PH_dim = length(EB)
    data_length = length(data_in)

    π_new = EB / data_length

    subint_new = ENij * inv(Diagonal(vec(EZ)))
    t_vec_new = vec(EN * inv(Diagonal(vec(EZ))))

    subint_new[diagind(subint_new)] .= 0.0

    t_ii = -subint_new * ones(PH_dim) - t_vec_new
    subint_new[diagind(subint_new)] .= t_ii

    return(subint_new, π_new)
end

function full_fitting_function(number_iterations, data_in, subint_guess, alpha_guess)
    T_new = subint_guess
    π_new = alpha_guess

    prev_logLik = -Inf
    iter = ProgressBar(1:number_iterations)
    for i in iter
        # old_T = T_new
        # old_π = π_new
        # prev_EB, prev_EZ, prev_EN, prev_ENij = PH_expectation_function(data_in, old_T, old_π)
        # prev_ENij
        EB, EZ, EN, ENij = PH_expectation_function(data_in, T_new, π_new)

        T_new, π_new = PH_update_function(data_in, EB, EZ, EN, ENij)
        logLik_new = likelihood_function(data_in, T_new, π_new)

        logLik_diff = prev_logLik - logLik_new
        set_description(iter, string(@sprintf("logLik: %.2f", logLik_new)))
        prev_logLik = logLik_new
        # println(logLik)
    end
    # π_new = π_new / sum(π_new)
    return(T_new, π_new)
end

function likelihood_function(data_in, subint_mat, alpha_vec)
    out_rates = -subint_mat * ones(size(subint_mat)[1])
    Lik = 0.0
    for i in eachindex(1:length(data_in))
        Lik += log((alpha_vec * exp(subint_mat * data_in[i]) * out_rates)[1])
    end
    return(Lik)
end



# EB, EZ, EN, ENij = PH_expectation_function(data_in, T, π)
# subint_new, π_new = PH_update_function(data_in, EB, EZ, EN, ENij)

dim_guess = 4

subint_guess = rand(Float64, (dim_guess, dim_guess)) * 1e-2
# subint_guess[1,3] = 0.0
# subint_guess[1,4] = 0.0
# subint_guess[2,4] = 0.0
# subint_guess[4,3] = 0.0
subint_guess[diagind(subint_guess)] .= 0.0
t_vec_guess = vec(rand(Float64, (1, dim_guess))) * 1e-1

t_ii = -subint_guess * ones(dim_guess) - t_vec_guess
subint_guess[diagind(subint_guess)] .= t_ii

alpha_guess = rand(Float64, (1, dim_guess))
# alpha_guess = alpha_guess / sum(alpha_guess)
alpha_guess = [1.0 0.0 0.0 0.0]

data_in = test_data[1:100]

# pre-run the function
full_fitting_function(1, data_in, subint_guess, alpha_guess)
T_new, π_new = full_fitting_function(1_000, data_in, subint_guess, alpha_guess)

# Compare distributions
fit_data = PH_dist_sampler(num_samples, T_new, π_new)

histogram(log.(test_data.+1.0))
histogram!(log.(fit_data.+1.0))