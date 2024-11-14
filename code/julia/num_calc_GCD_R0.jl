# Numerically calculate key outputs across parameters ranges

# Load packages
using IntervalSets
using Symbolics
using Serialization

# Load in output functions
GCD_func = Serialization.deserialize("GCD_func.jls");
N_offspring_func = Serialization.deserialize("N_offspring_func.jls");
R0_func = Serialization.deserialize("R0_func.jls");

# Define fixed parameters
(KJ, rhoJ, muJ, varPhi) = (3E8, (1/(12 * 1440)), (1 / (20 * 1440)), (300 / 1440))
(gV, gR, mu) = ((1/(1 * 1440)), (1/(1 * 1440)), (1/(20 * 1440)))
(bH, bB, eta, gH) = (1, 1, (1/(6 * 1440)), (1/(7 * 1440)))
(muH, KH) = (1/(365.25 * 65 * 1440), 1E5)

# Combine them for easy reference
L_vals = (KJ, rhoJ, muJ, varPhi)
V_vals = (gV, gR, mu)
Epi_vals = (bH, bB, eta, gH)
Host_vals = (muH, KH)


# Define parameter ranges

# Set resolution of variation
resolution = 10

# Probabilities
generic_prob_vec = range(0f0, 1f0, resolution)
sigma_vec = generic_prob_vec
pL_vec = generic_prob_vec
pP_vec = generic_prob_vec
pG_vec = generic_prob_vec

# Rates
# - These are defined by baseline values for durations, then inverted
function inverse_rate_range_function(baseline, resolution)
    # Get a large range of values defined from the inverse
    big_inv_vec = range(0f0, 10f0 * 1440, resolution+1)[2:end]
    # Get a smaller range of values centered around the baseline
    base_range_vec = range(0.5f0, 100f0, resolution) * baseline
    # Combine them 
    inv_vec = sort(unique([big_inv_vec; 1f0 ./ base_range_vec]), rev = true)

    return(1f0 ./ inv_vec)

end

sigma_baseline = 0.9
lQ_baseline = 1 / 480
lL_baseline = 0.1f0
lP_baseline = 0.2f0
lG_baseline = 1f0

lQ_vec = inverse_rate_range_function(lQ_baseline, resolution)
lL_vec = inverse_rate_range_function(lL_baseline, resolution)
lP_vec = inverse_rate_range_function(lP_baseline, resolution)
lG_vec = inverse_rate_range_function(lG_baseline, resolution)

# (sigma, lQ, lL, lP, lG, pQ, pL, pP, pG) = B_vars

# Set up grid of parameter combinations
using Base.Threads
using IterTools
parameter_grid = collect(IterTools.product(sigma_vec, pL_vec, pP_vec, pG_vec, lQ_vec, lL_vec, lP_vec, lG_vec))

# Calculate key outputs over the parameter combinations


# Prepare to store results
results = Vector{Float32}(undef, resolution^8::Int)

using ProgressBars
# Evaluate the function across all input combinations in parallel
for (i, (sigma, pL, pP, pG, lQ, lL, lP, lG)) in ProgressBar(enumerate(parameter_grid))

    B_vals = (sigma, pL, pP, pG, lQ, lL, lP, lG)

    results[i] = GCD_func(V_vals, B_vals)

end
