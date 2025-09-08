# Set up PRCC analysis
include("GCD_R0_num_calc.jl")

# Rates
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
# In terms of rates (for first four parameters)
const base_params_flighty = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.9, 1.0, 0.5, 0.5, 0.5, 50.0]
const base_params_persistent = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.66, 1.0, 0.7, 0.8, 0.9, 50.0]


function parameter_setup(baseline_vals, stretch_val) 
    ubs = (1f0 + stretch_val) .* baseline_vals
    lbs = max(0,(1f0 - stretch_val)) .* baseline_vals
    ubs[5:(end-1)] = [min(v,1) for v in ubs[5:(end-1)]]

    return lbs, ubs
end

# Create lower and upper bounds for the parameters
persistent_lbs, persistent_ubs = parameter_setup(base_params_persistent, 0.2)
flighty_lbs, flighty_ubs = parameter_setup(base_params_flighty, 0.2)

# Set up LHC sampling
using LatinHypercubeSampling

# Maximum variation parameter ranges
correction_term = 0.0
min_lbs = [1/(480.0), 1/(30.0), 1/(30.0), 1/(30.0), 0.01, 0.2, 0.2, 0.2, 0.2, 0.01]
max_ubs = [160/1440.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# Set number of LHC samples
# n_samples = 10_000_000::Int
n_samples = 100_000::Int

# Set up initial grid
using QuasiMonteCarlo
max_scaled_plan = QuasiMonteCarlo.sample(n_samples, min_lbs, max_ubs, LatinHypercubeSample())
flighty_scaled_plan = QuasiMonteCarlo.sample(n_samples, flighty_lbs, flighty_ubs, LatinHypercubeSample())
persistent_scaled_plan = QuasiMonteCarlo.sample(n_samples, persistent_lbs, persistent_ubs, LatinHypercubeSample())

# Set up grid of parameter combinations
using Base.Threads
using IterTools

function output_calc(LHS_samples)
    # Prepare to store results
    rows = Vector{NamedTuple}(undef, size(LHS_samples, 2))

    Threads.@threads for i in ProgressBar(1:size(LHS_samples, 2))
        params = LHS_samples[:,i]
        out = repnums_func(params[1:9])
        rows[i] = (
            lQ = params[1], lL = params[2], lP = params[3], lG = params[4],
            sigma = params[5], pQ = params[6], pL = params[7], pP = params[8], pG = params[9], dummy = params[10],
            N_offspring = out[1], RVH = out[2], RHV = out[3], R0 = out[4]
        )
    end
    return DataFrame(rows)
end

# Load output function
output_calc(QuasiMonteCarlo.sample(100, min_lbs, max_ubs, LatinHypercubeSample()))

# Get outputs for each type then join them
max_results = output_calc(max_scaled_plan)
flighty_results = output_calc(flighty_scaled_plan)
persistent_results = output_calc(persistent_scaled_plan)

# Calculate PRCC indices
using DataFrames
using Statistics
function PRCC_calc(input_results, output_names::Vector{String}, type::String)
    input_cols = names(input_results)[1:end-length(output_names)]
    results = DataFrame()
    for output_name in output_names
        rank_inputs = [sortperm(sortperm(input_results[!, col])) for col in input_cols]
        rank_output = sortperm(sortperm(input_results[!, output_name]))
        prccs = [cor(rank_inputs[i], rank_output) for i in 1:length(input_cols)]
        row = (; (Symbol(col) => prccs[i] for (i, col) in enumerate(input_cols))..., output = output_name, type = type)
        push!(results, row)
    end
    return results
end

output_names = ["N_offspring", "RVH", "RHV", "R0"]
max_PRCCs = PRCC_calc(max_results, output_names, "max")
flighty_PRCCs = PRCC_calc(flighty_results, output_names, "flighty")
persistent_PRCCs = PRCC_calc(persistent_results, output_names, "persistent")

all_PRCCs = vcat(max_PRCCs, flighty_PRCCs, persistent_PRCCs)

# Save outputs ----
# Parameter grid and outputs
using CodecZlib
open(joinpath(pwd(), "data", "julia_PRCC.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, all_PRCCs)
    close(gzip_io)
end
