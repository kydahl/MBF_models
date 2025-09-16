# Load in necessary output functions and constant parameters
include("Functions.jl")

# Define the flighty and persistent parameter sets [Note: the manuscript only uses the ``max'' parameter set]
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
const base_params_flighty = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.9, 1.0, 0.5, 0.5, 0.5, 50.0]
const base_params_persistent = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.66, 1.0, 0.7, 0.8, 0.9, 50.0]

function parameter_setup(baseline_vals, stretch_val) 
    ubs = (1f0 + stretch_val) .* baseline_vals
    lbs = max(0,(1f0 - stretch_val)) .* baseline_vals
    ubs[5:(end-1)] = [min(v,1) for v in ubs[5:(end-1)]]

    return lbs, ubs
end

# Flighty and persistent parameter ranges (lower and upper bounds)
flighty_lbs, flighty_ubs = parameter_setup(base_params_flighty, 0.2)
persistent_lbs, persistent_ubs = parameter_setup(base_params_persistent, 0.2)
# Maximum variation parameter ranges
min_lbs = [1/(480.0), 1/(30.0), 1/(30.0), 1/(30.0), 0.01, 0.2, 0.2, 0.2, 0.2, 0.01]
max_ubs = [160/1440.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# PRCC analysis ----

# Set up LHC sampling
using LatinHypercubeSampling

# Set number of LHC samples
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

# Save outputs
# Parameter grid and outputs
using CodecZlib
open(joinpath(dirname(pwd()), "data", "julia_PRCC.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, all_PRCCs)
    close(gzip_io)
end

# eFAST analysis ----

# Create function to calculate basic offspring number and basic reproduction number at the same time
function output_func(B_vals_in)
    repnums_func(B_vals_in)
end

# Significance test for parameters using eFAST indices
using Statistics
using GlobalSensitivity
using HypothesisTests
import DataFrames: groupby
# For NR times, calculate first and total order  eFAST indices for each parameter and save them to obtain a distribution of values for each parameter
function eFAST_func(lbs, ubs, sample_size, NR)
    # Give names for parameters, outputs, and index types
    param_names = ["lQ", "lL", "lP", "lG", "sigma", "pQ", "pL", "pP", "pG", "dummy"]
    output_names = ["N_offspring", "RVH", "RHV", "R0"]
    index_types = ["ST", "S1"]

    # Initialize DataFrame to hold results
    # Dataframe has columns for each parameter, each output, total and first order indices,
    eFAST_results = DataFrame(
        input = String[],
        output = String[],
        index_type = String[],
        value = Float64[],
        resample_num = Int[]
    )

    for resample_index in 1:NR # resample NR times
        # Calculate eFAST indices for the current resample
        eFAST_res = gsa(output_func, eFAST(num_harmonics = 11), [[lbs[i], ubs[i]] for i in 1:length(lbs)]; samples = sample_size)
        # Append results to the DataFrame
        for (i_out, outname) in enumerate(output_names)
            for (idx_type, arr) in zip(index_types, (eFAST_res.ST, eFAST_res.S1))
                for (i_param, pname) in enumerate(param_names)
                    push!(eFAST_results, (
                        input = pname,
                        output = outname,
                        index_type = idx_type,
                        value = arr[i_out, i_param],
                        resample_num = resample_index
                    ))
                end
            end
        end
    end

    # Calculate the mean and standard deviations of the first order (S1) and total order (ST) indices for each parameter and output
    mean_results = combine(
        groupby(eFAST_results, [:input, :output, :index_type]),
        :value => mean => :mean_value,
        :value => std => :std_value
    )

    # Determine the significance of each index value by conducting a two-sample t-test against the "dummy" parameter input
    # For each parameter (except "dummy"), output, and index_type, compare its distribution to "dummy"
    sig_results = DataFrame(
        input = String[],
        output = String[],
        index_type = String[],
        p_value = Float64[]
    )

    for pname in param_names[1:end-1] # exclude "dummy"
        for outname in output_names
            for idx_type in index_types
                vals_param = eFAST_results[(eFAST_results.input .== pname) .&& (eFAST_results.output .== outname) .&& (eFAST_results.index_type .== idx_type), :value]
                vals_dummy = eFAST_results[(eFAST_results.input .== "dummy") .&& (eFAST_results.output .== outname) .&& (eFAST_results.index_type .== idx_type), :value]
                # Test whether the parameter's values are significantly greater than the dummy parameter's values
                greater_sig_test = mannwhitney_onesided(vals_param, vals_dummy, alternative = :greater)

                push!(sig_results, (
                    input = pname,
                    output = outname,
                    index_type = idx_type,
                    p_value = greater_sig_test[2]
                ))
            end
        end
    end
    # Add the p_value results to the mean_results DataFrame
    mean_results = leftjoin(mean_results, sig_results, on=[:input, :output, :index_type])

    return mean_results
end

using Distributions

function mannwhitney_onesided(x::AbstractVector, y::AbstractVector; alternative::Symbol = :greater)
    test = MannWhitneyUTest(x, y)
    U = test.U
    nx, ny = length(x), length(y)

    # Exact null distribution if small, normal approx if large
    μ_U = nx * ny / 2
    σ_U = sqrt(nx * ny * (nx + ny + 1) / 12)

    if alternative == :greater
        # P(X > Y) means large U
        z = (U - μ_U - 0.5) / σ_U   # continuity correction
        p = 1 - cdf(Normal(), z)
    elseif alternative == :less
        # P(X < Y) means small U
        z = (U - μ_U + 0.5) / σ_U
        p = cdf(Normal(), z)
    else
        error("alternative must be :greater or :less")
    end

    return (U = U, pvalue = p, zscore = z)
end

test = eFAST_func(min_lbs, max_ubs, 1100, 10)


function get_gsa_results(sample_size, NR)

    # Define parameter sets and their names
    param_sets = [
        (persistent_lbs, persistent_ubs, "persistent"),
        (flighty_lbs, flighty_ubs, "flighty"),
        (min_lbs, max_ubs, "max"),
    ]

    # Collect results in a vector
    results_list = DataFrame[]

    for (lbs, ubs, set_name) in ProgressBar(param_sets)
        res = eFAST_func(lbs, ubs, sample_size, NR)
        res[!, :type] .= set_name
        push!(results_list, res)
    end

    # Concatenate all results into a single DataFrame
    all_results = vcat(results_list...)

    return all_results

end

get_gsa_results(1000, 10)

eFAST_results = get_gsa_results(11*10_000, 11*10)

# Save outputs to CSV
open(joinpath((pwd()), "data", "julia_eFAST.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, eFAST_results)
    close(gzip_io)
end
