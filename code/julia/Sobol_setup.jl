# Set up Sobol analysis
include("GCD_R0_num_calc.jl")

using GlobalSensitivity
using QuasiMonteCarlo
using DataFrames, CSV
using CodecZlib
# Define parameter ranges

# Rates at baseline
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)


# Function to set up parameter bounds around the baseline values
function parameter_setup(baseline_vals, stretch_val) 
    ubs = (1 + stretch_val) .* baseline_vals
    lbs = max(0,(1 - stretch_val)) .* baseline_vals
    ubs[5:(end-1)] = [min(v,1) for v in ubs[5:(end-1)]]

    return lbs, ubs
end

# In terms of rates (for first four parameters)
const base_params_flighty = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.9, 1.0, 0.5, 0.5, 0.5, 50.0]
const base_params_persistent = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.66, 1.0, 0.7, 0.8, 0.9, 50.0]
# In terms of durations (for first four parameters)
const base_invparams_flighty = [480.0, 10.0, 5.0, 1.0, 1.0 - 0.9, 1.0, 0.5, 0.5, 0.5, 50.0]
const base_invparams_persistent = [480.0, 10.0, 5.0, 1.0, 1.0 - 0.66, 1.0, 0.7, 0.8, 0.9, 50.0]

# Create lower and upper bounds for the parameters
persistent_lbs, persistent_ubs = parameter_setup(base_params_persistent, 0.2)
persistent_invlbs_temp, persistent_invubs_temp = parameter_setup(base_invparams_persistent, 0.2)
persistent_invlbs = copy(persistent_invlbs_temp)
persistent_invlbs[1:4] = 1.0 ./ persistent_invubs_temp[1:4] # change from durations back to rates
persistent_invubs = copy(persistent_invubs_temp)
persistent_invubs[1:4] = 1.0 ./ persistent_invlbs_temp[1:4] # change from durations back to rates
flighty_lbs, flighty_ubs = parameter_setup(base_params_flighty, 0.2)
flighty_invlbs_temp, flighty_invubs_temp = parameter_setup(base_invparams_flighty, 0.2)
flighty_invlbs = copy(flighty_invlbs_temp)
flighty_invlbs[1:4] = 1.0 ./ flighty_invubs_temp[1:4] # change from durations back to rates
flighty_invubs = copy(flighty_invubs_temp)
flighty_invubs[1:4] = 1.0 ./ flighty_invlbs_temp[1:4] # change from durations back to rates


# Define the absolute lower and upper bounds for the parameters

# In terms of rates (for first four parameters)
# adjusted lbs: lQ = 1/(480.0), all rates = 1/(30.0), all_probs = 0.2
min_lbs = [1/(480.0), 1/(30.0), 1/(30.0), 1/(30.0), 0.01, 0.2, 0.2, 0.2, 0.2, 0.01]
max_ubs = [160/1440.0, 2.0, 2.0, 2.0, 0.99, 0.8, 0.8, 0.8, 0.8, 0.99]

# Create function to calculate basic offspring number and basic reproduction number at the same time
function output_func(B_vals_in)
    repnums_func(B_vals_in)
end

# Significance test for parameters using eFAST indices
using Statistics

using HypothesisTests
# For NR times, calculate first and total order  eFAST indices for each parameter and save them to obtain a distribution of values for each parameter
function new_GSA_func(lbs, ubs, sample_size, NR)
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
        # A,B = QuasiMonteCarlo.generate_design_matrices(sample_size,lbs,ubs,SobolSample())
        # Sobol_res = gsa(output_func, Sobol(nboot = 100, conf_level = 0.95), A, B)
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

test = new_GSA_func(min_lbs, max_ubs, 1100, 10)


function new_get_gsa_results(sample_size, NR)

    # Define parameter sets and their names
    param_sets = [
        (persistent_lbs, persistent_ubs, "persistent"),
        (flighty_lbs, flighty_ubs, "flighty"),
        (min_lbs, max_ubs, "max"),
    ]

    # Collect results in a vector
    results_list = DataFrame[]

    for (lbs, ubs, set_name) in ProgressBar(param_sets)
        res = new_GSA_func(lbs, ubs, sample_size, NR)
        res[!, :type] .= set_name
        push!(results_list, res)
    end

    # Concatenate all results into a single DataFrame
    all_results = vcat(results_list...)

    return all_results

end

test = new_get_gsa_results(1000, 10)

new_eFAST_results_small = new_get_gsa_results(11*1_000, 11)
# new_eFAST_results = new_get_gsa_results(11*10_000, 11*10)

# Return all the unique entries in the type column of new_eFAST_results
# unique_types = unique(new_eFAST_results.type)

# Save outputs to CSV
open(joinpath((pwd()), "data", "new_eFAST_test_small.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, new_eFAST_results_small)
    close(gzip_io)
end


open(joinpath((pwd()), "data", "new_eFAST_test.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, new_eFAST_results)
    close(gzip_io)
end