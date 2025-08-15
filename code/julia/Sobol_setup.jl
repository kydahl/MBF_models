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
correction_term = 0.0
min_lbs = [1/((1/2)*1440.0), 1/(30.0), 1/(30.0), 1/(30.0), 0.0+correction_term, 0.0+correction_term, 0.0+correction_term, 0.0+correction_term, 0.0+correction_term, 0.0]
max_ubs = [160/1440.0, 2.0, 2.0, 2.0, 1.0-correction_term, 1.0-correction_term, 1.0-correction_term, 1.0-correction_term, 1.0-correction_term, 100.0]

# Create function to calculate basic offspring number and basic reproduction number at the same time
function output_func(B_vals_in)
    [N_offspring_func(B_vals_in), R0_func(B_vals_in)]
end

# function gsa_func(sample_size)
#     # Set number of samples
#     n_samples = sample_size #2^17::Int #10_000::Int

#     eFAST_persistent = gsa(output_func, eFAST(), [[persistent_lbs[i], persistent_ubs[i]] for i in 1:length(persistent_lbs)]; samples = n_samples)
#     eFAST_flighty = gsa(output_func, eFAST(), [[flighty_lbs[i], flighty_ubs[i]] for i in 1:length(flighty_lbs)]; samples = n_samples)
#     eFAST_max = gsa(output_func, eFAST(), [[min_lbs[i], max_ubs[i]] for i in 1:length(min_lbs)]; samples = n_samples)

#     eFAST_persistent_inv = gsa(output_func, eFAST(), [[persistent_invlbs[i], persistent_invubs[i]] for i in 1:length(persistent_invlbs)]; samples = n_samples)
#     eFAST_flighty_inv = gsa(output_func, eFAST(), [[flighty_invlbs[i], flighty_invubs[i]] for i in 1:length(flighty_invlbs)]; samples = n_samples)
#     eFAST_max_inv = gsa(output_func, eFAST(), [[invmin_lbs[i], invmax_ubs[i]] for i in 1:length(invmin_lbs)]; samples = n_samples)

#     # Create dataframe to hold all all_results

#     # Dataframe has columns for each parameter, each output, total and first order indices, and parameter set
#     param_names = ["lQ", "lL", "lP", "lG", "sigma", "pQ", "pL", "pP", "pG", "dummy"]
#     output_names = ["N_offspring", "R0"]
#     index_types = ["ST", "S1"]
#     paramset_types = ["persistent", "flighty", "max", "inv_persistent", "inv_flighty", "inv_max"]

#     all_results = DataFrame(
#         input = String[],
#         output = String[],
#         index_type = String[],
#         value = Float64[],
#         type = String[]
#     )

#     for (eFAST, paramset) in zip((eFAST_persistent, eFAST_flighty, eFAST_max, eFAST_persistent_inv, eFAST_flighty_inv, eFAST_max_inv), paramset_types)
#         for (i_out, outname) in enumerate(output_names)
#             for (idx_type, arr) in zip(index_types, (eFAST.ST, eFAST.S1))
#                 for (i_param, pname) in enumerate(param_names)
#                     push!(all_results, (
#                         input = pname,
#                         output = outname,
#                         index_type = idx_type,
#                         value = arr[i_out, i_param],
#                         type = paramset
#                     ))
#                 end
#             end
#         end
#     end
#     return(all_results)
# end

# # Create dataframe to add results across sample sizes
# samples_results = DataFrame(
#     input = String[],
#     output = String[],
#     index_type = String[],
#     value = Float64[],
#     type = String[],
#     sample_size = Int[]
# )

# using Base.Threads

# function get_gsa_results(sample_sizes)

#     # sample_sizes = broadcast(^, 2, power_range) #10^4::Int : 5*10^3::Int : 10^5::Int
#     results_list = Vector{DataFrame}(undef, length(sample_sizes))

#     for idx in ProgressBar(eachindex(sample_sizes))
#         i = sample_sizes[idx]
#         results = gsa_func(i)
#         results[!, :sample_size] .= i
#         results_list[idx] = results
#     end

#     for results in results_list
#         append!(samples_results, results)
#     end
#     # Save outputs to CSV
#     open(joinpath((pwd()), "data", "eFAST_test.csv"), "w") do io
#         gzip_io = GzipCompressorStream(io)
#         CSV.write(gzip_io, samples_results)
#         close(gzip_io)
#     end
# end

# sample_sizes = 10^5::Int : 5*10^4::Int : 10^6::Int
# get_gsa_results(10^2)
# get_gsa_results(sample_sizes) 

# Significance test for parameters using eFAST indices
using Statistics

using HypothesisTests
# For NR times, calculate first and total order  eFAST indices for each parameter and save them to obtain a distribution of values for each parameter
function new_GSA_func(lbs, ubs, sample_size, NR)
    # Give names for parameters, outputs, and index types
    param_names = ["lQ", "lL", "lP", "lG", "sigma", "pQ", "pL", "pP", "pG", "dummy"]
    output_names = ["N_offspring", "R0"]
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
        (persistent_invlbs, persistent_invubs, "inv_persistent"),
        (flighty_invlbs, flighty_invubs, "inv_flighty")
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

new_eFAST_results = new_get_gsa_results(11*10_000, 110)

# Return all the unique entries in the type column of new_eFAST_results
unique_types = unique(new_eFAST_results.type)

# Save outputs to CSV
open(joinpath((pwd()), "data", "new_eFAST_test.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, new_eFAST_results)
    close(gzip_io)
end