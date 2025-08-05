# Set up Sobol analysis
include("GCD_R0_num_calc.jl")

using GlobalSensitivity
using QuasiMonteCarlo
using DataFrames, CSV
using CodecZlib
# Define parameter ranges

# Rates at baseline
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
# In terms of rates (for first four parameters)
const base_params_flighty = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.9, 1.0, 0.5, 0.5, 0.5, 50.0]
const base_params_persistent = [1/480.0, 1/10.0, 1/5.0, 1/1.0, 1.0 - 0.66, 1.0, 0.7, 0.8, 0.9, 50.0]
# In terms of durations (for first four parameters)
const base_invparams_flighty = [480.0, 10.0, 5.0, 1.0, 1.0 - 0.9, 1.0, 0.5, 0.5, 0.5, 50.0]
const base_invparams_persistent = [480.0, 10.0, 5.0, 1.0, 1.0 - 0.66, 1.0, 0.7, 0.8, 0.9, 50.0]

# Function to set up parameter bounds around the baseline values
function parameter_setup(baseline_vals, stretch_val) 
    ubs = (1 + stretch_val) .* baseline_vals
    lbs = max(0,(1 - stretch_val)) .* baseline_vals
    ubs[5:(end-1)] = [min(v,1) for v in ubs[5:(end-1)]]

    return lbs, ubs
end

persistent_lbs, persistent_ubs = parameter_setup(base_params_persistent, 0.2)
persistent_invlbs, persistent_invubs = parameter_setup(base_invparams_persistent, 0.2)
flighty_lbs, flighty_ubs = parameter_setup(base_params_flighty, 0.2)
flighty_invlbs, flighty_invubs = parameter_setup(base_invparams_flighty, 0.2)


# Define the absolute lower and upper bounds for the parameters

# In terms of rates (for first four parameters)
min_lbs = [1/((1/2)*1440.0), 1/(30.0), 1/(30.0), 1/(30.0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
max_ubs = [160/1440.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0]
# In terms of durations (for first four parameters)
invmin_lbs = [1440/160.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
invmax_ubs = [(1/2)*1440.0, 30.0, 30.0, 30.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0]

# Create function to calculate basic offspring number and basic reproduction number at the same time
function output_func(B_vals_in)
    B_vals_in[1:4] = 1.0./B_vals_in[1:4] # change from durations back to rates
    [N_offspring_func(B_vals_in), R0_func(B_vals_in)]
end

function gsa_func(sample_size)
    # Set number of samples
    n_samples = sample_size #2^17::Int #10_000::Int

    # Set up QuasiMonteCarlo sampling
    # pad_val = round(Int, log2(n_samples)) + 4 # Pad to ensure enough samples
    # sampler = SobolSample(R = OwenScramble(base = 2, pad = 32))
    # sampler = SobolSample(R = OwenScramble(base = 2, pad = 32))

    # A_persistent, B_persistent = QuasiMonteCarlo.generate_design_matrices(n_samples, persistent_invlbs, persistent_invubs, sampler)
    # A_flighty, B_flighty = QuasiMonteCarlo.generate_design_matrices(n_samples, flighty_invlbs, flighty_invubs, sampler)
    # A_max, B_max = QuasiMonteCarlo.generate_design_matrices(n_samples, invmin_lbs, invmax_ubs, sampler)

    eFAST_persistent = gsa(output_func, eFAST(), [[persistent_invlbs[i], persistent_invubs[i]] for i in 1:length(persistent_invlbs)]; samples = n_samples)
    eFAST_flighty = gsa(output_func, eFAST(), [[flighty_invlbs[i], flighty_invubs[i]] for i in 1:length(flighty_invlbs)]; samples = n_samples)
    eFAST_max = gsa(output_func, eFAST(), [[invmin_lbs[i], invmax_ubs[i]] for i in 1:length(invmin_lbs)]; samples = n_samples)

    # # Get outputs for each type then join them
    # sobol_persistent = gsa(output_func, Sobol(), A_persistent, B_persistent)
    # sobol_flighty = gsa(output_func, Sobol(), A_flighty, B_flighty)
    # sobol_max = gsa(output_func, Sobol(), A_max, B_max)

    # Create dataframe to hold all all_results

    # Dataframe has columns for each parameter, each output, total and first order indices, and parameter set
    param_names = ["lQ", "lL", "lP", "lG", "sigma", "pQ", "pL", "pP", "pG", "dummy"]
    output_names = ["N_offspring", "R0"]
    index_types = ["ST", "S1"]
    paramset_types = ["persistent", "flighty", "max"]

    all_results = DataFrame(
        input = String[],
        output = String[],
        index_type = String[],
        value = Float64[],
        type = String[]
    )

    for (eFAST, paramset) in zip((eFAST_persistent, eFAST_flighty, eFAST_max), paramset_types)
        for (i_out, outname) in enumerate(output_names)
            for (idx_type, arr) in zip(index_types, (eFAST.ST, eFAST.S1))
                for (i_param, pname) in enumerate(param_names)
                    push!(all_results, (
                        input = pname,
                        output = outname,
                        index_type = idx_type,
                        value = arr[i_out, i_param],
                        type = paramset
                    ))
                end
            end
        end
    end
    return(all_results)
end

# Create dataframe to add results across sample sizes
samples_results = DataFrame(
        input = String[],
        output = String[],
        index_type = String[],
        value = Float64[],
        type = String[],
        sample_size = Int[]
    )

using Base.Threads

function get_gsa_results(sample_sizes)

    # sample_sizes = broadcast(^, 2, power_range) #10^4::Int : 5*10^3::Int : 10^5::Int
    results_list = Vector{DataFrame}(undef, length(sample_sizes))

    for idx in ProgressBar(eachindex(sample_sizes))
        i = sample_sizes[idx]
        results = gsa_func(i)
        results[!, :sample_size] .= i
        results_list[idx] = results
    end

    for results in results_list
        append!(samples_results, results)
    end
    # Save outputs to CSV
    open(joinpath((pwd()), "data", "eFAST_test.csv"), "w") do io
        gzip_io = GzipCompressorStream(io)
        CSV.write(gzip_io, samples_results)
        close(gzip_io)
    end
end

sample_sizes = 10^5::Int : 5*10^4::Int : 10^6::Int
get_gsa_results(10^2)
get_gsa_results(sample_sizes) #skipped 21 through 24
# # Save outputs ----
# # Parameter grid and outputs

# open(joinpath((pwd()), "data", "julia_eFAST.csv"), "w") do io
#     gzip_io = GzipCompressorStream(io)
#     CSV.write(gzip_io, all_results)
#     close(gzip_io)
# end