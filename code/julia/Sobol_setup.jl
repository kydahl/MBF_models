# Set up Sobol analysis
include("GCD_R0_num_calc.jl")

using GlobalSensitivity
using QuasiMonteCarlo
using DataFrames, CSV
using CodecZlib
# Define parameter ranges

# Rates at baseline
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
const base_params_flighty = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.9f0, 1.0f0, 0.5f0, 0.5f0, 0.5f0]

const base_params_persistent = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.66f0, 1.0f0, 0.7f0, 0.8f0, 0.9f0]
    
# Function to set up parameter bounds around the baseline values
function parameter_setup(baseline_vals, stretch_val) 
    ubs = (1f0 + stretch_val) .* baseline_vals
    lbs = max(0,(1f0 - stretch_val)) .* baseline_vals
    ubs[5:end] = [min(v,1) for v in ubs[5:end]]

    return lbs, ubs
end

persistent_lbs, persistent_ubs = parameter_setup(base_params_persistent, 0.1)
flighty_lbs, flighty_ubs = parameter_setup(base_params_flighty, 0.1)

# Set up LHC sampling
using LatinHypercubeSampling

# Define the absolute lower and upper bounds for the parameters
min_lbs = [1/(1*1440.0f0), 1/(1*1440.0f0), 1/(1*1440.0f0), 1/(1*1440.0f0), 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0]
max_ubs = [1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0]

# Create function to calculate basic offspring number and basic reproduction number at the same time
function output_func(B_vals_in)
    [N_offspring_func(B_vals_in), R0_func(B_vals_in)]
end

function gsa_func(sample_size)
    # Set number of samples
    n_samples = sample_size #2^17::Int #10_000::Int

    # # Set up QuasiMonteCarlo sampling
    # sampler = SobolSample(R = OwenScramble(base = 2, pad = 21))

    # A_persistent, B_persistent = QuasiMonteCarlo.generate_design_matrices(n_samples, persistent_lbs, persistent_ubs, sampler)
    # A_flighty, B_flighty = QuasiMonteCarlo.generate_design_matrices(n_samples, flighty_lbs, flighty_ubs, sampler)
    # A_max, B_max = QuasiMonteCarlo.generate_design_matrices(n_samples, min_lbs, max_ubs, sampler)
    
    eFAST_persistent = gsa(output_func, eFAST(), [[persistent_lbs[i], persistent_ubs[i]] for i in 1:length(persistent_lbs)]; 
                       samples = 100)
    eFAST_flighty = gsa(output_func, eFAST(), [[flighty_lbs[i], flighty_ubs[i]] for i in 1:length(flighty_lbs)]; 
                       samples = n_samples)
    eFAST_max = gsa(output_func, eFAST(), [[min_lbs[i], max_ubs[i]] for i in 1:length(min_lbs)]; 
                       samples = n_samples)

    # # Get outputs for each type then join them
    # sobol_persistent = gsa(output_func, Sobol(), A_persistent, B_persistent)
    # sobol_flighty = gsa(output_func, Sobol(), A_flighty, B_flighty)
    # sobol_max = gsa(output_func, Sobol(), A_max, B_max)

    # Create dataframe to hold all all_results

    # Dataframe has columns for each parameter, each output, total and first order indices, and parameter set
    param_names = ["lQ", "lL", "lP", "lG", "sigma", "pQ", "pL", "pP", "pG"]
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

sample_sizes = 100:100:10000
results_list = Vector{DataFrame}(undef, length(sample_sizes))

@threads for idx in ProgressBar(eachindex(sample_sizes))
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

# Save outputs ----
# Parameter grid and outputs

open(joinpath((pwd()), "data", "julia_sobol.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, all_results)
    close(gzip_io)
end