# Set up Sobol analysis
include("GCD_R0_num_calc.jl")

using GlobalSensitivity
# Define parameter ranges

# Rates at baseline
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
const base_params_flighty = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.9f0, 1.0f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0]

const base_params_persistent = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.66f0, 1.0f0, 0.7f0, 0.8f0, 0.9f0, 0.5f0]
    
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
min_lbs = [1/(3*1440.0f0), 1/(3*1440.0f0), 1/(3*1440.0f0), 1/(3*1440.0f0), 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0]
max_ubs = [60.0f0, 60.0f0, 60.0f0, 60.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 100.0f0]

# Set number of samples
n_samples = 2^15::Int #10_000::Int

# Create function to calculate basic offspring number and basic reproduction number at the same time
function output_func(B_vals_in)
    [N_offspring_func(B_vals_in), R0_func(B_vals_in)]
end

# Set up QuasiMonteCarlo sampling
sampler = SobolSample(R = OwenScramble(base = 2, pad = 19))


A_persistent, B_persistent = QuasiMonteCarlo.generate_design_matrices(n_samples, persistent_lbs, persistent_ubs, sampler)
A_flighty, B_flighty = QuasiMonteCarlo.generate_design_matrices(n_samples, flighty_lbs, flighty_ubs, sampler)
A_max, B_max = QuasiMonteCarlo.generate_design_matrices(n_samples, min_lbs, max_ubs, sampler)

# Get outputs for each type then join them
sobol_persistent = gsa(output_func, Sobol(), A_persistent, B_persistent)
sobol_flighty = gsa(output_func, Sobol(), A_flighty, B_flighty)
sobol_max = gsa(output_func, Sobol(), A_max, B_max)

# Create dataframe to hold all all_results
using DataFrames, CSV
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

for (sobol, paramset) in zip((sobol_persistent, sobol_flighty, sobol_max), paramset_types)
    for (i_out, outname) in enumerate(output_names)
        for (idx_type, arr) in zip(index_types, (sobol.ST, sobol.S1))
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

# Save outputs ----
# Parameter grid and outputs
using CodecZlib
open(joinpath((pwd()), "data", "julia_sobol.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, all_results)
    close(gzip_io)
end