# PRCC analysis of model outputs
using CSV
using CodecZlib
using DataFrames

# Load the gzipped CSV file
output_df = CSV.read(GzipDecompressorStream(open(joinpath(dirname(dirname(pwd())), "data", "julia_outputs.csv.gz"))), DataFrame)


# Quickly check monotonicity of outputs wrt each variable
using Plots

# downsampling
using StatsBase
using Measures
sample_indices = sample(1:n_samples, 10000, replace=false)
sampled_df = output_df[sample_indices,1:11]
