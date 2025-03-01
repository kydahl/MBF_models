# Set up PRCC analysis
include("GCD_R0_num_calc.jl")

# Rates
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
const base_params_flighty = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.9f0, 1.0f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0]

const base_params_persistent = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.66f0, 1.0f0, 0.7f0, 0.8f0, 0.9f0, 0.5f0]
    

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

# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in
min_lbs = [1/(3*1440.0f0), 1/(3*1440.0f0), 1/(3*1440.0f0), 1/(3*1440.0f0), 0.2f0, 0.2f0, 0.2f0, 0.2f0, 0.2f0, 0.0f0]
max_ubs = [60.0f0, 60.0f0, 60.0f0, 60.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 100.0f0]

# Set number of LHC samples
n_samples = 10_000::Int

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
    n_samples = size(LHS_samples)[2]
    GCD_results = Vector{Float64}(undef, n_samples::Int)
    N_offspring_results = Vector{Float64}(undef, n_samples::Int)
    R0_results = Vector{Float64}(undef, n_samples::Int)

    # Evaluate the function across all input combinations in parallel

    @threads for idx in ProgressBar(1:n_samples)#(i, (sigma, lQ, lL, lP, lG, pQ, pL, pP, pG)) in ProgressBar(enumerate(parameter_grid))
        # grid_point = LHS_samples[:,idx]
        # (lQ, lL, lP, lG, sigma, pL, pP, pG) = grid_point
        B_vals = LHS_samples[1:9,idx]
        # GCD values
        GCD_results[idx] = GCD_func(B_vals)
        # Basic offspring number values
        curr_N_offspring = N_offspring_func(B_vals)
        N_offspring_results[idx] = curr_N_offspring

        # R0 values
        # LVBEpiHost_vals = [LVB_vals; Epi_vals; Host_vals]
        R0_results[idx] = R0_func(B_vals)

    end
    scaled_plan_df = DataFrame(transpose(LHS_samples), [:lQ, :lL, :lP, :lG, :sigma, :pQ, :pL, :pP, :pG, :dummy])
    output_df = scaled_plan_df
    output_df[!,:GCD] = GCD_results
    output_df[!,:N_offspring] = N_offspring_results
    output_df[!,:R0] = R0_results

    return output_df
end

# Load output function
output_calc(QuasiMonteCarlo.sample(100, min_lbs, max_ubs, LatinHypercubeSample()))

# Get outputs for each type then join them
max_results = output_calc(max_scaled_plan)
# flighty_results = output_calc(flighty_scaled_plan)
# persistent_results = output_calc(persistent_scaled_plan)

max_results[!, :type] .= "max"
# flighty_results[!, :type] .= "flighty"
# persistent_results[!, :type] .= "persistent"

all_results = vcat(max_results) #, flighty_results, persistent_results)

# Save outputs ----
# Parameter grid and outputs
using CodecZlib
open(joinpath(dirname(dirname(pwd())), "data", "julia_outputs.csv.gz"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, all_results)
    close(gzip_io)
end

# using Plots
# # plot(all_results[!,:sigma],all_results[!,:R0], seriestype=:scatter, markersize =1, alpha = 0.5, color=:blue)
# using Measures
# input_names = names(all_results[!,1:8])
# output_names = names(all_results[!,9:11])

# # Iterate over output variables
# for output_name in output_names
#     output_label = join(["Log10 of", output_name], " ")  # Get the name of the output variable
#     # Create a grid of 8 plots for the current output variable
#     p = plot(layout=(2, 4), size = (1600, 1200), margin = 10mm)
#     subplot_index = 1
#     for input_name in input_names
#         # Extract data for the current pair of columns
#         x = all_results[!, input_name]
#         y = all_results[!, output_name]
#         # Transform y to base-10 logarithm
#         log_y = log10.(y)
#         # Add subplot for the current input variable
#         scatter!(x, log_y, subplot=subplot_index,
#                  xlabel=input_name, ylabel=output_label, markersize=1, legend = false)
#         subplot_index = subplot_index + 1
#     end

#     # Save the grid of plots
#     savefig(p, joinpath(dirname(dirname(pwd())), "figures", "grid_output_var$output_name.png"))
# end

# # Calculate PRCCs ----
# using StatsBase
# # Function to rank transform a DataFrame column by column
# function rank_transform_df(df::DataFrame)
#     rank_df = copy(df)  # Create a copy of the DataFrame
#     @threads for col in ProgressBar(names(df))
#         rank_df[!, col] = sortperm(sortperm(df[!, col]))  # Apply rank transform
#     end
#     return rank_df
# end


# # Rank transform inputs and outputs
# max_rank_df = rank_transform_df(max_results)
# flighty_rank_df = rank_transform_df(flighty_results)
# persistent_rank_df = rank_transform_df(persistent_results)

# function PRCC_calc(df::DataFrame)
#     input_names = names(df[!,1:8])
#     output_names = names(df[!,9:11])
#     # Compute correlations and store in a results DataFrame
#     results = DataFrame(Input = String[], Output = String[], Correlation = Float64[])

#     rank_df = rank_transform_df(df)
#     # Iterate over output variables
#     for output_name in output_names
#         for input_name in input_names
#             corr_value = cor(rank_df[!,input_name], rank_df[!,output_name])
#             println("Correlation between $input_name and $output_name: $corr_value")
#             push!(results, (string(input_name), string(output_name), corr_value))
#         end

#     end
# end

# # Save PRCC results
# CSV.write(joinpath(dirname(dirname(pwd())), "data", "PRCC_results.csv"), results)


# # Plot correlations for each output
# for output_name in output_names
#     # Filter the results for the current output
#     current_results = filter(row -> row.Output == string(output_name), results)
#     p = plot(size = (3200, 2400), margin = 10mm)

#     # Extract values for the bar plot
#     x_labels = current_results.Input
#     y_values = current_results.Correlation

#     # Create the bar plot
#     p = bar(
#         x_labels,
#         y_values,
#         title = "Correlations with $output_name",
#         xlabel = "Input Variables",
#         ylabel = "Correlation",
#         legend = false
#     )

#     # Display the plot
#     display(plot())
#     # Save the grid of plots
#     savefig(p, joinpath(dirname(dirname(pwd())), "figures", "PRCC_$output_name.png"))
# end
