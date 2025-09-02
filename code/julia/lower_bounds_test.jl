# Test for lower bound values for probabilities and rates so that basic offspring exceeds one

include("GCD_R0_num_calc.jl")

# range for lower bound of probabilities
p_min = 0.0
p_max = 1.0

# range for lower bound of rates
r_min = 30.0
r_max = 1/2.0

# range for lower bound of host-seeking rate
g_min = ((1/2)*1440.0)
g_max = ((1/24)*1440.0)

# define number of steps to vary each parameter
num_steps = 1_001

# set up data frame for all parameter combinations
using DataFrames

param_grid = DataFrame(
    p = repeat(collect(range(p_min, p_max, length=num_steps)), inner=num_steps*num_steps),
    r = repeat(collect(range(r_min, r_max, length=num_steps)), inner=num_steps, outer=num_steps),
    g = repeat(collect(range(g_min, g_max, length=num_steps)), outer=num_steps*num_steps),
    N_offspring = Vector{Float64}(undef, num_steps^3),
    check = Vector{Float64}(undef, num_steps^3)
)

# calculate N_offspring for each combination of parameters
using Base.Threads
@threads for i in ProgressBar(1:nrow(param_grid))
    B_vals_in = (1/param_grid[i, :g], 1/param_grid[i, :r], 1/param_grid[i, :r], 1/param_grid[i, :r], param_grid[i, :p], param_grid[i, :p], param_grid[i, :p], param_grid[i, :p], param_grid[i, :p])
    N_offspring = repnums_func(B_vals_in)[1]
    param_grid[i, :N_offspring] = N_offspring
    param_grid[i, :check] = N_offspring .>= 1.0
end

# save results
# Save outputs to CSV
using CodecZlib
open(joinpath((pwd()), "data", "lower_bounds_test.csv"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, param_grid)
    close(gzip_io)
end