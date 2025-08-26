using Distributions
using EMpht

data = rand(Exponential(1/10), 1_000)  # Generate some data to fit 
sample = EMpht.Sample(obs=data)        # Create an EMpht Sample object with this data
ph = empht(sample, p=5)                # Fit the data using p=5 phases

xGrid = range(0, 8, length=1_00)       # Create a grid to evaluate the density function over
fitPDFs = pdf.(ph, xGrid) 