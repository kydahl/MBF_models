# using Symbolics
# using Groebner
# using DynamicPolynomials
# using LinearAlgebra

# @polyvar a b c x1 x2 x3 x4

# system = [
#     a + b + c - x1,
#     a^2 + b^2 +c^2 - x2,
#     a^3 + b^3 + c^3 - x3,
#     a^4 + b^4 + c^4 - x4
# ]

# groebner(system)

# x = 3 * x1
# y = 3 * x2
# z = 3 * x3

# e1 = a + b + c - x
# e2 = a^2 + b^2 +c^2 - y
# e3 = a^3 + b^3 + c^3 - z


using HomotopyContinuation 

@var a b c # x1 x2 x3 x4
function EX(n, mu, sigma)
    # Calculate the nth moment of the inverse distribution of a lognormal distribution
    exp(-n * mu + (1//2) * n^2 * sigma^2)
end

mu = exp(2)
sigma = 0.5


EX1 = EX(1, mu, sigma)
EX2 = EX(2, mu, sigma)

x1 = 3 * EX(1, mu, sigma)
x2 = (3//2) * EX(2, mu, sigma)
x3 = (1//2) * EX(3, mu, sigma)
x4 = (1//8) * EX(4, mu, sigma)

f = System([
    a + b + c - x1,
    a^2 + b^2 +c^2 - x2,
    a^3 + b^3 + c^3 - x3,
    a^4 + b^4 + c^4 - x4
])

result = solve(f)

real_sols = real_solutions(result)