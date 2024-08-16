using Plots
using LinearAlgebra
using OptimalTransport
using Distances

# definining the circle shift map
function F(x)
    (x + θ) % 1
end

#defining the distance measure
function d(x::Real, y::Real)
    return min(((x % 1) - (y % 1))^2, ((x % 1) +1 - (y % 1))^2,  ((x % 1) -1 - (y % 1))^2)
end

# sampling some points
N = 1000
θ = sqrt(2)-1
x = collect(range(0, 1, length=N))
y = F.(x)

include("Single.jl")
include("SpectrumPlots.jl")

# using both sinkhorn and stochastic solver
@time T1 = EntropicTransferOperator(x, y, 5e-5, distance=d, solver=:sinkhorn)

@time T2 = EntropicTransferOperator(x, y, 5e-5, distance=d, solver=:stochastic)

#visualize results
plotspectrum(T1)

# alternative sampling from a dynamical system
# this is slightly non-snensical here as the system is not ergodic and we can for now only sample along a single trajectory 
# since DynamicalSystems is not designed to change initial point and rerun very well.
using DynamicalSystems
x0 = [0.5]
function F(x, p, t)
    SVector((x[1] + θ) % 1)
end

ds = DeterministicIteratedMap(F, x0, [])

T3 = EntropicTransferOperator(ds, 5e-5, N=1000, sample=300, lookahead=10)

plotspectrum(T3)

#using qudratic regularization

T4 = EntropicTransferOperator(x, y, 5, distance=d, solver=:quadratic)
plotspectrum(T4)
