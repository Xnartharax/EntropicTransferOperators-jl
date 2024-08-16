using StochasticOptimalTransport
using Random
using LogExpFunctions
using Distances

function ctransform(c, v::AbstractVector, x, ν, ε::Real)
    t = LogExpFunctions.logsumexp(
        (vᵢ - c(x, yᵢ)) / ε + log(νᵢ) for (vᵢ, yᵢ, νᵢ) in zip(v, ν.xs, ν.ps)
    )
    return -ε * (t)
end

function ctransform2(cost, μy, α::AbstractVector, ϵ::Real)
    N = length(α)
    t = LogExpFunctions.logsumexp.(
        (α[i] .- cost[i, :]) / ϵ .+ log(μy[i]) for i=1:N
    )
    return -ϵ * t
end

function solvestochastic(x::Matrix, y::Matrix, μx, μy, cost, ϵ;stepsize=nothing, dist=SqEuclidean())
    if stepsize === nothing
        stepsize= ϵ
    end
    N = length(μx)
    μ = StochasticOptimalTransport.DiscreteMeasure([x[i, :] for i=1:size(x, 1)], μx)
    ν = StochasticOptimalTransport.DiscreteMeasure([y[i, :] for i=1:size(y, 1)], μy)
    v = StochasticOptimalTransport.dual_v(Random.GLOBAL_RNG, dist, μ, ν, ϵ; maxiters=Int(1e6), stepsize=stepsize)
    u =  [ctransform(dist, v, x[i, :], ν, ϵ) for i=1:N]
    γ =  exp.((u .+ v' -cost)/ϵ)/N
    γ
end


function ∇h(cost, α, j, ϵ)
    tmp =(α .- cost[j, :])/ϵ
    return α - softmax(-tmp)
end
function entropicsag(cost, ϵ, s; iter=Int(1e6))
    # assumes symmetric i.e. dimensions of x and y are same and uniform weight
    N = size(cost, 1)
    α = ones(N)/N
    d = zeros(N)
    g = zeros(N, N)
    for k=1:iter
        j = Random.rand(1:N)
        d -= g[:, j]
        g[:, j] .= (1/N)*∇h(cost, α, j, ϵ)
        d += g[:, j]
        α += s*d
    end
    α
end
