using DynamicalSystems
using Distances
using SparseArrays
using OptimalTransport
using ArnoldiMethod
using Makie, GLMakie
using ProgressMeter
using Suppressor
include("StochasticFixes.jl")


struct EntropicTransferOperator
    discretization::StateSpaceSet
    γ::Union{Matrix, SparseMatrixCSC}
end

# TODO implement Base.show

# these convience wrappers basically only extract the discretization and evolution and pass the construction on
# TODO: implement other discretization options than trajectory

function EntropicTransferOperator(ds::DeterministicIteratedMap, ϵ::Real;
         N::Int=1000, sample::Int=1, lookahead::Int=1, kwargs...)
    traj = trajectory(ds, N*sample)[1]
    x = traj[1:sample:end-lookahead]
    y = traj[1+lookahead:sample:end]
    return EntropicTransferOperator(x, y, ϵ; kwargs...)
end

function EntropicTransferOperator(ds::CoupledODEs, ϵ::Real, Δt::Real; N=1000, sample=1, kwargs...)
    # TODO: more somphisticated sampling for example with Δt_sample and Δt_flow
    traj,t = trajectory(ds, N*Δt*sample, Δt=Δt)
    x = traj[1:sample:end-1]
    y = traj[2:sample:end]
    return EntropicTransferOperator(x, y, ϵ; kwargs...)
end
# this is useful if one-dimensional state spaces are required
Base.convert(T::Type{Float64}, s::SVector{1, Float64}) = s[1]

function EntropicTransferOperator(x::Union{AbstractArray{<:Real}, StateSpaceSet}, y::Union{AbstractArray{<:Real}, StateSpaceSet}, ϵ::Real; 
    makesparse=true, sparsitythresh=1e-3, distance=SqEuclidean(), solver=:sinkhorn, kwargs...)
    # intitialize transport problem
    
    if typeof(x) <: Matrix 
        M = size(x, 1)
        cost = pairwise(distance, x', y')
    else
        M = length(x)
        cost = pairwise(distance, x, y)
    end
    μx = ones(M)/M
    μy = ones(M)/M
    # solve optimal transport
    # TODO: figure out how to best pass the kwargs to OT solver
    if solver== :stochastic
        γ = solvestochastic(Matrix(x), Matrix(y), μx, μy, cost, ϵ, dist=distance, stepsize=max(ϵ, 1/N))
    elseif solver==:sinkhorn
        steps = Int(ceil(log(1/2, ϵ/maximum(cost))))
        γ =  sinkhorn(
            μy,
            μx,
            cost,
            ϵ,
            SinkhornEpsilonScaling(
                SinkhornStabilized();
                factor=0.5,
                steps=steps,
            ), maxiter=1000
        )*M
    elseif solver==:quadratic
        γ = quadreg(μx, μy, cost, ϵ)*M
    else
        error("No such solver technique $solver. Supported arguments are :sinkhorn, :stochastic and :quadratic.")
    end

    if makesparse   
        γ[γ .< maximum(γ)*sparsitythresh] .= 0
        γ = sparse(γ)
    end
    return EntropicTransferOperator(StateSpaceSet(x), γ)
end


function spectrum(T::EntropicTransferOperator; nev::Int=30,  which=:LR)
    vals, vecs = partialeigen(partialschur(T.γ; nev=nev, which=which)[1])
    return vals, vecs
end
function spectrum(A::Union{Matrix{<:Real}, SparseMatrixCSC}; nev::Int=30, which=:LR)
    vals, vecs = partialeigen(partialschur(A; nev=nev, which=which)[1])
    return vals, vecs
end

function realspectrum(T::EntropicTransferOperator; nev::Int=30)
    realspectrum(T.γ, nev=nev)
end


function realspectrum(A::Union{AbstractMatrix{<:Real}, SparseMatrixCSC}; nev::Int=30)
    vals, vecs = partialeigen(partialschur(A; nev=nev, which=:LR)[1])
    ind = isreal.(vals) .&& (real.(vals) .> 0)
    return real.(vals[ind]), real.(vecs[:, ind])
end

function sparsify(A; threshold=1e-4)
    B = SparseArrays.similar(A)
    ind = A .> maximum(A)*threshold
    B[ind] = A[ind]
    B
end
function sparse_sinkhorn_plan(solver; threshold=1e-3)
    OptimalTransport.absorb!(solver)
    return sparsify(Array(solver.cache.K), threshold=threshold)
end

