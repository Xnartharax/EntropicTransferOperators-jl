# TODO: Optimize this we only care about the spectrum here, memory usage is way to high
struct PersistentEntropicTransferOperator
    discretization::StateSpaceSet
    λs
    vecs
    ϵs::Vector{<:Real} 
end

function sparsify(A; threshold=1e-4)
    B = SparseArrays.similar(A)
    ind = A .> maximum(A)*threshold
    B[ind] = A[ind]
    B
end
function sparse_sinkhorn_plan(solver; threshold=1e-3)
    OptimalTransport.absorb!(solver)
    return sparsify(solver.cache.K, threshold=threshold)
end

function PersistentEntropicTransferOperator(x::StateSpaceSet, y::StateSpaceSet, ϵs::Vector{<:Real}; 
    makesparse=true, sparsitythresh=1e-3, distance=SqEuclidean(), nev =30, kwargs...)
    # a bit of technical fiddling with the library internals
    # basically an adaption of the sinkhorn function in 
    # https://github.com/JuliaOptimalTransport/OptimalTransport.jl/blob/master/src/entropic/sinkhorn_epsscaling.jl
    # Except I output the plan at every step
    M = length(x)
    cost = pairwise(distance, x,y)
    μx = ones(M)/M
    μy = ones(M)/M
    # ensure we start with the biggest epsilon
    sorting = sortperm(ϵs, rev=true)
    ϵs2 = ϵs[sorting]
    γs = []
    λs, vecs = [], []
    # initialize solver and perform Sinkhorn algorithm
    alg = SinkhornStabilized()
    solver = OptimalTransport.build_solver(μx, μy, cost, ϵs2[1], alg; kwargs...)
    OptimalTransport.solve!(solver)
    γ = OptimalTransport.sinkhorn_plan(solver)
    λ, vec = spectrum(γ*M)
    push!(λs, λ)
    push!(vecs, vec)
    # scaling over all ϵ using previous map as base
    for (i, εstep) in enumerate(ϵs2[2:end])
        solver = OptimalTransport.update_epsilon!(solver, εstep)
        OptimalTransport.solve!(solver)
        if makesparse
            γ = sparse_sinkhorn_plan(solver, threshold=sparsitythresh)
        else
            γ = OptimalTransport.sinkhorn_plan(solver)
        end
        λ, vec = spectrum(γ*M, nev=nev)
        push!(λs, λ)
        push!(vecs, vec)
    end
    # inverse sorting and save
    PersistentEntropicTransferOperator(x, λs[invperm(sorting)], vecs[invperm(sorting)], ϵs)
end




function PersistentEntropicTransferOperator(ds::DeterministicIteratedMap,  ϵs::Vector{<:Real};
    N::Int=1000, sample::Int=1, lookahead::Int=1, kwargs...)
    traj = trajectory(ds, N)
    x = traj[1:sample:end-lookahead]
    y = traj[1+lookahead:sample:end]
    return PeristentEntropicTransferOperator(x, y, ϵs; kwargs...)
end

function PersistentEntropicTransferOperator(ds::CoupledODEs, ϵs::Vector{<:Real}, Δt::Real; N=1000, sample=1, kwargs...)
    # TODO: more somphisticated sampling for example with Δt_sample and Δt_flow
    traj,t = trajectory(ds, N*Δt*sample, Δt=Δt)
    x = traj[1:sample:end-1]
    y = traj[2:sample:end]
    return PersistentEntropicTransferOperator(x, y, ϵs; kwargs...)
end


function realspectrum(Ts::PersistentEntropicTransferOperator)
    λs = []
    vecs = []
    for (λ, vec) in zip(Ts.λs, Ts.vecs)
        ind = isreal.(λ) .&& (real.(λ) .> 0)
        push!(λs, real.(λ[ind]))
        push!(vecs, vec[:, ind])
    end
    λs, vecs
end

function Base.getindex(T::PersistentEntropicTransferOperator, i)::EntropicTransferOperator
    return EntropicTransferOperator(T.ds, T.discretization, T.γs[i])
end

function Base.length(T::PersistentEntropicTransferOperator)
    return length(T.ϵs)
end

function get_vecs(xs, ys, vecs)
    vecs2 = []
    for (x, y) in zip(xs, ys)
        push!(vecs2, vecs[x][:, y])
    end
    vecs2
end

function average_vecs(vecs)
    avg = zeros(length(vecs[1]))
    signs = sign.(vecs[1])
    for vec in vecs
        if sum(sign.(vec) .!= signs) > length(vec)/2
            avg -= sign.(vec)
        else
            avg += sign.(vec)
        end
    end
    avg/length(vecs)
end

function average_vecs_weighted(vecs)
    avg = zeros(length(vecs[1]))
    signs = sign.(vecs[1])
    for (i, vec) in enumerate(vecs)
        if sum(sign.(vec) .!= signs) > length(vec)/2
            avg -= sign.(vec)/i
        else
            avg += sign.(vec)/i
        end
    end
    avg/sum(1/i for i=1:length(vecs))
end