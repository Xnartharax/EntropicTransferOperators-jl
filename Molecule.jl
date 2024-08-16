using NPZ
#data prep
subsamp = 10
lag = 10
dim = 30
# data from https://markovmodel.github.io/mdshare/ALA2
data = npzread("alanine-dipeptide-3x250ns-heavy-atom-positions.npz")
angledata = npzread("alanine-dipeptide-3x250ns-backbone-dihedrals.npz")
x = data["arr_0"][1:subsamp:end-lag+1, 1:dim]
y = data["arr_0"][lag:subsamp:end, 1:dim]
α = angledata["arr_0"][1:subsamp:end, :]




include("Single.jl")
include("Persistent.jl")
include("SpectrumPlots.jl")
include("Tracking.jl")

# entropic transfer operator
T = EntropicTransferOperator(x, y, 1e-2, solver=:stochastic, makesparse=false)

λ, vecs = spectrum(T, which=:LM, nev=100)

fig = plotspectrum(λ)

#visualize eigenvectors
i= 3
meas = real.(vecs[:, end-i])
l =  round(real(λ[end-i]); digits=2)
fig = Figure()
ax = Axis(fig[1, 1], title=L"Eigenvector for $\lambda_%$(i+1) = %$l$")
Makie.scatter!(ax, α[:, 1], α[:, 2], color=sign.(meas), markersize=5, )


# clustering
using Clustering
features = sign.(real.(vecs[:, end-3:end-1]))
R = kmeans(features', 4)
k1 = R.assignments .== 1
k2 = R.assignments .== 2
k3 = R.assignments .== 3
k4 = R.assignments .== 4


fig = Figure()
ax = Axis(fig[1,1], xlabel="ϕ", ylabel="ψ")
p1 =Makie.scatter!(ax,α[k1, 1], α[k1, 2], color=:red, markersize=2)
p2 =Makie.scatter!(ax,α[k2, 1], α[k2, 2], color=:blue, markersize=2)
p3 = Makie.scatter!(ax, α[k3, 1], α[k3, 2], color=:green, markersize=2)
p4 = Makie.scatter!(ax, α[k4, 1], α[k4, 2], color=:orange, markersize=2)

Legend(fig[1, 2], [p1, p2, p3], ["Cluster 1", "Cluster 2", "Cluster 3"])
