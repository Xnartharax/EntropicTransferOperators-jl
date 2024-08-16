using DynamicalSystems

ds = PredefinedDynamicalSystems.lorenz()

include("Single.jl")
include("Persistent.jl")
include("SpectrumPlots.jl")
include("Tracking.jl")

# using dynamical systems integration to directly sample from trajectory
T = EntropicTransferOperator(ds, 5, 0.1; N=2000, sample=10)

plotspectrum(T)

# getting 30 largest real eigenvalues
λ, vecs = spectrum(T, which=:LM)

#visualize
using LaTeXStrings
i = 4
l = round(real.(λ[end-i]), digits=2)
fig = Figure()
ax = Axis3(fig[1, 1]; viewmode=:stretch, title=L"\lambda_%$i = %$l")
Makie.scatter!(ax,Matrix(T.discretization), color=sign.(real.(vecs[:, end-i])))

# persistent part
ϵs = 10.0 .^ collect(range(-3, 1, length=50))
Ts = PersistentEntropicTransferOperator(ds, ϵs, 0.1; N=1500, sample=10, makesparse=false, maxiter=10000)


spectrum_barcode(Ts)

λs, vecs = realspectrum(Ts)
ls = greedyscan3(λs, vecs, tol=0.02, α=0.1)

spectrum_barcode_lines(λs, ϵs, ls)