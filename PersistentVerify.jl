include("EntropicTransferOperator.jl")
include("Persistent.jl")

N = 1000
θ = sqrt(2)-1

function F(x)
    (x + θ) % 1
end

x = collect(range(0, 1, length=N))
y = F.(x)
function d(x, y)
    return min(((x % 1) - (y % 1))^2, ((x % 1) +1 - (y % 1))^2,  ((x % 1) -1 - (y % 1))^2)
end
ϵs = 10 .^(-1* (range(1, 6, length=50)))
#define conversion routine
Base.convert(T::Type{Float64}, s::SVector{1, Float64}) = s[1]

T = PersistentEntropicTransferOperator(StateSpaceSet(x), StateSpaceSet(y),ϵs, distance=d, nev=1000)

expectedspectrum(ϵ) = [ exp(-ϵ*π^2*abs(k)^2)*exp(-2*π*im*k*θ) for k=-N/2+1:N/2]

λs2 = expectedspectrum.(ϵs)
λs1 = T.λs
using Plots
i = 30
Makie.scatter(real.(λs2[i]), imag.(λs2[i]), markersize=1)
Makie.scatter!(real.(λs1[i]), imag.(λs1[i]), markersize=1)
