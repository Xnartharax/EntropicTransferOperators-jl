n = 6
t = 4- 2.0^(-n+1)
function  F(x; t=4- 2.0^(-n+1))
    if x <=0.25
        return 2*x
    elseif  x < 1/2
        return t*(x-1/4)
    elseif  x < 3/4
        return t*(x-3/4)+1
    else
        return 2(x-1)+1
    end
end
N = 5000
x = collect(range(0, 1, length=N))
y = F.(x)


include("Single.jl")
include("Persistent.jl")
include("SpectrumPlots.jl")
include("Tracking.jl")

T = EntropicTransferOperator(x, y, 1e-5)

λ = eigvals(Matrix(T.γ))
Makie.scatter(real.(λ),imag.(λ))

Makie.lines!(0.5*cos.(0:0.001:2π), 0.5*sin.(0:0.001:2π), color=:red)
Makie.lines!([0.5+1/(5n), 0.5+1/(5n)], [-0.1, 0.1], color=:green)
Makie.lines!([0.5+1/(n^(2/3)), 0.5+1/(n^(2/3))], [-0.1, 0.1], color=:green)