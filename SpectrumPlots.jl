using Makie, GLMakie
using LinearAlgebra
include("Single.jl")
include("Persistent.jl")
function spectrum_barcode(Ts::PersistentEntropicTransferOperator)
    λs, _ = realspectrum(Ts)
    ϵs = Ts.ϵs
    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, limits=(nothing, (0.8, 1.015)), xlabel="ϵ", ylabel="λ")
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=5, color=:blue)
    end 
    f
end
function spectrum_barcode_line(Ts::PersistentEntropicTransferOperator, xs, ys)
    λs, _ = realspectrum(Ts)
    ϵs = Ts.ϵs
    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, limits=(nothing, (0.8, nothing)), xlabel="ϵ", ylabel="λ")
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=5, color=:blue)
    end 
    y = []
    for (i, λ) in zip(xs, λs[xs])
        push!(ys, λ[i])
    end
    Makie.lines(xs, y)
    f
end
function spectrum_barcode_line(λs, ϵs, xs, ys)
    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, limits=(nothing, (0.8, nothing)), xlabel="ϵ", ylabel="λ")
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=3, color=:blue)
    end 
    y = []
    x = []
    for (i, j) in zip(ys, xs)
        push!(y, λs[j][i])
        push!(x, ϵs[j])
    end
    print(y)
    Makie.lines!(ax, x, y)
    f
end
function spectrum_barcode!(ax, λs, ϵs)
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=3, color=:blue)
    end 
end
function spectrum_barcode(λs, ϵs)
    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, limits=(nothing, (0.8, 1.015)),xlabel="ϵ", ylabel="λ")
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=3, color=:blue)
    end 
    f
end


function spectrum_barcode_lines(λs, ϵs, lines)
    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, limits=(nothing, (0.8, 1.015)),xlabel="ϵ", ylabel="λ")
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=3, color=:blue)
    end 
    plots = []
    for (xs, ys) in lines
        y = []
        x = []
        for (i, j) in zip(ys, xs)
            push!(y, λs[j][i])
            push!(x, ϵs[j])
        end
        print(y)
        push!(plots, Makie.lines!(ax, x, y))
    end
    Legend(f[1, 2], plots, [L"λ_{%$i}" for i=1:length(plots)])
    f
end

function spectrum_barcode_lines!(ax, λs, ϵs, lines)
    for (i, ϵ) in enumerate(ϵs)
        y = λs[i]
        x = ones(length(y))*ϵ
        Makie.scatter!(ax, x, y, markersize=3, color=:blue)
    end 
    lplots = []
    for (xs, ys) in lines
        y = []
        x = []
        for (i, j) in zip(ys, xs)
            push!(y, λs[j][i])
            push!(x, ϵs[j])
        end
        pl = Makie.lines!(ax, x, y)
        push!(lplots, pl)
    end
    lplots
end

function plotspectrum(T::EntropicTransferOperator)
    λ = eigvals(Matrix(T.γ))
    plotspectrum(λ)
end

function plotspectrum(λ::Vector{<:Complex}; circ=true, markersize=3 , color=:blue)
    fig = Makie.Figure()
    ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im", aspect=1)
    if circ 
        θ= 0:0.01:2π
        x = cos.(θ)
        y = sin.(θ)
        poly!(collect(zip(x, y)); color=:lightblue)
    end
    plotspectrum!(ax, λ; markersize=markersize,color=color)
    fig
end
function plotspectrum!(ax, λ; markersize=3 , color=:blue)
    Makie.scatter!(ax, real.(λ), imag.(λ), markersize=markersize,color=color)
end