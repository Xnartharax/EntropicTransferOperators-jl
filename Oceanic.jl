using Rasters, NCDatasets

include("Single.jl")
include("Persistent.jl")
include("SpectrumPlots.jl")
include("Tracking.jl")

# load observational data (limited time range)
sst = Raster("sst.mnmean.nc"; name=:sst) 
# Data from https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc
replace_missing(sst, -9.96921f36)
# Taking only the latest high fidelity data and only the south-pacific
N = 600 # 50 years
bounds = X(28 .. 358-70), Y(-60 .. 20), Ti(2046-N:2046)
data1=Array(sst[bounds...])
mask = data1[:, :, 1] .!== missing
data = Float64.(data1[mask, :])
step = 12
x = data[:, 1:end-step]
y = data[:, step+1:end]

# computing transfer operator
T = EntropicTransferOperator(Matrix(x'), Matrix(y'), 500)

λ = eigvals(Matrix(T.γ))

plotspectrum(λ)