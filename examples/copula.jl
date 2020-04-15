
import HCubature
import PyPlot
import Random
import GSL
import LinearAlgebra
import SpecialFunctions
import Distributions

include("../src/drawRV.jl")
include("../src/data_structures.jl")
include("../src/matrix.jl")


PyPlot.close("all")
fig_num = 1

Random.seed!(25)

D = 2

ε = 1e-4
# xmin = ε*ones(2)
# xmax = [25.0; 1.0-ε]
# (val,err) = HCubature.hcubature(evalcopulapdf, xmin, xmax;
#                 norm = LinearAlgebra.norm, rtol = sqrt(eps(Float64)), atol = 0,
#                 maxevals = typemax(Int), initdiv=1)
# # to do: use semi-infinite domains by a change-of-variable.
# # https://github.com/stevengj/cubature/blob/master/README.md#infinite-intervals
# println("Integration value is ", val)

# verify beta cdf.

## specify marginals.
dist_beta = Distributions.Beta(4.0, 2.0)
ln_pdf_beta = xx->Distributions.logpdf(dist_beta, xx)
cdf_beta = xx->Distributions.cdf(dist_beta, xx)

dist_gamma = Distributions.Gamma(2.0, 1.0)
ln_pdf_gamma = xx->Distributions.logpdf(dist_gamma, xx)
cdf_gamma = xx->Distributions.cdf(dist_gamma, xx)

ln_pdfs = Vector{Function}(undef,D)
ln_pdfs[1] = ln_pdf_gamma
ln_pdfs[2] = ln_pdf_beta

sum_ln_pdf = xx->sum( ln_pdfs[d](xx[d]) for d = 1:D )

cdfs = Vector{Function}(undef,D)
cdfs[1] = cdf_gamma
cdfs[2] = cdf_beta

## put together with the copula.
off_diagonal_entries = zeros(Float64, fld(D*(D-1),2)) #rand(Float64, fld(D*(D-1),2))

f = xx->exp( evallncopulapdf(xx, off_diagonal_entries, cdfs, sum_ln_pdf, Val(:GaussianCopula)) )


xmin = ε*ones(D)
xmax = [1.0-ε; 25.0-ϵ]
(val,err) = HCubature.hcubature(f, xmin, xmax;
                norm = LinearAlgebra.norm, rtol = sqrt(eps(Float64)), atol = 0,
                maxevals = typemax(Int), initdiv=1)
# to do: use semi-infinite domains by a change-of-variable.
# https://github.com/stevengj/cubature/blob/master/README.md#infinite-intervals
println("Integration value is ", val)



# I am here.
@assert 1==2

# generate a heat map image
N_per_dim = 100
x1_coords = reverse(collect(LinRange(ε,25.0,N_per_dim))) # need to reverse it so origin is at bottom left corner.
x2_coords = collect(LinRange(ε,1.0-ε,N_per_dim))


Y = zeros(Float64, N_per_dim, N_per_dim)
for x1_ind = 1:N_per_dim
    for x2_ind = 1:N_per_dim
        x = [x1_coords[x1_ind]; x2_coords[x2_ind]]
        Y[x1_ind,x2_ind] = exp( evalcopulapdf(x, off_diagonal_entries, ) )
    end
end

# plot grayscale heat map
x1_mesh = repeat(x1_coords,1,N_per_dim)
x2_mesh = repeat(x2_coords',N_per_dim,1)
PyPlot.pcolormesh(x1_coords, x2_coords, Y, cmap = "Greens_r")
#PyPlot.imshow(Y, interpolation = "nearest", cmap="Greys_r")
PyPlot.title("copula pdf")
