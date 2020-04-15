## routines that relate to the generation of random objects.


function drawnormal(μ::T, σ::T)::T where T <: Real
    return randn()*σ + μ
end

# draw isotropic normal.
function drawnormal(μ::Vector{T}, σ::T)::Vector{T} where T <: Real
    return collect( randn()*σ + μ[d] for d = 1:length(μ) )
end


# returns integer vector of K elements, each randomly picked (with replacement) in 1:1:N
# this is drawing K IID realizations from a uniform PMF from with support 1:1:N
function drawfromuniformindexinterval(N::Int, K::Int)
    a = rand(K)*N;
    return int(ceil(a));
end

# generates a random positive definite matrix.
function generaterandomposdefmat(n)
    A = rand(n,n)
    out = A'*A
    return forcesymmetric(out)
end

function generaterandomsym(n)
    return extractsym(rand(n,n))
end

function generaterandomskew(n)
    return extractskew(rand(n,n))
end

function generaterandomortho(n)
    A = randn(n,n)
    B = A*A'
    s,Q = LinearAlgebra.eigen(B)
    return Q
end


function generateStiefel(n,p)
    Q = generaterandomortho(p)
    out = zeros(Float64, n,p)
    out[1:p,1:end] = Q
    return out
end

function generaterandomhankel(n)
    a = rand(2*n-1)
    out = Matrix{Float64}(undef,n,n)
    for j = 1:n
        for i = 1:n
            out[i,j] = a[i+j-1]
        end
    end
    return out
end

# function generateisotopysubgroup(n,p)
#     O = generaterandomortho(n-p)
#     out = zeros(Float64, n,n)
#     out[1:p,1:p] = eye(p)
#     out[(n-p):end, (n-p):end] = O
#     return out
# end

# generate an invertible nxn indefinite matrix, with k positive eigenvalues.
function generateindefmatrix(n,k=rand(2:n-1))
    @assert 2 <= k <= n

    U,S,V = LinearAlgebra.svd(rand(n,n))
    for i = 1:k
        S[i] = rand() # rand() is between [0,1]
    end
    for i = k+1:length(S)
        S[i] = -rand()
    end
    return U*LinearAlgebra.diagm(0=>S)*U'
end

function generateposdefmatrix(n)
    U,S,V = LinearAlgebra.svd(rand(n,n))
    S = rand(n)
    out = U*LinearAlgebra.diagm(0=>S)*U'

    return forcesymmetric(out)
end

function generatenegdefmatrix(n)
    U,S,V = LinearAlgebra.svd(rand(n,n))
    S = -rand(n)
    out = U*LinearAlgebra.diagm(0=>S)*U'

    return forcesymmetric(out)
end



### copula.

function probit(p::T)::T where T <: Real
    @assert zero(T) <= p <= one(T)
    return sqrt(2)*SpecialFunctions.erfinv(2*p-one(T))
end

function Gaussiancopula(u, Σ::Matrix{T}) where T <: Real
    N = length(u)
    @assert size(Σ) == (N,N)
    I_mat = LinearAlgebra.Diagonal{T}(LinearAlgebra.I, N)

    if all(isfinite.(u))
        return exp(-0.5*LinearAlgebra.dot(u,(inv(Σ)-I_mat)*u))/sqrt(LinearAlgebra.det(Σ))
    end

    return zero(T)
end

### TODO this should not use Distributions.jl
# function evalGaussiancopula(  x,
#                             R::Matrix{T},
#                             marginal_dists) where T <: Real
#     D = length(x)
#
#     PDF_eval_x = collect( Distributions.pdf(marginal_dists[d], x[d]) for d = 1:D )
#     probit_CDF_eval_x = collect( probit( Distributions.cdf(marginal_dists[d], x[d]) ) for d = 1:D )
#     CDF_eval_x = collect( Distributions.cdf(marginal_dists[d], x[d]) for d = 1:D )
#
#
#     #println(probit_CDF_eval_x)
#     out = Gaussiancopula(probit_CDF_eval_x,R)*prod(PDF_eval_x)
#
#     if isfinite(out) != true
#         println(prod(PDF_eval_x))
#         println(probit_CDF_eval_x)
#         println(CDF_eval_x)
#         println(x)
#     end
#     @assert isfinite(out) == true
#     return out
# end



### unused.

# algorithm 4.3 from Simulating Copulas
function drawfromGaussiancopula(inv_cdfs::Vector,
                                R::Matrix{T}) where T <: Real
    #
    D = size(R)[1]

    L = cholesky(R).L
    z = collect( randn() for d = 1:D )
    x = L*z
    for i = 1:D
        x[i] = inv_cdfs[i](inverseprobit(x[i]))
    end

    return x
end

function draw2D𝓝copulaγmarginals(  shape_parameters::Vector{T},
                                    scale_parameters::Vector{T},
                                    N::Int) where T <: Real
    #
    @assert length(shape_parameters) == length(scale_parameters)

    inv_cdfs = collect( xx->GSL.cdf_gamma_Pinv(xx, shape_parameters[d], scale_parameters[d]) for d = 1:D )
    X = collect( drawfromGaussiancopula(inv_cdfs, R) for n = 1:N )

    return X
end


####


# cdf( randn() ), cdf is that of the standard normal.
function drawstdnormalcdf()::Float64
    return 0.5*( 1.0 + SpecialFunctions.erf(randn()/sqrt(2)) )
end

function drawstdnormalcdf!(x::Vector{T}) where T <: Real
    for i = 1:length(x)
        x[i] = convert(T, drawstdnormalcdf())
    end

    return nothing
end
