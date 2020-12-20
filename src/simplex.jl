
### ℝ^{D-1} to D-dim simplex.

# Takes x ∈ (0,1)^D to somewhere within the interior of the D-dim simplex,
#   then we map it to the boundary of the (D+1)-dim simplex.
function unitcubetosimplexinterior(x::Vector{T})::Vector{T} where T
    D = length(x)

    s = Vector{T}(undef,D)
    for i = 1:D-1
        s[i] = x[i]*prod( one(T)-x[j] for j = i+1:D)
    end
    s[D] = x[D]

    return s
end

function unitcubetosimplex(x::Vector{T})::Vector{T} where T
    D = length(x)

    s = unitcubetosimplexinterior(x)
    push!(s, one(T) - sum(s[i] for i = 1:D))

    return s
end

function Euclideantounitcube(r::Vector{T})::Vector{T} where T
    D = length(r)

    x = Vector{T}(undef,D)
    for i = 1:D
        x[i] = realtounitinterval(r[i])
    end

    return x
end

function realtounitintervalV(x::T)::T where T
    return (one(T) + x/sqrt(one(T)+x^2))/2-0.5
end

function realtounitinterval2(x::T, a::T)::T where T
    return (one(T) + a*x/sqrt(one(T)+(a*x)^2))/2
end

function realtounitinterval(x::T)::T where T
    return (one(T) + x/sqrt(one(T)+x^2))/2
end
# function realtounitinterval(x::T)::T where T
#     return logicsticfunc(x)
# end

function derivativeofrealtounitinterval(x::T)::T where T
    return one(T)/(2*sqrt(one(T)+x^2)^3)
end

# package up.
function Euclideantosimplex(r::Vector{T})::Vector{T} where T
    x = Euclideantounitcube(r) # to the unit cube.
    s = unitcubetosimplex(x)

    return s
end

function Euclideantosimplexinterior(r::Vector{T})::Vector{T} where T
    x = Euclideantounitcube(r) # to the unit cube.
    s = unitcubetosimplexinterior(x)

    return s
end

### other direction: D-dim simplex to ℝ^{D-1}.
function simplextounitcube(s::Vector{T})::Vector{T} where T
    D = length(s) -1

    x = Vector{T}(undef,D)
    for i = 1:D-1
        x[i] = s[i]/(one(T) - sum( s[j] for j = i+1:D))
    end
    x[D] = s[D]

    return x
end

function simplexinteriortounitcube(s::Vector{T})::Vector{T} where T
    D = length(s)

    x = Vector{T}(undef,D)
    for i = 1:D-1
        x[i] = s[i]/(one(T) - sum( s[j] for j = i+1:D))
    end
    x[D] = s[D]

    return x
end

function unitcubetoEuclidean(x::Vector{T})::Vector{T} where T
    D = length(x)

    r = Vector{T}(undef,D)
    for i = 1:D
        r[i] = unitintervaltoreal(x[i])
    end

    return r
end

# inverse of realtounitinterval()
function unitintervaltoreal(y::T)::T where T
    u = 2*y-1
    return u/sqrt(one(T)-u^2)
end
# # test code.
# x = randn()
# h = realtounitinterval(x)
# x_rec = unitintervaltoreal(h)
# println("x = ", x)
# println("x_rec = ", x_rec)
# println("discrepancy = ", abs(x-x_rec))
# function unitintervaltoreal(y::T)::T where T
#     return invlogicsticfunc(y)
# end

# package up.
function simplextoEuclidean(s::Vector{T})::Vector{T} where T
    x = simplextounitcube(s) # to the unit cube.
    r = unitcubetoEuclidean(x)

    return r
end

function simplexinteriortoEuclidean(s::Vector{T})::Vector{T} where T
    x = simplexinteriortounitcube(s) # to the unit cube.
    r = unitcubetoEuclidean(x)

    return r
end




###### Jacobian of Euclideantosimplexinterior()

# derivative of realtounitinterval().
function ∂realtounitinterval(x::T)::T where T <: Real
    tmp = 2*sqrt(1+x^2)*(1+x^2)
    return one(T)/tmp
end

function ln∂realtounitinterval(x::T)::T where T <: Real
    tmp = log(2) + 3/2*log((1+x^2))
    return -tmp
end

# derivative of simplexinteriortounitcube().
function ∂simplexinteriortounitcube(u::Vector{T}, a::Int, b::Int) where T <: Real
    D = length(u)
    if a != b
        return -u[a]*prod( one(T)-u[j] for j = a+1:D if j != b)
    end

    if a == D
        return one(T)
    end

    return prod( one(T)-u[j] for j = a+1:D)
end

# ℝ^D to boundary of D-simplex.
function computeJacobiandiagforW( y::Vector{T}) where T <: Real
    D = length(y)

    # pre-compute.
    one_minus_u = collect( one(T) - realtounitinterval(y[i]) for i = 1:D )

    # compute Jacobian.
    J_diag = Vector{T}(undef,D)
    for i = 1:D-1

        factor1 = ∂realtounitinterval(y[i])
        factor2 = prod( one_minus_u[j] for j = i+1:D)

        J_diag[i] = factor1*factor2
    end
    J_diag[D] = ∂realtounitinterval(y[D])

    return J_diag
end

function computelogdetJacobianforW( y::Vector{T})::T where T <: Real
    D = length(y)

    # pre-compute.
    ln_one_minus_u = collect( log(one(T) - realtounitinterval(y[i])) for i = 1:D )

    # compute Jacobian.
    ln_det_J = zero(T)
    for i = 1:D-1

        ln_factor1 = ln∂realtounitinterval(y[i])
        ln_factor2 = sum( ln_one_minus_u[j] for j = i+1:D)

        update = ln_factor1 + ln_factor2
        ln_det_J += update
    end
    ln_det_J += ln∂realtounitinterval(y[D])

    return ln_det_J
end
