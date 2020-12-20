# unpacking and packing for optimization toolboxes.
function packbcW(b, c::Vector{T}, W) where T <: Real
    D = length(b[1])
    M = length(W)
    @assert length(b) == length(c) == M

    v = Vector{T}(undef,M+M*D+M*div(D*(D+1),2))

    st = 1
    fin = M
    v[st:fin] = c

    for m = 1:M
        v[fin+(m-1)*D+1:fin+m*D] = b[m]
    end
    st = fin+1
    fin = fin+M*D

    for m = 1:M
        L = cholesky(W[m]).L
        chol_c, chol_a = packCholeskyType(L)
        st = fin+1
        fin = fin + D
        v[st:fin] = chol_c

        st = fin+1
        fin = fin + div(D*(D+1),2)-D
        v[st:fin] = chol_a
    end

    return v
end

# forces c to be positive. Forces W to be posdef.
function unpackbcW(v::Vector{T}, M::Int, D::Int) where T <: Real

    st = 1
    fin = M
    c = v[st:fin]
    c = abs.(c) # force positive.

    b = collect( v[fin+(m-1)*D+1:fin+m*D] for m = 1:M )
    st = fin+1
    fin = fin+M*D

    W = Vector{Matrix{T}}(undef,M)
    for m = 1:M
        st = fin+1
        fin = fin + D
        chol_c = v[st:fin]
        chol_c = abs.(chol_c) # force positive.

        st = fin+1
        fin = fin + div(D*(D+1),2)-D
        chol_a = v[st:fin]

        L = parseCholeskyType(chol_c, chol_a)
        W[m] = L*L'
    end

    return b, c, W
end

# modifies L such that the lower triangular off-diagonal entries are replaced by a.
# The diagonal entries of L are unmodified.
function filloffdiagonalentriessym!(L::Matrix{T}, a::Vector{T})::Nothing where T <: Real
    N = size(L,1)
    @assert size(L,2) == size(L,1)

    p = 1
    for j = 1:N
        for i = j+1:N
            L[i,j] = a[p]
            p += 1
        end
    end
    @assert p-1 == round(Int, N*(N+1)/2-N) # debug.

    return nothing
end

function parseCholeskyType(c::Vector{T}, a::Vector{T})::Matrix{T} where T <: Real
    L = convert(Matrix{T}, LinearAlgebra.Diagonal(c))

    filloffdiagonalentriessym!(L, a)

    return L
end

# c is the diagonal, should be all positives.
function packCholeskyType(L)
    c = collect( L[i,i] for i = 1:size(L,1) )

    n = length(c)
    a = Vector{Float64}(undef, round(Int, n*(n+1)/2-n))

    p = 1
    for j = 1:n
        for i = j+1:n
            a[p] = L[i,j]
            p += 1
        end
    end
    @assert p-1 == round(Int, n*(n+1)/2-n) # debug.

    return c, a
end
