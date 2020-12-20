### routiens that exploit matrix structure for computational advantages.

"""
Assumes diagonal entries of L are all positive.
logdettriangular(L::AbstractMatrix{T})::T
"""
function logdettriangular(L::AbstractMatrix{T})::T where T <: Real
    return sum(log(L[i,i]) for i = 1:size(L)[1])
end


### routines for matrix transforms and operations.


"""
evalinverse!( Y::Matrix{T}.
                inv_E::Matrix{T},
                f::Vector{T},
                h::T)

Computes Y := inv(M), where M := [E f; f' h], where E is DxD, f is Dx1, h is 1x1.
Uses the Matrix inversion lemma; see Theorem 4.3.2 from Murphy's book.

Test:
M = Utilities.generaterandomposdefmat(3)
E = M[1:2,1:2]
f = vec(M[1:2,3])
h = M[3,3]
inv_E = inv(E)

Y = zeros(3,3)
evalinverse!( Y,
                inv_E,
                f,
                h)
println("inv(M)-Y = ", inv(M)-Y)
"""
function evalinverse!(  Y::Matrix{T},
                        inv_E::Matrix{T},
                        f::Vector{T},
                        h::T) where T <: Real
    #
    D = length(f)
    @assert size(inv_E) == (D,D)
    @assert size(Y) == (D+1,D+1)

    # compute the Schur complement of M w.r.t. E.
    schur = h - dot(f,inv_E*f)

    b = inv_E*f

    Y[1:D,1:D] = inv_E + (b*b')./schur
    Y[1:D, D+1] = -b ./ schur
    Y[D+1, 1:D] = -b ./ schur
    Y[D+1, D+1] = one(T)/schur

    return nothing
end

function evalinverse(   inv_E::Matrix{T},
                        f::Vector{T},
                        h::T)::Matrix{T} where T <: Real
    #
    Y = Matrix{T}(undef,size(inv_E))
    evalinverse!(Y, inv_E, f, h)

    return Y
end

function forcesymmetric(A::Matrix{T})::Matrix{T} where T <: Real
    return (A+A')./2
end

# replace the upper triangular entires with the lower ones.
function replaceupperbylowertriangular!(A::Matrix{T})::Nothing where T <: Real
    N = size(A,1)
    @assert size(A,2) == size(A,1)
    for j = 1:N
        for i = j+1:N
            A[j,i] = A[i,j]
        end
    end

    return nothing
end

function naivesqrtpsdmatrix(A)
    @assert LinearAlgebra.isposdef(A)
    s,Q = LinearAlgebra.eigen(A)

    s2 = sqrt.(s)

    B = Q'*LinearAlgebra.Diagonal(s2)*Q
    B = forcesymmetric(B)

    return B
end

function array2matrix(X::Array{Vector{T},L})::Matrix{T} where {T,L}

    N = length(X)
    D = length(X[1])

    out = Matrix{T}(undef,D,N)
    for n = 1:N
        out[:,n] = X[n]
    end

    return out
end


"""
    makeblockdiagonalmatrix(A::Vector{Matrix{T}})::Matrix{T}
"""
function makeblockdiagonalmatrix(A::Vector{Matrix{T}})::Matrix{T} where T <: Real
    K = length(A)

    N_rows = 0
    N_cols = 0

    for i = 1:K
        N_rows += size(A[i], 1)
        N_cols += size(A[i], 2)
    end

    out = zeros(T, N_rows, N_cols)
    st_rows = 1
    st_cols = 1

    for k = 1:K
        M,N = size(A[k])

        out[st_rows:st_rows+M-1,st_cols:st_cols+N-1] = A[k]

        st_rows += M
        st_cols += N
    end
    return out
end

"""
    invertblockmatrix!( out::Vector{Matrix{T}},
                    src::Vector{Matrix{T}})::Nothing
"""
function invertblockmatrix!(out::Vector{Matrix{T}},
                            src::Vector{Matrix{T}})::Nothing where T <: Real
    @assert length(out) == length(src)
    for i = 1:length(src)
        out[i] = inv(src[i])
    end
    return nothing
end

"""
    blockmatrixmatrixproduct(   A::Matrix{T},
                                B_set::Vector{Matrix{T}})::Matrix{T}

does A*makeblockdiagonalmatrix(B_set).

# Example:
```
D = 30
p = 30
B2_set = collect(randn(D,D) for i=1:p)
B2 = makeblockdiagonalmatrix( B2_set )
B = randn(D*p,D*p)
norm(B*B2 - blockmatrixmatrixproduct(B,B2_set))
```
"""
function blockmatrixmatrixproduct(  A::Matrix{T},
                                    B_set::Vector{Matrix{T}})::Matrix{T} where T <: Real
    K = length(B_set)

    M,N = size(B_set[1])
    @assert K*N == size(A,2)

    st = 1
    out = Matrix{T}(undef, size(A,1), K*M)
    for i = 1:length(B_set)
        # parse.
        st = (i-1)*size(B_set[i],1) + 1
        W = view(A, :, st:st+N-1)

        # matrix-vector product.
        out[:,st:st+N-1] = W*B_set[i]
    end

    return out
end


"""
    blockmatrixvectorproduct(   B_set::Vector{Matrix{T}},
                                v::Vector{T})::Vector{T}

does makeblockdiagonalmatrix(B_set)*v.

# Example:
```
D = 30
p = 30
B2_set = collect(randn(D,D) for i=1:p)
B2 = makeblockdiagonalmatrix( B2_set )
b = randn(D*p)
norm(B2*b - blockmatrixvectorproduct(B2_set,b))
```
"""
function blockmatrixvectorproduct(  B_set::Vector{Matrix{T}},
                                    v::Vector{T})::Vector{T} where T <: Real
    K = length(B_set)

    M,N = size(B_set[1])
    @assert K*N == length(v)

    st = 1
    out = Array{T}(K*M)
    for i = 1:length(B_set)
        # parse.
        st = (i-1)*size(B_set[i],2) + 1
        w = view(v, st:st+N-1)

        # matrix-vector product.
        out[st:st+N-1] = B_set[i]*w
    end

    return out
end


#### symmetric matrix efficient querying.

"""
    readuppertriangular(i::Int, j::Int, N::Int, v::Vector{T})::T

Query the (i-j)-th entry of an upper triangular matrix that is stored in v.

```
A = randn(7,7)

v = packageuppertriangle(A)

function packagemat(v, N)
    B = zeros(Float64, N, N)

    for j = 1:N
        for i = 1:j
            B[i,j] = readuppertriangular(i, j, N, v)
        end
    end

    return B
end

C = packagemat(v, size(A,1))
display(A-C)
```
"""
function readuppertriangular(i::Int, j::Int, N::Int, v::Vector{T})::T where T

    linear_ind = div(j*(j-1), 2) + i

    #Printf.@printf("(i,j) = (%d,%d), linear index = %d, offset = %d\n", i, j, linear_ind, offset)
    return v[linear_ind]
end

"""
    readsymmetric(i::Int, j::Int, N::Int, v::Vector{T})::T

Queries a symmetric matrix stored in a vector format.
"""
function readsymmetric(i::Int, j::Int, N::Int, v::Vector{T})::T where T
    if i < j
        return readuppertriangular(i, j, N, v)
    end

    return readuppertriangular(j, i, N, v)
end


"""
    writeuppertriangular!(v::Vector{T}, i::Int, j::Int, N::Int, val::T)::Nothing

Writes a value to an upper triangular matrix storage.
"""
function writeuppertriangular!(v::Vector{T}, i::Int, j::Int, N::Int, val::T)::Nothing where T

    linear_ind = div(j*(j-1), 2) + i
    v[linear_ind] = val

    return nothing
end

"""

```
D = 3
v = Vector{Float64}(undef, 6)
@time Utilities.writesymmetric!(v, 2, 3, 3, 1.23)

function fillv!(v, D)
    c = 0.0
    for j = 1:D
        for i = 1:j
            Utilities.writesymmetric!(v, i, j, D, c)
            c += 1
        end
    end

end

fillv!(v, D)
display(v)
```
"""
function writesymmetric!(v::Vector{T}, i::Int, j::Int, N::Int, val::T)::Nothing where T
    if i < j
        writeuppertriangular!(v, i, j, N, val)

        return nothing
    end

    writeuppertriangular!(v, j, i, N, val)

    return nothing
end

"""
    packageuppertriangle(L::Matrix{T})::Vector{T}

Extracts the upper triangular portion of a matrix. Includes diagonal entries.
"""
function packageuppertriangle(L::Matrix{T})::Vector{T} where T
    N = size(L,1)

    a = Vector{T}(undef, div(N*(N+1),2) )

    p = 1
    for j = 1:N
        for i = 1:j
            a[p] = L[i,j]
            p += 1
        end
    end

    return a
end
