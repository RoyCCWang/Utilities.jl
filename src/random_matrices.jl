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
