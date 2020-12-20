
function unpackdemo(D::Int = 2, M::Int = 3)

    W_array = collect( Utilities.generaterandomposdefmat(D) for m = 1:M )
    L_array = collect( cholesky(W_array[m]) for m = 1:M )

    b_array = collect( randn(D) for m = 1:M )
    c_array = rand(M)

    v = packbcW(b_array, c_array, W_array)

    b, c, W = unpackbcW(v,M)

    return v, b, c, W, b_array, c_array, W_array
end
