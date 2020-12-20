### routines for modifying/transforming, or searching a collection.

# use a to decide whether b should be removed.
function filterother(a::Vector{T}, b::Vector{T}, f::Function)::Vector{T} where T

    keep_inds = Vector{Int}(undef,0)

    for i = 1:length(a)
        if f(a[i])
            push!(keep_inds, i)
        end
    end

    return b[keep_inds]
end

function selectklargest(x, y, k::Int)

    sorted_ind = sortperm(x, rev = true)

    return y[sorted_ind[1:k]], sorted_ind[1:k]
end

# extract columns that correspond to
function extractselectcols( A::Matrix{T},
                            inclusion_array::BitArray{1}) where T <: Real

    N_obs = count(inclusion_array)
    out = Vector{Vector{T}}(undef,length(inclusion_array))

    counter = 1
    for j = 1:length(inclusion_array)

        if inclusion_array[j]
            out[counter] = A[:,j]
            counter += 1
        end

    end
    resize!(out,counter-1)
    @assert length(out) == N_obs

    return out
end

function finddiscrepancy(a::Vector{Int}, b::Vector{Int})
    out = Vector{Int}(undef,0)

    for i = 1:length(a)
        if a[i] != b[i]
            push!(out,i)
        end
    end

    return out
end


# discards the later half of a sequence, use the first half to make a symmetric sequence.
# [1,2,3,4,5] -> [1,2,0,2,1]
# [1,2,3,4,5,6] -> [1,2,3,3,2,1]
function forcesymmetric(h0::Vector{T}) where T <: Number
    h = copy(h0)

    mid_index = div(length(h0),2)

    if rem(length(h),2) == 1
        # odd number of entries.
        h[mid_index+2:end] = reverse(h[1:mid_index])
    else
        # even number of entries.
        h[mid_index+1:end] = reverse(h[1:mid_index])
    end

    return h
end

# [1,2,3,4,5] -> [1,2,0,-2,-1]
# [1,2,3,4,5,6] -> [1,2,3,-3,-2,-1]
function forcenegativesymmetric(h0::Vector{T}) where T <: Number
    h = copy(h0)

    mid_index = div(length(h0),2)

    if rem(length(h),2) == 1
        # odd number of entries.
        h[mid_index+2:end] = reverse(-h[1:mid_index])
        h[mid_index+1] = zero(T)
    else
        # even number of entries.
        h[mid_index+1:end] = reverse(-h[1:mid_index])
    end

    return h
end

# used for converting a given point x ∈ [a,b] to u ∈ [1:N].
function coordinatetolattice1D(x::T, a::T, b::T, N::Int)::T where T

    Δ = (b-a)/(N-1)
    return (x-a)/Δ +1
end

function coordinatetolattice1D(x, a::T, b::T, N::Int)::Vector{T} where T
    return collect( coordinatetolattice1D(x[i],a,b,N) for i = 1:length(x) )
end

function lattice1Dtocoordinate(k::T, a::T, b::T, N::Int)::T where T
    Δ = (b-a)/(N-1)
    return (k-1)*Δ +a
end

function lattice1Dtocoordinate(k, a::T, b::T, N::Int)::Vector{T} where T
    return collect( lattice1Dtocoordinate(k[i],a,b,N) for i = 1:length(k) )
end


# tensor product filtering for dimension D discrete signal.
function ranges2collection( x_ranges::Vector{LinRange{T}},
                            ::Val{D}) where {T,D}
    # set up.
    @assert !isempty(x_ranges)
    @assert length(x_ranges) == D
    N_array = collect( length(x_ranges[d]) for d = 1:D )
    N = prod(N_array)
    sz_N = tuple(N_array...)

    # Position.
    X_nD = Array{Vector{T},D}(undef,sz_N)
    for 𝑖 in CartesianIndices(sz_N)
        X_nD[𝑖] = Vector{T}(undef,D)

        for d = 1:D
            X_nD[𝑖][d] = x_ranges[d][𝑖[d]]
        end
    end

    return X_nD
end

"""
avgsumovergrid( x_ranges::Vector{LinRange{T}},
                p::Function,
                ::Val{D})::T

Approximates an integral.
"""
function avgsumovergrid( x_ranges::Vector{LinRange{T}},
                        p::Function,
                            ::Val{D})::T where {T,D}
    # set up.
    @assert !isempty(x_ranges)
    @assert length(x_ranges) == D
    N_array = collect( length(x_ranges[d]) for d = 1:D )
    sz_N = tuple(N_array...)
    N = prod(N_array)

    # pre-allocate.
    x = Vector{T}(undef,D)

    # sum
    out = zero(T)
    for 𝑖 in CartesianIndices(sz_N)
        # assumble input.
        for d = 1:D
            x[d] = x_ranges[d][𝑖[d]]
        end

        out += p(x)/N
    end

    return out
end
