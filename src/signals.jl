### related to signal analysis/synthesis, geometric transforms.

function rad2freq(r,fs)
    fn = fs/2
    #r/π = x/fn, x is in Hz.
    return fn*r/π
end

function freq2rad(x,fs)
    fn = fs/2

    return x*π/fn
end


# evaluates H(e^{j*ω}).
function computeDTFTviaformula(h,ω)
    N = length(h)
    return sum( h[n+1]*exp(-im*ω*n) for n = 0:N-1 )
end
# # test code
# ω = 99
# h = randn(5)
# Hn = h[1] + h[2]*exp(-im*ω) + h[3]*exp(-im*2*ω) + h[4]*exp(-im*3*ω) + h[5]*exp(-im*4*ω)
# Ht = computeDTFTviaformula(h,ω)
#
# println("discrepancy ", abs(Hn-Ht))

# length = 2N+1
function getdiracdeltaseq(N)
    out = zeros(Float64,2*N+1)
    out[N+1] = 1.0

    return out
end

# get highpass filter from a lowpass one.
function spectralinversion(h::Vector{T}) where T
    @assert isodd(length(h))

    out = collect( -h[i] for i = 1:length(h))
    out[div(length(h),2)+1] += one(T)

    return out
end

# assume x starts at time 0.
function upsample(x::Vector{T}, M::Int)::Vector{T} where T

    out = Vector{T}(undef, length(x)*M)
    #out = Vector{T}(undef, (length(x)-1)*M +1)
    fill!(out,zero(T))

    for n = 0:length(x)-1
        out[n*M+1] = x[n+1]
    end

    return out
end

function downsample(x::Vector{T}, M::Int)::Tuple{Vector{T},Vector{T}} where T

    x_padded = copy(x)
    Q = rem(length(x),M)
    if Q != 0
        #push!(out,zero(T))
        push!(x_padded,zeros(T,M-Q)...)
    end

    N = div(length(x_padded),M)

    out = Vector{T}(undef, N)

    for n = 0:N-1
        out[n+1] = x_padded[n*M+1]
    end

    return out, x_padded
end


function linearconv(  x::Vector{T},
                h::Vector{T}) where T

    out = Vector{T}(undef, length(x)+length(h)-1)
    for i = 1:length(x)+length(h)-1

        out[i] = sum( h[k]*getfinitesignalatindex(x,i-k+1) for k = 1:length(h) )
        # i-k+1 instead of i-k because Julia is 1-indexing.
    end

    return out
end

# remove the padding due to linear convolution. There will still be a bit of
#   shift probably due to phase-shift of system. This reasoning is untested.
# g is the filter.
# y is the signal.
function forcesamesizeafterlinearconv(g::Vector{T},y::Vector{T})::Vector{T} where T

    offset1 = div(length(g)-1,2)
    offset2 = length(g)-1-offset1

    return y[offset1:end-offset2]
end

# in progress. optimized for speed. eventually package pynogram without allocating N_chs worth of output.
# ignore entries of output that require out-of-bound h.
function fastlinearconv(a::Vector{T}, b::Vector{T})::Vector{T} where T
    m = length(a)
    n = length(b)

    out = zeros(T, m+n-1)
    @inbounds for j = 1:m
        for k = 1:n
            out[j+k-1] += a[j]*b[k]
        end
    end
    return out
end

function circconv(  x::Vector{T},
                h::Vector{T}) where T

    out = Vector{T}(undef, length(x)+length(h))
    for i = 1:length(x)+length(h)

        out[i] = sum( h[k]*x[mod(i-k,length(x))+1] for k = 1:length(h) )
    end

    return out
end
