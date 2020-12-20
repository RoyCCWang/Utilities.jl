### routines related to geometric transforms of a set of coordinates.


function shiftthenstretch(x,n::Int,T,N)
    return (x-n*T)/N
end

function shiftthencompress(x,n::Int,T,N)
    return (x-n*T)*N
end

# evaluates f(x/M)
function evalstretchedsignal(x, f::Function, M::Int)
    return f(x/M)
end

"""
    convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
"""
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function convertcompactdomain(x::Vector{T}, a::Vector{T}, b::Vector{T}, c::Vector{T}, d::Vector{T})::Vector{T} where T <: Real

    return collect( convertcompactdomain(x[i], a[i], b[i], c[i], d[i]) for i = 1:length(x) )
end
