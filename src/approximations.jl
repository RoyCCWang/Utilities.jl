### routines for approximating functions.

function smoothmin(x::Vector{T})::T where T <: Real
    return smoothextrema(x,convert(T,-999.0))
end

function smoothmax(x::Vector{T})::T where T <: Real
    return smoothextrema(x,convert(T,999.0))
end

function smoothextrema(x::Vector{T}, α::T) where T <: Real
    u = collect( log(x[i]) + x[i]*α for i = 1:length(x) )
    v = collect( x[i]*α for i = 1:length(x) )
    return exp( StatsFuns.logsumexp(u) - StatsFuns.logsumexp(v) )
    #return sum( x[i]*exp(x[i]*α) for i = 1:length(x))/sum( exp(x[i]*α) for i = 1:length(x))
end
# ## test code.
# α = 999.0
# D = 5
# x = randn(D)
# t = rand()
# y = t.*x
# smax_x = smoothextrema(x,α)
# smax_y = smoothextrema(y,α)
# println("smooth max of x is ", smax_x)
# println("max of x is ", maximum(x))
# println("smooth max of y is ", smax_y)
# println("max of y is ", maximum(y))
#
# smin_x = smoothextrema(x,-α)
# smin_y = smoothextrema(y,-α)
# println("smooth min of x is ", smin_x)
# println("min of x is ", minimum(x))
# println("smooth min of y is ", smin_y)
# println("min of y is ", minimum(y))
