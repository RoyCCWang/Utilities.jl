
# to do: fix this markdown.
"""
    setupcubicitp(ϕ::Array{T,D}, x_ranges::Vector{LinRange{T}}, amplification_factor::T)


Separable cubic spline interpolation on compact interval, with gain parameter.


Inputs: `x_ranges` is a vector of LinRange that specifies the sampling coordinates.
        `Φ is the set of sampling values associated with x_ranges.
        `amplification_factor` is the gain multiplier applied to the interpolated function.


Returns:    `f` is the interpolated function.
            `df` is the gradient of `f`.
            `d2f` is the Hessian of `f`. It throws an error if evaluated at any x in the compact interval that is specified by `x_ranges`.
"""
function setupcubicitp(ϕ::Array{T,D}, x_ranges::Vector{LinRange{T}}, amplification_factor::T) where {T,D}

    @assert D == length(x_ranges)
    N_array = collect( length(x_ranges[d]) for d = 1:D )

    itp_ϕ = Interpolations.interpolate(ϕ,
                Interpolations.BSpline(Interpolations.Cubic(
                    Interpolations.Flat(Interpolations.OnGrid()))))

    etp_ϕ = Interpolations.extrapolate(itp_ϕ, Interpolations.Line())
    #etp_ϕ = Interpolations.extrapolate(itp_ϕ, 0)

    st = collect( x_ranges[d][1] for d = 1:D )
    fin = collect( x_ranges[d][end] for d = 1:D )
    f = xx->etp_ϕ(interval2itpindex(xx,
                st,
                fin,
                N_array)...)*amplification_factor

    # chain rule for first derivatives.
    df = xx->( Interpolations.gradient(etp_ϕ,interval2itpindex(xx,
                st,
                fin,
                N_array)...) .* derivativeinterval2itpindex(st,fin,N_array) .*amplification_factor )

    # chain rule for second derivatives.
    d2f = xx->( Interpolations.hessian(etp_ϕ,interval2itpindex(xx,
                st,
                fin,
                N_array)...) .* (derivativeinterval2itpindex(st,fin,N_array).^2) .*amplification_factor )


    return f, df, d2f
end

function interval2itpindex(x::T, a::T, b::T, N::Int)::T where T <: Real
    return (x-a)/(b-a)*(N-1) + 1
end

# converts an x ∈ [a[d],b[d]]^D, x ∈ ℝ^D, to an index i ∈ [N]^D.
function interval2itpindex(x::Vector{T}, a::Vector{T}, b::Vector{T}, N::Vector{Int})::Vector{T} where T <: Real
    return collect( interval2itpindex(x[d], a[d], b[d], N[d]) for d = 1:length(N) )
end

# do for each entry in X.
function interval2itpindex(X::Vector{Vector{T}}, a::Vector{T}, b::Vector{T}, N::Vector{Int})::Vector{Vector{T}} where T <: Real
    return collect( collect( interval2itpindex(X[n][d], a[d], b[d], N[d]) for d = 1:length(N) ) for n = 1:length(X) )
end

# for HCubature.
function interval2itpindex(x::AbstractArray{T,D}, a::Vector{T}, b::Vector{T}, N::Vector{Int})::Vector{T} where {T,D}
    return collect( interval2itpindex(x[d], a[d], b[d], N[d]) for d = 1:length(N) )
end
# function interval2itpindex(x::StaticArrays.SArray{Tuple{D},T,1,D}, a::Vector{T}, b::Vector{T}, N::Vector{Int})::Vector{T} where {T,D}
#     return collect( interval2itpindex(x[d], a[d], b[d], N[d]) for d = 1:length(N) )
# end

### for AD.
function interval2itpindex(x::T, a, b, N::Int) where T <: Real
    return (x-a)/(b-a)*(N-1) + 1
end

# converts an x ∈ [a[d],b[d]]^D, x ∈ ℝ^D, to an index i ∈ [N]^D.
function interval2itpindex(x::Vector{T}, a, b, N::Vector{Int}) where T <: Real
    return collect( interval2itpindex(x[d], a[d], b[d], N[d]) for d = 1:length(N) )
end

# do for each entry in X.
function interval2itpindex(X::Vector{Vector{T}}, a, b, N::Vector{Int}) where T <: Real
    return collect( collect( interval2itpindex(X[n][d], a[d], b[d], N[d]) for d = 1:length(N) ) for n = 1:length(X) )
end

function derivativeinterval2itpindex(a::T, b::T, N::Int)::T where T <: Real
    return (N-1)/(b-a)
end

function derivativeinterval2itpindex(a::Vector{T}, b::Vector{T}, N::Vector{Int})::Vector{T} where T <: Real
    return collect( derivativeinterval2itpindex(a[d], b[d], N[d]) for d = 1:length(N) )
end




# # WIP: what is with the hessian?
# function setupcubicitpmutator(ϕ::Array{T,D}, x_ranges::Vector{LinRange{T}}, amplification_factor::T) where {T,D}
#
#     @assert D == length(x_ranges)
#     N_array = collect( length(x_ranges[d]) for d = 1:D )
#
#     itp_ϕ = Interpolations.interpolate(ϕ,
#                 Interpolations.BSpline(Interpolations.Cubic(
#                     Interpolations.Flat(Interpolations.OnGrid()))))
#
#     etp_ϕ = Interpolations.extrapolate(itp_ϕ, Interpolations.Line())
#     #etp_ϕ = Interpolations.extrapolate(itp_ϕ, 0)
#
#     st = collect( x_ranges[d][1] for d = 1:D )
#     fin = collect( x_ranges[d][end] for d = 1:D )
#     f = xx->etp_ϕ(interval2itpindex(xx,
#                 st,
#                 fin,
#                 N_array)...)*amplification_factor
#
#     # chain rule for first derivatives.
#     df = xx->( Interpolations.gradient(etp_ϕ,interval2itpindex(xx,
#                 st,
#                 fin,
#                 N_array)...) .* derivativeinterval2itpindex(st,fin,N_array) .*amplification_factor )
#
#     # chain rule for second derivatives.
#     d2f = xx->( Interpolations.hessian(etp_ϕ,interval2itpindex(xx,
#                 st,
#                 fin,
#                 N_array)...) .* (derivativeinterval2itpindex(st,fin,N_array).^2) .*amplification_factor )
#
#
#     return f, df, d2f
# end
#
#
# """
# Chain rule for first derivatives. Implements:
#
# df = xx->( Interpolations.gradient(etp_ϕ,interval2itpindex(xx,
#             st,
#             fin,
#             N_array)...) .* derivativeinterval2itpindex(st,fin,N_array) .*amplification_factor )
# """
# function firstderivativechainrulewrapper!( g::Vector{T},
#                                             x,
#                                             N_array::Vector{Int},
#                                             etp_ϕ,
#                                             a::T)::Nothing where T
#     #
#     @assert length(x) == length(g)
#
#     Interpolations.gradient!(g, etp_ϕ, interval2itpindex(  x,
#                                                         st,
#                                                         fin,
#                                                         N_array)...)
#
#     for d = 1:D
#         g[d] = g[d] * derivativeinterval2itpindex(a[d], b[d], N[d]) * a
#     end
#
#     return nothing
# end
