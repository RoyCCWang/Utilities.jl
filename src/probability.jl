# routines related to probability.

function getnormalizedensity(  f_tilde::Function,
                            x_min::Vector{T},
                            x_max::Vector{T},
                            max_integral_evals::Int = 10000 ) where T <: Real
    #
    D = length(x_min)

    # get normalizing constant.

    (Z, Z_err) = HCubature.hcubature(f_tilde, x_min, x_max;
                    norm = LinearAlgebra.norm, rtol = sqrt(eps(T)), atol = 0,
                    maxevals = max_integral_evals, initdiv = 1)

    # get normalized function.
    f = xx->(f_tilde(xx)/Z)

    return f, Z, Z_err
end

# infinite domain.
function getnormalizedensity(  f_tilde::Function,
                                D::Int,
                                max_integral_evals::Int = 10000 )
    #

    # get normalizing constant.
    h = getInfdomainintegrand(f_tilde, D)
    x_min = -ones(Float64, D)
    x_max = ones(Float64, D)

    (Z, Z_err) = HCubature.hcubature(h, x_min, x_max;
                    norm = LinearAlgebra.norm, rtol = sqrt(eps(Float64)), atol = 0,
                    maxevals = max_integral_evals, initdiv = 1)

    # get normalized function.
    f = xx->(f_tilde(xx)/Z)

    return f, Z, Z_err
end

function evalstdnormalcdf(x::T)::T where T
    return 0.5*( 1.0 + SpecialFunctions.erf(x/sqrt(2)) )
end

function evalstdnormalpdf(x::T)::T where T
    return exp(-0.5*x*x)/sqrt(2*π)
end


### Gaussian systems.

"""
    getMVNmarginalparams(   x2::Vector{T},
                            m1::Vector{T},
                            m2::Vector{T},
                            S11::Matrix{T},
                            S12::Matrix{T},
                            S22::Matrix{T})

S is covariance matrix.
Conditionals of joint Gaussian. Outputs m1_given_2, S1_given_2
"""
function getMVNmarginalparams(  x2::Vector{T},
                                m1::Vector{T},
                                m2::Vector{T},
                                S11::Matrix{T},
                                S12::Matrix{T},
                                S22::Matrix{T}) where T <: Real

   # GP posterior.
   S21 = S12'
   S1_given_2 = Utilities.forcesymmetric( S11 - S12*S22\S21 )
   m1_given_2 = m1 + S12*(S22\(x2-m2))

   return m1_given_2, S1_given_2
end
