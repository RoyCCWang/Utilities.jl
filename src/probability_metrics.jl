# integral probability metrics or distances.

# returns a reparameterization of f such that it can be integrated over an Inf domain.
# assumes the infinite domains are in the first D dimensions.
function getInfdomainintegrand(f::Function, D::Int) where T <: Real

    return tt->f( collect(tt[d]/(1-tt[d]^2) for d = 1:D) )*
                    prod( (1+tt[d]^2)/(1-tt[d]^2)^2 for d = 1:D )
end

# evaluates 𝕂𝕃(p||q) via numerical integration.
function evalKLintegral(p::Function,
                        q::Function,
                        x_min::Vector{T},
                        x_max::Vector{T},
                        max_integral_evals::Int = 10000) where T <: Real

    f = xx->(p(xx)*(log(p(xx))-log(q(xx))))

    D = length(x_min)

    (val, err) = HCubature.hcubature(f, x_min, x_max;
                    norm = LinearAlgebra.norm, rtol = sqrt(eps(Float64)), atol = 0,
                    maxevals = max_integral_evals, initdiv = 1)
    return val, err
end

# evaluates 𝕂𝕃(p||q) via numerical integration.
# assumes infinite domain for all dimensions.
function evalKLintegral(p::Function,
                        q::Function,
                        D::Int,
                        max_integral_evals::Int = 10000)

    f = xx->(p(xx)*(log(p(xx))-log(q(xx))))

    h = getInfdomainintegrand(f, D)
    x_min = -ones(Float64, D)
    x_max = ones(Float64, D)

    (val, err) = HCubature.hcubature(h, x_min, x_max;
                    norm = LinearAlgebra.norm, rtol = sqrt(eps(Float64)), atol = 0,
                    maxevals = max_integral_evals, initdiv = 1)
    return val, err
end
