

function generaterandomGMM(  mixture_weights::Vector{T},
                                D::Int) where T <: Real
    #
    N_components = length(mixture_weights)

    m_gen_array = Vector{Vector{Float64}}(undef,N_components)
    Σ_gen_array = Vector{Matrix{Float64}}(undef,N_components)
    for i = 1:N_components
        m_gen_array[i] = 3*i .*ones(T,D)

        S = randn(T,D,D)
        #S = Matrix{Float64}(LinearAlgebra.I,D,D)
        Σ_gen_array[i] = S'*S
        Σ_gen_array[i] = Σ_gen_array[i]./maximum(Σ_gen_array[i])
    end

    return getGMMdist(m_gen_array, Σ_gen_array, mixture_weights)
end

function getGMMdist(    μ_array::Vector{Vector{T}},
                        Σ_array::Vector{Matrix{T}},
                        mixture_weights::Vector{T}) where T <: Real
    #
    N_components = length(mixture_weights)

    dist_gen_array = collect( Distributions.MvNormal(μ_array[i], Σ_array[i]) for i = 1:N_components )
    dist_gen = Distributions.MixtureModel(dist_gen_array, mixture_weights)

    pdf = xx->sum( mixture_weights[n]*Distributions.pdf(dist_gen_array[n], xx) for n = 1:N_components )

    return pdf, dist_gen
end
