mutable struct CholeskyType{T}
    c::Vector{T} # diagonal entries.
    a::Vector{T} # offdiagonal entries.
end
