
"""
    logsumexp(x::Vector{Complex{T}})::Complex{T}

Log-sum-exp with choosing the shift to be the maximum absolute value for
the real and imaginary parts of x's components.

# Example:
```
a = randn(3) + randn(3) .* im
u = exp.(a)
RHS = logsumexp(u)
LHS = log(sum(exp.(u)))
println("RHS = ", RHS)
println("LHS = ", LHS)
```
"""
function logsumexp(x::Vector{Complex{T}})::Complex{T} where T <: Real

    # use the maximum magnitude for real and imaginary parts to get the shift, a.
    val_unused, 𝑖_real = findmax( real.(x) )
    val_unused, 𝑖_imag = findmax( imag.(x) )

    a::Complex{T} = real(x[𝑖_real]) + im*imag(𝑖_imag)

    # log-sum-exp.
    running_sum::Complex{T} = zero(Complex{T})
    for i = 1:length(x)
        running_sum += exp( x[i] - a )
    end

    return a + log(running_sum)
end
