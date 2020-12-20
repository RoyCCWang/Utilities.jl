### routines that relate to discrete math.

# check if n and m are coprime.
function iscoprime(n::Integer, m::Integer)
    if n == 0 && m == 0
        return false
    end

    if gcd(n, m) > 1
        false
    else
        true
    end
end

# this is just gcd() in Julia.
function EuclideanAlgotihm(a::Int, b::Int, MAX_ITERS::Int = 10000)

    iter = 1
    while iter < MAX_ITERS
        if a < b
            a,b = b,a
        end

        r::Int = rem(a,b)
        if r == 0
            return b, true
        end

        a = b
        b = r

        iter += 1
    end

    return 0, false
end
