### routines that deal with signal boundaries or windowing.


# frames in col-major-ordering.
# u contains the time stamps.
function extractframecenters(u, N_samples_per_frame, N_frames)

    fc = Vector{Float64}(undef,undef,N_frames)
    for j = 1:N_frames

        st = (j-1)*N_samples_per_frame +1
        fin = st+N_samples_per_frame - 1

        fc[j] = Statistics.mean(u[st:fin])
    end

    return fc
end

function extractframecenters(frames_time::Matrix{T})::Vector{T} where T

    N_frames = size(frames_time,2)

    fc = Vector{T}(undef,N_frames)
    for j = 1:N_frames
        fc[j] = Statistics.mean(frames_time[:,j])
    end

    return fc
end

# pad sequence y with zeros if length(y) < N
function padzeros(y::Vector{T},N)::Vector{T} where T
    if length(y) < N
        return [y; zeros(T,N-length(y))]
    end

    return y
end

function evalzeropaddedsignal(y::Vector{T}, n::Int)::T where T
    if 1 <= n <= length(y)
        return y[n]
    end

    return zero(T)
end
