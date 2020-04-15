

function plothistogram(fig_num::Int,
                        X::Vector{T},
                        n_bins::Int,
                        bounds::Vector{T},
                        use_bounds::Bool,
                        title_string::String,
                        normalize_counts_flag::Bool = false,
                        colour_code::String = "gray") where T <: Real
    #

    PyPlot.figure(fig_num)
    fig_num += 1

    if use_bounds
        PyPlot.plt.hist(X, n_bins, range = bounds, color = colour_code, density = normalize_counts_flag )
    else
        PyPlot.plt.hist(X, n_bins, color = colour_code, density = normalize_counts_flag )
    end
    PyPlot.title(title_string)

    return fig_num
end
