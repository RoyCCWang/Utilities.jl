# for 2D visualizations.


function visualizemeshgridpcolor(   x_ranges::Vector{LinRange{T}},
                                    Y::Matrix{T},
                                    marker_locations::Vector,
                                    marker_symbol::String,
                                    fig_num::Int,
                                    title_string::String,
                                    x1_title_string::String = "Dimension 1",
                                    x2_title_string::String = "Dimension 2") where T <: Real
    #
    @assert length(x_ranges) == 2
    x_coords = collect( collect(x_ranges[d]) for d = 1:2 )

    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.pcolormesh(x_coords[2], x_coords[1], Y, cmap = "Greens_r")
    PyPlot.xlabel(x2_title_string)
    PyPlot.ylabel(x1_title_string)
    PyPlot.title(title_string)

    for i = 1:length(marker_locations)
        pt = reverse(marker_locations[i])
        PyPlot.annotate(marker_symbol, xy=pt, xycoords="data")
    end

    PyPlot.plt.colorbar()

    return fig_num
end

function plot2Dhistogram(fig_num::Int,
                        X::Vector{Vector{T}},
                        n_bins::Int,
                        limit_a::Vector{T},
                        limit_b::Vector{T},
                        use_bounds::Bool,
                        title_string::String,
                        colour_code::String = "Greys",
                        use_color_bar::Bool = false)::Int where T <: Real

    PyPlot.figure(fig_num)
    fig_num += 1

    N_viz = length(X)
    p1 = collect(X[n][2] for n = 1:N_viz)
    p2 = collect(X[n][1] for n = 1:N_viz)

    bounds = [[limit_a[2], limit_b[2]], [limit_a[1], limit_b[1]]]

    if use_bounds
        PyPlot.plt.hist2d(p1, p2, n_bins, range = bounds, cmap=colour_code)
    else
        PyPlot.plt.hist2d(p1, p2, n_bins, cmap=colour_code)
    end
    PyPlot.plt.axis("equal")

    if use_color_bar
        PyPlot.plt.colorbar()
    end

    PyPlot.title(title_string)

    return fig_num
end
