### routines for visualizing DSP-related objects.

function visualizsignals(h_array, fig_num::Int, title_strings)::Int

    for j = 1:length(h_array)

        PyPlot.figure(fig_num)
        fig_num += 1
        PyPlot.plot(h_array[j])
        PyPlot.title(title_strings[j])
    end

    return fig_num
end

function visualizefbmagrsp(h_array, fig_num::Int, title_string_first_part = "adjusted h's spectrum")::Int

    for j = 1:length(h_array)
        title_string = Printf.@sprintf("%s, ch = %d", title_string_first_part, j)
        fig_num = plotmagnitudersp(h_array[j], fig_num, title_string)
    end

    return fig_num
end

function plotphasersp(h::Vector{Float64}, fig_num::Int, title_string::String = "Phase response")

    ω_set_fft, DFT_evals, ω_set, DTFT_evals = getfreqrsp(h)

    # visualize.
    phase_rsp_G = angle.(DTFT_evals)
    phase_rsp_fft = angle.(DFT_evals)


    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(ω_set_fft, phase_rsp_fft, ".", label = "DFT")
    PyPlot.plot(ω_set, phase_rsp_G, label = "DTFT")

    PyPlot.title(title_string)
    PyPlot.legend()

    return fig_num
end

function getfreqrsp(h::Vector{Float64}, resolution_multiple::Int = 20)
    N_samples = length(h)

    ω_set_fft = collect( LinRange(0,2*π-2*π/N_samples,N_samples))
    DFT_evals = fft(h)

    ω_set = collect( LinRange(0,2*π,N_samples*resolution_multiple) )
    DTFT_evals = collect( computeDTFTviaformula(h,ω_set[i]) for i = 1:length(ω_set) )

    return ω_set_fft, DFT_evals, ω_set, DTFT_evals
end

function plotmagnitudersp(h::Vector{Float64}, fig_num::Int, title_string::String = "Magnitude response")

    ω_set_fft, DFT_evals, ω_set, DTFT_evals = getfreqrsp(h)

    # visualize.
    mag_rsp_G = abs.(DTFT_evals)
    mag_rsp_fft = abs.(DFT_evals)


    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(ω_set_fft, mag_rsp_fft, ".", label = "DFT")
    PyPlot.plot(ω_set, mag_rsp_G, label = "DTFT")

    PyPlot.title(title_string)
    PyPlot.legend()

    return fig_num
end
