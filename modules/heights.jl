module Heights
    using Glob
    using FITSIO
    using ProgressMeter
    using CairoMakie, FileIO
    using PDFIO
    using Statistics
    using FFTW
    using DelimitedFiles

    
    include("functions.jl")
    include("tools.jl")
    include("plot.jl")
    include("data.jl")
    include("GaussianFit.jl")



function pos_angle(indir; low_fit=Dict(), high_fit=Dict(), threshold=0.2, savepath=nothing, show_plot=true)
    params = isfile(joinpath(indir, "params.json")) ? Tools.read_params(joinpath(indir, "params.json")) : Dict()

    low_file = joinpath(indir, "pulsar_low.txt")
    high_file = joinpath(indir, "pulsar_high.txt")

    has_low = isfile(low_file)
    has_high = isfile(high_file)

    if !has_low && !has_high
        error("Neither pulsar_low.txt nor pulsar_high.txt exists in $indir")
    end

    function _to_float(value)
        if value isa String
            return parse(Float64, value)
        elseif value isa Number
            return float(value)
        else
            return nothing
        end
    end

    function _lookup(params, name, freq)
        keys = ["$(name)_$(freq)", "$(freq)_$(name)", "$(name)$(freq)", "$(freq)$(name)"]
        for k in keys
            if haskey(params, k)
                return _to_float(params[k])
            end
        end
        return nothing
    end

    function _rvm_params(params, freq, override)
        alpha = get(override, "alpha", _lookup(params, "alpha", freq))
        beta = get(override, "beta", _lookup(params, "beta", freq))
        phi0 = get(override, "phi0", _lookup(params, "phi0", freq))
        phi0 = isnothing(phi0) ? get(override, "inflection", _lookup(params, "inflection", freq)) : phi0
        return Dict("alpha" => alpha, "beta" => beta, "phi0" => phi0)
    end

    function _make_profile(data)
        pulses, bins, pols = size(data)
        I = mean(data[:, :, 1], dims=1) |> vec
        Q = mean(data[:, :, 2], dims=1) |> vec
        U = mean(data[:, :, 3], dims=1) |> vec
        V = mean(data[:, :, 4], dims=1) |> vec
        L = sqrt.(Q .^ 2 .+ U .^ 2)
        psi = rad2deg.(0.5 .* atan.(U, Q))
        return I, L, V, psi
    end

    function _bin_longitude(bins, bin_st, bin_end)
        center = isnothing(bin_st) || isnothing(bin_end) ? (bins + 1) / 2 : (bin_st + bin_end) / 2
        width = 360.0 / bins
        return ((1:bins) .- center) .* width
    end

    low_params = _rvm_params(params, "low", low_fit)
    high_params = _rvm_params(params, "high", high_fit)

    fig = Figure(resolution = (1200, 700))
    axes = Dict()
    col = 1

    for (label, file, fit_params) in [("Low", low_file, low_params), ("High", high_file, high_params)]
        if !isfile(file)
            continue
        end

        data = Data.load_ascii_all(file)
        I, L, V, psi = _make_profile(data)
        bins = size(data, 2)
        phi = _bin_longitude(bins, get(params, "bin_st", nothing), get(params, "bin_end", nothing))

        if col == 1
            ax_ppa = Axis(fig[1, 1]; ylabel = "PA [deg]", title = "$(label) frequency: PPA and RVM")
            ax_prof = Axis(fig[2, 1]; xlabel = "Longitude [deg]", ylabel = "Flux density")
        else
            ax_ppa = Axis(fig[1, 2]; ylabel = "PA [deg]", title = "$(label) frequency: PPA and RVM")
            ax_prof = Axis(fig[2, 2]; xlabel = "Longitude [deg]", ylabel = "Flux density")
        end

        mask = (I .> 0) .& (L ./ max.(I, 1e-12) .>= threshold)
        scatter!(ax_ppa, phi[mask], psi[mask]; color = :black, markersize=4)

        if !isnothing(fit_params["alpha"]) && !isnothing(fit_params["beta"]) && !isnothing(fit_params["phi0"])
            phi_fit = range(minimum(phi), stop=maximum(phi), length=600)
            psi_fit = rvm_ppa(fit_params["alpha"], fit_params["beta"], phi_fit; phi0=fit_params["phi0"], deg=true)
            lines!(ax_ppa, phi_fit, psi_fit; color = :orange, linewidth = 2)
            lines!(ax_ppa, phi_fit, psi_fit .+ 90; color = :orange, linewidth = 2, linestyle = :dash)
            vlines!(ax_ppa, [fit_params["phi0"]]; color = :blue, linewidth = 2)
        elseif !isnothing(fit_params["phi0"])
            vlines!(ax_ppa, [fit_params["phi0"]]; color = :blue, linewidth = 2)
        end

        lines!(ax_prof, phi, I; color = :black, label = "Stokes I")
        lines!(ax_prof, phi, L; color = :red, label = "Linear L")
        lines!(ax_prof, phi, V; color = :blue, label = "Stokes V")
        axislegend(ax_prof; position = :rt)

        ax_ppa.xticksvisible = false
        ax_ppa.xgridvisible = false
        ax_ppa.ygridvisible = false
        ax_prof.xgridvisible = false
        ax_prof.ygridvisible = false

        col += 1
    end

    fig.layoutgap = 15
    if !isnothing(savepath)
        save(savepath, fig)
    end
    if show_plot
        display(fig)
    end

    return fig
end

"""
    rvm_ppa(alpha, beta, phi; phi0=0.0, deg=true)

Compute the Rotating Vector Model (RVM) position angle curve.

# Arguments
- `alpha`: magnetic inclination angle (deg or rad)
- `beta`: impact parameter (deg or rad)
- `phi`: longitude array (deg or rad)
- `phi0`: inflection point longitude (deg or rad, default 0)
- `deg`: if true, all angles are in degrees (default true)

# Returns
- Array of position angles (deg)
"""
function rvm_ppa(alpha, beta, phi; phi0=0.0, deg=true)
    if deg
        α = deg2rad(alpha)
        β = deg2rad(beta)
        φ = deg2rad.(phi)
        φ₀ = deg2rad(phi0)
    else
        α = alpha
        β = beta
        φ = phi
        φ₀ = phi0
    end
    ζ = α + β
    # RVM formula
    pa = atan.(sin(α) .* sin.(φ .- φ₀),
               sin(ζ) .* cos(α) .- cos(ζ) .* sin(α) .* cos.(φ .- φ₀))
    if deg
        return rad2deg.(pa)
    else
        return pa
    end
end

end