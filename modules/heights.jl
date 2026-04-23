module Heights
    using Glob
    using FITSIO
    using ProgressMeter
    using GLMakie, FileIO
    using PDFIO
    using Statistics
    using FFTW
    using DelimitedFiles
    using LinearAlgebra

    
    include("functions.jl")
    include("tools.jl")
    include("plot.jl")
    include("data.jl")
    include("GaussianFit.jl")



function pos_angle(indir; low_fit=Dict(), high_fit=Dict(), threshold=0.2, savepath=nothing, show_plot=true)
    return geometry_analysis(indir; threshold=threshold, savepath=savepath, show_plot=show_plot)
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


"""
    deopm_pa(pa_deg; tolerance=45.0)

Remove orthogonal polarization mode (OPM) jumps from a PA track by unwrapping 
along longitude with a 90° period. Walks bin-by-bin keeping a running "continuous" PA; 
whenever the next sample lies farther than `tolerance` degrees from it, shift the sample 
by ±90° (the OPM step) until it rejoins the track. The accumulated offset persists, 
so a whole OPM island is corrected end-to-end regardless of length. Finally the output 
is wrapped back to (-90°, 90°] — which also undoes any shifts of ±180° that correspond 
to the intrinsic mod-π PA wrap, leaving those pixels unchanged.
"""
function deopm_pa(pa_deg; tolerance=45.0)
    pa_out = copy(pa_deg)
    flipped = falses(length(pa_deg))
    offset = 0.0
    for i in 2:length(pa_deg)
        d = pa_out[i] + offset - pa_out[i-1]
        if abs(d) > tolerance
            step = 90.0 * sign(d)
            offset -= step
            flipped[i:end] .= .!flipped[i:end]
        end
        pa_out[i] += offset
    end
    pa_out = mod.(pa_out .+ 90, 180) .- 90
    return pa_out, flipped
end


"""
    fit_rvm(lon_deg, pa_deg, pa_err_deg; alpha_grid, beta_grid, phi0_grid)

Fit the Rotating Vector Model (RVM) to position angle data.

Grid search over α (inclination) and β (impact parameter), then for each
(α, β) pair solves linearly for PA0 and searches over φ0 (inflection longitude).

Returns NamedTuple with fields: alpha, beta, phi0, pa0 (all in degrees), chi2.
Returns nothing if fewer than 5 valid PA points.
"""
function fit_rvm(lon_deg, pa_deg, pa_err_deg;
                 alpha_grid=0:1:90, beta_grid=-10:0.5:10, phi0_grid=0:1:359)
    mask = .!(isnan.(pa_deg) .| isnan.(pa_err_deg))
    if count(mask) < 5
        return nothing
    end
    lon = lon_deg[mask]
    pa = pa_deg[mask]
    err = pa_err_deg[mask]
    w = 1 ./ err.^2

    # Adjust phi0_grid to be centered around the data if it's default
    if phi0_grid == 0:1:359 && minimum(lon) < 0
        phi0_grid = round(minimum(lon)):1:round(maximum(lon))
    end

    best = (chi2=Inf, alpha=NaN, beta=NaN, phi0=NaN, pa0=NaN)

    for alpha in alpha_grid
        for beta in beta_grid
            for phi0 in phi0_grid
                zeta = alpha + beta
                phi = deg2rad.(lon)
                phi0r = deg2rad(phi0)
                
                # RVM shape calculation
                num = sin(deg2rad(alpha)) .* sin.(phi .- phi0r)
                den = sin(deg2rad(zeta)) .* cos(deg2rad(alpha)) .- cos(deg2rad(zeta)) .* sin(deg2rad(alpha)) .* cos.(phi .- phi0r)
                model = rad2deg.(atan.(num, den))

                # Correct way to find PA0 for circular data (180-degree period)
                # We use the circular mean of the differences
                diffs = pa .- model
                # Transform to 2*theta to use standard circular mean (360-deg) for 180-deg data
                avg_rad = atan(sum(w .* sin.(deg2rad.(2 .* diffs))), 
                               sum(w .* cos.(deg2rad.(2 .* diffs))))
                pa0 = rad2deg(avg_rad) / 2.0

                # Calculate Chi2 using circular residuals (wrapped to [-90, 90])
                res = (pa .- model .- pa0)
                res_wrapped = mod.(res .+ 90, 180) .- 90
                chi2 = sum((res_wrapped ./ err).^2)

                if chi2 < best.chi2
                    best = (chi2=chi2, alpha=alpha, beta=beta, phi0=phi0, pa0=pa0)
                end
            end
        end
    end
    return (; alpha=best.alpha, beta=best.beta, phi0=best.phi0, pa0=best.pa0, chi2=best.chi2)
end


"""
    rvm_curve(params, lon_min_deg, lon_max_deg; npts=500)

Evaluate RVM curve over a dense longitude grid (degrees in, degrees out).
Returns (lon_dense, pa_rvm, pa_rvm_ortho) where pa_rvm_ortho is the 90° mode.
"""
function rvm_curve(params, lon_min_deg, lon_max_deg; npts=500)
    lon = range(lon_min_deg, lon_max_deg, length=npts)
    α = params.alpha
    β = params.beta
    φ0 = params.phi0
    pa0 = params.pa0
    ζ = α + β
    phi = deg2rad.(lon)
    phi0r = deg2rad(φ0)
    model = atan.(sin(deg2rad(α)) .* sin.(phi .- phi0r),
                  sin(deg2rad(ζ)) .* cos(deg2rad(α)) .-
                  cos(deg2rad(ζ)) .* sin(deg2rad(α)) .* cos.(phi .- phi0r))
    pa_deg = rad2deg.(model) .+ pa0
    pa_deg = mod.(pa_deg .+ 90, 180) .- 90
    pa_ortho_deg = mod.(pa_deg .+ 90 .+ 90, 180) .- 90
    return lon, pa_deg, pa_ortho_deg
end


"""
    pulse_width_fraction(profile, nbin; frac=0.1, on_st=nothing, on_end=nothing)

W_frac (e.g. W10): distance between the outermost longitude bins above
`frac * peak`, following Posselt et al. 2021 / Johnston et al. 2023.
If `on_st`/`on_end` are given, the peak and search are restricted to that
on-pulse window; this avoids noise spikes extending the width artificially.
"""
function pulse_width_fraction(profile, nbin; frac=0.1, on_st=nothing, on_end=nothing)
    if on_st !== nothing && on_end !== nothing
        prof = profile[on_st:on_end]
        offset = on_st - 1
    else
        prof = profile
        offset = 0
    end
    peak = maximum(prof)
    above = findall(x -> x > frac * peak, prof)
    if isempty(above)
        return 0.0
    end
    return (last(above) - first(above) + 1) * 360.0 / nbin
end


"""
    geometry_analysis(indir; threshold=0.2, savepath=nothing, show_plot=true)

Perform geometry analysis on pulsar data, including RVM fitting and α–β plane visualization.
"""

function _rvm_chi2_plane(lon_deg, pa_deg, pa_err_deg;
                         alpha_grid=0:2:180, beta_grid=-15:0.5:15, phi0_grid=0:2:359)
    lon = deg2rad.(lon_deg)
    err = pa_err_deg
    chi2_plane = fill(Inf, length(alpha_grid), length(beta_grid))

    for (i, alpha) in enumerate(alpha_grid)
        α = deg2rad(alpha)
        for (j, beta) in enumerate(beta_grid)
            ζ = α + deg2rad(beta)
            best_chi2 = Inf
            for phi0 in phi0_grid
                φ0 = deg2rad(phi0)
                num = sin(α) .* sin.(lon .- φ0)
                den = sin(ζ) .* cos(α) .- cos(ζ) .* sin(α) .* cos.(lon .- φ0)
                model = rad2deg.(atan.(num, den))
                diffs = pa_deg .- model
                avg_rad = atan(sum((1 ./ err.^2) .* sin.(deg2rad.(2 .* diffs))),
                               sum((1 ./ err.^2) .* cos.(deg2rad.(2 .* diffs))))
                pa0 = rad2deg(avg_rad) / 2.0
                res = pa_deg .- model .- pa0
                res_wrapped = mod.(res .+ 90, 180) .- 90
                chi2 = sum((res_wrapped ./ err).^2)
                if chi2 < best_chi2
                    best_chi2 = chi2
                end
            end
            chi2_plane[i, j] = best_chi2
        end
    end

    return alpha_grid, beta_grid, chi2_plane
end

function geometry_analysis(indir; threshold=0.2, savepath=nothing, show_plot=true)
    println("--- Starting Geometry Analysis (RVM Fitting) ---")
    params = isfile(joinpath(indir, "params.json")) ? Tools.read_params(joinpath(indir, "params.json")) : Dict()

    low_file = joinpath(indir, "pulsar_low.txt")
    high_file = joinpath(indir, "pulsar_high.txt")

    has_low = isfile(low_file)
    has_high = isfile(high_file)

    if !has_low && !has_high
        error("Neither pulsar_low.txt nor pulsar_high.txt exists in $indir")
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

    fig = Figure(resolution = (1400, 900))
    row = 1

    for (label, file) in [("Low", low_file), ("High", high_file)]
        if !isfile(file)
            continue
        end

        println("Processing $label frequency data...")
        data = Data.load_ascii_all(file)
        I, L, V, psi = _make_profile(data)
        bins = size(data, 2)
        phi = _bin_longitude(bins, get(params, "bin_st", nothing), get(params, "bin_end", nothing))

        mask = (I .> 0) .& (L ./ max.(I, 1e-12) .>= threshold)
        if sum(mask) < 5
            println("Not enough valid PA points for $label frequency after masking.")
            continue
        end

        pa_masked = psi[mask]
        lon_masked = phi[mask]
        pa_deopm, _ = deopm_pa(pa_masked)
        pa_err = fill(5.0, sum(mask))

        fit_result = fit_rvm(lon_masked, pa_deopm, pa_err;
                             alpha_grid=0:2:180, beta_grid=-15:0.5:15, phi0_grid=0:2:359)
        alpha_grid, beta_grid, chi2_plane = _rvm_chi2_plane(lon_masked, pa_deopm, pa_err;
                                                            alpha_grid=0:2:180, beta_grid=-15:0.5:15, phi0_grid=0:2:359)

        ax_ppa = Axis(fig[row, 1]; ylabel = "PA [deg]", title = "$(label) frequency: PA + RVM")
        ax_prof = Axis(fig[row+1, 1]; xlabel = "Longitude [deg]", ylabel = "Flux density")
        ax_ab = Axis(fig[row:(row+1), 2]; xlabel = "β [deg]", ylabel = "α [deg]", title = "χ² surface")

        scatter!(ax_ppa, lon_masked, pa_deopm; color = :black, markersize=4)
        if fit_result !== nothing
            lon_dense = range(minimum(phi), stop=maximum(phi), length=600)
            pa_rvm = rvm_ppa(fit_result.alpha, fit_result.beta, lon_dense; phi0=fit_result.phi0, deg=true) .+ fit_result.pa0
            pa_rvm = mod.(pa_rvm .+ 90, 180) .- 90
            lines!(ax_ppa, lon_dense, pa_rvm; color = :orange, linewidth = 2)
            lines!(ax_ppa, lon_dense, pa_rvm .+ 90; color = :orange, linewidth = 2, linestyle = :dash)
            vlines!(ax_ppa, [fit_result.phi0]; color = :blue, linewidth = 2)
            scatter!(ax_ab, [fit_result.beta], [fit_result.alpha]; color = :red, markersize=8)
            println("$(label): α=$(fit_result.alpha)°, β=$(fit_result.beta)°, φ0=$(fit_result.phi0)°, pa0=$(fit_result.pa0)°, χ²=$(fit_result.chi2)")
        end

        lines!(ax_prof, phi, I; color = :black, label = "Stokes I")
        lines!(ax_prof, phi, L; color = :red, label = "Linear L")
        lines!(ax_prof, phi, V; color = :blue, label = "Stokes V")
        axislegend(ax_prof; position = :rt)

        heatmap!(ax_ab, beta_grid, alpha_grid, chi2_plane; colormap = :viridis)
        Colorbar(fig[row:(row+1), 3], ax_ab, label = "χ²")

        ax_ppa.xticksvisible = false
        ax_ppa.xgridvisible = false
        ax_ppa.ygridvisible = false
        ax_prof.xgridvisible = false
        ax_prof.ygridvisible = false
        ax_ab.xgridvisible = true
        ax_ab.ygridvisible = true

        row += 2
    end

    fig.layoutgap = 15
    if !isnothing(savepath)
        save(savepath, fig)
        println("Plot successfully saved to: $savepath")
    end
    if show_plot
        screen = display(fig)
        if !isinteractive()
            wait(screen)
        end
    end

    println("--- Analysis Complete ---")
    return fig
end



"""
    position_angle(indir)

Wrapper function for position angle analysis (alias for geometry_analysis).
"""
function position_angle(indir; kwargs...)
    return geometry_analysis(indir; kwargs...)
end


end