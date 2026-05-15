module RVM
    using Statistics
    using LinearAlgebra
    using PyPlot # For plotting RVM results
    using Printf
    using CairoMakie, FileIO # For geometry_analysis chi2 maps

    # Use parent scope references to avoid redundant inclusions and recursion

    include("data.jl")
    include("plot.jl")
    include("tools.jl")
    include("heights.jl")

    # Helper function for extracting profile and PA (moved from Data.jl)
    function _extract_profile_and_pa(data, bin_st, bin_end, pulses_count)
        I_prof = vec(mean(data[:, bin_st:bin_end, 1], dims=1))
        Q_prof = vec(mean(data[:, bin_st:bin_end, 2], dims=1))
        U_prof = vec(mean(data[:, bin_st:bin_end, 3], dims=1))
        V_prof = vec(mean(data[:, bin_st:bin_end, 4], dims=1))
        L_prof = sqrt.(Q_prof.^2 .+ U_prof.^2)
        
        L_max = maximum(L_prof) # Used for SNR thresholding
        off_pulse_noise = std(data[:, 1:max(1, bin_st-1), 1])
        sigma_avg = off_pulse_noise / sqrt(pulses_count)

        pa_prof = [(L_prof[i] > 5.0*sigma_avg && L_prof[i] >= 0.3 * L_max) ? 0.5 * atan(U_prof[i], Q_prof[i]) * (180.0/pi) : NaN for i in 1:length(L_prof)]
        pa_err_prof = [(L_prof[i] > 5.0*sigma_avg && L_prof[i] >= 0.3 * L_max) ? 0.5 * sigma_avg / L_prof[i] * (180.0/pi) : NaN for i in 1:length(L_prof)]

        lon_prof = collect(range(-180, 180, length=length(I_prof)))

        return I_prof, L_prof, V_prof, pa_prof, pa_err_prof, lon_prof
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
        offset = 0.0
        for i in 2:length(pa_deg)
            if !isnan(pa_out[i]) && !isnan(pa_out[i-1])
                d = (pa_out[i] + offset) - pa_out[i-1]
                if abs(d) > tolerance
                    step = 90.0 * sign(d)
                    offset -= step
                end
            end
            if !isnan(pa_out[i])
                pa_out[i] += offset
            end
        end
        pa_out = mod.(pa_out .+ 90, 180) .- 90
        return pa_out
    end

    """
    fit_rvm(lon_deg, pa_deg, pa_err_deg; alpha_grid, beta_grid, phi0_grid)

    Fit the Rotating Vector Model (RVM) to position angle data.

    Grid search over α (inclination) and β (impact parameter), then for each
    (α, β) pair solves linearly for PA0 and searches over φ0 (inflection longitude).

        Returns NamedTuple with fields: alpha, beta, phi0, pa0, chi2, residuals.
    Returns nothing if fewer than 5 valid PA points.
    """
    function fit_rvm(lon_deg, pa_deg, pa_err_deg;
                     alpha_grid=0:1:90, beta_grid=-10:0.5:10, phi0_grid=0:1:359)
        mask = .!(isnan.(pa_deg) .| isnan.(pa_err_deg) .| (pa_err_deg .== 0.0))
        if count(mask) < 5
            return nothing
        end
        lon = lon_deg[mask]
        pa = pa_deg[mask]
        err = pa_err_deg[mask]
        w = 1 ./ err.^2

        # Adjust phi0_grid to be centered around the data if it's default
        if phi0_grid == 0:1:359 && minimum(lon) < 0
            phi0_grid = round(Int, minimum(lon)):1:round(Int, maximum(lon))
        end

        best = (chi2=Inf, alpha=NaN, beta=NaN, phi0=NaN, pa0=NaN)

        for alpha in alpha_grid
            for beta in beta_grid
                for phi0 in phi0_grid
                    # Calculate RVM curve for current alpha, beta, phi0
                    model = rvm_ppa(alpha, beta, lon; phi0=phi0, deg=true)

                    # Correct way to find PA0 for circular data (180-degree period)
                    diffs = pa .- model
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

        # Calculate final residuals for the best fit
        final_model = rvm_ppa(best.alpha, best.beta, lon; phi0=best.phi0, deg=true) .+ best.pa0
        final_res = mod.(pa .- final_model .+ 90, 180) .- 90

        return (; alpha=best.alpha, beta=best.beta, phi0=best.phi0, pa0=best.pa0, chi2=best.chi2, residuals=final_res)
    end

    """
    rvm_curve(params, lon_min_deg, lon_max_deg; npts=500)

    Evaluate RVM curve over a dense longitude grid (degrees in, degrees out).
    Returns (lon_dense, pa_rvm, pa_rvm_ortho) where pa_rvm_ortho is the 90° mode.
    """
    function rvm_curve(params, lon_min_deg, lon_max_deg; npts=500)
        lon = range(lon_min_deg, lon_max_deg, length=npts)
        pa_deg = rvm_ppa(params.alpha, params.beta, lon; phi0=params.phi0, deg=true) .+ params.pa0
        pa_deg = mod.(pa_deg .+ 90, 180) .- 90 # Wrap to (-90, 90]
        pa_ortho_deg = mod.(pa_deg .+ 90, 180) .- 90 # Orthogonal mode

        # Introduce NaNs for large jumps to prevent connecting them in plots
        for i in 2:length(pa_deg)
            if abs(pa_deg[i] - pa_deg[i-1]) > 90
                pa_deg[i-1] = NaN
            end
            if abs(pa_ortho_deg[i] - pa_ortho_deg[i-1]) > 90
                pa_ortho_deg[i-1] = NaN
            end
        end
        return lon, pa_deg, pa_ortho_deg
    end

    """
    _rvm_chi2_plane(lon_deg, pa_deg, pa_err_deg;
                         alpha_grid=0:2:180, beta_grid=-15:0.5:15, phi0_grid=0:2:359)

    Calculates the chi2 plane for RVM fitting.
    """
    function _rvm_chi2_plane(lon_deg, pa_deg, pa_err_deg;
                             alpha_grid=0:2:180, beta_grid=-15:0.5:15, phi0_grid=0:2:359)
        mask = .!(isnan.(pa_deg) .| isnan.(pa_err_deg) .| (pa_err_deg .== 0.0))
        if count(mask) < 5
            return alpha_grid, beta_grid, fill(Inf, length(alpha_grid), length(beta_grid))
        end
        lon = lon_deg[mask]
        pa = pa_deg[mask]
        err = pa_err_deg[mask]
        w = 1 ./ err.^2

        # Adjust phi0_grid to be centered around the data if it's default
        if phi0_grid == 0:2:359 && minimum(lon) < 0
            phi0_grid = round(Int, minimum(lon)):2:round(Int, maximum(lon))
        end

        chi2_plane = fill(Inf, length(alpha_grid), length(beta_grid))

        for (i, alpha) in enumerate(alpha_grid)
            for (j, beta) in enumerate(beta_grid)
                best_chi2_for_ab = Inf
                for phi0 in phi0_grid
                    # Calculate RVM curve for current alpha, beta, phi0
                    model = rvm_ppa(alpha, beta, lon; phi0=phi0, deg=true)

                    # Find optimal pa0
                    diffs = pa .- model
                    avg_rad = atan(sum(w .* sin.(deg2rad.(2 .* diffs))),
                                   sum(w .* cos.(deg2rad.(2 .* diffs))))
                    pa0 = rad2deg(avg_rad) / 2.0

                    # Calculate Chi2
                    res = (pa .- model .- pa0)
                    res_wrapped = mod.(res .+ 90, 180) .- 90
                    chi2 = sum((res_wrapped ./ err).^2)

                    if chi2 < best_chi2_for_ab
                        best_chi2_for_ab = chi2
                    end
                end
                chi2_plane[i, j] = best_chi2_for_ab
            end
        end
        return alpha_grid, beta_grid, chi2_plane
    end

    """
        calculate_rvm(X::Matrix{Float64}, t::Vector{Float64}; kernel_gamma::Float64=0.1, max_iter::Int=500)

    Implementacja algorytmu Relevance Vector Machine (uproszczony algorytm Tippinga).
    Służy do wyznaczenia wektorów relewantnych dla danych (alfa, beta).
    """
    function calculate_rvm(X::Matrix{Float64}, t::Vector{Float64}; kernel_gamma::Float64=0.1, max_iter::Int=500, tol::Float64=1e-5)
        N, D = size(X)

        # Macierz projektowa (Kernel Gaussowski)
        Phi = [exp(-kernel_gamma * sum(abs2, X[i,:] - X[j,:])) for i in 1:N, j in 1:N]

        alphas = fill(1e-3, N)
        beta_noise = 1.0 / (var(t) * 0.1)

        mu = zeros(N)
        active = collect(1:N)

        for iter in 1:max_iter
            Phi_a = Phi[:, active]
            A = Diagonal(alphas[active])

            # Obliczanie kowariancji i średniej wag (posterior)
            # Sigma = inv(beta_noise * Phi_a' * Phi_a + A)
            # Wykorzystujemy Cholesky'ego dla stabilności
            H = beta_noise * (Phi_a' * Phi_a) + A
            Sigma = inv(H)
            mu_a = beta_noise * Sigma * (Phi_a' * t)

            # Aktualizacja hiperparametrów gamma i alpha
            gamma = 1.0 .- alphas[active] .* diag(Sigma)
            alphas_new = gamma ./ (mu_a.^2)

            # Aktualizacja precyzji szumu
            res_sq = sum(abs2, t - Phi_a * mu_a)
            beta_noise = (N - sum(gamma)) / res_sq

            # Pruning (usuwanie nieistotnych wektorów)
            keep = findall(a -> a < 1e8, alphas_new)

            if length(keep) == length(active) && maximum(abs.(alphas_new .- alphas[active])) < tol
                active = active[keep]
                alphas[active] = alphas_new[keep]
                mu = mu_a[keep]
                break
            end

            active = active[keep]
            alphas[active] = alphas_new[keep]
            mu = mu_a[keep]
        end

        return (rv_indices=active, rv_coords=X[active, :], weights=mu, beta=beta_noise)
    end

    """
        calculate_pa_stats(pa_obs, pa_model; ndof_val=nothing)

    Oblicza błąd systematyczny, odchylenie standardowe i zredukowany chi2.
    """
    function calculate_pa_stats(pa_obs, pa_model; ndof_val=nothing)
        residuals = mod.(pa_obs .- pa_model .+ 90, 180) .- 90
        bias = mean(residuals)
        std_err = std(residuals)

        chi2 = sum(abs2, residuals ./ 5.0) # Zakładamy błąd pomiarowy 5 stopni jeśli nie podano

        ndof = isnothing(ndof_val) ? length(residuals) - 4 : ndof_val
        chi2_red = chi2 / ndof

        return (residuals=residuals, bias=bias, std_err=std_err, chi2_red=chi2_red, ndof=ndof)
    end

    """
    geometry_analysis(indir; chi2_red_max=10.0)

    Andrzej's approach: Processes low and high frequency data separately
    but sends them to a single plot for comparison.
    """
    function geometry_analysis(indir; chi2_red_max=10.0)
        p = Tools.read_params(joinpath(indir, "params.json"))

        # Load both frequency sets (assuming they are already converted to .txt by process_psrdata_16)
        l = Data.load_ascii_all(joinpath(indir, "pulsar_low.txt"))
        h = Data.load_ascii_all(joinpath(indir, "pulsar_high.txt"))

        bin_st, bin_end = p["bin_st"], p["bin_end"]
        pulses_l, bins_l = size(l, 1), size(l, 2)
        pulses_h, bins_h = size(h, 1), size(h, 2)

        # --- LOW FREQUENCY ANALYSIS BLOCK ---
        I_l, L_l, V_l, pa_l_raw, pa_err_l_raw, lon_l = _extract_profile_and_pa(l, bin_st, bin_end, pulses_l)
        pa_l = deopm_pa(pa_l_raw)

        res_l_pre = fit_rvm(lon_l, pa_l, pa_err_l_raw)
        scale_l = isnothing(res_l_pre) ? 1.0 : max(1.0, sqrt(res_l_pre.chi2 / (length(res_l_pre.residuals) - 4))) # Calculate reduced chi2
        res_l = fit_rvm(lon_l, pa_l, pa_err_l_raw .* scale_l)
        alphas_l, betas_l, chi2_l = _rvm_chi2_plane(lon_l, pa_l, pa_err_l_raw .* scale_l)

        # RVM Relevance Vector Machine
        rvm_res_l = calculate_rvm(hcat(lon_l, pa_l), pa_l; kernel_gamma=0.05)

        # --- HIGH FREQUENCY ANALYSIS BLOCK ---
        I_h, L_h, V_h, pa_h_raw, pa_err_h_raw, lon_h = _extract_profile_and_pa(h, bin_st, bin_end, pulses_h)
        pa_h = deopm_pa(pa_h_raw)

        res_h_pre = fit_rvm(lon_h, pa_h, pa_err_h_raw)
        scale_h = isnothing(res_h_pre) ? 1.0 : max(1.0, sqrt(res_h_pre.chi2 / (length(res_h_pre.residuals) - 4))) # Calculate reduced chi2
        res_h = fit_rvm(lon_h, pa_h, pa_err_h_raw .* scale_h)
        alphas_h, betas_h, chi2_h = _rvm_chi2_plane(lon_h, pa_h, pa_err_h_raw .* scale_h)

        # RVM Relevance Vector Machine
        rvm_res_h = calculate_rvm(hcat(lon_h, pa_h), pa_h; kernel_gamma=0.05)

        # Plotting
        Plot.geometry(chi2_l, chi2_h, alphas_l, betas_l, indir;
                      ndof_l=(isnothing(res_l) ? 1 : length(res_l.residuals) - 4),
                      ndof_h=(isnothing(res_h) ? 1 : length(res_h.residuals) - 4),
                      chi2_red_max=chi2_red_max, P_sec=get(p, "period", nothing),
                      fit_result_l=res_l, fit_result_h=res_h)

        if !isnothing(res_l)
            plot_rvm_results(lon_l, pa_l_raw, res_l, alphas_l, betas_l, chi2_l, rvm_res_l;
                                 name_mod="Low_RVM", outdir=indir)
        end

        if !isnothing(res_h)
            plot_rvm_results(lon_h, pa_h_raw, res_h, alphas_h, betas_h, chi2_h, rvm_res_h;
                                 name_mod="High_RVM", outdir=indir)
        end

        return (low_freq_result=res_l, high_freq_result=res_h)
    end

    """
        plot_rvm_results(lon, pa, fit_result, alpha_grid, beta_grid, chi2_plane, rvm_res; name_mod="RVM", outdir=".")

    Tworzy zaawansowany wykres RVM z elipsami chi2, wektorami relewantnymi i histogramem reszt.
    """
    function plot_rvm_results(lon, pa, fit_result, alpha_grid, beta_grid, chi2_plane, rvm_res; name_mod="RVM", outdir=".")
        PyPlot.figure(figsize=(12, 8))
        PyPlot.clf()

        # 1. Mapa konturowa Chi2 z elipsami i Relevance Vectors
        PyPlot.subplot2grid((2, 3), (0, 0), colspan=2)
        cp = PyPlot.contourf(beta_grid, alpha_grid, chi2_plane .- minimum(chi2_plane), levels=20, cmap="viridis")
        PyPlot.colorbar(cp, label=raw"$\Delta \chi^2$")

        # Elipsy ufności (dla 2 stopni swobody: 50%, 90%, 95%)
        # Progi: 1.39, 4.61, 5.99
        PyPlot.contour(beta_grid, alpha_grid, chi2_plane .- minimum(chi2_plane), levels=[1.39, 4.61, 5.99],
                colors=["white", "yellow", "red"], linestyles=["--", "-", "-"])

        # Relevance Vectors
        rv_idx = rvm_res.rv_indices
        PyPlot.scatter(rvm_res.rv_coords[:, 2], rvm_res.rv_coords[:, 1], c="cyan", marker="*", s=100, label="Relevance Vectors", zorder=5)

        # Punkt najlepszego dopasowania
        PyPlot.scatter([fit_result.beta], [fit_result.alpha], c="magenta", marker="x", s=100, label="Best Fit")

        PyPlot.xlabel(raw"Impact parameter $\beta$ ($^\circ$)")
        PyPlot.ylabel(raw"Magnetic inclination $\alpha$ ($^\circ$)")
        PyPlot.title(raw"$\chi^2$ Surface & Relevance Vectors")
        PyPlot.legend(fontsize=8)

        # 2. Histogram reszt i krzywa teoretyczna
        PyPlot.subplot2grid((2, 3), (0, 2))
        res = fit_result.residuals
        n, bins, _ = PyPlot.hist(res, bins=15, density=true, color="gray", alpha=0.7, label="Residuals")

        # Teoretyczna krzywa Gaussa (dla porwnania z chi2)
        mu_res, std_res = mean(res), std(res)
        x_theory = range(minimum(bins), maximum(bins), length=100)
        PyPlot.plot(x_theory, [1/(std_res * sqrt(2*pi)) * exp(-0.5*((x-mu_res)/std_res)^2) for x in x_theory],
             "r-", lw=2, label="Gaussian Fit")

        PyPlot.title("PA Residuals Distribution")
        PyPlot.xlabel("Resid. (deg)")
        PyPlot.legend(fontsize=8)

        # 3. Dopasowanie RVM do danych PA
        PyPlot.subplot2grid((2, 3), (1, 0), colspan=3)
        PyPlot.scatter(lon, pa, c="black", s=10, alpha=0.5)

        # Model curve
        lon_dense = range(minimum(lon)-10, maximum(lon)+10, length=500)
        alpha_r, beta_r, phi0_r = deg2rad(fit_result.alpha), deg2rad(fit_result.beta), deg2rad(fit_result.phi0)
        zeta_r = alpha_r + beta_r
        phi_dense_r = deg2rad.(lon_dense)
        pa_model = rad2deg.(atan.(sin(alpha_r) .* sin.(phi_dense_r .- phi0_r),
                           sin(zeta_r) .* cos(alpha_r) .- cos(zeta_r) .* sin(alpha_r) .* cos.(phi_dense_r .- phi0_r))) .+ fit_result.pa0

        PyPlot.plot(lon_dense, pa_model, "r-", lw=2, label="RVM Model")
        PyPlot.axvline(fit_result.phi0, color="blue", ls="--", alpha=0.5, label=raw"$\phi_0$")

        PyPlot.xlabel("Longitude (deg)")
        PyPlot.ylabel("PA (deg)")
        PyPlot.legend(fontsize=8)

        PyPlot.tight_layout()
        PyPlot.savefig(joinpath(outdir, "$(name_mod)_analysis.pdf"))
        println("RVM Plot saved to: $(name_mod)_analysis.pdf")
    end

end # module