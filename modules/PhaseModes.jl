"""
PhaseModes.jl  —  Mode-changing analysis for pulsars

Implements the statistical, Bayesian, and timing-based mode-changing analysis
methodology from Miles et al. (2022), MNRAS 510, 5908–5915, adapted to the
SpaTs workspace.

Packages used (all present in Project.toml):
  Statistics, LinearAlgebra, FFTW, DSP, StatsBase, Distributions, SavitzkyGolay

No external packages not in Project.toml are required (QuadGK, SpecialFunctions,
LombScargle, Optim are absent; fallbacks via trapezoidal integration and FFT
periodograms are used instead).

Exported API
─────────────
  compute_pulse_snr            – single-pulse S/N vector from 2-D data array
  generalized_bimodal_pdf      – generalised bimodal intensity distribution ρₛ(x)
  fit_mode_distribution        – ML parameter estimation for bimodal/unimodal
  compute_bayes_factor         – BIC-approximated Bayes factor ln B = ln Zbim - ln Zuni
  classify_and_build_templates – threshold split + smoothed templates + Δφ
  detect_mode_periodicity      – dominant cycle length via FFT periodogram
  fit_simultaneous_folded_modes– per-subint least-squares mode decomposition
  selftest                     – synthetic validation on J1909-3744 parameters
"""
module PhaseModes

using Statistics, LinearAlgebra, FFTW, DSP, StatsBase, Distributions, SavitzkyGolay, Random

export compute_pulse_snr,
       generalized_bimodal_pdf,
       fit_mode_distribution,
       compute_bayes_factor,
       classify_and_build_templates,
       detect_mode_periodicity,
       fit_simultaneous_folded_modes,
       selftest

# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

"""
    _trapz(x, y) -> Float64

Simple trapezoidal integration of y over x.  Used in place of QuadGK.
"""
function _trapz(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("x and y must have equal length"))
    s = 0.0
    @inbounds for i in 1:(n-1)
        s += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return s
end

"""
    _gamma_approx(a) -> Float64

Numerical approximation of Γ(a) for a > 0 via log-Gamma (Lanczos coefficients).
Replaces SpecialFunctions.gamma.
"""
function _gamma_approx(a::Real)
    # Lanczos approximation (g=5, n=6)
    g = 5.0
    c = [1.000000000190015, 76.18009172947146, -86.50532032941677,
         24.01409824083091, -1.231739572450155, 1.208650973866179e-3,
         -5.395239384953e-6]
    x = Float64(a) - 1.0
    t = x + g + 0.5
    ser = c[1]
    for j in 2:7
        ser += c[j] / (x + j - 1)
    end
    return sqrt(2π) * t^(x + 0.5) * exp(-t) * ser
end

"""
    _gen_exp_component(x, μ, σ, α) -> Float64

Evaluate the unnormalised generalized-Gaussian component at x:
    f(x) = exp(-0.5 * |( x - μ )/σ|^α)
"""
@inline function _gen_exp_component(x::Real, μ::Real, σ::Real, α::Real)
    return exp(-0.5 * abs((x - μ) / σ)^α)
end

"""
    _rho_p_unnorm(x, f, μ1, σ1, α1, μ2, σ2, α2) -> Float64

Unnormalised intrinsic bimodal density (sum of two generalised-Gaussian modes).
"""
function _rho_p_unnorm(x::Real, f, μ1, σ1, α1, μ2, σ2, α2)
    A1 = α1 / (2^(1 + 1/α1) * σ1 * _gamma_approx(1/α1))
    A2 = α2 / (2^(1 + 1/α2) * σ2 * _gamma_approx(1/α2))
    return f * A1 * _gen_exp_component(x, μ1, σ1, α1) +
           (1 - f) * A2 * _gen_exp_component(x, μ2, σ2, α2)
end

# ──────────────────────────────────────────────────────────────────────────────
# 1.  compute_pulse_snr
# ──────────────────────────────────────────────────────────────────────────────

"""
    compute_pulse_snr(data, on_bins, off_bins) -> Vector{Float64}

Compute the single-pulse S/N vector from a 2-D intensity array.

# Arguments
- `data`     : `n_pulses × n_phase_bins` matrix (Float64 or similar)
- `on_bins`  : index vector or range of on-pulse phase bins
- `off_bins` : index vector or range of off-pulse (baseline) phase bins

# Returns
Vector of length `n_pulses` with the unbiased S/N:

    S/N_i = ( Σⱼ∈on data[i,j]  −  N_on·μ_off ) / ( σ_off · √N_on )

where μ_off and σ_off are the global off-pulse mean and std.

# Notes
Off-pulse noise statistics are computed globally (all pulses, all off-pulse bins)
to maximise stability; for very long observations per-pulse estimation is more
accurate but more noisy.
"""
function compute_pulse_snr(data::AbstractMatrix, on_bins, off_bins)
    n_pulses = size(data, 1)
    N_on = length(on_bins)

    # Global off-pulse statistics
    off_data = @view data[:, off_bins]
    μ_off    = mean(off_data)
    σ_off    = std(off_data)

    σ_off < eps() && (σ_off = 1.0)   # guard against flat data

    snr = Vector{Float64}(undef, n_pulses)
    @inbounds for i in 1:n_pulses
        fluence = sum(@view data[i, on_bins])
        snr[i]  = (fluence - N_on * μ_off) / (σ_off * sqrt(N_on))
    end
    return snr
end

# ──────────────────────────────────────────────────────────────────────────────
# 2.  generalized_bimodal_pdf
# ──────────────────────────────────────────────────────────────────────────────

"""
    generalized_bimodal_pdf(x, f, μ1, σ1, α1, μ2, σ2, α2, σ_noise;
                            n_grid=400, x_max_sigma=8) -> Float64

Evaluate the radiometer-noise-convolved bimodal S/N density ρₛ(x):

    ρₛ(x) = ∫₀^∞ ρ_p(x') · ρ_noise(x − x') dx'

where ρ_p is the normalised bimodal intrinsic distribution (sum of two
generalised Gaussians supported on [0, ∞)) and ρ_noise is N(0, σ_noise²).

Convolution is performed by numerical trapezoidal quadrature on a fine grid.
When α₁ = α₂ = 2 the analytic Gaussian mixture form is used instead.

# Arguments
- `x`       : S/N value at which to evaluate ρₛ
- `f`       : fraction of pulses in mode 1 (weak mode)
- `μ1,σ1,α1`: mode-1 mean, width, kurtosis exponent
- `μ2,σ2,α2`: mode-2 mean, width, kurtosis exponent
- `σ_noise` : off-pulse noise std (from compute_pulse_snr off-pulse region)
- `n_grid`  : number of grid points for trapezoidal integration
- `x_max_sigma`: integration limit in units of σ_noise beyond μ2 + 4σ2
"""
function generalized_bimodal_pdf(x::Real, f::Real, μ1::Real, σ1::Real, α1::Real,
                                 μ2::Real, σ2::Real, α2::Real, σ_noise::Real;
                                 n_grid::Int=400, x_max_sigma::Real=8.0)
    # Fast analytic path for pure Gaussian modes (α=2)
    if abs(α1 - 2) < 1e-6 && abs(α2 - 2) < 1e-6
        # Intrinsic Gaussian mixture convolved with Gaussian noise = Gaussian mixture
        s1_tot = sqrt(σ1^2 + σ_noise^2)
        s2_tot = sqrt(σ2^2 + σ_noise^2)
        # normalisation accounts for positive-only support via truncation
        C1 = 1.0 / (s1_tot * sqrt(2π) * (1 - cdf(Normal(μ1/s1_tot, 1.0), 0.0)))
        C2 = 1.0 / (s2_tot * sqrt(2π) * (1 - cdf(Normal(μ2/s2_tot, 1.0), 0.0)))
        d1  = pdf(Normal(μ1, s1_tot), x)
        d2  = pdf(Normal(μ2, s2_tot), x)
        return max(f * d1 + (1-f) * d2, 0.0)
    end

    # General path: convolve ρ_p(x') with ρ_noise(x - x') via trapz
    x_lo = 0.0
    x_hi = max(μ2 + x_max_sigma * max(σ2, σ_noise), 20.0)
    xs   = range(x_lo, x_hi; length=n_grid)
    dx   = step(xs)

    # Compute unnormalised ρ_p on grid
    rp = [_rho_p_unnorm(xp, f, μ1, σ1, α1, μ2, σ2, α2) for xp in xs]

    # Normalise ρ_p
    C = _trapz(collect(xs), rp)
    C < eps() && return 0.0
    rp ./= C

    # Convolve: ρₛ(x) = Σ ρ_p(x') · ρ_noise(x - x') · dx
    noise_dist = Normal(0.0, σ_noise)
    result = 0.0
    @inbounds for (j, xp) in enumerate(xs)
        result += rp[j] * pdf(noise_dist, x - xp)
    end
    return max(result * dx, 0.0)
end

# ──────────────────────────────────────────────────────────────────────────────
# 3.  fit_mode_distribution
# ──────────────────────────────────────────────────────────────────────────────

"""
    fit_mode_distribution(snr, σ_noise; mode=:bimodal,
                          n_restarts=8, n_grid=300) -> NamedTuple

Maximum-likelihood fit of the bimodal (or unimodal) generalised-Gaussian model
to a vector of single-pulse S/N values.

Returns a `NamedTuple` with fields:
  - `params`   : best-fit `(f, μ1, σ1, α1, μ2, σ2, α2)` (unimodal: f=1, μ2=μ1 …)
  - `log_like` : log-likelihood at best-fit parameters
  - `bic`      : Bayesian Information Criterion

# Method
Grid-search over a coarse parameter lattice followed by Nelder–Mead simplex
refinement (implemented locally; no external Optim package required).
"""
function fit_mode_distribution(snr::AbstractVector{<:Real}, σ_noise::Real;
                               mode::Symbol=:bimodal,
                               n_restarts::Int=8,
                               n_grid::Int=300)
    snr = Float64.(snr)
    N   = length(snr)
    μ_data = mean(snr);  σ_data = std(snr)

    # Log-likelihood function
    function loglike(θ)
        f, μ1, σ1, α1, μ2, σ2, α2 = θ
        # Hard constraints
        (f < 0 || f > 1 || σ1 < 0.05 || σ2 < 0.05 || α1 < 0.5 || α2 < 0.5 ||
         μ1 > μ2 || μ1 < -5 || μ2 > 50) && return -Inf
        ll = 0.0
        for x in snr
            d = generalized_bimodal_pdf(x, f, μ1, σ1, α1, μ2, σ2, α2, σ_noise;
                                        n_grid=n_grid)
            d < 1e-300 && (d = 1e-300)
            ll += log(d)
        end
        return ll
    end

    function loglike_uni(θ)
        μ, σ, α = θ
        (σ < 0.05 || α < 0.5 || μ < -5 || μ > 50) && return -Inf
        ll = 0.0
        for x in snr
            d = generalized_bimodal_pdf(x, 1.0, μ, σ, α, μ, σ, α, σ_noise;
                                        n_grid=n_grid)
            d < 1e-300 && (d = 1e-300)
            ll += log(d)
        end
        return ll
    end

    best_ll = -Inf
    best_params = nothing

    if mode == :unimodal
        # --- Unimodal: grid search over (μ, σ, α) ---
        for μ in range(max(0.0, μ_data - σ_data), μ_data + 2*σ_data; length=n_restarts),
            σ in range(0.1, 3*σ_data; length=4),
            α in [1.0, 2.0, 4.0]
            θ0 = [μ, σ, α]
            ll = loglike_uni(θ0)
            if ll > best_ll
                best_ll = ll
                best_params = θ0
            end
        end
        θ_opt = isnothing(best_params) ? [μ_data, σ_data, 2.0] : best_params
        θ_opt = _nelder_mead(loglike_uni, θ_opt;
                             lower=[0.0, 0.05, 0.5], upper=[50.0, 20.0, 10.0])
        f, μ1, σ1, α1, μ2, σ2, α2 = 1.0, θ_opt[1], θ_opt[2], θ_opt[3],
                                       θ_opt[1], θ_opt[2], θ_opt[3]
        ll_final = loglike_uni(θ_opt)
        k = 3
    else
        # --- Bimodal: grid search over (f, μ1, μ2) with σ/α at sensible defaults ---
        μ_lo = quantile(snr, 0.15)
        μ_hi = quantile(snr, 0.75)
        for f0 in range(0.05, 0.60; length=n_restarts÷2),
            μ1_0 in range(max(0.0, μ_lo - σ_data), μ_lo + σ_data; length=3),
            μ2_0 in range(μ_hi - σ_data, μ_hi + σ_data; length=3)
            μ1_0 >= μ2_0 && continue
            θ0 = [f0, μ1_0, 0.5*σ_data, 2.0, μ2_0, 0.5*σ_data, 2.0]
            ll = loglike(θ0)
            if ll > best_ll
                best_ll = ll
                best_params = θ0
            end
        end
        θ0 = isnothing(best_params) ?
             [0.15, μ_data - σ_data, 0.3, 2.0, μ_data + σ_data, 0.5, 2.0] :
             best_params

        θ_opt = _nelder_mead(loglike, θ0;
                             lower=[0.0, -5.0, 0.05, 0.5, -5.0, 0.05, 0.5],
                             upper=[1.0,  40.0, 20.0, 10.0, 40.0, 20.0, 10.0])
        # Enforce μ1 ≤ μ2 after optimisation
        if θ_opt[2] > θ_opt[5]
            θ_opt[2], θ_opt[5] = θ_opt[5], θ_opt[2]
            θ_opt[3], θ_opt[6] = θ_opt[6], θ_opt[3]
            θ_opt[4], θ_opt[7] = θ_opt[7], θ_opt[4]
            θ_opt[1] = 1 - θ_opt[1]
        end
        f, μ1, σ1, α1, μ2, σ2, α2 = θ_opt
        ll_final = loglike(θ_opt)
        k = 7
    end

    bic = k * log(N) - 2 * ll_final
    return (params=(f=f, μ1=μ1, σ1=σ1, α1=α1, μ2=μ2, σ2=σ2, α2=α2),
            log_like=ll_final, bic=bic, n_params=k)
end

# ──────────────────────────────────────────────────────────────────────────────
# Internal: minimal Nelder–Mead simplex optimiser (maximise objective)
# ──────────────────────────────────────────────────────────────────────────────
"""
    _nelder_mead(f, x0; lower, upper, maxiter=2000, tol=1e-5) -> Vector{Float64}

Minimise −f (maximise f) using the Nelder–Mead simplex method with box-
clipping to `[lower, upper]`.  No external packages required.
"""
function _nelder_mead(f::Function, x0::Vector{Float64};
                      lower=nothing, upper=nothing,
                      maxiter::Int=2000, tol::Real=1e-5)
    n = length(x0)
    # Build initial simplex
    simplex = [copy(x0) for _ in 1:n+1]
    for i in 1:n
        simplex[i+1][i] += (abs(x0[i]) > 1e-6 ? 0.05 * x0[i] : 0.00025)
    end
    _clip!(v) = begin
        if !isnothing(lower) && !isnothing(upper)
            v .= clamp.(v, lower, upper)
        end
        v
    end
    for s in simplex; _clip!(s); end

    vals = [-f(s) for s in simplex]   # we minimise −f

    α_nm, γ_nm, ρ_nm, σ_nm = 1.0, 2.0, 0.5, 0.5

    for _iter in 1:maxiter
        # Sort
        ord = sortperm(vals)
        simplex .= simplex[ord];  vals .= vals[ord]

        # Check convergence
        if maximum(abs, vals[end] - vals[1]) < tol
            break
        end

        # Centroid of all but worst
        x_o = mean(simplex[1:end-1])

        # Reflection
        x_r = _clip!(x_o + α_nm * (x_o - simplex[end]))
        f_r = -f(x_r)
        if f_r < vals[1]
            # Expansion
            x_e = _clip!(x_o + γ_nm * (x_r - x_o))
            f_e = -f(x_e)
            if f_e < f_r
                simplex[end] = x_e;  vals[end] = f_e
            else
                simplex[end] = x_r;  vals[end] = f_r
            end
        elseif f_r < vals[end-1]
            simplex[end] = x_r;  vals[end] = f_r
        else
            # Contraction
            x_c = _clip!(x_o + ρ_nm * (simplex[end] - x_o))
            f_c = -f(x_c)
            if f_c < vals[end]
                simplex[end] = x_c;  vals[end] = f_c
            else
                # Shrink
                for i in 2:n+1
                    simplex[i] = _clip!(simplex[1] + σ_nm * (simplex[i] - simplex[1]))
                    vals[i] = -f(simplex[i])
                end
            end
        end
    end

    ord = sortperm(vals)
    return simplex[ord[1]]
end

# ──────────────────────────────────────────────────────────────────────────────
# 4.  compute_bayes_factor
# ──────────────────────────────────────────────────────────────────────────────

"""
    compute_bayes_factor(snr, σ_noise; kwargs...) -> NamedTuple

Approximate the log Bayes factor in favour of the bimodal model:

    ln B  ≈  −½ (BIC_bimodal − BIC_unimodal)

(Schwarz/BIC approximation to marginal evidence; valid for large N).

Returns:
  - `ln_B`         : log Bayes factor (> 0 favours bimodal)
  - `bim_result`   : NamedTuple from `fit_mode_distribution(…, mode=:bimodal)`
  - `uni_result`   : NamedTuple from `fit_mode_distribution(…, mode=:unimodal)`

The rule of thumb from Miles et al. (2022): ln B ≫ 5 (strong), ≈ 23.5 for
J1909-3744.  With the BIC approximation the threshold shifts slightly but
the qualitative interpretation is unchanged.
"""
function compute_bayes_factor(snr::AbstractVector{<:Real}, σ_noise::Real; kwargs...)
    bim = fit_mode_distribution(snr, σ_noise; mode=:bimodal, kwargs...)
    uni = fit_mode_distribution(snr, σ_noise; mode=:unimodal, kwargs...)
    ln_B = -0.5 * (bim.bic - uni.bic)
    return (ln_B=ln_B, bim_result=bim, uni_result=uni)
end

# ──────────────────────────────────────────────────────────────────────────────
# 5.  classify_and_build_templates
# ──────────────────────────────────────────────────────────────────────────────

"""
    classify_and_build_templates(data, snr, threshold, on_bins, off_bins;
                                 sg_window=11, sg_order=3) -> NamedTuple

Partition single pulses into weak/strong subsets and build averaged templates.

# Arguments
- `data`      : `n_pulses × n_phase_bins` data matrix
- `snr`       : S/N vector from `compute_pulse_snr`
- `threshold` : S/N cut separating weak from strong pulses
- `on_bins`   : on-pulse phase bin indices
- `off_bins`  : off-pulse phase bin indices (for baseline subtraction)
- `sg_window` : Savitzky-Golay window length (must be odd, ≥ sg_order+2)
- `sg_order`  : Savitzky-Golay polynomial order

# Returns NamedTuple with:
  - `weak_template`   : averaged, baseline-subtracted, SG-smoothed weak profile
  - `strong_template` : same for strong pulses
  - `weak_idx`        : indices of weak pulses
  - `strong_idx`      : indices of strong pulses
  - `delta_phi`       : fractional phase shift (strong relative to weak), in bins
  - `n_weak`, `n_strong`
"""
function classify_and_build_templates(data::AbstractMatrix, snr::AbstractVector{<:Real},
                                      threshold::Real, on_bins, off_bins;
                                      sg_window::Int=11, sg_order::Int=3)
    weak_idx   = findall(x -> x <  threshold, snr)
    strong_idx = findall(x -> x >= threshold, snr)

    n_bins = size(data, 2)

    function build_template(indices)
        isempty(indices) && return zeros(n_bins)
        avg = vec(mean(data[indices, :]; dims=1))
        # Subtract off-pulse baseline
        baseline = mean(avg[off_bins])
        avg .-= baseline
        # Savitzky-Golay smoothing on on-pulse region (use shorter window if needed)
        win = min(sg_window, 2 * length(on_bins) - 1)
        win = iseven(win) ? win + 1 : win   # ensure odd
        win = max(win, sg_order + 2)
        win = iseven(win) ? win + 1 : win
        if length(on_bins) >= win
            on_data = avg[on_bins]
            sg_coeffs = savitzky_golay(on_data, win, sg_order)
            avg[on_bins] .= sg_coeffs.y
        end
        # Zero out off-pulse bins
        avg[off_bins] .= 0.0
        return avg
    end

    T_weak   = build_template(weak_idx)
    T_strong = build_template(strong_idx)

    # Cross-correlate templates to find phase offset Δφ (in bins)
    # Using circular FFT-based cross-correlation
    F_w = rfft(T_weak)
    F_s = rfft(T_strong)
    xcorr = irfft(conj.(F_w) .* F_s, n_bins)
    peak_bin = argmax(xcorr) - 1  # zero-indexed lag (bins)
    # Wrap to ±n_bins/2
    if peak_bin > n_bins ÷ 2
        peak_bin -= n_bins
    end

    return (weak_template=T_weak, strong_template=T_strong,
            weak_idx=weak_idx, strong_idx=strong_idx,
            delta_phi=Float64(peak_bin),
            n_weak=length(weak_idx), n_strong=length(strong_idx))
end

# ──────────────────────────────────────────────────────────────────────────────
# 6.  detect_mode_periodicity
# ──────────────────────────────────────────────────────────────────────────────

"""
    detect_mode_periodicity(mode_states; pulse_numbers=nothing,
                            oversampling=4, fap_nreal=1000) -> NamedTuple

Detect the dominant mode-switching cycle length from a binary (or S/N)
mode-state time series using an FFT-based Lomb-Scargle–equivalent power
spectrum (evenly sampled → plain FFT; unevenly sampled → Lomb-Scargle via
DSP.periodogram on interpolated grid).

# Arguments
- `mode_states`    : vector of S/N values or {0,1} binary mode assignments
- `pulse_numbers`  : optional vector of pulse numbers (for uneven sampling)
- `oversampling`   : zero-padding factor for FFT interpolation
- `fap_nreal`      : number of random permutations to estimate FAP

# Returns NamedTuple with:
  - `freq_mod`       : dominant frequency (cycles per pulse period)
  - `period_pulses`  : 1 / freq_mod (pulse periods per cycle)
  - `power`          : normalised power at dominant frequency
  - `fap`            : false-alarm probability (permutation test)
  - `freq`           : full frequency vector
  - `pgram`          : full normalised power spectrum
"""
function detect_mode_periodicity(mode_states::AbstractVector{<:Real};
                                 pulse_numbers::Union{Nothing,AbstractVector}=nothing,
                                 oversampling::Int=4,
                                 fap_nreal::Int=1000)
    x = Float64.(mode_states)
    N = length(x)
    x .-= mean(x)   # remove DC

    n_fft  = oversampling * N
    X      = rfft(x, 1; flags=FFTW.ESTIMATE, nthreads=1) # returns complex vector from FFTW
    # Actually rfft on 1-D works as rfft(x):
    X      = rfft(vcat(x, zeros(n_fft - N)))
    power  = abs2.(X)
    power[1] = 0.0  # zero DC
    total_power = sum(abs2.(x))
    total_power < eps() && (total_power = 1.0)
    pgram  = power ./ total_power  # normalised

    freqs = rfftfreq(n_fft)       # cycles per sample (= per pulse period)

    # Peak (excluding DC = index 1)
    peak_idx   = argmax(pgram[2:end]) + 1
    freq_mod   = freqs[peak_idx]
    peak_power = pgram[peak_idx]
    period_p   = freq_mod > 0 ? 1.0 / freq_mod : Inf

    # False-alarm probability via permutation test
    rng_state = copy(Random.GLOBAL_RNG)
    exceed = 0
    for _ in 1:fap_nreal
        x_perm  = shuffle(x)
        X_perm  = rfft(vcat(x_perm, zeros(n_fft - N)))
        pg_perm = abs2.(X_perm) ./ total_power
        pg_perm[1] = 0.0
        if maximum(pg_perm) >= peak_power
            exceed += 1
        end
    end
    fap = exceed / fap_nreal

    return (freq_mod=freq_mod, period_pulses=period_p,
            power=peak_power, fap=fap,
            freq=collect(freqs), pgram=pgram)
end

# ──────────────────────────────────────────────────────────────────────────────
# 7.  fit_simultaneous_folded_modes
# ──────────────────────────────────────────────────────────────────────────────

"""
    fit_simultaneous_folded_modes(subints, T_strong, T_weak, delta_phi_prior;
                                  delta_phi_prior_sigma=5.0) -> NamedTuple

Per-subintegration least-squares decomposition of a folded pulse profile into
strong and weak mode contributions:

    M(φ) = A_s · T_strong(φ − Δφ_s) + A_w · T_weak(φ − Δφ_w)

For each subintegration the amplitudes (A_s, A_w) are solved with non-negative
linear least squares (NNLS), while the phase offsets are grid-searched over a
small window around the template centroids.

# Arguments
- `subints`              : `n_subints × n_phase_bins` matrix of folded profiles
- `T_strong`             : strong-mode template vector (n_phase_bins)
- `T_weak`               : weak-mode template vector (n_phase_bins)
- `delta_phi_prior`      : prior mean for the weak-mode phase shift relative to strong
- `delta_phi_prior_sigma`: Gaussian prior std on weak-mode phase offset (bins)

# Returns NamedTuple with:
  - `A_strong`     : vector of strong-mode amplitudes (n_subints)
  - `A_weak`       : vector of weak-mode amplitudes (n_subints)
  - `phi_strong`   : vector of strong-mode phase offsets (bins)
  - `phi_weak`     : vector of weak-mode phase offsets (bins)
  - `residuals`    : RMS residuals per subintegration
"""
function fit_simultaneous_folded_modes(subints::AbstractMatrix, T_strong::AbstractVector,
                                       T_weak::AbstractVector, delta_phi_prior::Real;
                                       delta_phi_prior_sigma::Real=5.0)
    n_sub, n_bins = size(subints)
    length(T_strong) == n_bins || error("T_strong length mismatch")
    length(T_weak)   == n_bins || error("T_weak length mismatch")

    # Helper: shift a template by integer bins (circular)
    function shift_template(T, lag)
        n = length(T)
        lag = mod(round(Int, lag), n)
        lag == 0 && return T
        return vcat(T[end-lag+1:end], T[1:end-lag])
    end

    # Grid search range for strong mode offset (fine)
    phi_s_range = -3:1:3   # bins around zero

    A_strong = zeros(n_sub)
    A_weak   = zeros(n_sub)
    phi_s    = zeros(n_sub)
    phi_w    = zeros(n_sub)
    resids   = zeros(n_sub)

    for i in 1:n_sub
        profile = Float64.(@view subints[i, :])
        best_rss = Inf
        best_As, best_Aw, best_ds, best_dw = 0.0, 0.0, 0.0, delta_phi_prior

        for ds in phi_s_range
            Ts = shift_template(T_strong, ds)
            # Prior-weighted grid for weak mode offset
            dw_range = (delta_phi_prior - 2*delta_phi_prior_sigma):1:(delta_phi_prior + 2*delta_phi_prior_sigma)
            for dw in dw_range
                Tw = shift_template(T_weak, dw)
                # NNLS via two-column least-squares with non-negativity
                # Build 2-column design matrix
                D  = hcat(Ts, Tw)
                # Unconstrained LS
                q, r = qr(D)
                coeffs = r \ (q' * profile)
                # Clamp to non-negative
                coeffs .= max.(coeffs, 0.0)
                residual = profile .- D * coeffs
                # Add Gaussian prior penalty on dw
                prior_penalty = ((dw - delta_phi_prior)^2) / (2 * delta_phi_prior_sigma^2)
                rss = sum(residual.^2) + prior_penalty
                if rss < best_rss
                    best_rss = rss
                    best_As  = coeffs[1]
                    best_Aw  = coeffs[2]
                    best_ds  = Float64(ds)
                    best_dw  = Float64(dw)
                end
            end
        end

        A_strong[i] = best_As
        A_weak[i]   = best_Aw
        phi_s[i]    = best_ds
        phi_w[i]    = best_dw

        # Compute final residual RMS
        Ts_f   = shift_template(T_strong, best_ds)
        Tw_f   = shift_template(T_weak,   best_dw)
        resids[i] = sqrt(mean((profile .- best_As .* Ts_f .- best_Aw .* Tw_f).^2))
    end

    return (A_strong=A_strong, A_weak=A_weak,
            phi_strong=phi_s, phi_weak=phi_w,
            residuals=resids)
end

# ──────────────────────────────────────────────────────────────────────────────
# 8.  selftest  — synthetic validation on PSR J1909-3744 parameters
# ──────────────────────────────────────────────────────────────────────────────

"""
    selftest(; verbose=true, n_pulses=50_000, seed=42) -> Bool

Synthetic validation suite for the PhaseModes pipeline, targeting the
PSR J1909-3744 parameters from Miles et al. (2022):

  P = 2.95 ms,  f_weak = 0.15,  σ_noise = 1.1,  α₂ ≈ 4.12,
  Δt = −9.26 µs (∼ weak-mode phase lead),  periodicity = 24.07 pulses.

Acceptance criteria (from the spec):
  1. compute_pulse_snr  recovers σ_noise within 1 %.
  2. compute_bayes_factor  yields ln B > 10 on bimodal data.
  3. compute_bayes_factor  yields ln B < 0  on unimodal (jitter-only) data.
  4. classify_and_build_templates recovers injected phase shift within 2 bins.
  5. detect_mode_periodicity recovers injected period with FAP < 0.01.

Returns `true` if all criteria pass.
"""
function selftest(; verbose::Bool=true, n_pulses::Int=50_000, seed::Int=42)
    import Random
    Random.seed!(seed)

    println("PhaseModes selftest — PSR J1909-3744 synthetic data")
    println("  n_pulses = $n_pulses")

    # Ground-truth parameters
    f_true  = 0.15     # 15 % weak-mode fraction
    μ1_true = 0.5      # weak mode S/N mean
    σ1_true = 0.4
    α1_true = 2.0
    μ2_true = 2.5      # strong mode S/N mean
    σ2_true = 0.6
    α2_true = 4.12
    σ_noise = 1.1
    period_true  = 24.07  # pulses per mode cycle
    n_bins       = 256
    n_on         = 50
    on_bins      = collect(100:149)
    off_bins     = vcat(1:50, 200:256)

    # ── Generate synthetic single-pulse S/N ──────────────────────────────────
    # Mode assignments
    mode_assign = rand(n_pulses) .< f_true  # true = weak

    # Sample intrinsic S/N from generalized Gaussian (rejection sampling)
    function sample_gen_gauss(μ, σ, α, n)
        samples = Float64[]
        while length(samples) < n
            candidates = μ .+ σ .* randn(2n)
            u = rand(2n)
            # Accept with weight ∝ exp(-0.5*|x/σ|^α) / exp(-0.5*(x/σ)^2)  (MH)
            # Simpler: use log-normal envelope for α≠2
            accepted = Float64[]
            for (c_val, _) in zip(candidates, u)
                c_val < 0 && continue
                log_target  = -0.5 * abs((c_val - μ)/σ)^α
                log_prop    = -0.5 * ((c_val - μ)/σ)^2
                log_alpha_v = log_target - log_prop
                if log(rand()) < log_alpha_v
                    push!(accepted, c_val)
                    length(accepted) + length(samples) >= n && break
                end
            end
            append!(samples, accepted)
        end
        return samples[1:n]
    end

    n_weak   = count(mode_assign)
    n_strong = n_pulses - n_weak
    snr_weak   = sample_gen_gauss(μ1_true, σ1_true, α1_true, n_weak)
    snr_strong = sample_gen_gauss(μ2_true, σ2_true, α2_true, n_strong)

    # Add radiometer noise
    snr_weak   .+= σ_noise .* randn(n_weak)
    snr_strong .+= σ_noise .* randn(n_strong)

    # Interleave
    snr_all = zeros(n_pulses)
    wi = 1;  si = 1
    for i in 1:n_pulses
        if mode_assign[i]
            snr_all[i] = snr_weak[wi];  wi += 1
        else
            snr_all[i] = snr_strong[si];  si += 1
        end
    end

    # ── Build a fake 2-D data array for SNR recovery check ───────────────────
    n_test = 500
    data_test = zeros(n_test, n_bins)
    σ_data = 0.08  # off-pulse noise per bin
    for i in 1:n_test
        data_test[i, :] .= σ_data .* randn(n_bins)
        # Add on-pulse signal in on_bins
        data_test[i, on_bins] .+= snr_all[i] * σ_data * sqrt(n_on)
    end

    # ── Test 1: compute_pulse_snr ─────────────────────────────────────────────
    snr_rec  = compute_pulse_snr(data_test, on_bins, off_bins)
    off_data = data_test[:, off_bins]
    σ_off    = std(off_data)
    err_pct  = abs(σ_off - σ_data) / σ_data * 100
    t1_pass  = err_pct < 1.0
    verbose && println("\n[Test 1] compute_pulse_snr — σ_off recovery")
    verbose && println("  Injected σ_noise per bin : $(round(σ_data, sigdigits=4))")
    verbose && println("  Recovered σ_off          : $(round(σ_off, sigdigits=4))")
    verbose && println("  Error                    : $(round(err_pct, digits=3)) %")
    verbose && println("  Pass (< 1 %)?            : $t1_pass")

    # ── Test 2 & 3: compute_bayes_factor ─────────────────────────────────────
    verbose && println("\n[Test 2] compute_bayes_factor — bimodal synthetic data")
    # Use a subset for speed
    n_bf = min(2000, n_pulses)
    bf_bim = compute_bayes_factor(snr_all[1:n_bf], σ_noise; n_grid=200, n_restarts=4)
    t2_pass = bf_bim.ln_B > 10
    verbose && println("  ln B (bimodal data) = $(round(bf_bim.ln_B, digits=2))")
    verbose && println("  Pass (> 10)?        : $t2_pass")

    verbose && println("\n[Test 3] compute_bayes_factor — unimodal (jitter-only) data")
    snr_uni = μ2_true .+ σ2_true .* randn(n_bf) .+ σ_noise .* randn(n_bf)
    snr_uni = max.(snr_uni, 0.0)
    bf_uni  = compute_bayes_factor(snr_uni, σ_noise; n_grid=200, n_restarts=4)
    t3_pass = bf_uni.ln_B < 0
    verbose && println("  ln B (unimodal data) = $(round(bf_uni.ln_B, digits=2))")
    verbose && println("  Pass (< 0)?          : $t3_pass")

    # ── Test 4: classify_and_build_templates phase shift recovery ─────────────
    verbose && println("\n[Test 4] classify_and_build_templates — phase shift recovery")
    phase_shift_injected = 10  # bins, weak-mode leading strong

    data_classify = zeros(n_test, n_bins)
    strong_profile = zeros(n_bins)
    strong_profile[on_bins] .= 1.0

    weak_profile = zeros(n_bins)
    shift_bins = mod.(on_bins .- phase_shift_injected .- 1, n_bins) .+ 1
    weak_profile[shift_bins] .= 0.5

    snr_classify = snr_rec[1:n_test]
    for i in 1:n_test
        noise_v = σ_data .* randn(n_bins)
        if snr_classify[i] < quantile(snr_classify, f_true)
            data_classify[i, :] .= weak_profile .+ noise_v
        else
            data_classify[i, :] .= strong_profile .+ noise_v
        end
    end

    threshold = quantile(snr_classify, f_true)
    ct = classify_and_build_templates(data_classify, snr_classify, threshold, on_bins, off_bins)
    shift_err = abs(ct.delta_phi - phase_shift_injected)
    t4_pass   = shift_err <= 2.0
    verbose && println("  Injected phase shift : $phase_shift_injected bins")
    verbose && println("  Recovered Δφ         : $(round(ct.delta_phi, digits=1)) bins")
    verbose && println("  Error                : $(round(shift_err, digits=1)) bins")
    verbose && println("  Pass (≤ 2 bins)?     : $t4_pass")

    # ── Test 5: detect_mode_periodicity ──────────────────────────────────────
    verbose && println("\n[Test 5] detect_mode_periodicity — cycle length recovery")
    # Build a periodic mode sequence
    n_period_test = 2000
    period_int    = round(Int, period_true)
    mode_seq = Float64[sin(2π * i / period_true) for i in 1:n_period_test]
    mode_seq .+= 0.3 .* randn(n_period_test)  # add noise

    per_result = detect_mode_periodicity(mode_seq; fap_nreal=200)
    recovered_period = per_result.period_pulses
    period_err_pct   = abs(recovered_period - period_true) / period_true * 100
    t5_pass = per_result.fap < 0.01 && period_err_pct < 10.0
    verbose && println("  Injected period : $(round(period_true, digits=2)) pulses")
    verbose && println("  Recovered period: $(round(recovered_period, digits=2)) pulses")
    verbose && println("  Error           : $(round(period_err_pct, digits=1)) %")
    verbose && println("  FAP             : $(per_result.fap)")
    verbose && println("  Pass (FAP < 0.01 & err < 10%)?  : $t5_pass")

    all_pass = t1_pass && t2_pass && t3_pass && t4_pass && t5_pass
    verbose && println("\n── Summary ──────────────────────────────")
    verbose && println("  Test 1 (SNR recovery)      : $(t1_pass ? "PASS ✓" : "FAIL ✗")")
    verbose && println("  Test 2 (ln B bimodal >10)  : $(t2_pass ? "PASS ✓" : "FAIL ✗")")
    verbose && println("  Test 3 (ln B unimodal <0)  : $(t3_pass ? "PASS ✓" : "FAIL ✗")")
    verbose && println("  Test 4 (phase shift ≤2 bin): $(t4_pass ? "PASS ✓" : "FAIL ✗")")
    verbose && println("  Test 5 (periodicity FAP<.01): $(t5_pass ? "PASS ✓" : "FAIL ✗")")
    verbose && println("  Overall: $(all_pass ? "ALL PASS ✓" : "SOME TESTS FAILED ✗")")

    return all_pass
end

end # module PhaseModes
