"""
P3Track — sliding longitude-resolved fluctuation spectra and the P3(t) track.

Same construction as Fig. 4 of Szary et al. (2022, ApJ 934, 23): the LRFS is
computed for a window of `window` consecutive pulses, shifted by `stride`
pulses at a time, and summed over on-pulse longitudes. P3 in every window
comes from a Gaussian fitted around the strongest peak of that spectrum.

Differences from the paper, all needed for short windows (16-32 pulses):

  * the window is tapered (periodic Hann) before the FFT — without it the
    DC term and slow intensity changes leak over the whole 0-0.5 range when
    L is short;
  * the window is zero-padded to `pad`·L points, so the peak is sampled
    finely enough for the Gaussian fit (padding interpolates, it does not
    add resolution: the main lobe of a pure sinusoid still has
    FWHM ≈ 1.44/L in frequency);
  * an off-pulse strip of the same width goes through the same pipeline,
    giving the radiometer-noise level of the spectrum in every window.

Resolution limit to keep in mind: frequencies below ≈ 2/L sit inside the
leakage of the DC term (Hann main-lobe half-width), so only P3 ≲ L/2 is
measurable, and ΔP3/P3 ≈ 1.44·P3/L is the resolving width of the feature
(the fitted centre can be more precise than that when the S/N is high).
"""
module P3Track

using FFTW
using Statistics
using LsqFit
using Random

export sliding_lrfs, p3_track, contrast_null, good_windows

"periodic Hann taper without zero end points (no pulse fully discarded at L = 16)"
hann_taper(L::Int) = sin.(π .* ((0:L-1) .+ 0.5) ./ L) .^ 2


"""
    off_pulse_bins(nbin, bin_st, bin_end; margin=0.1) -> UnitRange or Vector

Contiguous off-pulse strip of the on-pulse width, as far from the on-pulse
window as possible (circular in phase). Falls back to every off-pulse bin
when the strip does not fit.
"""
function off_pulse_bins(nbin::Int, bin_st::Int, bin_end::Int; margin::Real=0.1)
    w = bin_end - bin_st + 1
    m = round(Int, margin * w)
    # bins excluded: on-pulse ± margin (circular)
    excl = Set(mod1(b, nbin) for b in (bin_st - m):(bin_end + m))
    # strip centred half a turn away from the on-pulse centre
    c = round(Int, (bin_st + bin_end) / 2) + nbin ÷ 2
    strip = [mod1(b, nbin) for b in (c - w ÷ 2):(c - w ÷ 2 + w - 1)]
    if !any(in(excl), strip)
        return strip
    end
    return [b for b in 1:nbin if !(b in excl)]
end


"""
    sliding_lrfs(data, bin_st, bin_end; window=32, stride=1, pad=8,
                 off_bins=:auto) -> NamedTuple

Sliding longitude-averaged fluctuation spectra.

For every window start s (s = 1, 1+stride, …, N−window+1):

    X      = data[s:s+L−1, on] − column means      (static profile of the window)
    F(f,φ) = FFT_n( taper(n) · X(n,φ) ), zero-padded to pad·L
    P(f)   = Σ_φ |F(f,φ)|²

and the same for an off-pulse strip, scaled to the on-pulse number of bins.
Power is divided by Σ taper² so that P does not depend on L for white noise.

Arguments:
  data     – single pulses (N_pulses × N_bins)
  bin_st, bin_end – on-pulse window (1-indexed)
  window   – LRFS length L [pulses]
  stride   – step between window starts [pulses]
  pad      – zero-padding factor (nfft = pad·L)
  off_bins – :auto (`off_pulse_bins`), a bin collection, or nothing

Fields:
  freq      – frequencies [cycles/P], 0 … 0.5 (length nfft÷2+1)
  starts    – first pulse of each window
  centers   – window centre (starts + (L−1)/2), x-axis aligned with the pulse stack
  power     – on-pulse spectra, nwin × nf
  power_off – off-pulse spectra (same units), nwin × nf, or nothing
  noise     – per-window mean off-pulse power over f ≥ fmin (radiometer level)
  noise_std – per-window std of the off-pulse power over the same range
  fmin      – 2/L, lower edge of the usable range (DC leakage below)
  window, stride, pad, nfft, on_bins, off_bins
"""
function sliding_lrfs(data::AbstractMatrix, bin_st::Int, bin_end::Int;
                      window::Int=32, stride::Int=1, pad::Int=8, off_bins=:auto)
    N, nbin = size(data)
    L = window
    L ≤ N || error("window ($L) longer than the observation ($N pulses)")
    on = bin_st:bin_end
    off = off_bins === :auto ? off_pulse_bins(nbin, bin_st, bin_end) : off_bins
    nfft = pad * L
    nf = nfft ÷ 2 + 1
    freq = collect(0:nf-1) ./ nfft
    fmin = 2 / L
    use = freq .>= fmin

    tap = hann_taper(L)
    norm = sum(abs2, tap)
    starts = collect(1:stride:(N - L + 1))
    nwin = length(starts)

    power = zeros(nwin, nf)
    power_off = off === nothing ? nothing : zeros(nwin, nf)
    noise = fill(NaN, nwin)
    noise_std = fill(NaN, nwin)

    buf_on = zeros(nfft, length(on))
    plan_on = plan_rfft(buf_on, 1)
    if off !== nothing
        buf_off = zeros(nfft, length(off))
        plan_off = plan_rfft(buf_off, 1)
        scale_off = length(on) / length(off)
    end

    for (i, s) in enumerate(starts)
        rows = s:s+L-1
        X = data[rows, on]
        X .-= mean(X, dims=1)
        fill!(buf_on, 0.0)
        buf_on[1:L, :] .= tap .* X
        F = plan_on * buf_on
        power[i, :] .= vec(sum(abs2, F, dims=2)) ./ norm

        if off !== nothing
            Y = data[rows, off]
            Y .-= mean(Y, dims=1)
            fill!(buf_off, 0.0)
            buf_off[1:L, :] .= tap .* Y
            G = plan_off * buf_off
            power_off[i, :] .= vec(sum(abs2, G, dims=2)) ./ norm .* scale_off
            noise[i] = mean(power_off[i, use])
            noise_std[i] = std(power_off[i, use])
        end
    end

    return (freq=freq, starts=starts, centers=starts .+ (L - 1) / 2,
            power=power, power_off=power_off, noise=noise, noise_std=noise_std,
            fmin=fmin, window=L, stride=stride, pad=pad, nfft=nfft,
            on_bins=on, off_bins=off)
end


gauss_model(f, p) = p[1] .* exp.(-(f .- p[2]) .^ 2 ./ (2 .* p[3] .^ 2)) .+ p[4]


"""
    feature_peak(P, freq, srange, lo, hi, guard) -> (k, edge, contrast)

Index of the highest local maximum of P strictly inside `srange`, whether it
counts as an edge peak (no interior maximum, or within `guard` of `lo` /
`guard`/2 of `hi`), and its contrast P[k] / median(P[srange]).
"""
function feature_peak(P, freq, srange, lo, hi, guard)
    kmax = 0
    for k in srange[2:end-1]
        if P[k] > P[k-1] && P[k] ≥ P[k+1] && (kmax == 0 || P[k] > P[kmax])
            kmax = k
        end
    end
    k = kmax == 0 ? srange[argmax(view(P, srange))] : kmax
    edge = kmax == 0 || freq[k] < lo + guard || freq[k] > hi - guard / 2
    return k, edge, P[k] / median(view(P, srange))
end


"""
    p3_track(sl; frange=nothing, halfwidth=nothing) -> NamedTuple

P3 in every window of `sliding_lrfs` output `sl`: the highest local maximum
of P(f) strictly inside `frange` (default (fmin, 0.5)) — so a red continuum
rising towards fmin (null edges, slow intensity changes) cannot hide a
feature — refined with a Gaussian + constant
fitted over peak ± `halfwidth` (default 1/L, half the Hann main lobe).
Errors are 1σ from the fit covariance; a fit that fails or puts the centre
outside the fitted range falls back to the peak sample with NaN error.

Fields:
  f3, f3_err   – fitted feature frequency and 1σ [cycles/P]
  p3, p3_err   – 1/f3 and σ_f/f3² [P]
  fwhm         – fitted FWHM of the feature [cycles/P]
  peak         – P at the peak sample
  snr_off      – (peak − noise)/noise_std, against the off-pulse spectrum:
                 is the modulation above radiometer noise at all
  contrast     – peak / median on-pulse P over the search range: does the
                 feature stand out of the on-pulse fluctuation continuum
                 (jitter, energy variations); ~1–2 for a flat spectrum
  fit_ok       – Gaussian fit converged with the centre inside the range
  edge         – peak within 1/L of fmin (or 1/2L of the upper limit):
                 DC leakage or a feature outside the range, not a P3
  centers, starts, window – copied from `sl`
"""
function p3_track(sl; frange=nothing, halfwidth=nothing)
    freq = sl.freq
    L = sl.window
    lo, hi = frange === nothing ? (sl.fmin, 0.5) : frange
    lo = max(lo, sl.fmin)
    hw = halfwidth === nothing ? 1 / L : halfwidth
    srange = findall(f -> lo ≤ f ≤ hi, freq)
    isempty(srange) && error("empty search range ($lo, $hi) for window $L")

    nwin = size(sl.power, 1)
    f3 = fill(NaN, nwin); f3_err = fill(NaN, nwin); fwhm = fill(NaN, nwin)
    peak = fill(NaN, nwin); snr_off = fill(NaN, nwin); contrast = fill(NaN, nwin)
    fit_ok = falses(nwin)
    edge = falses(nwin)
    # peaks within half a main lobe of the search edges are leakage, not features
    guard = 1 / L

    for i in 1:nwin
        P = @view sl.power[i, :]
        k, edge[i], contrast[i] = feature_peak(P, freq, srange, lo, hi, guard)
        fpk = freq[k]
        peak[i] = P[k]
        if sl.power_off !== nothing
            snr_off[i] = (P[k] - sl.noise[i]) / sl.noise_std[i]
        end
        idx = findall(f -> abs(f - fpk) ≤ hw && f ≥ sl.fmin, freq)
        f3[i] = fpk
        length(idx) < 5 && continue
        x = freq[idx]; y = collect(P[idx])
        c0 = minimum(y)
        p0 = [P[k] - c0, fpk, 0.6 / L, c0]
        try
            fit = curve_fit(gauss_model, x, y, p0)
            p = fit.param
            if fit.converged && x[1] ≤ p[2] ≤ x[end] && p[1] > 0
                se = stderror(fit)
                f3[i] = p[2]
                f3_err[i] = se[2]
                fwhm[i] = 2.3548 * abs(p[3])
                fit_ok[i] = true
            end
        catch
            # keep the peak-sample estimate
        end
    end

    p3 = 1 ./ f3
    p3_err = f3_err ./ f3 .^ 2
    return (f3=f3, f3_err=f3_err, p3=p3, p3_err=p3_err, fwhm=fwhm, peak=peak,
            snr_off=snr_off, contrast=contrast, fit_ok=fit_ok, edge=edge,
            centers=sl.centers, starts=sl.starts, window=L, frange=(lo, hi))
end


"""
    contrast_null(data, sl; nshuffle=40, step=nothing, q=0.99, seed=1,
                  frange=nothing) -> NamedTuple

Per-window threshold for `p3_track(sl).contrast` from *local* pulse-order
shuffles: the L pulses of a window are permuted and the contrast of the
strongest interior peak recomputed. Shuffling keeps every single pulse
(brightness, spikes, jitter, nulls) and destroys only the ordering, i.e. any
periodicity — so the threshold answers "could these very pulses, in random
order, produce a feature this strong?". It compares the feature with the
pulsar's own on-pulse fluctuation continuum and stays meaningful for very
bright pulsars, where the off-pulse S/N is huge everywhere.

Local, not global: a global shuffle spreads bright episodes (mode changes,
J1825+0004 after pulse ~700) over every window and inflates the threshold
in the faint parts (J1825+0004, L = 64: 99% level 3.79 global vs 2.38 from
pulses 1–700 only). Shuffles are computed for windows starting every `step`
pulses (default L÷8) and pooled over anchors within ±L/2 of each window of
`sl`, giving ~8·nshuffle null values per window. Shuffled nulls put more
on/off edges into a window than the real data has, so the threshold is
conservative across null boundaries.

Fields: threshold (per window of `sl`), anchors, null (nanchor × nshuffle), q
"""
function contrast_null(data::AbstractMatrix, sl; nshuffle::Int=40, step=nothing,
                       q::Real=0.99, seed::Int=1, frange=nothing)
    N = size(data, 1)
    L = sl.window
    st = step === nothing ? max(1, L ÷ 8) : step
    on = sl.on_bins
    freq = sl.freq
    lo, hi = frange === nothing ? (sl.fmin, 0.5) : frange
    lo = max(lo, sl.fmin)
    srange = findall(f -> lo ≤ f ≤ hi, freq)
    guard = 1 / L

    tap = hann_taper(L)
    norm = sum(abs2, tap)
    buf = zeros(sl.nfft, length(on))
    plan = plan_rfft(buf, 1)
    P = zeros(length(freq))

    rng = Random.MersenneTwister(seed)
    anchors = collect(1:st:(N - L + 1))
    null = zeros(length(anchors), nshuffle)
    for (ia, a) in enumerate(anchors)
        X = data[a:a+L-1, on]
        X .-= mean(X, dims=1)
        for j in 1:nshuffle
            fill!(buf, 0.0)
            buf[1:L, :] .= tap .* X[Random.randperm(rng, L), :]
            F = plan * buf
            P .= vec(sum(abs2, F, dims=2)) ./ norm
            null[ia, j] = feature_peak(P, freq, srange, lo, hi, guard)[3]
        end
    end

    threshold = similar(sl.starts, Float64)
    for (i, s) in enumerate(sl.starts)
        near = findall(a -> abs(a - s) ≤ L ÷ 2, anchors)
        threshold[i] = quantile(vec(null[near, :]), q)
    end
    return (threshold=threshold, anchors=anchors, null=null, q=q)
end


"""
    good_windows(tr, threshold) -> BitVector

Windows with a usable P3: fit converged, peak not at the edge of the search
range, contrast above the shuffle threshold (scalar or per-window vector
from `contrast_null`).
"""
good_windows(tr, threshold) = tr.fit_ok .& .!tr.edge .& (tr.contrast .>= threshold)

end # module P3Track
