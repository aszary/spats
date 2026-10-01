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

export sliding_lrfs, p3_track, contrast_null, good_windows, window_length, p3_segments,
       p3_groups, merge_sections

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


"""
    window_length(p3; ncycles=4, lmin=16) -> Int

LRFS length for a pulsar with drift periodicity `p3`: max(lmin, ncycles·P3).
With fmin = 2/L and the edge guard of 1/L, P3 ≲ L/3 is measurable, so
ncycles = 4 leaves a margin for P3 wandering upwards.
"""
window_length(p3::Real; ncycles::Real=4, lmin::Int=16) = max(lmin, round(Int, ncycles * p3))


"""
    p3_segments(tr, good, nps; jump=0.25, maxgap=nothing, minwin=nothing, extend=0.5)
        -> Vector{NamedTuple}

Split the P3 track into continuous sections. Consecutive good windows i < j
(stride 1) belong to the same section when

    j − i ≤ maxgap   (default L÷2: a few rejected windows do not split it)
    |f3[j] − f3[i]| ≤ jump / L

With stride 1 the windows overlap by L−1 pulses, so a genuine slow change of
P3 — monotonic or not — moves f3 by ≪ 1/L between neighbours and stays in one
section; a jump of the feature (mode change, peak switching to another
feature or harmonic) is a discontinuity of order 1/L. The section therefore
needs no model of how P3 is allowed to vary; that is left to the folding.
Sections with fewer than `minwin` good windows (default max(3, L÷4)) are
dropped.

Pulse range: the window centres of the first and last window ("core"),
widened by `extend`·L on both sides — a window centred at c already carries
the modulation of pulses c ± L/2, so without widening every section would
lose half a window at each end (J0034-0721: ~100-P bursts read as 17–40-P
sections). Widening is clipped to the observation (1…`nps`) and, where two
sections would overlap, to the midpoint between their cores.

Fields per section: first, last (pulses), npulse, core_first, core_last,
win (window indices), p3_med, p3_first, p3_last, p3_min, p3_max,
dp3 (linear slope of P3 [P/P]), f3_med, contrast_med.
"""
function p3_segments(tr, good, nps::Int; jump::Real=0.25, maxgap=nothing, minwin=nothing,
                     extend::Real=0.5)
    runs = _runs(tr, good; jump=jump, maxgap=maxgap, minwin=minwin)
    return _sections(tr, runs, nps; extend=extend)
end

function _runs(tr, good; jump=0.25, maxgap=nothing, minwin=nothing)
    L = tr.window
    mg = maxgap === nothing ? max(1, L ÷ 2) : maxgap
    mw = minwin === nothing ? max(3, L ÷ 4) : minwin
    runs = Vector{Vector{Int}}()
    for j in findall(good)
        if !isempty(runs)
            i = runs[end][end]
            if j - i ≤ mg && abs(tr.f3[j] - tr.f3[i]) ≤ jump / L
                push!(runs[end], j)
                continue
            end
        end
        push!(runs, [j])
    end
    filter!(w -> length(w) ≥ mw, runs)
    return sort!(runs, by=w -> tr.centers[w[1]])
end

function _sections(tr, runs, nps; extend=0.5)
    L = tr.window
    cf = [tr.centers[w[1]] for w in runs]
    cl = [tr.centers[w[end]] for w in runs]
    ext = extend * L
    segs = NamedTuple[]
    for (k, w) in enumerate(runs)
        lo = cf[k] - ext
        hi = cl[k] + ext
        k > 1 && (lo = max(lo, (cl[k-1] + cf[k]) / 2))
        k < length(runs) && (hi = min(hi, (cl[k] + cf[k+1]) / 2))
        first = clamp(ceil(Int, lo), 1, nps)
        last = clamp(floor(Int, hi), 1, nps)
        # midpoint shared by two sections: give it to the earlier one
        k > 1 && !isempty(segs) && first ≤ segs[end].last && (first = segs[end].last + 1)
        p = tr.p3[w]
        c = tr.centers[w]
        dp3 = length(w) > 2 ? cov(c, p) / var(c) : 0.0
        push!(segs, (first=first, last=last, npulse=last - first + 1,
                     core_first=round(Int, cf[k]), core_last=round(Int, cl[k]), win=w,
                     p3_med=median(p), p3_first=p[1], p3_last=p[end],
                     p3_min=minimum(p), p3_max=maximum(p), dp3=dp3,
                     f3_med=median(tr.f3[w]), contrast_med=median(tr.contrast[w])))
    end
    return segs
end


"""
    p3_groups(segs, L; tol=1.0) -> Vector{Int}

Group label for every section: sections sorted by median f3 are merged while
neighbours differ by ≤ tol/L in frequency. tol = 1 is a resolution criterion:
the Hann main lobe has FWHM ≈ 1.44/L, so sections closer than ~1/L cannot be
told apart by this window and the difference is estimator scatter (J1825+0004,
L = 57: a 33-pulse section at P3 = 11.7 next to 14.5, Δf = 0.95/L, split off
with tol = 0.5). Sections of one regime separated by nulls or bad windows land in
one group and can be folded together; distinct regimes (J0034-0721 modes A/B/C:
f3 ≈ 0.08/0.15/0.25) stay apart. Labels are ordered by decreasing total
number of pulses (group 1 = dominant regime).
"""
function p3_groups(segs, L::Int; tol::Real=1.0)
    n = length(segs)
    n == 0 && return Int[]
    order = sortperm([s.f3_med for s in segs])
    raw = zeros(Int, n)
    g = 1
    raw[order[1]] = g
    for k in 2:n
        a, b = segs[order[k-1]], segs[order[k]]
        b.f3_med - a.f3_med > tol / L && (g += 1)
        raw[order[k]] = g
    end
    tot = [sum(segs[i].npulse for i in 1:n if raw[i] == gg) for gg in 1:g]
    rank = invperm(sortperm(tot, rev=true))
    return [rank[r] for r in raw]
end


"""
    merge_sections(tr, segs, groups, nps; extend=0.5) -> (segs, groups)

Merge sections that touch (next.first == prev.last + 1, i.e. split only at a
shared midpoint) and belong to the same group. With short windows the
estimator scatter of f3 between neighbouring windows can exceed `jump`/L
(J0820-1350, L = 19: 20 sections, all P3 ≈ 4.8); such splits are noise, not
regime changes, once grouping has put both sides together.
"""
function merge_sections(tr, segs, groups, nps::Int; extend::Real=0.5)
    isempty(segs) && return segs, groups
    order = sortperm([s.first for s in segs])
    runs = Vector{Vector{Int}}()
    grp = Int[]
    prev = 0
    for i in order
        if prev != 0 && groups[i] == groups[prev] && segs[i].first == segs[prev].last + 1
            append!(runs[end], segs[i].win)
        else
            push!(runs, copy(segs[i].win))
            push!(grp, groups[i])
        end
        prev = i
    end
    return _sections(tr, runs, nps; extend=extend), grp
end

end # module P3Track
