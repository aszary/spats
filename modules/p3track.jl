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
using Printf
using PyPlot
using Distributions

export sliding_lrfs, p3_track, contrast_null, good_windows, window_length, p3_segments,
       p3_groups, merge_sections, harmonic_groups, fundamental_track, select_groups,
       phase_fold, constant_fold, analyse, long_p3_pass, plot_track, plot_folds,
       plot_summary, harmonic_test, harmonic_power,
       template_significance, nyquist_pass, plot_nyquist

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
counts as an edge peak (no interior maximum, or within `guard` of `lo` — the
DC leakage), and its contrast P[k] / median(P[crange]). `crange` (default
`srange`) is the whole usable range f ≥ fmin in `p3_track`/`contrast_null`,
also when the search is restricted by `frange`: the median of a narrow
low-frequency range sits on the red continuum and hides a long-P3 feature
(J1825+0004 after ~715 at L = 228 in the second pass: no window passed). There is no guard at
the upper end: a feature at P3 ≈ 2.1 (f ≈ 0.48, J1001-5939) is real and only
has to be an interior maximum. A window of zapped pulses (P ≡ 0) returns
edge = true and contrast 0 instead of 0/0 = NaN, which broke the shuffle
quantile (J1524-5706, J1843-0211: pulses zeroed in the archive itself).
"""
function feature_peak(P, freq, srange, lo, hi, guard; crange=srange)
    kmax = 0
    for k in srange[2:end-1]
        if P[k] > P[k-1] && P[k] ≥ P[k+1] && (kmax == 0 || P[k] > P[kmax])
            kmax = k
        end
    end
    k = kmax == 0 ? srange[argmax(view(P, srange))] : kmax
    edge = kmax == 0 || freq[k] < lo + guard
    m = median(view(P, crange))
    # a window of zapped (all-zero) pulses has P ≡ 0: no feature, not 0/0
    m > 0 || return k, true, 0.0
    return k, edge, P[k] / m
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
  contrast     – peak / median on-pulse P over f ≥ fmin: does the
                 feature stand out of the on-pulse fluctuation continuum
                 (jitter, energy variations); ~1–2 for a flat spectrum
  fit_ok       – Gaussian fit converged with the centre inside the range
  edge         – no interior local maximum, or the peak within 1/L of fmin
                 (DC leakage), not a P3
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
    crange = findall(f -> f ≥ sl.fmin, freq)

    nwin = size(sl.power, 1)
    f3 = fill(NaN, nwin); f3_err = fill(NaN, nwin); fwhm = fill(NaN, nwin)
    peak = fill(NaN, nwin); snr_off = fill(NaN, nwin); contrast = fill(NaN, nwin)
    fit_ok = falses(nwin)
    edge = falses(nwin)
    # peaks within half a main lobe of the search edges are leakage, not features
    guard = 1 / L

    for i in 1:nwin
        P = @view sl.power[i, :]
        k, edge[i], contrast[i] = feature_peak(P, freq, srange, lo, hi, guard; crange=crange)
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
    crange = findall(f -> f ≥ sl.fmin, freq)
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
            null[ia, j] = feature_peak(P, freq, srange, lo, hi, guard; crange=crange)[3]
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


"""
    harmonic_groups(segs, groups, L; tol=1.0, check=nothing)
        -> (groups, harm::BitVector, tests)

Recognise groups that are the second harmonic of another, larger group
(|f3_a − 2·f3_b| ≤ tol/L, medians weighted by pulses) and relabel them to
that group: when the fundamental weakens, the strongest interior peak can
switch to 2·f3 (J0151-0635, L = 58: 40 pulses at P3 = 7.46 next to 14.3),
but the regime is the same. `harm[i]` marks sections whose track must be
halved in frequency (`fundamental_track`). Labels are re-ranked by size.

A frequency ratio of 2 alone cannot tell a harmonic from a separate mode
whose P3 happens to be half the other one. `check(sel)` — sel = section
indices of candidate group a — decides; it returns a NamedTuple with
`verdict` ∈ (:harmonic, :separate, :inconclusive) (see `harmonic_test`):
harmonic → joined, separate → kept as its own group, inconclusive → label
0 (sections to be discarded by the caller). Without `check` every 2:1
candidate is joined. `tests` lists (a, b, f_a, f_b, result).
"""
function harmonic_groups(segs, groups, L::Int; tol::Real=1.0, check=nothing)
    g = copy(groups)
    harm = falses(length(segs))
    tests = NamedTuple[]
    isempty(segs) && return g, harm, tests
    function stats(gs)
        ids = sort(unique(gs))
        npl = Dict(k => sum(segs[i].npulse for i in eachindex(segs) if gs[i] == k) for k in ids)
        f = Dict(k => sum(segs[i].f3_med * segs[i].npulse for i in eachindex(segs) if gs[i] == k) /
                      npl[k] for k in ids)
        return ids, npl, f
    end
    ids, npl, f = stats(g)
    for a in ids, b in ids
        a == b && continue
        if npl[b] > npl[a] && abs(f[a] - 2 * f[b]) ≤ tol / L
            sel = findall(==(a), g)
            isempty(sel) && continue
            res = check === nothing ? (verdict=:harmonic, is_harm=true) : check(sel)
            push!(tests, (a=a, b=b, f_a=f[a], f_b=f[b], result=res))
            if res.verdict == :harmonic
                for i in sel
                    g[i] = b
                    harm[i] = true
                end
            elseif res.verdict == :inconclusive
                g[sel] .= 0
            end
        end
    end
    ids, npl, _ = stats(g)
    order = sort(filter(!=(0), ids), by=k -> -npl[k])
    rank = Dict(k => r for (r, k) in enumerate(order))
    rank[0] = 0
    return [rank[k] for k in g], harm, tests
end


"""
    harmonic_power(F, cnt; h=1) -> Float64

Amplitude of the h-th Fourier component of a fold along P3 phase, RMS over
longitude, relative to the peak of the mean profile: h = 1 is the part of
the modulation that repeats once per fold cycle, h = 2 twice.
"""
function harmonic_power(F, cnt; h::Int=1)
    nb = size(F, 1)
    w = cnt ./ sum(cnt)
    prof = vec(sum(F .* w, dims=1))
    e = cis.(-2π * h .* ((1:nb) .- 0.5) ./ nb)
    c = vec(sum((F .- prof') .* (w .* e), dims=1))
    return sqrt(mean(abs2, c)) / maximum(prof)
end


"""
    harmonic_test(data, sl, tr, segs, sel; nshuffle=20, seed=1, nbins=8) -> NamedTuple

Is candidate group `sel` (sections whose feature sits at ≈ 2·f3 of another
group) a second harmonic, or a separate mode with half the P3? Its pulses
are demodulated and folded at *half* their tracked frequency (`phase_fold`
on a track halved in those sections, `nbins` phase bins), and the fold is
split along P3 phase into Fourier components (`harmonic_power`):

  h = 1 (once per fold cycle) — the fundamental at f3/2. A harmonic regime
        still carries its weaker fundamental, phase-locked to the pattern;
        a separate mode has nothing there.
  h = 2 — the tracked feature itself, folded twice per cycle; strong in
        both cases, so the fold's overall depth cannot decide (first
        version of this test: a synthetic separate mode at P3 = 4 next to
        P3 = 8 read depth 0.127 against a 0.096 shuffle maximum).

Verdict: :harmonic when h1 > max(h1 of `nshuffle` pulse-order shuffles)
(p ≲ 1/(nshuffle+1)); otherwise :separate if the group spans at least
`mincycles` cycles of the fundamental, else :inconclusive — a short group
has no power to show its fundamental (J0151-0635: 40 pulses ≈ 2.8 cycles;
synthetic harmonic with 42 pulses at high noise also failed), so absence of
h1 is not evidence of a separate mode there.
Synthetic check (`~/claude/work/scripts/p3track_harmonic_test.jl`, P3 = 8
regime + 400 pulses of either its 2nd-harmonic-dominated version or a
separate P3 = 4 mode): separate → :separate at noise 0.6 and 1.2;
harmonic → :harmonic at 0.6 and 1.2 (the latter barely: 0.084 vs 0.080).
Fields: verdict, is_harm, h1, h1_null_max, h1_null, h2, ncycles.
"""
function harmonic_test(data::AbstractMatrix, sl, tr, segs, sel; nshuffle::Int=20, seed::Int=1,
                       nbins::Int=8, mincycles::Real=10)
    mask = falses(length(segs)); mask[sel] .= true
    trh = fundamental_track(tr, segs, mask)
    tmp = zeros(Int, length(segs)); tmp[sel] .= 1
    fo = phase_fold(data, sl, trh, segs, tmp, 1; nshuffle=nshuffle, seed=seed, nbins=nbins,
                    stat=(F, c) -> harmonic_power(F, c; h=1))
    sig = fo.depth > maximum(fo.depth_null)
    ncyc = sum(segs[i].npulse * segs[i].f3_med / 2 for i in sel)
    verdict = sig ? :harmonic : (ncyc ≥ mincycles ? :separate : :inconclusive)
    return (verdict=verdict, is_harm=verdict == :harmonic, h1=fo.depth,
            h1_null_max=maximum(fo.depth_null), h1_null=fo.depth_null,
            h2=harmonic_power(fo.fold, fo.counts; h=2), ncycles=ncyc)
end


"""
    fundamental_track(tr, segs, harm) -> NamedTuple

Copy of the track with f3 halved (P3 doubled) in the windows of sections
flagged by `harmonic_groups`, so that every window carries the fundamental.
"""
function fundamental_track(tr, segs, harm)
    f3 = copy(tr.f3); f3_err = copy(tr.f3_err)
    for (i, s) in enumerate(segs)
        harm[i] || continue
        f3[s.win] ./= 2
        f3_err[s.win] ./= 2
    end
    return merge(tr, (f3=f3, f3_err=f3_err, p3=1 ./ f3, p3_err=f3_err ./ f3 .^ 2))
end


"""
    select_groups(segs, groups; ncycles=5) -> (segs, groups, dropped)

Keep only groups with at least `ncycles`·P3 pulses in total (P3 = pulse-
weighted median of the group); labels are re-ranked 1, 2, … by size.
`dropped` lists (old label, npulse, P3) of the removed groups.
"""
function select_groups(segs, groups; ncycles::Real=5)
    keep = Int[]
    dropped = Tuple{Int,Int,Float64}[]
    for g in sort(unique(groups))
        sel = findall(==(g), groups)
        npl = sum(segs[i].npulse for i in sel)
        p3g = group_p3(segs, sel)
        npl ≥ ncycles * p3g ? append!(keep, sel) : push!(dropped, (g, npl, p3g))
    end
    sort!(keep)
    old = groups[keep]
    ids = sort(unique(old))
    rank = Dict(k => r for (r, k) in enumerate(ids))
    return segs[keep], [rank[k] for k in old], dropped
end

"pulse-weighted P3 of the sections `sel`"
group_p3(segs, sel) = 1 / (sum(segs[i].f3_med * segs[i].npulse for i in sel) /
                           sum(segs[i].npulse for i in sel))


"linear interpolation of y(x) at x0, constant beyond the ends (x sorted)"
function _interp(x, y, x0)
    x0 ≤ x[1] && return y[1]
    x0 ≥ x[end] && return y[end]
    k = searchsortedlast(x, x0)
    t = (x0 - x[k]) / (x[k+1] - x[k])
    return (1 - t) * y[k] + t * y[k+1]
end


"""
    demodulate(data, on, pulses, f, L) -> Matrix{ComplexF64}

Complex modulation amplitude of every pulse n in `pulses` at its local
frequency f[j]: a Hann window of L pulses around n (shifted, not shrunk, at
the ends of the observation), column means removed, and

    Z(n, φ) = Σ_m taper(m) · X(m, φ) · e^{−2πi f (m − n)}.

The kernel is referenced to n itself, so for a modulation cos(θ(m) − ψ(φ))
arg Z(n, φ) = θ(n) − ψ(φ): the phase of the modulation *at that pulse*,
whatever the window position.

`loo` (leave-one-out, default): pulse n gets zero weight in its own window,
so its phase comes from its neighbours only. Without it the pulse's own
noise and jitter pull its phase towards wherever they best match the
template, and the fold of those very pulses shows structure even for
shuffled data (J0034-0721: shuffle-control depth 0.15–0.18 against 0.06 for
the constant-P3 fold).
"""
function demodulate(data::AbstractMatrix, on, pulses, f, L::Int; loo::Bool=true)
    N = size(data, 1)
    tap = hann_taper(L)
    Z = zeros(ComplexF64, length(pulses), length(on))
    for (j, n) in enumerate(pulses)
        s = clamp(n - L ÷ 2, 1, N - L + 1)
        X = data[s:s+L-1, on]
        X .-= mean(X, dims=1)
        k = tap .* cis.(-2π * f[j] .* ((s:s+L-1) .- n))
        loo && (k[n - s + 1] = 0)
        Z[j, :] .= vec(transpose(k) * X)
    end
    return Z
end


"""
    align_phases(Z; niter=10) -> (θ, T)

Modulation phase of every pulse against a common longitude template:
θ(n) = arg Σ_φ Z(n,φ)·conj T(φ), with T = mean_n Z(n,φ)·e^{−iθ(n)} iterated
from the strongest pulse. The template carries the longitude dependence of
the modulation phase (the drift band shape), so θ is one number per pulse,
consistent across sections separated by nulls or bad windows — that is what
lets separate sections of one group be folded together.
"""
function align_phases(Z::AbstractMatrix; niter::Int=10)
    T = Z[argmax(vec(sum(abs2, Z, dims=2))), :]
    θ = zeros(size(Z, 1))
    for _ in 1:niter
        θ .= angle.(Z * conj.(T))
        T = vec(mean(Z .* cis.(-θ), dims=1))
    end
    return θ, T
end


"""
    template_phase(T; frac=0.15, minrun=5, minpower=0.05) -> NamedTuple

Longitude dependence of the modulation phase, ψ(φ) = −arg T(φ), from a group
template (`align_phases`). This is what separates drift from amplitude
modulation in a fold: a drifting pattern has ψ changing steadily across the
emission (≈ W/P2 cycles), pure amplitude modulation has ψ flat within a
component (jumps of 1/2 cycle between components modulated in antiphase).

Bins with |T| ≥ frac·max|T| form contiguous runs (components). A run counts
only with ≥ `minrun` bins and ≥ `minpower` of the template power Σ|T|²
(frac = 0.2 with minpower = 0.1 cut off the weak trailing part of
J1825+0004's component — ~5% of the power, where ψ changes by ~¼ cycle — and
returned a confident "AM" from the flat main part alone; at 0.15/0.05 the
synthetic false-drift rate stays 0/80):
short low-amplitude runs are the overlap of components in antiphase (ψ steps
by ½ cycle where |T| cancels — a step, not a drift) or noise at the profile
edges, and with noise they read as strong gradients (synthetic antiphase AM
at high noise: a 4-bin run at |T| ≈ 0.3 gave Δψ = 0.62 ± 0.14).

Gradient per run from the amplitude-weighted phase increments,

    G_run = Σ_{j, j+1 ∈ run} conj(T_j)·T_{j+1},   slope = −arg(G_run)/2π  [cycles/bin],

no unwrapping (an unwrapped ψ random-walks through noisy bins), and a
½-cycle step where |T| is small weighs little. Δψ_run = slope × run length.
Runs are separate because the phase between separated components is defined
only mod 1 cycle, and drift may have opposite senses in them (bi-drifting).

Fields: psi (radians, unwrapped within mask runs, for plotting; NaN outside),
amp (|T|/max), mask, runs (counted runs), run_slope, run_dpsi (signed),
run_power (share of the counted runs' Σ|T|²),
dpsi = Σ|run_dpsi|, span = Σ run_dpsi, rms (|T|²-weighted RMS of ψ about the
run means [cycles]).
"""
function template_phase(T; frac::Real=0.15, minrun::Int=5, minpower::Real=0.05)
    amp = abs.(T) ./ maximum(abs.(T))
    mask = amp .>= frac
    psi = fill(NaN, length(T))
    raw = -angle.(T)
    ptot = sum(abs2, T)
    runs = UnitRange{Int}[]
    for r in true_runs(mask)
        psi[r] .= 2π .* unwrap_cycles(raw[r] ./ (2π))
        length(r) ≥ minrun && sum(abs2, T[r]) ≥ minpower * ptot && push!(runs, r)
    end
    run_slope = [run_gradient(T, r) for r in runs]
    run_dpsi = run_slope .* length.(runs)
    pr = [sum(abs2, T[r]) for r in runs]
    run_power = isempty(pr) ? Float64[] : pr ./ sum(pr)     # share of the counted runs' power
    ss = 0.0; sw = 0.0
    for r in runs
        w = amp[r] .^ 2; y = psi[r] ./ (2π)
        ym = sum(w .* y) / sum(w)
        ss += sum(w .* (y .- ym) .^ 2); sw += sum(w)
    end
    return (psi=psi, amp=amp, mask=mask, runs=runs, run_slope=run_slope, run_dpsi=run_dpsi,
            run_power=run_power, dpsi=sum(abs.(run_dpsi)), span=sum(run_dpsi),
            rms=sw > 0 ? sqrt(ss / sw) : NaN)
end

"""
    deep_dip(amp, w; depth=0.5) -> Bool

True when an interior bin of window `w` has amp < depth × the lower of the
maxima on its two sides — the cancellation point of two components
modulated in antiphase, where noise smears the ½-cycle phase step over
several bins (synthetic antiphase AM at high noise: 2/20 false partial
drifts without this check).
"""
function deep_dip(amp, w; depth::Real=0.5)
    for k in w.start+1:w.stop-1
        side = min(maximum(amp[w.start:k-1]), maximum(amp[k+1:w.stop]))
        amp[k] < depth * side && return true
    end
    return false
end

"amplitude-weighted phase gradient of ψ = −arg T over `r` [cycles/bin]"
run_gradient(T, r) = -angle(sum(conj(T[j]) * T[j+1] for j in r.start:r.stop-1)) / (2π)


"number of P3-phase bins: as `Functions.find_ybins` (2·P3, ≥ min_ppb pulses per bin, ≥ 4)"
fold_bins(p3, npulse; min_ppb=50) = max(4, min(floor(Int, 2 * p3), floor(Int, npulse / min_ppb)))


"fold rows `pulses` of `data[:, on]` into `nb` bins by phase ∈ [0,1); mean per bin"
function _fold(data, on, pulses, phase, nb)
    F = zeros(nb, length(on))
    cnt = zeros(Int, nb)
    for (j, n) in enumerate(pulses)
        b = clamp(floor(Int, phase[j] * nb) + 1, 1, nb)
        F[b, :] .+= data[n, on]
        cnt[b] += 1
    end
    return F ./ max.(cnt, 1), cnt
end


"""
    modulation_depth(F, cnt) -> Float64

RMS over longitude of the variance of the fold along P3 phase, as a fraction
of the squared peak of the mean profile: sqrt(mean_φ var_bins F(·,φ)) / max
profile. 0 for no modulation; it grows with how sharply the fold resolves the
pattern, so a fold with wrongly assigned phases smears and reads lower.
"""
function modulation_depth(F, cnt)
    w = cnt ./ sum(cnt)
    prof = vec(sum(F .* w, dims=1))
    v = vec(sum(w .* (F .- prof') .^ 2, dims=1))
    return sqrt(mean(v)) / maximum(prof)
end


"""
    phase_fold(data, sl, tr, segs, groups, g; nbins=nothing, niter=10,
               nshuffle=5, seed=1) -> NamedTuple

P3 fold of group `g` with compensation of variable P3. For every pulse of the
group's sections:

  1. local fundamental frequency f(n): `tr.f3` of the section's good windows
     interpolated over window centres (constant beyond the first/last one);
  2. complex modulation amplitude Z(n, φ) at f(n) (`demodulate`);
  3. phase θ(n) against the group template (`align_phases`).

Each pulse goes into bin floor(nb · (θ mod 2π)/2π). P3 never enters as a
constant: a monotonic or wandering P3 is followed through f(n), and phase
jumps (nulls, gaps between sections) are absorbed because θ is measured,
not integrated.

Self-alignment check: phases are measured on the data that are then folded,
so even noise aligned this way gives some fold structure. The whole chain
(steps 2–3 and the fold) is repeated on `nshuffle` pulse-order shuffles of
the group's pulses (same f(n), same pulses, order destroyed) and the
modulation depth of the real fold is reported against them. `stat` (default
`modulation_depth`) replaces the depth statistic for both, e.g.
`harmonic_power` in `harmonic_test`.

Fields: group, pulses, p3 (group P3), nb, fold, counts, phase (θ/2π mod 1),
theta (radians, per pulse), template, f (per pulse), depth, depth_null
(vector), coherence (mean over pulses of |Σ_φ Z conj T| / (‖Z‖·‖T‖):
how well single pulses match the template shape), tphase (`template_phase`:
longitude dependence of the modulation phase — the drift vs amplitude
modulation diagnostic), tsig (`template_significance`, block bootstrap over
the demodulations), tphase_fold / tsig_fold (the same from the fold's own
first harmonic, C(n,φ) = (I(n,φ) − ⟨I⟩)·e^{−iθ(n)}, bootstrapped pulse by
pulse — for short groups), verdict (tsig's, or tsig_fold's when the group
has < 5 bootstrap blocks), verdict_src (:block / :fold), bidrift,
p3_wander (std/median of the tracked P3 in the group's windows), phase_wander
(median rate at which the modulation phase leaves a constant-P3 clock
[cycles per 1000 pulses]), sections, on_bins.
"""
function phase_fold(data::AbstractMatrix, sl, tr, segs, groups, g::Int; nbins=nothing,
                    niter::Int=10, nshuffle::Int=5, seed::Int=1, stat=modulation_depth,
                    tp_kw=NamedTuple())
    on = sl.on_bins
    L = tr.window
    sel = findall(==(g), groups)
    pulses = Int[]
    f = Float64[]
    for i in sel
        s = segs[i]
        c = tr.centers[s.win]
        fw = tr.f3[s.win]
        for n in s.first:s.last
            push!(pulses, n)
            push!(f, _interp(c, fw, n))
        end
    end
    p3g = group_p3(segs, sel)
    nb = nbins === nothing ? fold_bins(p3g, length(pulses)) : nbins

    Z = demodulate(data, on, pulses, f, L)
    θ, T = align_phases(Z; niter=niter)
    phase = mod.(θ ./ (2π), 1.0)
    F, cnt = _fold(data, on, pulses, phase, nb)
    depth = stat(F, cnt)
    num = abs.(Z * conj.(T))
    den = sqrt.(vec(sum(abs2, Z, dims=2))) .* sqrt(sum(abs2, T))
    coherence = mean(num ./ max.(den, eps()))

    rng = Random.MersenneTwister(seed)
    depth_null = zeros(nshuffle)
    for k in 1:nshuffle
        d = copy(data)
        d[pulses, :] .= data[pulses[Random.randperm(rng, length(pulses))], :]
        Zs = demodulate(d, on, pulses, f, L)
        θs, _ = align_phases(Zs; niter=niter)
        Fs, cs = _fold(d, on, pulses, mod.(θs ./ (2π), 1.0), nb)
        depth_null[k] = stat(Fs, cs)
    end

    tp = template_phase(T; tp_kw...)
    # fold-based template: each pulse once, its own phase from the neighbours (leave-one-out),
    # so the pulses' contributions are independent and can be bootstrapped one by one —
    # a short group (few blocks of L/2 for the demodulation bootstrap) keeps its statistics
    X = data[pulses, on]
    C = (X .- mean(X, dims=1)) .* cis.(-θ)
    Tf = vec(mean(C, dims=1))
    tpf = template_phase(Tf; tp_kw...)
    tsf = template_significance(C, zeros(length(pulses)), Tf, tpf, L; seed=seed, block=1)
    tsb = template_significance(Z, θ, T, tp, L; seed=seed)

    # stability of the regime:
    # p3_wander — relative spread of the tracked P3 over the group's good windows (estimator
    #   scatter included; it is ≈ the resolution-limited jitter for a perfectly stable P3)
    # phase_wander — how fast the measured modulation phase leaves a constant-P3 clock:
    #   r(n) = θ(n)/2π − n/P3 within each section, median |r(n+Δ) − r(n)| over Δ = min(200,
    #   section/2) pulses, scaled to cycles per 1000 pulses (≈ 0.5 means a constant-P3 fold
    #   over 1000 pulses would be smeared by half a cycle)
    p3w = Float64[]
    for i in sel
        append!(p3w, tr.p3[segs[i].win])
    end
    p3_wander = length(p3w) > 1 ? std(p3w) / median(p3w) : NaN
    rates = Float64[]
    for i in sel
        sg = segs[i]
        idx = findall(n -> sg.first ≤ n ≤ sg.last, pulses)
        length(idx) < 20 && continue
        r = unwrap_cycles(θ[idx] ./ (2π) .- pulses[idx] ./ p3g)
        Δ = min(200, length(idx) ÷ 2)
        d = abs.(r[1+Δ:end] .- r[1:end-Δ])
        push!(rates, median(d) * 1000 / Δ)
    end
    phase_wander = isempty(rates) ? NaN : median(rates)
    # final verdict: the block bootstrap when the group has enough blocks, otherwise the
    # pulse-by-pulse fold estimate (synthetic short groups: no false drift, N = 40 drift 8/8 vs
    # 6/8; J1528-4109, 34 P: drift; short P3-only groups with z ≈ 11–13 from 2–3 blocks: z ≤ 1.9)
    fromfold = tsb.nblocks < 5
    return (group=g, pulses=pulses, p3=p3g, nb=nb, fold=F, counts=cnt, phase=phase, theta=θ,
            template=T, tphase=tp, tsig=tsb, tphase_fold=tpf, tsig_fold=tsf,
            verdict=fromfold ? tsf.verdict : tsb.verdict, verdict_src=fromfold ? :fold : :block,
            bidrift=fromfold ? tsf.bidrift : tsb.bidrift, p3_wander=p3_wander,
            phase_wander=phase_wander,
            f=f, depth=depth, depth_null=depth_null, coherence=coherence, on_bins=on,
            sections=sel)
end


"""
    template_significance(Z, θ, T, tp, L; nboot=300, block=nothing, seed=1,
                          dpsi_min=0.1, zdet=5.0, minblocks=5, lwin=5, maxfrac=0.5,
                          dipdepth=0.5, psimax=deg2rad(20), dpsi_am=0.25, zam=2.0,
                          minpow_drift=0.5)
        -> NamedTuple

Uncertainty and significance of the template phase gradient (`template_phase`
output `tp`) from a moving-block bootstrap over the group's pulses.

Why not the pulse-order shuffle used for the fold depth: it destroys the
modulation, the shuffled template is noise with random phase, and its Δψ
answers "is there modulation", not "is the modulation phase flat".

Bootstrap: blocks of `block` consecutive pulses of the group (default L÷2 —
the demodulation windows of pulses closer than ~L/2 share most of their
pulses, so single pulses are not independent) are drawn with replacement up
to the group size; T* = mean Z·e^{−iθ} over the drawn pulses (θ kept from
the full alignment) and the run gradients (`run_gradient`) are recomputed in
the *same* runs (components) as `tp`. σ of every run's Δψ is the bootstrap
standard deviation.

Under amplitude modulation every run's Δψ has expectation 0, so

    χ² = Σ_runs (Δψ_run / σ_run)²  ~  χ²(n_runs),

giving p and its one-sided normal equivalent z. A very bright AM pulsar can
show a tiny but formally significant gradient (profile asymmetry, slightly
offset components), so the verdict also needs a size:

  :drift        z ≥ `zdet` and Δψ ≥ `dpsi_min`, and the components that drift on
                their own (z_run ≥ `zdet`, |Δψ_run| ≥ `dpsi_min`) carry ≥
                `minpow_drift` = ½ of the template power — definition (b): the
                pattern moves through the dominant part of the emission. A
                significant gradient confined to weaker components is :partial
                (J1057-5226: flat dominant component, −0.19 ± 0.03 in the weaker
                one; J1543+0929: no component significant on its own)
  :partial      a global gradient that fails the power rule above, or some
                `lwin`-bin window inside the emission
                mask has a local phase change |Δψ_w| ≥ `dpsi_min` with
                z_w ≥ `zdet` (same bootstrap), spread over its bins — no
                single bin-to-bin increment carries more than `maxfrac` of
                the window's net change, no deep amplitude minimum inside
                the window (`deep_dip`), and every bin's phase known to
                ≤ `psimax` (bootstrap circular σ; noise-dominated bins at the
                profile edges jump by 60–100°/bin and, at high noise, gave
                2/20 false partial drifts in synthetic antiphase AM). A ½-cycle step where two components
                in antiphase cancel is one increment (or, with noise, a few
                around a deep |T| minimum) and fails that; a
                gradient confined to part of the emission passes (J1825+0004:
                flat phase over the bright peak, ~0.7 cycle change down the
                trailing flank, the same in all four time quarters — "partial
                drift", decision of 2026-10-01)
  :am           no gradient (z < `zam` = 2; at 3, J1511-5414 with a visible ramp
                of +0.09 ± 0.03, z = 2.7, read as AM), no partial window, and
                the upper limit Σ|Δψ_run| + 2·√(Σσ_run²) < `dpsi_am` = 0.25
                cycle — half the smallest firm drift seen (groups with
                verdict drift have Δψ ≈ 0.5–2.1). The first version required
                Σ(|Δψ_run| + 2σ_run) < 0.1, linear in the errors and far below
                any drift: no group of the 20-pulsar pilot reached it (J0709-5923
                Δψ = 0.04 with limit 0.17, J0849-6322 0.04 with 0.13).
  :inconclusive otherwise, or fewer than `minblocks` independent blocks
                (npulse / block) — too little data for a bootstrap

`dpsi_min` = 0.1 cycle is a working value, not calibrated physics.
Fields: verdict, z, p, chi2, nruns, sigma_run, z_run, dpsi_upper, nblocks, block,
partial (found, z, dpsi, bins, maxfrac, nwin — the best local window),
sigma_psi (per-bin phase error [rad]), pow_drift (power share of the
components drifting on their own), bidrift (two such components with
opposite drift senses).
"""
function template_significance(Z, θ, T, tp, L::Int; nboot::Int=300, block=nothing, seed::Int=1,
                               dpsi_min::Real=0.1, zdet::Real=5.0, minblocks::Real=5,
                               lwin::Int=5, maxfrac::Real=0.5, dipdepth::Real=0.5,
                               psimax::Real=deg2rad(20), dpsi_am::Real=0.25, zam::Real=2.0,
                               minpow_drift::Real=0.5)
    n = size(Z, 1)
    bl = block === nothing ? max(1, L ÷ 2) : block
    nblocks = n / bl
    runs = tp.runs
    nr = length(runs)
    # candidate windows for a local (partial) gradient: lwin bins inside any mask run
    wins = UnitRange{Int}[]
    for r in true_runs(tp.mask)
        for s0 in r.start:(r.stop - lwin + 1)
            push!(wins, s0:s0+lwin-1)
        end
    end
    nolocal = (found=false, z=NaN, dpsi=NaN, bins=0:-1, maxfrac=NaN, nwin=length(wins))
    if nr == 0 && isempty(wins)
        return (verdict=:inconclusive, z=NaN, p=NaN, chi2=NaN, nruns=0, sigma_run=Float64[],
                z_run=Float64[], dpsi_upper=NaN, nblocks=nblocks, block=bl, partial=nolocal,
                sigma_psi=fill(NaN, length(T)), pow_drift=0.0, bidrift=false)
    end
    Zr = Z .* cis.(-θ)
    rng = Random.MersenneTwister(seed)
    starts = 1:max(1, n - bl + 1)
    boot = zeros(nboot, nr)
    bootw = zeros(nboot, length(wins))
    cphase = zeros(ComplexF64, length(T))      # Σ e^{iδ}: per-bin phase scatter of T*
    idx = Int[]
    for b in 1:nboot
        empty!(idx)
        while length(idx) < n
            s0 = rand(rng, starts)
            append!(idx, s0:min(n, s0 + bl - 1))
        end
        resize!(idx, n)
        Tb = vec(mean(Zr[idx, :], dims=1))
        cphase .+= cis.(angle.(Tb .* conj.(T)))
        for (j, r) in enumerate(runs)
            boot[b, j] = run_gradient(Tb, r) * length(r)
        end
        for (j, w) in enumerate(wins)
            bootw[b, j] = run_gradient(Tb, w) * length(w)
        end
    end

    # global (whole components)
    if nr > 0
        σ = vec(std(boot, dims=1))
        zr = tp.run_dpsi ./ σ
        chi2 = sum(abs2, zr)
        p = ccdf(Chisq(nr), chi2)
        z = cquantile(Normal(), max(p, 1e-300))   # ≤ ~37; 1 − p would round to 1 below 1e-16
        upper = sum(abs.(tp.run_dpsi)) + 2 * sqrt(sum(abs2, σ))   # errors in quadrature
    else
        σ = Float64[]; zr = Float64[]; chi2 = NaN; p = NaN; z = NaN; upper = NaN
    end

    # per-bin phase error (circular σ of the bootstrap)
    σψ = sqrt.(-2 .* log.(clamp.(abs.(cphase ./ nboot), 1e-12, 1.0)))

    # local: strongest lwin-bin window with a significant, spread-out phase change
    partial = nolocal
    best = -Inf
    for (j, w) in enumerate(wins)
        d = run_gradient(T, w) * length(w)
        sw = std(view(bootw, :, j))
        zw = abs(d) / sw
        inc = [-angle(conj(T[k]) * T[k+1]) for k in w.start:w.stop-1]
        tot = abs(sum(inc))
        mf = tot > 0 ? maximum(abs.(inc)) / tot : Inf
        ok = zw ≥ zdet && abs(d) ≥ dpsi_min && mf ≤ maxfrac && !deep_dip(tp.amp, w; depth=dipdepth) &&
             maximum(σψ[w]) ≤ psimax
        if ok && zw > best
            best = zw
            partial = (found=true, z=zw, dpsi=d, bins=w, maxfrac=mf, nwin=length(wins))
        end
    end

    # power share of components that drift on their own (z_run ≥ zdet, |Δψ_run| ≥ dpsi_min)
    qual = nr > 0 ? (abs.(zr) .≥ zdet) .& (abs.(tp.run_dpsi) .≥ dpsi_min) : falses(0)
    pow_drift = nr > 0 ? sum(tp.run_power[qual]) : 0.0
    # bi-drifting: components drifting on their own with opposite senses (J1537-4912:
    # −0.19 ± 0.02 and +0.17 ± 0.03, the signs stable in all four time quarters)
    bidrift = nr > 0 && any(qual .& (tp.run_dpsi .> 0)) && any(qual .& (tp.run_dpsi .< 0))
    glob = nr > 0 && z ≥ zdet && tp.dpsi ≥ dpsi_min
    verdict = nblocks < minblocks ? :inconclusive :
              (glob && pow_drift ≥ minpow_drift) ? :drift :
              (glob || partial.found) ? :partial :
              (nr > 0 && z < zam && upper < dpsi_am) ? :am : :inconclusive
    return (verdict=verdict, z=z, p=p, chi2=chi2, nruns=nr, sigma_run=σ, z_run=zr,
            dpsi_upper=upper, nblocks=nblocks, block=bl, partial=partial, sigma_psi=σψ,
            pow_drift=pow_drift, bidrift=bidrift)
end


"""
    constant_fold(data, on, pulses, p3, nb) -> NamedTuple

Reference fold of the same pulses with one constant P3 and phase n/P3 from
the absolute pulse number — what `Tools.p3fold` does, restricted to the
group's pulses. Fields: fold, counts, depth.
"""
function constant_fold(data, on, pulses, p3, nb)
    F, cnt = _fold(data, on, pulses, mod.(pulses ./ p3, 1.0), nb)
    return (fold=F, counts=cnt, depth=modulation_depth(F, cnt))
end


"""
    analyse(data, bin_st, bin_end, p3; window=nothing, ncycles=5, nshuffle=5,
            second_pass=true) -> NamedTuple

Whole chain for one pulsar: window L = `window_length(p3)` (p3 from
params.json), sliding LRFS, P3 track, local shuffle threshold, sections,
groups (`p3_groups`), second harmonics joined to their fundamental
(`harmonic_groups`, `fundamental_track`), touching sections merged, groups
shorter than `ncycles`·P3 dropped, and a variable-P3 fold (`phase_fold`) plus
a constant-P3 reference fold (`constant_fold`) for every remaining group.

With `second_pass` the pulses left outside every group get a second look
with a longer window (`long_p3_pass`): a window chosen for P3 from
params.json only measures P3 ≲ L/3, so a second regime with a longer P3 is
invisible to it (J1825+0004 after the mode change at ~715: P3 ≈ 35–55 next
to 14.5; J0034-0721 mode A, P3 ~ 12, at L = 26).

Fields: L, sl, tr (fundamental track), threshold, good, segs, groups, harm
(before merging), harm_tests (`harmonic_groups`), dropped, folds, cfolds —
first pass; pass2 — the same
fields for the second pass plus p3_probe, probe (see `long_p3_pass`), or
nothing when there was no second pass; nyq — `nyquist_pass` result when
p3 ≤ `nyq_p3max` (2.2), run on the pulses outside the first-pass groups
before the second pass (which then also skips the Nyquist sections),
otherwise nothing.
"""
function analyse(data::AbstractMatrix, bin_st::Int, bin_end::Int, p3::Real; window=nothing,
                 ncycles::Real=5, nshuffle::Int=5, second_pass::Bool=true, tp_kw=NamedTuple(),
                 nyq_p3max::Real=2.2)
    N = size(data, 1)
    L = window === nothing ? window_length(p3) : window
    sl = sliding_lrfs(data, bin_st, bin_end; window=L)
    tr = p3_track(sl)
    thr = contrast_null(data, sl).threshold
    good = good_windows(tr, thr)
    segs = p3_segments(tr, good, N)
    groups = p3_groups(segs, L)
    groups, harm, htests = harmonic_groups(segs, groups, L;
        check=sel -> harmonic_test(data, sl, tr, segs, sel; nshuffle=max(20, nshuffle)))
    trf = fundamental_track(tr, segs, harm)
    keep0 = groups .!= 0               # inconclusive 2:1 candidates out
    segs, groups = segs[keep0], groups[keep0]
    segs, groups = merge_sections(trf, segs, groups, N)
    segs, groups, dropped = select_groups(segs, groups; ncycles=ncycles)
    folds = [phase_fold(data, sl, trf, segs, groups, g; nshuffle=nshuffle, tp_kw=tp_kw)
             for g in sort(unique(groups))]
    cfolds = [constant_fold(data, sl.on_bins, fo.pulses, fo.p3, fo.nb) for fo in folds]
    nyq = p3 ≤ nyq_p3max ? nyquist_pass(data, bin_st, bin_end, p3, segs; tp_kw=tp_kw) : nothing
    taken = nyq === nothing || !nyq.found ? segs :
            vcat(segs, [(first=first(r), last=last(r)) for r in nyq.sections])
    pass2 = second_pass ? long_p3_pass(data, bin_st, bin_end, L, taken;
                                       ncycles=ncycles, nshuffle=nshuffle, tp_kw=tp_kw) : nothing
    return (L=L, sl=sl, tr=trf, threshold=thr, good=good, segs=segs, groups=groups,
            harm=harm, harm_tests=htests, dropped=dropped, folds=folds, cfolds=cfolds, pass2=pass2,
            nyq=nyq)
end


"""
    nyquist_pass(data, bin_st, bin_end, p3, segs1; B=32, step=16, nperm=300, q=0.99,
                 seed=1, mincycles=2, nbins=8, tp_kw) -> NamedTuple

Drift detection at P3 ≈ 2 (f3 near the Nyquist frequency 0.5), where the
sliding LRFS fails: the feature at f3 and its mirror at 1 − f3 merge within
the Hann main lobe (±2/L) into one maximum exactly at f = 0.5, the edge of
the search range, and the window is rejected (`feature_peak` needs an
interior maximum). In batch v2/v3 no pulsar with P3 ≤ 2.13 (17) had a group
at its catalogue P3, the smallest group P3 was 2.15. Found and tested by the
session claude-ac (`~/claude/work/scripts/play/nyquist_*.jl`).

On pulses outside the first-pass sections (`segs1`):

1. Test B: blocks of `B` pulses (step `step`), folded at P3 = 2,
   A(φ) = Σ (−1)ⁿ xₙ(φ), S = Σ_φ |A|² against `nperm` shuffles of the
   block's pulses; blocks above the `q` quantile are significant and are
   merged into sections. (A sign change of A across the profile is what a
   drift looks like here; a uniform sign is alternating brightness.)
2. f3 from the incoherently summed periodogram of the sections
   (rectangular window, lobe ±1/M instead of Hann's ±2/M) on [0.40, 0.5];
   δ = 0.5 − f3.
3. f3 and its alias 0.5 + δ are separable only if the longest section M has
   M·2δ ≥ `mincycles`. If so: demodulation at f3 in each section, section
   phases aligned to the strongest one, template and its significance
   pulse by pulse (as `verdict_fold`) → drift / partial / am / inconclusive.
   The fold at the alias is the mirror image (Δψ with the opposite sign), so
   only |Δψ| and the relative signs of components (bi-drift) are meaningful;
   **the drift direction is unknown**. If not separable: verdict :nyquist —
   modulation at Nyquist without a phase measurement (at f = 0.5 exactly the
   demodulation is real and the template phase is 0 or π only).

J0846-3533 (claude-ac): sections 976 P, f3 = 0.4938 (P3 = 2.025, params
2.02), drift, |Δψ| = 0.61, z = 7.9 — no group in v2. J0943+2253: sections
≤ 96 P, f3 not separable from 0.5 → :nyquist.

Fields: found, sections, block_starts, block_z, block_sig, f3, delta, p3,
p3_alias, M, resolved, pulses, theta, verdict, tphase, tsig, fold, counts,
nb, depth, bidrift.
"""
function nyquist_pass(data::AbstractMatrix, bin_st::Int, bin_end::Int, p3::Real, segs1;
                      B::Int=32, step::Int=16, nperm::Int=300, q::Real=0.99, seed::Int=1,
                      mincycles::Real=2, nbins::Int=8, tp_kw=NamedTuple())
    N = size(data, 1)
    on = bin_st:bin_end
    free = free_pulses(segs1, N)
    sgn = [(-1.0)^n for n in 0:B-1]
    rng = Random.MersenneTwister(seed)
    starts = collect(1:step:max(1, N - B + 1))
    z = fill(NaN, length(starts)); sig = falses(length(starts))
    for (k, a) in enumerate(starts)
        a + B - 1 ≤ N || continue
        all(free[a:a+B-1]) || continue
        X = data[a:a+B-1, on]; X .-= mean(X, dims=1)
        all(iszero, X) && continue
        S = sum(abs2, sgn' * X)
        nul = [sum(abs2, sgn' * X[Random.randperm(rng, B), :]) for _ in 1:nperm]
        z[k] = (S - mean(nul)) / std(nul)
        sig[k] = S > quantile(nul, q)
    end
    secs = UnitRange{Int}[]
    for k in findall(sig)
        r = starts[k]:starts[k]+B-1
        if !isempty(secs) && first(r) ≤ last(secs[end]) + 1
            secs[end] = first(secs[end]):last(r)
        else
            push!(secs, r)
        end
    end
    empty = (found=false, sections=secs, block_starts=starts, block_z=z, block_sig=sig, f3=NaN,
             delta=NaN, p3=NaN, p3_alias=NaN, M=0, resolved=false, pulses=Int[], theta=Float64[],
             verdict=:none, tphase=nothing, tsig=nothing, fold=zeros(0, 0), counts=Int[], nb=0,
             depth=NaN, bidrift=false)
    isempty(secs) && return empty

    # 2. f3 from the sections' periodogram
    fgrid = collect(0.40:2e-5:0.5)
    Ptot = zeros(length(fgrid))
    for r in secs
        X = data[r, on]; X .-= mean(X, dims=1)
        n = collect(0:length(r)-1)
        for (i, f) in enumerate(fgrid)
            Ptot[i] += sum(abs2, transpose(cis.(-2π * f .* n)) * X)
        end
    end
    f3 = fgrid[argmax(Ptot)]
    δ = 0.5 - f3
    M = maximum(length.(secs))
    resolved = M * 2δ ≥ mincycles
    pulses = reduce(vcat, collect.(secs))
    if !resolved
        return merge(empty, (found=true, f3=f3, delta=δ, p3=1 / f3, p3_alias=1 / (0.5 + δ), M=M,
                             pulses=pulses, verdict=:nyquist))
    end

    # 3. demodulation at f3, section phases aligned to the strongest section
    Ts = [vec(mean((data[r, on] .- mean(data[r, on], dims=1)) .* cis.(-2π * f3 .* r), dims=1)) for r in secs]
    w = [length(r) * sum(abs2, t) for (r, t) in zip(secs, Ts)]
    ref = Ts[argmax(w)]
    rot = [angle(sum(t .* conj(ref))) for t in Ts]
    θ = reduce(vcat, [2π * f3 .* r .- rot[i] for (i, r) in enumerate(secs)])
    X = data[pulses, on]
    C = (X .- mean(X, dims=1)) .* cis.(-θ)
    Tf = vec(mean(C, dims=1))
    tp = template_phase(Tf; tp_kw...)
    ts = template_significance(C, zeros(length(pulses)), Tf, tp, B; seed=seed, block=1)
    F, cnt = _fold(data, on, pulses, mod.(θ ./ (2π), 1.0), nbins)
    return merge(empty, (found=true, f3=f3, delta=δ, p3=1 / f3, p3_alias=1 / (0.5 + δ), M=M,
                         resolved=true, pulses=pulses, theta=θ, verdict=ts.verdict, tphase=tp,
                         tsig=ts, fold=F, counts=cnt, nb=nbins, depth=modulation_depth(F, cnt),
                         bidrift=ts.bidrift))
end


"""
    plot_nyquist(data, nyq, outdir; nbin, name_mod, show_)

Nyquist path (`nyquist_pass`) for one pulsar: single pulses with the
sections of significant blocks, block z of test B, fold at f3 (direction of
drift unknown — the alias fold is its mirror image) and the template phase.
Writes `<name_mod>_nyquist.png` (+ .pdf with `pdf=true`).
"""
function plot_nyquist(data, nyq, outdir; nbin=size(data, 2), on=nothing, name_mod="pulsar",
                      pdf=true, show_=false)
    N = size(data, 1)
    on === nothing && error("on (on-pulse bins) required")
    lon = (collect(on) .- 1) .* 360.0 ./ nbin
    rc("font", size=7.)
    rc("axes", linewidth=0.5)
    rc("lines", linewidth=0.5)
    fig = figure(figsize=(7.0, 8.0))
    ax1 = fig.add_axes([0.10, 0.72, 0.86, 0.23])
    st = permutedims(data[:, on])
    imshow(st, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
           extent=(0.5, N + 0.5, lon[1], lon[end]), vmin=quantile(vec(st), 0.01),
           vmax=quantile(vec(st), 0.995))
    for r in nyq.sections
        ax1.axvspan(first(r) - 0.5, last(r) + 0.5, ymin=0.0, ymax=0.05, color="C3", lw=0)
    end
    ylabel("longitude (\$^\\circ\$)")
    tick_params(labelbottom=false)
    title(@sprintf("%s  Nyquist path: f\$_3\$ = %.5f (P\$_3\$ = %.4f or alias %.4f), M = %d, %s → %s",
                   name_mod, nyq.f3, nyq.p3, nyq.p3_alias, nyq.M,
                   nyq.resolved ? "resolved" : "aliases not resolved", nyq.verdict), fontsize=7)
    ax2 = fig.add_axes([0.10, 0.58, 0.86, 0.12], sharex=ax1)
    c = nyq.block_starts .+ 15.5
    plot(c, nyq.block_z, c="0.4", marker=".", ms=2)
    scatter(c[nyq.block_sig], nyq.block_z[nyq.block_sig], s=5, c="C3", zorder=3)
    xlim(0.5, N + 0.5); xlabel("pulse number"); ylabel("z (test B)")
    if nyq.resolved
        ax3 = fig.add_axes([0.10, 0.07, 0.40, 0.42])
        FF = vcat(nyq.fold, nyq.fold)
        imshow(FF, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
               extent=(lon[1], lon[end], 0, 2), vmax=quantile(vec(FF), 0.99))
        xlabel("longitude (\$^\\circ\$)"); ylabel("P\$_3\$ phase (cycles)")
        title(@sprintf("fold at f\$_3\$, %d P, %d bins, depth %.3f\n(alias fold = mirror image)",
                       length(nyq.pulses), nyq.nb, nyq.depth), fontsize=6)
        ax4 = fig.add_axes([0.58, 0.07, 0.34, 0.42])
        tp, ts = nyq.tphase, nyq.tsig
        plot(lon, tp.amp, c="lightgrey", lw=0.8); ylim(0, 1.05)
        xlabel("longitude (\$^\\circ\$)"); ylabel("|T| (norm.)")
        ax5 = ax4.twinx()
        ax5.plot(lon, rad2deg.(tp.psi), ".", ms=2, c="C3")
        ax5.set_ylabel("ψ (\$^\\circ\$), sign arbitrary")
        title(@sprintf("|Δψ| %.2f cyc (%s)\nz = %.1f → %s%s", tp.dpsi,
                       join([@sprintf("%+.2f±%.2f", a, b) for (a, b) in zip(tp.run_dpsi, ts.sigma_run)], ", "),
                       ts.z, ts.verdict, ts.bidrift ? ", bi-drift" : ""), fontsize=6)
    end
    savepath = joinpath(outdir, "$(name_mod)_nyquist.pdf")
    pdf && savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"), dpi=150)
    println(pdf ? savepath : replace(savepath, ".pdf" => ".png"))
    if show_
        PyPlot.show()
        println("Press Enter to close the figure.")
        readline(stdin; keep=false)
    end
    close()
end


"""
    prominence(y, k) -> Float64

Height of the local maximum y[k] above the higher of the two minima met
walking left and right from k until y exceeds y[k] (or the array ends).
"""
function prominence(y, k)
    lm = y[k]; i = k - 1
    while i ≥ 1 && y[i] ≤ y[k]
        lm = min(lm, y[i]); i -= 1
    end
    rm = y[k]; i = k + 1
    while i ≤ length(y) && y[i] ≤ y[k]
        rm = min(rm, y[i]); i += 1
    end
    return y[k] - max(lm, rm)
end


"""
    free_pulses(segs, N) -> BitVector

Pulses not covered by any section of `segs`.
"""
function free_pulses(segs, N::Int)
    free = trues(N)
    for s in segs
        free[s.first:s.last] .= false
    end
    return free
end

"contiguous runs of `true` in a BitVector, as UnitRanges"
function true_runs(m)
    runs = UnitRange{Int}[]
    i = 1
    while i ≤ length(m)
        if m[i]
            j = i
            while j < length(m) && m[j+1]
                j += 1
            end
            push!(runs, i:j)
            i = j + 1
        else
            i += 1
        end
    end
    return runs
end


"""
    long_p3_pass(data, bin_st, bin_end, L1, segs1; ncycles=5, nshuffle=5,
                 lprobe_max=256, min_free=nothing) -> NamedTuple or nothing

Second pass for regimes with a P3 longer than the first window L1 can
measure (P3 > L1/3), restricted to pulses outside the first-pass sections.

1. Probe: contrast spectra (each window divided by its median) of windows
   of Lp = min(`lprobe_max`, longest free stretch) pulses lying entirely in
   free stretches (step Lp/8), averaged; P3' = local maximum with
   3/Lp ≤ f < 3/L1 of the largest *prominence* (`prominence`) — the highest
   one is usually a bump on the red continuum next to 3/Lp (J1825+0004:
   P3' = 57, J2307+2225: 85). No P3' (or free stretches shorter than
   `min_free`, default 2·L1) → nothing.
2. Windows tried: min(window_length(P3'), longest free stretch) and the
   ladder L1·{2, 3, 4, 6, 8} up to the longest free stretch (powers of two
   alone skipped the working window for J1825+0004: 114, 228 but not ~160) — a wandering long
   P3 smears the probe spectrum, so the probe alone picks a window too short
   (J1825+0004: P3' = 21, L2 = 84, while the regime has P3 55 → 23). The
   window giving the most pulses in kept groups wins (`ladder` lists
   (L2, pulses) for all tried).
3. For each window: sliding LRFS at L2 on the whole observation, track searched only in
   f < 3/L1, local shuffle threshold over the same range; good windows must
   also have their centre on a free pulse and ≥ half of the window free.
4. Sections, groups, merging as in the first pass; every section is then
   cut to the longest free stretch inside it (no pulse is in two passes);
   groups shorter than `ncycles`·P3 dropped; folds as in the first pass.

Fields: as the first pass in `analyse` (L, sl, tr, threshold, good, segs,
groups, dropped, folds, cfolds) plus p3_probe, probe (Lp, freq, spectrum,
nwin, prominence), free (BitVector after the first pass), ladder.
"""
function long_p3_pass(data::AbstractMatrix, bin_st::Int, bin_end::Int, L1::Int, segs1;
                      ncycles::Real=5, nshuffle::Int=5, lprobe_max::Int=256, min_free=nothing,
                      tp_kw=NamedTuple())
    N = size(data, 1)
    free = free_pulses(segs1, N)
    runs = true_runs(free)
    isempty(runs) && return nothing
    longest = maximum(length.(runs))
    mf = min_free === nothing ? 2 * L1 : min_free
    longest < mf && return nothing

    # 1. probe spectrum over free stretches
    Lp = min(lprobe_max, longest)
    fhi = 3 / L1
    acc = nothing; nw = 0; freq = Float64[]
    for r in runs
        length(r) < Lp && continue
        sl = sliding_lrfs(data[r, :], bin_st, bin_end; window=Lp, stride=max(1, Lp ÷ 8),
                          off_bins=nothing)
        use = sl.freq .>= sl.fmin
        S = sl.power ./ [median(sl.power[i, use]) for i in 1:size(sl.power, 1)]
        acc = acc === nothing ? vec(sum(S, dims=1)) : acc .+ vec(sum(S, dims=1))
        nw += size(S, 1)
        freq = sl.freq
    end
    nw == 0 && return nothing
    spec = acc ./ nw
    cand = [k for k in 2:length(freq)-1 if 3 / Lp ≤ freq[k] < fhi &&
            spec[k] > spec[k-1] && spec[k] ≥ spec[k+1]]
    isempty(cand) && return nothing
    prom = [prominence(spec, k) for k in cand]
    kbest = cand[argmax(prom)]
    p3p = 1 / freq[kbest]
    probe = (Lp=Lp, freq=freq, spectrum=spec, nwin=nw, prominence=maximum(prom))

    # 2. windows to try: the probe's, and a ladder L1·{2,3,4,6,8} — a wandering
    #    long P3 smears the probe spectrum (J1825+0004 after ~715: P3 55 → 23,
    #    probe P3' = 21 → L2 = 84 sees only P3 ≤ 28), so no single P3' is reliable
    ladder = Int[]
    Lw = min(window_length(p3p), longest)
    Lw > L1 && push!(ladder, Lw)
    for k in (2, 3, 4, 6, 8)
        k * L1 ≤ longest && push!(ladder, k * L1)
    end
    ladder = sort(unique(filter(l -> L1 < l ≤ N ÷ 2, ladder)))
    isempty(ladder) && return nothing

    # 3–4. one pass per window, keep the one with most pulses in groups
    best = nothing
    tried = Tuple{Int,Int}[]
    for L2 in ladder
        r = _long_pass_at(data, bin_st, bin_end, L2, L1, free, N; ncycles=ncycles,
                          nshuffle=nshuffle, tp_kw=tp_kw)
        npl = isempty(r.segs) ? 0 : sum(s.npulse for s in r.segs)
        push!(tried, (L2, npl))
        if best === nothing || npl > best[2]
            best = (r, npl)
        end
    end
    r = best[1]
    return merge(r, (p3_probe=p3p, probe=probe, free=free, ladder=tried))
end

function _long_pass_at(data, bin_st, bin_end, L2, L1, free, N; ncycles=5, nshuffle=5,
                       tp_kw=NamedTuple())
    fhi = 3 / L1
    fr = (2 / L2, fhi)
    sl = sliding_lrfs(data, bin_st, bin_end; window=L2)
    tr = p3_track(sl; frange=fr)
    thr = contrast_null(data, sl; frange=fr).threshold
    good = good_windows(tr, thr)
    for (i, s) in enumerate(sl.starts)
        c = round(Int, sl.centers[i])
        good[i] &= free[c] && count(free[s:s+L2-1]) ≥ L2 / 2
    end
    segs = p3_segments(tr, good, N)
    groups = p3_groups(segs, L2)
    segs, groups = merge_sections(tr, segs, groups, N)
    keep = Int[]
    cut = NamedTuple[]
    for (i, s) in enumerate(segs)
        fr_runs = true_runs(free[s.first:s.last])
        isempty(fr_runs) && continue
        r = fr_runs[argmax(length.(fr_runs))] .+ (s.first - 1)
        push!(cut, merge(s, (first=r.start, last=r.stop, npulse=length(r))))
        push!(keep, i)
    end
    segs = cut; groups = groups[keep]
    dropped = Tuple{Int,Int,Float64}[]
    folds = NamedTuple[]; cfolds = NamedTuple[]
    if !isempty(segs)
        segs, groups, dropped = select_groups(segs, groups; ncycles=ncycles)
        folds = [phase_fold(data, sl, tr, segs, groups, g; nshuffle=nshuffle, tp_kw=tp_kw)
                 for g in sort(unique(groups))]
        cfolds = [constant_fold(data, sl.on_bins, fo.pulses, fo.p3, fo.nb) for fo in folds]
    end
    return (L=L2, sl=sl, tr=tr, threshold=thr, good=good, segs=segs, groups=groups,
            harm=falses(length(segs)), dropped=dropped, folds=folds, cfolds=cfolds)
end


"""
    plot_track(data, sl, tr, outdir; nbin, name_mod, good, threshold, snr_min,
               p3_ref, p3_lim, darkness, segs, groups, suffix, show_)

Sliding LRFS in the layout of Fig. 4 of Szary et al. (2022).

Panel 1: single pulses (on-pulse window), longitude vs pulse number.
Panel 2 (left): time-averaged fluctuation spectrum. Panel 2 (main):
spectra vs window centre, each divided by its own median over f ≥ fmin
(contrast); the dotted line marks fmin = 2/L, below which the DC leakage
sits. Panel 3: P3 per window. Black: windows in `good` (e.g. from
`good_windows`; default: fit converged and peak above `snr_min` off-pulse
σ); light grey: the rest. Dashed: `p3_ref`.

With `segs` and `groups` the windows of each section are coloured by group,
the pulse ranges of the sections are marked as a strip along the bottom of
panel 1, and every group is summarised (P3, sections, pulses) in panel 3.

The x-axis is the window *centre* (the paper uses the start pulse) so that
features line up with the single-pulse panel.

Writes `<name_mod>_sliding_lrfs_L<window><suffix>.pdf/.png`.
"""
function plot_track(data, sl, tr, outdir; nbin=size(data, 2), name_mod="pulsar",
                    good=nothing, threshold=nothing, snr_min=5.0, p3_ref=nothing,
                    p3_lim=nothing, darkness=0.995, segs=nothing, groups=nothing,
                    suffix="", pdf=true, show_=false)
    N = size(data, 1)
    on = sl.on_bins
    lon = (collect(on) .- 1) .* 360.0 ./ nbin
    L = sl.window

    rc("font", size=7.)
    rc("axes", linewidth=0.5)
    rc("lines", linewidth=0.5)
    fig = figure(figsize=(7.0, 6.0))

    # Panel 1: single pulses
    ax1 = fig.add_axes([0.20, 0.68, 0.77, 0.28])
    st = permutedims(data[:, on])
    imshow(st, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
           extent=(0.5, N + 0.5, lon[1], lon[end]),
           vmin=quantile(vec(st), 0.01), vmax=quantile(vec(st), darkness))
    ylabel("longitude (\$^\\circ\$)")
    tick_params(labelbottom=false)
    title(@sprintf("%s   L = %d   stride = %d", name_mod, L, sl.stride))

    # Panel 2: sliding spectra, each divided by its median over f ≥ fmin
    use = sl.freq .>= sl.fmin
    nrm = [median(sl.power[i, use]) for i in 1:size(sl.power, 1)]
    S = permutedims(sl.power ./ nrm)
    vmax = quantile(vec(S[use, :]), darkness)
    ax2 = fig.add_axes([0.20, 0.36, 0.77, 0.30], sharex=ax1)
    imshow(S, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
           extent=(sl.centers[1] - sl.stride / 2, sl.centers[end] + sl.stride / 2,
                   sl.freq[1], sl.freq[end]),
           vmin=0, vmax=vmax)
    axhline(y=sl.fmin, color="white", ls=":", lw=0.6)
    tick_params(labelleft=false, labelbottom=false)
    xlim(0.5, N + 0.5)

    ax2l = fig.add_axes([0.08, 0.36, 0.11, 0.30], sharey=ax2)
    plot(vec(mean(sl.power ./ nrm, dims=1)), sl.freq, color="grey")
    axhline(y=sl.fmin, color="black", ls=":", lw=0.6)
    ylim(sl.freq[1], sl.freq[end])
    xticks([])
    ylabel("frequency (1/P)")

    # Panel 3: P3 track
    ax3 = fig.add_axes([0.20, 0.07, 0.77, 0.27], sharex=ax1)
    if good === nothing
        good = tr.fit_ok .& (tr.snr_off .>= snr_min)
        crit = @sprintf("fit ok, S/N\$_{off}\$ ≥ %.0f", snr_min)
    else
        crit = threshold === nothing ? "good" :
               threshold isa Real ? @sprintf("fit ok, not at edge, contrast ≥ %.1f", threshold) :
               @sprintf("fit ok, not at edge, contrast ≥ local shuffle threshold (median %.1f)",
                        median(threshold))
    end
    bad = .!good
    any(bad) && plot(tr.centers[bad], tr.p3[bad], ".", ms=1.5, c="lightgrey", zorder=2)
    if segs !== nothing
        # sections coloured by group; good windows outside any section stay black
        insec = falses(length(good))
        for (si, sg) in enumerate(segs)
            col = "C$(mod(groups[si] - 1, 10))"
            w = sg.win
            insec[w] .= true
            errorbar(tr.centers[w], tr.p3[w], yerr=tr.p3_err[w], fmt=".", ms=1.5,
                     c=col, ecolor=col, elinewidth=0.3, capsize=0, zorder=3)
            ax1.axvspan(sg.first - 0.5, sg.last + 0.5, ymin=0.0, ymax=0.04, color=col, lw=0)
        end
        rest = good .& .!insec
        any(rest) && plot(tr.centers[rest], tr.p3[rest], ".", ms=1.5, c="black", zorder=3)
    elseif any(good)
        errorbar(tr.centers[good], tr.p3[good], yerr=tr.p3_err[good], fmt=".", ms=1.5,
                 c="black", ecolor="grey", elinewidth=0.3, capsize=0, zorder=3)
    end
    p3_ref === nothing || axhline(y=p3_ref, color="C3", ls="--", lw=0.6)
    if p3_lim !== nothing
        ylim(p3_lim...)
    elseif any(good)
        lo, hi = quantile(tr.p3[good], [0.02, 0.98])
        pad = 0.3 * (hi - lo) + 0.05 * hi
        ylim(max(0, lo - pad), hi + pad)
    end
    xlim(0.5, N + 0.5)
    xlabel("pulse number (window centre)")
    ylabel("P\$_3\$ (P)")
    minorticks_on()
    text(0.01, 0.95, @sprintf("%s: %d / %d", crit, count(good), length(good)),
         transform=ax3.transAxes, va="top", fontsize=6)
    if segs !== nothing
        for g in sort(unique(groups))
            sel = findall(==(g), groups)
            npl = sum(segs[i].npulse for i in sel)
            text(0.99, 0.95 - 0.08 * (g - 1),
                 @sprintf("group %d: P\$_3\$ ≈ %.2f, %d sections, %d pulses", g,
                          group_p3(segs, sel), length(sel), npl),
                 transform=ax3.transAxes, ha="right", va="top", fontsize=6,
                 color="C$(mod(g - 1, 10))")
        end
    end

    savepath = joinpath(outdir, "$(name_mod)_sliding_lrfs_L$(L)$(suffix).pdf")
    pdf && savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"), dpi=150)
    println(pdf ? savepath : replace(savepath, ".pdf" => ".png"))
    if show_
        PyPlot.show()
        println("Press Enter to close the figure.")
        readline(stdin; keep=false)
    end
    close()
end


"""
    plot_folds(res, outdir; nbin, name_mod, show_)

One column per group of `analyse` output `res`:
  row 1 — fold with variable-P3 compensation (`phase_fold`), two cycles;
  row 2 — constant-P3 fold of the same pulses (`constant_fold`), two cycles;
  row 3 — template |T| (grey) and phase ψ = −arg T (colour) vs longitude:
          a steady slope = drift, flat (or π steps) = amplitude modulation;
  row 4 — phase of every pulse against the constant-P3 phase,
          θ(n)/2π − n/P3 (unwrapped within sections): flat = constant P3,
          a slope or curvature = P3 wandering, steps between sections = phase
          jumps across nulls/gaps that a constant fold cannot follow.
Modulation depth (`modulation_depth`) is printed for both folds, with the
range from the shuffle control for the compensated one.

Writes `<name_mod>_p3fold_groups.pdf/.png`.
"""
function plot_folds(res, outdir; nbin=1024, name_mod="pulsar", darkness=0.99, pdf=true,
                    show_=false)
    ng = length(res.folds)
    ng == 0 && (println("no groups for $name_mod"); return)
    on = res.sl.on_bins
    lon = (collect(on) .- 1) .* 360.0 ./ nbin

    rc("font", size=7.)
    rc("axes", linewidth=0.5)
    rc("lines", linewidth=0.5)
    fig = figure(figsize=(2.6 * ng + 0.6, 9.0))
    for (k, (fo, cf)) in enumerate(zip(res.folds, res.cfolds))
        col = "C$(mod(fo.group - 1, 10))"
        for (row, F, lab) in ((1, fo.fold, "variable P\$_3\$"), (2, cf.fold, "constant P\$_3\$"))
            ax = subplot(4, ng, (row - 1) * ng + k)
            FF = vcat(F, F)
            imshow(FF, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
                   extent=(lon[1], lon[end], 0, 2), vmax=quantile(vec(FF), darkness))
            dep = row == 1 ?
                (isempty(fo.depth_null) ? @sprintf("%s  depth %.3f", lab, fo.depth) :
                 @sprintf("%s  depth %.3f (shuffle %.3f–%.3f)", lab, fo.depth,
                          minimum(fo.depth_null), maximum(fo.depth_null))) :
                @sprintf("%s  depth %.3f", lab, cf.depth)
            title(row == 1 ? @sprintf("group %d: P\$_3\$ ≈ %.2f, %d P, %d bins\n%s", fo.group,
                                      fo.p3, length(fo.pulses), fo.nb, dep) : dep,
                  fontsize=6, color=row == 1 ? col : "black")
            k == 1 && ylabel("P\$_3\$ phase (cycles)")
            row == 2 && xlabel("longitude (\$^\\circ\$)")
        end
        ax = subplot(4, ng, 2 * ng + k)
        tp = fo.tphase
        plot(lon, tp.amp, color="lightgrey", lw=0.8)
        ylim(0, 1.05)
        k == 1 && ylabel("|T| (norm.)")
        xlabel("longitude (\$^\\circ\$)")
        ax2 = ax.twinx()
        ax2.plot(lon, rad2deg.(tp.psi), ".", ms=2, c=col)
        ax2.set_ylabel("ψ = −arg T (\$^\\circ\$)")
        ts = fo.tsig
        title(@sprintf("Δψ %.2f cyc (%s)\nz = %.1f, upper %.2f → %s%s%s", tp.dpsi,
                       join([@sprintf("%+.2f±%.2f", d, e) for (d, e) in zip(tp.run_dpsi, ts.sigma_run)], ", "),
                       ts.z, ts.dpsi_upper, fo.verdict, fo.verdict_src == :fold ?
                       @sprintf(" (fold, z %.1f)", fo.tsig_fold.z) : "", fo.bidrift ? ", bi-drift" : ""),
              fontsize=6)
        ax = subplot(4, ng, 3 * ng + k)
        for i in fo.sections
            s = res.segs[i]
            idx = findall(n -> s.first ≤ n ≤ s.last, fo.pulses)
            n = fo.pulses[idx]
            r = unwrap_cycles(fo.theta[idx] ./ (2π) .- n ./ fo.p3)
            r .-= round(r[1])
            plot(n, r, ".", ms=1, c=col)
        end
        xlabel("pulse number")
        k == 1 && ylabel("θ/2π − n/P\$_3\$ (cycles)")
        minorticks_on()
    end
    suptitle(@sprintf("%s   L = %d", name_mod, res.L), fontsize=8)
    tight_layout(rect=(0, 0, 1, 0.96))
    savepath = joinpath(outdir, "$(name_mod)_p3fold_groups.pdf")
    pdf && savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"), dpi=150)
    println(pdf ? savepath : replace(savepath, ".pdf" => ".png"))
    if show_
        PyPlot.show()
        println("Press Enter to close the figure.")
        readline(stdin; keep=false)
    end
    close()
end

"""
    plot_summary(data, res, outdir; nbin, name_mod, p3_ref, darkness, show_)

All groups of both passes of `analyse` on one figure:
  panel 1 — single pulses with the pulse ranges of every group as strips
            (pass 1 along the bottom, pass 2 along the top);
  panel 2 — P3 of the windows of every section vs window centre, coloured by
            group (pass 1: dots, pass 2: crosses), each group labelled with
            its pass, L, P3, number of pulses and Δψ.
The per-pass figures (`plot_track`) show one window length each; this one
answers "which regimes were found, where, and are they drift or AM".

Writes `<name_mod>_p3track_summary.pdf/.png`.
"""
function plot_summary(data, res, outdir; nbin=size(data, 2), name_mod="pulsar", p3_ref=nothing,
                      darkness=0.995, pdf=true, show_=false)
    N = size(data, 1)
    on = res.sl.on_bins
    lon = (collect(on) .- 1) .* 360.0 ./ nbin
    passes = res.pass2 === nothing ? (res,) : (res, res.pass2)

    rc("font", size=7.)
    rc("axes", linewidth=0.5)
    rc("lines", linewidth=0.5)
    fig = figure(figsize=(7.0, 5.0))
    ax1 = fig.add_axes([0.10, 0.55, 0.87, 0.38])
    st = permutedims(data[:, on])
    imshow(st, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
           extent=(0.5, N + 0.5, lon[1], lon[end]),
           vmin=quantile(vec(st), 0.01), vmax=quantile(vec(st), darkness))
    ylabel("longitude (\$^\\circ\$)")
    tick_params(labelbottom=false)
    title(name_mod)
    ax2 = fig.add_axes([0.10, 0.09, 0.87, 0.44], sharex=ax1)
    ci = 0
    labels = String[]
    p3s = Float64[]
    for (ip, r) in enumerate(passes)
        for fo in r.folds
            col = "C$(mod(ci, 10))"
            ci += 1
            for i in fo.sections
                sg = r.segs[i]
                y0, y1 = ip == 1 ? (0.0, 0.04) : (0.96, 1.0)
                ax1.axvspan(sg.first - 0.5, sg.last + 0.5, ymin=y0, ymax=y1, color=col, lw=0)
                w = sg.win
                ax2.plot(r.tr.centers[w], r.tr.p3[w], ip == 1 ? "." : "x", ms=ip == 1 ? 1.5 : 2,
                         c=col)
                append!(p3s, r.tr.p3[w])
            end
            push!(labels, @sprintf("pass %d (L = %d): P\$_3\$ ≈ %.1f, %d P, Δψ = %.2f, z = %.1f → %s%s%s",
                                   ip, r.L, fo.p3, length(fo.pulses), fo.tphase.dpsi, fo.tsig.z,
                                   fo.verdict, fo.verdict_src == :fold ? " (fold)" : "",
                                   fo.bidrift ? ", bi-drift" : ""))
            ax2.text(0.99, 0.95 - 0.08 * (ci - 1), labels[end], transform=ax2.transAxes,
                     ha="right", va="top", fontsize=6, color=col)
        end
    end
    nq = hasproperty(res, :nyq) ? res.nyq : nothing
    if nq !== nothing && nq.found
        col = "C$(mod(ci, 10))"
        for r in nq.sections
            ax1.axvspan(first(r) - 0.5, last(r) + 0.5, ymin=0.48, ymax=0.52, color=col, lw=0)
            ax2.plot([first(r), last(r)], [nq.p3, nq.p3], "-", lw=2, c=col)
        end
        push!(p3s, nq.p3)
        lab = nq.resolved ? @sprintf("Nyquist path: P\$_3\$ ≈ %.3f (alias %.3f), %d P, |Δψ| = %.2f, z = %.1f → %s%s",
                                     nq.p3, nq.p3_alias, length(nq.pulses), nq.tphase.dpsi, nq.tsig.z,
                                     nq.verdict, nq.bidrift ? ", bi-drift" : "") :
              @sprintf("Nyquist path: %d P, aliases not resolved (M = %d) → %s", length(nq.pulses), nq.M, nq.verdict)
        ax2.text(0.99, 0.95 - 0.08 * ci, lab, transform=ax2.transAxes, ha="right", va="top", fontsize=6, color=col)
        ci += 1
    end
    p3_ref === nothing || ax2.axhline(y=p3_ref, color="grey", ls="--", lw=0.6)
    if !isempty(p3s)
        lo, hi = extrema(p3s)
        ax2.set_ylim(max(0, lo - 0.1 * (hi - lo) - 1), hi + 0.45 * (hi - lo) + 1)
    end
    ax2.set_xlim(0.5, N + 0.5)
    ax2.set_xlabel("pulse number (window centre)")
    ax2.set_ylabel("P\$_3\$ (P)")
    ax2.minorticks_on()
    savepath = joinpath(outdir, "$(name_mod)_p3track_summary.pdf")
    pdf && savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"), dpi=150)
    println(pdf ? savepath : replace(savepath, ".pdf" => ".png"))
    if show_
        PyPlot.show()
        println("Press Enter to close the figure.")
        readline(stdin; keep=false)
    end
    close()
end


"unwrap a phase sequence in cycles (jumps > 0.5 cycle folded back)"
function unwrap_cycles(x)
    y = copy(x)
    for i in 2:length(y)
        d = y[i] - y[i-1]
        y[i] -= round(d)
    end
    return y
end

end # module P3Track
