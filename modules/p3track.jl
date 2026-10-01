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

export sliding_lrfs, p3_track, contrast_null, good_windows, window_length, p3_segments,
       p3_groups, merge_sections, harmonic_groups, fundamental_track, select_groups,
       phase_fold, constant_fold, analyse, plot_track, plot_folds

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


"""
    harmonic_groups(segs, groups, L; tol=1.0) -> (groups, harm::BitVector)

Recognise groups that are the second harmonic of another, larger group
(|f3_a − 2·f3_b| ≤ tol/L, medians weighted by pulses) and relabel them to
that group: when the fundamental weakens, the strongest interior peak can
switch to 2·f3 (J0151-0635, L = 58: 40 pulses at P3 = 7.46 next to 14.3),
but the regime is the same. `harm[i]` marks sections whose track must be
halved in frequency (`fundamental_track`). Labels are re-ranked by size.
"""
function harmonic_groups(segs, groups, L::Int; tol::Real=1.0)
    g = copy(groups)
    harm = falses(length(segs))
    isempty(segs) && return g, harm
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
            for i in eachindex(segs)
                if g[i] == a
                    g[i] = b
                    harm[i] = true
                end
            end
        end
    end
    ids, npl, _ = stats(g)
    order = sort(ids, by=k -> -npl[k])
    rank = Dict(k => r for (r, k) in enumerate(order))
    return [rank[k] for k in g], harm
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
modulation depth of the real fold is reported against them.

Fields: group, pulses, p3 (group P3), nb, fold, counts, phase (θ/2π mod 1),
theta (radians, per pulse), template, f (per pulse), depth, depth_null
(vector), coherence (mean over pulses of |Σ_φ Z conj T| / (‖Z‖·‖T‖):
how well single pulses match the template shape), sections, on_bins.
"""
function phase_fold(data::AbstractMatrix, sl, tr, segs, groups, g::Int; nbins=nothing,
                    niter::Int=10, nshuffle::Int=5, seed::Int=1)
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
    depth = modulation_depth(F, cnt)
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
        depth_null[k] = modulation_depth(Fs, cs)
    end

    return (group=g, pulses=pulses, p3=p3g, nb=nb, fold=F, counts=cnt, phase=phase, theta=θ,
            template=T, f=f, depth=depth, depth_null=depth_null, coherence=coherence,
            on_bins=on, sections=sel)
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
    analyse(data, bin_st, bin_end, p3; window=nothing, ncycles=5, nshuffle=5)
        -> NamedTuple

Whole chain for one pulsar: window L = `window_length(p3)` (p3 from
params.json), sliding LRFS, P3 track, local shuffle threshold, sections,
groups (`p3_groups`), second harmonics joined to their fundamental
(`harmonic_groups`, `fundamental_track`), touching sections merged, groups
shorter than `ncycles`·P3 dropped, and a variable-P3 fold (`phase_fold`) plus
a constant-P3 reference fold (`constant_fold`) for every remaining group.

Fields: L, sl, tr (fundamental track), threshold, good, segs, groups, harm
(before merging), dropped, folds, cfolds.
"""
function analyse(data::AbstractMatrix, bin_st::Int, bin_end::Int, p3::Real; window=nothing,
                 ncycles::Real=5, nshuffle::Int=5)
    N = size(data, 1)
    L = window === nothing ? window_length(p3) : window
    sl = sliding_lrfs(data, bin_st, bin_end; window=L)
    tr = p3_track(sl)
    thr = contrast_null(data, sl).threshold
    good = good_windows(tr, thr)
    segs = p3_segments(tr, good, N)
    groups = p3_groups(segs, L)
    groups, harm = harmonic_groups(segs, groups, L)
    trf = fundamental_track(tr, segs, harm)
    segs, groups = merge_sections(trf, segs, groups, N)
    segs, groups, dropped = select_groups(segs, groups; ncycles=ncycles)
    folds = [phase_fold(data, sl, trf, segs, groups, g; nshuffle=nshuffle)
             for g in sort(unique(groups))]
    cfolds = [constant_fold(data, sl.on_bins, fo.pulses, fo.p3, fo.nb) for fo in folds]
    return (L=L, sl=sl, tr=trf, threshold=thr, good=good, segs=segs, groups=groups,
            harm=harm, dropped=dropped, folds=folds, cfolds=cfolds)
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
                    suffix="", show_=false)
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
    savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"), dpi=150)
    println(savepath)
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
  row 3 — phase of every pulse against the constant-P3 phase,
          θ(n)/2π − n/P3 (unwrapped within sections): flat = constant P3,
          a slope or curvature = P3 wandering, steps between sections = phase
          jumps across nulls/gaps that a constant fold cannot follow.
Modulation depth (`modulation_depth`) is printed for both folds, with the
range from the shuffle control for the compensated one.

Writes `<name_mod>_p3fold_groups.pdf/.png`.
"""
function plot_folds(res, outdir; nbin=1024, name_mod="pulsar", darkness=0.99, show_=false)
    ng = length(res.folds)
    ng == 0 && (println("no groups for $name_mod"); return)
    on = res.sl.on_bins
    lon = (collect(on) .- 1) .* 360.0 ./ nbin

    rc("font", size=7.)
    rc("axes", linewidth=0.5)
    rc("lines", linewidth=0.5)
    fig = figure(figsize=(2.6 * ng + 0.6, 7.0))
    for (k, (fo, cf)) in enumerate(zip(res.folds, res.cfolds))
        col = "C$(mod(fo.group - 1, 10))"
        for (row, F, lab) in ((1, fo.fold, "variable P\$_3\$"), (2, cf.fold, "constant P\$_3\$"))
            ax = subplot(3, ng, (row - 1) * ng + k)
            FF = vcat(F, F)
            imshow(FF, origin="lower", aspect="auto", cmap="viridis", interpolation="none",
                   extent=(lon[1], lon[end], 0, 2), vmax=quantile(vec(FF), darkness))
            dep = row == 1 ?
                @sprintf("%s  depth %.3f (shuffle %.3f–%.3f)", lab, fo.depth,
                         minimum(fo.depth_null), maximum(fo.depth_null)) :
                @sprintf("%s  depth %.3f", lab, cf.depth)
            title(row == 1 ? @sprintf("group %d: P\$_3\$ ≈ %.2f, %d P, %d bins\n%s", fo.group,
                                      fo.p3, length(fo.pulses), fo.nb, dep) : dep,
                  fontsize=6, color=row == 1 ? col : "black")
            k == 1 && ylabel("P\$_3\$ phase (cycles)")
            row == 2 && xlabel("longitude (\$^\\circ\$)")
        end
        ax = subplot(3, ng, 2 * ng + k)
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
    savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"), dpi=150)
    println(savepath)
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
