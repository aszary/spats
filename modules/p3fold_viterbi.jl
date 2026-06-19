module P3FoldViterbi

using Statistics

# Background / design rationale: p3fold-refine-notes.md §3.1-3.2.
#
# psrsalsa's `pfold` refine mode (foldP3 in fold.c) assigns phase block-by-block,
# greedily matching each block only against the sum of *previous* blocks
# (fold.c:158-176). A bad block can derail everything after it, and the
# correlation measure (template²·blockmap²) isn't normalized, so loud blocks
# dominate the phase choice.
#
# Here the phase is assigned per *pulse* (not per block) and globally via
# Viterbi: the whole observation is optimized at once, so one noisy pulse
# can't permanently knock later pulses off phase. The emission score is a
# normalized (zero-mean, unit-power) Pearson correlation, so faint and loud
# pulses are weighted by how well they match the template shape, not by their
# raw amplitude. The transition cost penalizes deviation from the phase
# increment implied by the nominal P3 but never forbids a jump outright, so
# genuine P3 changes can still be tracked (controlled by `continuity_weight`).


"""
    emission_corr(data, template, on) -> Matrix{Float64}

Per-pulse, per-phase-bin normalized correlation matrix (N_pulses × ybins).

`C[i, j]` is the zero-mean Pearson correlation between pulse `i`'s on-pulse
window (`on`) and phase-bin `j`'s on-pulse window in `template`. Pulses or
phase bins with zero variance (e.g. nulled pulses, still-empty template
bins) get correlation 0 against everything, so they neither win nor actively
mislead the Viterbi search.
"""
function emission_corr(data::AbstractMatrix, template::AbstractMatrix, on::UnitRange)
    ybins = size(template, 1)
    N = size(data, 1)

    tmpl_on = template[:, on]
    tmean = vec(mean(tmpl_on, dims=2))
    tc = tmpl_on .- tmean
    tnorm = [sqrt(sum(abs2, @view tc[j, :])) for j in 1:ybins]

    C = zeros(Float64, N, ybins)
    for i in 1:N
        d = @view data[i, on]
        dc = d .- mean(d)
        dn = sqrt(sum(abs2, dc))
        dn == 0 && continue
        for j in 1:ybins
            tnorm[j] == 0 && continue
            C[i, j] = sum(dc .* @view tc[j, :]) / (dn * tnorm[j])
        end
    end
    return C
end


"""
    circ_residual(x, m) -> Float64

Wrap `x` into (-m/2, m/2], i.e. the signed residual of `x` modulo `m`
closest to 0. Used to measure how far an actual phase step is from the
expected one, on a circle of `m` phase bins.
"""
function circ_residual(x::Real, m::Real)
    return mod(x + m / 2, m) - m / 2
end


"""
    viterbi(C, delta, continuity_weight) -> Vector{Int}

Globally optimal phase-bin path (1..ybins per pulse) given the emission
matrix `C` (N_pulses × ybins, see `emission_corr`).

Transition cost between consecutive pulses is
`continuity_weight * circ_residual(j - j', delta, ybins)^2`, i.e. a quadratic
penalty on deviating from the expected phase increment `delta = ybins / p3`.
`continuity_weight = 0` makes the search flat (any phase jump equally
"free", pure per-pulse template matching); larger values increasingly
enforce smooth, P3-like drift. This replaces foldP3's hard one-cycle-per-block
limit with a soft, tunable preference (notes §3.1).
"""
function viterbi(C::AbstractMatrix, delta::Real, continuity_weight::Real)
    N, ybins = size(C)
    cost = fill(Inf, N, ybins)
    backptr = zeros(Int, N, ybins)

    cost[1, :] .= -C[1, :]

    for i in 2:N
        for j in 1:ybins
            bestc = Inf
            bestp = 1
            for jp in 1:ybins
                resid = circ_residual((j - jp) - delta, ybins)
                trans = continuity_weight * resid^2
                c = cost[i-1, jp] + trans
                if c < bestc
                    bestc = c
                    bestp = jp
                end
            end
            cost[i, j] = bestc - C[i, j]
            backptr[i, j] = bestp
        end
    end

    path = zeros(Int, N)
    path[N] = argmin(@view cost[N, :])
    for i in N-1:-1:1
        path[i] = backptr[i+1, path[i+1]]
    end
    return path
end


"""
    build_template(data, path, ybins) -> Matrix{Float64}

Sum pulses into their assigned phase bin, ybins × N_bins (unnormalized sum,
same convention as `Tools.p3fold`).
"""
function build_template(data::AbstractMatrix, path::AbstractVector{Int}, ybins::Int)
    nb = size(data, 2)
    template = zeros(Float64, ybins, nb)
    for i in eachindex(path)
        template[path[i], :] .+= @view data[i, :]
    end
    return template
end


"""
    circ_unwrap_steps(path, ybins, delta) -> Vector{Float64}

Per-pulse phase increment `path[i+1] - path[i]`, unwrapped: each step is
resolved to the representative value closest to the expected nominal
increment `delta`, removing the ±ybins circular ambiguity. Length N-1.
"""
function circ_unwrap_steps(path::AbstractVector{Int}, ybins::Int, delta::Real)
    N = length(path)
    steps = zeros(Float64, N - 1)
    for i in 1:N-1
        raw = path[i+1] - path[i]
        steps[i] = delta + circ_residual(raw - delta, ybins)
    end
    return steps
end


"""
    instantaneous_p3(path, ybins, p3_nominal; window=20) -> Vector{Float64}

Smoothed, per-pulse instantaneous P3 [pulse periods] implied by the Viterbi
phase path — this is "the P3 actually used" for each pulse, as opposed to
the single nominal `p3` fed into `fold`.

The discrete per-pulse phase increments are unwrapped into a continuous
phase track (`circ_unwrap_steps`), then a local linear fit over a sliding
window of `window` pulses gives the local slope dphase/dpulse, from which
`P3 = ybins / slope`. Pulses too close to either end for a full window keep
`p3_nominal`.
"""
function instantaneous_p3(path::AbstractVector{Int}, ybins::Int, p3_nominal::Real; window::Int=20)
    delta = ybins / p3_nominal
    N = length(path)
    steps = circ_unwrap_steps(path, ybins, delta)
    cum = vcat(0.0, cumsum(steps))   # cum[i] = unwrapped phase relative to pulse 1

    p3_inst = fill(Float64(p3_nominal), N)
    halfw = max(1, window ÷ 2)
    for i in 1:N
        lo = max(1, i - halfw)
        hi = min(N, i + halfw)
        hi - lo < 2 && continue
        xs = lo:hi
        ys = @view cum[lo:hi]
        xbar = mean(xs)
        ybar = mean(ys)
        num = sum((x - xbar) * (y - ybar) for (x, y) in zip(xs, ys))
        den = sum((x - xbar)^2 for x in xs)
        slope = den == 0 ? delta : num / den
        p3_inst[i] = slope == 0 ? p3_nominal : ybins / slope
    end
    return p3_inst
end


"""
    fold(data, p3, bin_st, bin_end; ybins, n_iter, continuity_weight) -> NamedTuple

Refined P3-fold via per-pulse Viterbi phase assignment, intended as a
globally-optimized alternative to `pfold -p3fold`'s block-greedy refine
(p3fold-refine-notes.md §3).

Bootstraps from the same fixed-phase fold as `Tools.p3fold` (so `n_iter=0`
reproduces it exactly), then alternates:
  1. score every pulse against every phase bin of the current template
     (`emission_corr`),
  2. find the globally optimal phase path (`viterbi`),
  3. rebuild the template by summing pulses along that path,
for `n_iter` EM-like rounds.

Arguments:
  data     – single-pulse matrix (N_pulses × N_bins), real intensity
  p3       – nominal P3 [pulse periods P0]; only used to seed the bootstrap
             fold and to set the expected phase increment per pulse
             (`delta = ybins / p3`) for the continuity penalty
  bin_st, bin_end – on-pulse window (1-indexed) used for the correlation
             score; the folded output still spans all N_bins
  ybins    – number of P3-phase bins (states), default 10
  n_iter   – number of refit rounds, default 5
  continuity_weight – transition penalty strength; 0 = flat (any phase
             jump equally likely, pure template matching), default 0.2
  p3_window – smoothing window [pulses] for `p3_per_pulse`, default 20

Returns:
  folded       – ybins × N_bins matrix, the refined p3-fold
  phase        – Vector{Int}, assigned phase bin (1..ybins) per pulse
  confidence   – Vector{Float64}, correlation of each pulse against its
                 assigned phase bin (low ⇒ poor match, e.g. nulled pulse)
  margin       – Vector{Float64}, correlation gap between the best and
                 second-best phase bin per pulse (low ⇒ ambiguous phase,
                 notes §3.4)
  p3_per_pulse – Vector{Float64}, smoothed instantaneous P3 [pulse periods]
                 implied by `phase` (see `instantaneous_p3`) — the local P3
                 the Viterbi path is actually tracking at each pulse
"""
function fold(data::AbstractMatrix, p3::Real, bin_st::Int, bin_end::Int;
              ybins::Int=10, n_iter::Int=5, continuity_weight::Real=0.2,
              p3_window::Int=20)
    N = size(data, 1)
    on = bin_st:bin_end
    delta = ybins / p3

    path = [floor(Int, mod(i * delta, ybins)) + 1 for i in 1:N]
    template = build_template(data, path, ybins)

    C = zeros(Float64, N, ybins)
    for _ in 1:n_iter
        C = emission_corr(data, template, on)
        path = viterbi(C, delta, continuity_weight)
        template = build_template(data, path, ybins)
    end

    confidence = zeros(Float64, N)
    margin = zeros(Float64, N)
    if n_iter > 0
        for i in 1:N
            row = @view C[i, :]
            confidence[i] = row[path[i]]
            best, second = partialsort(row, 1:2, rev=true)
            margin[i] = best - second
        end
    end

    p3_per_pulse = instantaneous_p3(path, ybins, p3; window=p3_window)

    return (folded=template, phase=path, confidence=confidence, margin=margin,
            p3_per_pulse=p3_per_pulse)
end

end # module P3FoldViterbi
