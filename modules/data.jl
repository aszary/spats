module Data
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


    """
    Zero selected pulses (in place).

    `ranges` is the same `paz -w`-style string stored in `p["zaps"]`
    (e.g. `"120-135 640-705"` or a lone `"42"`), space-separated tokens of
    either a single 0-based subint number or a `start-end` 0-based,
    inclusive range — matching the `-Z 0-12` channel-zap syntax already
    used elsewhere in this file. `nothing` (no zaps defined) is a no-op.
    """
    function zap!(data; ranges=nothing)
        isnothing(ranges) && return data
        pulses = size(data, 1)
        for tok in split(ranges)
            bounds = split(tok, "-")
            lo = parse(Int, bounds[1]) + 1  # PSRCHIVE subints are 0-based, Julia rows are 1-based
            hi = length(bounds) == 1 ? lo : parse(Int, bounds[2]) + 1
            if lo < 1 || hi > pulses
                @warn "zap! range $tok (1-based $lo:$hi) outside data bounds 1:$pulses, clamping"
            end
            for z in max(lo, 1):min(hi, pulses)
                data[z, :] .= 0
            end
        end
    end


    """
        zaps_to_ranges(zaps) -> Vector{Tuple{Int,Int}}

    Parse a `p["zaps"]` string ("12 40-45", 0-based, inclusive) into 0-based
    (start, end) pairs. `nothing`/empty -> empty vector.
    """
    function zaps_to_ranges(zaps)
        (isnothing(zaps) || isempty(strip(zaps))) && return Tuple{Int,Int}[]
        out = Tuple{Int,Int}[]
        for tok in split(zaps)
            b = split(tok, "-")
            lo = parse(Int, b[1])
            push!(out, (lo, length(b) == 1 ? lo : parse(Int, b[2])))
        end
        return out
    end


    """
        pulses_to_zaps(rows) -> String

    1-based row numbers (Julia) -> `p["zaps"]` string (0-based PSRCHIVE
    subints), consecutive pulses merged into "a-b" ranges.
    """
    function pulses_to_zaps(rows)
        r = sort(unique(rows)) .- 1
        isempty(r) && return ""
        toks = String[]
        st = prev = r[1]
        for x in r[2:end]
            if x == prev + 1
                prev = x
                continue
            end
            push!(toks, st == prev ? "$st" : "$st-$prev")
            st = prev = x
        end
        push!(toks, st == prev ? "$st" : "$st-$prev")
        return join(toks, " ")
    end


    """
        find_bad_pulses(data, bin_st, bin_end; nsigma=5.0, margin=nothing)

    Automatic bad-pulse (RFI / broken data) detection, from the OFF-pulse
    region only — so nulls and bright pulses (real pulsar emission) are never
    rejected, only pulses where the noise itself is wrong:

      * `z_rms`  – off-pulse noise level vs its LOCAL running median
                   (`local_win` pulses) — sudden jumps (RFI bursts, dropouts)
                   flagged when |z| > nsigma; slow gain/scintillation drifts
                   are not (they raise signal and noise together)
      * `z_struct` – baseline structure: spread of off-pulse block means in
                   units of white noise (wavy / tilted / stepped baselines =
                   horizontal stripes in the image) — flagged when z > nsigma
      * `z_peak` – largest single off-pulse sample in units of that pulse's own
                   rms (narrow spikes) — flagged when z > nsigma AND the spike
                   is above the Gaussian extreme-value level (no false alarms
                   from the max of many noise samples)
      * all-zero pulses (already zapped / missing) are always flagged

    Baseline-structure and spike pulses are always cut (systematics). A
    noise-level outlier is only a SUSPECT: it is cut only if it adds more
    noise than signal (s/v < S/2V, see code), so bright pulses whose
    off-pulse merely looks noisier are kept (`kept_bright`).

    z = robust z-score, (x - median) / (1.4826·MAD), so the bad pulses
    themselves don't inflate the threshold. Flagged pulses separated by at
    most `fill_gap` good ones are merged into one bad region.

    Returns (bad = 1-based rows, z_rms, z_peak, z_struct, on_energy).
    """
    function find_bad_pulses(data, bin_st, bin_end; nsigma=5.0, margin=0, nblocks=8,
                             local_win=31, fill_gap=2)
        N, nb = size(data)
        off = vcat(1:max(0, bin_st - margin - 1), min(nb + 1, bin_end + margin + 1):nb)
        length(off) < 20 && (off = vcat(1:bin_st-1, bin_end+1:nb))
        off_rms  = [std(@view data[i, off]) for i in 1:N]
        off_peak = [maximum(abs, @view data[i, off]) / max(off_rms[i], eps()) for i in 1:N]
        on_energy = [sum(@view data[i, bin_st:bin_end]) for i in 1:N]
        # baseline structure: spread of the mean of `nblocks` off-pulse blocks,
        # in units of what white noise would give (≈1). Wavy / tilted /
        # stepped baselines (the horizontal stripes in the single-pulse image)
        # push it up. Ratio of two noise-scaled numbers -> gain-independent.
        bs = cld(length(off), nblocks)
        blocks = [off[k:min(k + bs - 1, end)] for k in 1:bs:length(off)]
        off_struct = [std([mean(@view data[i, b]) for b in blocks]) /
                      max(off_rms[i] / sqrt(bs), eps()) for i in 1:N]

        empty = off_rms .== 0
        good = .!empty
        function robust_z(x)
            z = fill(NaN, N)
            xs = x[good]
            isempty(xs) && return z
            med = median(xs)
            mad = 1.4826 * median(abs.(xs .- med))
            mad = mad > 0 ? mad : max(std(xs), eps())
            z[good] = (xs .- med) ./ mad
            return z
        end
        # noise level relative to its LOCAL running median: slow changes of
        # gain / scintillation raise noise AND signal together (seen on
        # J1110: whole 330-760 stretch elevated, flagging it threw away the
        # brightest pulses and LOWERED the profile S/N) — only sudden jumps
        # of the noise level are bad pulses
        lr = log.(max.(off_rms, eps()))
        gi = findall(good)
        h = local_win ÷ 2
        loc = zeros(N)
        for i in 1:N
            win = filter(j -> abs(j - i) <= h, gi)
            loc[i] = isempty(win) ? 0.0 : median(lr[win])
        end
        z_rms = robust_z(lr .- loc)
        z_peak = robust_z(off_peak)
        z_struct = robust_z(off_struct)

        # the max of many noise samples has a heavy right tail, so z_peak > nsigma
        # alone gives false alarms; also require the spike to exceed what pure
        # Gaussian noise practically never reaches anywhere in the observation
        # (extreme-value level for N·n_off samples, +1σ margin)
        k_spike = sqrt(2 * log(N * length(off))) + 1.0
        spike(i) = z_peak[i] > nsigma && off_peak[i] > k_spike
        # baseline structure and spikes are systematics, not white noise —
        # they don't average down, so they are always cut
        always = [spike(i) || z_struct[i] > nsigma for i in 1:N]
        suspect = [abs(z_rms[i]) > nsigma for i in 1:N]

        # S/N gate for noise-level suspects: cut only if the pulse brings more
        # noise than signal. Adding a pulse with signal s and variance v to
        # the sum raises S/N iff s/v > S/(2V) (S, V = totals). Bright pulses
        # whose off-pulse merely looks noisier (per-pulse digitiser scaling —
        # seen on real data: zapping them LOWERED S/N) are kept; pulses
        # carrying mostly noise are cut. Nulls are not suspects unless their
        # noise is abnormal, so they stay.
        n_on = bin_end - bin_st + 1
        sig = [on_energy[i] - n_on * mean(@view data[i, off]) for i in 1:N]
        noise_var = n_on .* off_rms .^ 2
        S, V = sum(sig[good]), sum(noise_var[good])
        worth(i) = V > 0 && noise_var[i] > 0 && sig[i] / noise_var[i] > 0.5 * S / V
        flag = [empty[i] || always[i] || (suspect[i] && !worth(i)) for i in 1:N]

        # bad stretches come in runs: close gaps of <= fill_gap pulses
        # between flagged pulses, so a bad region is cut as one block
        fl = findall(flag)
        for k in 1:length(fl)-1
            if 1 < fl[k+1] - fl[k] <= fill_gap + 1
                flag[fl[k]+1:fl[k+1]-1] .= true
            end
        end
        bad = findall(flag)
        return (bad=bad, z_rms=z_rms, z_peak=z_peak, z_struct=z_struct,
                on_energy=on_energy, nsigma=nsigma, sig=sig, noise_var=noise_var,
                kept_bright=findall(i -> suspect[i] && worth(i) && !always[i], 1:N))
    end


    """
        pulse_features(data, bin_st, bin_end) -> Vector{Pair{String,Vector{Float64}}}

    Many candidate per-pulse quality features (off-pulse = all bins outside
    bin_st:bin_end), computed so that hand labels (`label_bad_pulses` in
    spats.jl) can show which ones actually separate pulses a human would cut.
    Deliberately broad — most will turn out irrelevant:

      off_rms        off-pulse noise level
      off_rms_local  off_rms / running median of off_rms (31 pulses)
      off_mean       off-pulse mean (≈0 after debase)
      off_tilt       mean(right off) - mean(left off), in noise units
      off_struct     spread of 8 off-pulse block means / white-noise expectation
      off_lowfreq    fraction of off-pulse power in the lowest 5% of Fourier freqs
      off_ac1        lag-1 autocorrelation of off-pulse samples
      off_peak       max |off| / off_rms
      off_kurt       excess kurtosis of off-pulse samples
      off_skew       skewness of off-pulse samples
      off_mad_ratio  1.4826·MAD / rms (≈1 Gaussian; <1 heavy tails/outliers)
      off_nunique    fraction of distinct values (low = clipping / quantisation)
      on_energy      on-pulse sum
      on_snr         on-pulse sum / (off_rms·√n_on)
      on_peak        max on-pulse / off_rms
      row_rms_local  whole-row rms / running median of it
    """
    function pulse_features(data, bin_st, bin_end)
        N, nb = size(data)
        off = vcat(1:bin_st-1, bin_end+1:nb)
        on = bin_st:bin_end
        lft, rgt = 1:bin_st-1, bin_end+1:nb
        n_on = length(on)
        runmed(x, h=15) = [median(x[max(1, i - h):min(N, i + h)]) for i in 1:N]
        safe(a, b) = b > 0 ? a / b : NaN

        off_rms  = [std(@view data[i, off]) for i in 1:N]
        off_mean = [mean(@view data[i, off]) for i in 1:N]
        tilt = [(isempty(rgt) || isempty(lft)) ? NaN :
                safe(mean(@view data[i, rgt]) - mean(@view data[i, lft]), off_rms[i]) for i in 1:N]
        bs = cld(length(off), 8)
        blocks = [off[k:min(k + bs - 1, end)] for k in 1:bs:length(off)]
        struct_ = [safe(std([mean(@view data[i, b]) for b in blocks]), off_rms[i] / sqrt(bs)) for i in 1:N]
        nlow = max(1, round(Int, 0.05 * length(off) / 2))
        lowfreq = map(1:N) do i
            x = data[i, off] .- off_mean[i]
            P = abs2.(fft(x))[2:length(x)÷2+1]
            safe(sum(P[1:nlow]), sum(P))
        end
        ac1 = map(1:N) do i
            x = data[i, off] .- off_mean[i]
            safe(sum(x[1:end-1] .* x[2:end]), sum(abs2, x))
        end
        peak = [safe(maximum(abs, @view data[i, off]), off_rms[i]) for i in 1:N]
        moment(i, k) = (x = (data[i, off] .- off_mean[i]); mean(x .^ k))
        kurt = [safe(moment(i, 4), moment(i, 2)^2) - 3 for i in 1:N]
        skew = [safe(moment(i, 3), moment(i, 2)^1.5) for i in 1:N]
        madr = [safe(1.4826 * median(abs.(data[i, off] .- median(data[i, off]))), off_rms[i]) for i in 1:N]
        nuniq = [length(unique(@view data[i, off])) / length(off) for i in 1:N]
        on_energy = [sum(@view data[i, on]) - n_on * off_mean[i] for i in 1:N]
        row_rms = [std(@view data[i, :]) for i in 1:N]

        return ["off_rms" => off_rms,
                "off_rms_local" => off_rms ./ runmed(off_rms),
                "off_mean" => off_mean,
                "off_tilt" => tilt,
                "off_struct" => struct_,
                "off_lowfreq" => lowfreq,
                "off_ac1" => ac1,
                "off_peak" => peak,
                "off_kurt" => kurt,
                "off_skew" => skew,
                "off_mad_ratio" => madr,
                "off_nunique" => nuniq,
                "on_energy" => on_energy,
                "on_snr" => [safe(on_energy[i], off_rms[i] * sqrt(n_on)) for i in 1:N],
                "on_peak" => [safe(maximum(@view data[i, on]), off_rms[i]) for i in 1:N],
                "row_rms_local" => row_rms ./ runmed(row_rms)]
    end


    "Write `pulse_features` + label columns to a CSV (one row per pulse, 0-based pulse number)."
    function write_labels_csv(file, feats, labels::Pair...)
        N = length(last(first(feats)))
        mkpath(dirname(file))
        open(file, "w") do f
            println(f, join(vcat("pulse", first.(feats), [String(first(l)) for l in labels]), ","))
            for i in 1:N
                vals = vcat(string(i - 1), [string(round(v[i], sigdigits=6)) for v in last.(feats)],
                            [string(Int(last(l)[i])) for l in labels])
                println(f, join(vals, ","))
            end
        end
        println("Zapisano: $file")
    end


    """
        parse_paz_command(out) -> (subints, chans, cmd)

    Reads the `paz ...` line(s) that pazi prints on exit, e.g.
      paz -w "267 268 269" -e pazi pulsar.spCF
    Subints (0-based) from `-w "list"` / `-W "a b"`, channels from
    `-z "list"` / `-Z "a b"`; all paz lines are merged.
    """
    function parse_paz_command(out)
        cmds = [String(strip(l)) for l in split(out, '\n') if startswith(strip(l), "paz ")]
        subints, chans = Int[], Int[]
        arg(flag, c) = [something(m.captures[1], m.captures[2])
                        for m in eachmatch(Regex("$flag\\s+(?:\"([^\"]*)\"|(\\S+))"), c)]
        for c in cmds
            for s in arg("-w", c); append!(subints, parse.(Int, split(s))); end
            for s in arg("-z", c); append!(chans, parse.(Int, split(s))); end
            for s in arg("-W", c); a, b = parse.(Int, split(s)); append!(subints, a:b); end
            for s in arg("-Z", c); a, b = parse.(Int, split(s)); append!(chans, a:b); end
        end
        return (subints=sort(unique(subints)), chans=sort(unique(chans)),
                cmd=isempty(cmds) ? nothing : join(cmds, " ; "))
    end


    """
        run_pazi(archive) -> (subints, chans, cmd)

    Opens PSRCHIVE `pazi` on `archive`; zap pulses/channels there as usual
    and quit. The `paz` command pazi prints is captured and parsed — the
    archive itself is not modified by this function.
    """
    function run_pazi(archive)
        local out
        while true
            tmp = tempname()
            println("\npazi: zapuj impulsy (i ew. kanały) jak zwykle, potem wyjdź z pazi.\n" *
                    "Komenda paz zostanie odczytana automatycznie.")
            proc = run(pipeline(ignorestatus(`pazi $archive`), `tee $tmp`))
            out = read(tmp, String)
            rm(tmp, force=true)
            # pazi crashes e.g. when a click lands outside the subint range
            # (J0304: "isubint=1032 nsubint=1032") — before printing paz
            crashed = occursin("Error::", out) || any(p -> p.exitcode != 0, proc.processes)
            crashed || break
            print("pazi zakończył się błędem (np. kliknięcie poza zakresem) — zapy z tej sesji przepadły.\n" *
                  "Enter = otwórz pazi jeszcze raz, n = przejdź dalej bez zapów: ")
            strip(readline()) == "n" && break
        end
        z = parse_paz_command(out)
        println(isnothing(z.cmd) ? "pazi: nic nie wycięto." :
                "pazi: $(length(z.subints)) impulsów, $(length(z.chans)) kanałów  ←  $(z.cmd)")
        return z
    end


    "Features + label (+ auto) CSV for detector training, in repo input/zap_labels/<psr>.csv."
    function save_zap_labels(psr, data, bin_st, bin_end, rows)
        N = size(data, 1)
        label = falses(N); label[filter(i -> 1 <= i <= N, rows)] .= true
        auto = falses(N); auto[auto_zap(data, bin_st, bin_end).bad] .= true
        write_labels_csv(joinpath(@__DIR__, "..", "input", "zap_labels", "$psr.csv"),
                         pulse_features(data, bin_st, bin_end), "label" => label, "auto" => auto)
        tp, fp, fn = count(label .& auto), count(.!label .& auto), count(label .& .!auto)
        println("Ty: $(count(label))  automat: $(count(auto))  zgodne: $tp  " *
                "automat wyciął, Ty nie: $fp  Ty wyciąłeś, automat nie: $fn")
    end


    """
        offpulse_bins(data, bin_st, bin_end; nsig=5.0, pad=10)

    Off-pulse bins: outside the on-pulse window AND without emission in the
    mean profile (boxcar-smoothed, > `nsig`·σ above the median baseline),
    widened by `pad` bins. An interpulse or a wide component outside the
    window (J1057-5226) was counted as noise: profile S/N 21 instead of
    thousands and 1295 pulses cut as "noisy". Falls back to the whole
    outside if less than 20% would be left.
    """
    function offpulse_bins(data, bin_st, bin_end; nsig=5.0, pad=10)
        N, nb = size(data)
        off0 = vcat(1:bin_st-1, bin_end+1:nb)
        keep = [any(!iszero, @view data[i, :]) for i in 1:N]
        (count(keep) == 0 || length(off0) < 10) && return off0
        P = vec(mean(@view(data[keep, :]), dims=1))
        μ = median(P[off0])
        σ = 1.4826 * median(abs.(P[off0] .- μ))
        σ > 0 || return off0
        h = 4
        Ps = [mean(P[mod1(j + k, nb)] for k in -h:h) for j in 1:nb]
        em = (Ps .- μ) .> nsig * σ / sqrt(2h + 1)
        emx = falses(nb)
        for j in findall(em), k in -pad:pad; emx[mod1(j + k, nb)] = true; end
        off = filter(j -> !emx[j], off0)
        return length(off) >= 0.2 * length(off0) ? off : off0
    end


    """
        pulse_signal_noise(data, bin_st, bin_end) -> (s, v, good)

    Per pulse: on-pulse signal `s` (sum minus off-pulse baseline) and the
    variance `v` of that sum due to noise, measured EMPIRICALLY: the
    off-pulse is cut into non-overlapping windows as wide as the on-pulse
    window and the spread of their sums is taken. Unlike n_on·σ², this also
    includes correlated junk — baseline stripes, wiggles, contamination
    right next to the pulse — exactly what makes a stretch useless.
    Off-pulse windows carrying pulsar emission (leak of bright pulsars) are
    detected and excluded (see code; `nleak` in the result).
    With only a handful of windows per pulse `v` is noisy for ONE pulse but
    unbiased, so it is fine in sums over stretches; `vwhite` = n_on·σ² is
    the precise white-noise part, used to spot single noise jumps.
    Returns (s, v, vwhite, good, nwin); `good` = pulse not empty.
    """
    function pulse_signal_noise(data, bin_st, bin_end)
        N, nb = size(data)
        on = bin_st:bin_end
        n_on = length(on)
        off = offpulse_bins(data, bin_st, bin_end)
        segs = UnitRange{Int}[]
        for j in off
            (!isempty(segs) && last(segs[end]) == j - 1) ? (segs[end] = first(segs[end]):j) : push!(segs, j:j)
        end
        wins = UnitRange{Int}[]
        for sg in segs, j in 0:(length(sg) ÷ n_on)-1
            push!(wins, first(sg) + j * n_on : first(sg) + (j + 1) * n_on - 1)
        end
        good = [var(@view data[i, off]) > 0 for i in 1:N]
        base = [mean(@view data[i, off]) for i in 1:N]
        # MATCHED FILTER: on-pulse bins weighted by the mean profile (smoothed,
        # >= 0, mean weight 1). A plain sum over a wide on-pulse window
        # (J0820-4114: 452 bins, mostly empty edges) drowns a weaker pulse in
        # the noise of the empty bins, so it looked like a null and cutting it
        # LOWERED the S/N pazi shows (pazi uses its optimal width). With the
        # weights the empty edges hardly count — like pazi's optimal width.
        gi0 = findall(good)
        prof = isempty(gi0) ? zeros(n_on) :
               vec(mean(data[gi0, on] .- base[gi0], dims=1))
        hs = max(1, n_on ÷ 40)
        T = [max(0.0, mean(@view prof[max(1, j - hs):min(n_on, j + hs)])) for j in 1:n_on]
        T = sum(T) > 0 ? T .* (n_on / sum(T)) : ones(n_on)
        wsum(i, bins) = sum(T[j] * (data[i, bins[j]] - base[i]) for j in 1:n_on)
        s = [wsum(i, on) for i in 1:N]
        # Windows that contain PULSAR emission (bright pulsars leak outside
        # the on-pulse window: J0304, S/N 6392, had every pulse flagged)
        # must not count as noise: their sums follow the pulse brightness
        # (correlate with s) or sit at a steady non-zero level. Junk next to
        # the pulse (RFI, stripes) does neither, so it is still caught.
        W = [wsum(i, w) for i in 1:N, w in wins]
        gi = findall(good)
        ng = length(gi)
        leak = falses(length(wins))
        # rank (Spearman) correlation: a few very loud junk pulses (big s AND
        # big W) must not make a clean window look like it follows the pulse
        rank(x) = (r = similar(x); r[sortperm(x)] = 1:length(x); r)
        rs = rank(s[gi])
        for k in eachindex(wins)
            x = W[gi, k]
            md = median(x)
            spread = 1.4826 * median(abs.(x .- md))
            steady = spread > 0 && abs(md) > 4 * spread / sqrt(ng)
            follows = std(x) > 0 && std(s[gi]) > 0 && cor(rank(x), rs) > 4 / sqrt(ng)
            leak[k] = steady || follows
        end
        kw = findall(.!leak)
        clean_off = isempty(kw) ? off : vcat([collect(wins[k]) for k in kw]...)
        T2 = sum(abs2, T)
        vw = [T2 * var(@view data[i, clean_off]) for i in 1:N]
        v = length(kw) >= 3 ? [mean(abs2, @view W[i, kw]) for i in 1:N] : copy(vw)
        # lag-1 autocorrelation of the clean off-pulse (within windows):
        # stripes / baseline wiggles / slopes make neighbouring bins alike.
        # Works with any number of windows (J1750: wide on-pulse, only 2
        # windows -> v = white noise, blind to stripes; its stripes had
        # ac1 ~0.3 vs typical 0.02 ± 0.036)
        # The static off-pulse shape and emission leaking with the pulse
        # brightness (bright J0304-like wings: smooth, so ac1 would follow
        # the pulse) are regressed out first: x = d - base - c - s·m.
        segs_ac = isempty(kw) ? segs : wins[kw]
        ac_bins = vcat([collect(sg) for sg in segs_ac]...)
        R = [data[i, j] - base[i] for i in 1:N, j in ac_bins]
        sc = s .- mean(s[gi])
        ss = sum(abs2, sc[gi])
        for (c, j) in enumerate(ac_bins)
            col = @view R[:, c]
            col .-= mean(col[gi])
            ss > 0 && (col .-= sc .* (sum(sc[gi] .* col[gi]) / ss))
        end
        ac1 = zeros(N)
        for i in 1:N
            num = den = 0.0
            c0 = 0
            for sg in segs_ac
                L = length(sg)
                for j in 1:L-1
                    num += R[i, c0+j] * R[i, c0+j+1]
                end
                for j in 1:L
                    den += R[i, c0+j]^2
                end
                c0 += L
            end
            ac1[i] = den > 0 ? num / den : 0.0
        end
        return (s=s, v=v, vwhite=vw, good=good, nwin=length(kw), nleak=count(leak), ac1=ac1)
    end


    """
    Robust z of each pulse's log noise variance against its running median
    (`local_win` pulses): sudden noise jumps, not slow gain drifts. Doesn't
    depend on the tuned parameters, so `tune_autozap` computes it once.
    """
    function noise_jump_z(v, good; local_win=31)
        N = length(v)
        lv = log.(max.(v, eps()))
        gi = findall(good)
        h = local_win ÷ 2
        loc = [isempty(gi) ? 0.0 : median(lv[filter(j -> abs(j - i) <= h, gi)]) for i in 1:N]
        res = lv .- loc
        rg = res[good]
        mad = isempty(rg) ? 1.0 : max(1.4826 * median(abs.(rg .- median(rg))), eps())
        return res ./ mad
    end


    """
        auto_zap_core(sn; ...) -> NamedTuple

    Automatic JUNK zapping (built from hand-zapped examples and their
    failures on real data). It cuts junk, not merely weaker emission:

      1. single junk pulses: correlated off-pulse junk (stripes, wiggles,
         RFI next to the pulse) = window variance / white variance far above
         its TYPICAL value in this observation (robust z > `nsigma`; real
         data have correlated baseline structure in every pulse, so the
         white-noise expectation flagged ~half of all pulses — J0034, J0304)
         -> always cut; a sudden jump of the white noise (robust z vs a
         running median > `nsigma`) -> cut only if the pulse brings more
         noise than signal (bright pulses with digitiser-scaled noise stay);
      2. null / no-signal stretches: over a sliding window (singles left
         out) the signal is significantly (`zsig`·σ) below `f_null` of the
         typical on-state level (quantile `q_on` of the window efficiencies)
         — nulls and faded stretches, but NOT ordinary fainter emission
         (cutting that raised S/N on clean J0304 but is cherry-picking);
      3. noisy stretches: mean noise of the window > `noise_fac` × typical;
      4. stretches shorter than `min_run` are ignored (pulse-to-pulse
         modulation); kept islands < `min_keep` between cuts are cut; every
         cut region gets `pad` pulses on each side;
      5. safety: stretch cuts are dropped if they would lower the S/N.

    Window length adapts to the pulsar (average window S/N ≈ `target`).
    Returns (bad, single, stretch, reason, eff, E, w, snr_before,
    snr_after, why, n_single, n_cand). reason per pulse: 0 kept, 1 empty,
    2 single junk, 3 null/no signal, 4 noisy stretch, 5 margin/fill.
    """
    function auto_zap_core(sn; target=8.0, w=nothing, zsig=3.0, pad=1, min_keep=3, min_run=5,
                           nsigma=6.0, f_null=0.25, q_on=0.9, noise_fac=2.0, zdet=2.0,
                           f_rel=0.15, zrel=1.0, ac_sigma=4.5, min_run_fade=30,
                           local_win=31, noise_z=nothing, kw...)
        s, v, good = sn.s, sn.v, sn.good
        N = length(s)
        gi = findall(good)
        isempty(gi) && return (bad=Int[], single=Int[], stretch=Int[], reason=zeros(Int8, N),
                               eff=zeros(N), E=NaN, w=0, snr_before=NaN, snr_after=NaN,
                               why="brak danych", n_single=0, n_cand=0)
        robz(x) = (m = median(x[gi]); d = max(1.4826 * median(abs.(x[gi] .- m)), eps()); (x .- m) ./ d)

        # 1. single junk
        nz = isnothing(noise_z) ? noise_jump_z(sn.vwhite, good; local_win=local_win) : noise_z
        zc = sn.nwin >= 3 ? robz(log.(max.(v, eps()) ./ max.(sn.vwhite, eps()))) : zeros(N)
        S0, V0 = sum(s[gi]), sum(v[gi])
        # noise of ONE pulse: v from a few windows is a rough estimate (a
        # pulse with 25x white noise came out at 11x and was "worth"
        # keeping); the precise white noise, scaled to v, bounds it below
        kvw = median(v[gi] ./ max.(sn.vwhite[gi], eps()))
        vone(i) = max(v[i], kvw * sn.vwhite[i])
        worth(i) = V0 > 0 && vone(i) > 0 && s[i] / vone(i) > 0.5 * S0 / V0
        # stripes / wiggles seen as off-pulse autocorrelation far above its
        # typical value here (also when there are too few windows for zc)
        zac = hasproperty(sn, :ac1) ? robz(sn.ac1) : zeros(N)
        junk = good .& ((zc .> nsigma) .| (zac .> ac_sigma))
        jump = good .& (nz .> nsigma) .& .!junk .& [!worth(i) for i in 1:N]
        single = junk .| jump
        # single null pulse (an empty line in a steadily bright pulsar,
        # J1901): no signal at all (< zdet/2·σ; a real null is 0 ± 1σ) where
        # even the faintest 10% of pulses are clearly detected (>= 3·zdet·σ),
        # so it is not ordinary pulse-to-pulse modulation (lognormal σ=0.8:
        # ~0.02% false). Off for nulling/weak pulsars (there the stretch
        # test handles nulls).
        snr1 = [vone(i) > 0 ? s[i] / sqrt(vone(i)) : 0.0 for i in 1:N]
        rest = good .& .!single
        q10 = count(rest) > 0 ? quantile(snr1[rest], 0.1) : 0.0
        null1 = q10 >= 3 * zdet ? rest .& (snr1 .< zdet / 2) : falses(N)

        # 2./3. stretches, smoothed over the typical pulses
        ok = count(good .& .!single) > count(good) ÷ 2 ? (good .& .!single) : copy(good)
        if isnothing(w)
            S, V = sum(s[ok]), sum(v[ok])
            snr_pp = count(ok) > 0 && V > 0 ? S / sqrt(V) / sqrt(count(ok)) : 0.0
            w = snr_pp > 0 ? round(Int, (target / snr_pp)^2) : N ÷ 8
            w = clamp(isodd(w) ? w : w + 1, 3, max(3, N ÷ 8))
        end
        hw = w ÷ 2
        cs = vcat(0.0, cumsum(s .* ok)); cv = vcat(0.0, cumsum(v .* ok)); cn = vcat(0, cumsum(ok))
        lo_(i) = max(1, i - hw); hi_(i) = min(N, i + hw)
        Ss = [cs[hi_(i) + 1] - cs[lo_(i)] for i in 1:N]
        Vs = [cv[hi_(i) + 1] - cv[lo_(i)] for i in 1:N]
        Ns = [cn[hi_(i) + 1] - cn[lo_(i)] for i in 1:N]
        effw = [Vs[i] > 0 ? Ss[i] / Vs[i] : 0.0 for i in 1:N]
        e_on = quantile(effw[ok], q_on)
        # null = signal NOT DETECTED in the window (< zdet·σ) although typical
        # emission there would be clearly detected (e_on·V >= 2·zdet·σ) and
        # the window is below f_null of typical. A fainter-but-visible
        # stretch of a bright pulsar is not a null (that cut 32 pulses of a
        # clean J0304-like pulsar); real nulls and faded stretches are.
        # ...OR the window is clearly (`zrel`·σ) below `f_rel` (15%) of
        # typical: a faded stretch of a fairly bright pulsar is still weakly
        # "detected" (J1110-like: window S/N ~2 vs 22 when bright), but at
        # ~8% of normal it is junk; natural faint stretches stay well above 15%
        undet = [Ns[i] > 0 && Vs[i] > 0 && Ss[i] < zdet * sqrt(Vs[i]) &&
                 e_on * Vs[i] >= 2 * zdet * sqrt(Vs[i]) && Ss[i] < f_null * e_on * Vs[i] for i in 1:N]
        faded = [Ns[i] > 0 && Vs[i] > 0 && Ss[i] + zrel * sqrt(Vs[i]) < f_rel * e_on * Vs[i] for i in 1:N]
        # noisy: total noise OR white noise far above typical, over a SHORT
        # window of its own (noise needs no S/N; the signal window of a weak
        # pulsar can be 100+ pulses and would dilute a bad start). White
        # noise separately: with correlated baseline everywhere the total-
        # noise "typical" is inflated and a white-noise burst hides in it.
        hn = min_run
        lon(i) = max(1, i - hn); hin(i) = min(N, i + hn)
        cw = vcat(0.0, cumsum(sn.vwhite .* ok))
        Vn = [cv[hin(i) + 1] - cv[lon(i)] for i in 1:N]
        Wn = [cw[hin(i) + 1] - cw[lon(i)] for i in 1:N]
        Nn = [cn[hin(i) + 1] - cn[lon(i)] for i in 1:N]
        vtyp, wtyp = median(v[ok]), median(sn.vwhite[ok])
        # ...and brings proportionally more noise than signal (signal per
        # unit white noise < half of typical): a bright stretch whose noise
        # rose WITH its signal (digitiser scaling) is not junk
        # (each noise measure judged with its own signal/noise: stripes show
        # in the window noise V, a white-noise burst in the white noise W)
        Sn = [cs[hin(i) + 1] - cs[lon(i)] for i in 1:N]
        effv = [Vn[i] > 0 ? Sn[i] / Vn[i] : 0.0 for i in 1:N]
        effw_ = [Wn[i] > 0 ? Sn[i] / Wn[i] : 0.0 for i in 1:N]
        etv, etw = median(effv[ok]), median(effw_[ok])
        noisy = [Nn[i] > 0 &&
                 ((Vn[i] / (Nn[i] * vtyp) > noise_fac && effv[i] < 0.5 * etv) ||
                  (Wn[i] / (Nn[i] * wtyp) > noise_fac && effw_[i] < 0.5 * etw)) for i in 1:N]

        function runfilter(st, mr=min_run)
            st = copy(st); i = 1
            while i <= N
                if st[i]
                    j = i
                    while j < N && st[j+1]; j += 1; end
                    (j - i + 1) < mr && (st[i:j] .= false)
                    i = j + 1
                else
                    i += 1
                end
            end
            return st
        end
        function morph(mask)
            m = copy(mask); i = 1
            while i <= N
                if !m[i]
                    j = i
                    while j < N && !m[j+1]; j += 1; end
                    (i > 1 && j < N && (j - i + 1) < min_keep) && (m[i:j] .= true)
                    i = j + 1
                else
                    i += 1
                end
            end
            out = copy(m)
            for i in findall(m); out[max(1, i - pad):min(N, i + pad)] .= true; end
            return out
        end
        snr(mask) = (k = .!mask .& good; sv = sum(v[k]); sv > 0 ? sum(s[k]) / sqrt(sv) : NaN)

        # close short gaps inside a null region: with a long signal window,
        # neighbouring windows share most pulses, so a brief "detected" blip
        # in the middle of a faded stretch is a noise fluctuation (J1110-like)
        function closegaps(st, maxgap)
            st = copy(st)
            idx = findall(st)
            for k in 1:length(idx)-1
                g = idx[k+1] - idx[k] - 1
                0 < g <= maxgap && (st[idx[k]+1:idx[k+1]-1] .= true)
            end
            return st
        end
        # undetected (real nulls) from `min_run` pulses; merely FADED but still
        # detected stretches only when long (`min_run_fade`): in a bright,
        # strongly modulated pulsar (J1048-5832) short faint episodes of
        # real emission were cut as "nulls" all over (832 pulses) - that
        # raises the S/N but is cherry-picking; J1110-like fading is long
        nulls = runfilter(closegaps(undet .& ok, hw) .& good) .|
                runfilter(closegaps(faded .& ok, hw) .& good, max(min_run_fade, min_run))
        noisy = runfilter(noisy .& good)
        snr0 = snr(falses(N))
        # noisy stretches (already gated by signal/noise above) are cut
        # directly — the window-based S/N barely sees a white-noise burst
        # when baseline structure dominates the noise; nulls only if they
        # don't lower the S/N (safety)
        cur = .!good .| single .| noisy
        kept = any(noisy) ? ["zaszumione odcinki"] : String[]
        if any(nulls)
            if snr(morph(cur .| nulls)) >= snr(morph(cur))
                cur = cur .| nulls
                push!(kept, "nulle")
            else
                nulls .= false
            end
        end
        if any(null1)
            cur = cur .| null1
            push!(kept, "pojedyncze nulle")
        end
        stretch = nulls .| noisy
        best = morph(cur)
        why = isempty(kept) ? (any(single) ? "tylko pojedyncze śmieci" : "nic do wycięcia") :
              join(kept, " + ") * (any(single) ? " + pojedyncze śmieci" : "")

        reason = zeros(Int8, N)
        for i in findall(best)
            reason[i] = !good[i] ? 1 : single[i] ? 2 : (stretch[i] && noisy[i]) ? 4 :
                        ((stretch[i] && nulls[i]) || null1[i]) ? 3 : 5
        end
        return (bad=findall(best), single=findall(single), stretch=findall(stretch), reason=reason,
                eff=effw, E=f_null * e_on, w=w, snr_before=snr0, snr_after=snr(best), why=why,
                n_single=count(single), n_cand=count(stretch))
    end


    """
        auto_zap_core_scan(sn; ...)   — previous version, kept for comparison

    S/N-driven automatic zapping (built from hand-zapped examples: nulls,
    faded stretches, stripy/noisy stretches, bad starts — all of them are
    stretches that bring more noise than signal).

    Goal = the S/N of the summed pulses, i.e. what hand zapping optimises:

      1. single bad pulses are cut alone: a sudden jump of the white noise
         (robust z of log vwhite vs a running median > `nsigma`), or
         correlated junk (window variance / white variance above the χ²
         tail for `nwin` windows, ≈ 1 + nsigma·√(2/nwin)) — kept only if
         that doesn't lower the S/N;
      2. stretches: over a sliding window (singles excluded) each stretch
         has a signal-per-unit-noise efficiency; the cut level t is SCANNED
         (`n_thr` quantiles) and the one giving the highest S/N is kept.
         A stretch is cut only if even its signal + `zsig`·σ is below t
         (not just low by chance). Window length adapts to the pulsar: an
         average window has S/N ≈ `target` (bright pulsar -> short window,
         sharp null edges; weak -> long window), or fixed with `w`;
      3. stretches shorter than `min_run` pulses are not cut (ordinary
         pulse-to-pulse modulation); kept islands shorter than `min_keep`
         between cuts are cut too, and every cut region gets `pad` extra
         pulses on each side;
      4. stretch cutting is used only if it beats the singles-only result
         by at least `min_gain` (2%) in S/N.

    Returns (bad, single, stretch, eff, E = chosen cut level, w,
    snr_before, snr_after).
    """
    function auto_zap_core_scan(sn; target=8.0, w=nothing, zsig=3.0, pad=1, min_keep=3,
                           min_run=5, nsigma=6.0, min_gain=0.02, n_thr=60, local_win=31, noise_z=nothing)
        s, v, good = sn.s, sn.v, sn.good
        N = length(s)
        nz = isnothing(noise_z) ? noise_jump_z(sn.vwhite, good; local_win=local_win) : noise_z
        excess_thr = sn.nwin >= 3 ? 1 + nsigma * sqrt(2 / sn.nwin) : Inf
        single = .!good .| (nz .> nsigma) .| (sn.v ./ max.(sn.vwhite, eps()) .> excess_thr)

        # stretches are smoothed over the typical pulses (single suspects
        # left out, otherwise one stripe smears over a whole window); if the
        # suspects are most of the data, smooth over everything
        ok = count(good .& .!single) > count(good) ÷ 2 ? (good .& .!single) : copy(good)
        if isnothing(w)
            S, V = sum(s[ok]), sum(v[ok])
            snr_pp = count(ok) > 0 && V > 0 ? S / sqrt(V) / sqrt(count(ok)) : 0.0
            w = snr_pp > 0 ? round(Int, (target / snr_pp)^2) : N ÷ 8
            w = clamp(isodd(w) ? w : w + 1, 3, max(3, N ÷ 8))
        end
        hw = w ÷ 2
        # window sums via cumulative sums (O(N), so tuning over many
        # parameter sets stays fast)
        cs = vcat(0.0, cumsum(s .* ok))
        cv = vcat(0.0, cumsum(v .* ok))
        Ss = [cs[min(N, i + hw) + 1] - cs[max(1, i - hw)] for i in 1:N]
        Vs = [cv[min(N, i + hw) + 1] - cv[max(1, i - hw)] for i in 1:N]
        effw = [Vs[i] > 0 ? Ss[i] / Vs[i] : 0.0 for i in 1:N]

        # a few weak pulses in a row are ordinary pulse-to-pulse modulation,
        # not a null/bad stretch: only runs of >= min_run pulses are cut
        function runfilter(st)
            st = copy(st)
            i = 1
            while i <= N
                if st[i]
                    j = i
                    while j < N && st[j+1]; j += 1; end
                    (j - i + 1) < min_run && (st[i:j] .= false)
                    i = j + 1
                else
                    i += 1
                end
            end
            return st
        end

        function morph(mask)
            m = copy(mask)
            # short kept islands between two cut regions -> cut
            i = 1
            while i <= N
                if !m[i]
                    j = i
                    while j < N && !m[j+1]; j += 1; end
                    if i > 1 && j < N && (j - i + 1) < min_keep
                        m[i:j] .= true
                    end
                    i = j + 1
                else
                    i += 1
                end
            end
            out = copy(m)
            for i in findall(m)
                out[max(1, i - pad):min(N, i + pad)] .= true
            end
            return out
        end

        snr(mask) = (k = .!mask .& good; sv = sum(v[k]); sv > 0 ? sum(s[k]) / sqrt(sv) : NaN)
        snr0 = snr(falses(N))

        # reference: only single bad pulses (if that doesn't lower S/N),
        # otherwise only empty pulses
        ref, ref_single = morph(single), true
        snr(ref) >= snr0 || ((ref, ref_single) = (morph(.!good), false))
        ref_snr = snr(ref)

        # Stretches: scan the cut level t (signal per unit noise of a
        # stretch) and keep the one giving the HIGHEST S/N — the objective
        # behind hand zapping, computed. A fixed reference level failed on
        # J0034 (bright bursts + nulls): the "bright state" was so high that
        # nearly everything looked weak and nothing got cut. Guards against
        # chasing noise: a stretch must be below t even after +zsig·σ, at
        # least `min_run` long, and the gain over `ref` at least `min_gain`.
        # each cut level is tried with and without the single suspects, so a
        # bad batch of suspects can't veto good stretch cuts
        best, best_snr, t_best = ref, ref_snr, NaN
        stretch, with_single = falses(N), ref_single
        ncand = 0
        if count(ok) > 2
            for t in unique(quantile(effw[ok], range(0.02, 0.98, length=n_thr)))
                st = runfilter([ok[i] && Vs[i] > 0 && Ss[i] + zsig * sqrt(Vs[i]) < t * Vs[i] for i in 1:N])
                any(st) || continue
                ncand += 1
                for inc in (false, true)
                    m = morph(inc ? (single .| st) : (.!good .| st))
                    sm = snr(m)
                    if sm > best_snr
                        best, best_snr, t_best, stretch, with_single = m, sm, t, st, inc
                    end
                end
            end
        end
        if !(best_snr >= ref_snr * (1 + min_gain))
            best, t_best, stretch, with_single = ref, NaN, falses(N), ref_single
        end

        # why each pulse is cut (for the review plot):
        # 0 kept, 1 empty, 2 single junk, 3 weak/null stretch, 4 noisy stretch, 5 margin/fill
        vmed = median(v[good])
        hwv = [Vs[i] / max(1, count(ok[max(1, i - hw):min(N, i + hw)])) for i in 1:N]
        reason = zeros(Int8, N)
        for i in findall(best)
            reason[i] = !good[i] ? 1 :
                        (with_single && single[i]) ? 2 :
                        stretch[i] ? (hwv[i] > 2 * vmed ? 4 : 3) : 5
        end
        why = !isnan(t_best) ? "odcinki (próg $(round(t_best, sigdigits=3)))" * (with_single ? " + pojedyncze" : "") :
              ref_single && any(single) ? "tylko pojedyncze śmieci" :
              ncand == 0 ? "żaden odcinek nie jest wyraźnie słabszy" :
              "wycięcie odcinków nie podnosi S/N o ≥$(round(100min_gain))%"
        return (bad=findall(best), single=findall(single), stretch=findall(stretch), reason=reason,
                eff=effw, E=t_best, w=w, snr_before=snr0, snr_after=snr(best), why=why,
                n_single=count(single .& good), n_cand=ncand)
    end


    """
        profile_gate(data, bin_st, bin_end, reason) -> rows to restore

    Veto on the automat's cuts with the standard profile S/N (`profile_snr`):
    every contiguous cut run (reasons 2-5; empty rows don't matter) is put
    back if the integrated profile has a higher S/N with it, repeated until
    nothing changes. If the cuts still lower the S/N, all are put back.
    Only removes cuts, so it can't cherry-pick faint emission. Needed
    because on 199 pulsars the cuts alone lowered the profile S/N in more
    than half of them (median -0.4%, 12% of pulses cut vs Maciej's 1%):
    pulses with some noise but real signal are worth keeping.
    """
    function profile_gate(data, bin_st, bin_end, reason)
        N, nb = size(data)
        on, off = bin_st:bin_end, offpulse_bins(data, bin_st, bin_end)
        psnr(P) = (μ = mean(P[off]); σ = std(P[off]); σ > 0 ? sum(P[on] .- μ) / (σ * sqrt(length(on))) : -Inf)
        runs = UnitRange{Int}[]
        i = 1
        while i <= N
            if reason[i] >= 2
                j = i
                while j < N && reason[j+1] >= 2; j += 1; end
                push!(runs, i:j); i = j + 1
            else
                i += 1
            end
        end
        isempty(runs) && return Int[]
        Ptot = vec(sum(data, dims=1))
        Pr = [vec(sum(@view(data[r, :]), dims=1)) for r in runs]
        P = Ptot .- sum(Pr)
        cut = trues(length(runs))
        changed = true
        while changed
            changed = false
            s = psnr(P)
            for k in findall(cut)
                s2 = psnr(P .+ Pr[k])
                if s2 > s
                    P .+= Pr[k]; cut[k] = false; s = s2; changed = true
                end
            end
        end
        psnr(P) < psnr(Ptot) && (cut .= false)
        return reduce(vcat, [collect(runs[k]) for k in findall(.!cut)]; init=Int[])
    end


    """
    Convenience: `auto_zap_core` on a data matrix, then (`gate=true`) the
    profile-S/N veto `profile_gate` — vetoed pulses get reason 0.
    """
    function auto_zap(data, bin_st, bin_end; gate=true, kw...)
        sn = pulse_signal_noise(data, bin_st, bin_end)
        a = auto_zap_core(sn; kw...)
        gate || return a
        back = profile_gate(data, bin_st, bin_end, a.reason)
        isempty(back) && return a
        reason = copy(a.reason); reason[back] .= 0
        bad = findall(>(0), reason)
        return merge(a, (bad=bad, reason=reason, snr_after=zap_snr(sn, bad),
                         why=a.why * " (bramka S/N profilu przywróciła $(length(back)) imp.)"))
    end


    "S/N of the summed pulses (`pulse_signal_noise` result) with `rows` (1-based) removed."
    function zap_snr(sn, rows)
        k = copy(sn.good)
        k[filter(i -> 1 <= i <= length(k), rows)] .= false
        sv = sum(sn.v[k])
        return sv > 0 ? sum(sn.s[k]) / sqrt(sv) : NaN
    end


    "Pulse-level agreement of automatic vs hand zaps (1-based rows)."
    function zap_agreement(auto_rows, user_rows, N)
        a = falses(N); a[auto_rows] .= true
        u = falses(N); u[filter(i -> 1 <= i <= N, user_rows)] .= true
        tp, fp, fn = count(a .& u), count(a .& .!u), count(.!a .& u)
        prec = tp + fp > 0 ? tp / (tp + fp) : (count(u) == 0 ? 1.0 : 0.0)
        rec = tp + fn > 0 ? tp / (tp + fn) : 1.0
        f1 = prec + rec > 0 ? 2prec * rec / (prec + rec) : 0.0
        # both empty = perfect agreement
        (count(a) == 0 && count(u) == 0) && (f1 = 1.0)
        return (tp=tp, fp=fp, fn=fn, precision=prec, recall=rec, f1=f1)
    end


    "Hand labels (1-based rows) from input/zap_labels/<psr>.csv (`label` column), or nothing."
    function read_zap_labels(psr)
        f = joinpath(@__DIR__, "..", "input", "zap_labels", "$psr.csv")
        isfile(f) || return nothing
        lines = readlines(f)
        hdr = split(lines[1], ",")
        li = findfirst(==("label"), hdr)
        isnothing(li) && return nothing
        return [k - 1 for k in 2:length(lines) if split(lines[k], ",")[li] == "1"]
    end


    "Thrown by `review_zaps` when the user chooses to skip the pulsar (`s`)."
    struct SkipPulsar <: Exception end


    # 0-based "12 40-45" -> 1-based rows, clamped to 1:N
    zaps_to_rows(zaps, N) = sort(unique(filter(i -> 1 <= i <= N,
                                 [i + 1 for (a, b) in zaps_to_ranges(zaps) for i in a:b])))


    """
        select_pulses(data, r, rows, bin_st, bin_end, outdir; allow_skip=true, name_mod="pulsar")

    Interactive selection loop on the `Plot.zap_review` window (left before,
    right after). Every change redraws the window; returns the final 1-based
    rows. Commands (pulse numbers 0-based, as in zaps):

      Enter        accept          n   nothing / accept empty
      s            skip pulsar (throws `SkipPulsar`, only if `allow_skip`)
      +12 40-45    add             -40-45   remove
      12 40-45     replace all     u        undo last change
    """
    function select_pulses(data, r, rows, bin_st, bin_end, outdir; allow_skip=true, name_mod="pulsar")
        N = size(data, 1)
        rows = copy(rows)
        history = Vector{Vector{Int}}()
        while true
            Plot.zap_review(data, r, rows, bin_st, bin_end, outdir; show_=true, name_mod=name_mod)
            print("\nEnter = zatwierdź | n = nic | " * (allow_skip ? "s = pomiń pulsar | " : "") *
                  "+12 40-45 = dodaj | -40-45 = usuń | 12 40-45 = zastąp | u = cofnij:  ")
            ans = String(strip(readline()))
            Plot.close_all()
            isempty(ans) && break
            allow_skip && ans == "s" && throw(SkipPulsar())
            if ans == "n"
                rows = Int[]
                break
            end
            if ans == "u"
                isempty(history) ? println("Nie ma czego cofać.") : (rows = pop!(history))
            else
                try
                    new = if startswith(ans, "+")
                        sort(unique(vcat(rows, zaps_to_rows(ans[2:end], N))))
                    elseif startswith(ans, "-")
                        setdiff(rows, zaps_to_rows(ans[2:end], N))
                    else
                        zaps_to_rows(ans, N)
                    end
                    push!(history, rows)
                    rows = new
                catch
                    println("Nie rozumiem \"$ans\" — przykład: +12 40-45")
                    continue
                end
            end
            println("Teraz: $(length(rows)) impulsów  →  \"$(pulses_to_zaps(rows))\"")
        end
        return rows
    end

    """
        review_zaps(txtfile, p, params_file, outdir; mode=:auto_pazi, archive, ref_archive)

    Chooses the bad pulses of MY analysis and saves them in `p["zaps2"]`
    (+ `p["zap_chans2"]`). `p["zaps"]` (Maciej's / older runs) is never
    changed — only created empty if missing.

    Modes: `:auto_pazi` (default use) automat proposal + picture, then
    Enter / p (fix in pazi) / r (redo by hand) / n (none) / s (skip);
    `:auto` automat without asking; `:pazi` by hand in pazi; `:review` /
    old `find_bad_pulses` window.

    Comparison (`p["zap_cmp"]`): on the SAME unzapped data and on-pulse
    window, the standard profile S/N (`profile_snr`) and the automat's own
    measure for: no zaps / Maciej's zaps (pulses zapped in `ref_archive`,
    e.g. pazur/<psr>/pulsar.spCF, plus his `p["zaps"]`) / my zaps.
    """
    function review_zaps(txtfile, p, params_file, outdir; nsigma=5.0, mode=:review, archive=nothing,
                         ref_archive=nothing)
        data = load_ascii(txtfile)
        N = size(data, 1)
        bst, bend = p["bin_st"], p["bin_end"]
        haskey(p, "zaps") || (p["zaps"] = nothing)
        snr(rows) = profile_snr(data, bst, bend; drop=rows)
        r1(x) = round(x, digits=1)
        chans = Int[]

        if mode in (:auto_pazi, :auto)
            a = auto_zap(data, bst, bend)
            println("Automat: $(length(a.bad)) z $N impulsów ($(r1(100 * length(a.bad) / N))%)  →  " *
                    "\"$(pulses_to_zaps(a.bad))\"")
            kinds = ["pusty", "pojedynczy śmieć", "null/brak sygnału", "zaszumiony odcinek", "margines"]
            println("Powody: " * join(["$(kinds[k])=$(count(==(k), a.reason))" for k in 1:5], "  "))
            s0, s1 = snr(Int[]), snr(a.bad)
            println("S/N profilu: bez zapów=$(r1(s0))  po automacie=$(r1(s1))   " *
                    "(miara automatu: $(r1(a.snr_before)) → $(r1(a.snr_after)))")
            rows = copy(a.bad)
            if mode == :auto
                Plot.autozap_view(data, bst, bend, a, outdir; show_=false,
                                  title_="$(basename(rstrip(outdir, '/'))) — automat", snr_pazi=(s0, s1))
            else
                isnothing(archive) && error("review_zaps(mode=:auto_pazi) needs `archive`")
                viewer = Plot.autozap_view(data, bst, bend, a, outdir; show_=true, wait_=false,
                                           title_="$(basename(rstrip(outdir, '/'))) — automat", snr_pazi=(s0, s1))
                while true
                    print("Enter = zatwierdź | p = popraw w pazi (dane PO automacie) | r = od nowa ręcznie " *
                          "w pazi | n = bez zapów | s = pomiń pulsar: ")
                    ans = strip(readline())
                    if ans == "s"
                        Plot.close_png(viewer)
                        throw(SkipPulsar())
                    elseif ans == "n"
                        rows, chans = Int[], Int[]
                    elseif ans == "r"
                        z = run_pazi(archive)
                        rows, chans = [s + 1 for s in z.subints if 0 <= s < N], z.chans
                    elseif ans == "p"
                        preview = replace(archive, r"\.[^.]+$" => ".autozap.ar")
                        cp(archive, preview; force=true)
                        apply_zaps(preview, pulses_to_zaps(a.bad))
                        println("pazi pokazuje dane PO automacie — dozapuj, co jeszcze złe, i wyjdź.")
                        z = run_pazi(preview)
                        rm(preview, force=true)
                        rows = sort(unique(vcat(a.bad, [s + 1 for s in z.subints if 0 <= s < N])))
                        chans = z.chans
                    elseif ans != ""
                        println("Nie rozumiem \"$ans\".")
                        continue
                    end
                    break
                end
                Plot.close_png(viewer)
            end
        elseif mode == :pazi
            isnothing(archive) && error("review_zaps(mode=:pazi) needs `archive`")
            z = run_pazi(archive)
            rows, chans = [s + 1 for s in z.subints if 0 <= s < N], z.chans
        else
            r = find_bad_pulses(data, bst, bend; nsigma=nsigma)
            rows = copy(r.bad)
            println("Automatycznie (stary detektor): $(length(rows)) złych impulsów z $N  →  " *
                    "\"$(pulses_to_zaps(rows))\"")
            if mode == :review
                rows = select_pulses(data, r, rows, bst, bend, outdir)
            else
                Plot.zap_review(data, r, rows, bst, bend, outdir; show_=false)
            end
        end

        # final choice = training data for the detector too
        mode in (:auto_pazi, :pazi) && save_zap_labels(basename(rstrip(outdir, '/')), data, bst, bend, rows)
        zaps = pulses_to_zaps(rows)
        p["zaps2"] = isempty(zaps) ? nothing : zaps
        p["zap_chans2"] = isempty(chans) ? nothing : join(chans, " ")

        # comparison with Maciej's zapping (same data, same window)
        sn = pulse_signal_noise(data, bst, bend)
        cmp = Dict{String,Any}("n_pulses" => N, "n_mine" => length(rows),
                               "snr_none" => snr(Int[]), "snr_mine" => snr(rows),
                               "auto_none" => zap_snr(sn, Int[]), "auto_mine" => zap_snr(sn, rows))
        ref = isnothing(ref_archive) ? nothing : zapped_rows(ref_archive, bst, bend)
        mrows = nothing
        if !isnothing(ref) && ref.N == N
            mrows = sort(unique(vcat(ref.rows, zaps_to_rows(p["zaps"], N))))
            cmp["n_maciej"] = length(mrows)
            cmp["snr_maciej"] = snr(mrows)
            cmp["auto_maciej"] = zap_snr(sn, mrows)
            cmp["zaps_maciej"] = pulses_to_zaps(mrows)
        elseif !isnothing(ref)
            println("UWAGA: archiwum Macieja ma $(ref.N) impulsów, moje $N — nie porównuję.")
        end
        # pazi-like contrast before / automat / Maciej (pulsar_contrast.pdf)
        mode == :auto && Plot.zap_contrast(data, bst, bend, rows, haskey(cmp, "n_maciej") ? mrows : nothing,
                                           outdir; title_=basename(rstrip(outdir, '/')),
                                           off=offpulse_bins(data, bst, bend))
        p["zap_cmp"] = Dict(k => (v isa AbstractFloat && !isfinite(v)) ? nothing : v for (k, v) in cmp)
        Tools.save_params(params_file, p)
        println("Zapisane zaps2 = $(repr(p["zaps2"]))  zap_chans2 = $(repr(p["zap_chans2"]))   " *
                "(zaps Macieja/stare bez zmian: $(repr(p["zaps"])))")
        println("S/N profilu: bez zapów=$(r1(cmp["snr_none"]))  moje=$(r1(cmp["snr_mine"]))" *
                (haskey(cmp, "snr_maciej") ? "  Macieja=$(r1(cmp["snr_maciej"])) ($(cmp["n_maciej"]) imp.)" :
                                             "  Macieja: brak"))
        return p["zaps2"]
    end


    """
    Zero the weights of the `zaps` subints in a PSRCHIVE archive (in place),
    ranges as `-W "a b"`, single pulses as `-w "x y"`.
    """
    function apply_zaps(archive, zaps; chans=nothing)
        rngs = zaps_to_ranges(zaps)
        args = String[]
        # channels zapped in pazi (`p["zap_chans"]`, space-separated list)
        if !isnothing(chans) && !isempty(strip(chans))
            push!(args, "-z", chans)
            println("ZAPPING channels: $chans")
        end
        isempty(rngs) && isempty(args) && return
        singles = [string(a) for (a, b) in rngs if a == b]
        for (a, b) in rngs
            a == b && continue
            push!(args, "-W", "$a $b")
        end
        isempty(singles) || push!(args, "-w", join(singles, " "))
        println("ZAPPING subints: $zaps")
        run(pipeline(`paz -m $args $archive`, stderr="errs.txt"))
    end


    """
    S/N exactly as pazi / PSRCHIVE reports it ("Total S/N": `psrstat -c snr`
    on the fully scrunched profile, -j FTp), optionally after zapping
    `zaps` on a temporary copy. NaN if psrstat fails.
    """
    function psrchive_snr(archive; zaps=nothing)
        target = archive
        if !isnothing(zaps) && !isempty(zaps_to_ranges(zaps))
            target = replace(archive, r"\.[^.]+$" => ".snrtmp.ar")
            cp(archive, target; force=true)
            redirect_stdout(devnull) do
                apply_zaps(target, zaps)
            end
        end
        out = try
            read(pipeline(`psrstat -j FTp -c snr $target`; stderr=devnull), String)
        catch
            ""
        end
        target == archive || rm(target, force=true)
        m = match(r"snr=([-+0-9.eE]+)", out)
        return isnothing(m) ? NaN : parse(Float64, m.captures[1])
    end


    """
    Pulse zaps used by MY analysis: `p["zaps2"]` when it exists (also when
    it is `nothing` = checked, nothing cut), otherwise the old `p["zaps"]`
    (Maciej's / earlier runs). `p["zaps"]` itself is never overwritten.
    """
    zaps_of(p) = haskey(p, "zaps2") ? p["zaps2"] : get(p, "zaps", nothing)
    zap_chans_of(p) = haskey(p, "zaps2") ? get(p, "zap_chans2", nothing) : get(p, "zap_chans", nothing)


    """
        profile_snr(data, bin_st, bin_end; drop=Int[]) -> Float64

    Standard S/N of the integrated profile (Lorimer & Kramer): sum the
    pulses (all-zero rows = already zapped, and rows `drop` left out), then

        S/N = Σ_on (P − mean_off) / (σ_off · √n_on)

    with the fixed on-pulse window `bin_st:bin_end` and σ_off the rms of the
    integrated profile off the pulse. Same data + same window for every zap
    set, so the comparison between methods is fair (unlike pazi, which picks
    its own width).
    """
    function profile_snr(data, bin_st, bin_end; drop=Int[])
        N, nb = size(data)
        keep = [any(!iszero, @view data[i, :]) for i in 1:N]
        keep[filter(i -> 1 <= i <= N, collect(drop))] .= false
        count(keep) == 0 && return NaN
        P = vec(sum(@view(data[keep, :]), dims=1))
        on = bin_st:bin_end
        off = offpulse_bins(data, bin_st, bin_end)
        μ, σ = mean(P[off]), std(P[off])
        return σ > 0 ? sum(P[on] .- μ) / (σ * sqrt(length(on))) : NaN
    end


    """
    Pulses zapped in an EXISTING archive (e.g. Maciej's pazur/<psr>/pulsar.spCF,
    zapped directly in pazi): debase a temporary copy and return the 1-based
    rows that came out all-zero, plus the number of pulses. `nothing` if the
    archive can't be read.
    """
    function zapped_rows(archive, bin_st, bin_end)
        isfile(archive) || return nothing
        dir = mktempdir()
        try
            tmp = joinpath(dir, "ref.spCF")
            cp(archive, tmp)
            run(pipeline(`pmod -onpulse "$bin_st $bin_end" -device "/NULL" -debase $tmp`; stdout=devnull, stderr=devnull))
            txt = joinpath(dir, "ref.debase.txt")
            convert_psrfit_ascii(joinpath(dir, "ref.debase.gg"), txt)
            d = load_ascii(txt)
            return (rows=[i for i in axes(d, 1) if !any(!iszero, @view d[i, :])], N=size(d, 1))
        catch e
            println("Nie udało się odczytać zapów z $archive: $(sprint(showerror, e))")
            return nothing
        finally
            rm(dir; recursive=true, force=true)
        end
    end


    """
    Loads PSRCHIVE ASCII file
    """
    function load_ascii(infile)
        f = open(infile)
        lines = readlines(f)
        res = split(lines[1])
        pulses = parse(Int, res[6])
        bins = parse(Int, res[12])
        data = Array{Float64}(undef, pulses, bins)
        for i in 2: length(lines)
            res = split(lines[i])
            bin = parse(Int, res[3]) + 1
            pulse = parse(Int, res[1]) + 1
            inte = parse(Float64, res[4])
            data[pulse, bin] = inte
        end
        close(f)
        return data
    end


    """
        Loads PSRCHIVE ASCII file with all 4 polarizations.
        Returns a 3D array: (pulse, bin, pol)
    """
    function load_ascii_all(infile)
        f = open(infile)
        lines = readlines(f)
        
        header = split(lines[1])
        pulses = parse(Int, header[6])    # Nsub
        bins = parse(Int, header[12])     # Nbin
        pols = parse(Int, header[10])     # Npol

        data = Array{Float64}(undef, pulses, bins, pols)

        for i in eachindex(lines)
            startswith(lines[i], "F") && continue
            startswith(lines[i], "M") && continue
            res = split(lines[i])
            pulse = parse(Int, res[1]) + 1
            bin = parse(Int, res[3]) + 1
            for pol in 1:pols
                val = parse(Float64, res[3 + pol])
                data[pulse, bin, pol] = val
            end
        end

        close(f)
        return data
    end


    """
    Save PSRCHIVE ASCII file
    """
    function save_ascii(data, outfile)
        f = open(outfile, "w")
        pulse_num = size(data)[1]
        bin_num = size(data)[2]
        header = "File: filename Src: J1750-3503 Nsub: $pulse_num Nch: 1 Npol: 1 Nbin: $bin_num RMS: 0.0\n"
        write(f, header)

        for i in 1:pulse_num
            for j in 1:bin_num
                write(f, "$(i-1) 0 $(j-1) $(data[i,j])\n")
            end
        end
        close(f)
        return
    end


    """
    Loads data in FITS format
    """
    function load_fits(filename)
        f = FITS(filename)
        for hdu in f; println(typeof(hdu)); end
        he = read_header(f[5])
        println(he)
        close(f)
    end

    """
    Converts PSRFIT file to ASCII (using PSRCHIVE tools)
    """
    function convert_psrfit_ascii(infile, outfile)
        run(pipeline(`pdv -t -F $infile`, stdout="$outfile", stderr="errs.txt"))
        # change -t to -A to get frequancy information
        #@showprogress 1 for i in 1:pn  # psrchive indexing
        #end
    end


    """
    Uses PSRCHIVE to add .spCF files 
    
    indir: input directory
    outfile: output file name (full path)
    files: list of files to add (filenames only)
    txt: convert to txt?
    """
    function add_psrfiles(indir, outfile; files=nothing, sixteen=false)

        if files === nothing
            # Find all .spCF files in the input directory
            if sixteen == false
                files = filter(f -> endswith(f, ".spCF"), readdir(indir))
            else
                files = filter(f -> endswith(f, ".spCf16"), readdir(indir))
            end
        end
        
        # Sort files based on pulse numbers (e.g., 00000-00255, 00256-00511, etc.)
        sort!(files, by = f -> begin
            # Extract the pulse range from the filename (e.g., "2019-12-15-03:19:04_00000-00255.spCF")
            m = match(r"_(\d+)-(\d+)\.spC(?:F|f16)$", f)
            if isnothing(m)
                return typemax(Int)  # Files without proper format go to the end
            else
                return parse(Int, m.captures[1])  # Sort by the starting pulse number
            end
        end)
    
        if isempty(files)
            error("No .spCF or .spCf16 files found in directory: $indir")
        end
    
        println("Processing files in order:")
        for (i, f) in enumerate(files)
            println("$i. $f")
        end
    
        file_names = [joinpath(indir, file) for file in files]

        # connecting all files
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt")) # PSRCHIVE

    end


    function get_nsubint(outfile, params_file, params)
        # gets number of single pulses
        outstr = read(pipeline(`psrstat -c nsubint $outfile`, stderr="errs.txt"), String)
        nsubint = parse(Int, split(outstr, "=")[2])
        params["nsubint"] = nsubint
        println("Number of subintervals: $nsubint")
        Tools.save_params(params_file, params)
    end
 
    
    """
    Debase the data

    # Arguments
    """
    function debase(infile, params_file, params)


        if isnothing(params["bin_st"]) || isnothing(params["bin_end"])

            println("pmod -device \"/xw\" -debase $infile")
            run(pipeline(`pmod -device "/xw" -debase $infile`, `tee pmod_output.txt`))
            #run(pipeline(`pmod -device "/xw" -iformat PSRFITS -debase $outfile`, `tee pmod_output.txt`))

            # Read captured output to get bin_st and bin_end
            output = read("pmod_output.txt", String)
            rm("pmod_output.txt")  # cleanup

            # Extract onpulse values
            m = match(r"-onpulse '(.+) (.+)'", output)
            if !isnothing(m)
                bin_st, bin_end = parse.(Int, m.captures)
                # Check if onpulse region length is even
                region_length = bin_end - bin_st + 1
                if region_length % 2 != 0
                    println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                    bin_end -= 1
                    println("Adjusted onpulse range: $bin_st to $bin_end")
                end
                println("Found onpulse range: $bin_st to $bin_end")
                params["bin_st"] = bin_st
                params["bin_end"] = bin_end
                Tools.save_params(params_file, params)
            end
            # ASCII single pulses
            convert_psrfit_ascii(replace(infile, ".spCF"=>".debase.gg"), replace(infile, ".spCF"=>".debase.txt"))

        else
            run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $infile`))
            #run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -iformat PSRFITS -debase $outfile`))
            # ASCII single pulses
            convert_psrfit_ascii(replace(infile, ".spCF"=>".debase.gg"), replace(infile, ".spCF"=>".debase.txt"))
        end
    end


    """
    Calculate 2dfs and lrfs using PSRSALSA

    # ADD description!
    debased_file ??
    """
    function twodfs_lrfs(debased_file, params_file, p; detect=false)

        # Calculate 2dfs and lrfs
        run(pipeline(`pspec -w -oformat ASCII -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft $(p["nfft"]) -onpulse "$(p["bin_st"]) $(p["bin_end"])" $debased_file`,  stderr="errs.txt"))

        if p["p3"] == -1.0 || detect == true
            # Find P3
            run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))

            # Read captured output
            output = read("pspecDetect_output.txt", String)
            rm("pspecDetect_output.txt")  # cleanup

            # Extract P3 value from the last occurrence
            p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
            if !isempty(p3_matches)
                last_match = p3_matches[end]
                p3_value = parse(Float64, last_match.captures[1])
                p3_error = parse(Float64, last_match.captures[2])
                println("Found P3 = $p3_value ± $p3_error P0")
            else
                # no P3 found: type it by hand (e.g. read off the 2DFS/LRFS,
                # or from literature), or leave the pulsar marked "no P3"
                print("Nie znaleziono P3. Wpisz P3 ręcznie [P0] albo Enter = brak P3: ")
                ans = strip(readline())
                v = isempty(ans) ? nothing : tryparse(Float64, ans)
                if isnothing(v) || v <= 0
                    println("Brak P3 — pulsar oznaczony jako bez P3 (p3 = -1).")
                    p["p3"] = -1.0
                    Tools.save_params(params_file, p)
                    return
                end
                p3_value, p3_error = v, -1.0
            end

            ybins = Functions.find_ybins(p3_value)
            println("Number of ybins: $ybins")
            p["p3"] = p3_value
            p["p3_error"] = p3_error
            p["p3_ybins"] = ybins
            Tools.save_params(params_file, p)
        end
    end


    """
    Process data with my old routines
    """
    function process_data_andrzej(debased_file, outdir, p)

        println("DE $debased_file")
        # all data
        da0 = Data.load_ascii(replace(debased_file, ".gg"=>".txt"))  
        Plot.single(da0, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="all", show_=true)

        folded = Tools.p3fold(da0, p["p3"],  p["p3_ybins"])
        Plot.single(folded, outdir; darkness=0.97, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold", show_=true, repeat_num=4)

        println("Polarization cleaning threshold: ", p["clean_threshold"])
        # polarisation cleaning
        da = Data.load_ascii_all(replace(debased_file, ".gg"=>".txt"))  
        da2 = clean(da; threshold=p["clean_threshold"]) 
        Plot.single(da2, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="cleaned", show_=true)

        folded = Tools.p3fold(da2, p["p3"],  p["p3_ybins"])
        Plot.single(folded, outdir; darkness=0.97, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold_clean", show_=true, repeat_num=4)

        
    end



    """
    Process data with PSRCHIVE and PSRSALSA

    Bad pulses: after the first debase, `review_zaps` (`zap_mode=:pazi`,
    default: zap by hand in pazi, choice also saved as detector training
    data; `:review` my window; `:auto` automatic; `:off` keep
    `p["zaps"]`/`p["zap_chans"]`). Pulse and channel zaps are then applied
    with `paz` and the data debased again, so every later step (2DFS, P3,
    pfold, pulsar.debase.txt) works on the cleaned data.
    """
    function process_psrdata(indir, outdir; files=nothing, outfile="pulsar.spCF", params_file="params.json",
                             zap_mode=:pazi, zap_nsigma=5.0, ref_archive=nothing)

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        # ensure directory exists
        params_dir = dirname(params_file)
        if !isdir(params_dir)
            println("Directory $params_dir does not exist, creating it.")
            mkpath(params_dir)
        end

        # check if params_file exists if not creating default one
        if !isfile(params_file)
            println("File $params_file does not exist, creating default one.")
            p = Tools.default_params(params_file)
        else
            p = Tools.read_params(params_file)
        end

        # add all .spCF files 
        add_psrfiles(indir, outfile; files=files)

        # gets number of single pulses if needed
        if isnothing(p["nsubint"])
            get_nsubint(outfile, params_file, p)
        end

        # debase the data
        debase(outfile, params_file, p)

        # bad pulses: detect on the debased data, then zap + debase again
        # (outfile is rebuilt by psradd on every run, so zaps are re-applied)
        # (my zaps live in p["zaps2"]; p["zaps"] = Maciej's / older, untouched)
        zaps = zap_mode == :off ? zaps_of(p) :
               review_zaps(replace(outfile, ".spCF" => ".debase.txt"), p, params_file, outdir;
                           nsigma=zap_nsigma, mode=zap_mode, archive=outfile, ref_archive=ref_archive)
        chans = zap_chans_of(p)
        if !isnothing(zaps) || !isnothing(chans)
            apply_zaps(outfile, zaps; chans=chans)
            debase(outfile, params_file, p)
        end

        # Calculate 2dfs and lrfs
        twodfs_lrfs(debased_file, params_file, p; detect=false)

        # calculate p3-folded profile (impossible without P3)
        if !(p["p3"] isa Number) || p["p3"] <= 0
            println("Brak P3 — pomijam P3-fold.")
            return p, debased_file, outdir
        end
        println("P3-folding with:")
        # TODO experiment here
        println("pfold -p3fold_noonpulse -p3fold \"$(p["p3"]) $(p["p3_ybins"])\" -onpulse \"$(p["bin_st"]) $(p["bin_end"])\" -onpulsed \"/NULL\" -p3foldd \"/NULL\" -w -oformat ascii $debased_file")
        run(pipeline(`pfold  -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`,  stderr="errs.txt"))

        return p, debased_file, outdir
    end


    function plot_psrdata(data_dir, params)

        p = params

        d = load_ascii(joinpath(data_dir, "pulsar.debase.txt"))

        # Single pulses
        Plot.single(d, data_dir; darkness=0.5, bin_st= p["bin_st"], bin_end=p["bin_end"], start=p["pulse_start"], number= (p["pulse_end"] - p["pulse_start"]), name_mod="pulsar", show_=true)

        # 2DFS
        d2 = load_ascii(joinpath(data_dir, "pulsar.debase.1.2dfs"))
        Plot.twodfs(d2, data_dir, p; name_mod="pulsar", show_=true, average=Tools.average_profile(d))

        # LRFS
        d3 = load_ascii_all(joinpath(data_dir, "pulsar.debase.lrfs"))
        Plot.lrfs(d3, data_dir, p; name_mod="pulsar", show_=true)

        # P3-folded
        d5 = Data.load_ascii(joinpath(data_dir, "pulsar.debase.p3fold"))
        Plot.p3fold(d5, data_dir; start=1, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="pulsar", show_=true, repeat_num=4)
 
    end


    function process_all_data(outdir; base_root="/home/psr/data/new")
        # Get all subdirectories in base_root
        dirs = filter(isdir, readdir(base_root, join=true))

        if isempty(dirs)
            println("No data found in $base_root. Exiting...")
            return
        end
        
        for dir in dirs
            name = basename(dir)  # Extract directory name -> pulsar name
            data_dir = joinpath(dir, readdir(dir)[1]) # one level in
            out_dir = joinpath(outdir, name)

            println("\n\n"*"#"^79*"\n"*"#"^79*"\n"*"#"^79)
            println("Processing pulsar: ", name)
            println(out_dir)
            
            # TODO temporatory fix
            if isfile(joinpath(out_dir,"pulsar_single.pdf"))
                println("SKIPPING pulsar $name !!!!\n\n")
                continue 
            end

            # Creates out_dir if does not exists...
            if !isdir(out_dir)
                mkpath(out_dir)
                println("Created output directory: ", out_dir)
            end

            # TODO temporatory fix
            try
                # process data
                p = process_psrdata(data_dir, out_dir)
            
                # plot data
                plot_psrdata(out_dir, p)
            catch
                println("ERROR in pulsar $name !!!!\n\n")
            end
        end
    end

    function combine_pngs(outdir)
        paths = filter(isdir, readdir(outdir, join=true))
        pulsars = basename.(paths)

        pdfs = []

        # create pdf for all pulsars
        for (i,pulsar) in enumerate(pulsars)

            fig = Figure(size = (500, 800))
            Label(fig[0, :], text = pulsar, fontsize = 20, tellwidth = false)

            ax1 = Axis(fig[1, 1])
            ax2 = Axis(fig[1, 2])
            ax3 = Axis(fig[2, 1])
            ax4 = Axis(fig[2, 2])
            ax5 = Axis(fig[3, 1])
            ax6 = Axis(fig[3, 2])


            # load images
            try
                img1 = load(joinpath(paths[i], "pulsar_single.png"))
                image!(ax1, rotr90(img1))
            catch e
                println("No pulsar_single.png for $pulsar")
            end
            try
                img2 = load(joinpath(paths[i], "pulsar_lrfs.png"))
                image!(ax2, rotr90(img2))
            catch e
                println("No pulsar_lrfs.png for $pulsar")
            end
            try
                img3 = load(joinpath(paths[i], "pulsar_2dfs.png"))
                image!(ax3, rotr90(img3))
            catch e
                println("No pulsar_2dfs.png for $pulsar")
            end
            try
                img4 = load(joinpath(paths[i], "pulsar_p3fold.png"))
                image!(ax4, rotr90(img4))
            catch e
                println("No pulsar_p3fold.png for $pulsar")
            end
            try
                img5 = load(joinpath(paths[i], "all_single.png"))
                image!(ax5, rotr90(img5))
            catch e
                println("No  all_single.png for $pulsar")
            end
            try
                img6 = load(joinpath(paths[i], "p3fold_single.png"))
                image!(ax6, rotr90(img6))
            catch e
                println("No p3fold_single.png for $pulsar")
            end


            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
            hidedecorations!(ax4)
            hidedecorations!(ax5)
            hidedecorations!(ax6)


            filename = joinpath(outdir, "page_$i.pdf")
            save(filename, fig)
            push!(pdfs, basename(filename))
        end

        # combine pdfs
        cmd = vcat(["pdfunite"], pdfs, ["all.pdf"])
        write("convert.txt", join(cmd, " "))
        println(join(cmd, " "))

    end


    function combine_pngs_to_pdf(out_dir::String)

        # Collect and sort PNG files recursively
        png_paths = sort([
            joinpath(root, file)
            for (root, _, files) in walkdir(out_dir)
            for file in files
            if endswith(lowercase(file), ".png")
        ])

        isempty(png_paths) && return

        pulsar_name = basename(out_dir)
        output_pdf_dir = "/home/psr/output/"
        mkpath(output_pdf_dir)

        pages = Figure[]
        pdf_paths = String[]
        images_per_page = 4
        fig=nothing

        for (i, path) in enumerate(png_paths)
            page_index = div(i - 1, images_per_page) + 1
            position_in_page = (i - 1) % images_per_page

            if position_in_page == 0
                fig = Figure(size = (400, 400))
            end

            row = div(position_in_page, 2) + 1
            col = mod(position_in_page, 2) + 1

            ax = Axis(fig[row, col];
                xticksvisible = false, 
                yticksvisible = false,
                xgridvisible = false, 
                ygridvisible = false,
                leftspinevisible = false, 
                rightspinevisible = false,
                topspinevisible = false, 
                bottomspinevisible = false,
                xlabelvisible = false, 
                ylabelvisible = false,
                titlevisible = false,
            )

            image!(ax, load(path))

            if (i % images_per_page == 0) || (i == length(png_paths))
                push!(pages, fig)
                pdf_full_path = joinpath(output_pdf_dir, "$(pulsar_name)_page_$(page_index).pdf")
                pdf_name = "$(pulsar_name)_page_$(page_index).pdf"
                CairoMakie.save(pdf_full_path, fig)
                push!(pdf_paths, pdf_name)
            end
        end

        # Combine all PDF pages into one file using pdfunite
        combined_pdf_path = joinpath(output_pdf_dir, "combined.pdf")
        cmd_args = vcat(["pdfunite"], pdf_paths, [combined_pdf_path])
        println("Running command: ", join(cmd_args, " "))
        write("convert.txt", join(cmd_args, " "))
        #run(`$(cmd_args...)`)
        #println("Combined PDF created at: $combined_pdf_path")
    end


    """
        clean(data; threshold=0.7)

    Filter polarization data based on the ratio of linear polarization to intensity.

    # Arguments
    - `data`: 3D array with dimensions (num, bins, npol) containing polarization data
    - `data[:, :, 1]`: Stokes I intensity
    - `data[:, :, 2]`: Stokes Q parameter  
    - `data[:, :, 3]`: Stokes U parameter
    - `threshold`: filtering threshold (default 0.7) - minimum value for L/I ratio

    # Returns
    - `data_clean`: 2D array containing only intensity I for points meeting the L/I > threshold criterion

    # Description
    The function calculates linear polarization L = √(Q² + U²) for each data point and retains 
    only those points where the ratio of linear polarization to intensity exceeds the given threshold.
    """
    function clean(data; threshold=0.7)
        # four polarisation data
        num, bins, npol = size(data)
        data_clean = zeros(num, bins)
        for i in 1:num
            for j in 1:bins 
                I = data[i, j, 1]
                Q = data[i, j, 2]
                U = data[i, j, 3]
                L = sqrt(Q^2 + U^2)
                if L / I > threshold
                    data_clean[i, j] = I
                end
            end
        end
        return data_clean
    end


    """
    just some quick cleaning..
    """
    function remove_folders(dirname)
        # Get all subdirectories in dirname
        dirs = filter(isdir, readdir(dirname, join=true))
        for dir in dirs
            # Get all files in the subdirectory
            files = filter(isfile, readdir(dir, join=true))
            filename = basename(files[1])
            if length(files)==1 && filename=="params.json"
                println(basename(dir))
                rm(dir, recursive=true)
            end
        end
    end

    """
    some more cleaning
    """
    function remove_notinteresting(interesting_filename, data_dir)

        keep_dirs = Set(readlines(interesting_filename))
        all_dirs = filter(isdir, readdir(data_dir, join=true))

        for dir in all_dirs
            dir_name = basename(dir)
            if dir_name ∉ keep_dirs
                rm(dir, recursive=true, force=true)
            end
        end
        
    end


    """
    Process multifrequency data with PSRCHIVE and PSRSALSA
    """
    function process_psrdata_16(indir, outdir; files=nothing, outfile="pulsar.spCf16", params_file="params.json")

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)

        # if no catalog
        if !isdir(outdir)
            mkdir(outdir)
        end

        # check if params_file exists if not creating default one
        if !isfile(params_file)
            println("File $params_file does not exist, creating default one.")
            p = Tools.default_params(params_file)
        else
            p = Tools.read_params(params_file)
        end

        # add all .spCf16 files 
        add_psrfiles(indir, outfile; files=files, sixteen=true)

        if !isnothing(p["zaps"])
            println("ZAPPING subints starts")
            run(pipeline(`paz -m -w "$(p["zaps"])" $outfile`)) # zero weights  # this is the way
            println("ZAPPING subints ends")
        end
        #run(pipeline(`paz -m -W "49 194" -W "532 580" -W "880 982" $outfile`)) # zero weights  # this is the way

        # divide to two frequencies
        low_filename, high_filename, mid_filename = multifrequency_split(outfile)

        # gets number of single pulses if needed
        if isnothing(p["nsubint"])
            get_nsubint(low_filename, params_file, p)
        end

        # debase the data
        debase_16(low_filename, high_filename, mid_filename, params_file, p)

        # Calculate 2dfs and lrfs => finds P3
        low_debase = replace(low_filename, ".low"=>"_low.debase.gg")
        twodfs_lrfs(low_debase, params_file, p; detect=false)
        high_debase = replace(high_filename, ".high"=>"_high.debase.gg")
        twodfs_lrfs(high_debase, params_file, p; detect=false)

        # single pulses - low
        da_lo = Data.load_ascii(replace(low_filename, ".low"=>"_low_debase.txt"))
        Plot.single(da_lo, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="low", show_=true)

        # single pulses - high
        da_hi = Data.load_ascii(replace(high_filename, ".high"=>"_high_debase.txt"))
        Plot.single(da_hi, outdir; darkness=0.7, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="high", show_=true)

        # low freq p3-fold with cleaning
        #=
        da_low = Data.load_ascii_all(replace(low_filename, ".low"=>"_low_debase.txt"))
        da_low = clean(da_low; threshold=p["clean_threshold"]) 
        folded_low = Tools.p3fold(da_low, p["p3"],  p["p3_ybins"])
        Plot.single(folded_low, outdir; darkness=0.9, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold_low", show_=true, repeat_num=4)
        =#

        # high freq p3-fold with cleaning
        #=
        da_high = Data.load_ascii_all(replace(high_filename, ".high"=>"_high_debase.txt"))
        da_high = clean(da_high; threshold=p["clean_threshold"]) 
        folded_high = Tools.p3fold(da_high, p["p3"],  p["p3_ybins"])
        Plot.single(folded_high, outdir; darkness=0.9, number=nothing, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="p3fold_high", show_=true, repeat_num=4)
        =#

        # PSRSALSA p3folding low freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low", show_=true, repeat_num=4)
        # norefine
        #println("pfold -p3fold_norefine -p3fold \"$(p["p3"]) $(p["p3_ybins"])\" -onpulse \"$(p["bin_st"]) $(p["bin_end"])\" -onpulsed \"/NULL\" -p3foldd \"/NULL\" -w -oformat ascii $low_debase")
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_norefine", show_=true, repeat_num=4)
 
        # PSRSALSA p3folding high freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        #run(pipeline(`pfold -p3fold_nritt 50 -p3fold_cpb 50 -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high", show_=true, repeat_num=4)
        # norefine
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, outdir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_norefine", show_=true, repeat_num=4)

        return p, nothing, outdir

    end


    function multifrequency_split(spCf16_file)
        println(spCf16_file)

        low = replace(spCf16_file, ".spCf16"=>".low")
        low_txt = replace(spCf16_file, ".spCf16"=>"_low.txt")
        high = replace(spCf16_file, ".spCf16"=>".high")
        high_txt = replace(spCf16_file, ".spCf16"=>"_high.txt")
        mid = replace(spCf16_file, ".spCf16"=>".mid")
        mid_txt = replace(spCf16_file, ".spCf16"=>"_mid.txt")

        #run(pipeline(`paz -Z 0-7 -e high $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 8-15 -e low $spCf16_file`,stderr="errs.txt"))
        
        # 3 channels  
        run(pipeline(`paz -Z 0-12 -e high $spCf16_file`,stderr="errs.txt"))
        run(pipeline(`paz -Z 3-15 -e low $spCf16_file`,stderr="errs.txt"))
        run(pipeline(`paz -Z 0-6 -Z 9-15 -e mid $spCf16_file`,stderr="errs.txt")) # 2 channels?

        # 6 channels  
        #run(pipeline(`paz -Z 0-9 -e high $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 6-15 -e low $spCf16_file`,stderr="errs.txt"))
        #run(pipeline(`paz -Z 0-4 -Z 11-15 -e mid $spCf16_file`,stderr="errs.txt"))

        # fscrunch
        run(pipeline(`pam -F -e high $high`,stderr="errs.txt"))
        run(pipeline(`pam -F -e low $low`,stderr="errs.txt"))
        run(pipeline(`pam -F -e mid $mid`,stderr="errs.txt"))

        run(pipeline(`pdv -A -F $high`, stdout=high_txt, stderr="errs.txt"))
        run(pipeline(`pdv -A -F $low`, stdout=low_txt, stderr="errs.txt"))
        run(pipeline(`pdv -A -F $mid`, stdout=mid_txt, stderr="errs.txt"))
        #println("done")
        # change -t to -A to get frequancy information
        return low, high, mid

    end


    """
    Debase the data

    # Arguments
    """
    function debase_16(low, high, mid, params_file, params)

        # finding bin_st and bin_end based on low frequency file
        if isnothing(params["bin_st"]) || isnothing(params["bin_end"])

            println("pmod -device \"/xw\" -debase $low")
            run(pipeline(`pmod -device "/xw" -debase $low`, `tee pmod_output.txt`))
            #run(pipeline(`pmod -device "/xw" -iformat PSRFITS -debase $outfile`, `tee pmod_output.txt`))

            # Read captured output to get bin_st and bin_end
            output = read("pmod_output.txt", String)
            rm("pmod_output.txt")  # cleanup

            # Extract onpulse values
            m = match(r"-onpulse '(.+) (.+)'", output)
            if !isnothing(m)
                bin_st, bin_end = parse.(Int, m.captures)
                # Check if onpulse region length is even
                region_length = bin_end - bin_st + 1
                if region_length % 2 != 0
                    println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                    bin_end -= 1
                    println("Adjusted onpulse range: $bin_st to $bin_end")
                end
                println("Found onpulse range: $bin_st to $bin_end")
                params["bin_st"] = bin_st
                params["bin_end"] = bin_end
                Tools.save_params(params_file, params)
            end
            # ASCII single pulses
            convert_psrfit_ascii(replace(low, ".low"=>".debase.gg"), replace(low, ".low"=>"_low_debase.txt"))

        else
            run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $low`))
            #run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -iformat PSRFITS -debase $outfile`))
            # ASCII single pulses
            convert_psrfit_ascii(replace(low, ".low"=>".debase.gg"), replace(low, ".low"=>"_low_debase.txt"))
        end

        # changing low frequency psrfit filename
        mv(replace(low, ".low"=>".debase.gg"), replace(low, "pulsar.low"=>"pulsar_low.debase.gg"), force=true)

        # high frequency runs
        run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $high`))
        convert_psrfit_ascii(replace(high, ".high"=>".debase.gg"), replace(high, ".high"=>"_high_debase.txt"))
        # changing high frequency psrfit filename
        mv(replace(high, ".high"=>".debase.gg"), replace(high, "pulsar.high"=>"pulsar_high.debase.gg"), force=true)

        # mid frequency runs
        run(pipeline(`pmod -onpulse "$(params["bin_st"]) $(params["bin_end"])" -device "/NULL" -debase $mid`))
        convert_psrfit_ascii(replace(mid, ".mid"=>".debase.gg"), replace(mid, ".mid"=>"_mid_debase.txt"))
        # changing high frequency psrfit filename
        mv(replace(mid, ".mid"=>".debase.gg"), replace(mid, "pulsar.mid"=>"pulsar_mid.debase.gg"), force=true)

    end


    """
    Analyse p3folds

    # Arguments
    - `indir`: data input directory
    - `type`: type of p3folds "refine" or "norefine" filenames are based on that
 
    """
    function analyse_p3folds_16(indir, type)

        # parameters file 
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filename end
        low = joinpath(indir, "pulsar_low.debase.p3fold_" * type)
        high = joinpath(indir, "pulsar_high.debase.p3fold_" * type)

        #println(low)
        l = Data.load_ascii(low)
        Plot.p3fold(l, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_$(type)_analyse", show_=true, repeat_num=1)

        #println(high)
        h = Data.load_ascii(high)
        Plot.p3fold(h, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_$(type)_analyse", show_=true, repeat_num=1)
 
        #diff = normalize_01(l) .- normalize_01(h)
        diff = normalize_per_pulse(l) .- normalize_per_pulse(h)
        Plot.p3fold(diff, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_high_$(type)_analyse", show_=true, repeat_num=1)

        nl = normalize_per_pulse(l)
        nh = normalize_per_pulse(h)
        shifts, corr_values = find_shift_per_pulse(nl, nh)

        println("Mean shift: $(mean(shifts)) bins")
        println("Median shift: $(median(shifts)) bins")

        profile_idx = 1
        p1 = nl[profile_idx, :]
        p2 = nh[profile_idx, :]
        bins = 1:length(p1)

        writedlm(joinpath(indir,"p_low.txt"), hcat(bins, p1), ' ')
        writedlm(joinpath(indir,"p_high.txt"), hcat(bins, p2), ' ')


        println(size(nl))

        Plot.analyse_p3folds(nl, nh, p)
        Plot.analyse_p3folds2(nl, nh, p)

    end


    """
    Analyse p3folds

    # Arguments
    - `indir`: data input directory
    - `type`: type of p3folds "refine" or "norefine" filenames are based on that
 
    """
    function analyse_p3folds_16_new(indir, type; n_comp=3, psr=nothing, save_summary_dir=nothing)

        # parameters file
        p = Tools.read_params(joinpath(indir, "params.json"))
        psr = something(psr, psr_from_dir(indir))

        # low, high frequancy filename end
        low = joinpath(indir, "pulsar_low.debase.p3fold_" * type)
        high = joinpath(indir, "pulsar_high.debase.p3fold_" * type)

        #println(low)
        l = Data.load_ascii(low)
        #Plot.p3fold(l, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_$(type)_analyse", show_=true, repeat_num=1)

        #println(high)
        h = Data.load_ascii(high)
        #Plot.p3fold(h, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_$(type)_analyse", show_=true, repeat_num=1)
 
        #diff = normalize_01(l) .- normalize_01(h)
        diff = normalize_per_pulse(l) .- normalize_per_pulse(h)
        #Plot.p3fold(diff, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_high_$(type)_analyse", show_=true, repeat_num=1)

        nl = normalize_per_pulse(l)
        nh = normalize_per_pulse(h)
        #=
        profile_idx = 1
        p1 = nl[profile_idx, :]
        p2 = nh[profile_idx, :]
        bins = 1:length(p1)

        writedlm(joinpath(indir,"p_low.txt"), hcat(bins, p1), ' ')
        writedlm(joinpath(indir,"p_high.txt"), hcat(bins, p2), ' ')
        =#
        # TODO work here...

        #Plot.analyse_p3folds3(nl, nh, p, n_comp) # old
        return Plot.analyse_p3folds4(nl, nh, p, n_comp; psr=psr, save_summary_dir=save_summary_dir) # with 's' to skip added


        #println(size(nl))


    end


    """
    Offset analysis on single-pulse data.

    Groups single pulses (pulsar_low.txt / pulsar_high.txt) into blocks of
    npulse pulses each and averages every block into one mean profile.
    The last block may contain fewer pulses if npulse does not divide evenly.
    Then fits Gaussian components to measure the frequency-dependent offset.

    # Arguments
    - `indir`:  data directory (must contain params.json, pulsar_low.txt, pulsar_high.txt)
    - `n_comp`: number of Gaussian components to fit (default 3)
    - `npulse`: number of single pulses per averaged profile (default 150)
    """
    function analyse_average_offset(indir; n_comp=3, npulse=150, psr=nothing)

        p = Tools.read_params(joinpath(indir, "params.json"))
        psr = something(psr, psr_from_dir(indir))

        l_all = Data.load_ascii_all(joinpath(indir, "pulsar_low.txt"))
        h_all = Data.load_ascii_all(joinpath(indir, "pulsar_high.txt"))
        l = l_all[:, :, 1]
        h = h_all[:, :, 1]

        n_pulses = size(l, 1)
        bins     = size(l, 2)
        @assert size(h, 1) == n_pulses "Low and high pulse counts differ ($(n_pulses) vs $(size(h,1)))"

        n_profiles = cld(n_pulses, npulse)

        nl_avg = zeros(n_profiles, bins)
        nh_avg = zeros(n_profiles, bins)

        for i in 1:n_profiles
            st = (i - 1) * npulse + 1
            en = min(i * npulse, n_pulses)
            prof_l = vec(mean(l[st:en, :], dims=1))
            prof_h = vec(mean(h[st:en, :], dims=1))
            prof_l .-= minimum(prof_l)
            prof_h .-= minimum(prof_h)
            nl_avg[i, :] = prof_l / maximum(prof_l)
            nh_avg[i, :] = prof_h / maximum(prof_h)
        end

        println("Created $n_profiles averaged profiles from $n_pulses pulses (npulse = $npulse, last block = $(n_pulses - (n_profiles-1)*npulse))")

        Plot.analyse_average_offset(nl_avg, nh_avg, p, n_comp; npulse=npulse, n_pulses=n_pulses, psr=psr)

    end


    """
    Pulsar name out of a data directory such as ".../J1842-0359_uhf_16" — the
    directories carry a suffix (_16, _16c, _uhf_16) that the catalogue and
    input/offsets.csv do not, so only the PSRJ part is kept. Returns "" when
    the directory does not look like a pulsar at all, which leaves the
    separations file untouched rather than writing a nonsense row.
    """
    function psr_from_dir(dir)
        m = match(r"J\d{4}[+-]\d{2,4}", basename(rstrip(dir, '/')))
        return isnothing(m) ? "" : m.match
    end


    function normalize_01(data)
        min_val = minimum(data)
        max_val = maximum(data)
        return (data .- min_val) ./ (max_val - min_val)
    end

    function normalize_per_pulse(data)
        normalized = similar(data)
        for i in 1:size(data, 1)
            pulse = data[i, :]
            min_val = minimum(pulse)
            max_val = maximum(pulse)
            normalized[i, :] = (pulse .- min_val) ./ (max_val - min_val)
        end
        return normalized
    end


    function find_shift_per_pulse(data1, data2)
        n_pulses = size(data1, 1)
        shifts = zeros(Int, n_pulses)
        correlations = zeros(n_pulses)
        
        for i in 1:n_pulses
            # Cross-correlation
            ccf = ifft(fft(data1[i, :]) .* conj.(fft(data2[i, :])))
            ccf = real.(ccf)
            
            # Znajdź maksimum
            max_idx = argmax(ccf)
            correlations[i] = ccf[max_idx]
            
            # Przelicz na przesunięcie (z uwzględnieniem wraparound)
            n_bins = length(ccf)
            shifts[i] = max_idx - 1
            if shifts[i] > n_bins ÷ 2
                shifts[i] = shifts[i] - n_bins
            end
        end
        
        return shifts, correlations
    end


    """
    Process multifrequency data with PSRCHIVE and PSRSALSA single pulse analysis
    """
    function process_psrdata_single(indir, outdir; files=nothing, outfile="pulsar.spCf16", params_file="params.json")

        # full paths
        params_file = joinpath(outdir, params_file)
        outfile = joinpath(outdir, outfile)

        # if no catalog
        if !isdir(outdir)
            mkdir(outdir)
        end

        # check if params_file exists if not creating default one
        if !isfile(params_file)
            println("File $params_file does not exist, creating default one.")
            p = Tools.default_params(params_file)
        else
            p = Tools.read_params(params_file)
        end

        # speeding up
        # add all .spCf16 files 
        add_psrfiles(indir, outfile; files=files, sixteen=true)

        # divide to two frequencies
        low_filename, high_filename, mid_filename = multifrequency_split(outfile)

        # gets number of single pulses if needed
        if isnothing(p["nsubint"])
            get_nsubint(low_filename, params_file, p)
        end

        # debase the data
        debase_16(low_filename, high_filename, mid_filename, params_file, p)

        # speeding ends here

        # comment those three if speeding up is off
        #low_filename = "/home/psr/output/J2139+2242_single/pulsar.low"
        #high_filename = "/home/psr/output/J2139+2242_single/pulsar.high"
        #mid_filename = "/home/psr/output/J2139+2242_single/pulsar.mid"

        low_debase = replace(low_filename, ".low"=>"_low.debase.gg")
        Tools.generate_snr(low_debase) 
        low_snrfile = low_debase * ".snr.txt"

        twodfs_lrfs(low_debase, params_file, p; detect=false) # just to detect p2

        # single pulses low
        da_lo = Data.load_ascii(replace(low_filename, ".low"=>"_low_debase.txt"))

        high_debase = replace(high_filename, ".high"=>"_high.debase.gg")
        Tools.generate_snr(high_debase) 
        high_snrfile = high_debase * ".snr.txt"
 
        # single pulses high 
        da_hi = Data.load_ascii(replace(high_filename, ".high"=>"_high_debase.txt"))

        # single pulses mid 
        da_mi = Data.load_ascii(replace(mid_filename, ".mid"=>"_mid_debase.txt"))

        Tools.track_subpulses_snr3(da_lo, da_hi, da_mi, 10, low_snrfile, on_st=p["bin_st"], on_end=p["bin_end"])
        
        Plot.analyse_single(da_lo, da_hi, da_mi, p; pulse_start=250)

        # TODO single pulse detection work starts here
        #Tools.track_subpulses_manual(da_mi, on_st=p["bin_st"], on_end=p["bin_end"])

        return 
    end

    function p3fold_psrdata(indir)

        # parameters file 
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filenames
        low = joinpath(indir, "pulsar.low")
        high = joinpath(indir, "pulsar.high")

        low_debase = replace(low, ".low"=>"_low.debase.gg")
        high_debase = replace(high, ".high"=>"_high.debase.gg")
        
        #-p3fold_norefine -p3fold_smooth 3 -p3fold_nritt 1 -p3fold_cpb 1

        
        # PSRSALSA p3folding low freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low", show_=true, repeat_num=4)
        # norefine
        #println("pfold -p3fold_norefine -p3fold \"$(p["p3"]) $(p["p3_ybins"])\" -onpulse \"$(p["bin_st"]) $(p["bin_end"])\" -onpulsed \"/NULL\" -p3foldd \"/NULL\" -w -oformat ascii $low_debase")
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $low_debase`,  stderr="errs.txt"))
        mv(replace(low_debase, ".gg"=>".p3fold"), replace(low_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(low_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_low_norefine", show_=true, repeat_num=4)
 
        # PSRSALSA p3folding high freq
        run(pipeline(`pfold -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_refine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_refine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high", show_=true, repeat_num=4)
        # norefine
        run(pipeline(`pfold -p3fold_norefine -p3fold "$(p["p3"]) $(p["p3_ybins"])" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $high_debase`,  stderr="errs.txt"))
        mv(replace(high_debase, ".gg"=>".p3fold"), replace(high_debase, ".gg"=>".p3fold_norefine"), force=true)
        d5 = Data.load_ascii(replace(high_debase, ".gg"=>".p3fold_norefine"))
        Plot.p3fold(d5, indir; start=1, bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=0.9, name_mod="pulsar_high_norefine", show_=true, repeat_num=4)

        #=
        =#




    end



    """
    Remove orthogonal polarization mode (OPM) jumps from a PA track by
    unwrapping along longitude with a 90° period. Walks bin-by-bin keeping
    a running "continuous" PA; whenever the next sample lies farther than
    `tolerance` degrees from it, shift the sample by ±90° (the OPM step)
    until it rejoins the track. The accumulated offset persists, so a whole
    OPM island is corrected end-to-end regardless of length. Finally the
    output is wrapped back to (-90°, 90°] — which also undoes any shifts
    of ±180° that correspond to the intrinsic mod-π PA wrap, leaving those
    pixels unchanged.
    """
    function deopm_pa(pa_deg; tolerance=45.0)
        pa_out = copy(pa_deg)
        n = length(pa_out)
        last_unwrapped = NaN
        offset = 0.0
        flipped = 0
        for i in 1:n
            isnan(pa_out[i]) && continue
            pa_in = pa_out[i]
            if isnan(last_unwrapped)
                last_unwrapped = pa_in
                continue
            end
            shifted = pa_in + offset
            d = shifted - last_unwrapped
            while d > tolerance
                offset -= 90; shifted -= 90; d -= 90
            end
            while d < -tolerance
                offset += 90; shifted += 90; d += 90
            end
            pa_out[i] = mod(shifted + 90, 180) - 90
            abs(pa_out[i] - pa_in) > 0.1 && (flipped += 1)
            last_unwrapped = shifted
        end
        return pa_out, flipped
    end


    """
    Fit the Rotating Vector Model (RVM) to position angle data.

    Grid search over α (inclination) and β (impact parameter), then for each
    (α, β) pair solves linearly for PA0 and searches over φ0 (inflection longitude).

    Returns NamedTuple with fields: alpha, beta, phi0, pa0 (all in degrees), chi2.
    Returns nothing if fewer than 5 valid PA points.
    """
    function fit_rvm(lon_deg, pa_deg, pa_err_deg;
                     return_map=false,
                     n_alpha=85, n_beta=40, n_phi0=60,
                     alpha_range=(5.0, 175.0), beta_range=(-20.0, 20.0))
        mask = .!isnan.(pa_deg) .& .!isnan.(pa_err_deg) .& (pa_err_deg .> 0)
        if sum(mask) < 5
            return nothing
        end

        lon_rad = lon_deg[mask]     .* (π / 180)
        pa_rad  = pa_deg[mask]      .* (π / 180)
        sigma   = pa_err_deg[mask]  .* (π / 180)
        inv_var = 1.0 ./ (sigma .^ 2)

        alphas = range(alpha_range[1], alpha_range[2], length=n_alpha) .* (π / 180)
        betas  = range(beta_range[1],  beta_range[2],  length=n_beta)  .* (π / 180)
        # Search φ0 over the full longitude range
        lon_min, lon_max = extrema(lon_rad)
        margin = (lon_max - lon_min) * 0.5
        phi0s = range(lon_min - margin, lon_max + margin, length=n_phi0)

        best_chi2  = Inf
        best_alpha = 0.0
        best_beta  = 0.0
        best_phi0  = 0.0
        best_pa0   = 0.0

        na = length(alphas)
        nb = length(betas)
        chi2_map = return_map ? fill(Inf, na, nb) : nothing

        for (ia, alpha) in enumerate(alphas)
            ca = cos(alpha); sa = sin(alpha)
            for (ib, beta) in enumerate(betas)
                zeta = alpha + beta
                cz = cos(zeta);  sz = sin(zeta)
                best_ab  = Inf
                best_phi0_ab = phi0s[1]
                best_pa0_ab  = 0.0
                for phi0 in phi0s
                    dphi = lon_rad .- phi0
                    # Komesaroff (1970) sign convention (matches publication).
                    rvm_shape = atan.(.-sa .* sin.(dphi),
                                      sz .* ca .- cz .* sa .* cos.(dphi))

                    # PA is defined mod π. Use a weighted circular mean on
                    # 2·diff to get PA0, then wrap residuals to (-π/2, π/2].
                    diffs = pa_rad .- rvm_shape
                    Sx = sum(sin.(2 .* diffs) .* inv_var)
                    Cx = sum(cos.(2 .* diffs) .* inv_var)
                    pa0 = atan(Sx, Cx) / 2
                    residuals = mod.(diffs .- pa0 .+ π/2, π) .- π/2
                    chi2 = sum((residuals .^ 2) .* inv_var)

                    if chi2 < best_ab
                        best_ab      = chi2
                        best_phi0_ab = phi0
                        best_pa0_ab  = pa0
                    end
                end
                if return_map
                    chi2_map[ia, ib] = best_ab
                end
                if best_ab < best_chi2
                    best_chi2  = best_ab
                    best_alpha = alpha
                    best_beta  = beta
                    best_phi0  = best_phi0_ab
                    best_pa0   = best_pa0_ab
                end
            end
        end

        ndof = max(sum(mask) - 4, 1)
        result = (alpha    = best_alpha * (180/π),
                  beta     = best_beta  * (180/π),
                  phi0     = best_phi0  * (180/π),
                  pa0      = best_pa0   * (180/π),
                  chi2     = best_chi2,
                  chi2_red = best_chi2 / ndof,
                  ndof     = ndof)
        if return_map
            return result, chi2_map, collect(alphas .* (180/π)), collect(betas .* (180/π))
        else
            return result
        end
    end


    """
    Evaluate RVM curve over a dense longitude grid (degrees in, degrees out).
    Returns (lon_dense, pa_rvm, pa_rvm_ortho) where pa_rvm_ortho is the 90° mode.
    """
    function rvm_curve(params, lon_min_deg, lon_max_deg; npts=500)
        lon = collect(range(lon_min_deg, lon_max_deg, length=npts))
        lon_r = lon .* (π / 180)
        phi0  = params.phi0 * (π / 180)
        pa0   = params.pa0  * (π / 180)
        alpha = params.alpha * (π / 180)
        beta  = params.beta  * (π / 180)
        zeta  = alpha + beta

        dphi = lon_r .- phi0
        # Komesaroff (1970) sign convention (matches publication).
        pa_r = pa0 .+ atan.(.-sin(alpha) .* sin.(dphi),
                             sin(zeta) .* cos(alpha) .- cos(zeta) .* sin(alpha) .* cos.(dphi))
        pa_deg       = mod.(pa_r .* (180/π) .+ 90, 180) .- 90
        pa_ortho_deg = mod.(pa_deg .+ 90, 180) .- 90

        # Break line at wrap jumps so renderers don't draw vertical segments
        idx = eachindex(pa_deg)
        for i in Iterators.drop(idx, 1)
            if abs(pa_deg[i] - pa_deg[i-1]) > 90
                pa_deg[i-1] = NaN
            end
            if abs(pa_ortho_deg[i] - pa_ortho_deg[i-1]) > 90
                pa_ortho_deg[i-1] = NaN
            end
        end

        return lon, pa_deg, pa_ortho_deg
    end


    function position_angle(indir)
        # parameters file
        p = Tools.read_params(joinpath(indir, "params.json"))

        # low, high frequancy filenames
        low = joinpath(indir, "pulsar.low")
        high = joinpath(indir, "pulsar.high")

        # txt files
        lt = joinpath(indir, "pulsar_low.txt") # no dabese here needed?
        ht = joinpath(indir, "pulsar_high.txt") # no dabese here needed?

        l = Data.load_ascii_all(lt)
        h = Data.load_ascii_all(ht)

        bin_st  = p["bin_st"]
        bin_end = p["bin_end"]

        pulses_l, bins_l, _ = size(l)
        pulses_h, bins_h, _ = size(h)

        # Longitude axis (centered at 0)
        db_l = (bin_end + 1) - bin_st
        dl_l = 360.0 * db_l / bins_l
        lon_l = collect(range(-dl_l/2.0, dl_l/2.0, length=db_l))

        db_h = (bin_end + 1) - bin_st
        dl_h = 360.0 * db_h / bins_h
        lon_h = collect(range(-dl_h/2.0, dl_h/2.0, length=db_h))

        # Average Stokes profiles over on-pulse window
        I_l = vec(mean(l[:, bin_st:bin_end, 1], dims=1))
        Q_l = vec(mean(l[:, bin_st:bin_end, 2], dims=1))
        U_l = vec(mean(l[:, bin_st:bin_end, 3], dims=1))
        V_l = vec(mean(l[:, bin_st:bin_end, 4], dims=1))
        Lin_l = sqrt.(Q_l.^2 .+ U_l.^2)

        I_h = vec(mean(h[:, bin_st:bin_end, 1], dims=1))
        Q_h = vec(mean(h[:, bin_st:bin_end, 2], dims=1))
        U_h = vec(mean(h[:, bin_st:bin_end, 3], dims=1))
        V_h = vec(mean(h[:, bin_st:bin_end, 4], dims=1))
        Lin_h = sqrt.(Q_h.^2 .+ U_h.^2)

        # PA masked where L < 5σ on averaged Q, U (Johnston et al. 2023)
        off_noise_l = std(l[:, 1:bin_st-1, 1])
        off_noise_h = std(h[:, 1:bin_st-1, 1])
        sigma_avg_l = off_noise_l / sqrt(pulses_l)
        sigma_avg_h = off_noise_h / sqrt(pulses_h)
        thresh_l = 5.0 * sigma_avg_l
        thresh_h = 5.0 * sigma_avg_h

        pa_l = [Lin_l[i] > thresh_l ? 0.5 * atan(U_l[i], Q_l[i]) * (180.0/pi) : NaN for i in 1:db_l]
        pa_h = [Lin_h[i] > thresh_h ? 0.5 * atan(U_h[i], Q_h[i]) * (180.0/pi) : NaN for i in 1:db_h]

        # PA errors: σ_PA = 0.5 * σ_noise_avg / L  (in degrees)
        pa_err_l = [Lin_l[i] > thresh_l ? 0.5 * sigma_avg_l / Lin_l[i] * (180.0/pi) : NaN for i in 1:db_l]
        pa_err_h = [Lin_h[i] > thresh_h ? 0.5 * sigma_avg_h / Lin_h[i] * (180.0/pi) : NaN for i in 1:db_h]

        # Detect and undo orthogonal polarization mode jumps before RVM fit
        pa_l, flipped_l = deopm_pa(pa_l)
        pa_h, flipped_h = deopm_pa(pa_h)
        println("De-OPM: flipped $flipped_l (low) and $flipped_h (high) bins")

        # Select OPM branch (publication convention): global +90° shift
        pa_shift = 90.0
        pa_l = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_l]
        pa_h = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_h]

        # Fit RVM to low-frequency PA (more points typically)
        println("Fitting RVM (low frequency)...")
        rvm_params_l = fit_rvm(lon_l, pa_l, pa_err_l)
        if !isnothing(rvm_params_l)
            println("  α = $(round(rvm_params_l.alpha, digits=1))°, " *
                    "β = $(round(rvm_params_l.beta, digits=1))°, " *
                    "φ₀ = $(round(rvm_params_l.phi0, digits=2))°, " *
                    "PA₀ = $(round(rvm_params_l.pa0, digits=1))°, " *
                    "χ²/ndof = $(round(rvm_params_l.chi2_red, digits=2))")
            lon_rvm_l, pa_rvm_l, pa_rvm_l_ortho = rvm_curve(rvm_params_l,
                                                              minimum(lon_l), maximum(lon_l))
        else
            println("  Not enough PA points for RVM fit (low)")
            lon_rvm_l = pa_rvm_l = pa_rvm_l_ortho = nothing
        end

        println("Fitting RVM (high frequency)...")
        rvm_params_h = fit_rvm(lon_h, pa_h, pa_err_h)
        if !isnothing(rvm_params_h)
            println("  α = $(round(rvm_params_h.alpha, digits=1))°, " *
                    "β = $(round(rvm_params_h.beta, digits=1))°, " *
                    "φ₀ = $(round(rvm_params_h.phi0, digits=2))°, " *
                    "PA₀ = $(round(rvm_params_h.pa0, digits=1))°, " *
                    "χ²/ndof = $(round(rvm_params_h.chi2_red, digits=2))")
            lon_rvm_h, pa_rvm_h, pa_rvm_h_ortho = rvm_curve(rvm_params_h,
                                                              minimum(lon_h), maximum(lon_h))
        else
            println("  Not enough PA points for RVM fit (high)")
            lon_rvm_h = pa_rvm_h = pa_rvm_h_ortho = nothing
        end

        phi0_l = isnothing(rvm_params_l) ? nothing : rvm_params_l.phi0
        phi0_h = isnothing(rvm_params_h) ? nothing : rvm_params_h.phi0

        Plot.position_angle(lon_l, pa_l, pa_err_l, I_l, Lin_l, V_l,
                            lon_h, pa_h, pa_err_h, I_h, Lin_h, V_h,
                            indir; show_=true,
                            lon_rvm_l=lon_rvm_l, pa_rvm_l=pa_rvm_l,
                            pa_rvm_l_ortho=pa_rvm_l_ortho, phi0_l=phi0_l,
                            lon_rvm_h=lon_rvm_h, pa_rvm_h=pa_rvm_h,
                            pa_rvm_h_ortho=pa_rvm_h_ortho, phi0_h=phi0_h)
    end

    function geometry_analysis(indir; chi2_red_max=10.0)

        # parameters file
        p = Tools.read_params(joinpath(indir, "params.json"))

        lt = joinpath(indir, "pulsar_low.txt")
        ht = joinpath(indir, "pulsar_high.txt")

        l = load_ascii_all(lt)
        h = load_ascii_all(ht)

        bin_st  = p["bin_st"]
        bin_end = p["bin_end"]

        pulses_l, bins_l, _ = size(l)
        pulses_h, bins_h, _ = size(h)

        db_l = (bin_end + 1) - bin_st
        lon_l = collect(range(-360.0 * db_l / bins_l / 2, 360.0 * db_l / bins_l / 2, length=db_l))
        db_h = (bin_end + 1) - bin_st
        lon_h = collect(range(-360.0 * db_h / bins_h / 2, 360.0 * db_h / bins_h / 2, length=db_h))

        Q_l   = vec(mean(l[:, bin_st:bin_end, 2], dims=1))
        U_l   = vec(mean(l[:, bin_st:bin_end, 3], dims=1))
        Lin_l = sqrt.(Q_l.^2 .+ U_l.^2)
        Q_h   = vec(mean(h[:, bin_st:bin_end, 2], dims=1))
        U_h   = vec(mean(h[:, bin_st:bin_end, 3], dims=1))
        Lin_h = sqrt.(Q_h.^2 .+ U_h.^2)

        # Mean Stokes I for W10 estimation
        I_l_full = vec(mean(l[:, :, 1], dims=1))
        I_h_full = vec(mean(h[:, :, 1], dims=1))

        # sigma_avg: statistically correct noise on the mean profile, used both
        # for detection threshold (5σ on L) and for PA error (0.5·σ/L). Combined
        # with the χ²_red-display mode in Plot.geometry (cutoff 10), this gives
        # the Johnston+ 2023 Fig. 1 convention: minimum χ²_red ≈ 1 at best fit,
        # the "χ² > 10" region masked as blank.
        sigma_avg_l = std(l[:, 1:bin_st-1, 1]) / sqrt(pulses_l)
        sigma_avg_h = std(h[:, 1:bin_st-1, 1]) / sqrt(pulses_h)
        thresh_l = 5.0 * sigma_avg_l
        thresh_h = 5.0 * sigma_avg_h

        pa_l     = [Lin_l[i] > thresh_l ? 0.5 * atan(U_l[i], Q_l[i]) * (180.0/π) : NaN for i in 1:db_l]
        pa_err_l = [Lin_l[i] > thresh_l ? 0.5 * sigma_avg_l / Lin_l[i] * (180.0/π) : NaN for i in 1:db_l]
        pa_h     = [Lin_h[i] > thresh_h ? 0.5 * atan(U_h[i], Q_h[i]) * (180.0/π) : NaN for i in 1:db_h]
        pa_err_h = [Lin_h[i] > thresh_h ? 0.5 * sigma_avg_h / Lin_h[i] * (180.0/π) : NaN for i in 1:db_h]

        # Match the preprocessing used in position_angle: undo OPM jumps,
        # then select the same OPM branch (global +90° shift).
        pa_l, _ = deopm_pa(pa_l)
        pa_h, _ = deopm_pa(pa_h)
        pa_shift = 90.0
        pa_l = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_l]
        pa_h = [isnan(x) ? x : mod(x + pa_shift + 90, 180) - 90 for x in pa_h]

        # Two-pass with floored rescaling: first fit with σ_avg, then inflate
        # errors by √χ²_red (if > 1) so χ²_red_min = 1 by construction. Floor at
        # 1 prevents deflating σ when the model fits below statistical noise
        # (e.g., sparse high-freq PA points) — there the region stays wide.
        println("Fitting RVM (low frequency, pass 1)...")
        res_l_pre = fit_rvm(lon_l, pa_l, pa_err_l;
            n_alpha=171, n_beta=161, n_phi0=200)
        scale_l = max(1.0, sqrt(res_l_pre.chi2_red))
        println("  χ²_red=$(round(res_l_pre.chi2_red,digits=2))  → σ scale factor = $(round(scale_l,digits=2))")

        println("Fitting RVM chi² map (low frequency, pass 2)...")
        res_l, chi2_l, alphas_deg, betas_deg = fit_rvm(lon_l, pa_l, pa_err_l .* scale_l;
            return_map=true, n_alpha=171, n_beta=161, n_phi0=200)
        println("  best: α=$(round(res_l.alpha,digits=1))° β=$(round(res_l.beta,digits=1))° χ²/ndof=$(round(res_l.chi2_red,digits=4))")

        println("Fitting RVM (high frequency, pass 1)...")
        res_h_pre = fit_rvm(lon_h, pa_h, pa_err_h;
            n_alpha=171, n_beta=161, n_phi0=200)
        scale_h = max(1.0, sqrt(res_h_pre.chi2_red))
        println("  χ²_red=$(round(res_h_pre.chi2_red,digits=2))  → σ scale factor = $(round(scale_h,digits=2))")

        println("Fitting RVM chi² map (high frequency, pass 2)...")
        res_h, chi2_h, _, _ = fit_rvm(lon_h, pa_h, pa_err_h .* scale_h;
            return_map=true, n_alpha=171, n_beta=161, n_phi0=200)
        println("  best: α=$(round(res_h.alpha,digits=1))° β=$(round(res_h.beta,digits=1))° χ²/ndof=$(round(res_h.chi2_red,digits=4))")

        # W10 (pulse width at 10 % of peak flux), cf. Johnston et al. 2023, Eq. (3)
        nbin  = p["nbin"]
        W10_l = pulse_width_fraction(I_l_full, nbin; frac=0.1,
                                     on_st=bin_st, on_end=bin_end)
        W10_h = pulse_width_fraction(I_h_full, nbin; frac=0.1,
                                     on_st=bin_st, on_end=bin_end)
        P_sec = p["period"]
        println("P = $P_sec s   W10: low=$(round(W10_l, digits=2))°  high=$(round(W10_h, digits=2))°")

        # Diagnostic: χ²_red percentiles — helps see if a seemingly "missing"
        # region is just above the cutoff or genuinely a bad fit.
        function _stats(chi2_map, ndof)
            red = chi2_map[isfinite.(chi2_map)] ./ ndof
            sort!(red)
            n = length(red)
            q(p) = red[clamp(round(Int, p*n), 1, n)]
            (min=red[1], p50=q(0.5), p90=q(0.9), p95=q(0.95), p99=q(0.99), max=red[end])
        end
        st_l = _stats(chi2_l, res_l.ndof)
        st_h = _stats(chi2_h, res_h.ndof)
        println("χ²_red low:  min=$(round(st_l.min,digits=2))  p50=$(round(st_l.p50,digits=2))  p90=$(round(st_l.p90,digits=2))  p95=$(round(st_l.p95,digits=2))  p99=$(round(st_l.p99,digits=2))  max=$(round(st_l.max,digits=2))")
        println("χ²_red high: min=$(round(st_h.min,digits=2))  p50=$(round(st_h.p50,digits=2))  p90=$(round(st_h.p90,digits=2))  p95=$(round(st_h.p95,digits=2))  p99=$(round(st_h.p99,digits=2))  max=$(round(st_h.max,digits=2))")

        Plot.geometry(chi2_l, chi2_h, alphas_deg, betas_deg, indir;
                      show_=true, P_sec=P_sec, W_deg_l=W10_l, W_deg_h=W10_h,
                      ndof_l=res_l.ndof, ndof_h=res_h.ndof,
                      chi2_red_max=chi2_red_max)
    end

    """
    W_frac (e.g. W10): distance between the outermost longitude bins above
    `frac * peak`, following Posselt et al. 2021 / Johnston et al. 2023.
    If `on_st`/`on_end` are given, the peak and search are restricted to that
    on-pulse window; this avoids noise spikes extending the width artificially.
    """
    function pulse_width_fraction(profile, nbin; frac=0.1,
                                   on_st=nothing, on_end=nothing)
        st = isnothing(on_st)  ? 1              : on_st
        en = isnothing(on_end) ? length(profile) : on_end
        seg = @view profile[st:en]
        threshold = frac * maximum(seg)
        above = findall(>=(threshold), seg)
        isempty(above) && return NaN
        return (last(above) - first(above) + 1) * 360.0 / nbin
    end

end # module
