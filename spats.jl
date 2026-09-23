module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")
    include("modules/phase_modulation.jl")
    include("modules/p3fold_viterbi.jl")
    include("modules/relations.jl")


    function test(outdir)
        d = Data.load_ascii("input/1.txt")
        Plot.single(d, outdir; darkness=0.3, number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.average(d, outdir; number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.lrfs_obsolete(d, outdir; darkness=0.1, start=1, name_mod="1", bin_st=500, bin_end=530, show_=true)
    end

    
    function test2(outdir)    
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00000-00255.spCF", outdir*"1.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00256-00511.spCF", outdir*"2.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00512-00767.spCF", outdir*"3.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00768-01029.spCF", outdir*"4.txt")  
        data1 = Data.load_ascii(outdir*"1.txt")
        data2 = Data.load_ascii(outdir*"2.txt")
        data3 = Data.load_ascii(outdir*"3.txt")
        data4 = Data.load_ascii(outdir*"4.txt")
        data = vcat(data1, data2, data3, data4)
        Plot.single(data, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        Plot.lrfs_obsolete(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        folded = Tools.p3fold(data, 20, 40)
        Plot.single(folded, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319_p3fold", show_=true)
    end


    function test3(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        Data.convert_psrfit_ascii(joinpath(outdir, "pulsar.debase.gg"), joinpath(outdir, "pulsar.debase.txt"))

        outfile="pulsar.spCF"
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        p = Tools.read_params(joinpath(outdir, "params.json"))
        d4 = Data.load_ascii_all(joinpath(outdir, "pulsar.debase.txt"))
        d1 = Data.clean(d4; threshold=0.0031)
        Plot.single(d1, outdir; darkness=0.7, number=100, bin_st=p["bin_st"], bin_end=p["bin_end"], start=210, name_mod="pulsar", show_=true)
        Plot.lrfs_obsolete(d1, outdir; darkness=0.3, start=210, name_mod="pulsar", bin_st=p["bin_st"], bin_end=p["bin_end"], show_=true)
        Data.twodfs_lrfs(debased_file, outdir, p)
        lrfs_file = replace(debased_file, "gg"=>"lrfs")
        data = Data.load_ascii_all(lrfs_file)
        Plot.lrfs(data, outdir, p; darkness=0.3, name_mod="pulsar", show_=true)
        twodfs_file = replace(debased_file, "gg"=>"1.2dfs")
        data_2dfs = Data.load_ascii_all(twodfs_file)
        Plot.twodfs(data_2dfs, outdir, p; darkness=0.3, name_mod="pulsar", show_=true)
        folded = Data.load_ascii(joinpath(outdir, "pulsar.debase.p3fold"))
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="pulsar", show_=true, repeat_num=4)        

    end

    """
    For Sardinia poster
    """
    function J1539_6322_Sard(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        Data.convert_psrfit_ascii(joinpath(outdir, "pulsar.debase.gg"), joinpath(outdir, "pulsar.debase.txt"))
        outfile="pulsar.spCF"
        outfile = joinpath(outdir, outfile)
        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        p = Tools.read_params(joinpath(outdir, "params.json"))
        d4 = Data.load_ascii_all(joinpath(outdir, "pulsar.debase.txt"))
        d1 = Data.clean(d4; threshold=0.001)
        Plot.single(d1, outdir; darkness=0.9, number=150, bin_st=p["bin_st"], bin_end=p["bin_end"], start=1, name_mod="pulsar", show_=true)
    end

    function process_psrdata(indir, outdir)
        p, debased_file, outdir = Data.process_psrdata(indir, outdir)
        Data.plot_psrdata(outdir, p)
        Data.process_data_andrzej(debased_file, outdir, p)
    end

    function process_psrdata_16(indir, outdir)
        p, debased_file, outdir = Data.process_psrdata_16(indir, outdir)
    end

    function process_psrdata_single(indir, outdir)
        Data.process_psrdata_single(indir, outdir)
    end


    """
    Phase-drift vs amplitude-modulation test on already-processed data.

    Reads pulsar.debase.txt and params.json from `outdir`, computes the
    coherent phase-slope statistic at f3 = 1/P3, compares against an
    amplitude-modulation null distribution, and saves a 3-panel PDF/PNG.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      phase_modulation(vpmout*"J1110-5637")
    """
    function phase_modulation(outdir; nreal=6000, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        p3_error = haskey(p, "p3_error") ? Float64(p["p3_error"]) : 0.0
        result = PhaseDrift.drift_test(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            p3_error=p3_error, nreal=nreal)
        println("Feature SNR:  $(round(result.snr, digits=1))")
        println("Slope:        $(round(result.slope, digits=4)) rad/bin  " *
                "($(round(rad2deg(result.slope), digits=2)) °/bin)")
        println("Significance: $(round(result.significance, digits=1)) σ")
        Plot.phase_drift(result, outdir, Int(p["nbin"]);
                         name_mod="pulsar", show_=show_)
        return result
    end


    """
    Phase-stability variant of `phase_modulation` — same data, same
    `PhaseDrift.drift_test`, same top two panels, but the bottom panel shows
    the local phase gradient dψ/dφ(φ) instead of the null-slope histogram,
    and a reduced χ² quantifies it.

    The coherent slope of `phase_modulation` collapses the whole profile into
    one number, so it cannot tell a systematic drift from a phase jump: pure
    amplitude modulation with a node produces a 180° step in ψ, which averages
    to a slope near zero and reads as "no drift". Here that case shows up as
    an isolated spike in dψ/dφ and a huge χ²_red.

      χ²_red ~ 1  – phase changes systematically (constant gradient fits)
      χ²_red ≫ 1  – phase jumps with longitude

    `snr_min` (default 3) sets which bins are trusted in that χ²; below it the
    Gaussian phase-error approximation σ_ψ ≈ σ_off/|L| breaks down. Writes
    `pulsar_phase_stability.pdf/.png`, so `phase_modulation` output is kept.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      phase_modulation2(vpmout*"J1110-5637")
    """
    function phase_modulation2(outdir; nreal=6000, snr_min=3.0, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        p3_error = haskey(p, "p3_error") ? Float64(p["p3_error"]) : 0.0
        result = PhaseDrift.drift_test(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            p3_error=p3_error, nreal=nreal, snr_min=snr_min)
        println("Feature SNR:  $(round(result.snr, digits=1))")
        println("Slope:        $(round(result.slope, digits=4)) rad/bin  " *
                "($(round(rad2deg(result.slope), digits=2)) °/bin)")
        println("Significance: $(round(result.significance, digits=1)) σ")
        if isnan(result.chi2_red)
            println("Stability:    too few high-S/N bins for χ² " *
                    "(only $(result.chi2_n) usable increments)")
        else
            println("Stability:    χ²_red = $(round(result.chi2_red, digits=2)) " *
                    "(dof $(result.chi2_dof)) — " *
                    (result.chi2_red < 2 ? "systematic drift" : "phase jumps with longitude"))
        end
        Plot.phase_stability(result, outdir, Int(p["nbin"]);
                             name_mod="pulsar", show_=show_)
        return result
    end


    """
    Windowed phase-stability test — `phase_modulation2` for pulsars whose
    global LRFS bin is empty.

    `phase_modulation2` asks the right question (is the phase gradient
    systematic, or is it a jump?) but reads it off the single global FFT bin,
    which P3 wobble empties: the f3 feature smears over Δk ≈ k·ΔP3/P3 bins
    with k = N/P3, so at short P3 in a long observation nothing is left. For
    J2053-7200 (P3 = 3.06, k ≈ 340, wobble ±1.3% ⇒ ~9 bins) it reports 0.5σ
    and a χ² from 8 usable increments, while `phase_modulation3` sees the same
    modulation at 80σ in short windows.

    This variant measures the same increments from the *windowed* LRFS:
    dψ/dφ(φ) = arg Σ_b conj(L_b[φ])·L_b[φ+1], summing the pairwise products
    over windows instead of over longitude (`phase_modulation3` does the
    latter and gets slope(t)). The products are invariant to each window's
    absolute phase, so the sum survives P3 wobble. Goodness of fit of a
    constant gradient is Monte-Carlo calibrated against flat-phase surrogates:

      χ²_red ~ 1  – constant gradient fits → genuine drift
      χ²_red ≫ 1  – phase jumps with longitude → not a drift

    Interpret it only when the modulation is actually detected — run
    `phase_modulation3` first and check its significance. The window sum is
    coherent, so for a drifter that reverses (J1750-3503) restrict the range
    to one episode with `pulse_st`/`pulse_end`, otherwise the episodes cancel.
    See `PhaseDrift.drift_test_profile` for the synthetic calibration,
    including the one regime that stays genuinely degenerate.

    Writes `pulsar_phase_stability_windowed.pdf/.png`, so the
    `phase_modulation2` output is kept.

    Typical call after phase_modulation3 has shown a significant drift:
      phase_modulation3(vpmout*"J2053-7200")
      phase_modulation2a(vpmout*"J2053-7200")
    """
    function phase_modulation2a(outdir; window=32, stride=1, nreal=500,
                                pulse_st=nothing, pulse_end=nothing, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        result = PhaseDrift.drift_test_profile(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            window=window, stride=stride, nreal=nreal,
            pulse_st=pulse_st, pulse_end=pulse_end)
        ptxt = result.p_value == 0 ? "p < $(round(1/nreal, sigdigits=1))" :
                                     "p = $(round(result.p_value, sigdigits=2))"
        println("Median window SNR: $(round(result.snr_med, digits=2))")
        println("Windows summed:    $(result.nwin) " *
                "(pulses $(result.pulse_range[1])-$(result.pulse_range[2]))")
        println("Gradient:          $(round(rad2deg(result.slope), digits=2)) °/bin " *
                "($(result.slope > 0 ? "positive" : "negative") drift sense)")
        println("Stability:         χ²_red = $(round(result.chi2_red, digits=2)) " *
                "($ptxt) — " *
                (result.chi2_red < 2 ? "systematic drift" : "phase jumps with longitude"))
        Plot.phase_stability_windowed(result, outdir, Int(p["nbin"]);
                                      name_mod="pulsar", show_=show_)
        return result
    end


    """
    Sliding-window, reversal-tolerant drift test — the detector that
    `phase_modulation`/`phase_modulation2` cannot be for drifters that switch
    drift direction (e.g. J1750-3503): their single coherent slope averages
    episodes of opposite sign to ~zero. Here the coherent slope statistic is
    evaluated in a short window slid pulse-by-pulse and the per-window drift
    quadratures |Im g_b| are combined incoherently, so episodes of + and −
    drift add instead of cancelling; the slope(t) panel shows the reversals
    directly. Significance comes from flat-phase surrogates whose noise is
    bootstrapped from the pulsar's own off-pulse region, shared across
    overlapping windows (see `PhaseDrift.drift_test_sliding` for the full
    construction and the caveats).

    `window` should be about half the shortest expected drift episode
    (J1750-3503: negative episodes 28±4 P → default 16). Scanning window
    ∈ {8,16,32,64,128} is a useful diagnostic, but quote the significance at
    the pre-chosen default or apply a trials correction. `stride=window`
    turns it into the disjoint-block variant.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1750-3503/.../", vpmout*"J1750-3503")
      phase_modulation3(vpmout*"J1750-3503")
    """
    function phase_modulation3(outdir; window=32, stride=1, nreal=1000, sig_min=3.0, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        result = PhaseDrift.drift_test_sliding(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            window=window, stride=stride, nreal=nreal, sig_min=sig_min)
        ptxt = result.p_value == 0 ? "p < $(round(1/nreal, sigdigits=1))" :
                                     "p = $(round(result.p_value, sigdigits=2))"
        npos = count(result.detected .& (result.slope .> 0))
        nneg = count(result.detected .& (result.slope .< 0))
        println("Median window SNR: $(round(result.snr_med, digits=2))")
        println("Drift detection:   $(round(result.significance, digits=1)) σ ($ptxt)")
        println("Windows ≥ $(sig_min)σ:    $(npos + nneg) of $(length(result.slope)) " *
                "($npos positive, $nneg negative drift)")
        if npos > 0 && nneg > 0
            println("                   both drift senses present — direction reverses")
        end
        Plot.phase_drift_sliding(result, outdir, Int(p["nbin"]);
                                 name_mod="pulsar", show_=show_)
        return result
    end


    """
    Run `phase_modulation3` on a list of pulsars (passed as a Vector of names or a file path).
    Automatically resolves output directory structure: `vpmout * name * "_16"`.
    """
    function phase_modulation3_list(vpmout, psr_input;
                                    window=32, stride=1, nreal=1000, sig_min=3.0, show_=false)
        names = String[]
        if psr_input isa AbstractString
            path = normpath(psr_input)
            if !isfile(path) && isfile(joinpath(@__DIR__, psr_input))
                path = normpath(joinpath(@__DIR__, psr_input))
            end
            isfile(path) || error("Pulsar list file not found: $psr_input")
            for line in eachline(path)
                s = strip(line)
                (isempty(s) || startswith(s, "#")) && continue
                push!(names, String(first(split(s))))
            end
        elseif psr_input isa AbstractVector
            for item in psr_input
                push!(names, String(item))
            end
        else
            error("Invalid psr_input: expected String file path or Vector of pulsar names")
        end

        results = Dict{String, Any}()

        for (i, name) in enumerate(names)
            outdir = vpmout * name * "_16"
            if !isdir(outdir)
                outdir = joinpath(vpmout, name * "_16")
            end
            if !isdir(outdir)
                outdir = vpmout * name
            end

            if !isdir(outdir)
                @warn "[$i/$(length(names))] Output dir not found for $name ($outdir), skipping."
                continue
            end

            println("\n[$i/$(length(names))] Running phase_modulation3 for PSR $name ($outdir)...")
            try
                res = phase_modulation3(outdir; window=window, stride=stride,
                                       nreal=nreal, sig_min=sig_min, show_=show_)
                results[name] = res
            catch e
                @warn "Failed for $name: $e"
            end
        end

        println("\nBatch phase_modulation3 complete: $(length(results)) / $(length(names)) processed.")
        return results
    end


    """
    Globally-optimized P3-fold (per-pulse Viterbi phase assignment), as an
    alternative to the `pfold -p3fold` refine used elsewhere in this file
    (e.g. `Data.process_psrdata_16` / `Data.p3fold_psrdata`). See
    p3fold-refine-notes.md for the rationale.

    Reads pulsar.debase.txt and params.json from `outdir` (same convention
    as `phase_modulation`), runs `P3FoldViterbi.fold`, plots the refined
    p3-fold, and reports per-pulse confidence/margin diagnostics.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      p3fold_refine(vpmout*"J1110-5637")
    """
    function p3fold_refine(outdir; ybins=nothing, n_iter=5, continuity_weight=0.2, darkness=1.0, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        yb   = isnothing(ybins) ? Int(p["p3_ybins"]) : ybins
        p3   = Float64(p["p3"])
        result = P3FoldViterbi.fold(
            data, p3, Int(p["bin_st"]), Int(p["bin_end"]);
            ybins=yb, n_iter=n_iter, continuity_weight=continuity_weight)
        println("Mean confidence: $(round(sum(result.confidence)/length(result.confidence), digits=3))")
        println("Mean margin:     $(round(sum(result.margin)/length(result.margin), digits=3))")
        folded_const = Tools.p3fold(data, p3, yb)
        Plot.p3fold_compare(result.folded, folded_const, result.p3_per_pulse, p3, outdir;
                            bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=darkness,
                            name_mod="pulsar_viterbi", show_=show_, repeat_num=4)
        return result
    end


    """
    Coherent, matched-filter P3-fold (`P3FoldViterbi.coherent_fold`) — for
    pulsars where the per-pulse modulation depth is too weak for any blind
    per-pulse matching (`p3fold_refine`/Viterbi, or `pfold`'s own block
    refine) to resolve, but where `phase_modulation`/`drift_test` still
    detects a highly significant *aggregate* coherent phase drift. See
    p3fold-refine-notes.md §3.3 and the conversation that motivated this.

    Reads pulsar.debase.txt and params.json from `outdir` (same convention
    as `phase_modulation`/`p3fold_refine`).

    Typical call after process_psrdata, once `phase_modulation` has
    confirmed a significant coherent drift:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      phase_modulation(vpmout*"J1110-5637")
      p3fold_coherent(vpmout*"J1110-5637")
    """
    function p3fold_coherent(outdir; ybins=nothing, lowpass_cutoff=1/300, filter_order=6, n_groups=4, darkness=1.0, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        yb   = isnothing(ybins) ? Int(p["p3_ybins"]) : ybins
        p3   = Float64(p["p3"])
        result = P3FoldViterbi.coherent_fold_jackknife(
            data, p3, Int(p["bin_st"]), Int(p["bin_end"]);
            ybins=yb, lowpass_cutoff=lowpass_cutoff, filter_order=filter_order, n_groups=n_groups)
        println("Matched-filter SNR: $(round(result.snr, digits=1))")
        folded_const = Tools.p3fold(data, p3, yb)
        intensity, _ = Tools.intensity_pulses(data[:, Int(p["bin_st"]):Int(p["bin_end"])])
        Plot.p3fold_compare(result.folded, folded_const, result.p3_per_pulse, p3, outdir;
                            bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=darkness,
                            name_mod="pulsar_coherent", show_=show_, repeat_num=4,
                            label="coherent fold", p3_per_pulse_err=result.p3_per_pulse_err,
                            intensity=intensity) # , p3_ylim=(-150, 150))
        return result
    end


    """
    Controlled comparison of `Tools.p3fold` (the "constant P3" naive fold
    used inside `p3fold_coherent`/`p3fold_refine`) against PSRSALSA's own
    `pfold -p3fold_norefine`, run on the *same* archive (`pulsar.debase.gg`)
    that `pulsar.debase.txt` was dumped from — unlike the low/high
    frequency-split norefine comparison in `Data.process_psrdata_16`
    ([data.jl:711-721]), which uses a different (half-bandwidth) archive and
    can't isolate whether a "looks much worse" gap is data or algorithm.

    Both sides see: the same pulses, the same zaps (PSRSALSA via its own
    `-w`, ours via `Data.zap!`), the same nominal P3/ybins/on-pulse window.
    If they still don't match after this, the difference is genuinely
    algorithmic (PSRSALSA's norefine fold vs `Tools.p3fold`'s modulo-P3
    sum), not a data/preprocessing mismatch.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      p3fold_norefine_compare(vpmout*"J1110-5637")
    """
    function p3fold_norefine_compare(outdir; ybins=nothing, darkness=1.0, show_=true)
        p  = Tools.read_params(joinpath(outdir, "params.json"))
        p3 = Float64(p["p3"])
        yb = isnothing(ybins) ? Int(p["p3_ybins"]) : ybins

        debase_gg = joinpath(outdir, "pulsar.debase.gg")
        run(pipeline(`pfold -p3fold_norefine -p3fold "$p3 $yb" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debase_gg`, stderr="errs.txt"))
        mv(replace(debase_gg, ".gg"=>".p3fold"), replace(debase_gg, ".gg"=>".p3fold_norefine"), force=true)
        folded_norefine = Data.load_ascii(replace(debase_gg, ".gg"=>".p3fold_norefine"))

        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        folded_const = Tools.p3fold(data, p3, yb)

        Plot.p3fold_compare(folded_norefine, folded_const, fill(p3, size(data, 1)), p3, outdir;
                            bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=darkness,
                            name_mod="pulsar_norefine_compare", show_=show_, repeat_num=4,
                            label="psrsalsa norefine")
        return (folded_norefine=folded_norefine, folded_const=folded_const)
    end


    """
    Automatically process all pulsars found in `dataroot`.
    For each pulsar directory, picks the lexicographically first observation
    subdirectory and runs process_psrdata_16 + analyse_p3folds_16_new.
    Errors for individual pulsars are caught so the loop continues.

    Arguments:
      dataroot  – parent directory containing JXXXX-XXXX subdirs
      vpmout    – output directory prefix (same convention as main())
      n_comp    – number of components passed to analyse_p3folds_16_new
    """
    function analyse_all(dataroot="/home/psr/data/new/", vpmout="/home/psr/output/")
        isdir(dataroot) || error("dataroot not found: $dataroot")

        p3_file = joinpath(@__DIR__, "input/drift_pulsars_P3.txt")
        p3_map = Dict{String,String}()
        if isfile(p3_file)
            for line in eachline(p3_file)
                parts = split(strip(line))
                length(parts) >= 2 && (p3_map[parts[1]] = parts[2])
            end
        end

        psr_dirs = sort(filter(d -> isdir(joinpath(dataroot, d)), readdir(dataroot)))
        isempty(psr_dirs) && (@warn "No pulsar directories found in $dataroot"; return)

        for psr in psr_dirs
            psr_path = joinpath(dataroot, psr)
            obs_dirs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
            if isempty(obs_dirs)
                @warn "No observation subdirectory for $psr, skipping"
                continue
            end
            indir  = joinpath(psr_path, obs_dirs[1]) * "/"
            outdir = vpmout * psr * "_16"
            if isdir(outdir)
                println("=== Skipping $psr (outdir exists) ===")
                continue
            end
            p3_str = get(p3_map, psr, "unknown")
            println("=== Processing $psr (obs: $(obs_dirs[1])), P3 = $p3_str ===")
            try
                # parse P3 value and error from "6.6(2)" notation and write to params.json
                m = match(r"^(\d+(?:\.\d+)?)\((\d+)\)$", p3_str)
                if !isnothing(m)
                    val_str = m.captures[1]
                    err_digits = parse(Int, m.captures[2])
                    p3_value = parse(Float64, val_str)
                    dot_pos = findfirst('.', val_str)
                    decimal_places = isnothing(dot_pos) ? 0 : length(val_str) - dot_pos
                    p3_error = err_digits * 10.0^(-decimal_places)
                    p3_ybins = Data.Functions.find_ybins(p3_value)
                    isdir(outdir) || mkdir(outdir)
                    params_file = joinpath(outdir, "params.json")
                    p = isfile(params_file) ? Tools.read_params(params_file) : Tools.default_params(params_file)
                    p["p3"] = p3_value
                    p["p3_error"] = p3_error
                    p["p3_ybins"] = p3_ybins
                    Tools.save_params(params_file, p)
                end
                process_psrdata_16(indir, outdir)
                print("n_comp for $psr [default=2]: ")
                n_comp_input = strip(readline())
                n_comp = isempty(n_comp_input) ? 2 : parse(Int, n_comp_input)
                Data.analyse_p3folds_16_new(outdir, "norefine"; n_comp=n_comp)
            catch e
                @warn "Failed for $psr: $e"
            end
            println("Analysis finished for PSR $psr")
            print("Continue? [Enter/y/yes = next, q = quit]: ")
            input = strip(readline())
            input == "q" && break
        end
    end












    function analyse_separations_todo(vpmout; csv_file=normpath(joinpath(@__DIR__, "..", "input", "separations_todo.csv")), type="norefine")
    if !isfile(csv_file)
        for alt in [normpath(joinpath(pwd(), "input", "separations_todo.csv")),
                    normpath(joinpath(@__DIR__, "input", "separations_todo.csv")),
                    "input/separations_todo.csv"]
            if isfile(alt)
                csv_file = alt
                break
            end
        end
    end
    isfile(csv_file) || error("separations_todo.csv not found: $csv_file")
    for (i, line) in enumerate(eachline(csv_file))
        i == 1 && continue  # header
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        fields = split(s, ',')
        length(fields) < 2 && continue
        
        name = String(strip(fields[1]))
        n_comp_val = tryparse(Int, strip(fields[2]))
        n_comp = isnothing(n_comp_val) ? 2 : n_comp_val
        outdir = vpmout * name * "_16"
        if !isdir(outdir)
            @warn "Output directory for $name not found: $outdir, skipping"
            continue
        end
        println("=== Running analyse_p3folds_16_new for $name (n_comp = $n_comp) ===")
        try
            Data.analyse_p3folds_16_new(outdir, type; n_comp=n_comp)
        catch e
            @warn "Failed for $name: $e"
        end
    end
end










"""
Phase-drift vs amplitude-modulation test on already-processed data.

    Reads pulsar.debase.txt and params.json from `outdir`, computes the
    coherent phase-slope statistic at f3 = 1/P3, compares against an
    amplitude-modulation null distribution, and saves a 3-panel PDF/PNG.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      phase_modulation(vpmout*"J1110-5637")
    """
    function phase_modulation(outdir; nreal=6000, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        p3_error = haskey(p, "p3_error") ? Float64(p["p3_error"]) : 0.0
        result = PhaseDrift.drift_test(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            p3_error=p3_error, nreal=nreal)
        println("Feature SNR:  $(round(result.snr, digits=1))")
        println("Slope:        $(round(result.slope, digits=4)) rad/bin  " *
                "($(round(rad2deg(result.slope), digits=2)) °/bin)")
        println("Significance: $(round(result.significance, digits=1)) σ")
        Plot.phase_drift(result, outdir, Int(p["nbin"]);
                         name_mod="pulsar", show_=show_)
        return result
    end


    """
    Phase-stability variant of `phase_modulation` — same data, same
    `PhaseDrift.drift_test`, same top two panels, but the bottom panel shows
    the local phase gradient dψ/dφ(φ) instead of the null-slope histogram,
    and a reduced χ² quantifies it.

    The coherent slope of `phase_modulation` collapses the whole profile into
    one number, so it cannot tell a systematic drift from a phase jump: pure
    amplitude modulation with a node produces a 180° step in ψ, which averages
    to a slope near zero and reads as "no drift". Here that case shows up as
    an isolated spike in dψ/dφ and a huge χ²_red.

      χ²_red ~ 1  – phase changes systematically (constant gradient fits)
      χ²_red ≫ 1  – phase jumps with longitude

    `snr_min` (default 3) sets which bins are trusted in that χ²; below it the
    Gaussian phase-error approximation σ_ψ ≈ σ_off/|L| breaks down. Writes
    `pulsar_phase_stability.pdf/.png`, so `phase_modulation` output is kept.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      phase_modulation2(vpmout*"J1110-5637")
    """
    function phase_modulation2(outdir; nreal=6000, snr_min=3.0, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        p3_error = haskey(p, "p3_error") ? Float64(p["p3_error"]) : 0.0
        result = PhaseDrift.drift_test(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            p3_error=p3_error, nreal=nreal, snr_min=snr_min)
        println("Feature SNR:  $(round(result.snr, digits=1))")
        println("Slope:        $(round(result.slope, digits=4)) rad/bin  " *
                "($(round(rad2deg(result.slope), digits=2)) °/bin)")
        println("Significance: $(round(result.significance, digits=1)) σ")
        if isnan(result.chi2_red)
            println("Stability:    too few high-S/N bins for χ² " *
                    "(only $(result.chi2_n) usable increments)")
        else
            println("Stability:    χ²_red = $(round(result.chi2_red, digits=2)) " *
                    "(dof $(result.chi2_dof)) — " *
                    (result.chi2_red < 2 ? "systematic drift" : "phase jumps with longitude"))
        end
        Plot.phase_stability(result, outdir, Int(p["nbin"]);
                             name_mod="pulsar", show_=show_)
        return result
    end


    """
    Windowed phase-stability test — `phase_modulation2` for pulsars whose
    global LRFS bin is empty.

    `phase_modulation2` asks the right question (is the phase gradient
    systematic, or is it a jump?) but reads it off the single global FFT bin,
    which P3 wobble empties: the f3 feature smears over Δk ≈ k·ΔP3/P3 bins
    with k = N/P3, so at short P3 in a long observation nothing is left. For
    J2053-7200 (P3 = 3.06, k ≈ 340, wobble ±1.3% ⇒ ~9 bins) it reports 0.5σ
    and a χ² from 8 usable increments, while `phase_modulation3` sees the same
    modulation at 80σ in short windows.

    This variant measures the same increments from the *windowed* LRFS:
    dψ/dφ(φ) = arg Σ_b conj(L_b[φ])·L_b[φ+1], summing the pairwise products
    over windows instead of over longitude (`phase_modulation3` does the
    latter and gets slope(t)). The products are invariant to each window's
    absolute phase, so the sum survives P3 wobble. Goodness of fit of a
    constant gradient is Monte-Carlo calibrated against flat-phase surrogates:

      χ²_red ~ 1  – constant gradient fits → genuine drift
      χ²_red ≫ 1  – phase jumps with longitude → not a drift

    Interpret it only when the modulation is actually detected — run
    `phase_modulation3` first and check its significance. The window sum is
    coherent, so for a drifter that reverses (J1750-3503) restrict the range
    to one episode with `pulse_st`/`pulse_end`, otherwise the episodes cancel.
    See `PhaseDrift.drift_test_profile` for the synthetic calibration,
    including the one regime that stays genuinely degenerate.

    Writes `pulsar_phase_stability_windowed.pdf/.png`, so the
    `phase_modulation2` output is kept.

    Typical call after phase_modulation3 has shown a significant drift:
      phase_modulation3(vpmout*"J2053-7200")
      phase_modulation2a(vpmout*"J2053-7200")
    """
    function phase_modulation2a(outdir; window=32, stride=1, nreal=500,
                                pulse_st=nothing, pulse_end=nothing, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        result = PhaseDrift.drift_test_profile(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            window=window, stride=stride, nreal=nreal,
            pulse_st=pulse_st, pulse_end=pulse_end)
        ptxt = result.p_value == 0 ? "p < $(round(1/nreal, sigdigits=1))" :
                                     "p = $(round(result.p_value, sigdigits=2))"
        println("Median window SNR: $(round(result.snr_med, digits=2))")
        println("Windows summed:    $(result.nwin) " *
                "(pulses $(result.pulse_range[1])-$(result.pulse_range[2]))")
        println("Gradient:          $(round(rad2deg(result.slope), digits=2)) °/bin " *
                "($(result.slope > 0 ? "positive" : "negative") drift sense)")
        println("Stability:         χ²_red = $(round(result.chi2_red, digits=2)) " *
                "($ptxt) — " *
                (result.chi2_red < 2 ? "systematic drift" : "phase jumps with longitude"))
        Plot.phase_stability_windowed(result, outdir, Int(p["nbin"]);
                                      name_mod="pulsar", show_=show_)
        return result
    end


    """
    Sliding-window, reversal-tolerant drift test — the detector that
    `phase_modulation`/`phase_modulation2` cannot be for drifters that switch
    drift direction (e.g. J1750-3503): their single coherent slope averages
    episodes of opposite sign to ~zero. Here the coherent slope statistic is
    evaluated in a short window slid pulse-by-pulse and the per-window drift
    quadratures |Im g_b| are combined incoherently, so episodes of + and −
    drift add instead of cancelling; the slope(t) panel shows the reversals
    directly. Significance comes from flat-phase surrogates whose noise is
    bootstrapped from the pulsar's own off-pulse region, shared across
    overlapping windows (see `PhaseDrift.drift_test_sliding` for the full
    construction and the caveats).

    `window` should be about half the shortest expected drift episode
    (J1750-3503: negative episodes 28±4 P → default 16). Scanning window
    ∈ {8,16,32,64,128} is a useful diagnostic, but quote the significance at
    the pre-chosen default or apply a trials correction. `stride=window`
    turns it into the disjoint-block variant.

    Typical call after process_psrdata:
      process_psrdata("/home/psr/data/new/J1750-3503/.../", vpmout*"J1750-3503")
      phase_modulation3(vpmout*"J1750-3503")
    """
    function phase_modulation3(outdir; window=32, stride=1, nreal=1000, sig_min=3.0, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=haskey(p, "zaps") ? p["zaps"] : nothing)
        result = PhaseDrift.drift_test_sliding(
            data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
            window=window, stride=stride, nreal=nreal, sig_min=sig_min)
        ptxt = result.p_value == 0 ? "p < $(round(1/nreal, sigdigits=1))" :
                                     "p = $(round(result.p_value, sigdigits=2))"
        npos = count(result.detected .& (result.slope .> 0))
        nneg = count(result.detected .& (result.slope .< 0))
        println("Median window SNR: $(round(result.snr_med, digits=2))")
        println("Drift detection:   $(round(result.significance, digits=1)) σ ($ptxt)")
        println("Windows ≥ $(sig_min)σ:    $(npos + nneg) of $(length(result.slope)) " *
                "($npos positive, $nneg negative drift)")
        if npos > 0 && nneg > 0
            println("                   both drift senses present — direction reverses")
        end
        Plot.phase_drift_sliding(result, outdir, Int(p["nbin"]);
                                 name_mod="pulsar", show_=show_)
        return result
    end























"""
Run `phase_modulation3` on a list of pulsars (passed as a Vector of names or a file path).

function phase_modulation3_list(vpmout::String, psrs::Union{Vector{String}, String};
                                window=32, stride=1, nreal=1000, sig_min=3.0, show_=false)
    
    names = if psrs isa String
        isfile(psrs) || error("Pulsar list file not found: $psrs")
        [String(first(split(strip(line)))) for line in eachline(psrs)
         if !isempty(strip(line)) && !startswith(strip(line), "#")]
    else
        psrs
    end

    results = Dict{String, Any}()

    for (i, name) in enumerate(names)
        outdir = vpmout * name
        if !isdir(outdir)
            @warn "[$i/$(length(names))] Output dir not found for $name ($outdir), skipping."
            continue
        end

        println("[$i/$(length(names))] Processing PSR $name...")
        try
            res = phase_modulation3(outdir; window=window, stride=stride,
                                   nreal=nreal, sig_min=sig_min, show_=show_)
            results[name] = res
        catch e
            @warn "Failed for $name: $e"
        end
    end

    return results
end

"""












    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"

        # PSR J0034-0721
        #process_psrdata_16("/home/psr/data/new/J0034-0721/2019-10-18-22:29:51/", vpmout*"J0034-0721_16")
        #Data.analyse_p3folds_16(vpmout*"J0034-0721_16", "norefine")
        #process_psrdata(vpmout*"J0034-0721", vpmout*"J0034-0721") # P. nice
        #Data.analyse_p3folds_16_new(vpmout*"J0034-0721_16", "norefine"; n_comp=2)
        #Data.analyse_p3folds_16_new(vpmout*"J0034-0721_16", "refine"; n_comp=2)
        #Data.position_angle(vpmout*"J0034-0721_16")
        #Data.geometry_analysis(vpmout*"J0034-0721_16")
        #test(vpmout)
        #test2(vpmout)
        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)

        #Tools.clean_all(vpmout)
        #analyse_all()

        #Data.analyse_p3folds_16_new(vpmout*"J1539-6322_16", "norefine", n_comp=2)

        #Data.analyse_separations_todo(vpmout; csv_file=joinpath(@__DIR__, "..", "input", "separations_todo.csv"), type="norefine")















        # P-Pdot diagram based on the ATNF catalogue (input/psrcat.db)
        #Plot.ppdot("output")
        # the same diagram with the component offsets (input/offsets.csv)
        #Plot.ppdot_offsets("output")
        
        #Plot.ppdot_w50("output"; mode=:w50, name_mod="w50_ms")

        # Parameter relations distinguishing Drifting vs P3-only pulsars
        #Relations.plot_all_relations("output")

        phase_modulation3_list("output", "input/p3only_pulsars_P3.txt")



    end

end # module

SpaTs.main()

println("Bye")