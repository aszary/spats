module SpaTs
    using ArgParse
    using Glob
    using JSON
    using Statistics

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")
    include("modules/phase_modulation.jl")
    include("modules/p3fold_viterbi.jl")
    include("modules/offset_analysis.jl")


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

    function process_psrdata(indir, outdir; zap_mode=:pazi, zap_nsigma=5.0)
        p, debased_file, outdir = Data.process_psrdata(indir, outdir; zap_mode=zap_mode, zap_nsigma=zap_nsigma)
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
        Data.zap!(data; ranges=Data.zaps_of(p))
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
        Data.zap!(data; ranges=Data.zaps_of(p))
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
        Data.zap!(data; ranges=Data.zaps_of(p))
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
        Data.zap!(data; ranges=Data.zaps_of(p))
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
        Data.zap!(data; ranges=Data.zaps_of(p))
        p3   = Float64(p["p3"])
        # same as p3fold_coherent: p["p3_ybins"] can be stale (J2053 aliasing)
        yb   = isnothing(ybins) ? Data.Functions.find_ybins(p3, size(data, 1)) : ybins
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
    function p3fold_coherent(outdir; ybins=nothing, lowpass_cutoff=1/300, filter_order=6, n_groups=4, p3_window=60, bin_search=2, darkness=1.0, show_=true, name_mod="pulsar_coherent")
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=Data.zaps_of(p))
        p3   = Float64(p["p3"])
        # p["p3_ybins"] bywa nieaktualne/za duze wzgledem faktycznej liczby
        # impulsow na cykl P3 (np. P3≈3 z ybins=10 daje aliasing) — licz
        # wlasciwa wartosc z find_ybins zamiast ufac bezkrytycznie plikowi
        yb   = isnothing(ybins) ? Data.Functions.find_ybins(p3, size(data, 1)) : ybins
        result = P3FoldViterbi.coherent_fold_jackknife(
            data, p3, Int(p["bin_st"]), Int(p["bin_end"]);
            ybins=yb, lowpass_cutoff=lowpass_cutoff, filter_order=filter_order, n_groups=n_groups,
            p3_window=p3_window, bin_search=bin_search)
        println("Matched-filter SNR: $(round(result.snr, digits=1))")
        folded_const = Tools.p3fold(data, p3, yb)
        intensity, _ = Tools.intensity_pulses(data[:, Int(p["bin_st"]):Int(p["bin_end"])])
        Plot.p3fold_compare(result.folded, folded_const, result.p3_per_pulse, p3, outdir;
                            bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=darkness,
                            name_mod=name_mod, show_=show_, repeat_num=4,
                            label="coherent fold", p3_per_pulse_err=result.p3_per_pulse_err,
                            intensity=intensity) # , p3_ylim=(-150, 150))
        # diagnostic for WHY a given peak/dip on P3(t) should be trusted:
        # overlay the n_groups independent longitude-subband tracks against
        # the combined result — consistent across most subbands = real wobble
        Plot.subband_consistency(result.group_p3, result.p3_per_pulse, outdir;
                                 p3_nominal=p3, name_mod=name_mod, show_=show_)
        return result
    end


    """
    `p3fold_coherent` with a curvature-penalized (Whittaker/Hodrick-Prescott
    style) smoother instead of a fixed Butterworth lowpass — no
    `lowpass_cutoff`/`filter_order` to pick or scan for at all, and it
    penalises how much the P3(t) track is allowed to *bend* rather than
    which Fourier frequencies survive, matching the physical prior that P3
    changes gradually (no jumps). The smoothing strength is itself picked
    automatically per pulsar (discrepancy principle against the noise level
    estimated from the baseband's own spectrum). See
    `P3FoldViterbi.whittaker_smooth_auto` for the derivation.

    If the automatic smoothing feels too aggressive (or too loose), pass
    `auto_strength` — `<1` allows more bend (less smoothing) than the
    automatic noise estimate calls for, `>1` forces more.

    Typical call — no cutoff scan needed first:
      process_psrdata("/home/psr/data/new/J1110-5637/.../", vpmout*"J1110-5637")
      p3fold_coherent_auto(vpmout*"J1110-5637")
      p3fold_coherent_auto(vpmout*"J1110-5637", auto_strength=0.5)  # less smoothing
    """
    function p3fold_coherent_auto(outdir; ybins=nothing, n_groups=4, p3_window=60,
                                  auto_strength=1.0, darkness=1.0, show_=true,
                                  name_mod="pulsar_coherent_auto")
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=Data.zaps_of(p))
        p3   = Float64(p["p3"])
        yb   = isnothing(ybins) ? Data.Functions.find_ybins(p3, size(data, 1)) : ybins
        result = P3FoldViterbi.coherent_fold_jackknife(
            data, p3, Int(p["bin_st"]), Int(p["bin_end"]);
            ybins=yb, n_groups=n_groups, p3_window=p3_window, auto=true,
            auto_strength=auto_strength)
        println("Matched-filter SNR: $(round(result.snr, digits=1))")
        folded_const = Tools.p3fold(data, p3, yb)
        intensity, _ = Tools.intensity_pulses(data[:, Int(p["bin_st"]):Int(p["bin_end"])])
        Plot.p3fold_compare(result.folded, folded_const, result.p3_per_pulse, p3, outdir;
                            bin_st=p["bin_st"], bin_end=p["bin_end"], darkness=darkness,
                            name_mod=name_mod, show_=show_, repeat_num=4,
                            label="coherent fold (auto)", p3_per_pulse_err=result.p3_per_pulse_err,
                            intensity=intensity)
        return result
    end


    """
    Diagnostic run before `p3fold_coherent` — scans several `lowpass_cutoff`
    values and plots how track "trustworthiness" (subband consistency) and
    noise level trade off, so the cutoff for `p3fold_coherent` is picked from
    numbers instead of by re-running the fold with different values by hand.
    See `P3FoldViterbi.scan_lowpass_cutoff` for the metric definitions.

    Typical call, before settling on a cutoff for a given pulsar:
      process_psrdata("/home/psr/data/new/J1750-3503/.../", vpmout*"J1750-3503")
      p3fold_cutoff_scan(vpmout*"J1750-3503")
      # inspect the plot, then:
      p3fold_coherent(vpmout*"J1750-3503", lowpass_cutoff=<chosen value>)
    """
    function p3fold_cutoff_scan(outdir; cutoffs=nothing, min_reliability=0.5,
                                filter_order=6, n_groups=8, ybins=nothing, chosen=nothing, show_=true)
        p    = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=Data.zaps_of(p))
        p3 = Float64(p["p3"])
        N  = size(data, 1)
        yb = isnothing(ybins) ? Data.Functions.find_ybins(p3, N) : ybins
        # Grid now scales with P3 instead of one fixed range for every
        # pulsar: a cutoff of X cycles/pulse corresponds to a timescale of
        # 1/X pulses, so bound that timescale physically —
        #   fastest meaningful: ~2 P3 cycles (can't define a slope/wobble
        #     faster than that; least 1 cycle just to have samples)
        #   slowest meaningful: ~N/3 pulses (anything slower than a third
        #     of the whole observation isn't "variability", it's a constant
        #     offset) — capped so very short observations still get a
        #     sensible spread
        if isnothing(cutoffs)
            lo_timescale = max(3.0, 2 * p3)
            hi_timescale = max(lo_timescale * 5, N / 3)
            cutoffs = 1.0 ./ (10 .^ range(log10(hi_timescale), log10(lo_timescale), length=20))
        end
        results = P3FoldViterbi.scan_lowpass_cutoff(
            data, p3, Int(p["bin_st"]), Int(p["bin_end"]);
            cutoffs=cutoffs, filter_order=filter_order, n_groups=n_groups, ybins=yb)
        for r in results
            if hasproperty(r, :valid) && !r.valid
                println("cutoff=$(round(r.cutoff, sigdigits=3))  ODRZUCONY " *
                        "(P3=2π/nachylenie eksplodowało gdzieś na tym cutoffie — " *
                        "cały punkt pominięty, nie tylko wygładzony)")
            else
                println("cutoff=$(round(r.cutoff, sigdigits=3))  " *
                        "consistency=$(round(r.consistency, digits=3))  " *
                        "p3_std=$(round(r.p3_std, digits=3))  " *
                        "signal_std=$(round(r.signal_std, digits=3))  " *
                        "R=$(round(r.reliability, digits=2))" *
                        (r.reliability >= min_reliability ? "" : " (szum>sygnał)"))
            end
        end
        # snr is identical across all rows by design (computed pre-filter, so
        # lowpass_cutoff can't affect it) — it's a sanity reference, not part
        # of the scan itself. Pure noise gives ~1.41 for this statistic
        # (different from PhaseDrift.drift_test's own ~1.25 floor — the two
        # use different magnitude statistics); snr < 3 means the f3 feature
        # itself is marginal, independent of anything the cutoff scan found.
        if !isempty(results)
            snr0 = results[1].snr
            floor_ratio = snr0 / sqrt(2)
            println("snr=$(round(snr0, digits=1)) jest takie samo dla każdego cutoffu " *
                    "(liczone przed filtrowaniem) — czysty szum daje tu ~1.41, więc to " *
                    "$(round(floor_ratio, digits=1))× powyżej szumu" *
                    (snr0 < 3 ? " — SŁABY sygnał, traktuj fold ostrożnie." : "."))
        end
        # The grid can have several comparable local peaks in signal_std (seen
        # in practice: two nearly-tied peaks a decade apart), so refining only
        # the single global argmax risks locking onto whichever one happened
        # to be a hair higher due to grid noise. Instead: find every local
        # maximum on the grid, refine *each* one continuously (golden-section,
        # within its own immediate neighbours), then take whichever refined
        # result actually scores best — a real comparison after polishing,
        # not a bet on the raw grid alone.
        #
        # Selection criterion: signal_std grows with cutoff by construction
        # (a wider filter always lets more signal AND noise through), so
        # maximising it alone always drifts to the top edge of the grid
        # whenever consistency is near noise level (seen: J1110, J2053 picked
        # the noisiest cutoff, J2053 even beyond the grid). So only cutoffs
        # whose combined track is still majority signal (reliability R >=
        # min_reliability, see cutoff_metrics) compete on signal_std. If no
        # cutoff reaches that, P3(t) variability isn't reliably detected at
        # all — then take the most reproducible track (max consistency).
        # Refinement brackets stay inside the scanned grid.
        valid = filter(r -> isfinite(r.consistency) && isfinite(r.signal_std), results)
        suggested = nothing
        refined = nothing
        refine_trace = Tuple{Float64,Float64}[]
        if !isempty(valid)
            sorted = sort(valid, by = r -> r.cutoff)
            sorted_cutoffs = [r.cutoff for r in sorted]
            n = length(sorted)
            use_signal = any(r -> r.reliability >= min_reliability, sorted)
            objective = use_signal ? :signal_std : :consistency
            key(r) = use_signal ? (r.reliability >= min_reliability ? r.signal_std : -Inf) :
                                  r.consistency
            if !use_signal
                println("UWAGA: żaden cutoff nie ma R >= $min_reliability — zmienność P3(t) " *
                        "nie jest wiarygodnie wykryta (szum dominuje). Wybieram cutoff " *
                        "o najwyższej spójności; traktuj P3(t) jako ~stałe.")
            end
            keys_ = [key(r) for r in sorted]
            is_local_max(i) = isfinite(keys_[i]) &&
                              (i == 1 || keys_[i] >= keys_[i-1]) &&
                              (i == n || keys_[i] >= keys_[i+1])
            candidate_idxs = [i for i in 1:n if is_local_max(i)]
            isempty(candidate_idxs) && (candidate_idxs = [argmax(keys_)])

            best_key = -Inf
            for idx in candidate_idxs
                # grid point itself is a candidate (golden-section never
                # evaluates the bracket ends)
                if keys_[idx] > best_key
                    best_key, refined = keys_[idx], merge(sorted[idx], (trace=Tuple{Float64,Float64}[],))
                end
                n == 1 && continue
                lo = sorted_cutoffs[max(1, idx-1)]
                hi = sorted_cutoffs[min(n, idx+1)]
                r = P3FoldViterbi.refine_cutoff(
                    data, p3, Int(p["bin_st"]), Int(p["bin_end"]), lo, hi;
                    filter_order=filter_order, n_groups=n_groups, ybins=yb,
                    objective=objective, min_reliability=min_reliability)
                append!(refine_trace, r.trace)
                if r.valid && key(r) > best_key
                    best_key, refined = key(r), r
                end
            end
            suggested = refined.cutoff
            println("Wybrany cutoff=$(round(suggested, sigdigits=3)) " *
                    "(skala ~$(round(1/suggested, digits=0)) impulsów)  " *
                    "consistency=$(round(refined.consistency, digits=3))  " *
                    "R=$(round(refined.reliability, digits=2))  " *
                    "kryterium=$(use_signal ? "max signal_std przy R>=$min_reliability" : "max consistency (fallback)")")
        end
        # only valid (non-blown-up) cutoffs get drawn at all — a rejected
        # cutoff doesn't even get a point on the scan plots, per the whole
        # candidate being thrown out rather than patched
        Plot.lowpass_cutoff_scan(valid, outdir; name_mod="pulsar",
                                 chosen=something(chosen, suggested), chosen_metrics=refined,
                                 refine_trace=refine_trace,
                                 # r at which R = min_reliability (inverse Spearman-Brown)
                                 r_threshold=min_reliability / (n_groups - (n_groups - 1) * min_reliability),
                                 show_=show_)
        # drop cutoffs whose P3(t) track is basically identical to a
        # neighbour's before plotting the track overlay — several grid
        # points (typically the tightest cutoffs, all collapsing to ≈flat
        # lines at nominal P3) add nothing and just clutter the legend
        tracks_for_plot = P3FoldViterbi.dedupe_redundant_cutoffs(valid)
        Plot.p3_tracks_vs_cutoff(tracks_for_plot, outdir; name_mod="pulsar",
                                 highlight_cutoff=suggested, p3_nominal=p3, show_=show_)
        # diagnostic for WHY results differ across cutoffs: a 2D map (cutoff
        # x pulse number) — a feature present across many rows is real, one
        # confined to a single row is a filter artifact specific to that cutoff
        Plot.p3_cutoff_heatmap(valid, outdir; p3_nominal=p3, name_mod="pulsar", show_=show_)
        return (results=results, suggested_cutoff=suggested)
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
        # same rule as p3fold_coherent (p["p3_ybins"] can be stale); used for
        # BOTH sides of the comparison, so they stay identical
        yb = isnothing(ybins) ? Data.Functions.find_ybins(p3) : ybins

        debase_gg = joinpath(outdir, "pulsar.debase.gg")
        run(pipeline(`pfold -p3fold_norefine -p3fold "$p3 $yb" -onpulse "$(p["bin_st"]) $(p["bin_end"])" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debase_gg`, stderr="errs.txt"))
        mv(replace(debase_gg, ".gg"=>".p3fold"), replace(debase_gg, ".gg"=>".p3fold_norefine"), force=true)
        folded_norefine = Data.load_ascii(replace(debase_gg, ".gg"=>".p3fold_norefine"))

        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        Data.zap!(data; ranges=Data.zaps_of(p))
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


    """
    Minimal preparation of `vpmout/<pulsar>/pulsar.debase.txt` from
    `dataroot/<pulsar>/<first obs>/`: psradd + debase only (no P3 search,
    no 2DFS/pfold, no zapping). Without params.json, pmod asks for the
    on-pulse window interactively, as in normal processing.
    """
    function prepare_debase(psr; vpmout="/home/psr/output/", dataroot="/home/psr/data/new/")
        psr_path = joinpath(dataroot, psr)
        obs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
        isempty(obs) && error("brak katalogu obserwacji w $psr_path")
        indir, outdir = joinpath(psr_path, obs[1]) * "/", vpmout * psr
        mkpath(outdir)
        params_file = joinpath(outdir, "params.json")
        p = isfile(params_file) ? Tools.read_params(params_file) : Tools.default_params(params_file)
        outfile = joinpath(outdir, "pulsar.spCF")
        println("Przygotowuję $psr (obs: $(obs[1])): psradd + debase")
        Data.add_psrfiles(indir, outfile)
        isnothing(get(p, "nsubint", nothing)) && Data.get_nsubint(outfile, params_file, p)
        Data.debase(outfile, params_file, p)
    end


    """
    Tryb OZNACZANIA złych impulsów — do opracowania/dostrojenia detektora na
    podstawie Twoich decyzji. Nic nie zapuje i nie zmienia params.json.

    Otwiera `pazi pulsar.spCF` — zapujesz tak jak zwykle i wychodzisz z pazi.
    Komenda `paz -w "..."`, którą pazi wypisuje, jest odczytywana
    automatycznie. Zapisuje `input/zap_labels/<pulsar>.csv`: ~16 cech każdego
    impulsu + `label` (Twoje z pazi) + `auto` (obecny detektor), oraz
    `input/zap_labels/<pulsar>.paz.txt` z oryginalną komendą (impulsy i
    kanały). Na koniec porównanie z automatem.

      label_bad_pulses("J0034-0721")

    Archiwum jest zawsze budowane od nowa z surowych danych (`prepare_debase`:
    psradd + debase), żeby wcześniejsze zapy nie były już w nim ukryte.
    """
    function label_bad_pulses(psr; vpmout="/home/psr/output/", dataroot="/home/psr/data/new/")
        outdir = vpmout * psr
        archive = joinpath(outdir, "pulsar.spCF")
        # always rebuild from raw data: earlier runs may have zapped
        # pulsar.spCF in place (paz -m), which would hide those pulses in pazi
        prepare_debase(psr; vpmout=vpmout, dataroot=dataroot)
        p = Tools.read_params(joinpath(outdir, "params.json"))
        data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
        N = size(data, 1)
        bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
        println("\n=== OZNACZANIE: $psr ($N impulsów, on-pulse $bst-$bend) ===")
        z = Data.run_pazi(archive)
        rows = [s + 1 for s in z.subints if 0 <= s < N]
        Data.save_zap_labels(psr, data, bst, bend, rows)
        cmdfile = joinpath(@__DIR__, "input", "zap_labels", "$psr.paz.txt")
        open(f -> println(f, something(z.cmd, "")), cmdfile, "w")
        println("Twoje zaps: \"$(Data.pulses_to_zaps(rows))\"" *
                (isempty(z.chans) ? "" : "   kanały: $(join(z.chans, " "))"))
        return z
    end


    """
    Strojenie automatycznego zapowania (`Data.auto_zap_core`) na Twoich
    etykietach z pazi (`input/zap_labels/<pulsar>.csv`).

    Dla każdego pulsara z etykietami i danymi (`vpmout/<psr>/pulsar.debase.txt`
    — po oznaczaniu to są dane BEZ zapów) sprawdza 729 kombinacji parametrów
    i liczy zgodność z Tobą (F1: 1 = identycznie, 0 = nic wspólnego) oraz
    S/N: bez zapów / Twoje / automat. Wypisuje najlepsze ustawienia,
    zapisuje pełną tabelę `input/zap_labels/tuning.csv` i dla najlepszych
    parametrów wykres `vpmout/<psr>/pulsar_zap_compare.png` (czerwone = Ty,
    niebieskie = automat).

      tune_autozap(["J0034-0721", "J0459-0210", "J0624-0424", ...])
    """
    function tune_autozap(psrs; vpmout="/home/psr/output/")
        sets = []
        for psr in psrs
            rows = Data.read_zap_labels(psr)
            outdir = vpmout * psr
            tf, pf = joinpath(outdir, "pulsar.debase.txt"), joinpath(outdir, "params.json")
            if isnothing(rows) || !isfile(tf) || !isfile(pf)
                println("SKIP $psr — brak etykiet (input/zap_labels/$psr.csv) albo danych")
                continue
            end
            p = Tools.read_params(pf)
            data = Data.load_ascii(tf)
            bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
            sn = Data.pulse_signal_noise(data, bst, bend)
            push!(sets, (psr=psr, outdir=outdir, data=data, bst=bst, bend=bend, sn=sn,
                         nz=Data.noise_jump_z(sn.vwhite, sn.good), rows=rows, N=size(data, 1)))
            println("$psr: $(size(data, 1)) impulsów, Twoje zapy: $(length(rows))")
        end
        isempty(sets) && (println("Brak danych do strojenia."); return nothing)

        snr(sn, rows) = (k = copy(sn.good); k[filter(i -> 1 <= i <= length(k), rows)] .= false;
                         sv = sum(sn.v[k]); sv > 0 ? sum(sn.s[k]) / sqrt(sv) : NaN)
        grid = Iterators.product((5.0, 8.0, 12.0), (2.0, 3.0, 4.0), (0, 1, 2),
                                 (3, 5, 8), (5.0, 6.0, 8.0), (0.01, 0.02, 0.05))
        results = []
        for (target, zsig, pad, min_run, nsigma, min_gain) in grid
            kw = (target=target, zsig=zsig, pad=pad, min_run=min_run, nsigma=nsigma, min_gain=min_gain)
            f1s, ratios = Float64[], Float64[]
            for st in sets
                a = Data.auto_zap_core(st.sn; noise_z=st.nz, kw...)
                push!(f1s, Data.zap_agreement(a.bad, st.rows, st.N).f1)
                push!(ratios, a.snr_after / snr(st.sn, st.rows))
            end
            push!(results, (kw=kw, f1=mean(f1s), minf1=minimum(f1s), snr_ratio=mean(ratios), f1s=f1s))
        end
        sort!(results, by=r -> (-r.f1, -r.minf1))

        names_ = [st.psr for st in sets]
        csv = joinpath(@__DIR__, "input", "zap_labels", "tuning.csv")
        open(csv, "w") do f
            println(f, join(vcat(["target", "zsig", "pad", "min_run", "nsigma", "min_gain",
                                  "mean_f1", "min_f1", "snr_auto_over_user"], names_), ","))
            for r in results
                println(f, join(vcat(collect(values(r.kw)), round(r.f1, digits=4), round(r.minf1, digits=4),
                                     round(r.snr_ratio, digits=4), round.(r.f1s, digits=3)), ","))
            end
        end
        println("\nPełna tabela: $csv\n\nNajlepsze 10 (średnie F1 po pulsarach):")
        for r in results[1:min(10, end)]
            println("  F1=$(round(r.f1, digits=3)) (min $(round(r.minf1, digits=2)))  " *
                    "S/N auto/Twoje=$(round(r.snr_ratio, digits=3))  $(r.kw)")
        end

        best = results[1].kw
        default = (target=8.0, zsig=3.0, pad=1, min_run=5, nsigma=6.0, min_gain=0.02)
        for (label, kw) in (("DOMYŚLNE", default), ("NAJLEPSZE", best))
            println("\n$label $kw")
            for st in sets
                a = Data.auto_zap_core(st.sn; noise_z=st.nz, kw...)
                m = Data.zap_agreement(a.bad, st.rows, st.N)
                println("  $(rpad(st.psr, 12)) Ty=$(lpad(length(st.rows), 4)) auto=$(lpad(length(a.bad), 4))  " *
                        "F1=$(round(m.f1, digits=2)) prec=$(round(m.precision, digits=2)) " *
                        "rec=$(round(m.recall, digits=2))  S/N bez=$(round(a.snr_before, digits=1)) " *
                        "Ty=$(round(snr(st.sn, st.rows), digits=1)) auto=$(round(a.snr_after, digits=1))  okno=$(a.w)")
                if label == "NAJLEPSZE"
                    Plot.zap_compare(st.data, st.bst, st.bend, a, st.rows, st.outdir;
                                     title_="$(st.psr): F1=$(round(m.f1, digits=2))  " *
                                            "czerwone = Ty ($(length(st.rows))), niebieskie = automat ($(length(a.bad)))")
                end
            end
        end
        return results
    end


    """
    Ślepy test automatu: przerabia pulsary od nowa, automat (DOMYŚLNE
    parametry, bez wiedzy o Twoich zapach) wycina, a dopiero potem porównanie
    z Twoimi etykietami (input/zap_labels/<psr>.csv — nie są zmieniane).

    Dla każdego pulsara: dane od zera (psradd + debase) → automat → wydruk
    (zgodność, S/N bez / Twoje / automat) → `pulsar_zap_compare.png` →
    pazi NA DANYCH PO AUTOMACIE (dozapuj, jeśli coś zostało — zapisane jako
    przeoczone) → Twoja ocena. Wyniki: input/zap_labels/blind_test.csv.

      blind_test_autozap(["J0034-0721", "J0459-0210", ...])
    """
    function blind_test_autozap(psrs; vpmout="/home/psr/output/", dataroot="/home/psr/data/new/",
                                interactive=true, show_plots=!interactive)
        out = joinpath(@__DIR__, "input", "zap_labels", "blind_test.csv")
        mkpath(dirname(out))
        summary = String[]
        isfile(out) || open(f -> println(f, "pulsar,n_pulses,n_auto,n_user,f1,precision,recall," *
                                           "snr_none,snr_user,snr_auto,n_missed_in_pazi,verdict,auto_zaps"),
                            out, "w")
        verdicts = Dict("1" => "lepiej_niz_ja", "2" => "tak_samo", "3" => "gorzej")
        f(x) = isnan(x) ? "" : string(round(x, digits=1))
        for psr in psrs
            outdir = vpmout * psr
            println("\n===== ŚLEPY TEST: $psr =====")
            try
                prepare_debase(psr; vpmout=vpmout, dataroot=dataroot)
                p = Tools.read_params(joinpath(outdir, "params.json"))
                data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
                N = size(data, 1)
                bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
                sn = Data.pulse_signal_noise(data, bst, bend)
                a = Data.auto_zap_core(sn)
                user = Data.read_zap_labels(psr)
                user_s = isnothing(user) ? "brak etykiet" : "$(length(user))"
                m = isnothing(user) ? nothing : Data.zap_agreement(a.bad, user, N)
                println("Automat: $(length(a.bad)) z $N  →  \"$(Data.pulses_to_zaps(a.bad))\"")
                println("Ty: $user_s" * (isnothing(user) ? "" : "  →  \"$(Data.pulses_to_zaps(user))\""))
                isnothing(m) || println("Zgodność: F1=$(round(m.f1, digits=2))  wspólne=$(m.tp)  " *
                                        "tylko automat=$(m.fp)  tylko Ty=$(m.fn)")
                snr_u = isnothing(user) ? NaN : Data.zap_snr(sn, user)
                # standard S/N of the integrated profile (same window), before and
                # after the automat's zaps
                sp0 = Data.profile_snr(data, bst, bend)
                sp1 = Data.profile_snr(data, bst, bend; drop=a.bad)
                println("S/N profilu: bez zapów=$(f(sp0))  po automacie=$(f(sp1))  " *
                        "($(isnan(sp0) || isnan(sp1) ? "?" : string(round(100 * (sp1 / sp0 - 1), digits=1)) * "%"))")
                println("   (miara automatu, liczy też falowanie tła: $(f(a.snr_before)) → $(f(a.snr_after))" *
                        (isnothing(user) ? "" : ", Twoje=$(f(snr_u))") * ")")
                println("Decyzja: $(a.why)   [okno=$(a.w), podejrzane pojedyncze=$(a.n_single), " *
                        "progi z odcinkami=$(a.n_cand), okna off: $(sn.nwin) użyte / $(sn.nleak) z emisją pulsara]")
                kinds = ["pusty", "pojedynczy śmieć", "null/brak sygnału", "zaszumiony odcinek", "margines"]
                println("Powody: " * join(["$(kinds[k])=$(count(==(k), a.reason))" for k in 1:5], "  "))
                push!(summary, "$(rpad(psr, 12)) wycięte=$(lpad(length(a.bad), 4)) z $(lpad(N, 4))  " *
                               "S/N profilu: $(f(sp0)) → $(f(sp1))   (miara automatu: " *
                               "$(f(a.snr_before)) → $(f(a.snr_after)))")
                # main view: original (left) vs what the automat cuts and why (right)
                Plot.autozap_view(data, bst, bend, a, outdir; show_=show_plots, title_="$psr — automat",
                                  snr_pazi=(sp0, sp1))
                # side file: your zaps vs automat (saved only)
                Plot.zap_compare(data, bst, bend, a, something(user, Int[]), outdir; show_=false,
                                 title_="$psr: czerwone = Ty ($user_s), niebieskie = automat ($(length(a.bad)))")

                missed, v = Int[], ""
                if interactive
                    # look at the automatic result in pazi
                    archive = joinpath(outdir, "pulsar.spCF")
                    preview = joinpath(outdir, "pulsar.autozap.ar")
                    cp(archive, preview; force=true)
                    Data.apply_zaps(preview, Data.pulses_to_zaps(a.bad))
                    println("pazi pokazuje dane PO automacie. Dozapuj, jeśli coś zostało (zapisane jako " *
                            "przeoczone), i wyjdź.")
                    z = Data.run_pazi(preview)
                    missed = [s + 1 for s in z.subints if 0 <= s < N && !(s + 1 in a.bad)]
                    rm(preview, force=true)
                    print("Ocena automatu: 1 = lepiej niż ja, 2 = tak samo dobrze, 3 = gorzej, Enter = bez oceny: ")
                    v = get(verdicts, strip(readline()), "")
                end
                open(out, "a") do io
                    println(io, join([psr, N, length(a.bad), isnothing(user) ? "" : length(user),
                                      isnothing(m) ? "" : round(m.f1, digits=3),
                                      isnothing(m) ? "" : round(m.precision, digits=3),
                                      isnothing(m) ? "" : round(m.recall, digits=3),
                                      f(a.snr_before), f(snr_u), f(a.snr_after), length(missed), v,
                                      "\"$(Data.pulses_to_zaps(a.bad))\""], ","))
                end
                println("Zapisano do $out")
            catch e
                e isa InterruptException && rethrow()
                println("BŁĄD dla $psr: $(sprint(showerror, e))")
            end
        end
        println("\n===== PODSUMOWANIE ŚLEPEGO TESTU =====")
        foreach(println, summary)
        println("Wykresy: $(vpmout)<pulsar>/pulsar_autozap.png (co wycina automat), " *
                "pulsar_zap_compare.png (Ty vs automat)   Tabela: $out")
    end


    """
    Eksport danych pojedynczych impulsów (bez zapów) do debugowania automatu
    lokalnie: `input/zap_data/<psr>.f16.gz` (macierz N×nbin, Float16,
    wierszami) + `<psr>.json` (N, nbin, bin_st, bin_end). ~1.5–2 MB/pulsar.

      export_zap_data(["J0034-0721", "J0820-1350", ...])
    """
    function export_zap_data(psrs; vpmout="/home/psr/output/", dataroot="/home/psr/data/new/")
        dir = joinpath(@__DIR__, "input", "zap_data")
        mkpath(dir)
        for psr in psrs
            try
                outdir = vpmout * psr
                isfile(joinpath(outdir, "pulsar.debase.txt")) || prepare_debase(psr; vpmout=vpmout, dataroot=dataroot)
                p = Tools.read_params(joinpath(outdir, "params.json"))
                data = Data.load_ascii(joinpath(outdir, "pulsar.debase.txt"))
                N, nb = size(data)
                # Float16 tops out at 65504 -> store data/scale, scale in the json
                mx = maximum(abs, filter(isfinite, data); init=0.0)
                scale = mx > 0 ? mx / 30000 : 1.0
                raw = joinpath(dir, "$psr.f16")
                open(raw, "w") do io
                    write(io, Float16.(permutedims(data) ./ scale))   # row-major: pulse after pulse
                end
                run(`gzip -f $raw`)
                open(joinpath(dir, "$psr.json"), "w") do io
                    JSON.print(io, Dict("N" => N, "nbin" => nb, "bin_st" => p["bin_st"],
                                        "bin_end" => p["bin_end"], "scale" => scale))
                end
                println("$psr: $(N)×$(nb) → $(raw).gz ($(round(filesize(raw * ".gz") / 1e6, digits=1)) MB)")
            catch e
                e isa InterruptException && rethrow()
                println("BŁĄD dla $psr: $(sprint(showerror, e))")
            end
        end
        println("\nWyślij: git add input/zap_data input/zap_labels && git commit -m \"dane do zapow\" && git push")
    end


    # matched-filter SNR of the current pulsar.debase.txt (with its zaps),
    # or `nothing` when P3/on-pulse aren't known yet
    function coherent_snr(outdir)
        pf = joinpath(outdir, "params.json")
        tf = joinpath(outdir, "pulsar.debase.txt")
        (isfile(pf) && isfile(tf)) || return nothing
        p = Tools.read_params(pf)
        (get(p, "p3", -1.0) isa Number && p["p3"] > 0 && !isnothing(get(p, "bin_st", nothing))) || return nothing
        data = Data.load_ascii(tf)
        Data.zap!(data; ranges=Data.zaps_of(p))
        return P3FoldViterbi.coherent_fold(data, Float64(p["p3"]), Int(p["bin_st"]), Int(p["bin_end"]);
                                           warn_weak=false).snr
    end


    """
    Przerabia od nowa WSZYSTKIE pulsary z `dataroot` moją metodą zapowania
    do `vpmout` (domyślnie pazur2). Dla każdego: psradd → debase → automat
    proponuje złe impulsy + obrazek → Enter / p (popraw w pazi) / r (od nowa
    ręcznie) / n (bez zapów) / s (pomiń) → paz → debase → 2DFS/LRFS/P3-fold.

    Params: moje zapy w `"zaps2"` / `"zap_chans2"`; `"zaps"` nietknięte
    (przy nowym pulsarze przepisane z params Macieja albo puste). Z
    `ref_root` (pazur, Maciej — tylko odczyt) brane jest okno on-pulse (jeśli
    ustawione) i jego zapy (impulsy z zerową wagą w jego pulsar.spCF +
    jego "zaps") — tylko dla pulsarów z `"zapped": true`. Porównanie S/N
    profilu (te same dane i okno: bez zapów / moje / Macieja) zapisywane w
    params `"zap_cmp"`; zbiorczo: `compare_zap_methods()`.

      reprocess_all_zap()                                  # wszystkie
      reprocess_all_zap(only=["J1750-3503", "J2053-7200"]) # wybrane
      reprocess_all_zap(zap_mode=:auto)                    # automat bez pytań

    Wiersze są DOPISYWANE do `vpmout/zap_summary.csv`, a pulsary, które już
    tam są (status inny niż "blad"), są pomijane — po przerwaniu wystarczy
    uruchomić ponownie, żeby kontynuować. Żeby przerobić pulsar jeszcze raz:
    usuń jego wiersz albo `skip_done=false`.

    Status w tabeli: ok / brak_P3 / pominiety (s) / blad.
    """
    function reprocess_all_zap(; dataroot="/home/psr/data/new/", vpmout="/home/psr/data/OUTPUT/pazur2/",
                               ref_root="/home/psr/data/OUTPUT/pazur/",
                               only=nothing, zap_mode=:auto_pazi, zap_nsigma=5.0,
                               reverse=false, half=false, skip_done=true, warn_pct=10.0)
        endswith(vpmout, "/") || (vpmout *= "/")
        mkpath(vpmout)
        psr_dirs = sort(filter(d -> isdir(joinpath(dataroot, d)), readdir(dataroot)); rev=reverse)
        isnothing(only) || (psr_dirs = filter(in(Set(only)), psr_dirs))
        if half && !isempty(psr_dirs)
            n_all = length(psr_dirs)
            psr_dirs = psr_dirs[1:cld(n_all, 2)]
            println("Połowa listy: $(length(psr_dirs)) z $n_all pulsarów " *
                    "($(first(psr_dirs)) … $(last(psr_dirs)))")
        end
        println("Pulsarów do przerobienia: $(length(psr_dirs))   wyniki: $vpmout   " *
                "porównanie z: $(something(ref_root, "-"))")
        summary_file = joinpath(vpmout, "zap_summary.csv")
        header = "pulsar,status,n_pulses,n_zapped,pct_zapped,snr_none,snr_mine,snr_maciej,n_maciej,zaps2"
        done = Set{String}()
        if isfile(summary_file)
            for line in Iterators.drop(eachline(summary_file), 1)
                parts = split(line, ",")
                length(parts) >= 2 && parts[2] != "blad" && push!(done, parts[1])
            end
        else
            open(f -> println(f, header), summary_file, "w")
        end
        rows = String[]
        fmt(x) = isnothing(x) ? "" : x isa AbstractFloat ? string(round(x, digits=2)) : string(x)
        # dopisanie po każdym pulsarze — nic nie ginie przy przerwaniu
        record(row) = (push!(rows, row); open(f -> println(f, row), summary_file, "a"))

        for psr in psr_dirs
            if skip_done && psr in done
                println("SKIP $psr — już w $summary_file")
                continue
            end
            psr_path = joinpath(dataroot, psr)
            obs_dirs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
            if isempty(obs_dirs)
                @warn "Brak katalogu obserwacji dla $psr, pomijam"
                continue
            end
            indir  = joinpath(psr_path, obs_dirs[1]) * "/"
            outdir = vpmout * psr
            # Maciej's results (read only): on-pulse window + his zaps
            ref_dir = isnothing(ref_root) ? nothing : joinpath(ref_root, psr)
            ref_params = isnothing(ref_dir) ? nothing : joinpath(ref_dir, "params.json")
            ref_archive = isnothing(ref_dir) ? nothing : joinpath(ref_dir, "pulsar.spCF")
            println("\n===== $psr  (obs: $(obs_dirs[1]))" *
                    (!isnothing(ref_params) && isfile(ref_params) ? "  [jest wynik Macieja]" : "") * " =====")
            try
                pf = joinpath(outdir, "params.json")
                if !isfile(pf)
                    # new params; from Maciej's take ONLY the on-pulse window
                    # (if set) and his "zaps" as a record — his p3 etc. are
                    # not inherited (a placeholder p3 would switch off the
                    # P3 detection); "zaps" stays untouched, mine go to zaps2
                    mkpath(outdir)
                    p0 = Tools.default_params(pf)
                    if !isnothing(ref_params) && isfile(ref_params)
                        rp = Tools.read_params(ref_params)
                        for k in ("bin_st", "bin_end")
                            isnothing(get(rp, k, nothing)) || (p0[k] = rp[k])
                        end
                        p0["zaps"] = get(rp, "zaps", nothing)
                    end
                    Tools.save_params(pf, p0)
                end
                # compare only with pulsars Maciej marked as done ("zapped": true)
                ref_done = !isnothing(ref_params) && isfile(ref_params) &&
                           get(Tools.read_params(ref_params), "zapped", false) == true
                ra = (ref_done && isfile(ref_archive)) ? ref_archive : nothing
                ref_done && isnothing(ra) && println("UWAGA: Maciej ma \"zapped\": true, ale brak $ref_archive")
                p, _, _ = Data.process_psrdata(indir, outdir; zap_mode=zap_mode, zap_nsigma=zap_nsigma,
                                               ref_archive=ra)
                zaps = Data.zaps_of(p)
                nz = sum((b - a + 1 for (a, b) in Data.zaps_to_ranges(zaps)); init=0)
                np = something(get(p, "nsubint", nothing), 0)
                pct = np > 0 ? 100 * nz / np : 0.0
                c = get(p, "zap_cmp", Dict())
                status = (p["p3"] isa Number && p["p3"] > 0) ? "ok" : "brak_P3"
                println("$psr [$status]: odrzucono $nz/$np impulsów ($(round(pct, digits=1))%), " *
                        "S/N profilu: bez zapów $(fmt(get(c, "snr_none", nothing))) → moje " *
                        "$(fmt(get(c, "snr_mine", nothing)))" *
                        (haskey(c, "snr_maciej") ? ", Macieja $(fmt(c["snr_maciej"]))" : ""))
                pct > warn_pct && println("UWAGA: odrzucono ponad $(warn_pct)% impulsów — szum może " *
                                          "siedzieć w jednym kanale częstotliwości; rozważ ścieżkę _16 " *
                                          "z zapowaniem kanałów (paz -z) zamiast wyrzucania impulsów.")
                record("$psr,$status,$np,$nz,$(round(pct, digits=2)),$(fmt(get(c, "snr_none", nothing)))," *
                       "$(fmt(get(c, "snr_mine", nothing))),$(fmt(get(c, "snr_maciej", nothing)))," *
                       "$(fmt(get(c, "n_maciej", nothing))),\"$(something(zaps, ""))\"")
            catch e
                # Ctrl+C stops the whole run (not just this pulsar); the
                # interrupted pulsar is not recorded, so a rerun redoes it
                e isa InterruptException && rethrow()
                if e isa Data.SkipPulsar
                    println("$psr — pominięty na życzenie.")
                    record("$psr,pominiety,,,,,,,,\"\"")
                else
                    @warn "Błąd dla $psr: $e"
                    msg = replace(sprint(showerror, e), '"' => '\'', '\n' => ' ')
                    record("$psr,blad,,,,,,,,\"$(first(msg, 200))\"")
                end
            end
            if zap_mode in (:review, :pazi)
                print("Dalej? [Enter = następny, q = koniec]: ")
                strip(readline()) == "q" && break
            end
        end
        println("\nPodsumowanie (ta sesja), pełna tabela: $summary_file")
        println(header)
        foreach(println, rows)
        return rows
    end


    """
    W PEŁNI AUTOMATYCZNE wycięcie złych impulsów we wszystkich pulsarach,
    bez pytań i bez okien; struktura wyników jak u Macieja:
    `out/<psr>/params.json` + `out/<psr>/pulsar.spCF` (z wyciętymi
    impulsami, waga 0) — pełną analizę (2DFS/LRFS/P3-fold) można potem
    puścić normalnie na `out`.

    Dla każdego pulsara z `dataroot`:
      1. params = KOPIA params Macieja (`ref_root/<psr>`), jeśli są; inaczej
         domyślne. Brak okna on-pulse → z `p3_root/<psr>_16` (Czarek7); P3:
         Macieja, jeśli > 0, inaczej z `p3_root/<psr>_16`.
         Bez P3 albo bez okna on-pulse pulsar NIE jest teraz przerabiany
         (żaden folder nie powstaje) — trafia do `out/do_recznie.csv`
         (pulsar, powód, obserwacja) do zrobienia osobno na końcu.
      2. psradd → debase → automat (`review_zaps(mode=:auto)`): moje zapy
         w `"zaps2"`, `"zaps"` Macieja bez zmian, porównanie S/N profilu z
         jego zapami w `"zap_cmp"`, obrazki `pulsar_autozap.pdf` i `pulsar_contrast.pdf`.
      3. paz: wycięcie `zaps2` w pulsar.spCF; pliki pośrednie (debase) usunięte.

    Postęp w `out/zap_summary.csv` (ponowne uruchomienie kontynuuje).

      autozap_all()
      autozap_all(only=["J0034-0721"])
    """
    function autozap_all(; dataroot="/home/psr/data/new/", out="/home/psr/data/OUTPUT/pazur2/",
                         ref_root="/home/psr/data/OUTPUT/pazur/", p3_root="/home/psr/data/OUTPUT/Czarek7/",
                         only=nothing, skip_done=true, keep_png=true)
        endswith(out, "/") || (out *= "/")
        mkpath(out)
        psrs = sort(filter(d -> isdir(joinpath(dataroot, d)), readdir(dataroot)))
        isnothing(only) || (psrs = filter(in(Set(only)), psrs))
        summary_file = joinpath(out, "zap_summary.csv")
        header = "pulsar,status,n_pulses,n_zapped,pct_zapped,snr_none,snr_mine,snr_maciej,n_maciej,p3,p3_from,n_rfi,n_null,zaps2"
        done = Set{String}()
        if isfile(summary_file)
            for line in Iterators.drop(eachline(summary_file), 1)
                parts = split(line, ",")
                length(parts) >= 2 && parts[2] != "blad" && push!(done, parts[1])
            end
        else
            open(f -> println(f, header), summary_file, "w")
        end
        fmt(x) = isnothing(x) ? "" : x isa AbstractFloat ? string(round(x, digits=2)) : string(x)
        record(row) = open(f -> println(f, row), summary_file, "a")
        valid(x) = x isa Number && x > 0
        readp(f) = isfile(f) ? (try Tools.read_params(f) catch; nothing end) : nothing
        n_ok = n_skip = n_err = 0
        manual_file = joinpath(out, "do_recznie.csv")
        tmp_params = joinpath(mktempdir(), "params.json")
        manual = String[]
        println("AUTOMAT: $(length(psrs)) pulsarów → $out   (params/zapy Macieja: $ref_root, P3/on-pulse zapas: $p3_root)")

        for (k, psr) in enumerate(psrs)
            if skip_done && psr in done
                continue
            end
            psr_path = joinpath(dataroot, psr)
            obs = sort(filter(d -> isdir(joinpath(psr_path, d)), readdir(psr_path)))
            if isempty(obs)
                record("$psr,brak_danych,,,,,,,,,,,,\"\"")
                continue
            end
            indir, outdir = joinpath(psr_path, obs[1]) * "/", out * psr
            println("\n[$k/$(length(psrs))] ===== $psr =====")
            try
                pf = joinpath(outdir, "params.json")
                rp = readp(joinpath(ref_root, psr, "params.json"))
                cp16 = something(readp(joinpath(p3_root, psr * "_16", "params.json")),
                                 readp(joinpath(p3_root, psr, "params.json")), Dict())
                # 1. params: copy of Maciej's (same structure), else default
                #    (nothing is written before the checks below pass)
                p = !isnothing(rp) ? deepcopy(rp) : Tools.default_params(tmp_params)
                haskey(p, "zaps") || (p["zaps"] = nothing)
                for key in ("bin_st", "bin_end")
                    isnothing(get(p, key, nothing)) && (p[key] = get(cp16, key, nothing))
                end
                p3_from = valid(get(p, "p3", nothing)) && !isnothing(rp) ? "Maciej" : ""
                if p3_from == "" && valid(get(cp16, "p3", nothing))
                    p["p3"], p3_from = cp16["p3"], "Czarek7"
                    haskey(cp16, "p3_error") && (p["p3_error"] = cp16["p3_error"])
                end
                p3_from == "" && (p["p3"] = get(p, "p3", -1.0); valid(p["p3"]) || (p["p3"] = -1.0))
                p["nsubint"] = nothing          # recount on the rebuilt archive
                # no P3 / no on-pulse window anywhere -> not processed now,
                # listed in out/do_recznie.csv to be done separately at the end
                why = String[]
                (isnothing(p["bin_st"]) || isnothing(p["bin_end"])) && push!(why, "brak_onpulse")
                valid(p["p3"]) || push!(why, "brak_P3")
                if !isempty(why)
                    st = join(why, "+")
                    println("Pomijam ($st) — dopisany do $(manual_file) do zrobienia na końcu.")
                    record("$psr,$st,,,,,,,,$(fmt(p["p3"])),$p3_from,,,\"\"")
                    isfile(manual_file) || open(f -> println(f, "pulsar,powod,obs"), manual_file, "w")
                    open(f -> println(f, "$psr,$st,$(obs[1])"), manual_file, "a")
                    push!(manual, psr)
                    n_skip += 1
                    continue
                end
                mkpath(outdir)
                Tools.save_params(pf, p)

                # 2. archive + automat (no questions, no windows)
                outfile = joinpath(outdir, "pulsar.spCF")
                Data.add_psrfiles(indir, outfile)
                Data.get_nsubint(outfile, pf, p)
                Data.debase(outfile, pf, p)
                ref_archive = joinpath(ref_root, psr, "pulsar.spCF")
                ra = (!isnothing(rp) && get(rp, "zapped", false) == true && isfile(ref_archive)) ? ref_archive : nothing
                zaps = Data.review_zaps(joinpath(outdir, "pulsar.debase.txt"), p, pf, outdir;
                                        mode=:auto, archive=outfile, ref_archive=ra)
                p["zapped2"] = true
                Tools.save_params(pf, p)

                # 3. cut them in pulsar.spCF (zero weights, like Maciej's) + clean up
                isnothing(zaps) || Data.apply_zaps(outfile, zaps)
                for f in ("pulsar.debase.gg", "pulsar.debase.txt", "pulsar_zap_review.png")
                    rm(joinpath(outdir, f), force=true)
                end
                keep_png || foreach(f -> rm(joinpath(outdir, f), force=true), ("pulsar_autozap.pdf", "pulsar_contrast.pdf"))

                c = get(p, "zap_cmp", Dict())
                np, nz = get(c, "n_pulses", 0), get(c, "n_mine", 0)
                pct = np > 0 ? 100 * nz / np : 0.0
                record("$psr,ok,$np,$nz,$(round(pct, digits=2)),$(fmt(get(c, "snr_none", nothing)))," *
                       "$(fmt(get(c, "snr_mine", nothing))),$(fmt(get(c, "snr_maciej", nothing)))," *
                       "$(fmt(get(c, "n_maciej", nothing))),$(fmt(p["p3"])),$p3_from,$(fmt(get(c, "n_rfi", nothing)))," *
                       "$(fmt(get(c, "n_null", nothing))),\"$(something(zaps, ""))\"")
                n_ok += 1
            catch e
                e isa InterruptException && rethrow()
                n_err += 1
                msg = replace(sprint(showerror, e), '"' => '\'', '\n' => ' ')
                println("BŁĄD dla $psr: $(first(msg, 300))")
                record("$psr,blad,,,,,,,,,,,,\"$(first(msg, 200))\"")
            end
        end
        println("\nAUTOMAT gotowy: ok=$n_ok  pominięte=$n_skip  błędy=$n_err   tabela: $summary_file")
        isempty(manual) || println("Do zrobienia osobno (brak P3 / on-pulse), $(length(manual)) w tej sesji: " *
                                   "$(join(manual, ", "))\nPełna lista: $manual_file")
    end


    """
    Końcowe porównanie: moje zapowanie (pazur2, `zaps2`) vs Macieja (pazur).

    Czyta `p["zap_cmp"]` z `mine/<psr>/params.json` — liczone przy
    przerabianiu na TYCH SAMYCH danych i tym samym oknie on-pulse:
    standardowe S/N profilu (bez zapów / moje / Macieja) i miara automatu.
    Zapisuje do `outdir` (poza folderami wyników):
      zap_comparison.csv  — tabela per pulsar,
      zap_comparison.pdf  — S/N moje vs Macieja, zysk S/N, % wyciętych.

      compare_zap_methods()
    """
    function compare_zap_methods(; mine="/home/psr/data/OUTPUT/pazur2/",
                                 outdir="/home/psr/data/OUTPUT/zap_comparison/", show_=false)
        mkpath(outdir)
        T = []
        for psr in sort(readdir(mine))
            pf = joinpath(mine, psr, "params.json")
            isfile(pf) || continue
            c = get(Tools.read_params(pf), "zap_cmp", nothing)
            (isnothing(c) || !haskey(c, "snr_maciej") || isnothing(c["snr_maciej"])) && continue
            g(k) = (v = get(c, k, nothing); isnothing(v) ? NaN : Float64(v))
            push!(T, (psr=psr, N=Int(c["n_pulses"]), n_mine=Int(c["n_mine"]), n_maciej=Int(c["n_maciej"]),
                      none=g("snr_none"), mine=g("snr_mine"), maciej=g("snr_maciej"),
                      a_none=g("auto_none"), a_mine=g("auto_mine"), a_maciej=g("auto_maciej")))
        end
        isempty(T) && (println("Brak pulsarów z porównaniem (zap_cmp z wynikiem Macieja) w $mine"); return nothing)

        csv = joinpath(outdir, "zap_comparison.csv")
        open(csv, "w") do f
            println(f, "pulsar,n_pulses,n_mine,n_maciej,snr_none,snr_mine,snr_maciej,gain_mine_pct," *
                       "gain_maciej_pct,better,auto_none,auto_mine,auto_maciej")
            for t in T
                gm, gj = 100 * (t.mine / t.none - 1), 100 * (t.maciej / t.none - 1)
                better = abs(t.mine - t.maciej) <= 0.005 * t.none ? "remis" : t.mine > t.maciej ? "moja" : "Macieja"
                println(f, join([t.psr, t.N, t.n_mine, t.n_maciej, round(t.none, digits=2), round(t.mine, digits=2),
                                 round(t.maciej, digits=2), round(gm, digits=2), round(gj, digits=2), better,
                                 round(t.a_none, digits=2), round(t.a_mine, digits=2), round(t.a_maciej, digits=2)], ","))
            end
        end

        gm = [100 * (t.mine / t.none - 1) for t in T]
        gj = [100 * (t.maciej / t.none - 1) for t in T]
        tol = [0.5 for _ in T]   # remis: różnica < 0.5% S/N bez zapów
        nb = count(gm .- gj .> tol); nw = count(gj .- gm .> tol); ne = length(T) - nb - nw
        println("\n===== MOJA METODA vs MACIEJA ($(length(T)) pulsarów, standardowe S/N profilu) =====")
        println("lepsza moja: $nb   lepsza Macieja: $nw   remis (<0.5%): $ne")
        println("średni zysk S/N: moja $(round(mean(gm), digits=1))%   Macieja $(round(mean(gj), digits=1))%   " *
                "(mediana: $(round(median(gm), digits=1))% / $(round(median(gj), digits=1))%)")
        println("średnio wycięte: moja $(round(mean(t.n_mine / t.N for t in T) * 100, digits=1))%   " *
                "Macieja $(round(mean(t.n_maciej / t.N for t in T) * 100, digits=1))%")
        png = Plot.zap_methods_comparison(T, outdir; show_=show_)
        println("Tabela: $csv\nWykres: $png")
        return T
    end

    """
        autozap_test(psrs; out="/home/psr/data/OUTPUT/pazur2_test/", ...)

    Test before the full run: `autozap_all` on a few representative pulsars
    into a SEPARATE folder (pazur2 stays clean for the full run), the
    comparison with Maciej (`out/zap_comparison/`), then the pictures one by
    one (Enter = next): what the automat cut and why, and Maciej vs automat.
    """
    function autozap_test(psrs; out="/home/psr/data/OUTPUT/pazur2_test/", dataroot="/home/psr/data/new/",
                          ref_root="/home/psr/data/OUTPUT/pazur/", p3_root="/home/psr/data/OUTPUT/Czarek7/",
                          process=true, show_=true)
        out = rstrip(out, '/') * "/"
        process && autozap_all(dataroot=dataroot, out=out, ref_root=ref_root, p3_root=p3_root, only=psrs,
                           skip_done=false)
        T = compare_zap_methods(mine=out, outdir=out * "zap_comparison/", show_=false)
        r1(x) = isnothing(x) ? "-" : string(round(x, digits=1))
        println("
===== TEST: podsumowanie =====")
        println(rpad("pulsar", 13), lpad("impulsy", 8), lpad("auto", 7), lpad("auto%", 7), lpad("Maciej", 8),
                lpad("S/N bez", 10), lpad("auto", 9), lpad("Maciej", 9))
        for psr in psrs
            pf = joinpath(out, psr, "params.json")
            isfile(pf) || (println(rpad(psr, 13), "  (brak wyniku - zob. do_recznie.csv / zap_summary.csv)"); continue)
            c = get(Tools.read_params(pf), "zap_cmp", Dict())
            N, nm = get(c, "n_pulses", 0), get(c, "n_mine", 0)
            println(rpad(psr, 13), lpad(N, 8), lpad(nm, 7), lpad(r1(N > 0 ? 100 * nm / N : nothing), 7),
                    lpad(string(get(c, "n_maciej", "-")), 8), lpad(r1(get(c, "snr_none", nothing)), 10),
                    lpad(r1(get(c, "snr_mine", nothing)), 9), lpad(r1(get(c, "snr_maciej", nothing)), 9))
        end
        show_ || return T
        println("
Obrazki po kolei (Enter = następny, s + Enter = koniec oglądania).")
        for psr in psrs, f in ("pulsar_autozap.pdf", "pulsar_contrast.pdf")
            png = joinpath(out, psr, f)
            isfile(png) || continue
            println("$psr: $f")
            proc = Plot.open_png(png)
            print("Enter = następny: ")
            ans = strip(readline())
            Plot.close_png(proc)
            ans == "s" && break
        end
        cpng = joinpath(out, "zap_comparison", "zap_comparison.pdf")
        isfile(cpng) && Plot.show_png(cpng; prompt="Wykres zbiorczy - Enter = koniec: ")
        return T
    end

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

        # PSR J0134-2937
        #process_psrdata_16("/home/psr/data/new/J0134-2937/2019-10-18-22:46:59/", vpmout*"J0134-2937_16")
        #Data.analyse_p3folds_16_new(vpmout*"J0134-2937_16", "norefine"; n_comp=2)


        # PSR J0151-0635
        #process_psrdata_16("/home/psr/data/new/J0151-0635/2020-04-13-10:00:14/", vpmout*"J0151-0635_16")
        #Data.analyse_p3folds_16(vpmout*"J0151-0635_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J0151-0635_16", "norefine"; n_comp=2)

        # PSR J0255-5304
        #process_psrdata_16("/home/psr/data/new/J0255-5304/2019-10-18-23:29:57/", vpmout*"J0255-5304_16")
        #Data.analyse_p3folds_16(vpmout*"J0255-5304_16", "norefine")

        # PSR J0421-0345
        #process_psrdata_16("/home/psr/data/new/J0421-0345/2019-10-27-21:58:44/", vpmout*"J0421-0345_16")
        #Data.analyse_p3folds_16(vpmout*"J0421-0345_16", "norefine")
        #process_psrdata(vpmout*"J0421-0345", vpmout*"J0421-0345") #

        # PSR J0630-2834
        #process_psrdata_16("/home/psr/data/new/J0630-2834/2019-11-18-21:05:35/", vpmout*"J0630-2834_16")
        #Data.analyse_p3folds_16(vpmout*"J0630-2834_16", "norefine")
        #process_psrdata(vpmout*"J0630-2834", vpmout*"J0630-2834") # not nice

        # PSR J0820-1350
        #process_psrdata("/home/psr/data/uhf/J0820-1350/2020-07-11-10:30:12/816/single/", vpmout*"J0820-1350_uhf")
        #process_psrdata_16("/home/psr/data/new/J0820-1350/2020-01-11-01:05:56/", vpmout*"J0820-1350_16")
        #Data.analyse_p3folds_16(vpmout*"J0820-1350_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J0820-1350_16", "norefine", n_comp=2)
        #process_psrdata_single("/home/psr/data/new/J0820-1350/2020-01-11-01:05:56/", vpmout*"J0820-1350_single") # very bad no signal at least at the beginning
        #process_psrdata(vpmout*"J0820-1350", vpmout*"J0820-1350") # nice

        # PSR J0823+0159
        #process_psrdata("/home/psr/data/uhf/J0823+0159/2023-02-27-20:38:58/816/single", vpmout*"J0823+0159_uhf")
        #process_psrdata_16("/home/psr/data/new/J0823+0159/2019-11-05-23:33:48/", vpmout*"J0823+0159_16")
        #process_psrdata_single("/home/psr/data/new/J0823+0159/2019-11-05-23:33:48/", vpmout*"J0823+0159_single")
        #process_psrdata(vpmout*"J0823+0159", vpmout*"J0823+0159") #  nice single pulses, but p3folds no

        # PSR J0856-6137
        #process_psrdata_16("/home/psr/data/new/J0856-6137/2020-01-04-19:54:45/", vpmout*"J0856-6137_16")
        #Data.analyse_p3folds_16(vpmout*"J0856-6137_16", "norefine")

        # PSR J0904-7459
        #process_psrdata_16("/home/psr/data/new/J0904-7459/2019-10-19-05:41:13", vpmout*"J0904-7459_16")
        #Data.analyse_p3folds_16(vpmout*"J0904-7459_16", "norefine")

        # PSR J0934-5249
        #process_psrdata_16("/home/psr/data/new/J0934-5249/2019-11-05-23:07:57/", vpmout*"J0934-5249_16")
        #Data.analyse_p3folds_16(vpmout*"J0934-5249_16", "norefine")

        # PSR J0959-4809
        #process_psrdata_16("/home/psr/data/new/J0959-4809/2019-10-19-06:23:58/", vpmout*"J0959-4809_16")
        #Data.analyse_p3folds_16(vpmout*"J0959-4809_16", "norefine")

        # PSR J1034-3224
        #process_psrdata_16("/home/psr/data/new/J1034-3224/2019-10-19-07:31:26/", vpmout*"J1034-3224_16")
        #Data.analyse_p3folds_16(vpmout*"J1034-3224_16", "norefine")
        #process_psrdata(vpmout*"J1034-3224", vpmout*"J1034-3224") # single bad, p3folds bad, P3 stable

        # PSR J1048-5832
        #process_psrdata_16("/home/psr/data/new/J1048-5832/2020-08-29-13:00:11/", vpmout*"J1048-5832_16")
        #Data.analyse_p3folds_16(vpmout*"J1048-5832_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1048-5832_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1048-5832_16")
        #Data.geometry_analysis(vpmout*"J1048-5832_16")

        # PSR J1110-5637
        #process_psrdata("/home/psr/data/new/J1110-5637/2019-10-19-08:18:28/", vpmout*"J1110-5637")
        #process_psrdata_16("/home/psr/data/new/J1110-5637/2019-10-19-08:18:28/", vpmout*"J1110-5637_16")
        #Data.analyse_p3folds_16(vpmout*"J1110-5637_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1110-5637_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1110-5637_16")
        #Data.geometry_analysis(vpmout*"J1110-5637_16")
        #phase_modulation(vpmout*"J1110-5637")
        #p3fold_refine(vpmout*"J1110-5637", continuity_weight=0.05)
        #p3fold_coherent(vpmout*"J1110-5637")

        # PSR J1114-6100
        #process_psrdata_16("/home/psr/data/new/J1114-6100/2019-10-19-08:30:30/", vpmout*"J1114-6100_16")
        #Data.analyse_p3folds_16(vpmout*"J1114-6100_16c", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1114-6100_16c", "norefine")

        # PSR J1133-6250
        #process_psrdata_16("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_16")
        #process_psrdata_single("/home/psr/data/new/J1133-6250/2019-11-06-00:46:43/", vpmout*"J1133-6250_single") # some offset? probably not
        #process_psrdata(vpmout*"J1133-6250", vpmout*"J1133-6250") # single not stable, P. p3fold (10 ybins)
        #phase_modulation(vpmout*"J1133-6250")
        #p3fold_coherent(vpmout*"J1133-6250")

        # PSR J1137-6700
        #process_psrdata_16("/home/psr/data/new/J1137-6700/2019-11-06-00:36:02/", vpmout*"J1137-6700_16")
        #Data.analyse_p3folds_16_new(vpmout*"J1137-6700_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1137-6700_16")

        # PSR J1232-4742
        #process_psrdata_16("/home/psr/data/new/J1232-4742/2020-04-11-22:21:18/", vpmout*"J1232-4742_16")
        #Data.analyse_p3folds_16(vpmout*"J1232-4742_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1232-4742_16", "norefine")

        # PSR J1312-5402
        #process_psrdata_16("/home/psr/data/new/J1312-5402/2020-03-30-19:11:06/", vpmout*"J1312-5402_16")
        #Data.analyse_p3folds_16(vpmout*"J1312-5402_16", "norefine")

        # PSR J1345-6115
        #Data.analyse_p3folds_16_new(vpmout*"J1345-6115_16", "norefine"; n_comp=2)

        # PSR J1414-6802
        #process_psrdata_16("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_16") # start diging there...
        #process_psrdata("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_test") # start diging there...
        #process_psrdata_single("/home/psr/data/new/J1414-6802/2020-07-24-16:17:58/", vpmout*"J1414-6802_single") # weak signal no offset
        #process_psrdata(vpmout*"J1414-6802", vpmout*"J1414-6802") # start diging there...

        # PSR J1512-5431
        #process_psrdata_16("/home/psr/data/new/J1512-5431/2020-04-11-23:16:42/", vpmout*"J1512-5431_16")
        #Data.analyse_p3folds_16(vpmout*"J1512-5431_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1512-5431_16", "norefine", n_comp=2)

        # PSR J1514-4834
        #process_psrdata_16("/home/psr/data/new/J1514-4834/2019-12-16-08:22:16/", vpmout*"J1514-4834_16")
        #process_psrdata(vpmout*"J1514-4834", vpmout*"J1514-4834") #

        # PSR J1524-5706
        #Data.analyse_p3folds_16(vpmout*"JJ1524-5706_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1524-5706_16", "norefine")

        # PSR J1534-4428
        #process_psrdata_16("/home/psr/data/new/J1534-4428/2020-01-03-01:03:02/", vpmout*"J1534-4428_16")
        #Data.analyse_p3folds_16(vpmout*"J1534-4428_16", "norefine")

        # PSR J1555-3134
        #process_psrdata_16("/home/psr/data/new/J1555-3134/2019-10-31-06:34:39/", vpmout*"J1555-3134_16")
        #Data.analyse_p3folds_16(vpmout*"J1555-3134_16", "norefine")
        #Data.analyse_p3folds_16(vpmout*"J1555-3134_16", "refine")
        #process_psrdata_single("/home/psr/data/new/J1555-3134/2019-10-31-06:34:39/", vpmout*"J1555-3134_single") # no offset
        #process_psrdata(vpmout*"J1555-3134", vpmout*"J1555-3134") # single not stable, P. p3fold (10 ybins)

        # PSR J1651-5222
        #process_psrdata_16("/home/psr/data/new/J1651-5222/2019-10-31-15:25:01/", vpmout*"J1651-5222_16")
        #Data.analyse_p3folds_16(vpmout*"J1651-5222_16", "norefine")

        # PSR J1651-7642
        #Data.analyse_p3folds_16_new(vpmout*"J1651-7642_16", "norefine"; n_comp=2)
        #Data.position_angle(vpmout*"J1651-7642_16")

        # PSR J1703-1846
        #process_psrdata_16("/home/psr/data/new/J1703-1846/2019-11-05-16:26:36/", vpmout*"J1703-1846_16")
        #Data.analyse_p3folds_16(vpmout*"J1703-1846_16", "norefine")


        # PSR J1741-0840
        #process_psrdata_16("/home/psr/data/new/J1741-0840/2019-11-05-16:44:10/", vpmout*"J1741-0840_16")
        #Data.analyse_p3folds_16(vpmout*"J1741-0840_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1741-0840_16", "norefine")

        # PSR J1744-1610
        #test3("/home/psr/data/new/J1744-1610/2020-06-12-20:00:13/", vpmout)

        # PSR J1807+0756
        #process_psrdata_16("/home/psr/data/new/J1807+0756/2020-06-07-22:03:05/", vpmout*"J1807+0756_16")
        #Data.analyse_p3folds_16(vpmout*"J1807+0756_16", "refine")

        # PSR J1834-0010
        #process_psrdata(vpmout*"J1834-0010", vpmout*"J1834-0010") # P. nice

        # PSR J1842-0359
        #process_psrdata_16("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_16")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_16", "refine")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1842-0359_16", "norefine", n_comp=4)
        #process_psrdata("/home/psr/data/uhf/J1842-0359/2020-07-25-19:58:19/816/single/", vpmout*"J1842-0359_uhf")
        #process_psrdata_16("/home/psr/data/uhf/J1842-0359/2020-07-25-19:58:19/816/single/", vpmout*"J1842-0359_uhf_16")
        #Data.analyse_p3folds_16(vpmout*"J1842-0359_uhf_16", "norefine")
        #process_psrdata_single("/home/psr/data/new/J1842-0359/2019-11-05-18:03:43/", vpmout*"J1842-0359_single")
        #process_psrdata(vpmout*"J1842-0359", vpmout*"J1842-0359") #  nice single pulses, nice P. p3fold
        #Data.position_angle(vpmout*"J1842-0359_16")
        #Data.geometry_analysis(vpmout*"J1842-0359_16")

        # PSR J1857+0057
        #process_psrdata_16("/home/psr/data/new/J1857+0057/2020-02-03-06:00:48/", vpmout*"J1857+0057_16")
        #Data.analyse_p3folds_16(vpmout*"J1857+0057_16", "norefine")

        # PSR J1900-7951
        #process_psrdata_16("/home/psr/data/new/J1900-7951/2020-05-31-22:21:41/", vpmout*"J1900-7951_16")
        #Data.analyse_p3folds_16(vpmout*"J1900-7951_16", "norefine")

        # PSR J1901-0906
        #process_psrdata_16("/home/psr/data/new/J1901-0906/2019-12-16-13:57:56/", vpmout*"J1901-0906_16")    # P3=3
        #Data.analyse_p3folds_16(vpmout*"J1901-0906_16", "norefine")
        #process_psrdata_16("/home/psr/data/new/J1901-0906/2019-12-16-13:57:56/", vpmout*"J1901-0906_16b")   # P3=5
        #Data.analyse_p3folds_16(vpmout*"J1901-0906_16b", "norefine")

        # PSR J1903+2225
        #process_psrdata_16("/home/psr/data/new/J1903+2225/2020-08-20-20:53:16/", vpmout*"J1903+2225_16")
        #Data.analyse_p3folds_16_new(vpmout*"J1903+2225_16", "norefine"; n_comp=1)

        # PSR J1904-1224
        #process_psrdata_16("/home/psr/data/new/J1904-1224/2020-04-14-05:57:45", vpmout*"J1904-1224_16")
        #Data.analyse_p3folds_16(vpmout*"J1904-1224_16", "norefine")

        # PSR J1907+0731
        #process_psrdata(vpmout*"J1907+0731", vpmout*"J1907+0731") #

        # PSR J1919+0134
        #test3("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #process_psrdata_16("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout*"J1919+0134_16")
        #Data.analyse_p3folds_16(vpmout*"J1919+0134_16", "norefine")

        # PSR J1921+2153
        #process_psrdata_16("/home/psr/data/new/J1921+2153/2020-08-08-19:17:21/", vpmout*"J1921+2153_16")
        #Data.analyse_p3folds_16(vpmout*"J1921+2153_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J1921+2153_16", "norefine", n_comp=2)
        #process_psrdata(vpmout*"J1921+2153", vpmout*"J1921+2153") # not nice

        # PSR J2046+1540
        #process_psrdata_16("/home/psr/data/new/J2046+1540/2019-11-26-17:18:47/", vpmout*"J2046+1540_16")
        #Data.analyse_p3folds_16_new(vpmout*"J2046+1540_16", "norefine", n_comp=2)
        #Data.position_angle(vpmout*"J2046+1540_16")


        # PSR J2053-7200
        #process_psrdata("/home/psr/data/new/J2053-7200/2019-11-26-18:06:44/", vpmout*"J2053-7200")
        #process_psrdata_16("/home/psr/data/new/J2053-7200/2019-11-26-18:06:44/", vpmout*"J2053-7200_16")
        #Data.analyse_p3folds_16(vpmout*"J2053-7200_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J2053-7200_16", "norefine", n_comp=2)
        #Data.position_angle(vpmout*"J2053-7200_16")
        #Data.geometry_analysis(vpmout*"J2053-7200_16")
        #phase_modulation(vpmout*"J2053-7200")
        #phase_modulation2(vpmout*"J2053-7200")
        #phase_modulation2a(vpmout*"J2053-7200")
        #phase_modulation3(vpmout*"J2053-7200")
        #p3fold_coherent(vpmout*"J2053-7200")

        # PSR J2139+2242
        #process_psrdata("/home/psr/data/new/J2139+2242/2020-09-07-20:01:41/", vpmout*"J2139+2242")
        #process_psrdata_16("/home/psr/data/new/J2139+2242/2020-09-07-20:01:41/", vpmout*"J2139+2242_16")
        #Data.analyse_p3folds_16(vpmout*"J2139+2242_16", "norefine")
        #Data.analyse_p3folds_16_new(vpmout*"J2139+2242_16", "norefine", n_comp=2)
        #phase_modulation(vpmout*"J2139+2242")
        #p3fold_coherent(vpmout*"J2139+2242")
        #p3fold_norefine_compare(vpmout*"J2139+2242"; darkness=1.0)


        # PSR J1224-6407
        #process_psrdata_16("/home/psr/data/new/J1224-6407/2019-10-28-03:47:09/", vpmout*"J1224-6407_16")
        #Data.analyse_p3folds_16_new(vpmout*"J1224-6407_16", "norefine", n_comp=2)

        # PSR J1703-4442
        #process_psrdata_16("/home/psr/data/new/J1703-4442/2020-03-11-02:22:24/", vpmout*"J1703-4442_16")
        #Data.analyse_p3folds_16_new(vpmout*"J1703-4442_16", "norefine", n_comp=2)


        # PSR J0304+1932 # AVERAGE
        #process_psrdata_16("/home/psr/data/new/J0304+1932/2019-10-18-23:05:44/", vpmout*"J0304+1932_16")
        #Data.analyse_p3folds_16_new(vpmout*"J0304+1932_16", "norefine", n_comp=2)
        #Data.analyse_average_offset(vpmout*"J0304+1932_16";n_comp=2, npulse=150)

        # PSR J1429-5935
        #process_psrdata("/home/psr/data/new/J1429-5935/2020-02-14-04:21:11", vpmout*"J1429-5935")
        #phase_modulation(vpmout*"J1429-5935")
        #p3fold_coherent(vpmout*"J1429-5935")


        # PSR J1750-3503
        #process_psrdata("/home/psr/data/new/J1750-3503/2019-12-14-14:22:12/", vpmout*"J1750-3503")
        #phase_modulation(vpmout*"J1750-3503")
        #phase_modulation2(vpmout*"J1750-3503")
        #phase_modulation3(vpmout*"J1750-3503")
        #p3fold_coherent(vpmout*"J1750-3503", lowpass_cutoff=1/300)

        # TEST: skan lowpass_cutoff + finalny coherent fold dla kilku pulsarów
        # po kolei. Dane już przetworzone (process_psrdata było wcześniej
        # odpalone) dla wszystkich czterech. Uwaga na J2053-7200: krótkie
        # P3≈3.06, mocno "wobbluje" — opisany wcześniej w kodzie jako trudny
        # przypadek dla metod opartych na pojedynczym globalnym binie FFT
        # (drift_test daje tu tylko 0.5σ, mimo że phase_modulation3 w
        # krótkich oknach widzi modulację na 80σ) — spodziewaj się niższego
        # SNR niż dla J1539-6322, to oczekiwane, nie błąd. J1133-6250 ma
        # wcześniejszą notatkę "single not stable" — może też wypaść słabo.
        # J2053-7200: P3≈3.06 jest tak krótkie, że p["p3_ybins"] z pliku
        # (często=10) dawało artefakt aliasingu (sztuczny "szachownicowy"
        # wzór w foldzie) — p3fold_cutoff_scan/p3fold_coherent liczą teraz
        # wlasciwa wartosc automatycznie przez Functions.find_ybins, zamiast
        # slepo ufac plikowi, wiec nie trzeba juz nic wymuszac recznie.
        # TEST zadan 1-3 (dedupe cutoffow, siatka zalezna od P3, linia
        # odniesienia p3_nominal): odpal skan na dwoch roznych pulsarach,
        # zeby zobaczyc czy nowa P3-zalezna siatka wyglada sensownie dla
        # dwoch bardzo roznych wartosci P3, i czy dedupe/linia dzialaja
        # TEST KONCOWY (wersja ostateczna): kilku zroznicowanych pulsarow,
        # pelny pipeline - skan cutoffu (P3-zalezna siatka, dedupe, mapa
        # heatmap) + finalny fold (spojny p3_window=60, nowy bin_search,
        # spojnosc podpasm). Dobor celowo zroznicowany:
        #   J1750-3503 - slaby sygnal (SNR~3), znany odwracajacy sie dryf
        #   J1539-6322 - mocny sygnal (SNR~26), kontrast dla powyzszego
        #   J1110-5637 - "podrecznikowy" czysty przypadek z dokumentacji
        #   J2053-7200 - bardzo krotkie P3≈3, trudny przypadek dla ybins
        #
        # ETAP 1 — OZNACZANIE zlych impulsow recznie (do dostrojenia
        # detektora). Wyniki: input/zap_labels/<pulsar>.csv -> git push.
        # Dopisz tu pulsary, ktore chcesz oznaczyc (rozne typy: czyste,
        # zaszumione, z RFI, z nullami). Brak pulsar.debase.txt -> najpierw
        # psradd + debase (prepare_debase). Ctrl+C przerywa calosc; juz
        # zapisane CSV zostaja. Pulsar z gotowym CSV jest pomijany.
        label_psrs = ["J0034-0721", "J0304+1932", "J0459-0210", "J0624-0424", "J0820-1350",
                      "J2346-0609", "J1110-5637", "J1750-3503", "J1539-6322"]
        # (wylaczone — oznaczanie reczne w pazi juz niepotrzebne)
        #for psr in label_psrs
        #    isfile(joinpath(@__DIR__, "input", "zap_labels", "$psr.csv")) && continue
        #    label_bad_pulses(psr; vpmout=vpmout)
        #end

        # ETAP 2 — SLEPY TEST automatu (domyslne parametry, bez wiedzy o
        # Twoich zapach), porownanie z etykietami + ocena w pazi:
        # input/zap_labels/blind_test.csv, <psr>/pulsar_zap_compare.png
        # tylko automat: bez pazi i bez pytan, wykres po kazdym pulsarze
        # (Enter = nastepny), na koncu tabela podsumowania
        # pierwsze 5 znanych (z J0820) + 4 nowe (sprawdzenie na pulsarach, ktorych
        # automat "nie widzial" przy strojeniu)
        test_psrs = vcat(label_psrs[1:5], ["J1903+0135", "J1901-0906", "J1946+1805", "J1555-3134"])
        #blind_test_autozap(["J0955-5304"]; vpmout=vpmout, interactive=false)
        #blind_test_autozap(test_psrs; vpmout=vpmout, interactive=false)

        # eksport danych (bez zapow) do debugowania automatu lokalnie
        #export_zap_data(label_psrs; vpmout=vpmout)

        # ETAP 2b — strojenie parametrow na etykietach (po slepym tescie)
        #tune_autozap(label_psrs; vpmout=vpmout)

        # ETAP 3 — PRZERABIANIE WSZYSTKICH pulsarow moja metoda do pazur2
        # (pusty, moj). Moje zapy -> params "zaps2"; "zaps" Macieja nietkniete.
        # Z pazur (Maciej, tylko odczyt) brane: okno on-pulse + jego zapy do
        # porownania (S/N profilu na tych samych danych). Dla kazdego
        # pulsara: obrazek, potem Enter = zatwierdz | p = popraw w pazi |
        # r = od nowa recznie | n = bez zapow | s = pomin. Ctrl+C przerywa;
        # ponowne uruchomienie kontynuuje (pazur2/zap_summary.csv).
        # (wersja z pytaniami przy kazdym pulsarze — na pozniej)
        #reprocess_all_zap(vpmout="/home/psr/data/OUTPUT/pazur2/", ref_root="/home/psr/data/OUTPUT/pazur/",
        #                  zap_mode=:auto_pazi)

        # ETAP 3a — W PELNI AUTOMATYCZNIE (bez pytan i okien): params kopiowane
        # od Macieja (+ on-pulse/P3 z Czarek7 _16, gdy brak), automat wycina
        # zle impulsy -> zaps2, pulsar.spCF z wycietymi impulsami; struktura
        # jak w pazur. Ctrl+C przerywa, ponowne uruchomienie kontynuuje.
        # TEST przed pelnym przebiegiem: reprezentatywne pulsary do
        # /home/psr/data/OUTPUT/pazur2_test/ + porownanie z Maciejem + obrazki
        #autozap_test(["J0034-0721",   # silny nulling
        #             "J0304+1932",   # czysty, jasny (nic nie ciac)
        #             "J0459-0210",   # nulle / wygaszenia
        #             "J0820-4114",   # szerokie okno on-pulse
        #             "J1046-5813",   # Maciej cial pojedyncze impulsy
        #             "J1048-5832",   # jasny, silna modulacja (wczesniej za duzo ciete)
        #             "J1049-5833",   # dlugie odcinki bez sygnalu
        #             "J1054-6452",   # Maciej cial 31 impulsow
        #             "J1056-6258",   # bardzo jasny, czysty
        #             "J1057-5226",   # interpuls
        #             "J1110-5637",   # wygaszony srodek
        #             "J1539-6322",   # zaszumiony koniec
        #             "J1750-3503"])  # paski RFI

        # ETAP 3a — PELNY PRZEBIEG (test z lokalna regula RFI zaakceptowany: 6:0:5 vs Maciej):
        # params kopiowane od Macieja (+ on-pulse/P3 z Czarek7 _16, gdy brak),
        # automat wycina zle impulsy -> zaps2; struktura jak w pazur.
        # Ctrl+C przerywa, ponowne uruchomienie kontynuuje.
        autozap_all(dataroot="/home/psr/data/new/", out="/home/psr/data/OUTPUT/pazur2/",
                    ref_root="/home/psr/data/OUTPUT/pazur/", p3_root="/home/psr/data/OUTPUT/Czarek7/")

        # ETAP 4 — porownanie moja vs Macieja (tabela + wykres poza
        # folderami wynikow: /home/psr/data/OUTPUT/zap_comparison/)
        compare_zap_methods(mine="/home/psr/data/OUTPUT/pazur2/",
                            outdir="/home/psr/data/OUTPUT/zap_comparison/", show_=false)

        # test skanu cutoffu + coherent fold na 4 pulsarach (wylaczony)
        #for psr_test in ["J1750-3503", "J1539-6322", "J1110-5637", "J2053-7200"]
        #    println("\n===== $psr_test =====")
        #    try
        #        scan_result = p3fold_cutoff_scan(vpmout*psr_test)
        #        if !isnothing(scan_result.suggested_cutoff)
        #            println(">>> OSTATECZNIE WYBRANY cutoff dla $psr_test: $(scan_result.suggested_cutoff)")
        #            p3fold_coherent(vpmout*psr_test, lowpass_cutoff=scan_result.suggested_cutoff)
        #        end
        #    catch e
        #        println("ERROR dla $psr_test: $e")
        #    end
        #end


        #test(vpmout)
        #test2(vpmout)
        #Data.process_all_data(vpmout)
        #Data.combine_pngs_to_pdf(vpmout)
        #Data.combine_pngs(vpmout)
        #Data.remove_folders(vpmout)
        #Data.remove_notinteresting("input/pulsars_interesting.txt", vpmout)

        #Tools.clean_all(vpmout)
        #analyse_all()

        # pętla po reprocess_list.csv wyniesiona do reprocess_half_list() —
        # wywołaj ją ręcznie z REPL-a kiedy będziesz chciał przerabiać listę
    end

    """
    Pętla po `input/reprocess_list.csv` (tylko pierwsza połowa listy) —
    interaktywnie pyta o n_comp per pulsar, woła `Data.analyse_p3folds_16_new`,
    pomija pulsary już zapisane w `input/separations.csv` (poza `force_redo`).
    Wywołaj ręcznie z REPL-a, np.:
      reprocess_half_list()
    """
    function reprocess_half_list(; vpmout="/home/psr/output/")
        todo_file = joinpath(@__DIR__, "input", "reprocess_list.csv")
        sep_file  = joinpath(@__DIR__, "input", "separations.csv")
        summary_dir = joinpath(@__DIR__, "output", "offset_summaries")

        # pulsary do wymuszonego przerobienia od nowa (nadpiszą wiersz w separations.csv)
        force_redo = Set{String}(["J0812-3905", "J0952-3839", "J1034-3224"])

        # pulsary już przetworzone (są w separations.csv)
        done_set = Set{String}()
        if isfile(sep_file)
            for line in eachline(sep_file)
                parts = split(strip(line), ",")
                length(parts) >= 1 && push!(done_set, strip(parts[1]))
            end
        end

        todo_lines = readlines(todo_file)[2:end]  # skip header
        half_n = cld(length(todo_lines), 2)
        todo_lines = todo_lines[1:half_n]
        println("Przerabiam pierwszą połowę listy: $half_n pulsarów.")

        skipped = String[]
        errored = String[]
        processed = String[]

        # awaryjna kopia wyników w pamięci — na wypadek gdyby zapis do separations.csv zawiódł
        fallback_rows = String[]

        for line in todo_lines
            parts = split(strip(line), ",")
            length(parts) < 2 && continue
            psr_name    = strip(parts[1])
            n_comp_default = parse(Int, strip(parts[2]))

            if psr_name in done_set && !(psr_name in force_redo)
                println("SKIP $psr_name — już przetworzone")
                push!(skipped, psr_name)
                continue
            end

            indir = vpmout * psr_name * "_16"
            if !isdir(indir)
                println("SKIP $psr_name — brak folderu $indir")
                push!(skipped, psr_name)
                continue
            end

            # --- pokaż uśredniony profil i zapytaj o liczbę komponentów ---
            n_comp = n_comp_default
            try
                p_json   = Tools.read_params(joinpath(indir, "params.json"))
                low_file = joinpath(indir, "pulsar_low.debase.p3fold_norefine")
                l        = Data.load_ascii(low_file)
                nl       = Data.normalize_per_pulse(l)
                avg_prof = vec(sum(nl, dims=1))
                bst = something(get(p_json, "bin_st", nothing), 1)
                bnd = something(get(p_json, "bin_end", nothing), length(avg_prof))

                figure(figsize=(8, 4))
                plot(bst:bnd, avg_prof[bst:bnd], color="steelblue", lw=1.5)
                title("$psr_name — ile komponentów widać?")
                xlabel("bin")
                tight_layout()
                show()

                print("Ile komponentów? (Enter = $n_comp_default): ")
                inp = strip(readline(stdin; keep=false))
                close("all")
                if !isempty(inp)
                    n_comp = parse(Int, inp)
                end
            catch e
                println("Nie udało się pokazać podglądu ($e), używam n_comp=$n_comp_default")
            end

            println("\n=== $psr_name (n_comp=$n_comp) ===")
            try
                res = Data.analyse_p3folds_16_new(indir, "norefine", n_comp=n_comp, save_summary_dir=summary_dir)
                if !isnothing(res) && !isnothing(res.row)
                    push!(fallback_rows, res.row)
                end
                push!(processed, psr_name)
            catch e
                println("ERROR $psr_name: $e")
                println("Pomijam i idę dalej.")
                push!(errored, psr_name)
            end
        end

        println("\n\n===== PODSUMOWANIE =====")
        println("Przetworzone ($(length(processed))): ", join(processed, ", "))
        println("Pominięte ($(length(skipped))): ", join(skipped, ", "))
        println("Błędy ($(length(errored))): ", join(errored, ", "))

        # awaryjny wydruk CSV — na wypadek gdyby zapis do separations.csv gdzieś zawiódł
        if !isempty(fallback_rows)
            println("\n===== AWARYJNA KOPIA separations.csv (skopiuj ręcznie w razie potrzeby) =====")
            println("Nazwa,ncomp,lon 1,lon 1 err,lon 2,lon 2 err,lon 3,lon 3 err,lon 4,lon 4 err,sep,sep err")
            for r in fallback_rows
                println(r)
            end
            println("===== KONIEC KOPII =====")
        end
    end

end # module

SpaTs.main()

println("Bye")
