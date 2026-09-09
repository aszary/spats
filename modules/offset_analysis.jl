module OffsetAnalysis

using Statistics
using Printf
using PyPlot
PyPlot.matplotlib.use("Tkagg")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
const I_NS    = 1e45          # moment of inertia [g cm^2]
const YR_S    = 3.156e7       # seconds per year
const B_CONST = 3.2e19        # B = B_CONST * sqrt(P * Pdot) [G]

# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

"""
Read psrcat.db and return a Dict: name => NamedTuple of all available params.
Derives tau_c, B, Edot from P0/P1.
"""
function read_psrcat_full(filename)
    isfile(filename) || error("psrcat.db not found: $filename")

    out   = Dict{String, Dict{String,Float64}}()
    rec   = Dict{String,Float64}()
    jname = ""
    bname = ""

    scalar_keys = Set(["P0","P1","F0","F1","DM","W50","S1400","DIST_DM","DIST_DM1","RM"])

    function flush_rec!()
        name = isempty(jname) ? bname : jname
        isempty(name) && return
        r = copy(rec)
        # derive P0 from F0 if missing
        if !haskey(r, "P0") && haskey(r, "F0") && r["F0"] > 0
            r["P0"] = 1.0 / r["F0"]
        end
        # derive P1 from F1 if missing
        if !haskey(r, "P1") && haskey(r, "F1") && haskey(r, "P0")
            r["P1"] = -r["F1"] / r["F0"]^2
        end
        # compute derived quantities
        p  = get(r, "P0", NaN)
        pd = get(r, "P1", NaN)
        if isfinite(p) && isfinite(pd) && p > 0 && pd > 0
            r["TAU_C"]  = p / (2 * pd) / YR_S          # characteristic age [yr]
            r["B"]      = B_CONST * sqrt(p * pd)        # surface B field [G]
            r["EDOT"]   = 4π^2 * I_NS * pd / p^3       # spin-down luminosity [erg/s]
            r["LOG_P"]  = log10(p)
            r["LOG_PD"] = log10(pd)
            r["LOG_TAU"] = log10(r["TAU_C"])
            r["LOG_B"]   = log10(r["B"])
            r["LOG_EDOT"] = log10(r["EDOT"])
        end
        if haskey(r, "DM")
            r["LOG_DM"] = r["DM"] > 0 ? log10(r["DM"]) : NaN
        end
        if haskey(r, "W50")
            r["LOG_W50"] = r["W50"] > 0 ? log10(r["W50"]) : NaN
        end
        if haskey(r, "S1400")
            r["LOG_S1400"] = r["S1400"] > 0 ? log10(r["S1400"]) : NaN
        end
        out[name] = r
    end

    for line in eachline(filename)
        if startswith(line, "@-")
            flush_rec!()
            empty!(rec)
            jname = ""
            bname = ""
            continue
        end
        (isempty(strip(line)) || startswith(line, "#")) && continue
        parts = split(line)
        length(parts) < 2 && continue
        key, val = parts[1], parts[2]
        if key == "PSRJ"
            isempty(jname) && (jname = val)
        elseif key == "PSRB"
            isempty(bname) && (bname = val)
        elseif key in scalar_keys && !haskey(rec, key)
            v = tryparse(Float64, val)
            isnothing(v) || (rec[key] = v)
        end
    end
    flush_rec!()
    return out
end

"""
Read offsets.csv. Returns Dict: name => (value, err, ncomp, grade).
For multi-component pulsars returns separation change (last - first offset),
for single-component pulsars returns the raw shift.
"""
function read_offsets(filename)
    out = Dict{String, NamedTuple{(:value,:err,:ncomp,:grade,:offsets,:errs),
                                   Tuple{Float64,Float64,Int,Int,Vector{Float64},Vector{Float64}}}}()
    isfile(filename) || (@warn "offsets.csv not found: $filename"; return out)
    for (nline, line) in enumerate(eachline(filename))
        nline == 1 && continue
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        f = split(s, ',')
        length(f) < 10 && continue
        name  = String(strip(f[1]))
        offs  = Float64[]
        errs  = Float64[]
        for c in 1:4
            v = tryparse(Float64, strip(f[2c]))
            e = tryparse(Float64, strip(f[2c+1]))
            (isnothing(v) || isnothing(e)) && continue
            push!(offs, v)
            push!(errs, e)
        end
        isempty(offs) && continue
        grade = something(tryparse(Int, strip(f[10])), 0)
        val, err = length(offs) > 1 ?
            (offs[end] - offs[1], sqrt(errs[end]^2 + errs[1]^2)) :
            (offs[1], errs[1])
        out[name] = (value=val, err=err, ncomp=length(offs), grade=grade,
                     offsets=offs, errs=errs)
    end
    return out
end

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

"""Spearman rank correlation coefficient."""
function spearman_r(x, y)
    n = length(x)
    rx = invperm(sortperm(x)) |> float
    ry = invperm(sortperm(y)) |> float
    return cor(rx, ry)
end

"""Two-tailed p-value for Spearman r given n samples."""
function spearman_pval(r, n)
    n <= 2 && return 1.0
    t  = r * sqrt((n - 2) / max(1 - r^2, 1e-15))
    # approximate: use normal for large n
    z  = t / sqrt(1 + t^2 / (n - 2))
    pv = 2 * (1 - _stdnorm_cdf(abs(z)))
    return clamp(pv, 0.0, 1.0)
end

function _stdnorm_cdf(z)
    # Abramowitz & Stegun approximation
    t = 1 / (1 + 0.2316419 * abs(z))
    c = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 +
            t * (-1.821255978 + t * 1.330274429))))
    p = 1 - (1/sqrt(2π)) * exp(-0.5*z^2) * c
    return z >= 0 ? p : 1 - p
end

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

"""
Analyse correlations between pulsar parameters (from psrcat.db) and the
measured high-low component offsets (from offsets.csv).

For each parameter a scatter plot is shown with:
  - data points colour-coded by number of components
  - Spearman rank correlation coefficient r_s and p-value
  - simple linear regression line

A final summary panel ranks all parameters by |r_s|.

Arguments:
  offset_file   path to offsets.csv
  catalogue     path to psrcat.db
  outdir        directory for saved figures
  min_grade     minimum quality grade to include (1-10)
  min_sigma     minimum offset significance (offset/err) to include
"""
function analyse_offset_correlations(;
        offset_file = normpath(joinpath(@__DIR__, "..", "input", "offsets.csv")),
        catalogue   = normpath(joinpath(@__DIR__, "..", "input", "psrcat.db")),
        outdir      = normpath(joinpath(@__DIR__, "..", "output")),
        min_grade   = 6,
        min_sigma   = 2.0)

    println("\n══════════════════════════════════════════════")
    println("  OFFSET vs PULSAR PARAMETER CORRELATION ANALYSIS")
    println("══════════════════════════════════════════════")

    # --- load data ---
    cat  = read_psrcat_full(catalogue)
    offs = read_offsets(offset_file)
    @printf("Catalogue: %d pulsars\n", length(cat))
    @printf("Offsets:   %d pulsars\n", length(offs))

    # --- filter offsets ---
    good = Dict(n => d for (n, d) in offs
                if d.grade >= min_grade && abs(d.value) / max(d.err, 1e-10) >= min_sigma)
    @printf("After grade≥%d and SNR≥%.1f: %d pulsars\n", min_grade, min_sigma, length(good))

    # --- match with catalogue ---
    matched_names = sort(collect(intersect(keys(good), keys(cat))))
    @printf("Matched with catalogue: %d pulsars\n\n", length(matched_names))
    isempty(matched_names) && (@warn "No matches — check pulsar name format"; return)

    # parameters to analyse (key in psrcat dict, label, log-scale for x-axis)
    params = [
        ("LOG_P",    "Period P [s] (log₁₀)",             true),
        ("LOG_PD",   "Period derivative Ṗ [s/s] (log₁₀)",true),
        ("LOG_TAU",  "Characteristic age τ_c [yr] (log₁₀)",true),
        ("LOG_B",    "Surface B field [G] (log₁₀)",       true),
        ("LOG_EDOT", "Spin-down luminosity Ė [erg/s] (log₁₀)",true),
        ("DM",       "Dispersion Measure DM [pc/cm³]",     false),
        ("LOG_DM",   "DM [pc/cm³] (log₁₀)",               true),
        ("W50",      "Pulse width W₅₀ [ms]",               false),
        ("LOG_W50",  "Pulse width W₅₀ [ms] (log₁₀)",      true),
        ("S1400",    "Flux S₁₄₀₀ [mJy]",                  false),
        ("DIST_DM",  "Distance (DM) [kpc]",                false),
    ]

    # collect all offset values
    colors_by_ncomp = Dict(1=>"#1976D2", 2=>"#E65100", 3=>"#388E3C", 4=>"#7B1FA2")

    summary = Tuple{String,Float64,Float64,Int}[]  # (label, r_s, pval, n)

    mkpath(outdir)

    for (key, label, _) in params
        # extract parameter values for matched pulsars
        xvals = Float64[]
        yvals = Float64[]
        nvals = Int[]
        nnames = String[]

        for n in matched_names
            v = get(cat[n], key, NaN)
            isfinite(v) || continue
            push!(xvals, v)
            push!(yvals, good[n].value)
            push!(nvals, good[n].ncomp)
            push!(nnames, n)
        end

        length(xvals) < 5 && continue

        rs   = spearman_r(xvals, yvals)
        pval = spearman_pval(rs, length(xvals))
        push!(summary, (label, rs, pval, length(xvals)))

        @printf("%-45s  r_s=%+.3f  p=%.3f  n=%d\n", label, rs, pval, length(xvals))

        # --- scatter plot ---
        figure(figsize=(7, 5))

        # group by ncomp for colour
        for nc in sort(unique(nvals))
            idx = nvals .== nc
            col = get(colors_by_ncomp, nc, "gray")
            scatter(xvals[idx], yvals[idx],
                    c=col, s=40, alpha=0.8, zorder=3,
                    label="$nc component$(nc>1 ? "s" : "")")
        end

        # linear regression line
        _plot_regression!(xvals, yvals)

        axhline(0.0, color="gray", lw=0.8, ls="--", zorder=1)
        xlabel(label, fontsize=10)
        ylabel("Offset high − low (°)", fontsize=10)
        sig_str = pval < 0.001 ? "p<0.001" : @sprintf("p=%.3f", pval)
        title(@sprintf("r_s = %+.3f  %s  (n=%d)", rs, sig_str, length(xvals)),
              fontsize=10)
        legend(fontsize=8, loc="best")
        tight_layout()

        safe_key = replace(key, "/" => "_", " " => "_")
        savefig(joinpath(outdir, "offset_corr_$(safe_key).pdf"))
        savefig(joinpath(outdir, "offset_corr_$(safe_key).png"), dpi=150)
        show()
        println("  Press Enter for next plot, 'q' to skip remaining...")
        inp = readline(stdin; keep=false)
        close("all")
        lowercase(strip(inp)) == "q" && break
    end

    # --- summary bar chart ---
    _plot_summary(summary, outdir)

    # --- component separation vs P ---
    _plot_separation_vs_params(matched_names, good, cat, outdir)

    println("\nDone. Figures saved to $outdir")
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _plot_regression!(x, y)
    isempty(x) && return
    mx, my = mean(x), mean(y)
    b = sum((x .- mx) .* (y .- my)) / max(sum((x .- mx).^2), 1e-30)
    a = my - b * mx
    xs = range(minimum(x), maximum(x), length=50)
    plot(collect(xs), a .+ b .* collect(xs),
         color="black", lw=1.2, ls="-", alpha=0.6, zorder=2, label="_nolegend_")
end

function _plot_summary(summary, outdir)
    isempty(summary) && return

    # sort by |r_s|
    sort!(summary, by=x -> abs(x[2]), rev=true)
    labels = [s[1] for s in summary]
    rs     = [s[2] for s in summary]
    pvals  = [s[3] for s in summary]

    colors = [p < 0.05 ? (r > 0 ? "#E53935" : "#1E88E5") :
              (r > 0 ? "#EF9A9A" : "#90CAF9")
              for (r, p) in zip(rs, pvals)]

    figure(figsize=(8, max(4, 0.45 * length(labels))))
    barh(1:length(labels), rs, color=colors, edgecolor="black", linewidth=0.5)
    axvline(0, color="black", lw=0.8)
    axvline( 0.3, color="gray", lw=0.7, ls="--", alpha=0.5)
    axvline(-0.3, color="gray", lw=0.7, ls="--", alpha=0.5)
    yticks(1:length(labels), labels, fontsize=8)
    xlabel("Spearman r_s", fontsize=10)
    title("Offset correlations — ranked by |r_s|\n(filled = p<0.05, red=positive, blue=negative)", fontsize=9)
    xlim(-1, 1)
    tight_layout()
    savefig(joinpath(outdir, "offset_corr_summary.pdf"))
    savefig(joinpath(outdir, "offset_corr_summary.png"), dpi=150)
    show()
    println("Summary plot — Press Enter to continue.")
    readline(stdin; keep=false)
    close("all")
end

"""
Extra plot: for multi-component pulsars, show individual component offsets
as a function of longitude, one panel per parameter bin.
"""
function _plot_separation_vs_params(matched_names, good, cat, outdir)
    # scatter: total separation change vs P (log scale)
    multi = [n for n in matched_names if good[n].ncomp >= 2]
    isempty(multi) && return

    figure(figsize=(7, 5))
    ps_all  = [get(get(cat, n, Dict()), "LOG_P",   NaN) for n in multi]
    sep_all = [good[n].value for n in multi]
    err_all = [good[n].err   for n in multi]
    tau_all = [get(get(cat, n, Dict()), "LOG_TAU", NaN) for n in multi]

    finite = isfinite.(ps_all) .& isfinite.(tau_all)
    ps_f, sep_f, err_f, tau_f = ps_all[finite], sep_all[finite], err_all[finite], tau_all[finite]

    isempty(ps_f) && (close("all"); return)

    # colour by log(tau_c)
    sc = scatter(ps_f, sep_f, c=tau_f,
                 cmap="plasma", s=50, alpha=0.85, zorder=3)
    errorbar(ps_f, sep_f, yerr=err_f,
             fmt="none", ecolor="gray", alpha=0.5, capsize=2, zorder=2)
    colorbar(sc, label="log₁₀(τ_c / yr)")
    axhline(0, color="gray", lw=0.8, ls="--")
    _plot_regression!(ps_f, sep_f)

    rs = spearman_r(ps_f, sep_f)
    pval = spearman_pval(rs, length(ps_f))
    sig_str = pval < 0.001 ? "p<0.001" : @sprintf("p=%.3f", pval)
    xlabel("log₁₀(P / s)", fontsize=10)
    ylabel("Separation change (°)\n[offset_last − offset_first]", fontsize=10)
    title(@sprintf("Multi-component pulsars: separation vs P\nr_s=%+.3f  %s  n=%d",
                   rs, sig_str, length(ps_f)), fontsize=10)
    tight_layout()
    savefig(joinpath(outdir, "offset_separation_vs_P.pdf"))
    savefig(joinpath(outdir, "offset_separation_vs_P.png"), dpi=150)
    show()
    println("Press Enter to close.")
    readline(stdin; keep=false)
    close("all")
end

end  # module OffsetAnalysis
