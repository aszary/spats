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

    # parameters to analyse: (key, short label, use_log)
    params = [
        ("LOG_P",    "Period P [s] (log₁₀)",                  true),
        ("LOG_PD",   "Period derivative Ṗ [s/s] (log₁₀)",     true),
        ("LOG_TAU",  "Characteristic age τ_c [yr] (log₁₀)",   true),
        ("LOG_B",    "Surface B field [G] (log₁₀)",            true),
        ("LOG_EDOT", "Spin-down luminosity Ė [erg/s] (log₁₀)", true),
        ("DM",       "Dispersion Measure DM [pc/cm³]",          false),
        ("LOG_DM",   "DM [pc/cm³] (log₁₀)",                   true),
        ("W50",      "Pulse width W₅₀ [ms]",                   false),
        ("LOG_W50",  "Pulse width W₅₀ [ms] (log₁₀)",          true),
        ("S1400",    "Flux S₁₄₀₀ [mJy]",                      false),
        ("DIST_DM",  "Distance (DM) [kpc]",                    false),
    ]

    # human-readable interpretation hint for each parameter
    param_hints = Dict(
        "LOG_P"    => "Dłuższy okres → pulsar wolniej się obraca, inny tryb emisji?\nKorelacja sugeruje czy aberracja/retardacja zależy od P.",
        "LOG_PD"   => "Wyższe Ṗ → silniejszy spin-down, młodszy pulsar.\nKorelacja z Ṗ może wskazywać na związek z energią emisji.",
        "LOG_TAU"  => "Starsze pulsary (duże τ_c) vs młode.\nCzy offset zmienia się w czasie życia pulsara?",
        "LOG_B"    => "Silniejsze pole B → inny mechanizm emisji?\nModel A/R przewiduje zależność od geometrii magnetosfery.",
        "LOG_EDOT" => "Pulsary tracące więcej energii często mają inny kształt profilu.\nKorelacja może wskazywać na wpływ wiatru pulsarowego.",
        "DM"       => "DM ≈ całka z gęstości elektronów wzdłuż linii widzenia.\nBezpośredni wpływ na offset jest mało prawdopodobny,\nale DM koreluje z odległością i ze strumieniem.",
        "LOG_DM"   => "To samo co DM ale w skali log — lepiej widać duży zakres wartości.",
        "W50"      => "Szerszy impuls → bardziej rozległa emisja lub geometria pod większym kątem.\nMoże korelować z offset jeśli geometria decyduje o przesunięciu.",
        "LOG_W50"  => "To samo co W₅₀ ale w skali log.",
        "S1400"    => "Strumień na 1400 MHz — jaśniejsze pulsary są bliżej lub silniej emitują.\nSłaba korelacja z offsetem byłaby zaskoczeniem (to powinien być efekt geometryczny).",
        "DIST_DM"  => "Odległość oszacowana z DM.\nBezpośrednio nie powinna wpływać na offset — to test na systematykę.",
    )

    # collect all offset values
    colors_by_ncomp = Dict(1=>"#1976D2", 2=>"#E65100", 3=>"#388E3C", 4=>"#7B1FA2")

    summary        = Tuple{String,Float64,Float64,Int}[]       # (label, r_s, pval, n)
    summary_by_nc  = Tuple{String,Int,Float64,Float64,Int}[]  # (label, nc, r_s, pval, n)

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

        # --- overall correlation ---
        rs   = spearman_r(xvals, yvals)
        pval = spearman_pval(rs, length(xvals))
        push!(summary, (label, rs, pval, length(xvals)))

        # --- per-ncomp correlations for summary ---
        for (nc, idx_nc) in [(1, nvals.==1), (2, nvals.==2), (99, nvals.>=3)]
            xg, yg = xvals[idx_nc], yvals[idx_nc]
            length(xg) < 4 && continue
            rsg = spearman_r(xg, yg)
            pvg = spearman_pval(rsg, length(xg))
            push!(summary_by_nc, (label, nc, rsg, pvg, length(xg)))
        end

        # --- print ---
        @printf("%-45s  r_s=%+.3f  p=%.3f  n=%d\n", label, rs, pval, length(xvals))

        # --- scatter: fixed 2×2 grid, groups: all / 1c / 2c / 3+c ---
        idx1  = nvals .== 1
        idx2  = nvals .== 2
        idx3p = nvals .>= 3
        panel_groups = [
            (xvals,       yvals,       "Wszystkie",    "#555555"),
            (xvals[idx1], yvals[idx1], "1 komponent",  colors_by_ncomp[1]),
            (xvals[idx2], yvals[idx2], "2 komponenty", colors_by_ncomp[2]),
            (xvals[idx3p],yvals[idx3p],"3+ komponenty",colors_by_ncomp[3]),
        ]

        figure(figsize=(10, 7))
        for (pi, (xg, yg, glabel, gcol)) in enumerate(panel_groups)
            subplot(2, 2, pi)
            if length(xg) >= 2
                scatter(xg, yg, c=gcol, s=35, alpha=0.8, zorder=3)
                _plot_best_fit_ax!(gca(), xg, yg)
            end
            axhline(0, color="gray", lw=0.8, ls="--")
            xlabel(label, fontsize=8)
            ylabel("Offset (°)", fontsize=8)
            if length(xg) >= 4
                rsg = spearman_r(xg, yg)
                pvg = spearman_pval(rsg, length(xg))
                sstr = pvg < 0.001 ? "p<0.001" : @sprintf("p=%.3f", pvg)
                istr = pvg < 0.05 ? " ★" : ""
                title("$glabel  (n=$(length(xg)))\nr_s=$(round(rsg,digits=3))  $sstr$istr",
                      fontsize=8)
            else
                title("$glabel  (n=$(length(xg)))\nza mało danych", fontsize=8)
            end
        end
        suptitle(label, fontsize=10, fontweight="bold")
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

    # --- 2×2 summary: all / 1-comp / 2-comp / 3+comp ---
    _plot_summary_grid(summary, summary_by_nc, outdir)

    # --- component separation vs P ---
    _plot_separation_vs_params(matched_names, good, cat, outdir)

    println("\nDone. Figures saved to $outdir")
    _print_legend()
end

function _print_legend()
    println("""
\n╔══════════════════════════════════════════════════════════════════════╗
║              JAK CZYTAĆ WYKRESY KORELACJI                           ║
╠══════════════════════════════════════════════════════════════════════╣
║ OŚ Y (na każdym wykresie):                                          ║
║   Offset high−low [°] = o ile stopni profil pulsara przesuwa się    ║
║   między 1023 MHz a 1523 MHz.                                        ║
║   >0 → profil przy wyższej częstotliwości jest przesunięty w prawo  ║
║   <0 → profil przy wyższej częstotliwości jest przesunięty w lewo   ║
║   Dla pulsarów z 2+ składowymi: zmiana separacji (ostatnia−pierwsza)║
╠══════════════════════════════════════════════════════════════════════╣
║ PARAMETRY (oś X):                                                   ║
║                                                                     ║
║  P (okres)     Jak szybko obraca się pulsar. Dłuższy okres = starszy║
║                lub słabiej wyhamowany. Model A/R: offset ~ 1/P²     ║
║                → spodziewana korelacja ujemna z P.                  ║
║                                                                     ║
║  Ṗ (dP/dt)    Jak szybko pulsar traci energię obrotową.             ║
║                Ṗ duże = pulsar młody i energetyczny.                ║
║                                                                     ║
║  τ_c = P/(2Ṗ) Wiek charakterystyczny [lata]. Duże τ_c = stary.     ║
║                Jeśli korelacja z τ_c → offset zmienia się z wiekiem.║
║                                                                     ║
║  B = 3.2e19√(PṖ)  Pole magnetyczne powierzchni [Gauss].            ║
║                Duże B → silniejsza magnetosfera, inna geometria.    ║
║                                                                     ║
║  Ė = 4π²IṖ/P³ Świecistość spin-down [erg/s]. Proxy energii emisji. ║
║                                                                     ║
║  DM            Miara dyspersji [pc/cm³] = całka z gęstości e⁻.     ║
║                Koreluje z odległością. Nie powinna wpływać na offset.║
║                Korelacja byłaby artefaktem.                         ║
║                                                                     ║
║  W₅₀           Szerokość impulsu przy 50% maksimum [ms].           ║
║                Szerszy profil → emisja z większego obszaru albo     ║
║                obserwacja pod dużym kątem do osi magnetycznej.      ║
║                                                                     ║
║  S₁₄₀₀         Strumień radiowy [mJy]. Proxy jasności / odległości. ║
║                Korelacja z offsetem byłaby podejrzana (selekcja?).  ║
║                                                                     ║
║  Odległość     Z DM. Nie powinna korelować — test systematyki.      ║
╠══════════════════════════════════════════════════════════════════════╣
║ STATYSTYKA:                                                         ║
║  r_s ∈ [-1,1]  Korelacja Spearmana (rangowa, odporna na outliery).  ║
║  |r_s| > 0.5  → silna korelacja                                     ║
║  |r_s| 0.3-0.5 → umiarkowana                                        ║
║  |r_s| < 0.3  → słaba lub brak                                      ║
║  p < 0.05     → korelacja statystycznie ISTOTNA (mało prawdopodobna ║
║                 przy braku zależności)                               ║
║  p > 0.05     → nie możemy odrzucić hipotezy że to przypadek        ║
╠══════════════════════════════════════════════════════════════════════╣
║ MODELE DOPASOWANIA (czarna linia = najlepszy wg AIC):               ║
║  linear     y = a + b·x          → prosta zależność liniowa         ║
║  quadratic  y = a + b·x + c·x²  → minimum/maksimum w środku zakresu║
║  log        y = a + b·log(x)    → efekt nasycenia (szybki wzrost,   ║
║                                    potem plateau)                   ║
║  power      |y| ~ x^b           → potęgowa (jak w astronomii B~P^α) ║
║  AIC mniejszy = lepszy model (uwzględnia liczbę parametrów)         ║
╚══════════════════════════════════════════════════════════════════════╝
""")
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
Fit multiple models to (x, y) data, pick the best by AIC, plot it.
Returns (best_model_name, AIC_best).
Models tried:
  - linear:      y = a + b*x
  - quadratic:   y = a + b*x + c*x²
  - logarithmic: y = a + b*log(x)   (only when all x > 0)
  - power law:   log|y| = a + b*log(x)  (only when all x>0 and all |y|>0)
"""
function _plot_best_fit!(x, y)
    isempty(x) && return ("none", Inf)
    n = length(x)
    xs_plot = collect(range(minimum(x), maximum(x), length=200))

    # AICc = corrected AIC for small samples: penalises extra params harder
    function aicc(rss, k)
        rss <= 0 && return Inf
        aic_val = n * log(rss / n) + 2 * k
        # small-sample correction term
        denom = n - k - 1
        denom <= 0 && return Inf
        aic_val + 2 * k * (k + 1) / denom
    end

    function ols_rss_ypred(A)
        coef   = A \ y
        ypred  = A * coef
        rss    = sum((y .- ypred).^2)
        return coef, rss, ypred
    end

    results = Tuple{String, Float64, Vector{Float64}, Function}[]

    # --- linear ---
    A_lin = hcat(ones(n), x)
    c_lin, rss_lin, _ = ols_rss_ypred(A_lin)
    push!(results, ("linear", aicc(rss_lin, 2), c_lin,
          xs -> c_lin[1] .+ c_lin[2] .* xs))

    # --- quadratic ---
    A_qua = hcat(ones(n), x, x.^2)
    c_qua, rss_qua, _ = ols_rss_ypred(A_qua)
    push!(results, ("quadratic", aicc(rss_qua, 3), c_qua,
          xs -> c_qua[1] .+ c_qua[2] .* xs .+ c_qua[3] .* xs.^2))

    # --- logarithmic (needs all x > 0) ---
    if all(x .> 0)
        lx = log.(x)
        A_log = hcat(ones(n), lx)
        c_log, rss_log, _ = ols_rss_ypred(A_log)
        push!(results, ("log", aicc(rss_log, 2), c_log,
              xs -> c_log[1] .+ c_log[2] .* log.(max.(xs, 1e-300))))
    end

    # --- power law: fit in log-log, but compute RSS back in y-space ---
    if all(x .> 0) && all(abs.(y) .> 0)
        signs = sign.(y)
        lx    = log.(x)
        ly    = log.(abs.(y))
        A_pw  = hcat(ones(n), lx)
        c_pw  = A_pw \ ly
        ypred_pw = signs .* exp.(A_pw * c_pw)   # back to y-space
        rss_pw   = sum((y .- ypred_pw).^2)       # RSS in y-space
        push!(results, ("power", aicc(rss_pw, 2), c_pw,
              xs -> signs[1] .* exp.(c_pw[1] .+ c_pw[2] .* log.(max.(xs, 1e-300)))))
    end

    # pick best (lowest AIC)
    best = argmin([r[2] for r in results])
    bname = results[best][1]
    baic  = results[best][2]

    # plot all models faintly, best prominently
    style_map = Dict("linear"=>"--", "quadratic"=>"-.", "log"=>":", "power"=>"--")
    for (i, (mname, _, _, mfun)) in enumerate(results)
        ys = mfun(xs_plot)
        lw = i == best ? 1.8 : 0.8
        al = i == best ? 0.85 : 0.3
        col = i == best ? "black" : "gray"
        lab = i == best ? "$mname (best, AIC=$(round(Int,baic)))" : "_nolegend_"
        plot(xs_plot, ys, color=col, lw=lw, ls=get(style_map, mname, "-"),
             alpha=al, zorder=2, label=lab)
    end

    return bname, baic
end

"""Combined 2×2 summary: all pulsars + 1-comp + 2-comp + 3+comp."""
function _plot_summary_grid(summary, summary_by_nc, outdir)
    isempty(summary) && return

    panels = [
        ("Wszystkie",    [(s[1],s[2],s[3]) for s in summary]),
        ("1 komponent",  [(s[1],s[3],s[4]) for s in summary_by_nc if s[2]==1]),
        ("2 komponenty", [(s[1],s[3],s[4]) for s in summary_by_nc if s[2]==2]),
        ("3+ komponenty",[(s[1],s[3],s[4]) for s in summary_by_nc if s[2]==99]),
    ]

    # short labels for y-axis (strip units/log info)
    short_label(l) = replace(l, r" \(log₁₀\)" => " (log)",
                                r" \[.*?\]"    => "",
                                r"Dispersion Measure " => "")

    figure(figsize=(10, 7))
    for (pi, (ptitle, rows)) in enumerate(panels)
        subplot(2, 2, pi)
        isempty(rows) && (title("$ptitle\n(brak danych)", fontsize=8); continue)

        sort!(rows, by=r -> abs(r[2]), rev=true)
        labs  = [short_label(r[1]) for r in rows]
        rs    = [r[2] for r in rows]
        pvs   = [r[3] for r in rows]
        n_bar = length(labs)

        bar_colors = [isnan(p) ? "#cccccc" :
                      p < 0.05 ? (r > 0 ? "#E53935" : "#1E88E5") :
                                  (r > 0 ? "#FFCDD2" : "#BBDEFB")
                      for (r, p) in zip(rs, pvs)]

        ax = gca()
        ax.barh(1:n_bar, rs, color=bar_colors, edgecolor="black", linewidth=0.4)
        ax.axvline(0,    color="black", lw=0.8)
        ax.axvline( 0.3, color="gray",  lw=0.6, ls="--", alpha=0.5)
        ax.axvline(-0.3, color="gray",  lw=0.6, ls="--", alpha=0.5)
        ax.set_yticks(1:n_bar)
        ax.set_yticklabels(labs, fontsize=7)
        ax.set_xlim(-1, 1)
        ax.set_xlabel("Spearman r_s", fontsize=8)
        ax.set_title("$ptitle\n(czerwony/niebieski = p<0.05, jasny = nieistotne)", fontsize=8)
    end

    tight_layout()
    savefig(joinpath(outdir, "offset_corr_summary.pdf"))
    savefig(joinpath(outdir, "offset_corr_summary.png"), dpi=150)
    show()
    println("Summary — Press Enter to continue.")
    readline(stdin; keep=false)
    close("all")
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
    _plot_best_fit!(ps_f, sep_f)

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

"""Simple linear regression drawn on given Axes."""
function _linear_fit_ax!(ax, x, y)
    length(x) < 2 && return
    mx, my = mean(x), mean(y)
    b = sum((x .- mx) .* (y .- my)) / max(sum((x .- mx).^2), 1e-30)
    a = my - b * mx
    xs = collect(range(minimum(x), maximum(x), length=100))
    ax.plot(xs, a .+ b .* xs, color="black", lw=1.4, ls="--", alpha=0.7, zorder=2)
end

"""Best-fit model selection drawn on given matplotlib Axes object."""
function _plot_best_fit_ax!(ax, x, y)
    isempty(x) && return
    n = length(x)
    n < 3 && return

    xs_plot = collect(range(minimum(x), maximum(x), length=200))

    function aicc(rss, k)
        rss <= 0 && return Inf
        d = n - k - 1
        d <= 0 && return Inf
        n * log(rss / n) + 2 * k + 2 * k * (k + 1) / d
    end

    ols_c(A) = (c = A \ y; (c, sum((y .- A*c).^2)))

    cands = Tuple{String, Float64, Function}[]

    A_lin = hcat(ones(n), x)
    c, r = ols_c(A_lin)
    let cv = copy(c); push!(cands, ("linear", aicc(r,2), xs -> cv[1] .+ cv[2].*xs)); end

    A_qua = hcat(ones(n), x, x.^2)
    c, r = ols_c(A_qua)
    let cv = copy(c); push!(cands, ("quadratic", aicc(r,3), xs -> cv[1] .+ cv[2].*xs .+ cv[3].*xs.^2)); end

    if all(x .> 0)
        lx = log.(x); A_log = hcat(ones(n), lx)
        c, r = ols_c(A_log)
        let cv = copy(c); push!(cands, ("log", aicc(r,2), xs -> cv[1] .+ cv[2].*log.(max.(xs,1e-300)))); end
    end

    if all(x .> 0) && all(abs.(y) .> 0)
        signs = sign.(y); lx = log.(x); A_pw = hcat(ones(n), lx)
        c_pw = A_pw \ log.(abs.(y))
        ypred_pw = signs .* exp.(A_pw * c_pw)
        r_pw = sum((y .- ypred_pw).^2)
        let cpw = copy(c_pw), sg = signs[1]
            push!(cands, ("power", aicc(r_pw,2),
                  xs -> sg .* exp.(cpw[1] .+ cpw[2].*log.(max.(xs,1e-300)))))
        end
    end

    best = argmin([c[2] for c in cands])
    bname, _, bfun = cands[best]
    ys = bfun(xs_plot)

    # clip plotted line to data y-range (±20% padding) so bad extrapolation stays invisible
    ylo = minimum(y) - 0.2 * (maximum(y) - minimum(y))
    yhi = maximum(y) + 0.2 * (maximum(y) - minimum(y))
    ys_clipped = clamp.(ys, ylo, yhi)

    ax.plot(xs_plot, ys_clipped, color="black", lw=1.5, ls="--", alpha=0.7,
            label="$bname (best)", zorder=2)
end

"""Summary bar chart of correlations broken down by number of components."""
function _plot_summary_by_nc(summary_by_nc, outdir)
    isempty(summary_by_nc) && return

    all_nc = sort(unique([s[2] for s in summary_by_nc]))
    all_labels = unique([s[1] for s in summary_by_nc])

    nL = length(all_labels)
    nC = length(all_nc)
    nc_label = Dict(1=>"1 komponent", 2=>"2 komponenty", 99=>"3+ komponenty")

    # stacked rows so all groups fit on screen
    figure(figsize=(9, 3.5 * nC))
    for (j, nc) in enumerate(all_nc)
        subplot(nC, 1, j)
        rows = [(s[1], s[3], s[4]) for s in summary_by_nc if s[2] == nc]
        sort!(rows, by=r -> abs(r[2]), rev=true)
        labs = [r[1] for r in rows]
        rs   = [r[2] for r in rows]
        pvs  = [r[3] for r in rows]
        cols = [p < 0.05 ? (r > 0 ? "#E53935" : "#1E88E5") :
                            (r > 0 ? "#EF9A9A" : "#90CAF9")
                for (r,p) in zip(rs, pvs)]
        barh(1:length(labs), rs, color=cols, edgecolor="black", linewidth=0.4)
        axvline(0,    color="black", lw=0.8)
        axvline( 0.3, color="gray",  lw=0.6, ls="--", alpha=0.5)
        axvline(-0.3, color="gray",  lw=0.6, ls="--", alpha=0.5)
        yticks(1:length(labs), labs, fontsize=7)
        xlim(-1, 1)
        xlabel("Spearman r_s", fontsize=8)
        title("$(get(nc_label, nc, "$nc komp."))   (czerwony=p<0.05 dodatni, niebieski=p<0.05 ujemny)", fontsize=8)
    end

    tight_layout()
    savefig(joinpath(outdir, "offset_corr_summary_by_ncomp.pdf"))
    savefig(joinpath(outdir, "offset_corr_summary_by_ncomp.png"), dpi=150)
    show()
    println("Summary per ncomp — Press Enter to continue.")
    readline(stdin; keep=false)
    close("all")
end

end  # module OffsetAnalysis
