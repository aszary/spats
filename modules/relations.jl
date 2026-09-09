module Relations

using PyPlot
using Printf

# Pulsars renamed in the catalogue: J name translation mapping
const PSR_RENAMED = Dict("J1402-5124" => "J1402-5021")

struct PulsarRecord
    name::String
    params::Dict{String,Float64}
end

"""
Read pulsar names from a text file or CSV into a Set{String}.
"""
function read_pulsar_set(filename; is_csv=false)
    names = Set{String}()
    isfile(filename) || return names
    for (i, line) in enumerate(eachline(filename))
        is_csv && i == 1 && continue  # header
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        name = is_csv ? String(strip(split(s, ',')[1])) : String(first(split(s)))
        haskey(PSR_RENAMED, name) && (name = PSR_RENAMED[name])
        push!(names, name)
    end
    return names
end

"""
Read full ATNF psrcat database and derive all required parameter features.
"""
function read_psrcat_extended(catalogue)
    isfile(catalogue) || error("psrcat database not found: $catalogue")

    records = PulsarRecord[]
    rec = Dict{String,Float64}()
    jname, bname = "", ""

    yr_s = 3.15576e7
    inertia = 1e45

    function _flush_record(name, r)
        f0 = get(r, "F0", NaN)
        p0 = get(r, "P0", NaN)
        if !isfinite(f0) && isfinite(p0) && p0 > 0
            f0 = 1 / p0
        elseif !isfinite(p0) && isfinite(f0) && f0 > 0
            p0 = 1 / f0
        end
        p1 = get(r, "P1", NaN)
        if !isfinite(p1)
            f1 = get(r, "F1", NaN)
            isfinite(f1) && isfinite(f0) && f0 > 0 && (p1 = -f1 / f0^2)
        end

        (isfinite(p0) && isfinite(p1) && p0 > 0) || return

        dict = Dict{String,Float64}()
        dict["LOG_P"] = log10(p0)
        dict["P"] = p0

        if p1 > 0
            dict["LOG_PD"] = log10(p1)
            dict["PD"] = p1

            tau_c = p0 / (2 * p1 * yr_s)
            dict["TAU"] = tau_c
            dict["LOG_TAU"] = tau_c > 0 ? log10(tau_c) : NaN

            b_surf = 3.2e19 * sqrt(p0 * p1)
            dict["B"] = b_surf
            dict["LOG_B"] = b_surf > 0 ? log10(b_surf) : NaN

            edot = (4 * pi^2 * inertia * p1) / (p0^3)
            dict["EDOT"] = edot
            dict["LOG_EDOT"] = edot > 0 ? log10(edot) : NaN
        end

        dm = get(r, "DM", NaN)
        if isfinite(dm) && dm > 0
            dict["DM"] = dm
            dict["LOG_DM"] = log10(dm)
        end

        w50 = get(r, "W50", NaN)
        if isfinite(w50) && w50 > 0
            dict["W50"] = w50
            dict["LOG_W50"] = log10(w50)
        end

        s1400 = get(r, "S1400", NaN)
        if isfinite(s1400) && s1400 > 0
            dict["S1400"] = s1400
        end

        dist = get(r, "DIST_DM", get(r, "DIST", NaN))
        if isfinite(dist) && dist > 0
            dict["DIST_DM"] = dist
        end

        push!(records, PulsarRecord(name, dict))
    end

    for line in eachline(catalogue)
        if startswith(line, "@-")
            _flush_record(isempty(jname) ? bname : jname, rec)
            empty!(rec)
            jname, bname = "", ""
            continue
        end
        (isempty(strip(line)) || startswith(line, "#")) && continue
        parts = split(line)
        length(parts) < 2 && continue
        key, val = parts[1], parts[2]
        if key == "PSRJ" && isempty(jname)
            jname = val
        elseif key == "PSRB" && isempty(bname)
            bname = val
        elseif key in ("P0", "P1", "F0", "F1", "DM", "W50", "S1400", "DIST_DM", "DIST") && !haskey(rec, key)
            v = tryparse(Float64, val)
            isnothing(v) || (rec[key] = v)
        end
    end
    _flush_record(isempty(jname) ? bname : jname, rec)

    return records
end

const PARAM_INFO = Dict(
    "LOG_P"    => ("Period \$P\$ [s]", true),
    "LOG_PD"   => ("Period derivative \$\\dot{P}\$ [s s\$^{-1}\$]", true),
    "LOG_TAU"  => ("Characteristic age \$\\tau_c\$ [yr]", true),
    "LOG_B"    => ("Surface magnetic field \$B\$ [G]", true),
    "LOG_EDOT" => ("Spin-down luminosity \$\\dot{E}\$ [erg s\$^{-1}\$]", true),
    "DM"       => ("Dispersion Measure DM [pc cm\$^{-3}\$]", false),
    "LOG_DM"   => ("Dispersion Measure DM [pc cm\$^{-3}\$]", true),
    "W50"      => ("Pulse width \$W_{50}\$ [ms]", false),
    "LOG_W50"  => ("Pulse width \$W_{50}\$ [ms]", true),
    "S1400"    => ("Flux \$S_{1400}\$ [mJy]", false),
    "DIST_DM"  => ("Distance \$d_{\\text{DM}}\$ [kpc]", false),
)

"""
Plot scatter relation between two parameters `px` and `py`, differentiating
between Drifting, P3-only, and Background ATNF pulsars.
"""
function plot_relation(records, px::String, py::String, drifting_set, p3only_set, outdir;
                       name_mod="$(lowercase(px))_vs_$(lowercase(py))", show_=false)

    label_x, is_log_x = get(PARAM_INFO, px, (px, false))
    label_y, is_log_y = get(PARAM_INFO, py, (py, false))

    bg_x, bg_y = Float64[], Float64[]
    drift_x, drift_y = Float64[], Float64[]
    p3_x, p3_y = Float64[], Float64[]

    for rec in records
        haskey(rec.params, px) && haskey(rec.params, py) || continue
        vx = rec.params[px]
        vy = rec.params[py]
        (isfinite(vx) && isfinite(vy)) || continue

        if rec.name in p3only_set
            push!(p3_x, vx)
            push!(p3_y, vy)
        elseif rec.name in drifting_set
            push!(drift_x, vx)
            push!(drift_y, vy)
        else
            push!(bg_x, vx)
            push!(bg_y, vy)
        end
    end

    rc("font", size=10.)
    rc("axes", linewidth=0.7)
    rc("lines", linewidth=0.7)

    figure(figsize=(6.5, 5.2))
    ax = gca()

    if is_log_x
        ax.set_xscale("log")
        bg_x_plot = 10 .^ bg_x
        drift_x_plot = 10 .^ drift_x
        p3_x_plot = 10 .^ p3_x
    else
        bg_x_plot = bg_x
        drift_x_plot = drift_x
        p3_x_plot = p3_x
    end

    if is_log_y
        ax.set_yscale("log")
        bg_y_plot = 10 .^ bg_y
        drift_y_plot = 10 .^ drift_y
        p3_y_plot = 10 .^ p3_y
    else
        bg_y_plot = bg_y
        drift_y_plot = drift_y
        p3_y_plot = p3_y
    end

    plot(bg_x_plot, bg_y_plot, ".", ms=2.5, c="0.75", alpha=0.35, label="ATNF background", zorder=1)
    plot(drift_x_plot, drift_y_plot, "o", ms=4.5, c="tab:blue", alpha=0.85, mec="black", mew=0.4, label="Drifting", zorder=3)
    plot(p3_x_plot, p3_y_plot, "^", ms=5.0, c="tab:red", alpha=0.85, mec="black", mew=0.4, label="P3-only", zorder=4)

    xlabel(label_x)
    ylabel(label_y)
    minorticks_on()
    legend(fontsize=8, loc="best", framealpha=0.9)
    grid(true, which="both", ls=":", lw=0.4, alpha=0.5)

    savepath = joinpath(outdir, "relation_$(name_mod).pdf")
    savefig(savepath)
    savefig(replace(savepath, ".pdf" => ".png"))
    println("Saved relation plot: $savepath (Drifting: $(length(drift_x)), P3-only: $(length(p3_x)), BG: $(length(bg_x)))")

    if show_
        PyPlot.show()
        println("Press Enter to close the figure.")
        readline(stdin; keep=false)
    end
    close()
end

"""
Generates key parameter relation figures.
"""
function plot_all_relations(outdir;
                            catalogue=normpath(joinpath(@__DIR__, "..", "input", "psrcat.db")),
                            offsets=normpath(joinpath(@__DIR__, "..", "input", "offsets.csv")),
                            offsets_p3only=normpath(joinpath(@__DIR__, "..", "input", "offsets_p3only.csv")),
                            drift_list=normpath(joinpath(@__DIR__, "..", "input", "drift_pulsars_P3.txt")),
                            show_=false)

    mkpath(outdir)
    records = read_psrcat_extended(catalogue)

    drifting_set = read_pulsar_set(offsets, is_csv=true)
    union!(drifting_set, read_pulsar_set(drift_list, is_csv=false))

    p3only_set = read_pulsar_set(offsets_p3only, is_csv=true)

    pairs = [
        ("LOG_P", "LOG_PD"),
        ("LOG_P", "LOG_W50"),
        ("LOG_EDOT", "LOG_W50"),
        ("LOG_TAU", "LOG_W50"),
        ("LOG_B", "LOG_W50"),
        ("LOG_P", "DM"),
        ("LOG_EDOT", "S1400"),
        ("LOG_TAU", "LOG_EDOT")
    ]

    for (px, py) in pairs
        plot_relation(records, px, py, drifting_set, p3only_set, outdir; show_=show_)
    end
end

end # module Relations
