module LrfsBatch

using DelimitedFiles
using JSON
using PyPlot
using Printf

# Try to safely resolve the Data module (either from parent SpaTs or standalone)
try
    import ..Data
catch
    try
        import .Data
    catch
        try
            import Main.Data
        catch
            try
                include("data.jl")
                import .Data
            catch
                # Data module not available
            end
        end
    end
end

include("LrfsDiagnostics.jl")
using .LrfsDiagnostics

export batch_analyze_lrfs, plot_lrfs_results

"""
    _read_debase_ascii(debase_file::String) -> Matrix{Float64}

Fast, non-allocating streaming reader for PSRCHIVE ASCII debase files.
Avoids allocating millions of temporary SubString vectors.
"""
function _read_debase_ascii(debase_file::String)
    open(debase_file, "r") do io
        line1 = readline(io)
        isempty(line1) && return zeros(Float64, 0, 0)
        h = split(line1)
        length(h) < 12 && return zeros(Float64, 0, 0)
        n_pulses = parse(Int, h[6])
        n_bins   = parse(Int, h[12])
        
        X = zeros(Float64, n_pulses, n_bins)
        for line in eachline(io)
            itr = eachsplit(line)
            item1 = iterate(itr)
            isnothing(item1) && continue
            item2 = iterate(itr, item1[2])
            isnothing(item2) && continue
            item3 = iterate(itr, item2[2])
            isnothing(item3) && continue
            item4 = iterate(itr, item3[2])
            isnothing(item4) && continue
            
            p = parse(Int, item1[1]) + 1
            b = parse(Int, item3[1]) + 1
            v = parse(Float64, item4[1])
            if 1 <= p <= n_pulses && 1 <= b <= n_bins
                @inbounds X[p, b] = v
            end
        end
        return X
    end
end

"""
    batch_analyze_lrfs(vpmout::String, list_file::String, out_csv::String="lrfs_classifications.csv"; force_debase::Bool=false)

Reads a list of pulsars from `list_file`, finds their output directories in `vpmout`,
loads full-range single-pulse data (generating `pulsar.debase.txt` directly from `pulsar.spCf16`
if necessary, ignoring high/low frequency splits), runs LRFS phase tracking,
and saves the results to `out_csv`.
"""
function batch_analyze_lrfs(vpmout::String, list_file::String, out_csv::String="lrfs_classifications.csv"; force_debase::Bool=false)
    if !isfile(list_file)
        error("Pulsar list file not found: $list_file")
    end

    # Parse both pulsar names and known P3 values from the list file.
    # Expected format per line: "JNAME P3value(error)"  e.g. "J0601-0527 2.041(7)"
    # The error suffix (e.g. "(7)") is stripped before parsing.
    names   = String[]
    p3_dict = Dict{String, Float64}()   # name => P3 in pulses (nothing if missing/unparseable)

    for line in eachline(list_file)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        parts = split(s)
        psr_name = String(parts[1])
        push!(names, psr_name)
        if length(parts) >= 2
            # Strip parenthesised error: "2.041(7)" → "2.041"
            raw_p3 = replace(String(parts[2]), r"\(.*\)" => "")
            p3_parsed = tryparse(Float64, raw_p3)
            if !isnothing(p3_parsed) && p3_parsed > 1.0
                p3_dict[psr_name] = p3_parsed
            end
        end
    end

    println("Found $(length(names)) pulsars in list ($(length(p3_dict)) with known P3).")

    results = []

    for (i, name) in enumerate(names)
        known_p3 = get(p3_dict, name, nothing)   # Float64 or nothing
        try
            # Directory resolution fallback (check all 4 possible variations)
            candidates = [
                joinpath(vpmout, name * "_16"),
                vpmout * name * "_16",
                joinpath(vpmout, name),
                vpmout * name
            ]
            outdir = ""
            for cand in candidates
                if isdir(cand)
                    outdir = cand
                    break
                end
            end
            
            if isempty(outdir)
                @warn "[$i/$(length(names))] Skipping $name — Output dir not found."
                continue
            end
            
            params_file = joinpath(outdir, "params.json")
            
            # Check for full-range debase file, or create it directly from pulsar.spCf16
            # Note: We deliberately exclude pulsar_high_debase.txt and pulsar_low_debase.txt
            # to ensure analysis covers the full frequency range rather than split sub-bands.
            debase_target = joinpath(outdir, "pulsar.debase.txt")
            debase_file = ""

            # Check if pulsar.spCf16 exists
            spcf16_file = ""
            for cand in ["pulsar.spCf16", "pulsar.spCF16", "pulsar.spcf16"]
                cand_path = joinpath(outdir, cand)
                if isfile(cand_path)
                    spcf16_file = cand_path
                    break
                end
            end

            if !isempty(spcf16_file) && isfile(params_file)
                # If pulsar.debase.txt does not exist or force_debase is requested, create it directly from spCf16
                if !isfile(debase_target) || force_debase
                    if isdefined(@__MODULE__, :Data) && isdefined(Data, :make_fullrange_debase)
                        println("[$i/$(length(names))] Creating full-range pulsar.debase.txt from $(basename(spcf16_file)) for $name ...")
                        try
                            result = Data.make_fullrange_debase(outdir; spCf16_file=basename(spcf16_file), force=force_debase)
                            if !isnothing(result) && isfile(result)
                                debase_file = result
                            end
                        catch debase_err
                            @warn "[$i/$(length(names))] Failed to create debase file from spCf16 for $name: $debase_err"
                        end
                    else
                        @warn "[$i/$(length(names))] Cannot create debase file for $name: Data.make_fullrange_debase not available"
                    end
                else
                    debase_file = debase_target
                end
            elseif isfile(debase_target)
                debase_file = debase_target
            end

            if !isfile(params_file) || isempty(debase_file)
                missing_arr = []
                if !isfile(params_file); push!(missing_arr, "params.json"); end
                if isempty(debase_file); push!(missing_arr, "full-range pulsar.debase.txt (no spCf16 found to generate from; high/low files ignored)"); end
                
                missing_str = join(missing_arr, ", ")
                @warn "[$i/$(length(names))] Skipping $name — required full-range files missing in $outdir: $missing_str"
                continue
            end
            
            # Load params to get bin_st and bin_end
            params = JSON.parsefile(params_file)
            bin_st = get(params, "bin_st", 1)
            bin_end = get(params, "bin_end", nothing)
            
            # Fast streaming load of data matrix from PSRCHIVE ASCII dump
            X = _read_debase_ascii(debase_file)
            if size(X, 1) == 0 || size(X, 2) == 0
                @warn "[$i/$(length(names))] Skipping $name — failed to parse debase file $debase_file."
                continue
            end
            
            if bin_end === nothing
                bin_end = size(X, 2)
            end
            
            println("[$i/$(length(names))] Analyzing $name (P3 source: $(isnothing(known_p3) ? "auto-detect" : "known=$(known_p3)")) ...")
            res = lrfs_phase_track(X, bin_st, bin_end; known_p3=known_p3)

            push!(results, (name,
                            isnothing(known_p3) ? NaN : known_p3,  # P3 from list
                            res.p3_pulses,                          # P3 actually used (quantised)
                            res.phase_slope,
                            string(res.p3_source),
                            string(res.classification)))

        catch e
            @warn "[$i/$(length(names))] Failed to analyze $name: $e — skipping."
            continue
        end
    end

    # Save to CSV
    open(out_csv, "w") do io
        write(io, "Name,P3_Known,P3_Used,Phase_Slope,P3_Source,Classification\n")
        for r in results
            write(io, "$(r[1]),$(r[2]),$(r[3]),$(r[4]),$(r[5]),$(r[6])\n")
        end
    end
    println("\nSaved results to $out_csv")
    return out_csv
end

"""
    plot_lrfs_results(csv_file::String, out_plot::String="lrfs_chart.png")

Reads the CSV generated by `batch_analyze_lrfs` and creates a scatter plot 
to visually separate amplitude-modulated pulsars from drifting pulsars 
based solely on the LRFS phase tracking.
"""
function plot_lrfs_results(csv_file::String, out_plot::String="lrfs_chart.png")
    if !isfile(csv_file)
        @warn "CSV file not found: $csv_file"
        return
    end
    data, header = readdlm(csv_file, ',', header=true)
    if isempty(data)
        @warn "No successful results in $csv_file to plot."
        return
    end
    
    names   = data[:, 1]
    p3_known = [tryparse(Float64, string(v)) for v in data[:, 2]]
    p3_used  = convert(Vector{Float64}, data[:, 3])
    slopes   = convert(Vector{Float64}, data[:, 4])
    # col 5 = P3_Source, col 6 = Classification
    classes  = data[:, 6]

    # Use known P3 for x-axis where available, fall back to LRFS-detected P3
    p3_plot = [(!isnothing(p3_known[i]) && !isnan(p3_known[i])) ? p3_known[i] : p3_used[i]
               for i in 1:length(names)]
    
    # Setup plot
    PyPlot.figure(figsize=(10, 8))
    
    # Separate by class for coloring
    idx_am = classes .== "amplitude_modulation"
    idx_drift = classes .== "pure_drift"
    idx_undef = classes .== "undetermined"
    
    # Scatter: X = P3 (Pulses), Y = Phase Slope (rad/bin)
    PyPlot.scatter(p3_plot[idx_am], abs.(slopes[idx_am]), color="red", label="Amplitude Modulation (Flat phase)", alpha=0.7, s=60)
    PyPlot.scatter(p3_plot[idx_drift], abs.(slopes[idx_drift]), color="blue", label="Subpulse Drift (Slanted phase)", alpha=0.7, s=60)
    
    if any(idx_undef)
        PyPlot.scatter(p3_plot[idx_undef], abs.(slopes[idx_undef]), color="gray", label="Undetermined (Noise)", alpha=0.7, s=40, marker="x")
    end
    
    # Add labels for strong drifters (slope > 0.15) to identify them easily
    for i in 1:length(names)
        if classes[i] == "pure_drift" && abs(slopes[i]) > 0.15
            PyPlot.annotate(names[i], (p3_plot[i], abs(slopes[i])), fontsize=8, alpha=0.6)
        end
    end
    
    PyPlot.xlabel("Dominant P3 Period (pulses)")
    PyPlot.ylabel("|LRFS Phase Slope| (radians/bin)")
    PyPlot.title("Pulsar Classification via LRFS Phase Tracking")
    
    # Threshold line
    PyPlot.axhline(0.05, color="gray", linestyle="--", alpha=0.5, label="Classification Threshold (0.05)")
    
    # Log scale on X axis often looks better for periods
    # PyPlot.xscale("log") 
    
    PyPlot.legend()
    PyPlot.grid(true, alpha=0.3)
    
    PyPlot.savefig(out_plot, dpi=300, bbox_inches="tight")
    PyPlot.close()
    
    println("Saved classification chart to $out_plot")
end

end # module
