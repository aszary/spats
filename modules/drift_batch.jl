module DriftBatch

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

# Include DriftDiagnostics
include("DriftDiagnostics.jl")
using .DriftDiagnostics

export batch_analyze_drift, plot_drift_results

"""
    batch_analyze_drift(vpmout::String, list_file::String, out_csv::String)

Reads a list of pulsars from `list_file`, finds their output directories in `vpmout`,
loads the single-pulse data and parameters, and runs `analyze_drift`.
Saves the results to `out_csv`.
"""
function batch_analyze_drift(vpmout::String, list_file::String, out_csv::String="drift_classifications.csv")
    if !isfile(list_file)
        error("Pulsar list file not found: $list_file")
    end

    names = String[]
    for line in eachline(list_file)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        push!(names, String(first(split(s))))
    end
    
    println("Found $(length(names)) pulsars in list.")
    
    results = []
    
    for (i, name) in enumerate(names)
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
                @warn "[$i/$(length(names))] Skipping $name — Output dir not found. Checked: $(candidates[1]) etc."
                continue
            end
            
            params_file = joinpath(outdir, "params.json")
            
            # Check for different debase file naming conventions
            debase_file = ""
            for name_variant in ["pulsar.debase.txt", "pulsar_high_debase.txt", "pulsar_low_debase.txt"]
                cand_debase = joinpath(outdir, name_variant)
                if isfile(cand_debase)
                    debase_file = cand_debase
                    break
                end
            end
            
            # If no debase file found but spCf16 exists, create pulsar.debase.txt on the fly
            if isempty(debase_file) && isfile(joinpath(outdir, "pulsar.spCf16")) && isfile(params_file)
                if isdefined(@__MODULE__, :Data) && isdefined(Data, :make_fullrange_debase)
                    println("[$i/$(length(names))] Creating pulsar.debase.txt for $name ...")
                    try
                        result = Data.make_fullrange_debase(outdir)
                        if !isnothing(result) && isfile(result)
                            debase_file = result
                        end
                    catch debase_err
                        @warn "[$i/$(length(names))] Failed to create debase file for $name: $debase_err"
                    end
                else
                    @warn "[$i/$(length(names))] Cannot create debase file for $name: Data.make_fullrange_debase not available"
                end
            end
            
            if !isfile(params_file) || isempty(debase_file)
                missing_arr = []
                if !isfile(params_file); push!(missing_arr, "params.json"); end
                if isempty(debase_file); push!(missing_arr, "any *debase.txt (and no spCf16 to generate from)"); end
                
                missing_str = join(missing_arr, ", ")
                @warn "[$i/$(length(names))] Skipping $name — required files missing in $outdir: $missing_str"
                continue
            end
            
            # Load params to get bin_st and bin_end
            params = JSON.parsefile(params_file)
            bin_st = get(params, "bin_st", 1)
            bin_end = get(params, "bin_end", nothing)
            
            # Load data matrix from PSRCHIVE ASCII dump
            lines = readlines(debase_file)
            if isempty(lines)
                @warn "[$i/$(length(names))] Skipping $name — debase file $debase_file is empty."
                continue
            end
            header = split(lines[1])
            n_pulses = parse(Int, header[6])
            n_bins = parse(Int, header[12])
            
            X = zeros(Float64, n_pulses, n_bins)
            for j in 2:length(lines)
                res = split(lines[j])
                pulse = parse(Int, res[1]) + 1
                bin = parse(Int, res[3]) + 1
                X[pulse, bin] = parse(Float64, res[4])
            end
            
            if bin_end === nothing
                bin_end = n_bins
            end
            
            println("[$i/$(length(names))] Analyzing $name ...")
            res = analyze_drift(X, bin_st, bin_end)
            
            push!(results, (name, res.travel_sig, res.asymmetry_ratio, res.phase_gradient, res.mean_peak_lag, res.score, res.classification))
            
        catch e
            @warn "[$i/$(length(names))] Failed to analyze $name: $e — skipping."
            continue
        end
    end
    
    # Save to CSV
    open(out_csv, "w") do io
        write(io, "Name,Travel_Sig,Asymmetry_Ratio,Phase_Gradient,Mean_Peak_Lag,Score,Classification\n")
        for r in results
            write(io, "$(r[1]),$(r[2]),$(r[3]),$(r[4]),$(r[5]),$(r[6]),$(r[7])\n")
        end
    end
    println("\nSaved results to $out_csv")
    return out_csv
end

"""
    plot_drift_results(csv_file::String, out_plot::String="drift_plot.png")

Reads the CSV generated by `batch_analyze_drift` and creates a scatter plot 
to visually separate amplitude-modulated pulsars from drifting pulsars.
"""
function plot_drift_results(csv_file::String, out_plot::String="drift_plot.png")
    if !isfile(csv_file)
        @warn "CSV file not found: $csv_file"
        return
    end
    data, header = readdlm(csv_file, ',', header=true)
    if isempty(data)
        @warn "No successful results in $csv_file to plot."
        return
    end
    
    names = data[:, 1]
    asym = convert(Vector{Float64}, data[:, 3])
    grad = convert(Vector{Float64}, data[:, 4])
    lags = convert(Vector{Float64}, data[:, 5])
    classes = data[:, 7]
    
    # Setup plot
    PyPlot.figure(figsize=(10, 8))
    
    # Separate by class for coloring
    idx_am = classes .== "amplitude_modulation"
    idx_drift = classes .== "pure_drift"
    idx_bi = classes .== "bi_drift"
    idx_undef = classes .== "undetermined"
    
    # Scatter: X = Asymmetry (2DFS), Y = Phase Gradient (Cross-Spec)
    PyPlot.scatter(abs.(asym[idx_am]), abs.(grad[idx_am]), color="red", label="Amplitude Modulation (P3-only)", alpha=0.7, s=60)
    PyPlot.scatter(abs.(asym[idx_drift]), abs.(grad[idx_drift]), color="blue", label="Pure Drift", alpha=0.7, s=60)
    PyPlot.scatter(abs.(asym[idx_bi]), abs.(grad[idx_bi]), color="green", label="Bi-Drift", alpha=0.7, s=60)
    PyPlot.scatter(abs.(asym[idx_undef]), abs.(grad[idx_undef]), color="gray", label="Undetermined", alpha=0.7, s=40, marker="x")
    
    # Add labels for extreme cases (optional, could get crowded)
    for i in 1:length(names)
        if classes[i] == "pure_drift" && (abs(asym[i]) > 0.3 || abs(grad[i]) > 0.3)
            PyPlot.annotate(names[i], (abs(asym[i]), abs(grad[i])), fontsize=8, alpha=0.6)
        end
    end
    
    PyPlot.xlabel("|2DFS Asymmetry Ratio|")
    PyPlot.ylabel("|Cross-Spectrum Phase Gradient|")
    PyPlot.title("Subpulse Modulation Diagnostics")
    
    PyPlot.axhline(0.1, color="gray", linestyle="--", alpha=0.5)
    PyPlot.axvline(0.1, color="gray", linestyle="--", alpha=0.5)
    
    PyPlot.legend()
    PyPlot.grid(true, alpha=0.3)
    
    PyPlot.savefig(out_plot, dpi=300, bbox_inches="tight")
    PyPlot.close()
    
    println("Saved classification chart to $out_plot")
end

end # module
