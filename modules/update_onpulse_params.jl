"""
update_onpulse_params.jl

Modifies on-pulse regions in `params.json` for all pulsars based on `input/onpulsefixing.csv`:
- Parses MP components (C1-C4):
    - Sets `bin_st` to the minimum bin across all MP components (rounded to integer).
    - Sets `bin_end` to the maximum bin across all MP components (rounded to integer).
- If IP (Interpulse) components exist (C1-C4):
    - Sets `IP_start` to the minimum bin across all IP components (rounded to integer).
    - Sets `IP_end` to the maximum bin across all IP components (rounded to integer).
- Creates a backup `params_pre_reonpuls.json` before overwriting `params.json`.
- Can search across multiple base directories (e.g. `/home/psr/output`, `./output`, etc.).
"""

module OnpulseFixer

using JSON

export parse_onpulse_csv, update_params_files, update_single_params!

function parse_float_safe(val)
    if ismissing(val) || isnothing(val)
        return nothing
    end
    s = strip(string(val))
    if isempty(s)
        return nothing
    end
    s_clean = replace(s, "," => ".")
    try
        v = parse(Float64, s_clean)
        return v != 0.0 ? v : nothing
    catch
        return nothing
    end
end

"""
    _parse_csv_line(line::AbstractString) -> Vector{String}

RFC4180-compliant CSV line parser handling quoted fields containing commas.
Uses pure Julia standard string operations without external package dependencies.
"""
function _parse_csv_line(line::AbstractString)
    fields = String[]
    current = IOBuffer()
    in_quotes = false
    i = 1
    chars = collect(line)
    n = length(chars)

    while i <= n
        c = chars[i]
        if c == '"'
            if in_quotes && i < n && chars[i+1] == '"'
                # Escaped quote: ""
                write(current, '"')
                i += 1
            else
                in_quotes = !in_quotes
            end
        elseif c == ',' && !in_quotes
            push!(fields, String(take!(current)))
        else
            write(current, c)
        end
        i += 1
    end
    push!(fields, String(take!(current)))
    return fields
end

"""
    parse_onpulse_csv(csv_path::String) -> Dict{String, Dict{String, Any}}

Reads onpulsefixing.csv and aggregates min/max bins for MP and IP using standard library only.
"""
function parse_onpulse_csv(csv_path::String="input/onpulsefixing.csv")
    if !isfile(csv_path)
        error("CSV file not found at: $csv_path")
    end

    results = Dict{String, Dict{String, Any}}()

    # Define component column index triplets for bins (1-based index in the CSV):
    # C1: col 16 (region), col 17 (bin_start), col 18 (bin_end)
    # C2: col 19 (region), col 20 (bin_start), col 21 (bin_end)
    # C3: col 22 (region), col 23 (bin_start), col 24 (bin_end)
    # C4: col 25 (region), col 26 (bin_start), col 27 (bin_end)
    comp_cols = [
        (16, 17, 18),
        (19, 20, 21),
        (22, 23, 24),
        (25, 26, 27)
    ]

    open(csv_path, "r") do io
        # Read and discard header
        header_line = readline(io)
        for line in eachline(io)
            stripped_line = strip(line)
            isempty(stripped_line) && continue
            row = _parse_csv_line(stripped_line)
            length(row) < 1 && continue

            jname = strip(row[1])
            isempty(jname) && continue

            mp_starts = Float64[]
            mp_ends   = Float64[]
            ip_starts = Float64[]
            ip_ends   = Float64[]

            for (reg_idx, st_idx, en_idx) in comp_cols
                if reg_idx > length(row)
                    continue
                end
                reg = uppercase(strip(row[reg_idx]))
                st = st_idx <= length(row) ? parse_float_safe(row[st_idx]) : nothing
                en = en_idx <= length(row) ? parse_float_safe(row[en_idx]) : nothing

                if reg == "MP"
                    !isnothing(st) && push!(mp_starts, st)
                    !isnothing(en) && push!(mp_ends, en)
                elseif reg == "IP"
                    !isnothing(st) && push!(ip_starts, st)
                    !isnothing(en) && push!(ip_ends, en)
                end
            end

            entry = Dict{String, Any}()
            if !isempty(mp_starts) && !isempty(mp_ends)
                entry["bin_st"]  = round(Int, minimum(mp_starts))
                entry["bin_end"] = round(Int, maximum(mp_ends))
            end

            if !isempty(ip_starts) && !isempty(ip_ends)
                entry["IP_start"] = round(Int, minimum(ip_starts))
                entry["IP_end"]   = round(Int, maximum(ip_ends))
            end

            results[jname] = entry
        end
    end

    return results
end

"""
    update_single_params!(params_path::String, onpulse_info::Dict{String, Any}; dry_run::Bool=false)

Updates a single params.json file, saving a backup as params_pre_reonpuls.json.
"""
function update_single_params!(params_path::String, onpulse_info::Dict{String, Any}; dry_run::Bool=false)
    if !isfile(params_path)
        return false
    end

    # Read existing params
    params = open(params_path, "r") do f
        JSON.parse(f)
    end

    backup_path = joinpath(dirname(params_path), "params_pre_reonpuls.json")

    # Update MP values
    if haskey(onpulse_info, "bin_st") && haskey(onpulse_info, "bin_end")
        params["bin_st"]  = onpulse_info["bin_st"]
        params["bin_end"] = onpulse_info["bin_end"]
    end

    # Update or add IP values
    if haskey(onpulse_info, "IP_start") && haskey(onpulse_info, "IP_end")
        params["IP_start"] = onpulse_info["IP_start"]
        params["IP_end"]   = onpulse_info["IP_end"]
    end

    if !dry_run
        # Create backup if not already present or refresh it
        cp(params_path, backup_path; force=true)

        # Write updated params
        open(params_path, "w") do f
            JSON.print(f, params, 4)
        end
    end

    return true
end

"""
    update_params_files(search_dirs::Vector{String}=["/home/psr/output", "output"];
                        csv_path::String="input/onpulsefixing.csv",
                        dry_run::Bool=false)

Finds and updates params.json files in all target directories for each pulsar in the CSV.
"""
function update_params_files(search_dirs::Vector{String}=["/home/psr/output", "output"];
                            csv_path::String="input/onpulsefixing.csv",
                            dry_run::Bool=false)
    println("Loading on-pulse regions from $csv_path ...")
    onpulse_data = parse_onpulse_csv(csv_path)
    println("Loaded $(length(onpulse_data)) pulsar definitions.")

    updated_count = 0
    not_found = String[]

    for (jname, info) in onpulse_data
        found = false
        # Possible folder name variants: Jname, Jname_16, Jname_32
        folder_variants = [jname, "$(jname)_16", "$(jname)_32"]

        for sdir in search_dirs
            if !isdir(sdir)
                continue
            end
            for fvar in folder_variants
                candidate_dir = joinpath(sdir, fvar)
                pfile = joinpath(candidate_dir, "params.json")
                if isfile(pfile)
                    update_single_params!(pfile, info; dry_run=dry_run)
                    updated_count += 1
                    found = true
                    println("  [OK] Updated $pfile -> bin_st=$(get(info, "bin_st", "N/A")), bin_end=$(get(info, "bin_end", "N/A"))" *
                            (haskey(info, "IP_start") ? ", IP_start=$(info["IP_start"]), IP_end=$(info["IP_end"])" : ""))
                end
            end
        end

        if !found
            push!(not_found, jname)
        end
    end

    println("\nSummary:")
    println("  Updated params.json files: $updated_count")
    println("  Pulsars not located in search dirs: $(length(not_found))")
    return (updated=updated_count, not_found=not_found)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    search_paths = length(ARGS) > 0 ? ARGS : ["/home/psr/output", "output"]
    OnpulseFixer.update_params_files(search_paths)
end
