module Data
    using Glob
    using FITSIO
    using ProgressMeter
    using Dates

    include("functions.jl")

    """
    Zero selected pulses
    """
    function zap!(data; ranges=nothing)
        pulses, bins = size(data)
        zaps = []
        for rn in ranges
            if length(rn) == 1
                push!(zaps, rn)
            else
                zaps = vcat(zaps, collect(rn[1]:rn[2]))
            end
            for z in zaps
                for i in 1:bins
                    data[z,i] = 0
                end
            end
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
        run(pipeline(`pdv -t -F -p $infile`, stdout="$outfile", stderr="errs.txt"))
        # change -t to -A to get frequancy information
        #@showprogress 1 for i in 1:pn  # psrchive indexing
        #end
    end

    """
    Process data with PSRCHIVE and PSRSALSA
    """
    function process_psrdata(indir, outdir; outfile="pulsar.spCF", files=nothing)
        # Tworzenie nowego folderu w outdir z datą i godziną
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
        new_outdir = joinpath(outdir, "output_$timestamp")
        mkpath(new_outdir)  # Tworzenie folderu

        if files === nothing
            # Znajdź wszystkie pliki .spCF w katalogu wejściowym
            files = filter(f -> endswith(f, ".spCF"), readdir(indir))
        end
        
        # Sortowanie plików na podstawie numerów impulsów
        sort!(files, by = f -> begin
            m = match(r"_(\d+)-(\d+)\.spCF$", f)
            if isnothing(m)
                return typemax(Int)  # Pliki, które nie pasują do formatu, idą na koniec
            else
                return parse(Int, m.captures[1])  # Sortowanie po numerze początkowego impulsu
            end
        end)
    
        if isempty(files)
            error("No .spCF files found in directory: $indir")
        end
    
        println("Processing files in order:")
        for (i, f) in enumerate(files)
            println("$i. $f")
        end
    
        file_names = [joinpath(indir, file) for file in files]

        outfile = joinpath(new_outdir, outfile)  # Wskazujemy nowy folder

        # Łączenie wszystkich plików
        run(pipeline(`psradd $file_names -o $outfile`, stderr="errs.txt"))

        # Debasing danych
        run(pipeline(`pmod -device "/xw" -debase $outfile`, `tee pmod_output.txt`))
        output = read("pmod_output.txt", String)
        rm("pmod_output.txt")  # Cleanup

        # Extract onpulse values
        m = match(r"-onpulse '(\d+) (\d+)'", output)
        if !isnothing(m)
            bin_st, bin_end = parse.(Int, m.captures)
            region_length = bin_end - bin_st + 1
            if region_length % 2 != 0
                println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                bin_end -= 1
                println("Adjusted onpulse range: $bin_st to $bin_end")
            end
            println("Found onpulse range: $bin_st to $bin_end")
        end

        debased_file = replace(outfile, ".spCF" => ".debase.gg")

        # Calculate 2dfs and lrfs
        run(pipeline(`pspec -w -2dfs -lrfs -profd "/NULL" -onpulsed "/NULL" -2dfsd "/NULL" -lrfsd "/NULL" -nfft 256 -onpulse "$(bin_st) $(bin_end)" $debased_file`, stderr="errs.txt"))

        # Find P3
        run(pipeline(`pspecDetect -v -device "/xw" $debased_file`, `tee pspecDetect_output.txt`))
        output = read("pspecDetect_output.txt", String)
        rm("pspecDetect_output.txt")  # Cleanup

        p3_matches = collect(eachmatch(r"P3\[P0\]\s*=\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)", output))
        if !isempty(p3_matches)
            last_match = p3_matches[end]
            p3_value = parse(Float64, last_match.captures[1])
            p3_error = parse(Float64, last_match.captures[2])
            println("Found P3 = $p3_value ± $p3_error P0")
        end
        ybins = Functions.find_ybins(p3_value)
        println("Number of ybins: $ybins")

        run(pipeline(`pfold  -p3fold "$p3_value $ybins" -onpulse "$bin_st $bin_end" -onpulsed "/NULL" -p3foldd "/NULL" -w -oformat ascii $debased_file`, stderr="errs.txt"))

        # Przenoszenie plików do nowo utworzonego folderu
        for file in readdir(outdir)
            if !isdir(joinpath(outdir, file))  # Upewnij się, że to plik, a nie folder
                mv(joinpath(outdir, file), new_outdir)  # Przenieś pliki do nowego folderu
            end
        end

        return bin_st - 20, bin_end + 20
    end

end  # module

