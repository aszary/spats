module SpaTs
    using ArgParse
    using Glob
    using JSON
    using FITSIO
    using PyPlot
    using PyCall
    const mcolors = pyimport("matplotlib.colors")
    using FFTW
    using Statistics


    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


            
 """
Renders and saves a 2DFS plot from the file pulsar.debase.1.2dfs using PyPlot.


Arguments:
- outdir: The output directory where the pulsar data is stored.
- pulsar_name: The name of the pulsar to create the plot for.
- show_plot: A boolean flag to decide whether to display the plot (default: true).
"""
function plot_2dfs(outdir::String, pulsar_name::String; show_plot::Bool=true)


    filepath = joinpath(outdir, "pulsar.debase.1.2dfs")


    if !isfile(filepath)
        println("File does not exist: $filepath")
        return
    end

    println("Inspecting FITS file: $filepath")

    f = nothing
    data = nothing
    try
        f = FITS(filepath)


        for (i, hdu) in enumerate(f)
            println("HDU $i:")
            println("  Type: ", typeof(hdu))

            try
                if hdu isa FITSIO.ImageHDU
                    img = read(hdu)
                    if ndims(img) == 2
                        println("  -> Found 2D image of size ", size(img))
                        data = img
                        break
                    else
                        println("  -> Not a 2D image (ndims=$(ndims(img)))")
                    end


                elseif hdu isa FITSIO.TableHDU
                    names = FITSIO.colnames(hdu)
                    println("  Columns: ", names)


                    for name in names
                        col_data = read(hdu, name)
                        println("    Column '$name' -> type: ", typeof(col_data), ", size: ", size(col_data))


                        if isa(col_data, AbstractArray) && ndims(col_data) == 2 && eltype(col_data) <: Number
                            println("  ✅ Found 2D numeric column '$name' in HDU $i")
                            data = col_data
                            break
                        end
                    end
                end
            catch e
                println("  -> Failed to read HDU $i: $e")
            end


            if data !== nothing
                break
            end
        end


        if data === nothing
            println("❌ No suitable 2D data found in FITS file.")
            close(f)
            return
        end


        close(f)


        n_p3, n_p2 = size(data)
        p3_range = range(0, stop=0.5, length=n_p3)
        p2_range = range(160, stop=200, length=n_p2)  # Pulse longitude in deg, ograniczony do 160–200


        fig, ax = subplots()
        im = ax.imshow(data';
            extent=[160, 200, 0, 0.5],  # Zamieniono osie x i y
            origin="lower",
            aspect="auto",
            cmap="gray_r",  # Dodano przecinek
            vmin=0,  # Minimalna wartość skali kolorów
            vmax=50000  # Maksymalna wartość skali kolorów
        )


        ax.set_xlabel("Pulse longitude (deg)")  # Etykieta osi x
        ax.set_ylabel("Fluctuation frequency (P/P3)")  # Etykieta osi y
        ax.set_title("2DFS – $pulsar_name")
        colorbar(im, ax=ax, label="Power")


        savepath = joinpath(outdir, "2dfs_" * pulsar_name * ".png")
        savefig(savepath)
        println("✅ 2DFS plot saved to: $savepath")


        if show_plot
            show()
        else
            close(fig)
        end


    catch e
        println("❌ Error handling FITS file: $e")
        if f !== nothing
            close(f)
        end
    end
end



function process_psrdata(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="test", show_=true, repeat_num=4)
        
    end

function Plot_2dfs_zmiany(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    f = FITS(filepath)
    data = read(f[4], "DATA")
    close(f)

    if data === nothing
        println("❌ No suitable 2D data found in 2DFS.")
        return
    end

    n_bins, n_subints = size(data)
    data = data'  # Transponujemy: (subints, bins)

    # === Zakresy osi ===
    pulse_longitude = LinRange(160.0, 200.0, n_bins)    # Oś X: długość impulsu
    fluct_freq = LinRange(0.0, 0.5, n_subints)          # Oś Y: częstotliwość fluktuacji (P/P3)

    # === Parametry mapy kolorów ===
    vmin = 0.0
    vmax = 0.07

    fig, ax = subplots(figsize=(7, 6))

    im = ax.imshow(data;
        extent=[160, 200, 0.0, 0.5],
        origin="lower",
        aspect="auto",
        cmap="gray",
        vmin=vmin,
        vmax=vmax
    )

    # === Opisy osi i tytuł ===
    ax.set_xlabel("Pulse longitude (deg)")
    ax.set_ylabel("Fluctuation frequency (P/P₃)")
    ax.set_title("2DFS – $pulsar_name")

    # === Pasek kolorów ===
    cbar = colorbar(im, ax=ax)
    cbar.set_label("Power")
    cbar.set_ticks([0.0, 0.02, 0.04, 0.06])

    # === Zapis i/lub wyświetlenie ===
    savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
    savefig(savepath, dpi=150)
    println("✅ 2DFS plot saved to: $savepath")

    if show_plot
        show()
    else
        close(fig)
    end
end



function Plot_2dfs_zmiany_pplot(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("🔍 Reading 2DFS from: $filepath")
    f = FITS(filepath)
    data = read(f[4], "DATA")
    close(f)

    data = data'  # transpozycja

    println("Dane min: ", minimum(data))
    println("Dane max: ", maximum(data))
    println("Czy są NaN: ", any(isnan, data))
    println("Czy są Inf: ", any(isinf, data))

    # Konwersja do Float64, by można było przypisać 1e-10
    data_clean = Float64.(data)

    data_clean[isnan.(data_clean)] .= 1e-10
    data_clean[data_clean .<= 0] .= 1e-10

    NBIN = 512
    left_bin = 160
    right_bin = 200
    n_region = right_bin - left_bin + 1

    p2min = -NBIN / 2
    p2max = p2min + NBIN * n_region / NBIN

    p3min = 0.0
    p3max = 0.5

    fig, ax = subplots(figsize=(7,6))

    im = ax.imshow(data_clean;
        extent=[p2min, p2max, p3min, p3max],
        origin="lower",
        aspect="auto",
        cmap="gray",
        norm=matplotlib[:colors][:LogNorm](vmin=1e-4, vmax=0.07)
    )

    ax.set_xlabel("Pulse longitude (deg)")
    ax.set_ylabel("Fluctuation frequency (P/P₃)")
    ax.set_title("2DFS – $pulsar_name")

    cbar = colorbar(im, ax=ax)
    cbar.set_label("Power")
    cbar.set_ticks([0.0, 0.02, 0.04, 0.06])

    savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
    savefig(savepath, dpi=150)
    println("✅ 2DFS plot saved to: $savepath")

    if show_plot
        show()
    else
        close(fig)
    end
end







function detailed_check_fits(outdir::String, pulsar_name::String)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end
    
    println("🔍 Opening FITS file: $filepath")
    f = FITS(filepath)

    println("\n=== Nagłówki HDU ===")
    for (i, hdu) in enumerate(f)
        println("\n--- HDU #$i ---")
        hdr = read_header(hdu)
        println(hdr)  # wypisujemy nagłówek jako tekst
    end

    println("\n🔍 Reading DATA from HDU #4")
    data = read(f[4], "DATA")
    close(f)

    println("\n=== Podstawowe informacje o danych ===")
    println("Rozmiar oryginalnych danych: ", size(data))
    println("Typ danych: ", eltype(data))
    println("Minimum: ", minimum(data))
    println("Maksimum: ", maximum(data))
    println("Średnia: ", mean(Float64.(data)))
    println("Mediana: ", median(Float64.(data)))
    println("Liczba NaN: ", count(isnan, data))
    println("Liczba Inf: ", count(isinf, data))
    println("Liczba zer: ", count(x -> x == 0, data))
    println("Liczba wartości <= 0: ", count(x -> x <= 0, data))

    nrows = min(5, size(data,1))
    ncols = min(10, size(data,2))
    println("\nPróbka danych (pierwsze $nrows wierszy i $ncols kolumn):")
    for r in 1:nrows
        println(data[r, 1:ncols])
    end

    data_t = data'
    println("\nPo transpozycji rozmiar: ", size(data_t))
    println("\nPróbka danych po transpozycji (pierwsze $nrows wierszy i $ncols kolumn):")
    for r in 1:nrows
        println(data_t[r, 1:ncols])
    end

    println("\n✅ Sprawdzenie zakończone.")
end






function Plot_2dfs_simple(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("🔍 Reading 2DFS from: $filepath")
    f = FITS(filepath)
    data = read(f[4], "DATA")
    close(f)

    data = data'  # transpozycja
    dataf = Float64.(data)

    NBIN = 512

    # Dopasuj indeksy do rozmiaru danych
    left_bin = 1
    right_bin = size(dataf, 2)

    n_region = right_bin - left_bin + 1

    dataf_region = dataf[:, left_bin:right_bin]

    # Zamiana zer i ujemnych na minimalną dodatnią
    dataf_region[dataf_region .<= 0] .= 1e-10

    vmin = minimum(dataf_region)
    vmax = maximum(dataf_region)

    p2min = -NBIN / 2
    p2max = p2min + n_region  # dokładnie n_region na osi x

    p3min = 0.0
    p3max = 0.5

    fig, ax = subplots(figsize=(7,6))

    im = ax.imshow(dataf_region;
        extent=[p2min, p2max, p3min, p3max],
        origin="lower",
        aspect="auto",
        cmap="gray_r",
        vmin=vmin,
        vmax=vmax
    )

    ax.set_xlabel("Pulse longitude (deg)")
    ax.set_ylabel("Fluctuation frequency (P/P₃)")
    ax.set_title("2DFS – $pulsar_name")

    cbar = colorbar(im, ax=ax)
    cbar.set_label("Power")
    cbar.set_ticks([vmin, (vmin+vmax)/2, vmax])

    savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
    savefig(savepath, dpi=150)
    println("✅ 2DFS plot saved to: $savepath")

    if show_plot
        show()
    else
        close(fig)
    end
end



function Plot_2dfs2(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    lines = readlines(filepath)
    data = Float64[]
    for line in lines
        # Pomijamy puste linie
        if isempty(strip(line))
            continue
        end
        vals = split(strip(line))
        # Sprawdzamy, czy linia zawiera same liczby
        can_parse_all = all(v -> try parse(Float64, v); true catch; false end, vals)
        if can_parse_all
            for v in vals
                push!(data, parse(Float64, v))
            end
        else
            # Pomijamy linie, które nie da się sparsować do liczb
            continue
        end
    end

    # Podaj wymiary zgodnie z tym, czego się spodziewasz w pliku
    nrows = 129
    ncols = 92

    if length(data) != nrows*ncols
        println("⚠️ Warning: Unexpected data length: $(length(data)), expected $(nrows*ncols)")
        return
    end

    arr = reshape(data, ncols, nrows)'
    x = range(160, stop=200, length=ncols)
    y = range(0, stop=0.5, length=nrows)

    figure()
    imshow(arr, extent=(minimum(x), maximum(x), minimum(y), maximum(y)),
           aspect="auto", origin="lower", cmap="viridis")
    colorbar(label="Intensity")
    xlabel("Pulse longitude (deg)")
    ylabel("Fluctuation frequency (P/P3)")
    title("2DFS plot for $pulsar_name")
    if show_plot
        show()
    end
end



function Plot_ostateczny(outdir::String, pulsar_name::String, period::Float64; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("🔍 Reading 2DFS from: $filepath")
    f = FITS(filepath)

    # Nie pobieramy okresu z pliku tylko bierzemy argument period
    println("✅ Używany okres (period) z argumentu: ", period)

    # Bezpiecznie czytamy dane
    data_col = nothing
    try
        data_col = read(f[4], "DATA")
    catch e
        println("❌ Błąd przy odczycie kolumny DATA: ", e)
        close(f)
        return
    end

    close(f)

    # Sprawdzenie, czy data_col jest OK
    if data_col === nothing
        println("❌ data_col nie został poprawnie wczytany.")
        return
    end

    data = data_col'  # transpozycja
    dataf = Float64.(data)

    n_y, n_x = size(dataf)
    pulse_longitudes = range(0, stop=360, length=n_x)
    fluct_freqs = range(0, stop=0.5, length=n_y)
    P_over_P3 = period ./ fluct_freqs
    P_over_P3[1] = NaN

    x_indices = findall(x -> 160 <= x <= 200, pulse_longitudes)
    if isempty(x_indices)
        println("❌ No pulse longitude values found in the range 160–200°.")
        return
    end

    cropped_data = dataf[:, x_indices]
    cropped_longitudes = pulse_longitudes[x_indices]

    fig, ax = subplots()
    im = ax.imshow(cropped_data;
        origin="lower",
        aspect="auto",
        extent=[minimum(cropped_longitudes), maximum(cropped_longitudes),
                minimum(P_over_P3[2:end]), maximum(P_over_P3[2:end])],
        cmap="gray",
        vmin=1.0,
        vmax=200
    )
    im.set_clim(0.0, 0.07)

    ax.set_xlabel("Pulse longitude (degrees)")
    ax.set_ylabel("Fluctuation frequency (P/P3)")
    ax.set_title("2DFS: $pulsar_name")

    colorbar(im, ax=ax, label="Intensity")

    if show_plot
        show()
    else
        savefig(joinpath(outdir, pulsar_name, "$pulsar_name-2dfs.png"))
    end

    println("✅ Wygenerowano wykres dla $pulsar_name.")
end



function Plot_ostateczny_simple(outdir::String, pulsar_name::String, period::Float64; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("🔍 Reading 2DFS from: $filepath")
    f = FITS(filepath)
    
    data_col = nothing
    try
        data_col = read(f[4], "DATA")
    catch e
        println("❌ Error reading DATA column: ", e)
        close(f)
        return
    end

    close(f)

    if data_col === nothing
        println("❌ DATA column not loaded correctly")
        return
    end

    data = data_col'
    dataf = Float64.(data)

    n_y, n_x = size(dataf)
    pulse_longitudes = range(0, stop=360, length=n_x)
    fluct_freqs = range(0, stop=0.5, length=n_y)
    P_over_P3 = period ./ fluct_freqs
    P_over_P3[1] = NaN

    fig, ax = subplots()
    im = ax.imshow(dataf;
        origin="lower",
        aspect="auto",
        extent=[minimum(pulse_longitudes), maximum(pulse_longitudes),
                minimum(P_over_P3[2:end]), maximum(P_over_P3[2:end])],
        cmap="gray_r",
        norm=mcolors.LogNorm(vmin=1e-4, vmax=1.0)
    )

    #im.set_clim(0.0, 0.07)


    ax.set_xlabel("Pulse longitude (degrees)")
    ax.set_ylabel("Fluctuation frequency (P/P3)")
    ax.set_title("2DFS: $pulsar_name")

    colorbar(im, ax=ax, label="Intensity")

    if show_plot
        show()
    else
        savefig(joinpath(outdir, pulsar_name, "$pulsar_name-2dfs.png"))
    end

    println("✅ Wykres wygenerowany dla $pulsar_name.")
end






function check_period_column(fits_path::String)
    f = FITS(fits_path)
    try
        hdr = read_header(f[4])
        tf = get(hdr, "TFIELDS", 0)
        for i in 1:tf
            key = "TTYPE$i"
            if haskey(hdr, key) && hdr[key] == "PERIOD"
                # Odczytaj kolumnę PERIOD
                period_data = read(f[4], "PERIOD")
                println("Wartości w kolumnie PERIOD:")
                println(period_data)
                close(f)
                return true
            end
        end
        close(f)
        return false
    catch e
        close(f)
        rethrow(e)
    end
end





 function test(outdir)
        d = Data.load_ascii("input/1.txt")
        Plot.single(d, outdir; darkness=0.3, number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.average(d, outdir; number=256, bin_st=400, bin_end=600, start=1, name_mod="1", show_=true)
        Plot.lrfs(d, outdir; darkness=0.1, start=1, name_mod="1", bin_st=500, bin_end=530, show_=true)
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
        Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        folded = Tools.p3fold(data, 20, 40)
        Plot.single(folded, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319_p3fold", show_=true)
    end


    function process_psrdata(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="test", show_=true, repeat_num=4)
    end
    
    
    
    function main()
        # output directory for local run
        localout = "output"
        # output directory for VPM
        vpmout = "/home/psr/output/"
        indir = "/home/psr/data/"

        plot_2dfs("/home/psr/output/J0738-4042", "J0738-4042", show_plot=true)
        #process_psrdata("/home/psr/data/new/J0738-4042/2019-11-18-21:45:21/", vpmout)
        #period_value = 1.6039332673478421  # Twoja znana wartość

        #Plot_ostateczny_simple("/home/psr/output", "J1919+0134", period_value; show_plot=true)

        #process_psrdata("/home/psr/data/new/J0924-5302/2019-12-15-01:28:41/", vpmout)

        #Plot_2dfs_zmiany("/home/psr/output", "J1919+0134", show_plot=true)
        #Plot_2DFS_from_pspec("/home/psr/output", "J1919+0134", show_plot=true)
        #Plot_2dfs_simple("/home/psr/output", "J1919+0134", show_plot=true)
        #Check_2dfs_file("/home/psr/output", "J1919+0134")
        #detailed_check_fits("/home/psr/output", "J1919+0134")

    
    end

    function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table! s begin
            "--indir", "-i"
                help = "input directory"
                default = "input"
            "--outdir", "-o"
                help = "output directory"
                default = "output"
            "--plot", "-p"
                help = "plots to create"
                default = []
                nargs = '*'
        end
        return parse_args(s)
    end

end # module

SpaTs.main()

println("Bye")
