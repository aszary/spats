module SpaTs
    using ArgParse
    using Glob
    using JSON
    using FITSIO
    using PyPlot
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
            cmap="gray",  # Dodano przecinek
            vmin=0,  # Minimalna wartość skali kolorów
            vmax=5000  # Maksymalna wartość skali kolorów
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










   

    
    
    
    function main()
        # output directory for local run
        localout = "output"
        # output directory for VPM
        vpmout = "/home/psr/output/"
        indir = "/home/psr/data/"

        #plot_2dfs("/home/psr/output", "J1919+0134", show_plot=true)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        
        #Plot_2dfs_zmiany("/home/psr/output", "J1919+0134", show_plot=true)
        #Plot_2DFS_from_pspec("/home/psr/output", "J1919+0134", show_plot=true)
        Plot_2dfs_zmiany_pplot("/home/psr/output", "J1919+0134", show_plot=true)

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
