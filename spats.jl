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
    data = data'

    # Skalowanie jak w pplot.c
    data_min = minimum(data)
    data_max = maximum(data)
    scale = 1.0
    scale2 = 0.0
    new_min = data_min + (data_max - data_min) * scale2
    new_max = data_min + (data_max - data_min) / scale
    if new_max <= new_min
        new_min, new_max = data_min, data_max
    end

    profile_P3 = sum(data, dims=2)[:,1]
    profile_P2 = sum(data, dims=1)[1,:]

    fig = figure(figsize=(8, 8))
    gs = matplotlib[:gridspec][:GridSpec](2, 2,
        width_ratios=[1, 4], height_ratios=[1, 4],
        wspace=0.05, hspace=0.05)

    axMain = fig.add_subplot(gs[1, 1])
    axLeft = fig.add_subplot(gs[1, 0], sharey=axMain)
    axTop = fig.add_subplot(gs[0, 1], sharex=axMain)

    axTop.xaxis.set_tick_params(labelbottom=false)
    axLeft.yaxis.set_tick_params(labelleft=false)

    p2_min = -n_bins / 2
    p2_max = n_bins / 2 - 1
    p3_min = 0.0
    p3_max = 0.5

    x_vals = LinRange(p2_min, p2_max, n_bins)
    y_vals = LinRange(p3_min, p3_max, n_subints)

    im = axMain.imshow(data;
        cmap="gray_r",
        norm=matplotlib[:colors][:Normalize](vmin=new_min, vmax=new_max),
        origin="lower",
        extent=[p2_min, p2_max, p3_min, p3_max],
        aspect="auto"
    )

    colorbar(im, ax=axMain, label="Power", shrink=0.9)

    #axLeft.plot(profile_P3, y_vals, color="black", lw=1.5)
    #axLeft.set_xlabel("Power")

    #axTop.fill_between(x_vals, 0, profile_P2, facecolor="lightgray", edgecolor="black", alpha=0.5)
    #axTop.plot(x_vals, profile_P2, color="black", lw=1.5)
    #axTop.set_ylabel("Power")

    axMain.set_xlabel("fluctuation frequency (cycles/period)")
    axMain.set_ylabel("fluctuation frequency (cycles/period)")
    fig.suptitle("2DFS – $pulsar_name")

    savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
    savefig(savepath)
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
        
        Plot_2dfs_zmiany("/home/psr/output", "J1919+0134", show_plot=true)


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
