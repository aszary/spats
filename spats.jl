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
                            println("  Found 2D numeric column '$name' in HDU $i")
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
            println(" No suitable 2D data found in FITS file.")
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
        println(" 2DFS plot saved to: $savepath")


        if show_plot
            show()
        else
            close(fig)
        end


    catch e
        println(" Error handling FITS file: $e")
        if f !== nothing
            close(f)
        end
    end
end



function read_2dfs_file(filename::String)
    println("Reading FITS file: $filename")
    
    # Check if file exists
    if !isfile(filename)
        error("File does not exist: $filename")
    end
    
    # Open the FITS file
    fits = FITS(filename)
    
    try
        # Print header information for debugging
        header = read_header(fits[1])
        println("FITS Header Information:")
        
        # Try to get dimensions from header
        NrBins = get(header, "NAXIS1", nothing)
        NrSubints = get(header, "NAXIS2", nothing)
        
        if NrBins === nothing || NrSubints === nothing
            println("Warning: Could not find NAXIS1/NAXIS2 in header")
            println("Header contents:")
            for (key, value) in header
                println("$key: $value")
            end
            
            # Try to infer dimensions from data
            data = read(fits[1])
            if data !== nothing
                println("Data shape: ", size(data))
                NrBins, NrSubints = size(data)
            else
                error("Could not read data from FITS file")
            end
        else
            println("NAXIS1: $NrBins")
            println("NAXIS2: $NrSubints")
            
            # Read the data
            data = read(fits[1])
        end
        
        # Try to get other parameters from header
        f2_min = get(header, "F2_MIN", nothing)
        f2_max = get(header, "F2_MAX", nothing)
        f3_min = get(header, "F3_MIN", nothing)
        f3_max = get(header, "F3_MAX", nothing)
        
        if f2_min === nothing || f2_max === nothing || f3_min === nothing || f3_max === nothing
            println("Warning: Could not find some required parameters in header")
        end
        
        return (data, NrBins, NrSubints, f2_min, f2_max, f3_min, f3_max)
    catch e
        println("Error reading FITS file:")
        println(e)
        println("Header contents:")
        if haskey(locals(), :header)
            for (key, value) in header
                println("$key: $value")
            end
        end
        rethrow()
    finally
        # Close the file
        close(fits)
    end
end

    
    
    function main()
        # output directory for local run
        localout = "output"
        # output directory for VPM
        vpmout = "/home/psr/output/"
        indir = "/home/psr/data/"

        data, NrBins, NrSubints, f2_min, f2_max, f3_min, f3_max = read_2dfs_file(vpmout * "/pulsar.debase.1.2dfs") 

        println("Number of bins: $NrBins")
        println("Number of sub-integrations: $NrSubints")
        println("f2_min: $f2_min")
        println("f2_max: $f2_max")
        println("f3_min: $f3_min")
        println("f3_max: $f3_max")

        #plot_2dfs("/home/psr/output", "J1919+0134", show_plot=true)



    end


end # module

SpaTs.main()

println("Bye")
