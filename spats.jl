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



function analyze_fits_file(filename::String)
    println("Analyzing FITS file: $filename")
    
    # Check if file exists
    if !isfile(filename)
        println("Error: File does not exist")
        return
    end
    
    # Open the FITS file
    fits = FITS(filename)
    
    try
        # Get number of HDUs
        println("\nNumber of HDUs: ", length(fits))
        
        # Analyze each HDU
        for (i, hdu) in enumerate(fits)
            println("\n=== HDU $i ===")
            
            # Print header information
            header = read_header(hdu)
            println("\nHeader Keywords:")
            for (key, value) in header
                println("$key: $value")
            end
            
            # Print data information if available
            try
                data = read(hdu)
                println("\nData Information:")
                println("Data type: ", typeof(data))
                println("Data size: ", size(data))
                println("Data shape: ", size(data))
                println("Data element type: ", eltype(data))
                
                # Print some statistics if numeric data
                if isa(data, AbstractArray) && eltype(data) <: Number
                    println("\nData Statistics:")
                    println("Min value: ", minimum(data))
                    println("Max value: ", maximum(data))
                    println("Mean value: ", mean(data))
                    println("Standard deviation: ", std(data))
                end
            catch e
                println("\nNo data in this HDU or error reading data:")
                println(e)
            end
        end
    catch e
        println("\nError analyzing FITS file:")
        println(e)
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

        analyze_fits_file(vpmout * "/pulsar.debase.1.2dfs") 

        #plot_2dfs("/home/psr/output", "J1919+0134", show_plot=true)



    end


end # module

SpaTs.main()

println("Bye")
