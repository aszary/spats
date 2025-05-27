module SpaTs
    using ArgParse
    using Glob
    using JSON
    #=using CairoMakie
    using FileIO
    using PDFIO
    using ImageMagick
=#
    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


   
   #= function J0820Mac(outdir)
        Data.convert_psrfit_ascii("/home/psr/data/J0820-1350/2020-01-11-01:05:56_00768-01055.spCF" , "/home/psr/output/2.txt")
        data = Data.load_ascii("/home/psr/output/2.txt")



        Plot.single(data, outdir, darkness=0.5, bin_st=1 , bin_end=1024, number=nothing, name_mod="J0820Mac", show_=true)
        Plot.lrfs(data, outdir, darkness=0.1, start=1,  bin_st=1, bin_end=1024, name_mod="J0820Mac", change_fftphase=false, show_=true)
        Plot.average(data, outdir, bin_st=1, bin_end=1024, number=nothing, name_mod="J0820Mac", show_=true)
    end
=#





    # TODO użyj JSON, czytanie do słownika d
    # TODO do poprawy, przeczytać pmod_output.txt
    # TODO zakres pulsów do rysowania (pobrany z JSON)
   

#=

function J0820Mac(outdir)
    params_file = "params.json"
    output_txt = joinpath(outdir, "J0924-5814.txt")
    debased_file = replace(output_txt, ".txt" => ".debase.gg")

    # === Load or initialize JSON parameter library ===
    p = isfile(params_file) ? Tools.read_params(params_file) : Tools.default_params(params_file)

    # === Ensure output directory exists ===
    if !isdir(outdir)
        mkpath(outdir)
        println("Created output directory: ", outdir)
    end

    # === Convert PSRFIT file to ASCII format ===
    input_file = "/home/psr/data/new/J0924-5814/2019-10-19-05:58:44/2019-10-19-05:58:44_00512-00767.spCF"
    
    if !isfile(output_txt)
        println("Converting PSRFIT to ASCII...")
        Data.convert_psrfit_ascii(input_file, output_txt)
    else
        println("Found existing ASCII file.")
    end

    # === Debase the data using pmod ===
    if !isfile(debased_file)
        println("Running pmod to debase...")
        run(pipeline(`pmod -device "/xw" -debase $output_txt`, `tee pmod_output.txt`))

        # Read the captured output from pmod
        output = read("pmod_output.txt", String)
       # Cleanup the log file

        # Extract onpulse range from the output
        m = match(r"-onpulse '(\d+) (\d+)'", output)
        if !isnothing(m)
            bin_st, bin_end = parse.(Int, m.captures)
            # Ensure the onpulse region length is even
            region_length = bin_end - bin_st + 1
            if region_length % 2 != 0
                println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                bin_end -= 1
                println("Adjusted onpulse range: $bin_st to $bin_end")
            end
            println("Found onpulse range: $bin_st to $bin_end")
            p["bin_st"] = bin_st
            p["bin_end"] = bin_end
        end
    else
        println("Using existing debased file.")
    end

    # === Save updated params back to JSON ===
    Tools.save_params(params_file, p)
    println("Parameters updated and saved to $params_file")

    # === Load data and plot ===
    data = Data.load_ascii(output_txt)

    Plot.single(data, outdir;darkness=0.5,bin_st=p["bin_st"],bin_end=p["bin_end"],start=p["pulse_start"],number=p["number"],name_mod="J0820Mac",show_=true)
    Plot.lrfs(data, outdir;darkness=0.1,start=p["pulse_start"],bin_st=p["bin_st"],bin_end=p["bin_end"],name_mod="J0820Mac",change_fftphase=false,show_=true)
    Plot.average(data, outdir;bin_st=p["bin_st"],bin_end=p["bin_end"],number=nothing,name_mod="J0820Mac",show_=true)
end

=#


















    
    function process_psrfit_files(base_dir::String, output_dir::String; name_mod::Union{String, Nothing}=nothing)
        # Step 1: Extract base directory 
        base_name = basename(base_dir)
       
   

        # Step 2: Create output subdirectory
        output_subdir = joinpath(output_dir, base_name)
        println(output_subdir)
        if !isdir(output_subdir)
            mkpath(output_subdir)
            println("Created output directory: ", output_subdir)
        end

        # Step 3 create JSON file
        params_file = joinpath(output_subdir, "p.json")



        #=
        #Step 3.5 Check if already processed
        if isdir(output_subdir) && isfile(joinpath(output_subdir, "p.json"))
            println("Skipping already processed catalogue: $base_name")
            return
        end
        =#




        println(params_file)
        
        p = isfile(params_file) ? Tools.read_params(params_file) : Tools.default_params(params_file)


        # Step 4: Find second-level catalogue
        second_catalogue = joinpath(base_dir, readdir(base_dir)[1])
        println("Second catalogue found: ", second_catalogue)
    
        # Step 5: Find .spCF files
        spcf_files = filter(f -> occursin("spCF", f), readdir(second_catalogue, join=true))
        Data.sort!(spcf_files)

        converted_txt_files = String[]

        #step 6: combining .spCF into one 
        output_file = joinpath(output_subdir, "converted.spCF")
        file_names = [joinpath(base_name, file) for file in spcf_files]

        Data.sort!(spcf_files)

        run(pipeline(`psradd $file_names -o $output_file`, stderr="errs.txt"))
        out_txt=replace(output_file ,".spCF" => ".txt")

        # Step 7: Convert spCF -> ascii
        Data.convert_psrfit_ascii(output_file, out_txt)

        rm(output_file) #cleanup .spCF file


        

        # Step 8: Debase the combined ASCII file using pmod
       # debased_file = replace(out_txt, ".txt" => ".debase.gg")
        #run(pipeline(`pmod -device "/xw" -debase $out_txt`, `tee pmod_output.txt`))
    
        # Read captured output and cleanup
        # output = read("pmod_output.txt", String)
        # rm("pmod_output.txt")  # Cleanup 

        # Step 9: Load combined data
        combined_data = Data.load_ascii(out_txt)

        #output_debase = read(debased_file, String)

        # Extract onpulse values
        #m = match(r"-onpulse\s+'(\d+)\s+(\d+)'", output_debase)
        #=p["bin_st"] = bin_st
        p["bin_end"] = bin_end=#

        if isfile(joinpath(output_subdir, "p.json"))
          println("File exitst")
          
        else
         m = match(r"-onpulse '(\d+) (\d+)'", output)
            if !isnothing(m)
                bin_st, bin_end = parse.(Int, m.captures)
              #  p["bin_st"] = bin_st
               #p["bin_end"] = bin_end

                # Ensure onpulse region length is even
                region_length = bin_end - bin_st + 1
                if region_length % 2 != 0
                  println("Warning: Onpulse region length ($region_length) is not even. Adjusting bin_end to make it even.")
                   bin_end -= 1
                   println("Adjusted onpulse range: $bin_st to $bin_end")
                end
               p["bin_st"] = bin_st
               p["bin_end"] = bin_end
               println("Found onpulse range: $bin_st to $bin_end")
            
             end
             Tools.save_params(params_file, p)
              println("Parameters updated and saved to $params_file")
        end
           #= Tools.save_params(params_file, p)
            println("Parameters updated and saved to $params_file")
            =#
    
        # Step 10: Plot
        Plot.single(combined_data, output_subdir; darkness=0.5, bin_st= p["bin_st"], bin_end= p["bin_end"], start= p["pulse_start"], number= (p["pulse_end"] - p["pulse_start"]), name_mod=name_mod, show_=false)

        Plot.lrfs(combined_data, output_subdir; darkness=0.1, start= p["pulse_start"], bin_st= p["bin_st"], bin_end= p["bin_end"], name_mod=name_mod, change_fftphase=false, show_=false)

        #Plot.average(combined_data, output_subdir; bin_st=p["bin_st"], bin_end=p["bin_end"], number=(p["pulse_end"]-p["pulse_start"]), name_mod=name_mod, show_=false)
        
        #=
        Plot.single(combined_data, output_subdir, darkness=0.5, bin_st=bin_st, bin_end=bin_end, number=nothing, name_mod=name_mod, show_=false)
        Plot.lrfs(combined_data, output_subdir, darkness=0.1, start=1, bin_st=bin_st, bin_end=bin_end, name_mod=name_mod, change_fftphase=false, show_=false)
        Plot.average(combined_data, output_subdir,bin_st=bin_st, bin_end=bin_end, number=nothing, name_mod=name_mod, show_=false)
        =#
    end
    


    function process_all_catalogues(output_dir::String, base_root::String="/home/psr/data/new")
        # Get all subdirectories in base_root
        catalogues = filter(isdir, readdir(base_root, join=true))
    
        if isempty(catalogues)
            println("No catalogues found in $base_root. Exiting...")
            return
        end
    
        for catalogue in catalogues
            base_name = basename(catalogue)  # Extract directory name
            println("Processing catalogue: ", base_name)
            process_psrfit_files(catalogue, output_dir, name_mod=base_name * "Mac")
        end
    end
    #=
    function combine_pdfs(output_dir::String)
    # Find all .pdf files in subdirectories
    pdf_paths = String[]
    for (root, _, files) in walkdir(output_dir)
        for file in files
            if endswith(file, "single.pdf")
                push!(pdf_paths, joinpath(root, file))
            end
        end
    end

    # Sort for consistent ordering
    sort!(pdf_paths)
    combined_pdf_path = joinpath(output_dir, "singles.pdf")

    if isempty(pdf_paths)
        println("No PDF files found in subdirectories of $output_dir.")
        return
    end

    # Use pdfunite to combine PDFs

    #Zapisz to .txt i daj paths do wszystkich 
    try
        println(`pdfunite $(pdf_paths...) $combined_pdf_path`)
       # run(`pdfunite $(pdf_paths...) $combined_pdf_path`)
        #println(`pdfunite $(pdf_paths...) $combined_pdf_path`)
        #println("Combined PDF created at: $combined_pdf_path")
    catch e
        println("Error combining PDFs: ", e)
    end
end
=#



#=
function combine_pdfs(output_dir::String)
    # Find all .pdf files in subdirectories
    pdf_paths = String[]
    for (root, _, files) in walkdir(output_dir)
        for file in files
            if endswith(file, "single.pdf")
                push!(pdf_paths, joinpath(root, file))
            end
        end
    end

    # Sort for consistent ordering
    sort!(pdf_paths)

    if isempty(pdf_paths)
        println("No PDF files found in subdirectories of $output_dir.")
        return
    end

    # Extract pulsar name
    pulsar_name = splitpath(output_dir)[end]

    # True combined path (where you actually want to save)
    combined_pdf_path = joinpath("/home/aszary/output/Maciej", pulsar_name * ".pdf ")

    # Create a version of the paths just for printing (fake Maciej-style paths)
    fake_pdf_paths = [replace(path, "/home/psr/output" => "/home/aszary/output/Maciej") for path in pdf_paths]

    try
        println("Simulated command:")
        println("pdfunite " * join(fake_pdf_paths, " ") * " " * combined_pdf_path)
        # Actual execution (if you want)
        # run(`pdfunite $(pdf_paths...) $combined_pdf_path`)
    catch e
        println("Error combining PDFs: ", e)
    end
end
=#

#=
function combine_pdfs(output_dir::String)
    # Find all .pdf files in subdirectories
    pdf_paths = String[]
    for (root, _, files) in walkdir(output_dir)
        for file in files
            if endswith(file, "single.pdf")
                push!(pdf_paths, joinpath(root, file))
            end
        end
    end



function render_pdf_to_image(path::String)
    img = ImageMagick.load(path)  # Will only work if ImageMagick CLI is working
    return convert(Matrix{RGB{N0f8}}, img)
end


    # Sort for consistent ordering
    sort!(pdf_paths)

    if isempty(pdf_paths)
        println("No PDF files found in subdirectories of $output_dir.")
        return
    end

    # Extract pulsar name
    pulsar_name = splitpath(output_dir)[end]

    # True combined path (where you actually want to save)
    combined_pdf_path = joinpath("/home/aszary/output/Maciej", pulsar_name * ".pdf")

    # Create a version of the paths just for printing (fake Maciej-style paths)
    fake_pdf_paths = [replace(path, "/home/psr/output" => "/home/aszary/output/Maciej") for path in pdf_paths]

    # NEW: Plotting and rendering
    rendered_pages = String[]
    fig = nothing
    page_count = 0


 # Rendering PDFs to pages with Makie
    rendered_pages = String[]
    fig = nothing
    page_count = 0

    for (i, path) in enumerate(pdf_paths)
        if (i - 1) % 4 == 0
            if fig !== nothing
                save("/tmp/page_$page_count.pdf", fig)
                push!(rendered_pages, "/tmp/page_$page_count.pdf")
                page_count += 1
            end
            fig = Figure(size = (800, 800))
        end

        row = div((i - 1), 2) % 2 + 1
        col = (i - 1) % 2 + 1
        ax = Axis(fig[row, col])

        try
            #image_data = load(path)
            image_data = render_pdf_to_image(path)
  # First page preview
            image!(ax, image_data)
            ax.title = splitpath(path)[end]
        catch e
            ax.title = "Error loading: $(basename(path))"
        end
    end

    # Save the last page
    if fig !== nothing
        save("/tmp/page_$page_count.pdf", fig)
        push!(rendered_pages, "/tmp/page_$page_count.pdf")
    end

    # Combine all rendered pages into one
    try
        println("Simulated command:")
        println("pdfunite " * join(rendered_pages, " ") * " " * combined_pdf_path)
        run(`pdfunite $(rendered_pages...) $combined_pdf_path`)
    catch e
        println("Error combining PDFs: ", e)
    end
end
=#








       # Run processing for all catalogues
    function J0034Mac(output_dir)
        process_all_catalogues(output_dir, "/home/psr/data/new")
       # combine_pdfs("/home/psr/output")
    end 


    function main()
        # output directory for local run
        localout = "output"
        # output directory for VPM
        vpmout = "/home/psr/output/"



        #test(vpmout)
        #J0820Mac(vpmout)
        J0034Mac(vpmout)
        #mkieth()
        #J1651()
        #J1705()
        #B0320()
        #J1750_remote()
        #J1750_local()
        #J1750_paper()
        #J1750_paper2()
        #J1750_average()
        #J1750_modeled()
        #J1750_calculations()

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