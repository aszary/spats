"""
Renders and saves a 2DFS plot from the file pulsar.debase.1.2dfs using PyPlot.


Arguments:
- outdir: The output directory where the pulsar data is stored.
- pulsar_name: The name of the pulsar to create the plot for.
- show_plot: A boolean flag to decide whether to display the plot (default: true).
"""
function plot_2dfs(outdir::String, pulsar_name::String; show_plot::Bool=true)


    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")


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


        savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
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





function plot_2dfs_zmiany(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("File does not exist: $filepath")
        return
    end

    println("Inspecting FITS file: $filepath")

    f = nothing
    data = nothing

    try
        f = FITS(filepath)
        hdu = f[4]  # Directly accessing HDU 4, since we know the structure
        data = read(hdu, "DATA")  # Reading the "DATA" column
        close(f)

        if data === nothing
            println("❌ No suitable 2D data found.")
            return
        end

        # === Extract size and coordinate ranges ===
        n_p3, n_p2 = size(data)
        p3_range = range(0, stop=0.5, length=n_p3)
        p2_range = range(160, stop=200, length=n_p2)

        # === Compute column sums for left panel ===
        col_sums = sum(data, dims=1)[1, :]  # Collapse each column into a single value

        # === Create side-by-side plots ===
        fig, axs = subplots(1, 2, figsize=(10, 5), width_ratios=[1, 4])

        # Left panel: column sum profile (rotated vertically)
        axs[1].plot(-col_sums, p2_range)
        axs[1].invert_xaxis()  # Flip horizontally so it appears on the left
        axs[1].set_ylabel("Pulse longitude (deg)")
        axs[1].set_xlabel("Power sum")
        axs[1].grid(true)

        # Right panel: original 2DFS image
        im = axs[2].imshow(data';
            extent=[160, 200, 0, 0.5],
            origin="lower",
            aspect="auto",
            cmap="gray",
            vmin=0,
            vmax=5000
        )

        axs[2].set_xlabel("Pulse longitude (deg)")
        axs[2].set_ylabel("Fluctuation frequency (P/P3)")
        axs[2].set_title("2DFS – $pulsar_name")
        colorbar(im, ax=axs[2], label="Power")

        # Save figure
        savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
        savefig(savepath)
        println("✅ 2DFS plot with left panel saved to: $savepath")

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



function plot_2dfs_zmiany2(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("File does not exist: $filepath")
        return
    end

    println("Inspecting FITS file: $filepath")

    f = nothing
    data = nothing

    try
        f = FITS(filepath)
        hdu = f[4]  # Directly access HDU 4
        data = read(hdu, "DATA")
        close(f)

        if data === nothing
            println("❌ No suitable 2D data found.")
            return
        end

        # === Dimensions and axis ranges ===
        n_p3, n_p2 = size(data)
        p3_range = range(0, stop=0.5, length=n_p3)
        p2_range = range(160, stop=200, length=n_p2)

        # === Left panel: sum over columns instead of rows ===
        col_sums = sum(data, dims=1)[1, :]  # Collapse along rows → sum each column

        # === Subplots setup ===
        fig, axs = subplots(1, 2, figsize=(10, 5), width_ratios=[1, 4])

        # Left panel: column sum profile (vertical)
        axs[1].plot(col_sums, p2_range)
        axs[1].invert_xaxis()  # Flip to keep panel visually on the left
        axs[1].set_xscale("log")
        axs[1].set_ylabel("Pulse longitude (deg)")
        axs[1].set_xlabel("Column sum (Power)")
        axs[1].grid(true)

        # Right panel: 2DFS image
        im = axs[2].imshow(data';
            extent=[160, 200, 0, 0.5],
            origin="lower",
            aspect="auto",
            cmap="gray",
            vmin=0,
            vmax=5000
        )

        axs[2].set_xlabel("Pulse longitude (deg)")
        axs[2].set_ylabel("Fluctuation frequency (P/P3)")
        axs[2].set_title("2DFS – $pulsar_name")
        colorbar(im, ax=axs[2], label="Power")

        # Save output
        savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
        savefig(savepath)
        println("✅ 2DFS plot with left column-sum panel saved to: $savepath")

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
    





function plot_lrfs(outdir::String, pulsar_name::String; show_plot::Bool=true)
filepath = joinpath(outdir, pulsar_name, "pulsar.debase.lrfs")

if !isfile(filepath)
    println("File does not exist: $filepath")
    return
end

println("Inspecting FITS file: $filepath")

f = nothing
data = nothing

try
    f = FITS(filepath)

    # Iteracja przez dostępne HDU w pliku FITS
    for (i, hdu) in enumerate(f)
        println("HDU $i:")
        println("  Type: ", typeof(hdu))

        try
            # Próbujemy znaleźć dane w różnych typach HDU
            if hdu isa FITSIO.ImageHDU
                img = read(hdu)
                if ndims(img) == 2  # Sprawdzamy, czy to jest obraz 2D
                    println("  -> Found 2D image of size ", size(img))
                    data = img
                    break
                else
                    println("  -> Not a 2D image (ndims=$(ndims(img)))")
                end

            elseif hdu isa FITSIO.TableHDU
                names = FITSIO.colnames(hdu)
                println("  Columns: ", names)

                # Przeglądamy kolumny, by znaleźć odpowiednie dane
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

        # Jeśli dane zostały znalezione, przerywamy iterację
        if data !== nothing
            break
        end
    end

    # Jeśli nie znaleziono danych, zwracamy błąd
    if data === nothing
        println("❌ No suitable LRFS data found.")
        close(f)
        return
    end

    # Zamykanie pliku FITS
    close(f)

    # === Rozmiar danych ===
    n_freq, n_long = size(data)

    # Automatyczne dopasowanie zakresu na podstawie rozmiaru danych
    freq_range = range(0, stop=0.5, length=n_freq)  # Frequencies based on data size
    long_range = range(0, stop=360, length=n_long)  # Pulse longitude range based on data size

    # === Tworzenie wykresu ===
    fig, ax = subplots()
    im = ax.imshow(data';
        extent=[0, 360, 0, 0.5],  # [x_min, x_max, y_min, y_max]
        origin="lower",
        aspect="auto",
        cmap="gray",
        vmin=0,
        vmax=maximum(data)  # Dynamic range based on data values
    )

    ax.set_xlabel("Pulse longitude (deg)")
    ax.set_ylabel("Fluctuation frequency (P/P3)")
    ax.set_title("LRFS – $pulsar_name")
    colorbar(im, ax=ax, label="Power")

    # === Zapis wykresu ===
    savepath = joinpath(outdir, pulsar_name, "lrfs_" * pulsar_name * ".png")
    savefig(savepath)
    println("✅ LRFS plot saved to: $savepath")

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














function inspect_fits(filepath::String)

    if !isfile(filepath)
        println("File does not exist: $filepath")
        return
    end

    println("Inspecting FITS file: $filepath")
    
    try
        f = FITS(filepath)

        for (i, hdu) in enumerate(f)
            println("HDU $i:")
            println("  Type: ", typeof(hdu))
            try
                if hdu isa FITSIO.ImageHDU
                    d = read(hdu)
                    println("  -> ImageHDU with dims: ", size(d))
                elseif hdu isa FITSIO.TableHDU
                    println("  -> TableHDU with columns: ", colnames(hdu))
                else
                    println("  -> Unknown HDU type.")
                end
            catch e
                println("  -> Failed to read HDU $i: $e")
            end
        end

        close(f)
    catch e
        println("Error opening FITS file: $e")
    end
end















function plot_lrfs22(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.lrfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading LRFS from: $filepath")
    
    f = nothing
    data = nothing

    try
        f = FITS(filepath)
        hdu = f[4]  # Możliwe że 2 albo inny — trzeba sprawdzić!
        data = read(hdu, "DATA")
        close(f)

        if data === nothing
            println("❌ No suitable 2D data found in LRFS.")
            return
        end

        n_freq, n_long = size(data)
        freq_range = range(0, stop=0.5, length=n_freq)
        long_range = range(0, stop=360*(n_long-1)/n_long, length=n_long)

        fig, ax = subplots()
        im = ax.imshow(data';
            extent=[0, 360, 0, 0.5],
            origin="lower",
            aspect="auto",
            cmap="gray",
            vmin=0,
            vmax=maximum(data)
        )

        ax.set_xlabel("Pulse phase (deg)")
        ax.set_ylabel("Fluctuation frequency (P/P3)")
        ax.set_title("LRFS – $pulsar_name")
        colorbar(im, ax=ax, label="Power")

        savepath = joinpath(outdir, pulsar_name, "lrfs_" * pulsar_name * ".png")
        savefig(savepath)
        println("✅ LRFS plot saved to: $savepath")

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





function plot2dfs333(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    # Sprawdzenie, czy plik istnieje
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    f = nothing
    data = nothing

    try
        # Otwarcie pliku FITS
        f = FITS(filepath)
        hdu = f[4]  # Możliwe że 2DFS jest w HDU 4, ale warto to zweryfikować
        data = read(hdu, "DATA")
        
        # Sprawdźmy wymiary danych (powinny być 2D)
        if ndims(data) != 2
            println("❌ Data is not 2D as expected.")
            return
        end

        close(f)

        # Jeśli dane są puste, zakończ działanie funkcji
        if data === nothing
            println("❌ No suitable 2D data found in 2DFS.")
            return
        end

        # Zakładając, że P3 to częstotliwość, a P2 to pulsacja, dostosujemy odpowiednie zakresy
        nP3, nP2 = size(data)
        P3 = LinRange(0.0, 0.5, nP3)  # Zakładając, że P3 jest od 0 do 0.5
        P2 = LinRange(-0.5, 0.5, nP2)  # Zakładając, że P2 jest od -0.5 do 0.5

        # Tworzenie wykresu z jednym panelem (środkowym)
        fig = figure(figsize=(8, 6))

        # Główny panel (2DFS)
        axMain = fig.add_subplot(111)

        # Rysowanie wykresu 2DFS
        img = axMain.imshow(data', origin="lower", aspect="auto", cmap="gray_r", 
                            extent=[first(P2), last(P2), first(P3), last(P3)])
        colorbar(img, ax=axMain, label="Power", shrink=0.9)

        axMain.set_xlabel("P2 (cycles per period)")
        axMain.set_ylabel("P3 (cycles per period)")
        fig.suptitle("2DFS – $pulsar_name")

        # Zapis do pliku PNG
        savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
        savefig(savepath, dpi=300)
        println("✅ 2DFS plot saved to: $savepath")

        # Wyświetlenie wykresu, jeśli show_plot jest ustawione na true
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










function inspect_fits22(filepath::String)
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("🔎 Inspecting FITS file: $filepath")

    f = FITS(filepath)
    for (i, hdu) in enumerate(f)
        println("\n📂 HDU $i:")
        println("  Type: ", typeof(hdu))

        try
            if hdu isa FITSIO.ImageHDU
                img = read(hdu)
                println("  - Image dimensions: ", size(img))
                println("  - Image element type: ", eltype(img))

            elseif hdu isa FITSIO.TableHDU
                names = FITSIO.colnames(hdu)
                println("  - Table columns: ", names)

                for name in names
                    col_data = read(hdu, name)
                    println("    Column '$name': type ", eltype(col_data), ", size ", size(col_data))
                end

            else
                println("  - Unknown HDU type.")
            end
        catch e
            println("  ⚠️ Error reading HDU $i: $e")
        end
    end
    close(f)
end




    

function plot_correct_2dfs(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    try
        # Open FITS
        f = FITS(filepath, "r")
        hdu = f[4]
        data_raw = read(hdu, "DATA")
        
        # Check dimensions
        if ndims(data_raw) != 2
            println("❌ Data is not 2D as expected.")
            close(f)
            return
        end

        # Try to read scaling if available
        header = read_header(hdu)
        dat_scl = tryparse(Float64, get(header, "DAT_SCL", "1.0"))
        dat_offs = tryparse(Float64, get(header, "DAT_OFFS", "0.0"))
        period = tryparse(Float64, get(header, "PERIOD", "1.0"))  # seconds

        data = data_raw .* dat_scl .+ dat_offs
        
        close(f)

        # Prepare for 2D FFT
        n_pulses, n_bins = size(data)
        println("ℹ️ Data shape: $n_pulses pulses × $n_bins bins")

        # 2D FFT over pulse axis (time)
        F = fft(data, 1)  # FFT along pulses
        F = fftshift(F, 1)  # Center zero frequency vertically
        power = abs.(F).^2

        # Create axis scales
        freq = fftshift(fftfreq(n_pulses, 1))  # P/P3 axis
        pulse_long = LinRange(0, 360, n_bins + 1)[1:end-1]  # Pulse longitude in degrees

        # Plot using PyPlot
        figure(figsize=(8, 6))
        ax = gca()
        im = ax.imshow(power, aspect="auto", cmap="Greys", extent=[pulse_long[1], pulse_long[end], freq[1], freq[end]], vmin=0, vmax=quantile(vec(power), 0.95))
        ax.set_xlabel("Pulse longitude [deg]")
        ax.set_ylabel("Fluctuation frequency [cpp]")
        ax.set_title("2DFS – $pulsar_name")
        ax.set_ylim(0, 0.5)
        colorbar(im, ax=ax, label="Power")

        savepath = joinpath(outdir, pulsar_name, "corrected_2dfs_" * pulsar_name * ".png")
        savefig(savepath)

        println("✅ 2DFS saved to: $savepath")
        
        if show_plot
            display()
        end

    catch e
        println("❌ Error while handling FITS file: $e")
    end
end



# Define Hanning window
function hanning(N::Int)
    if N <= 1
        return ones(N)
    else
        return 0.5 .- 0.5 * cos.(2pi .* (0:(N-1)) ./ (N-1))
    end
end



    function plot_2dfs_koncowy(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    try
        # Open FITS
        f = FITS(filepath, "r")
        hdu = f[4]  # always 4th HDU for data

        # Read necessary fields
        data_raw = read(hdu, "DATA")
        dat_scl = tryparse(Float64, get(read_header(hdu), "DAT_SCL", "1.0"))
        dat_offs = tryparse(Float64, get(read_header(hdu), "DAT_OFFS", "0.0"))
        period = tryparse(Float64, get(read_header(hdu), "PERIOD", "1.0"))  # seconds

        close(f)

        # Scale data
        data = data_raw .* dat_scl .+ dat_offs

        # Check dimensions
        n_pulses, n_bins = size(data)
        println("ℹ️ Data shape: $n_pulses pulses × $n_bins bins")

        # Detrend: Remove mean from each pulse profile
        for i in 1:n_pulses
            data[i, :] .-= mean(data[i, :])
        end

        # Apply 2D window (Hanning window) to reduce spectral leakage
        win_pulses = hanning(n_pulses)
        win_bins = hanning(n_bins)
        window = win_pulses * win_bins'
        data_windowed = data .* window

        # Perform 2D FFT
        F = fftshift(fft(fft(data_windowed, 1), 2))

        # Power spectrum
        power = abs.(F).^2

        # Normalize (optional, matches the Song et al. paper appearance better)
        power ./= maximum(power)

        # Create axis: frequency scales
        pulse_freq = fftshift(fftfreq(n_pulses, 1))  # "cycles per pulse"
        bin_freq = fftshift(fftfreq(n_bins, 1))      # "cycles per bin"

        # Convert pulse longitude (bins) into degrees
        pulse_long = LinRange(0, 360, n_bins + 1)[1:end-1]

        # Plot 2DFS
        figure(figsize=(8, 6), dpi=300)
        ax = gca()
        extent = [pulse_long[1], pulse_long[end], pulse_freq[1], pulse_freq[end]]

        im = ax.imshow(power, aspect="auto", cmap="Greys", extent=extent, origin="lower",
                       vmin=0, vmax=quantile(vec(power), 0.95), interpolation="none")

        ax.set_xlabel("Pulse longitude [deg]")
        ax.set_ylabel("Fluctuation frequency (cycles per period, cpp)")
        ax.set_title("2DFS – $pulsar_name")
        ax.set_ylim(0, 0.5)  # Standard limit as in Song et al.
        ax.set_xlim(160, 200)


        colorbar(im, ax=ax, label="Normalized Power")

        # Save figure
        savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
        savefig(savepath)

        println("✅ 2DFS saved to: $savepath")

        if show_plot
            # Nic nie trzeba, wykres i tak się pokaże
        end

    catch e
        println("❌ Error while handling FITS file: $e")
    end
end






function plot_2dfs_koncowy2(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    try
        # Open FITS file
        f = FITS(filepath, "r")
        hdu = f[4]  # always 4th HDU for data

        # Read the raw data, check if it's an appropriate type, and force conversion to Float64
        data_raw = read(hdu, "DATA")  # Read raw data
        
        # Ensure the data is in correct floating point format (Float64 or Float32)
        if typeof(data_raw) <: AbstractArray && eltype(data_raw) <: Integer
            println("⚠️ Data is in integer format, converting to Float64...")
            data_raw = Float64.(data_raw)
        end

        close(f)

        # Check dimensions
        n_pulses, n_bins = size(data_raw)
        println("ℹ️ Data shape: $n_pulses pulses × $n_bins bins")

        # Detrend: Remove mean from each pulse profile
        for i in 1:n_pulses
            data_raw[i, :] .-= mean(data_raw[i, :])
        end

        # Apply 2D window (Hanning window)
        win_pulses = 0.5 .* (1 .- cos.(2π .* (0:n_pulses-1) ./ (n_pulses-1)))
        win_bins = 0.5 .* (1 .- cos.(2π .* (0:n_bins-1) ./ (n_bins-1)))
        window = win_pulses * win_bins'
        data_windowed = data_raw .* window

        # Perform 2D FFT
        F = fftshift(fft(fft(data_windowed, 1), 2))

        # Power spectrum
        power = abs.(F).^2

        # Normalize power
        power ./= maximum(power)

        # Frequency axes
        pulse_freq = fftshift(fftfreq(n_pulses, 1))  # cycles per pulse
        bin_freq = fftshift(fftfreq(n_bins, 1))      # cycles per bin

        # Longitude (bin) to degrees
        pulse_long = LinRange(0, 360, n_bins + 1)[1:end-1]

        # Plot
        figure(figsize=(8, 6))
        ax = gca()
        extent = [pulse_long[1], pulse_long[end], pulse_freq[1], pulse_freq[end]]

        im = ax.imshow(power, aspect="auto", cmap="Greys", extent=extent, origin="lower",
                       vmin=0, vmax=quantile(vec(power), 0.95), interpolation="none")

        ax.set_xlabel("Pulse longitude [deg]")
        ax.set_ylabel("Fluctuation frequency (cycles per period, cpp)")
        ax.set_title("2DFS – $pulsar_name")
        ax.set_xlim(160, 200)
        ax.set_ylim(0, 0.5)

        colorbar(im, ax=ax, label="Normalized Power")

        # Save figure
        savepath = joinpath(outdir, pulsar_name, "2dfs_" * pulsar_name * ".png")
        savefig(savepath, dpi=300)

        println("✅ 2DFS saved to: $savepath")

        if show_plot
            # Display automatically handled by the environment
        end

    catch e
        println("❌ Error while handling FITS file: $e")
    end
end




function plot_simple_2dfs_fixed(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    try
        f = FITS(filepath, "r")
        hdu = f[4]
        data = read(hdu, "DATA")
        close(f)

        n_pulses, n_bins = size(data)
        println("ℹ️ Data shape: $n_pulses pulses × $n_bins bins")

        # 2D FFT
        F = fftshift(fft(fft(data, 1), 2))
        power = abs.(F).^2

        # Osie
        pulse_long = LinRange(0, 360, n_bins + 1)[1:end-1]
        pulse_freq = fftshift(fftfreq(n_pulses, 1))

        # Wybieramy tylko fragment 160–200 stopni
        mask = (pulse_long .>= 160) .& (pulse_long .<= 200)
        power_crop = power[:, mask]
        pulse_long_crop = pulse_long[mask]

        extent = [pulse_long_crop[1], pulse_long_crop[end], pulse_freq[1], pulse_freq[end]]

        # Rysujemy tylko wycinek
        fig = PyPlot.figure(figsize=(8, 6))
        ax = PyPlot.gca()
        im = ax.imshow(power_crop, aspect="auto", cmap="Greys", extent=extent, origin="lower",
                       vmin=0, vmax=quantile(vec(power), 0.95), interpolation="none")

        ax.set_xlabel("Pulse longitude [deg]")
        ax.set_ylabel("Fluctuation frequency (P/P3)")
        ax.set_title("2DFS – $pulsar_name")
        ax.set_xlim(160, 200)
        ax.set_ylim(0, 0.5)

        PyPlot.colorbar(im, ax=ax, label="Power")

        savepath = joinpath(outdir, pulsar_name, "simple_2dfs_" * pulsar_name * ".png")
        PyPlot.savefig(savepath, dpi=300)

        println("✅ 2DFS saved to: $savepath")

        if show_plot
            PyPlot.show()
        end

    catch e
        println("❌ Error while handling FITS file: $e")
    end
end





function plot_2dfs_pulse_longitude(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    try
        f = FITS(filepath)
        hdu = f[4]
        data = read(hdu, "DATA")
        close(f)

        if ndims(data) != 2
            println("❌ Data is not 2D.")
            return
        end

        n_pulses, n_bins = size(data)
        println("ℹ️ Data shape: $n_pulses pulses × $n_bins bins")

        # Perform FFT
        F = fftshift(fft(fft(data, 1), 2))
        power = abs.(F).^2

        # Axes
        pulse_freq = fftshift(fftfreq(n_pulses, 1))     # fluctuation frequency (P/P3)
        bin_freq = fftshift(fftfreq(n_bins, 1))         # fluctuation frequency (P/P2)
        pulse_long = LinRange(0, 360, n_bins + 1)[1:end-1]  # Pulse longitude in degrees

        # Crop longitude 160–200 deg
        mask = (pulse_long .>= 160) .& (pulse_long .<= 200)
        power_crop = power[:, mask]
        pulse_long_crop = pulse_long[mask]

        extent = [pulse_long_crop[1], pulse_long_crop[end], pulse_freq[1], pulse_freq[end]]

        fig = PyPlot.figure(figsize=(8, 6))
        ax = PyPlot.gca()
        im = ax.imshow(power_crop, aspect="auto", cmap="Greys", extent=extent, origin="lower",
                       vmin=0, vmax=quantile(vec(power), 0.95), interpolation="none")

        ax.set_xlabel("Pulse longitude [deg]")
        ax.set_ylabel("Fluctuation frequency (P/P3)")
        ax.set_xlim(160, 200)
        ax.set_ylim(0, 0.5)

        PyPlot.colorbar(im, ax=ax, label="Power")
        ax.set_title("2DFS – $pulsar_name")

        savepath = joinpath(outdir, pulsar_name, "2dfs_pulse_long_" * pulsar_name * ".png")
        PyPlot.savefig(savepath, dpi=300)
        println("✅ 2DFS pulse longitude plot saved to: $savepath")

        if show_plot
            PyPlot.show()
        else
            PyPlot.close(fig)
        end

    catch e
        println("❌ Error while handling FITS file: $e")
    end
end


function plot2dfs_p2_frequency(outdir::String, pulsar_name::String; show_plot::Bool=true)
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("✅ Reading 2DFS from: $filepath")

    f = nothing
    data = nothing

    try
        f = FITS(filepath)
        hdu = f[4]
        data = read(hdu, "DATA")
        close(f)

        if ndims(data) != 2
            println("❌ Data is not 2D as expected.")
            return
        end

        n_pulses, n_bins = size(data)
        println("ℹ️ Data shape: $n_pulses pulses × $n_bins bins")

        # FFT
        F = fftshift(fft(fft(data, 1), 2))
        power = abs.(F).^2

        # Axes
        pulse_freq = fftshift(fftfreq(n_pulses, 1))  # cycles per pulse (P/P3)
        bin_freq = fftshift(fftfreq(n_bins, 1))      # cycles per bin (P/P2)

        # Fluctuation freq (P/P2) przeskalowana na zakres -100 do 100
        bin_freq_scaled = bin_freq .* 100

        extent = [bin_freq_scaled[1], bin_freq_scaled[end], pulse_freq[1], pulse_freq[end]]

        # Plot
        fig = PyPlot.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        img = ax.imshow(power, origin="lower", aspect="auto", cmap="gray_r",
                        extent=extent, vmin=0, vmax=quantile(vec(power), 0.95),
                        interpolation="none")
        
        ax.set_xlabel("Fluctuation frequency (P/P2)")
        ax.set_ylabel("Fluctuation frequency (P/P3)")
        ax.set_xlim(-100, 100)
        ax.set_ylim(0, 0.5)
        fig.suptitle("2DFS – $pulsar_name")

        PyPlot.colorbar(img, ax=ax, label="Power", shrink=0.9)

        savepath = joinpath(outdir, pulsar_name, "2dfs_p2freq_" * pulsar_name * ".png")
        PyPlot.savefig(savepath, dpi=300)

        println("✅ 2DFS (P2 freq) plot saved to: $savepath")

        if show_plot
            PyPlot.show()
        else
            PyPlot.close(fig)
        end

    catch e
        println("❌ Error while handling FITS file: $e")
        if f !== nothing
            close(f)
        end
    end
end





function plot_2dfs_kon(outdir::String, pulsar_name::String; show_plot::Bool=true)

    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("File does not exist: $filepath")
        return
    end
    

    println("File does exist and we can read 2DFS from : $filepath")

    f = nothing
    data = nothing
    try
        f = FITS(filepath, "r")
        hdu = f[4]
    header = read_header(hdu)  # <-- dodaj to
    data = read(hdu, "DATA")

        #check dimensions
        if ndims(data) != 2
            println("Data is not 2D")
            close(f)
            return
        end

        dat_scl = tryparse(Float64, get(header, "DAT_SCL", "1.0"))
        dat_offs = tryparse(Float64, get(header, "DAT_OFFS", "0.0"))
        period = tryparse(Float64, get(header, "PERIOD", "1.0"))  # seconds

    

        close(f)

        n_pulses, n_bins = size(data)
        println("Data shape: $n_pulses pulses × $n_bins bins")

        # Pulse longitude array
        pulse_long = LinRange(0, 360, n_bins+1)[1:end-1]

        # Frequency array (P/P3)
        freq = fftshift(fftfreq(n_pulses, 1))

        # Nowa zmienna - moc sygnału
        power = abs2.(data)

        # Wytnij fragment pulse longitude 160° - 200°
        mask = (pulse_long .>= 160) .& (pulse_long .<= 200)
        pulse_long_crop = pulse_long[mask]
        power_crop = power[:, mask]  # Uwaga: zachowaj wszystkie częstotliwości (wiersze), ale tylko wybrane kolumny!

        # Plotowanie
        figure(figsize=(8, 6))
        ax = gca()
        im = ax.imshow(
            power_crop, 
            aspect="auto", 
            cmap="Greys",
            extent=[pulse_long_crop[1], pulse_long_crop[end], freq[1], freq[end]],
            vmin=0, vmax=quantile(vec(power_crop), 0.95)
        )
        ax.set_xlabel("Pulse Longitude [deg]")
        ax.set_ylabel("P/P₃ [cpp]")
        ax.set_title("2DFS – $pulsar_name (160°–200°)")
        ax.set_ylim(0, 0.5)  # np. ograniczenie Y do sensownego zakresu
        colorbar(im, ax=ax, label="Power")

        savepath = joinpath(outdir, pulsar_name, "corrected_2dfs_" * pulsar_name * ".png")
        savefig(savepath)

        println("✅ 2DFS saved to: $savepath")

        if show_plot
            show()
        end

        

    catch e
        println("❌ Error while handling FITS file: $e")
    end
end


function inspect_fits223(filepath::String)
    if !isfile(filepath)
        println("❌ File does not exist: $filepath")
        return
    end

    println("🔎 Inspecting FITS file: $filepath")

    f = FITS(filepath)
    for (i, hdu) in enumerate(f)
        println("\n📂 HDU $i:")
        println("  Type: ", typeof(hdu))

        try
            if hdu isa FITSIO.ImageHDU
                img = read(hdu)
                println("  - Image dimensions: ", size(img))
                println("  - Image element type: ", eltype(img))

            elseif hdu isa FITSIO.TableHDU
                names = FITSIO.colnames(hdu)
                println("  - Table columns: ", names)

                for name in names
                    col_data = read(hdu, name)
                    println("    Column '$name': type ", eltype(col_data), ", size ", size(col_data))
                end

            else
                println("  - Unknown HDU type.")
            end
        catch e
            println("  ⚠️ Error reading HDU $i: $e")
        end
    end
    close(f)
end


function plot_2dfs_ostateczne(outdir, pulsar_name; title_text="2DFS Plot", cmap="gray", show_plot=false)
    # Generowanie ścieżki do pliku
    filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")

    if !isfile(filepath)
        println("File does not exist: $filepath")
        return
    end

    println("File exists and we can read 2DFS from: $filepath")

    f = nothing
    data = nothing
    try
        # Otwórz plik FITS
        f = FITS(filepath, "r")
        hdu = f[4]  # Czwarty HDU zawiera dane
        
        # Odczyt nagłówka i danych
        header = read_header(hdu)
        data = read(hdu, "DATA")

        # Sprawdzenie wymiarów danych
        if ndims(data) != 2
            println("Data is not 2D")
            close(f)
            return
        end

        # Odczytywanie wartości z nagłówka
        dat_scl = tryparse(Float64, get(header, "DAT_SCL", "1.0"))
        dat_offs = tryparse(Float64, get(header, "DAT_OFFS", "0.0"))
        period = tryparse(Float64, get(header, "PERIOD", "1.0"))  # Czas w sekundach
        
        # Skala danych (skalowanie, przesunięcie)
        data = (data .- dat_offs) .* dat_scl

        # Przygotowanie danych dla osi
        P2_vals = LinRange(-1.0, 1.0, size(data, 2))  # Przykładowe wartości dla 1/P2
        P3_vals = LinRange(0.0, 0.5, size(data, 1))  # Przykładowe wartości dla 1/P3

        # Tworzenie wykresu
        fig, ax = subplots(figsize=(8,6))

        # Rysowanie obrazu 2DFS
        cax = ax.imshow(data,
            extent=[minimum(P2_vals), maximum(P2_vals), minimum(P3_vals), maximum(P3_vals)],
            aspect="auto",
            origin="lower",
            cmap=cmap
        )

        # Dodanie paska kolorów
        colorbar(cax, ax=ax, label="Spectral Power")

        # Ustawienia osi
        ax.set_ylabel("Fluctuation Frequency (1/P₃) [cpp]")
        ax.set_xlabel("Fluctuation Frequency (1/P₂) [cpp]")

        # Linia pomocnicza w 1/P₂ = 0
        ax.axvline(x=0, color="white", linestyle="--", linewidth=1.5, label="1/P₂ = 0")

        # Symetria - odbicie w poziomie (jeśli to ma sens)
        vertical_profile = sum(data, dims=1)[:]  # Integracja w pionie (po 1/P₃)
        mirrored_profile = reverse(vertical_profile)  # Odbicie
        P2_centered = LinRange(minimum(P2_vals), maximum(P2_vals), length=length(vertical_profile))

        ax2 = ax.twinx()
        ax2.plot(P2_centered, mirrored_profile, color="cyan", linestyle="--", alpha=0.6, label="Mirrored Profile")
        ax2.set_yticks([])
        ax2.set_ylim(ax.get_ylim())

        # Tytuł wykresu
        ax.set_title(title_text)
        ax.legend(loc="upper right")
        tight_layout()
        
        # Zamykanie pliku FITS
        close(f)

        # Pokazanie wykresu, jeśli 'show_plot' jest ustawione na true
        if show_plot
            show()
        end

    catch e
        println("Error: $e")
        close(f)
    end
end




function wykres2dfs(folder_path, outdir, pulsar_name; show_plot=false)

    # Funkcja do odczytu danych z pliku FITS
    function read_fits_data(file_path, hdu_number, column_name)
        fits_file = FITSIO.FITS(file_path)               # Otwórz plik FITS
        hdu = fits_file[hdu_number]                      # Wybierz odpowiedni HDU
        data_column = read(hdu, column_name)             # Odczytaj dane z wybranej kolumny
        FITSIO.close(fits_file)                          # Zamknij plik FITS
        return data_column
    end
    
    # Funkcja do obliczania 2DFS
    function compute_2dfs(data)
        # Przekształcenie danych na Complex{Float64}
        data_complex = Complex{Float64}.(data)
    
        # Wykonanie FFT na pierwszym wymiarze
        data_fft = FFTW.fft(data_complex, 1)  # Wykonaj FFT wzdłuż pierwszego wymiaru (np. wiersze)
    
        # Wykonanie FFT na drugim wymiarze
        data_2dfs = FFTW.fft(data_fft, 2)  # Wykonaj FFT wzdłuż drugiego wymiaru (np. kolumny)
        
        return data_2dfs
    end
    
    
    
    # Funkcja do generowania wykresu 2DFS z użyciem plot_2dfs_ostateczne
    function plot_2dfs_ostateczne(outdir, pulsar_name; title_text="2DFS Plot", cmap="gray", show_plot=false)
    
        # Generowanie ścieżki do pliku
        filepath = joinpath(outdir, pulsar_name, "pulsar.debase.1.2dfs")
    
        if !isfile(filepath)
            println("File does not exist: $filepath")
            return
        end
    
        println("File exists and we can read 2DFS from: $filepath")
    
        f = nothing
        data = nothing
        try
            # Otwórz plik FITS
            f = FITS(filepath, "r")
            hdu = f[4]  # Czwarty HDU zawiera dane
            
            # Odczyt nagłówka i danych
            header = read_header(hdu)
            data = read(hdu, "DATA")
    
            # Sprawdzenie wymiarów danych
            if ndims(data) != 2
                println("Data is not 2D")
                close(f)
                return
            end
    
            # Odczytywanie wartości z nagłówka
            dat_scl = tryparse(Float64, get(header, "DAT_SCL", "1.0"))
            dat_offs = tryparse(Float64, get(header, "DAT_OFFS", "0.0"))
            period = tryparse(Float64, get(header, "PERIOD", "1.0"))  # Czas w sekundach
    
            # Skala danych (skalowanie, przesunięcie)
            data = (data .- dat_offs) .* dat_scl
    
            # Przygotowanie danych dla osi
            P2_vals = LinRange(-1.0, 1.0, size(data, 2))  # Przykładowe wartości dla 1/P2
            P3_vals = LinRange(0.0, 0.5, size(data, 1))  # Przykładowe wartości dla 1/P3
    
            # Tworzenie wykresu
            fig, ax = subplots(figsize=(8,6))
    
            # Rysowanie obrazu 2DFS
            cax = ax.imshow(data,
                extent=[minimum(P2_vals), maximum(P2_vals), minimum(P3_vals), maximum(P3_vals)],
                aspect="auto",
                origin="lower",
                cmap=cmap
            )
    
            # Dodanie paska kolorów
            colorbar(cax, ax=ax, label="Spectral Power")
    
            # Ustawienia osi
            ax.set_ylabel("Fluctuation Frequency (1/P₃) [cpp]")
            ax.set_xlabel("Fluctuation Frequency (1/P₂) [cpp]")
    
            # Linia pomocnicza w 1/P₂ = 0
            ax.axvline(x=0, color="white", linestyle="--", linewidth=1.5, label="1/P₂ = 0")
    
            # Symetria - odbicie w poziomie (jeśli to ma sens)
            vertical_profile = sum(data, dims=1)[:]  # Integracja w pionie (po 1/P₃)
            mirrored_profile = reverse(vertical_profile)  # Odbicie
            P2_centered = LinRange(minimum(P2_vals), maximum(P2_vals), length=length(vertical_profile))
    
            ax2 = ax.twinx()
            ax2.plot(P2_centered, mirrored_profile, color="cyan", linestyle="--", alpha=0.6, label="Mirrored Profile")
            ax2.set_yticks([])
            ax2.set_ylim(ax.get_ylim())
    
            # Tytuł wykresu
            ax.set_title(title_text)
            ax.legend(loc="upper right")
            tight_layout()
            
            # Zamykanie pliku FITS
            close(f)
    
            # Pokazanie wykresu, jeśli 'show_plot' jest ustawione na true
            if show_plot
                show()
            end
    
        catch e
            println("Error: $e")
            close(f)
        end
    end

    # --- Główna część programu: przetwarzanie pliku FITS ---
    
    # Bez potrzeby przetwarzania wielu plików, przechodzimy od razu do jednego
    file_path = joinpath(folder_path, "pulsar.debase.1.2dfs")
    if isfile(file_path)
        println("Przetwarzam plik: ", file_path)
        data = read_fits_data(file_path, 4, "DATA")  # Odczytaj dane z kolumny 'DATA' w HDU 4
        data_2dfs = compute_2dfs(data)                 # Oblicz 2DFS
        plot_2dfs_ostateczne(outdir, pulsar_name; title_text="2DFS for $pulsar_name", cmap="gray", show_plot=true)  # Wywołaj funkcję plot_2dfs_ostateczne
    else
        println("Plik pulsar.debase.1.2dfs nie istnieje w ścieżce: $file_path")
    end

    println("Przetwarzanie pliku zakończone.")
end



function check_and_plot_2dfs(file_path::String, col_p2::String, col_p3::String)
    FITS(file_path) do f
        hdu = f[4]  # 4. HDU
        
        if !(hdu isa FITSIO.TableHDU)
            println("Błąd: HDU 4 nie jest tabelą.")
            return
        end

        # Wczytaj dane kolumn
        p2_values = read(hdu, col_p2)
        p3_values = read(hdu, col_p3)

        # Wypisz rozmiary
        println("Rozmiar kolumny $col_p2: ", size(p2_values))
        println("Rozmiar kolumny $col_p3: ", size(p3_values))

        # Oblicz odwrotności
        inv_p2 = 1.0 .÷ p2_values
        inv_p3 = 1.0 .÷ p3_values

        # Wykres
        plot(inv_p2, inv_p3, label="1/$col_p2 vs 1/$col_p3", linewidth=2)
        xlabel("1/$col_p2")
        ylabel("1/$col_p3")
        legend()
        grid(true)
        title("Wykres odwrotności")
        show()

        # Wypisz pierwsze 10 wartości
        println("Pierwsze 10 wartości 1/$col_p2: ", inv_p2[1:10])
        println("Pierwsze 10 wartości 1/$col_p3: ", inv_p3[1:10])

        return inv_p2, inv_p3
    end
end
