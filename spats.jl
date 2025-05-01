module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


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
        #Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        #Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        #folded = Tools.p3fold(data, 20, 40)
        #Plot.single(folded, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319_p3fold", show_=true)
    end


    function load_ascii2(infile::String)
        raw_lines = readlines(infile)
        data = []
        for line in raw_lines
            parts = split(strip(line))
            if length(parts) >= 6
                row = tryparse.(Float64, parts[1:6])
                if all(!isnothing, row)
                    push!(data, row)
                end
            end
        end
        return data
    end
    


    function repuls(vpmout::String, num_files::Int)
        all_data = []
    
        # Regex do sprawdzania liczby: dopuszcza liczby dodatnie, ujemne, naukowe (1e-3 itp.)
        number_regex = r"^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$"
    
        for i in 1:num_files
            infile = "$(vpmout)$(i).txt"
            outfile = "$(vpmout)$(i)_zmiany.txt"
    
            raw_lines = readlines(infile)
    
            # Dokładne sprawdzanie i parsowanie tylko poprawnych linii
            data = []
            for line in raw_lines
                parts = split(strip(line))
                if length(parts) >= 6 && all(x -> occursin(number_regex, x), parts[1:6])
                    push!(data, parse.(Float64, parts[1:6]))
                end
            end
    
            # Jeśli brak danych, pomiń ten plik
            if isempty(data)
                println("Plik $(infile) nie zawiera prawidłowych danych – pominięty.")
                continue
            end
    
            # Oblicz moc sygnału z kolumn 5 i 6
            magnitudes = [sqrt(row[5]^2 + row[6]^2) for row in data]
            max_val = maximum(magnitudes)
            cap = 0.1 * max_val
    
            # Zmodyfikowane dane (nie usuwamy kolumn)
            modified_data = [
                (row[1], row[2], row[3], magnitude > cap ? row[4] : 0.0, row[5], row[6])
                for (row, magnitude) in zip(data, magnitudes)
            ]
    
            # Zapis do nowego pliku z formatowaniem liczb bez .0 jeśli niepotrzebne
            open(outfile, "w") do io
                for row in modified_data
                    formatted_row = [x == floor(x) ? Int(x) : x for x in row]
                    println(io, join(formatted_row, " "))
                end
            end
    
            # Dołącz do całości
            append!(all_data, [collect(row) for row in modified_data])
        end
    
        # Jeśli zebrano jakiekolwiek dane, zrób wykres
        if !isempty(all_data)
            # Debugowanie: sprawdzenie danych przed generowaniem wykresu
            println("Dane zebrane do wykresu:")
            println(all_data[1:min(5, end)])  # Wypisz pierwsze 5 wierszy
    
            # Zbieramy wszystkie dane z plików wynikowych
            data = []
            for i in 1:num_files
                zmianyfile = "$(vpmout)$(i)_zmiany.txt"
                if isfile(zmianyfile)
                    append!(data, load_ascii2(zmianyfile))
                end
            end
    
            # Debugowanie: sprawdzenie danych przed generowaniem wykresu
            println("Dane przekazane do wykresu:")
            println(data[1:min(5, end)])
    
            # Użycie Plot.single do wygenerowania wykresu
            Plot.single(data, vpmout; darkness=0.5, number=nothing, start=1, name_mod="J1319", show_=true)
        else
            println("Nie znaleziono żadnych prawidłowych danych do wykresu.")
        end
    end
    
   
    
    
    
    
    
    



    function process_psrdata(indir, outdir)
        p = Data.process_psrdata(indir, outdir)
        folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        Plot.p3fold(folded, outdir; start=3, bin_st=p["bin_st"]-20, bin_end=p["bin_end"]+20, name_mod="test", show_=true, repeat_num=4)
        
    end


    function main()
        # output directory for VPM
        vpmout = "/home/psr/output/"
        num_files = 4  # lub inna liczba plików, które chcesz przetworzyć
        repuls(vpmout, num_files)
        #test(vpmout)
        #test2(vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)

    end

end # module

SpaTs.main()

println("Bye")