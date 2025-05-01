module SpaTs
    using ArgParse
    using Glob
    using JSON
    using Printf

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



    function repuls(vpmout::String, num_files::Int)
        all_data = []

        # Regex do sprawdzania liczby: dopuszcza liczby dodatnie, ujemne, naukowe (1e-3 itp.)
        number_regex = r"^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$"

        for i in 1:num_files
            infile = "$(vpmout)$(i).txt"               # Wczytujemy plik np. 1.txt, 2.txt, itd.
            outfile = "$(vpmout)$(i)_zmiany.txt"       # Tworzymy nowy plik zmodyfikowany

            raw_lines = readlines(infile)

            data = []
            for line in raw_lines
                parts = split(strip(line))
                if length(parts) >= 7 && all(x -> occursin(number_regex, x), parts[1:7])
                    push!(data, parse.(Float64, parts[1:7]))  # Parsujemy kolumny 1–7
                end
            end

            if isempty(data)
                println("Plik $(infile) nie zawiera prawidłowych danych – pominięty.")
                continue
            end

            # Oblicz magnitudę sygnału z kolumn 5 i 6
            magnitudes = [sqrt(row[5]^2 + row[6]^2) for row in data]
            max_val = maximum(magnitudes)
            cap = 0 * max_val  # 10% z maksymalnej magnitudy

            # Modyfikujemy TYLKO kolumnę 4 — inne pozostają bez zmian
            modified_data = [
                begin
                    new_row = copy(row)
                    new_row[4] = magnitude > cap ? row[4] : 0.0
                    new_row
                end
                for (row, magnitude) in zip(data, magnitudes)
            ]

            # Zapis do nowego pliku
            open(outfile, "w") do io
                for row in modified_data
                    # Formatujemy zapis, aby nie zmieniać formatu liczb (zachowujemy oryginalny format)
                    formatted_row = [x == floor(x) ? string(Int(x)) : @sprintf("%.6f", x) for x in row]
                    println(io, join(formatted_row, " "))
                end
            end

            # Zbieramy dane do wspólnej tabeli
            append!(all_data, [collect(row) for row in modified_data])
        end

        # Ładowanie danych z plików zmodyfikowanych
        data1 = Data.load_ascii(vpmout * "1_zmiany.txt")
        data2 = Data.load_ascii(vpmout * "2_zmiany.txt")
        data3 = Data.load_ascii(vpmout * "3_zmiany.txt")
        data4 = Data.load_ascii(vpmout * "4_zmiany.txt")

        # Łączenie danych
        combined_data = vcat(data1, data2, data3, data4)

        # Wykres
        if !isempty(all_data)
            Plot.single(combined_data, vpmout; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
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