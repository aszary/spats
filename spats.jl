module SpaTs
    using ArgParse
    using Glob
    using JSON
    using FITSIO
    using PyPlot
    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")


    function J0820(args)
        # J0820
        #Data.convert_psrfit_ascii("input/J0820.high", "input/high.txt") # does not work without TEMPO?
        #Data.convert_psrfit_ascii("input/J0820.low", "input/low.txt")

        data = Data.load_ascii("input/low.txt")
        data2 = Data.load_ascii("input/high.txt")

        Plot.single0(data, args["outdir"]; bin_st=400, bin_end=600)

        Plot.single(data, args["outdir"]; bin_st=400, bin_end=600, name_mod="low")
        Plot.single(data2, args["outdir"]; bin_st=400, bin_end=600, name_mod="high")

        Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1, bin_st=400, bin_end=600, name_mod="low")
        Plot.lrfs(data2, args["outdir"]; darkness=0.1, start=1, bin_st=400, bin_end=600, name_mod="high")

        p3data = Data.load_ascii("input/low.debase.p3fold")
        Plot.p3fold(data, args["outdir"]; start=1, bin_st=400, bin_end=600, name_mod="low", darkness=1.0)

        # OLD
        #Plot.offset_subtract(data, data2, "output"; bin_st=400, bin_end=600, name_mod="sub", darkness=1.0, repeat_num=1)  # it may look nice? Nope! but why? try subtract p3folds

        #Plot.average(data, args["outdir"]; bin_st=40, bin_end=90, number=nothing, name_mod="low")
        #Plot.average(data3, args["outdir"]; bin_st=40, bin_end=90, number=nothing, name_mod="high4")
        #Plot.average2(data, data3, args["outdir"]; bin_st=40, bin_end=90, number=nothing, name_mod="lowhigh4")
        #Plot.offset_subtract(data, data3, "output"; bin_st=40, bin_end=90, name_mod="B0320_p3fold", darkness=1.0, repeat_num=4)  # it may look nice? Nope! but why? try subtract p3folds


        #data = Data.load_ascii("input/low.p3fold")
        #data2 = Data.load_ascii("input/mid.p3fold")
        #data3 = Data.load_ascii("input/high4.p3fold")
        #Plot.offset(data, data3, "output/"; bin_st=50, bin_end=80, name_mod="low_high", darkness=1.0)
        #Plot.offset_points(data, data3, "output/"; bin_st=50, bin_end=80, name_mod="db_low_high")
        #Plot.offset_points3(data, data2, data3, "output/"; bin_st=50, bin_end=80, name_mod="low_mid_high4")
        #Plot.crosscorplot(data, data2, "/home/szary/work/B0320/"; bin_st=50, bin_end=80, name_mod="low_high4")
        #Plot.offset(data, data2, args["outdir"]; bin_st=50, bin_end=80, name_mod="low_high4", darkness=1.0)

        #Plot.single0(data, args["outdir"]; number=250, bin_st=50, bin_end=100)
        #Plot.single0(data, args["outdir"])
        #Plot.single(data, args["outdir"]; darkness=0.5, number=nothing, name_mod="all")
        #Plot.single(data, args["outdir"]; start=1, number=17, bin_st=50, bin_end=80, name_mod="low_p3fold", darkness=1.0)
        #Plot.p3fold(data, args["outdir"]; start=1, bin_st=50, bin_end=80, name_mod="high4", darkness=1.0)
        #Plot.single(data, args["outdir"]; start=1, number=512, bin_st=30, bin_end=95, name_mod="p3fold")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1, number=512, bin_st=50, bin_end=80, name_mod="low")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1000, number=512, bin_st=40, bin_end=85, name_mod="2")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1250, number=512, bin_st=40, bin_end=85, name_mod="3")
        #Plot.p3_evolution(data, args["outdir"]; darkness=0.1, bin_st=40, bin_end=85, name_mod="1", number=128)
        #Plot.p3_evolution(data, args["outdir"]; darkness=0.1, bin_st=40, bin_end=85, name_mod="2", number=256)
    end



    function mkieth()

        data = []
        push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/c_1.txt"))
        push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/c_101.txt"))
        push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/c_201.txt"))
        #push!(data, Data.load_ascii("/home/szary/work/J1651-4246/mkieth/b_301.txt"))
        println(size(data))

        d = Array{Float64}(undef, 100 * (size(data)[1]-1) + size(data[end])[1], size(data[end])[2])

        for (i,da) in enumerate(data)
            pulses, bins = size(da)
            for j in 1:pulses
                for k in i:bins
                    d[100*(i-1)+j, k] = da[j, k]
                end
            end
        end
        println(size(d))
        Plot.single(d, args["outdir"]; darkness=0.5, number=nothing, name_mod="J1651-4246")
        Plot.lrfs(d, args["outdir"]; darkness=0.1, start=1, name_mod="J1651-4246", change_fftphase=false)

    end


    function J1651()
        data = Data.load_ascii("/home/szary/work/MeerTime/J1651/1700.txt")
        Plot.single(data, args["outdir"]; start=1, number=nothing, darkness=0.5, name_mod="J1651")
        Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1,  bin_st=300, bin_end=700, name_mod="J1651", change_fftphase=false)
    end


    function J1705()
        data = Data.load_ascii("/fred/oz005/users/aszary/search/J1705-1906/2019-04-03T02:19:47/1444.5/grand.txt")
        Plot.single(data, args["outdir"]; darkness=0.5, number=nothing, name_mod="J1705-1906")
        Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1,  bin_st=400, bin_end=500, name_mod="J1705-1906", change_fftphase=false)
    end


    function B0320()
        # B0320
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.high4", "input/high4.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.low", "input/low.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.pazi", "input/full_full.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/add.ar.mid", "input/mid.txt")
        #Data.convert_psrfit_ascii("/data/szary/B0320/core/rawvoltages/SAP0/BEAM0/B0320+39_L570031_SAP0_BEAM0.paz.fscr.pdmp.AR", "input/core/")

        #data = Data.load_ascii("input/low.txt")
        #data3 = Data.load_ascii("input/high4.txt")

        data = Data.load_ascii("input/low.p3fold")
        data3 = Data.load_ascii("input/high4.p3fold")

        #Plot.average(data, args["outdir"]; bin_st=40, bin_end=90, number=nothing, name_mod="low")
        #Plot.average(data3, args["outdir"]; bin_st=40, bin_end=90, number=nothing, name_mod="high4")
        #Plot.average2(data, data3, args["outdir"]; bin_st=40, bin_end=90, number=nothing, name_mod="lowhigh4")
        Plot.offset_subtract(data, data3, "output"; bin_st=40, bin_end=90, name_mod="B0320_p3fold", darkness=1.0, repeat_num=4)  # it may look nice? Nope! but why? try subtract p3folds


        #data = Data.load_ascii("input/low.p3fold")
        #data2 = Data.load_ascii("input/mid.p3fold")
        #data3 = Data.load_ascii("input/high4.p3fold")
        #Plot.offset(data, data3, "output/"; bin_st=50, bin_end=80, name_mod="low_high", darkness=1.0)
        #Plot.offset_points(data, data3, "output/"; bin_st=50, bin_end=80, name_mod="db_low_high")
        #Plot.offset_points3(data, data2, data3, "output/"; bin_st=50, bin_end=80, name_mod="low_mid_high4")
        #Plot.crosscorplot(data, data2, "/home/szary/work/B0320/"; bin_st=50, bin_end=80, name_mod="low_high4")
        #Plot.offset(data, data2, args["outdir"]; bin_st=50, bin_end=80, name_mod="low_high4", darkness=1.0)

        #Plot.single0(data, args["outdir"]; number=250, bin_st=50, bin_end=100)
        #Plot.single0(data, args["outdir"])
        #Plot.single(data, args["outdir"]; darkness=0.5, number=nothing, name_mod="all")
        #Plot.single(data, args["outdir"]; start=1, number=17, bin_st=50, bin_end=80, name_mod="low_p3fold", darkness=1.0)
        #Plot.p3fold(data, args["outdir"]; start=1, bin_st=50, bin_end=80, name_mod="high4", darkness=1.0)
        #Plot.single(data, args["outdir"]; start=1, number=512, bin_st=30, bin_end=95, name_mod="p3fold")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1, number=512, bin_st=50, bin_end=80, name_mod="low")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1000, number=512, bin_st=40, bin_end=85, name_mod="2")
        #Plot.lrfs(data, args["outdir"]; darkness=0.1, start=1250, number=512, bin_st=40, bin_end=85, name_mod="3")
        #Plot.p3_evolution(data, args["outdir"]; darkness=0.1, bin_st=40, bin_end=85, name_mod="1", number=128)
        #Plot.p3_evolution(data, args["outdir"]; darkness=0.1, bin_st=40, bin_end=85, name_mod="2", number=256)
    end


    function J1750_remote()

        # ozStar data
        ozdir1 = "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-09-29-12:41:18/1284/single/"
        ozdir2 = "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/"
        ozdir3 = "/fred/oz005/users/aszary/search_processed/J1750-3503/2020-02-24-03:46:43/1284/single/"
        ozdir4 = "/fred/oz005/users/aszary/search_processed/J1750-3503/2020-03-29-04:11:50/1284/single/"
        ozdir5 = "/fred/oz005/users/aszary/search_processed/J1750-3503/2020-05-07-23:14:22/1284/single/"
        ozdir6 = "/fred/oz005/users/aszary/search_processed/J1750-3503/2020-05-30-22:04:58/1284/single/"


        # first session
        #Data.convert_psrfit_ascii("/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-09-29-12:41:18/1284/single/2019-09-29-12:41:18_00000-00294.spCF", "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-09-29-12:41:18/1284/single/2019-09-29-12:41:18_00000-00294.txt")
        #data = Data.load_ascii("/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-09-29-12:41:18/1284/single/2019-09-29-12:41:18_00000-00294.txt")
        #Plot.single(data, "output"; start=1, number=nothing, name_mod="J1750_1")
        #Plot.single(data, "output"; start=1, number=nothing, bin_st=450, bin_end=700, name_mod="J1750_1_zoom")
        #Plot.lrfs(data, "output"; darkness=0.3, start=1, name_mod="J1750_1", bin_st=450, bin_end=700, change_fftphase=false)
        #Plot.p3_evolution(data, "output"; darkness=0.3, bin_st=450, bin_end=700, name_mod="J1750_1_128", number=128)


        # second session
        #Data.convert_psrfit_ascii("/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00000-00255.spCF", "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00000-00255.txt")
        #Data.convert_psrfit_ascii("/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00256-00511.spCF", "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00256-00511.txt")
        #Data.convert_psrfit_ascii("/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00512-00767.spCF", "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00512-00767.txt")
        #Data.convert_psrfit_ascii("/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00768-01030.spCF", "/fred/oz005/users/aszary/search_processed/J1750-3503_old/2019-12-14-14:22:12/1284/single/2019-12-14-14:22:12_00768-01030.txt")

        #=
        data1 = Data.load_ascii("$(ozdir2)2019-12-14-14:22:12_00000-00255.txt")
        data2 = Data.load_ascii("$(ozdir2)2019-12-14-14:22:12_00256-00511.txt")
        data3 = Data.load_ascii("$(ozdir2)2019-12-14-14:22:12_00512-00767.txt")
        data = vcat(data1, data2, data3, data4)
        data4 = Data.load_ascii("$(ozdir2)2019-12-14-14:22:12_00768-01030.txt")
        #Data.save_ascii(data, "output/2019-12-14-14:22:12_00000-01030.txt")
        #println(size(data))
        #Plot.single(data, "output"; start=1, number=nothing, name_mod="J1750_2")
        #Plot.lrfs(data, "output"; darkness=0.3, start=1, name_mod="J1750_2", bin_st=450, bin_end=700, change_fftphase=false)
        =#

        # third session
        #Data.convert_psrfit_ascii("$(ozdir3)2020-02-24-03:46:43_00000-00255.spCF", "$(ozdir3)2020-02-24-03:46:43_00000-00255.txt")
        #Data.convert_psrfit_ascii("$(ozdir3)2020-02-24-03:46:43_00256-00445.spCF", "$(ozdir3)2020-02-24-03:46:43_00256-00445.txt")

        # fourth session
        #Data.convert_psrfit_ascii("$(ozdir4)2020-03-29-04:11:50_00000-00255.spCF", "$(ozdir4)2020-03-29-04:11:50_00000-00255.txt")
        #Data.convert_psrfit_ascii("$(ozdir4)2020-03-29-04:11:50_00256-00449.spCF", "$(ozdir4)2020-03-29-04:11:50_00256-00449.txt")

        # fith session
        #Data.convert_psrfit_ascii("$(ozdir5)2020-05-07-23:14:22_00000-00255.spCF", "$(ozdir5)2020-05-07-23:14:22_00000-00255.txt")
        #Data.convert_psrfit_ascii("$(ozdir5)2020-05-07-23:14:22_00256-00446.spCF", "$(ozdir5)2020-05-07-23:14:22_00256-00446.txt")

        # sixth session
        Data.convert_psrfit_ascii("$(ozdir6)2020-05-30-22:04:58_00000-00255.spCF", "$(ozdir6)2020-05-30-22:04:58_00000-00255.txt")
        Data.convert_psrfit_ascii("$(ozdir6)2020-05-30-22:04:58_00256-00446.spCF", "$(ozdir6)2020-05-30-22:04:58_00256-00446.txt")


        # seventh session in the future
        #Data.convert_psrfit_ascii("$(ozdir7)", "$(ozdir7).txt")
        #Data.convert_psrfit_ascii("$(ozdir7)", "$(ozdir7).txt")

    end


    function J1750_local()
        # local data
        data_dir = "/home/szary/work/MeerTime/J1750/data/"
        outdir = "/home/szary/work/MeerTime/J1750/"
        outdir2 = "/home/szary/work/MeerTime/J1750/tracks2/"

        # first session
        #data = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #Plot.single(data, outdir; start=725, number=125, bin_st=450, bin_end=700, name_mod="J1750_2_part")
        #Plot.lrfs(data, outdir; darkness=0.3, start=725, number=125, name_mod="J1750_2_part", bin_st=450, bin_end=700, change_fftphase=false)
        #Plot.p3_evolution(data, outdir; darkness=0.5, bin_st=450, bin_end=700, name_mod="J1750_1_128", number=128)

        #data = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294_refine.p3fold")
        #Plot.p3fold(data, outdir; start=1, number=46, bin_st=nothing, bin_end=nothing, name_mod="J1750_1", darkness=1.0, cmap="viridis")
        #Plot.p3fold(data, outdir; start=1, number=46, bin_st=450, bin_end=700, name_mod="J1750_1_refine", darkness=1.0, cmap="viridis")

        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294_low.p3fold")
        #data2 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294_high.p3fold")
        #Plot.offset_subtract(data1, data2, outdir; bin_st=450, bin_end=700, name_mod="J1750_1_p3fold", darkness=1.0, repeat_num=3)  # it may look nice? Nope! but why? try subtract p3folds


        # second session
        #data = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #Plot.single(data, outdir; start=725, number=125, bin_st=450, bin_end=700, name_mod="J1750_2_part")
        #Plot.single(data, outdir; start=1, number=300, bin_st=450, bin_end=700, name_mod="J1750_2_part1")
        #Plot.lrfs(data, outdir; darkness=0.3, start=725, number=125, name_mod="J1750_2_part", bin_st=450, bin_end=700, change_fftphase=false)
        #Plot.p3_evolution(data, outdir; darkness=0.5, bin_st=450, bin_end=700, name_mod="J1750_2_64", number=64)

        #data = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030_refine.p3fold")
        #Plot.p3fold(data, outdir; start=1, number=nothing, bin_st=450, bin_end=700, name_mod="J1750_2_refine", darkness=1.0, cmap="viridis")

        #data = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.p3fold")
        #Plot.p3fold(data, outdir; start=1, number=nothing, bin_st=450, bin_end=700, name_mod="J1750_2", darkness=1.0, cmap="viridis")


        # third session
        #data1 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00255.txt")
        #data2 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00256-00445.txt")
        #data = vcat(data1, data2)
        #Data.save_ascii(data, "$(data_dir)2020-02-24-03:46:43_00000-00445.txt")

        #data = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        #Data.zap!(data; ranges=[241, 335, 334])
        #Plot.single(data, outdir; start=1, number=nothing, bin_st=415, bin_end=620, name_mod="J1750_3", darkness=0.3)
        #Plot.single(data, outdir; start=1, number=nothing, bin_st=nothing, bin_end=nothing, name_mod="J1750_3", darkness=0.7)
        #Plot.lrfs(data, outdir; darkness=0.3, start=1, number=256, name_mod="J1750_4_part", bin_st=400, bin_end=700, change_fftphase=false)
        #Plot.p3_evolution(data, outdir; darkness=0.5, bin_st=400, bin_end=700, name_mod="J1750_3_128", number=128)
        #Plot.p3_evolution(data, outdir; darkness=0.5, bin_st=400, bin_end=700, name_mod="J1750_3_128", number=128)


        # both
        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #Plot.average2(data1, data2, outdir, bin_st=450, bin_end=700, number=nothing, name_mod="J1750")


        # offsets
        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294_low.txt")
        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #data2 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294_high.txt")
        #data2 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #Plot.average2(data1, data2, outdir, bin_st=450, bin_end=700, number=nothing, name_mod="J1750_low_high")
        #Plot.offset(data1, data2, outdir; bin_st=450, bin_end=700, name_mod="J1750_1", darkness=0.7, repeat_num=1)
        #Plot.offset_subtract(data1, data2, outdir; bin_st=450, bin_end=700, name_mod="J1750_1", darkness=0.7, repeat_num=1)
        #Plot.offset_points(data1, data2, outdir; bin_st=450, bin_end=700, name_mod="J1750_1", repeat_num=1)


        # P2 estimate
        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        #data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        #data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        #Plot.single(data1, outdir; start=1, number=nothing, bin_st=nothing, bin_end=nothing, name_mod="J1750_test", darkness=0.3)
        #Tools.p2_estimate(data1; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=8) # 18
        #Tools.p2_estimate(data2; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=4, template_num=37) # 18
        #Tools.p2_estimate(data3; on_st=350, on_end=650, off_st=50, off_end=350, thresh=3.5, win=15) #  no! # to noisy!
        #Tools.p2_estimate(data4; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=12) # 17.6
        #Tools.p2_estimate(data5; on_st=450, on_end=700, off_st=100, off_end=350, thresh=3.5, win=16, template_num=nothing) # 18


        # fourth session
        #data1 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00255.txt")
        #data2 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00256-00449.txt")
        #data = vcat(data1, data2)
        #Data.save_ascii(data, "$(data_dir)2020-03-29-04:11:50_00000-00449.txt")

        #data = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        #Plot.single(data, outdir; start=1, number=nothing, bin_st=400, bin_end=700, name_mod="J1750_4", darkness=0.3)
        #Plot.single(data, outdir; start=1, number=nothing, bin_st=nothing, bin_end=nothing, name_mod="J1750_4", darkness=0.3)
        #Plot.p3_evolution(data, outdir; darkness=0.5, bin_st=400, bin_end=700, name_mod="J1750_4", number=128)
        #Plot.lrfs(data, outdir; darkness=0.3, start=1, number=256, name_mod="J1750_4_part", bin_st=400, bin_end=700, change_fftphase=false)


        # fifth session
        #data1 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00255.txt")
        #data2 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00256-00446.txt")
        #data = vcat(data1, data2)
        #Data.save_ascii(data, "$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        #data = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        #Plot.single(data, outdir; start=1, number=nothing, bin_st=400, bin_end=700, name_mod="J1750_5", darkness=0.3)
        #Plot.single(data, outdir; start=1, number=nothing, bin_st=nothing, bin_end=nothing, name_mod="J1750_5", darkness=0.3)
        #Plot.p3_evolution(data, outdir; darkness=0.5, bin_st=400, bin_end=700, name_mod="J1750_5", number=128)
        #Plot.lrfs(data, outdir; darkness=0.3, start=1, number=256, name_mod="J1750_5", bin_st=400, bin_end=700, change_fftphase=false)


        # average profiles - need to align, check X-ray project
        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        #data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        #data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        #Plot.averageX([data1, data2, data3, data4, data5], outdir, bin_st=400, bin_end=700, number=nothing, name_mod="J1750")

        # track subpulses
        #data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        #data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        #data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        #data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        #Data.zap!(data3; ranges=[241, 335, 334])
        #Plot.single(data1, outdir; start=1, number=10, bin_st=400, bin_end=700, name_mod="short", darkness=0.6)

        #peaks = Tools.track_subpulses(data1, 18.0, thresh=2.1, thresh2=0.7, on_st=500, on_end=650)
        #Plot.subpulses(data1, outdir, peaks; start=1, number=294, bin_st=500, bin_end=650, name_mod="J1750_1_new", darkness=0.6)
        #peaks = Tools.track_subpulses(data2, 18.0, thresh=2.2, thresh2=0.8, on_st=500, on_end=650)
        #Plot.subpulses(data2, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_2", darkness=0.6)
        #peaks = Tools.track_subpulses(data3, 18.0, thresh=2.2, thresh2=0.8, on_st=425, on_end=600)
        #Plot.subpulses(data3, outdir, peaks; start=1, number=256, bin_st=425, bin_end=600, name_mod="J1750_3", darkness=0.6)
        #peaks = Tools.track_subpulses(data4, 18.0, thresh=2.1, thresh2=0.8, on_st=500, on_end=650)
        #Plot.subpulses(data4, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_4", darkness=0.6)

        # HERE
        #p2, template = Tools.p2_estimate(data2; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=4, template_num=37) # 18
        #peaks = Tools.track_subpulses(data2, p2, thresh=2.1, thresh2=0.7, on_st=500, on_end=650)
        #peaks = Tools.track_subpulses_template(data2, template, thresh=2.1, thresh2=0.7, on_st=500, on_end=650)
        #Plot.subpulses(data2, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_2_temp", darkness=0.6)

        #peaks = Tools.track_subpulses(data5, 18, thresh=2.1, thresh2=0.7, on_st=450, on_end=650)
        #Plot.subpulses(data5, outdir, peaks; start=1, number=256, bin_st=450, bin_end=650, name_mod="J1750_5", darkness=0.6)
        #Plot.subpulses(data5, outdir, peaks; start=1, number=256, bin_st=450, bin_end=650, name_mod="J1750_5", darkness=0.6)

        # convolution + template test here
        #x = collect(1:100)
        #data = Tools.gauss(x, [1, 30, 7, 0])
        #template = Tools.gauss(x, [0.5, 60, 3, 0])
        #peaks = Tools.track_subpulses_template(data, template, thresh=2.1, thresh2=0.7, on_st=500, on_end=650)
        # m1 + m2 - (m3+1) = 0

        # template generation
        #data = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #p3data = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.p3fold")
        #template, pulses = Tools.template(p3data, 18; on_st=450, on_end=700, off_st=100, off_end=350, dbins=[-0.5727272727272728]) # dbins=LinRange(-0.55, -0.6, 100)
        #println(length(template))
        #peaks = Tools.track_subpulses_template(data, template, thresh=2.1, thresh2=0.7, on_st=500, on_end=650)
        #Plot.single(pulses, outdir; start=1, number=nothing, name_mod="template", darkness=0.6, show_=true)

        # subpulse tracking (Second Session)
        #data = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        #peaks = Tools.track_subpulses(data, 18, thresh=0.4, thresh2=0.5, on_st=450, on_end=700) # better
        #peaks = Tools.track_subpulses_smoothing(data; lambda=1000.0) # good
        #tracks = Tools.group_tracks_obsolete(peaks, 18.0)
        #Plot.subpulses(data, outdir, peaks; start=1, number=1031, bin_st=450, bin_end=700, name_mod="subpulses", darkness=0.6)
        #Plot.tracks(data, outdir, tracks, peaks; start=1, number=1031, bin_st=450, bin_end=700, name_mod="tracks", darkness=0.6)
        ##Plot.group_tracks(data, outdir2, peaks; start=1, number=1031, bin_st=450, bin_end=700, name_mod="2", darkness=0.6)
        #Plot.show_tracks(data, outdir2, peaks; start=1, number=1031, bin_st=450, bin_end=700, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis(outdir2; win=100, start=1, number=1031, bin_st=450, bin_end=700, name_mod="2", darkness=0.6)
        ##Plot.tracks_analysis2(outdir2; win=30, start=1, number=1031, bin_st=450, bin_end=700, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis3(outdir2; lambda=10000.0, start=1, number=1031, bin_st=450, bin_end=700, name_mod="2", darkness=0.6)

        # subpulse tracking (First Session)
        data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        peaks = Tools.track_subpulses(data1, 18, thresh=0.4, thresh2=0.5, on_st=450, on_end=700) # better
        ##Plot.group_tracks(data1, "$(outdir)tracks1", peaks; start=1, number=295, bin_st=450, bin_end=700, name_mod="1", darkness=0.6)
        #Plot.show_tracks(data1, "$(outdir)tracks1", peaks; start=1, number=295, bin_st=450, bin_end=700, name_mod="1", darkness=0.6)
        #Plot.tracks_analysis("$(outdir)tracks1"; win=50, start=1, bin_st=450, bin_end=700, name_mod="1", darkness=0.6)
        #Plot.tracks_analysis2("$(outdir)tracks1"; win=50, start=1, bin_st=450, bin_end=700, name_mod="1", darkness=0.6)
        Plot.tracks_analysis3("$(outdir)tracks1"; lambda=10000.0, start=1, number=1031, bin_st=450, bin_end=700, name_mod="1", darkness=0.6)


        # subpulse tracking (Third Session)
        #data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        #Data.zap!(data3; ranges=[241, 335, 334])
        #peaks = Tools.track_subpulses(data3, 18, thresh=0.4, thresh2=0.5, on_st=350, on_end=650) # better
        ##Plot.group_tracks(data3, "$(outdir)tracks3", peaks; start=1, number=446, bin_st=400, bin_end=650, name_mod="3", darkness=0.6)
        #Plot.tracks_analysis("$(outdir)tracks3"; win=50, start=1, bin_st=450, bin_end=700, name_mod="3", darkness=0.6)
        #Plot.tracks_analysis2("$(outdir)tracks3"; win=50, start=1, bin_st=450, bin_end=700, name_mod="3", darkness=0.6)
        #Plot.tracks_analysis3("$(outdir)tracks3"; lambda=1000.0, start=1, bin_st=450, bin_end=700, name_mod="3", darkness=0.6)


        # TODO subpulse tracking (Fourth Session)
        #data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        #peaks = Tools.track_subpulses(data4, 18, thresh=0.4, thresh2=0.5, on_st=350, on_end=650) # better
        #Plot.group_tracks(data4, "$(outdir)tracks4", peaks; start=1, bin_st=400, bin_end=650, name_mod="4", darkness=0.6)
        #Plot.tracks_analysis("$(outdir)tracks4"; win=50, start=1, bin_st=450, bin_end=700, name_mod="4", darkness=0.6)
        #Plot.tracks_analysis2("$(outdir)tracks4"; win=50, start=1, bin_st=450, bin_end=700, name_mod="4", darkness=0.6)
        #Plot.tracks_analysis3("$(outdir)tracks4"; lambda=1000.0, start=1, bin_st=450, bin_end=700, name_mod="4", darkness=0.6)

        # sixth session data connection
        #data1 = Data.load_ascii("$(data_dir)2020-05-30-22:04:58_00000-00255.txt")
        #data2 = Data.load_ascii("$(data_dir)2020-05-30-22:04:58_00256-00446.txt")
        #data = vcat(data1, data2)
        #Data.save_ascii(data, "$(data_dir)2020-05-30-22:04:58_00000-00446.txt")
        #data6 = Data.load_ascii("$(data_dir)2020-05-30-22:04:58_00000-00446.txt")
        #Plot.single(data6, outdir; start=1, number=nothing, name_mod="6_", darkness=0.6, show_=true)
    end

    function J1750_paper()

        data_dir = "/home/szary/work/MeerTime/J1750/new_data/"
        old_dir = "/home/szary/work/MeerTime/J1750/data/"
        data_remote = "/home/aszary/J1750/new_data/"
        data_patrick = "/home/szary/work/MeerTime/J1750/data_patrick/"

        #outdir = "/home/szary/work/MeerTime/J1750/"
        outdir = "/home/szary/work/MeerTime/J1750/old_plots/"

        #Data.convert_psrfit_ascii("$(data_patrick)/20190929_124118.debase.hp", "$(data_patrick)/2019-09-29-12:41:18_00000-00294.txt")
        # convert to txt
        #=
        Data.convert_psrfit_ascii("$(data_remote)/2019-09-29-12:41:18_00000-00294.spCF", "$(data_remote)/2019-09-29-12:41:18_00000-00294.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2019-12-14-14:22:12_00000-00255.spCF", "$(data_remote)/2019-12-14-14:22:12_00000-00255.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2019-12-14-14:22:12_00256-00511.spCF", "$(data_remote)/2019-12-14-14:22:12_00256-00511.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2019-12-14-14:22:12_00512-00767.spCF", "$(data_remote)/2019-12-14-14:22:12_00512-00767.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2019-12-14-14:22:12_00768-01030.spCF", "$(data_remote)/2019-12-14-14:22:12_00768-01030.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-02-24-03:46:43_00000-00255.spCF", "$(data_remote)/2020-02-24-03:46:43_00000-00255.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-02-24-03:46:43_00256-00445.spCF", "$(data_remote)/2020-02-24-03:46:43_00256-00445.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-03-29-04:11:50_00000-00255.spCF", "$(data_remote)/2020-03-29-04:11:50_00000-00255.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-03-29-04:11:50_00256-00449.spCF", "$(data_remote)/2020-03-29-04:11:50_00256-00449.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-05-07-23:14:22_00000-00255.spCF", "$(data_remote)/2020-05-07-23:14:22_00000-00255.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-05-07-23:14:22_00256-00446.spCF", "$(data_remote)/2020-05-07-23:14:22_00256-00446.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-05-30-22:04:58_00000-00255.spCF", "$(data_remote)/2020-05-30-22:04:58_00000-00255.txt")
        Data.convert_psrfit_ascii("$(data_remote)/2020-05-30-22:04:58_00256-00446.spCF", "$(data_remote)/2020-05-30-22:04:58_00256-00446.txt")
        =#

        # connect data
        #=
        data1 = Data.load_ascii("$(data_remote)2019-12-14-14:22:12_00000-00255.txt")
        data2 = Data.load_ascii("$(data_remote)2019-12-14-14:22:12_00256-00511.txt")
        data3 = Data.load_ascii("$(data_remote)2019-12-14-14:22:12_00512-00767.txt")
        data4 = Data.load_ascii("$(data_remote)2019-12-14-14:22:12_00768-01030.txt")
        data = vcat(data1, data2, data3, data4)
        Data.save_ascii(data, "$(data_remote)2019-12-14-14:22:12_00000-01030.txt")
        data1 = Data.load_ascii("$(data_remote)2020-02-24-03:46:43_00000-00255.txt")
        data2 = Data.load_ascii("$(data_remote)2020-02-24-03:46:43_00256-00445.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_remote)2020-02-24-03:46:43_00000-00445.txt")
        data1 = Data.load_ascii("$(data_remote)2020-03-29-04:11:50_00000-00255.txt")
        data2 = Data.load_ascii("$(data_remote)2020-03-29-04:11:50_00256-00449.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_remote)2020-03-29-04:11:50_00000-00449.txt")
        data1 = Data.load_ascii("$(data_remote)2020-05-07-23:14:22_00000-00255.txt")
        data2 = Data.load_ascii("$(data_remote)2020-05-07-23:14:22_00256-00446.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_remote)2020-05-07-23:14:22_00000-00446.txt")
        data1 = Data.load_ascii("$(data_remote)2020-05-30-22:04:58_00000-00255.txt")
        data2 = Data.load_ascii("$(data_remote)2020-05-30-22:04:58_00256-00446.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_remote)2020-05-30-22:04:58_00000-00446.txt")
        =#

        # Loading data
        data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        Data.zap!(data3; ranges=[241, 335, 334])
        data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        data6 = Data.load_ascii("$(data_dir)2020-05-30-22:04:58_00000-00446.txt")

        #data1p = Data.load_ascii("$(data_patrick)2019-09-29-12:41:18_00000-00294.txt")
        #Plot.single_J1750(data1p, outdir; start=1, number=nothing, name_mod="1p", darkness=0.7, show_=true, panel="a")

        # single pulse plots
        #=
        data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        Plot.single_J1750(data1, outdir; start=1, number=nothing, name_mod="1", darkness=0.7, show_=true, panel="a")
        data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        Plot.single_J1750(data2, outdir; start=1, number=nothing, name_mod="2", darkness=0.7, show_=true, panel="b")
        data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        Data.zap!(data3; ranges=[241, 335, 334])
        Plot.single_J1750(data3, outdir; start=1, number=nothing, name_mod="3", darkness=0.7, show_=true, panel="c")
        data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        Plot.single_J1750(data4, outdir; start=1, number=nothing, name_mod="4", darkness=0.7, show_=true, panel="d")
        data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        Plot.single_J1750(data5, outdir; start=1, number=nothing, name_mod="5", darkness=0.7, show_=true, panel="e")
        data6 = Data.load_ascii("$(data_dir)2020-05-30-22:04:58_00000-00446.txt")
        Plot.single_J1750(data6, outdir; start=1, number=nothing, name_mod="6", darkness=0.7, show_=true, panel="f")
        =#

        # get information
        #=
        files = Glob.glob("*_pipeline_info.json", data_dir)
        for file in files
            println("\n$file")
            js = JSON.parsefile(file)
            #println(js)
            println(length(split(js["input_data"]["header"]["ANTENNAE"], ",")))
            println(js["input_data"]["header"]["SCHEDULE_BLOCK_ID"])
            #println(js["snr_per_pulse"])
            println(js["snr"])
            println(js["input_data"]["header"]["TSAMP"])
            npulses = [295, 1031, 446, 450, 447, 447, 446]
            period = 0.6840110402274631
            duration = [period * np for np in npulses]
            println(duration)
        end
        =#

        # p3_evolution plots
        #=
        data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294.txt")
        Plot.p3_evolution_J1750(data1, outdir; panel="a", step=1, darkness=1.0, bin_st=450, bin_end=700, name_mod="1", number=128, verbose=true)
        readline(stdin)
        data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt")
        Plot.p3_evolution_J1750(data2, outdir; panel="b", step=1, darkness=1.0, bin_st=450, bin_end=700, name_mod="2", number=128, verbose=true)
        readline(stdin)
        data3 = Data.load_ascii("$(data_dir)2020-02-24-03:46:43_00000-00445.txt")
        Data.zap!(data3; ranges=[241, 335, 334])
        Plot.p3_evolution_J1750(data3, outdir; panel="c", step=1, darkness=1.0, bin_st=450, bin_end=700, name_mod="3", number=128, verbose=true)
        readline(stdin)
        data4 = Data.load_ascii("$(data_dir)2020-03-29-04:11:50_00000-00449.txt")
        Plot.p3_evolution_J1750(data4, outdir; panel="d", step=1, darkness=1.0, bin_st=450, bin_end=700, name_mod="4", number=128, verbose=true)
        readline(stdin)
        data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        Plot.p3_evolution_J1750(data5, outdir; panel="e", step=1, darkness=1.0, bin_st=450, bin_end=700, name_mod="5", number=128, verbose=true)
        readline(stdin)
        data6 = Data.load_ascii("$(data_dir)2020-05-30-22:04:58_00000-00446.txt")
        Plot.p3_evolution_J1750(data6, outdir; panel="f", step=1, darkness=1.0, bin_st=450, bin_end=700, name_mod="6", number=128, verbose=true)
        =#

        #p3fold plots
        #p3data1 = Data.load_ascii("$(data_dir)2019-09-29-12:41:18_00000-00294_norefine.p3fold")
        #p3data2 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.p3fold")
        #Plot.p3fold(p3data1, outdir; start=1, number=nothing, bin_st=350, bin_end=674, name_mod="1", darkness=1.0, cmap="viridis")
        #Plot.p3fold(p3data2, outdir; start=1, number=nothing, bin_st=350, bin_end=674, name_mod="2", darkness=1.0, cmap="viridis")
        #Plot.p3fold_two(p3data1, p3data2, outdir; start=1, number=nothing, bin_st=350, bin_end=674, name_mod="12", darkness=1.0, cmap="viridis")
        # test
        #p3data3 = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030_b.p3fold")
        #Plot.p3fold_two(p3data1, p3data3, outdir; start=1, number=nothing, bin_st=350, bin_end=674, name_mod="12b", darkness=1.0, cmap="viridis")
        # peaks detection
        #peaks = Tools.track_subpulses(p3data1, 18, thresh=0.4, thresh2=0.5, on_st=350, on_end=674)
        #Plot.group_tracks(p3data1, "$(data_dir)tracks1", peaks; start=1, number=22, bin_st=350, bin_end=674, name_mod="1p3_", darkness=0.6)
        #peaks2 = Tools.track_subpulses(p3data2, 18, thresh=0.4, thresh2=0.5, on_st=350, on_end=674)
        #Plot.group_tracks(p3data2, "$(data_dir)tracks2", peaks2; start=1, number=22, bin_st=350, bin_end=674, name_mod="2p3_", darkness=0.6)
        #Plot.p3fold_twotracks(p3data1, p3data2, "$(data_dir)tracks", outdir; start=1, number=nothing, bin_st=350, bin_end=674, name_mod="12", darkness=1.0, cmap="viridis")

        # P2 estimate
        #win = 12
        #Tools.p2_estimate(data1; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        #Tools.p2_estimate(data2; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        #Tools.p2_estimate(data3; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win) #  no! # too noisy!
        #Tools.p2_estimate(data4; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        #Tools.p2_estimate(data5; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win) # too noisy!
        #Tools.p2_estimate(data6; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)


        # track subpulses
        #peaks = Tools.track_subpulses(data1, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data1, "$(outdir)tracks1", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="1", darkness=0.6)
        #Plot.tracks_analysis3("$(outdir)tracks1"; lambda=10000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="1", darkness=0.6)
        #peaks = Tools.track_subpulses(data2, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data2, "$(outdir)tracks2", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis3("$(outdir)tracks2"; lambda=10000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #peaks = Tools.track_subpulses(data3, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data3, "$(outdir)tracks3", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="3", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks3"; lambda=10000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="3", darkness=0.6)
        #peaks = Tools.track_subpulses(data4, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data4, "$(outdir)tracks4", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="4", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks4"; lambda=1000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="4", darkness=0.6)
        #peaks = Tools.track_subpulses(data5, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data5, "$(outdir)tracks5", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="5", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks5"; lambda=1000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="5", darkness=0.6)
        #peaks = Tools.track_subpulses(data6, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data6, "$(outdir)tracks6", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="6", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks6"; lambda=1000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="6", darkness=0.6)
        #Plot.single_J1750(data3, outdir; start=100, number=200, bin_st=350, bin_end=650, name_mod="3test", darkness=0.4, show_=true, panel="tt")

        # some checks
        #Plot.tracks_analysis("$outdir/tracks2"; win=40, start=1, number=1031, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis2("$outdir/tracks2"; win=40, start=1, number=1031, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis3("$outdir/tracks2"; lambda=1000.0, start=1, number=1031, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis4("$outdir/tracks2"; start=1, number=1031, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)

        # the Plot - old version
        #Plot.driftrate_J1750(outdir; lambda=0.0003, show_=true)
        #Plot.driftrate_analysis_J1750(outdir; lambda=1000.0, show_=true)

        # average profiles
        #Plot.average_J1750([data1, data2, data3, data4, data5, data6], outdir; lambda=1000.0, bin_st=350, bin_end=650, name_mod="123456", show_=true)

        Plot.driftdirection_J1750([data1, data2, data3, data4, data5, data6], outdir; lambda=1000.0, bin_st=350, bin_end=650, name_mod="123456", show_=false)

    end


    function J1750_paper2()
        data_dir = "/home/szary/work/MeerTime/J1750/new_data/"
        data_patrick = "/home/szary/work/MeerTime/J1750/data_patrick/"
        outdir = "/home/szary/work/MeerTime/J1750/"

        # convert to txt
        #=
        Data.convert_psrfit_ascii("$(data_patrick)/20190929_124118.debase.hp", "$(data_patrick)/20190929_124118.debase.hp.txt")
        Data.convert_psrfit_ascii("$(data_patrick)/20191214_142212.debase.hp", "$(data_patrick)/20191214_142212.debase.hp.txt")
        Data.convert_psrfit_ascii("$(data_patrick)/20200224_034643.debase.hp", "$(data_patrick)/20200224_034643.debase.hp.txt")
        Data.convert_psrfit_ascii("$(data_patrick)/20200329_041150.debase.hp", "$(data_patrick)/20200329_041150.debase.hp.txt")
        Data.convert_psrfit_ascii("$(data_patrick)/20200507_231422.debase.hp", "$(data_patrick)/20200507_231422.debase.hp.txt")
        Data.convert_psrfit_ascii("$(data_patrick)/20200530_220458.debase.hp", "$(data_patrick)/20200530_220458.debase.hp.txt")
        Data.convert_psrfit_ascii("$(data_patrick)/20200625_212452.debase.hp", "$(data_patrick)/20200625_212452.debase.hp.txt")
        =#

        data1 = Data.load_ascii("$(data_patrick)20190929_124118.debase.hp.txt")
        data2 = Data.load_ascii("$(data_patrick)20191214_142212.debase.hp.txt")
        #data3b = Data.load_ascii("$(data_dir)2019-12-14-14:22:12_00000-01030.txt") # old for tests
        data3 = Data.load_ascii("$(data_patrick)20200224_034643.debase.hp.txt")
        Data.zap!(data3; ranges=[241, 335, 334])
        data4 = Data.load_ascii("$(data_patrick)20200329_041150.debase.hp.txt")
        data5 = Data.load_ascii("$(data_patrick)20200507_231422.debase.hp.txt")
        data6 = Data.load_ascii("$(data_patrick)20200530_220458.debase.hp.txt")
        data7 = Data.load_ascii("$(data_patrick)20200625_212452.debase.hp.txt")
        datas = [data1, data2, data3, data4, data5, data6, data7]


        #Plot.single_J1750(data2, outdir; start=1, number=nothing, name_mod="2p", darkness=0.7, show_=true, panel="b")

        # single pulse plots # not used!
        Plot.single_J1750(data1, outdir; start=1, number=nothing, name_mod="1p", darkness=0.5, show_=true, panel="")
        Plot.lrfs(data1, outdir; darkness=0.3, start=1, number=256, name_mod="J1750_1", bin_st=350, bin_end=650, change_fftphase=false)
        return
        #=
        Plot.single_J1750(data2, outdir; start=1, number=nothing, name_mod="2p", darkness=0.7, show_=true, panel="b")
        Plot.single_J1750(data3, outdir; start=1, number=nothing, name_mod="3p", darkness=0.7, show_=true, panel="c")
        Plot.single_J1750(data4, outdir; start=1, number=nothing, name_mod="4p", darkness=0.7, show_=true, panel="d")
        Plot.single_J1750(data5, outdir; start=1, number=nothing, name_mod="5p", darkness=0.7, show_=true, panel="e")
        Plot.single_J1750(data6, outdir; start=1, number=nothing, name_mod="6p", darkness=0.7, show_=true, panel="f")
        Plot.single_J1750(data7, outdir; start=1, number=nothing, name_mod="7p", darkness=0.7, show_=true, panel="g")
        =#

        # P2 estimate
        #=
        win = 8
        Tools.p2_estimate(data1; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        Tools.p2_estimate(data2; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        #Tools.p2_estimate(data3; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win) #  no! # too noisy!
        Tools.p2_estimate(data4; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        #Tools.p2_estimate(data5; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win) # too noisy!
        Tools.p2_estimate(data6; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        Tools.p2_estimate(data7; on_st=350, on_end=675, off_st=25, off_end=350, thresh=3, win=win)
        =#

        # calculate SNR files
        data_files = ["$(data_patrick)20190929_124118.debase.hp", "$(data_patrick)20191214_142212.debase.hp", "$(data_patrick)20200224_034643.debase.hp", "$(data_patrick)20200329_041150.debase.hp", "$(data_patrick)20200507_231422.debase.hp", "$(data_patrick)20200530_220458.debase.hp", "$(data_patrick)20200625_212452.debase.hp"]

        #Tools.generate_snr("input/single_pulses.gg") # for the modeled data
        #=
        for f in data_files
            Tools.generate_snr(f)
        end
        =#
        thresh, thresh2 = 5, 0.1
        #Tools.generate_snr(data_files[3])
        #peaks = Tools.track_subpulses_snr(data3, 18, data_files[3]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        #=
        peaks = Tools.track_subpulses_snr2(data1, 18, data_files[1]*".snr.txt", thresh=thresh, on_st=350, on_end=650)
        Plot.group_tracks(data1, "$(outdir)tracks/1", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="1", darkness=0.5)
        peaks = Tools.track_subpulses_snr(data2, 18, data_files[2]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        Plot.group_tracks(data2, "$(outdir)tracks/2", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="2", darkness=0.5)
        peaks = Tools.track_subpulses_snr(data3, 18, data_files[3]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        Plot.group_tracks(data3, "$(outdir)tracks/3", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="3", darkness=0.6)
        peaks = Tools.track_subpulses_snr(data4, 18, data_files[4]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        Plot.group_tracks(data4, "$(outdir)tracks/4", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="4", darkness=0.6)
        peaks = Tools.track_subpulses_snr(data5, 18, data_files[5]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        Plot.group_tracks(data5, "$(outdir)tracks/5", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="5", darkness=0.6)
        peaks = Tools.track_subpulses_snr(data6, 18, data_files[6]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        Plot.group_tracks(data6, "$(outdir)tracks/6", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="6", darkness=0.6)
        peaks = Tools.track_subpulses_snr(data7, 18, data_files[7]*".snr.txt", thresh=thresh, thresh2=thresh2, on_st=350, on_end=650)
        Plot.group_tracks(data7, "$(outdir)tracks/7", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="7", darkness=0.6)
        =#
        #=
        =#

        #Plot.test_track_subpulses_snr(data2, outdir, 18, data_files[2]*".snr.txt", thresh2=thresh2, bin_st=350, bin_end=650, name_mod="2") # OLD OBSOLETE
        Plot.test_track_subpulses_snr_new(data2, outdir, 18, data_files[2]*".snr.txt", bin_st=350, bin_end=650, name_mod="2") # in the paper - version 3

        # track subpulses # OLD # OBSOLOTE
        #peaks = Tools.track_subpulses(data1, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data1, "$(outdir)tracks/1", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="1", darkness=0.6)
        #Plot.tracks_analysis3("$(outdir)tracks/1"; lambda=10000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="1", darkness=0.6)
        #peaks = Tools.track_subpulses(data2, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data2, "$(outdir)tracks/2", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #Plot.tracks_analysis3("$(outdir)tracks/2"; lambda=10000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="2", darkness=0.6)
        #peaks = Tools.track_subpulses(data3, 15, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data3, "$(outdir)tracks/3", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="3", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks/3"; lambda=10000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="3", darkness=0.6)
        #peaks = Tools.track_subpulses(data4, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data4, "$(outdir)tracks/4", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="4", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks/4"; lambda=1000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="4", darkness=0.6)
        #peaks = Tools.track_subpulses(data5, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data5, "$(outdir)tracks/5", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="5", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks/5"; lambda=1000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="5", darkness=0.6)
        #peaks = Tools.track_subpulses(data6, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data6, "$(outdir)tracks/6", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="6", darkness=0.5)
        #Plot.tracks_analysis3("$(outdir)tracks/6"; lambda=1000.0, start=1, number=nothing, bin_st=350, bin_end=650, name_mod="6", darkness=0.6)
        #peaks = Tools.track_subpulses(data7, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(data7, "$(outdir)tracks/7", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="7", darkness=0.5)
        #Tools.convert_tracks("$outdir/tracks/") # converts from bins to longitude! # IMPORTANT # do it once # do not use it anymore! group_tracks improved

        # P3 grand-average NOPE! check Figure 2
        #=
        Plot.lrfs(data3, outdir; darkness=0.3, start=1, number=256, name_mod="J1750_2", bin_st=350, bin_end=650, change_fftphase=false)
        slices, pulses = random_slices(datas, 256)
        Plot.lrfs_average_J1750(slices, pulses, outdir; darkness=1.0, bin_st=350, bin_end=650, name_mod="1234567", number=128, verbose=true)
        println(size(slices))
        println(typeof(pulses[1]))
        for (i, k) in pulses[1]
            println("$i $k")
        end
        =#

        # New single pulses plot
        #Plot.singlepulses_J1750(datas, outdir; name_mod="1234567", darkness=0.5, show_=false)

        # New P3 evolution plot
        #Plot.p3evolutions_J1750(datas, outdir, 1, 128; bin_st=350, bin_end=650, name_mod="1234567", darkness=0.7, show_=true)
        # get P3 values
         #Plot.p3_evolution_J1750(data7, outdir; panel="g", step=1, darkness=1.0, bin_st=350, bin_end=650, name_mod="7", number=128, verbose=true)
         #Plot.p3_evolution_J1750(data2, outdir; panel="b", step=1, darkness=1.0, bin_st=350, bin_end=650, name_mod="2", number=128, verbose=true)

        #The P3-folded plot
        #p3data1 = Data.load_ascii("$(data_patrick)20190929_124118.debase.hp.p3fold")
        #p3data1a = Data.load_ascii("$(data_patrick)20190929_124118.debase.hp.norefine.p3fold")
        #p3data1b = Data.load_ascii("$(data_patrick)20190929_124118.debase.hp.refine.p3fold")
        #p3data2 = Data.load_ascii("$(data_patrick)20191214_142212.debase.hp.p3fold")
        #Plot.p3fold(p3data1, outdir; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="1", darkness=1.0, cmap="viridis")
        #Plot.p3fold_two(p3data1, p3data2, outdir; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="12", darkness=0.7, cmap="viridis")
        #peaks = Tools.track_subpulses(p3data1, 18, thresh=0.4, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(p3data1, "$(data_patrick)tracks/1", peaks; start=1, number=44, bin_st=350, bin_end=650, name_mod="1p3_", darkness=0.7)
        #peaks = Tools.track_subpulses(p3data2, 18, thresh=0.4, thresh2=0.5, on_st=350, on_end=650)
        #Plot.group_tracks(p3data2, "$(data_patrick)tracks/2", peaks; start=1, number=44, bin_st=350, bin_end=650, name_mod="2p3_", darkness=0.7)
        #Plot.p3fold_twotracks(p3data1, p3data2, "$(data_patrick)tracks", outdir; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="12", darkness=0.7, cmap="viridis", show_=true)

        # The Plot
        #Plot.driftrate_J1750_2(outdir; spar=0.6, show_=true) # not too good, but I tried...
        #Plot.driftrate_J1750_2(outdir; lambda=200.0, show_=true) # in the paper
        #Plot.driftrate_J1750_2(outdir; lambda=100.0, show_=true) # tests
        #Plot.driftrate_analysis_J1750_2(outdir; lambda=200.0, show_=false, datas=[data1, data2, data3, data4, data5, data6, data7]) # values in the paper
        #Plot.driftrate_J1750_3(outdir; lambda=200.0, show_=true) # segmented fits # not used?

        # Timescales
        #Plot.driftdirection_J1750_2([data1, data2, data3, data4, data5, data6], outdir; lambda=1000.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true)
        # Timescales - new
        #Plot.driftdirection_J1750_3([data1, data2, data3, data4, data5, data6], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true) # in paper up to version 2
        #Plot.driftdirection_J1750_3b([data1, data2, data3, data4, data5, data6], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true) # in the paper version 3
        #Plot.driftdirection_J1750_4([data1, data2, data3, data4, data5, data6], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true) # in the paper version 4

        #Plot.driftdirection_analysis_J1750_3(outdir; lambda=200.0, show_=false) # TODO
        #Plot.driftdirection_analysis_J1750_4(outdir; lambda=200.0, show_=true) # changes vs longitude
        #Plot.p2_analysis_J1750(outdir; lambda=200.0, show_=true) # p2 stability

        # average profiles
        #Plot.average_J1750([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true) # old too many ranges
        #Plot.average_J1750_2([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true)
        #Plot.average_J1750_3([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true) # in the paper

        # DO NOT USE THESE
        #Plot.average_J1750_4([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true)
        #Plot.average_J1750_stability([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=nothing, bin_end=nothing, name_mod="1234567", show_=true) # 350 650
    end


    function J1750_average()

        outdir = "/home/szary/work/MeerTime/J1750/"
        data_dir2 = "/home/szary/work/MeerTime/J1750/new_data2/"

        # connect data
        #=
        data1 = Data.load_ascii("$(data_dir2)2019-12-14-14:22:12_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2019-12-14-14:22:12_00256-00511.txt")
        data3 = Data.load_ascii("$(data_dir2)2019-12-14-14:22:12_00512-00767.txt")
        data4 = Data.load_ascii("$(data_dir2)2019-12-14-14:22:12_00768-01030.txt")
        data = vcat(data1, data2, data3, data4)
        Data.save_ascii(data, "$(data_dir2)2019-12-14-14:22:12_00000-01030.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-02-24-03:46:43_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-02-24-03:46:43_00256-00445.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-02-24-03:46:43_00000-00445.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-03-29-04:11:50_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-03-29-04:11:50_00256-00449.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-03-29-04:11:50_00000-00449.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-05-07-23:14:22_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-05-07-23:14:22_00256-00446.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-05-07-23:14:22_00000-00446.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-05-30-22:04:58_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-05-30-22:04:58_00256-00446.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-05-30-22:04:58_00000-00446.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-06-25-21:24:52_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-06-25-21:24:52_00256-00445.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-06-25-21:24:52_00000-00445.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-07-25-17:53:20_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-07-25-17:53:20_00256-00441.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-07-25-17:53:20_00000-00441.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-08-23-16:06:06_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-08-23-16:06:06_00256-00444.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-08-23-16:06:06_00000-00444.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-09-19-14:14:25_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-09-19-14:14:25_00256-00445.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-09-19-14:14:25_00000-00445.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-11-16-12:47:10_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-11-16-12:47:10_00256-00445.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-11-16-12:47:10_00000-00445.txt")
        data1 = Data.load_ascii("$(data_dir2)2020-12-12-08:36:03_00000-00255.txt")
        data2 = Data.load_ascii("$(data_dir2)2020-12-12-08:36:03_00256-00445.txt")
        data = vcat(data1, data2)
        Data.save_ascii(data, "$(data_dir2)2020-12-12-08:36:03_00000-00445.txt")
        =#

        data1 = Data.load_ascii("$(data_dir2)2019-09-29-12:41:18_00000-00294.txt")
        data2 = Data.load_ascii("$(data_dir2)2019-12-14-14:22:12_00000-01030.txt")
        data3 = Data.load_ascii("$(data_dir2)2020-02-24-03:46:43_00000-00445.txt")
        data4 = Data.load_ascii("$(data_dir2)2020-03-29-04:11:50_00000-00449.txt")
        data5 = Data.load_ascii("$(data_dir2)2020-05-07-23:14:22_00000-00446.txt")
        data6 = Data.load_ascii("$(data_dir2)2020-05-30-22:04:58_00000-00446.txt")
        data7 = Data.load_ascii("$(data_dir2)2020-06-25-21:24:52_00000-00445.txt")
        data8 = Data.load_ascii("$(data_dir2)2020-07-25-17:53:20_00000-00441.txt")
        data9 = Data.load_ascii("$(data_dir2)2020-08-23-16:06:06_00000-00444.txt")
        data10 = Data.load_ascii("$(data_dir2)2020-09-19-14:14:25_00000-00445.txt")
        data11 = Data.load_ascii("$(data_dir2)2020-11-16-12:47:10_00000-00445.txt")
        data12 = Data.load_ascii("$(data_dir2)2020-12-12-08:36:03_00000-00445.txt")

        datas = [data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12]
        #datas = [data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data12]
        # average profile for almost whole data
        # peaks detection
        #=
        peaks = Tools.track_subpulses(data8, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        Plot.group_tracks(data8, "$(outdir)tracks/8", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="8", darkness=0.6)
        peaks = Tools.track_subpulses(data9, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        Plot.group_tracks(data9, "$(outdir)tracks/9", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="9", darkness=0.6)
        peaks = Tools.track_subpulses(data10, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        Plot.group_tracks(data10, "$(outdir)tracks/10", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="10", darkness=0.6)
        peaks = Tools.track_subpulses(data11, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        Plot.group_tracks(data11, "$(outdir)tracks/11", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="11", darkness=0.6)
        peaks = Tools.track_subpulses(data12, 18, thresh=0.5, thresh2=0.5, on_st=350, on_end=650)
        Plot.group_tracks(data12, "$(outdir)tracks/12", peaks; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="12", darkness=0.6)
        =#

        #Plot.average_J1750_stability([data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12], outdir; lambda=200.0, bin_st=nothing, bin_end=nothing, name_mod="12", show_=true) # 350 650
        #Plot.average_J1750_2(datas, outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="12", show_=true)
        #Plot.average_J1750_2(datas, outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="12", show_=true)
        #Plot.average_J1750_3(datas, outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="12", show_=true)
        Plot.average_J1750_3(datas, outdir; lambda=200.0, bin_st=390, bin_end=600, name_mod="12", show_=true)
    end


    function J1750_modeled()

        data_patrick = "/home/szary/work/MeerTime/J1750/data_patrick/"
        data_modeled= "/home/szary/work/MeerTime/J1750/modeled/"
        outdir = "/home/szary/work/MeerTime/J1750/modeled/"

        nsp = 13

        data1 = Data.load_ascii("$(data_modeled)single_pulses_6_212.ascii")
        data1m = Data.load_ascii("$(data_modeled)single_pulses_218.ascii")
        data1mnew = Data.load_ascii("$(data_modeled)single_pulses_13_247.ascii")
        data1mnewv3 = Data.load_ascii("$(data_modeled)single_pulses_13_257.ascii")
        data1new = Data.load_ascii("$(data_modeled)single_pulses_13_245.ascii") # Łokej
        #data2 = Data.load_ascii("$(data_modeled)single_pulses_noiseless.ascii")
        data1p = Data.load_ascii("$(data_patrick)20190929_124118.debase.hp.txt")
        data2p = Data.load_ascii("$(data_patrick)20191214_142212.debase.hp.txt")

        #peaks = Tools.track_subpulses(data1, 18, thresh=0.5, thresh2=0.5, on_st=600, on_end=900, off_st=200, off_end=500)
        #Plot.group_tracks(data1, "$(outdir)tracks/", peaks; start=1, number=nothing, bin_st=600, bin_end=900, name_mod="1", darkness=0.6)
        #Plot.driftrate_J1750_2(outdir; lambda=200.0, show_=true)

        #Plot.single(data1, outdir; darkness=0.5, bin_st=330, bin_end=400, start=1, number=size(data1p)[1], name_mod="model_$(nsp)_", show_=false)
        #Plot.single(data1, outdir; darkness=0.5, bin_st=285, bin_end=320, start=1, number=size(data1p)[1], name_mod="model_$(nsp)_", show_=false)
        #Plot.lrfs(data1, outdir; darkness=0.5, bin_st=330, bin_end=400, start=1, name_mod="model_$(nsp)_", change_fftphase=false)

        #Plot.single(data1, outdir; darkness=0.8, bin_st=600, bin_end=680, start=1, number=size(data1p)[1], name_mod="model_$(nsp)_", show_=false)
        #Plot.lrfs(data1, outdir; darkness=0.5, bin_st=600, bin_end=680, start=1, name_mod="model_$(nsp)_", change_fftphase=false)

        #Plot.single(data1p, outdir; darkness=0.5, bin_st=425, bin_end=565, start=1, number=nothing, name_mod="obs")
        #Plot.lrfs(data1p, outdir; darkness=0.5, bin_st=445, bin_end=545, start=1, name_mod="obs", change_fftphase=false)

        #Plot.lrfses([data1p, data1], outdir; darkness=[0.5, 0.87], bin_st=[400, 270], bin_end=[600, 470], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)")

        # not anymore in the paper alpha > 90
        #nsp = 6
        #Plot.singles([data1p, data1], outdir; darkness=[0.5, 0.87], bin_st=[400, 270], bin_end=[600, 470], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)", show_=false)
        #Plot.lrfses([data1p, data1], outdir; darkness=[0.5, 0.3], bin_st=[420, 290], bin_end=[570, 440], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)")
        #Plot.modeledpulses_J1750([data2p, data1m], outdir; darkness=[0.7, 1.0], name_mod="obs_model_218", show_=false)

        # in the paper! # alpha < 90
        #nsp = 13
        #Plot.singles([data1p, data1new], outdir; darkness=[0.5, 0.87], bin_st=[300, 245], bin_end=[700, 645], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)", show_=false)
        #Plot.singles([data1p, data1new], outdir; darkness=[0.5, 0.87], bin_st=[400, 345], bin_end=[600, 545], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)", show_=false) # this one (f10.pdf)
        #Plot.lrfses([data1p, data1new], outdir; darkness=[0.5, 0.3], bin_st=[420, 380], bin_end=[570, 510], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)")

        #Plot.modeledpulses_J1750([data2p, data1mnew], outdir; darkness=[0.7, 1.0], name_mod="obs_model_247", show_=false) # f11.pdf v.1-2
        Plot.modeledpulses_J1750([data2p, data1mnewv3], outdir; darkness=[0.7, 1.0], name_mod="obs_model_257", show_=false) # f11.pdf v.3

    end


    function random_slices(datas, len)
        slices = []
        pulses = []
        for i in 1:size(datas)[1]
            #println(size(datas[i]))
            indexes = collect(1:size(datas[i])[1]) # pulses
            trys = 0
            max_try = 100
            while (length(indexes) > len) && (trys < max_try)
                mx_ind = length(indexes) - len
                ind = rand(1:mx_ind, 1)[1]
                j = indexes[ind]
                k = indexes[ind+len]
                #println(abs(k-j))
                if abs(k-j) == len
                    splice!(indexes, ind:ind+len-1) #  remove pulses
                    #println("$j, $k, $(j-k)")
                    #println(indexes)
                    push!(slices, datas[i][j:k-1, :]) # collect data
                    push!(pulses, Dict(i=>(j, k-1))) # collect data
                    #println(size(slices[1]))
                end
                #return
                trys += 1
                #println(trys)
            end
            #println(indexes)
        end
        return slices, pulses
    end

    function J1750_calculations()
        # P3 Rahul
        edot0 = (2.3 + 0.2) * 1e32 # err 0.2
        edot = 4.7e30
        p3 = (edot / edot0) ^ (-0.6 - 0.1) # err 0.1
        println("Rahul P3: ", p3)

        # RS1975
        edot1 = 4e31
        nsp = 13
        p3_rs = (5.6 / nsp) * (edot / edot1) ^ 0.5
        println("RS 1975 P3: ", p3_rs)

        # aliasing
        n = 7
        p3_obs = 43.5 # err 0.4
        p3_tr = 1 / (n + 1 / (p3_obs))  # Gupta formula p3_obs in P_3
        println("True P3: ", p3_tr)

    end

    function test(outdir)

        d = Data.load_ascii("input/1.txt")
        Plot.single(d, outdir; darkness=0.3, number=256, bin_st=200, bin_end=800, start=1, name_mod="1", show_=true)
        Plot.average(d, outdir; number=256, bin_st=200, bin_end=800, start=1, name_mod="1", show_=true)
        Plot.lrfs(d, outdir; darkness=0.1, start=1, name_mod="1", bin_st=500, bin_end=530, show_=true)

    end

    function J1319(outdir)    
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00000-00255.spCF", outdir*"1.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00256-00511.spCF", outdir*"2.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00512-00767.spCF", outdir*"3.txt")  
        Data.convert_psrfit_ascii("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/2019-12-15-03:19:04_00768-01029.spCF", outdir*"4.txt")  
        data1 = Data.load_ascii(outdir*"1.txt")
        data2 = Data.load_ascii(outdir*"2.txt")
        data3 = Data.load_ascii(outdir*"3.txt")
        data4 = Data.load_ascii(outdir*"4.txt")
        data = vcat(data1, data2, data3, data4)
        #Plot.single(data, outdir; darkness=0.5, number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        #Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J1319", show_=true)
        #Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J1319", bin_st=400, bin_end=600, show_=true)
        #Plot.p3fold(data, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
    end

    function J0151(outdir)
        Data.convert_psrfit_ascii("/home/psr/data/new/J0151-0635/2020-04-13-10:00:14/2020-04-13-10:00:14_00000-00255.spCF", outdir*"1.txt") 
        Data.convert_psrfit_ascii("/home/psr/data/new/J0151-0635/2020-04-13-10:00:14/2020-04-13-10:00:14_00256-00511.spCF", outdir*"2.txt") 
        Data.convert_psrfit_ascii("/home/psr/data/new/J0151-0635/2020-04-13-10:00:14/2020-04-13-10:00:14_00512-00767.spCF", outdir*"3.txt") 
        Data.convert_psrfit_ascii("/home/psr/data/new/J0151-0635/2020-04-13-10:00:14/2020-04-13-10:00:14_00768-01038.spCF", outdir*"4.txt")
        data1 = Data.load_ascii(outdir*"1.txt")
        data2 = Data.load_ascii(outdir*"2.txt")
        data3 = Data.load_ascii(outdir*"3.txt")
        data4 = Data.load_ascii(outdir*"4.txt")
        data = vcat(data1, data2, data3, data4)
        #Plot.single(data, outdir; darkness=0.5, number=150, bin_st=400, bin_end=600, start=1, name_mod="J0151", show_=true)
        #Plot.average(data, outdir; number=nothing, bin_st=400, bin_end=600, start=1, name_mod="J0151", show_=true)
        #Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="J0151", bin_st=400, bin_end=600, show_=true)
        Plot.p3fold(data, outdir; start=1, number=nothing, repeat_num=4, cmap="inferno", bin_st=nothing, bin_end=nothing, darkness=0.5, name_mod="0")
    end


    function J1319_psrchive(indir, outdir)
        Data.process_psrchive(indir, outdir, ["2019-12-15-03:19:04_00000-00255.spCF", "2019-12-15-03:19:04_00256-00511.spCF", "2019-12-15-03:19:04_00512-00767.spCF", "2019-12-15-03:19:04_00768-01029.spCF"], "J1319-6105.spCF")
    end



    function process_psrfit_files(base_dir::String, output_dir::String; name_mod::Union{String, Nothing}=nothing)
        # Step 1: Extract base directory name
        base_name = basename(base_dir)
    
        # Step 2: Auto-generate name_mod if not provided
        name_mod = isnothing(name_mod) ? base_name * "Mac" : name_mod
    
        # Step 3: Create output subdirectory
        output_subdir = joinpath(output_dir, base_name)
        if !isdir(output_subdir)
            mkpath(output_subdir)
            println("Created output directory: ", output_subdir)
        end
    
        # Step 4: Find second-level catalogue
        second_catalogue = joinpath(base_dir, readdir(base_dir)[1])
        println("Second catalogue found: ", second_catalogue)
    
        # Step 5: Find .spCF files
        spcf_files = filter(f -> occursin("spCF", f), readdir(second_catalogue, join=true))
        converted_txt_files = String[]
    
        try
            # Step 6: Combining .spCF into one 
            output_file = joinpath(output_subdir, "converted.spCF")
            file_names = [joinpath(base_name, file) for file in spcf_files]
            run(pipeline(`psradd $file_names -o $output_file`, stderr="errs.txt"))
            out_txt = replace(output_file, ".spCF" => ".txt")
    
            # Step 7: Convert spCF -> ascii
            Data.convert_psrfit_ascii(output_file, out_txt)
    
            # Step 8: Load combined data
            combined_data = Data.load_ascii(out_txt)
        
            # Step 9: Plot
            Plot.single(combined_data, output_subdir, darkness=0.5, bin_st=1, bin_end=1024, number=nothing, name_mod=name_mod, show_=false)
            Plot.lrfs(combined_data, output_subdir, darkness=0.1, start=1, bin_st=1, bin_end=1024, name_mod=name_mod, change_fftphase=false, show_=false)
            Plot.average(combined_data, output_subdir, bin_st=1, bin_end=1024, number=nothing, name_mod=name_mod, show_=false)
            
            # Step 10: P3 Fold
            folded = Tools.p3fold(combined_data, 4.81, 24)
            Plot.p3fold(folded, output_subdir, start=3, bin_st=470, bin_end=550, name_mod=name_mod, show_=false, repeat_num=4)
        catch e
            println("Error encountered: ", e)
            println("Skipping this dataset and continuing...")
        end
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
    
    # Run processing for all catalogues
    function J0034Mac(output_dir)
        process_all_catalogues(output_dir, "/home/psr/data/new")
    end
    
    function process_psrdata(indir, outdir)
        bin_st, bin_end, pulsar_name = Data.process_psrdata(indir, outdir)
        folded = Data.load_ascii(joinpath(outdir, pulsar_name, "pulsar.debase.p3fold"))
        Plot.p3fold(folded, joinpath(outdir, pulsar_name); start=3, bin_st=bin_st, bin_end=bin_end, name_mod="test", show_=true, repeat_num=4)
        
    end

    function J1750_psrdata(indir, outdir)
        bin_st, bin_end = Data.process_psrdata(indir, outdir; files= ["J1750-3503_MeerKAT_2019-09-29.debase.hp"])
        folded = Data.load_ascii(outdir*"/pulsar.debase.p3fold")
        Plot.p3fold(folded, outdir; start=3, bin_st=bin_st, bin_end=bin_end, name_mod="test", show_=true, repeat_num=4)

    end

    function fold_test(indir, outdir)
        # TODO fix this!
        Data.convert_psrfit_ascii("J1750-3503_MeerKAT_2019-09-29.debase.hp", outdir*"1.txt")
        data = Data.load_ascii(outdir*"1.txt")
        #Plot.single(data, outdir; darkness=0.5, number=nothing, bin_st=470, bin_end=550, start=1, name_mod="test", show_=true)
        #Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="test", bin_st=470, bin_end=550, show_=true)
        folded = Tools.p3fold(data, 4.81, 24)
        #println(size(folded))
        Plot.single(data, outdir; darkness=0.5, number=100, bin_st=470, bin_end=550, start=1, name_mod="test", show_=true)
        Plot.lrfs(data, outdir; darkness=0.1, start=1, name_mod="test", bin_st=470, bin_end=550, show_=true)
        Plot.p3fold(folded, outdir; start=3, bin_st=470, bin_end=550, name_mod="test", show_=true, repeat_num=4)
    end

    function print_lrfs_header_from_folder(folder::String)
    
        # Resolve the ~ to the full absolute path
        folder = abspath(expanduser(folder))
    
        # Construct the full path to the pulsar.debase.lrfs file
        fits_path = joinpath(folder, "pulsar.debase.lrfs")
    
        # Check if the file exists
        if !isfile(fits_path)
            println("File does not exist: $fits_path")
            return
        end
    
        println("Opening FITS file: $fits_path")
    
        # Open the FITS file
        f = FITS(fits_path)
    
        # Loop through each HDU (Header/Data Unit) and print its header
        for (i, hdu) in enumerate(f)
            println("=== HEADER HDU $i ===")
            try
                if typeof(hdu) <: ImageHDU
                    # For ImageHDU, get the header with `hdu.header`
                    header = read_header(hdu)
                elseif typeof(hdu) <: TableHDU
                    # For TableHDU, get the header with `hdu.header`
                    header = read_header(hdu)
                else
                    println("Unsupported HDU type: $(typeof(hdu))")
                    continue
                end
    
                # Print the header keys and values
                for (key, value) in header
                    println("  $key = $value")
                end
            catch e
                println("Error reading header for HDU $i: $e")
            end
            println()
        end
    
        close(f)
    end
    

    
    



 
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
            im = ax.imshow(data;
                extent=[160, 200, 0, 0.5],  # Zamieniono osie x i y
                origin="lower",
                aspect="auto",
                cmap="gray",  # Dodano przecinek
                #vmin=0,  # Minimalna wartość skali kolorów
                #vmax=0.06  # Maksymalna wartość skali kolorów
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
    
    
    
    
    
    
    


    
    
    
    
    function print_first_10_lines(filepath::String)
        println("Reading first 10 lines from FITS file: $filepath")
        try
            f = open(filepath, "r")  # Open the FITS file in read mode
            for i in 1:10
                line = readline(f)  # Read a single line
                println(line)  # Print the line
            end
            close(f)  # Close the file after reading
        catch e
            println("Error reading FITS file: $e")  # Handle any errors
        end
    end
    
 
    
    
    
    
    

    


    



    
    
    
    function main()
        # output directory for local run
        localout = "output"
        # output directory for VPM
        vpmout = "/home/psr/output/"
        indir = "/home/psr/data/"

        #J1319(vpmout)
        #process_psrdata("/home/psr/data/new/J1319-6105/2019-12-15-03:19:04/", vpmout)
        #process_psrdata("/home/psr/data/new/J1919+0134/2020-02-02-11:45:29/", vpmout)
        #process_psrdata("/home/psr/data/new/J1057-5226/2019-06-21-15:37:29", vpmout)
        #print_lrfs_header_from_folder("~/output/J1919+0134")

        plot_2dfs("/home/psr/output", "J1057-5226", show_plot=true)
        #inspect_fits("/home/psr/output/J1919+0134/pulsar.debase.1.2dfs")
        #print_first_10_lines("/home/psr/output/J1057-5226/pulsar.debase.1.2dfs")


        #J1750_psrdata(indir, vpmout)
        #fold_test(indir, vpmout)
        #test(vpmout)
        #J0820(args)
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
