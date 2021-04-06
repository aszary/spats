module SpaTs
    using ArgParse
    using Glob
    using JSON

    include("modules/data.jl")
    include("modules/plot.jl")
    include("modules/tools.jl")

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

        # single pulse plots # not used!
        #=
        Plot.single_J1750(data1, outdir; start=1, number=nothing, name_mod="1p", darkness=0.7, show_=true, panel="a")
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

        # track subpulses
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
        #Plot.singlepulses_J1750(datas, outdir; name_mod="1234567", darkness=0.5, show_=true)

        # New P3 evolution plot
        #Plot.p3evolutions_J1750(datas, outdir, 1, 128; name_mod="1234567", darkness=0.7, show_=true)
        # get P3 values
        # Plot.p3_evolution_J1750(data7, outdir; panel="g", step=1, darkness=1.0, bin_st=350, bin_end=650, name_mod="7", number=128, verbose=true)

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
        #Plot.p3fold_twotracks(p3data1, p3data2, "$(data_patrick)tracks", outdir; start=1, number=nothing, bin_st=350, bin_end=650, name_mod="12", darkness=0.7, cmap="viridis")

        # The Plot
        #Plot.driftrate_J1750_2(outdir; spar=0.6, show_=true) # not too good, but I tried...
        #Plot.driftrate_J1750_2(outdir; lambda=200.0, show_=true) # in the paper
        #Plot.driftrate_analysis_J1750_2(outdir; lambda=200.0, show_=false)
        #Plot.driftrate_J1750_3(outdir; lambda=200.0, show_=true) # segmented fits # not used?
        #Plot.driftrate_analysis_J1750_3(outdir; lambda=200.0, show_=false) # TODO

        # Timescales # too messy switch to drift rate?
        #Plot.driftdirection_J1750_2([data1, data2, data3, data4, data5, data6], outdir; lambda=1000.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true)
        # Timescales - new
        Plot.driftdirection_J1750_3([data1, data2, data3, data4, data5, data6], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=false)

        # average profiles
        #Plot.average_J1750([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true) # old too many ranges
        #Plot.average_J1750_2([data1, data2, data3, data4, data5, data6, data7], outdir; lambda=200.0, bin_st=350, bin_end=650, name_mod="1234567", show_=true)


    end


    function J1750_modeled()

        data_patrick = "/home/szary/work/MeerTime/J1750/data_patrick/"
        data_modeled= "/home/szary/work/MeerTime/J1750/modeled/"
        outdir = "/home/szary/work/MeerTime/J1750/modeled/"

        nsp = 6

        data1 = Data.load_ascii("$(data_modeled)single_pulses_6_212.ascii")
        data1m = Data.load_ascii("$(data_modeled)single_pulses_218.ascii")
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

        # in the paper!
        #Plot.singles([data1p, data1], outdir; darkness=[0.5, 0.87], bin_st=[400, 270], bin_end=[600, 470], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)", show_=false)
        #Plot.lrfses([data1p, data1], outdir; darkness=[0.5, 0.3], bin_st=[420, 290], bin_end=[570, 440], start=1, number=size(data1p)[1], name_mod="obs_model_$(nsp)")

        Plot.modeledpulses_J1750([data2p, data1m], outdir; darkness=[0.7, 1.0], name_mod="obs_model_218", show_=false)

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


    function main()
        args = parse_commandline()
        for (arg, val) in args
            println("  $arg  =>  $val")
        end

        #mkieth()
        #J1651()
        #J1705()
        #B0320()
        #J1750_remote()
        #J1750_local()
        #J1750_paper()
        J1750_paper2()
        #J1750_modeled()

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
