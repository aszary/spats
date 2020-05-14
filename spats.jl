module SpaTs
    using ArgParse

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
        Data.convert_psrfit_ascii("$(ozdir5)2020-05-07-23:14:22_00000-00255.spCF", "$(ozdir5)2020-05-07-23:14:22_00000-00255.txt")
        Data.convert_psrfit_ascii("$(ozdir5)2020-05-07-23:14:22_00256-00446.spCF", "$(ozdir5)2020-05-07-23:14:22_00256-00446.txt")


    end


    function J1750_local()
        # local data
        data_dir = "/home/szary/work/MeerTime/J1750/data/"
        outdir = "/home/szary/work/MeerTime/J1750/"


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
        #Plot.single(data1, outdir; start=1, number=nothing, bin_st=nothing, bin_end=nothing, name_mod="J1750_test", darkness=0.3)
        #Tools.p2_estimate(data1; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=8) # 18
        #Tools.p2_estimate(data2; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=4, template_num=37) # 18
        #Tools.p2_estimate(data3; on_st=350, on_end=650, off_st=50, off_end=350, thresh=3.5, win=15) #  no! # to noisy!
        #Tools.p2_estimate(data4; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=12) # 17.6

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
        data5 = Data.load_ascii("$(data_dir)2020-05-07-23:14:22_00000-00446.txt")
        #Data.zap!(data3; ranges=[241, 335, 334])
        #Plot.single(data1, outdir; start=1, number=10, bin_st=400, bin_end=700, name_mod="short", darkness=0.6)

        #peaks = Tools.track_subpulses(data1, 18.0, thresh=2.1, thresh2=0.7, on_st=500, on_end=650)
        #Plot.tracks(data1, outdir, peaks; start=1, number=294, bin_st=500, bin_end=650, name_mod="J1750_1_new", darkness=0.6)
        #peaks = Tools.track_subpulses(data2, 18.0, thresh=2.2, thresh2=0.8, on_st=500, on_end=650)
        #Plot.tracks(data2, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_2", darkness=0.6)
        #peaks = Tools.track_subpulses(data3, 18.0, thresh=2.2, thresh2=0.8, on_st=425, on_end=600)
        #Plot.tracks(data3, outdir, peaks; start=1, number=256, bin_st=425, bin_end=600, name_mod="J1750_3", darkness=0.6)
        #peaks = Tools.track_subpulses(data4, 18.0, thresh=2.1, thresh2=0.8, on_st=500, on_end=650)
        #Plot.tracks(data4, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_4", darkness=0.6)

        #p2, template = Tools.p2_estimate(data2; on_st=450, on_end=700, off_st=100, off_end=350, thresh=5, win=4, template_num=37) # 18
        #peaks = Tools.track_subpulses(data2, p2, thresh=2.1, thresh2=0.7, on_st=500, on_end=650, template=nothing)
        #Plot.tracks(data2, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_2_new", darkness=0.6)

        peaks = Tools.track_subpulses(data5, 18, thresh=2.1, thresh2=0.7, on_st=500, on_end=650, template=nothing)
        Plot.tracks(data5, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_5", darkness=0.6)
        Plot.tracks(data5, outdir, peaks; start=1, number=256, bin_st=500, bin_end=650, name_mod="J1750_5", darkness=0.6)


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
        J1750_local()

    end

    function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table s begin
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
