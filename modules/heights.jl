module Heights
    using Glob
    using FITSIO
    using ProgressMeter
    using CairoMakie, FileIO
    using PDFIO
    using Statistics
    using FFTW
    using DelimitedFiles

    
    include("functions.jl")
    include("tools.jl")
    include("plot.jl")
    include("data.jl")
    include("GaussianFit.jl")



