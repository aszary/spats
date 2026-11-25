# Unit tests for Basu & Mitra (2018) Fig 4 peak tracking implementation
using Test
using Statistics
using FFTW

# Include required modules if not already defined
if !isdefined(@__MODULE__, :Plot)
    include("plot.jl")
end

function test_basu_fig4_synthetic()
    println("Running test_basu_fig4_synthetic...")
    # Synthetic parameters
    N_periods = 1000
    N_bins = 64
    bin_st = 20
    bin_end = 45
    n_on = bin_end - bin_st + 1

    # Frequencies matching Basu & Mitra (2018) J1555-3134
    fp1 = 0.057
    fp2 = 0.098

    # Create synthetic data matrix: pulses x bins
    data = zeros(Float64, N_periods, N_bins)

    # In on-pulse region, add two modulated sinusoidal signals:
    # Bin 25 has Peak 1, Bin 35 has Peak 2, Bin 40 has both, Bin 10 (off-pulse) has pure noise
    for n in 1:N_periods
        # Peak 1 on bin 25
        data[n, 25] += 5.0 * cos(2π * fp1 * n)
        # Peak 2 on bin 35
        data[n, 35] += 4.0 * cos(2π * fp2 * n)
        # Both on bin 40
        data[n, 40] += 3.0 * cos(2π * fp1 * n) + 2.5 * cos(2π * fp2 * n)
    end

    # Add Gaussian noise across all bins (σ = 0.3)
    rng_noise = 0.3 .* randn(N_periods, N_bins)
    data .+= rng_noise

    peaks = [
        (name="Peak 1", fp=0.057, fwhm=0.027, color="red"),
        (name="Peak 2", fp=0.098, fwhm=0.025, color="green")
    ]

    # Run sliding FFT peak extraction
    res = Plot.find_sliding_lrfs_peaks(data, bin_st, bin_end;
                                       nfft=256, shift=50, nsigma=3.0, peaks=peaks)

    amps = res.amps
    valid = res.valid
    b_idx_25 = 25 - bin_st + 1
    b_idx_35 = 35 - bin_st + 1

    # 1. Check detection of Peak 1 in bin 25
    @test count(valid[1, :, b_idx_25]) > 0.8 * size(amps, 2)
    # Peak 1 amplitude should be recovered significantly above zero
    mean_amp_p1 = mean(amps[1, valid[1, :, b_idx_25], b_idx_25])
    @test mean_amp_p1 > 3.0

    # 2. Check detection of Peak 2 in bin 35
    @test count(valid[2, :, b_idx_35]) > 0.8 * size(amps, 2)
    mean_amp_p2 = mean(amps[2, valid[2, :, b_idx_35], b_idx_35])
    @test mean_amp_p2 > 2.5

    # 3. Check 3-sigma noise rejection: on a bin with pure noise (e.g. bin 21 inside on-pulse)
    b_idx_noise = 21 - bin_st + 1
    # Pure noise should have 0 or very few false positive detections (< 10%)
    false_pos_1 = count(valid[1, :, b_idx_noise]) / size(amps, 2)
    false_pos_2 = count(valid[2, :, b_idx_noise]) / size(amps, 2)
    @test false_pos_1 < 0.10
    @test false_pos_2 < 0.10

    println("All test_basu_fig4_synthetic tests passed successfully!")
    return true
end

function test_find_frequency_peaks()
    println("Running test_find_frequency_peaks...")
    # Generate a frequency grid and a spectrum with two known peaks
    freq = collect(range(0.001, 0.5, length=500))
    # Peak 1 at index 100 (~0.1 cy/P), Peak 2 at index 250 (~0.25 cy/P)
    spec = 0.05 .* randn(length(freq)) .+ 0.2
    spec .+= 1.5 .* exp.(-0.5 .* ((freq .- freq[100]) ./ 0.005).^2)
    spec .+= 0.8 .* exp.(-0.5 .* ((freq .- freq[250]) ./ 0.005).^2)

    pks = Tools.find_frequency_peaks(spec, freq; min_prominence=0.05, min_snr=3.0, max_peaks=5)
    @test length(pks) >= 2
    @test any(abs(pk - 100) <= 2 for pk in pks)
    @test any(abs(pk - 250) <= 2 for pk in pks)
    println("test_find_frequency_peaks passed successfully!")
    return true
end

function test_determine_optimum_nfft()
    println("Running test_determine_optimum_nfft...")
    # Case 1: Strong modulation with many significant spectral points at length 512
    N = 1024
    M = 32
    data_strong = randn(N, M) .* 0.2
    # Add strong broadband or multi-harmonic periodic modulation to bin 16
    for k in 1:6
        data_strong[:, 16] .+= 1.5 .* cos.(2π * (0.05 * k) .* (1:N))
    end
    opt_nfft_1 = Tools.determine_optimum_nfft(data_strong; min_nfft=32, min_sig_fraction=0.15)
    @test opt_nfft_1 >= 64
    @test ispow2(opt_nfft_1)

    # Case 2: Pure noise (no bin has >= 25% significant points at high nfft) -> progressively reduced to min_nfft
    data_noise = randn(N, M)
    opt_nfft_noise = Tools.determine_optimum_nfft(data_noise; min_nfft=32, min_sig_fraction=0.25)
    @test opt_nfft_noise == 32
    @test ispow2(opt_nfft_noise)

    println("test_determine_optimum_nfft passed successfully!")
    return true
end

function test_welch_lrfs_and_diagnostic()
    println("Running test_welch_lrfs_and_diagnostic...")
    N_pulses = 1000
    N_bins = 20
    data = 0.2 .* randn(N_pulses, N_bins)
    # Add sinusoidal modulation at 0.1 cy/P in bin 10
    fp_true = 0.10
    for n in 1:N_pulses
        data[n, 10] += 3.0 * cos(2π * fp_true * n)
    end

    res32 = Tools.welch_lrfs_1d(data, 32; overlap_fraction=0.5)
    res128 = Tools.welch_lrfs_1d(data, 128; overlap_fraction=0.5)

    @test abs(res32.fp - fp_true) < 0.03
    @test abs(res128.fp - fp_true) < 0.015
    # Window with N=128 should have smaller FWHM than N=32
    @test res128.fwhm < res32.fwhm

    all_res = Tools.analyze_fft_lengths(data; N_list=[32, 64, 128])
    @test length(all_res) == 3
    for r in all_res
        @test r.theoretical_fwhm ≈ 1.44 / r.N
    end

    println("test_welch_lrfs_and_diagnostic passed successfully!")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_basu_fig4_synthetic()
    test_find_frequency_peaks()
    test_determine_optimum_nfft()
    test_welch_lrfs_and_diagnostic()
end
