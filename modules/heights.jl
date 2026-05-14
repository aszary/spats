module Heights
    using Statistics
    using LsqFit # For spectral fitting in R2F
    using PhysicalConstants.CODATA2018 # For speed of light (used in calculate_emission_height_rvm)
    #using GaussianFit # To interface with GaussianFit.jl
    using CairoMakie, FileIO # For plotting in geometry_analysis
    #using .RVM # Access RVM module

    include("functions.jl")
    include("tools.jl")
    include("plot.jl")
    include("GaussianFit.jl") # GaussianFit.jl is still needed for Heights.jl
    include ("RVM.jl")

    # Define physical constants
    const C_LIGHT_KMS = PhysicalConstants.CODATA2018.c0 / 1000 # Speed of light in km/s

    export calculate_emission_height_rvm, calculate_beam_opening_angle, calculate_emission_height_r2f, fit_r2f_spectrum
    export fwhm_from_sigma, width_at_10_percent_from_sigma, fwhm_deg_from_sigma_bins, width_at_10_percent_deg_from_sigma_bins

    """
        fwhm_from_sigma(sigma)

    Calculates Full Width at Half Maximum (FWHM) from Gaussian standard deviation (sigma).
    """
    fwhm_from_sigma(sigma::Real) = 2 * sqrt(2 * log(2)) * sigma

    """
        width_at_10_percent_from_sigma(sigma)

    Calculates the width at 10% of the peak amplitude from Gaussian standard deviation (sigma).
    """
    width_at_10_percent_from_sigma(sigma::Real) = 2 * sqrt(2 * log(10)) * sigma

    """
        fwhm_deg_from_sigma_bins(sigma_bins, nbin)

    Calculates Full Width at Half Maximum (FWHM) in degrees from Gaussian standard deviation (sigma) in bins.
    """
    fwhm_deg_from_sigma_bins(sigma_bins::Real, nbin::Int) = fwhm_from_sigma(sigma_bins) / nbin * 360.0

    """
        width_at_10_percent_deg_from_sigma_bins(sigma_bins, nbin)

    Calculates the width at 10% of the peak amplitude in degrees from Gaussian standard deviation (sigma) in bins.
    """
    width_at_10_percent_deg_from_sigma_bins(sigma_bins::Real, nbin::Int) = width_at_10_percent_from_sigma(sigma_bins) / nbin * 360.0

    """
        calculate_emission_height_rvm(phi_PPA, phi_c, P_seconds; filter_polarized_samples=false)

    Calculates the absolute emission height (h_em) based on relativistic aberration and retardation (A/R) effects.

    # Arguments
    - `phi_PPA`: Polarization Position Angle (PPA) inflection point in degrees.
    - `phi_c`: Geometric center of the intensity profile in degrees.
    - `P_seconds`: Pulsar period in seconds.
    - `filter_polarized_samples`: A boolean indicating whether to apply filtering for highly linearly polarized samples.
                                  (Note: Actual filtering logic would depend on input polarization data,
                                  which is not directly handled here. This flag is a placeholder
                                  to indicate the *option` for such a feature.)

    # Returns
    - `h_em`: Emission height in km.
    """
    function calculate_emission_height_rvm(phi_PPA::Real, phi_c::Real, P_seconds::Real;
                                           filter_polarized_samples::Bool=false)
        # The filtering logic itself would depend on having access to polarization data
        # and a PPA traverse. For this function, we assume phi_PPA is already derived
        # potentially with such filtering applied.
        # The `filter_polarized_samples` flag serves as a reminder for future expansion
        # or to indicate that the input phi_PPA should be pre-processed if this is true.

        delta_phi_deg = phi_PPA - phi_c
        h_em = (C_LIGHT_KMS / 4) * (delta_phi_deg / 360.0) * P_seconds
        return h_em
    end

    """
        calculate_beam_opening_angle(W_deg, alpha_deg, beta_deg)

    Calculates the beam opening angle (ρ) using the pulse width and viewing geometry.

    # Arguments
    - `W_deg`: Pulse width in degrees. This can be derived from GaussianFit.jl's sigma.
    - `alpha_deg`: Magnetic inclination angle in degrees.
    - `beta_deg`: Impact parameter in degrees.

    # Returns
    - `rho_deg`: Beam opening angle in degrees.
    """
    function calculate_beam_opening_angle(W_deg::Real, alpha_deg::Real, beta_deg::Real)
        alpha_rad = deg2rad(alpha_deg)
        beta_rad  = deg2rad(beta_deg)
        W_rad     = deg2rad(W_deg)

        sin2_rho_half = sin(alpha_rad) * sin(alpha_rad + beta_rad) * sin(W_rad / 4)^2 + sin(beta_rad / 2)^2
        # Ensure argument to asin is within [-1, 1] due to potential floating point inaccuracies
        sin2_rho_half = max(0.0, min(1.0, sin2_rho_half))
        rho_half_rad = asin(sqrt(sin2_rho_half))
        rho_deg = rad2deg(2 * rho_half_rad)
        return rho_deg
    end

    """
        calculate_emission_height_r2f(rho_nu_deg, P_seconds; s::Real=0.6)

    Calculates the geometric emission height (h_nu) at a specific frequency.

    # Arguments
    - `rho_nu_deg`: Beam opening angle at the given frequency in degrees.
    - `P_seconds`: Pulsar period in seconds.
    - `s`: Mapping parameter (typically 0.6 for conal emission).
           Note: The exact way to incorporate 's' into the formula is not specified in the prompt.
           This function calculates h_nu based on the primary formula. If 's' is meant to
           modify the height, it should be applied as a scaling factor externally or
           a specific formula for its integration needs to be provided.

    # Returns
    - `h_nu`: Emission height in km.
    """
    function calculate_emission_height_r2f(rho_nu_deg::Real, P_seconds::Real; s::Real=0.6)
        # The formula is h_nu = (85.9 / rho_nu)^2 * (2πcP)
        # 85.9 is assumed to be in degrees.
        h_nu = (85.9 / rho_nu_deg)^2 * (2 * pi * C_LIGHT_KMS * P_seconds)
        return h_nu
    end

    """
        fit_r2f_spectrum(frequencies_MHz, heights_km)

    Fits the resulting heights across multiple frequencies to the Thorsett relation:
    h(f) = h_0 + (f/f_h)^b.

    # Arguments
    - `frequencies_MHz`: Vector of frequencies in MHz.
    - `heights_km`: Vector of emission heights in km corresponding to `frequencies_MHz`.

    # Returns
    - A NamedTuple containing the fitted parameters `h_0`, `f_h`, `b`, and the fit result.
    """
    function fit_r2f_spectrum(frequencies_MHz::AbstractVector{<:Real}, heights_km::AbstractVector{<:Real})
        if length(frequencies_MHz) != length(heights_km)
            error("Frequencies and heights vectors must have the same length.")
        end
        if length(frequencies_MHz) < 3
            error("At least 3 data points are required for fitting h(f) = h_0 + (f/f_h)^b.")
        end

        # Model function for h(f) = h_0 + (f/f_h)^b
        # Parameters: p = [h_0, f_h, b]
        model(f, p) = p[1] .+ (f ./ p[2]).^p[3]

        # Initial parameter guess
        h0_p0 = minimum(heights_km) * 0.5
        fh_p0 = mean(frequencies_MHz)
        b_p0  = -2/3

        p0 = [h0_p0, fh_p0, b_p0]

        # Bounds for parameters
        lower_bounds = [0.0, minimum(frequencies_MHz) * 0.1, -5.0]
        upper_bounds = [maximum(heights_km) * 2, maximum(frequencies_MHz) * 10, 0.0]

        fit = curve_fit(model, frequencies_MHz, heights_km, p0; lower=lower_bounds, upper=upper_bounds)

        params = coef(fit)
        h_0_fit = params[1]
        f_h_fit = params[2]
        b_fit   = params[3]

        return (h_0 = h_0_fit, f_h = f_h_fit, b = b_fit, fit_result = fit)
    end

    """
        pos_angle(indir)

    Wrapper function for position angle analysis (alias for geometry_analysis).
    """
    function pos_angle(indir; kwargs...)
        return RVM.geometry_analysis(indir; kwargs...)
    end

end