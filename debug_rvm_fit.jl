include("modules/data.jl")
include("modules/tools.jl")

# Simulate test data
n_pulses = 100
n_bins = 512

# Create simple RVM-like PA data
alpha_true = 30.0
zeta_true = 50.0
phi0_true = -5.0
PA0_true = 45.0

lon_test = collect(range(-20.0, 20.0, length=50))

# Generate synthetic PA data with some noise
pa_model = Tools.rvm_model(lon_test, [PA0_true, alpha_true, zeta_true, phi0_true])
pa_model_wrapped = mod.(pa_model .+ 90.0, 180.0) .- 90.0

pa_data = pa_model_wrapped .+ randn(length(lon_test)) .* 3.0  # add noise
pa_err = fill(3.0, length(lon_test))
mask = trues(length(lon_test))

println("True params: PA0=$PA0_true, α=$alpha_true, ζ=$zeta_true, φ0=$phi0_true")
println("PA data: ", round.(pa_data, digits=1))
println()

# Test chi2_grid
alphas, betas, chi2_map, best_p = Tools.chi2_grid(lon_test, pa_data, pa_err, mask)
println("Chi2 grid found: PA0=$(round(best_p[1], digits=1)), α=$(round(best_p[2], digits=1)), ζ=$(round(best_p[3], digits=1)), φ0=$(round(best_p[4], digits=1))")

# Test fit_rvm (should refine it)
result = Tools.fit_rvm(lon_test, pa_data, pa_err, mask; p0=best_p)
println("After fit_rvm: PA0=$(round(result.PA0, digits=1)), α=$(round(result.alpha, digits=1)), ζ=$(round(result.zeta, digits=1)), φ0=$(round(result.phi0, digits=1))")
println("Chi2_red=$(round(result.chi2_red, digits=2)), rms=$(round(result.rms_deg, digits=2))°")

# Check fitted curve
pa_fitted = Tools.rvm_model(lon_test, [result.PA0, result.alpha, result.zeta, result.phi0])
pa_fitted_wrapped = mod.(pa_fitted .+ 90.0, 180.0) .- 90.0

println("\nFitted PA:  ", round.(pa_fitted_wrapped, digits=1))
println("Observed PA:", round.(pa_data, digits=1))
println("Residuals:  ", round.(pa_fitted_wrapped .- pa_data, digits=1))
