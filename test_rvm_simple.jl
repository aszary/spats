# Quick test: synthetic RVM data fitting

# Simple synthetic RVM
function rvm_model_test(phi, p)
    PA0, alpha, zeta, phi0 = p
    a  = deg2rad(alpha)
    z  = deg2rad(zeta)
    dp = deg2rad.(phi .- phi0)
    num = sin(a) .* sin.(dp)
    den = sin(z) .* cos(a) .- cos(z) .* sin(a) .* cos.(dp)
    return PA0 .+ rad2deg.(atan.(num, den))
end

# Test data
lon_test = collect(range(-20.0, 20.0, length=50))
p_true = [45.0, 30.0, 50.0, -5.0]

pa_true = rvm_model_test(lon_test, p_true)
pa_wrapped = mod.(pa_true .+ 90.0, 180.0) .- 90.0

# Add noise
pa_obs = pa_wrapped .+ randn(length(lon_test)) .* 2.0
pa_err = fill(2.0, length(lon_test))

println("True params: PA0=$(p_true[1]), α=$(p_true[2]), ζ=$(p_true[3]), φ0=$(p_true[4])")
println("PA observed (first 10): ", round.(pa_obs[1:10], digits=1))
println("PA true (first 10):    ", round.(pa_wrapped[1:10], digits=1))

# Test unwrapping discontinuity handling
pa_jump = [80.0, -85.0]
pa_jump_wrapped = mod.(pa_jump .+ 90.0, 180.0) .- 90.0
println("\nWrap test: $pa_jump -> $pa_jump_wrapped")
println("Diff: ", pa_jump_wrapped[2] - pa_jump_wrapped[1], " ° (should be ~-175 wrapping smoothly)")
