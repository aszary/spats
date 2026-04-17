# Quick sanity check: PA values from RVM model

PA0, alpha, zeta, phi0 = 45.0, 30.0, 50.0, -5.0

lon = [-10.0, -5.0, 0.0, 5.0, 10.0]

# Manual calculation (same as rvm_model)
a  = deg2rad(alpha)
z  = deg2rad(zeta)
dp_rad = deg2rad.(lon .- phi0)

num = sin(a) .* sin.(dp_rad)
den = sin(z) .* cos(a) .- cos(z) .* sin(a) .* cos.(dp_rad)

pa_unwrapped = PA0 .+ rad2deg.(atan.(num, den))
pa_wrapped = mod.(pa_unwrapped .+ 90.0, 180.0) .- 90.0

println("PA unwrapped: $pa_unwrapped")
println("PA wrapped: $pa_wrapped")

# Expected smooth curve, no kinks
if any(abs.(diff(pa_wrapped)) .> 30)
    println("WARNING: Large jumps in PA! Wrapping issue?")
else
    println("OK: PA curve is smooth")
end
