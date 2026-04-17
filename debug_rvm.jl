include("modules/tools.jl")

# Simple test: does rvm_model match chi2_grid formula?

# Test params
PA0, alpha, zeta, phi0 = 45.0, 30.0, 50.0, -5.0
p = [PA0, alpha, zeta, phi0]

# Test longitudes
lon = [-10.0, -5.0, 0.0, 5.0, 10.0]

# Using rvm_model directly
a  = deg2rad(alpha)
z  = deg2rad(zeta)
dp = deg2rad.(lon .- phi0)
num = sin(a) .* sin.(dp)
den = sin(z) .* cos(a) .- cos(z) .* sin(a) .* cos.(dp)
pa_manual = PA0 .+ rad2deg.(atan.(num, den))

# Using rvm_model function
pa_rvm = Tools.rvm_model(lon, p)

println("Manual RVM calculation:")
println(pa_manual)
println("\nUsing rvm_model function:")
println(pa_rvm)
println("\nDifference:")
println(pa_manual .- pa_rvm)

# Wrap both to [-90, 90]
pa_manual_wrapped = mod.(pa_manual .+ 90.0, 180.0) .- 90.0
pa_rvm_wrapped = mod.(pa_rvm .+ 90.0, 180.0) .- 90.0

println("\n=== WRAPPED ===")
println("Manual wrapped:")
println(pa_manual_wrapped)
println("\nrvm_model wrapped:")
println(pa_rvm_wrapped)
println("\nWrapped difference:")
println(pa_manual_wrapped .- pa_rvm_wrapped)
