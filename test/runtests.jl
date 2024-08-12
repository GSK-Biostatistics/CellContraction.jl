using CellContraction
using Test

# I have hardcoded here the simulated biomarker values from the control model, and taken +/- 5% of these values to account for numerical differences across different package installations.
# If, after installation, your control model simulation produces a Delta(biomarkers) above or below this 5% threshold, maybe something in the installation went wrong.
baseline = abs.(Float64[-88.4257, 32.8727, 193.9, 271.81, 77.91, 345.483, -0.0200801, 71.6222, 462.881, 391.259, 157.99, 308.78, 36.97, 168.09, 22.7545, 113.92, 223.6, 0.247252, -0.140091])
lower_bound = 0.95 .* baseline 
upper_bound = 1.05 .* baseline

@testset "CellContraction" begin
    # check that we can run the control simulation to steady-state (1000 beats) and that both AP, CaT and AT curves fully repolarise
    @test sum(runSimulCompBio(ntuple(x->1.0, Val(14)), 1000)[end-2:end]) == 3
    # check that control simulation biomarker values are within expected range (+/- 5% of baseline)
    @test (abs.(runSimulCompBio(ntuple(x->1.0, Val(14)), 1000)[1:end-3]) > lower_bound) && (abs.(runSimulCompBio(ntuple(x->1.0, Val(14)), 1000)[1:end-3]) < upper_bound)
end
