module Solver
export solveEP

using ForceImport
@force using ..Initializer
@force using ..Model

using DifferentialEquations
using Sundials


function solveEP(n_beats, p)
	stim_amp = p[66]
	stim_period = p[68]
	stim_times = collect(0:stim_period:n_beats*stim_period-1)
	anti_stim_times = collect(0:stim_period:n_beats*stim_period-1).+1
	
	function affect1!(integrator)
		integrator.p[81] = stim_amp
	end
	cb1 = PresetTimeCallback(stim_times, affect1!)
	function affect2!(integrator)
		integrator.p[81] = 0.0
	end
	cb2 = PresetTimeCallback(anti_stim_times, affect2!)
	cb = CallbackSet(cb1, cb2)
	
	u0 = initStates()
	tspan = (0.0, n_beats*stim_period)
	prob = ODEProblem(computeRates!, u0, tspan, p)

	alg = CVODE_BDF()
	sol = solve(prob, alg, maxiters=10000000, callback=cb)

	return sol
end


end  # module
