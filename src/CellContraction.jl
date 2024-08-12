module CellContraction

using ForceImport, Reexport

include("Currents/Currents.jl")        # Functions implementing different ionic currents
include("Initializer/Initializer.jl")  # ODE initial guess + parameters
include("Model/Model.jl")              # The actual EP + Contr model
include("Solver/Solver.jl")            # ODE solver
include("Utils/Utils.jl")              # Utility functions
include("Runner/Runner.jl")            # Powerful functions to run the model and return solution + different properties
include("Drug/Drug.jl")                # Drug module implementing the pore-block model
include("DoseResponse/Bayesian.jl")    # Dose-response analysis using Bayesian modelling
include("DoseResponse/Median.jl")      # Dose-response analysis using fit to median values
include("PerturbResponse/Mech.jl")     # Perturbation-response analysis for mechanisms

# "Merging" of the modules
@force    using .Currents
@reexport using .Currents
@force    using .Initializer
@reexport using .Initializer
@force    using .Model
@reexport using .Model
@force    using .Utils
@reexport using .Utils
@force    using .Solver
@reexport using .Solver
@force    using .Runner
@reexport using .Runner
@force    using .Drug
@reexport using .Drug
@force    using .Bayesian
@reexport using .Bayesian
@force    using .Median
@reexport using .Median
@force    using .Mech
@reexport using .Mech

end  # module
