# CellContraction.jl
A Julia implementation of the human ventricular cell electro-mechanics [Margara et al. 2021](https://doi.org/10.1016/j.pbiomolbio.2020.06.007) model, ported from [Virtual Assay](https://www.cs.ox.ac.uk/ccs/virtual-assay/). The model implements the electro-mechanical coupling of the [ToR-ORd](https://doi.org/10.7554/eLife.48890) electrophysiology model and the [Land](https://doi.org/10.1016/j.yjmcc.2017.03.008) active contraction model.

---
## Information

**Status**: `Completed`

**Type**: `Company project`

**Development year**: `2023-2024`

**Author**: Stefano Longobardi

---
## Getting Started

To get a copy of the project, make sure to have [git](https://git-scm.com) installed on your machine. Then, assuming you have correctly enstablished SSH connection to GitHub, type in your terminal:

```
git clone git@github.com:GSK-Biostatistics/CellContraction.jl.git
```

---
## Prerequisites
To be able to run the project, you need to have [Julia](https://julialang.org/downloads/) (>=1.10.0) installed on your machine.

---
## Installation
Go to the package directory and enter the Julia REPL, then enter the Pkg REPL by pressing ] from the Julia REPL
```
$ cd CellContraction.jl/
$ julia
julia> ]
```
In the Pkg REPL, activate `CellContraction.jl` project environment from the current working directory and instantiate it
```
(@v1.10) pkg> activate .
(CellContraction) pkg> instantiate
```
The last command will 'resolve' the environment to install all the latest versions of the dependencies compatible with the project. A `Manifest.toml` file will be generated as a result of this operation.

Run tests to check that the installed package is working properly
```
(CellContraction) pkg> test
```
This should result in 2 out of 2 total tests passed.

---
## Usage
You can run `CellConctraction.jl` both locally (commonly for control model simulations) and on the HPC (commonly for population simulations).
>Note: by default, all simulations are run for 1000 beats to simulate a steady-state condition. This behaviour can be changed internally.

### Local execution
In the `run/locally/` folder, we have provided 4 main scripts that showcase key functionalities. File names are quite self-explanatory:
1. `run_contr.jl` - Run control model simulations.
2. `run_db_contr.jl` - Run drug block simulations for the control model. This employs a simple pore-block model to calculate the percentage of remaining active current for each ion channel that is blocked by specifying IC50 and Hill coefficient values.
3. `run_mech_contr.jl` - Run mechanism perturbation simulations for the control model. With 'mechanism' we mean any parameter that we have exposed at the user level, allowing to modify its behavior by scaling it with a factor of choice.
4. `run_mech_contr_bonus.jl` - Allows to additionally visualise the effect that blocking ion channels has onto respective current time courses (pure ion channel block).
> Note: all these scripts utilise input data from the `run/locally/data/` folder. To customise your simulation (for example, to add information about a new compound or to run different concentrations or a different mechanism), please adapt these input files using them as templates.

### HPC execution
In the `run/hpc/` folder, we have provided 3 main scripts to perform the exact same operations as described above, but this time in parallel. This is very useful when we need to simulate the full population of models and not just the control model. The scripts can be executed by submitting custom batch scripts to the `Slurm` workload manager. We have provided a utility bash script that automatically generates batch scripts from available templates to be readily submitted to Slurm. The utility is called `generate_slurm_from_template.sh` and can be used as follows:
```
[stefano@server]$ cd run/hpc/
[stefano@server hpc]$ chmod +x generate_slurm_from_template.sh
```
Example 1: baseline population simulation
```
[stefano@server hpc]$ ./generate_slurm_from_template.sh \
    --julia="/home/stefano/.juliaup/bin/julia" \
    --project="/home/stefano/CellContraction.jl" \
    --script="/home/stefano/CellContraction.jl/run/hpc/scripts/run_pop.jl" \
    --nworkers="32" \
    --template="/home/stefano/CellContraction.jl/run/hpc/templates/pop.slurm"
```
Example 2: drug block population simulation (28 compounds)
```
[stefano@server hpc]$ ./generate_slurm_from_template.sh \
    --julia="/home/stefano/.juliaup/bin/julia" \
    --project="/home/stefano/CellContraction.jl" \
    --script="/home/stefano/CellContraction.jl/run/hpc/scripts/run_db_pop.jl" \
    --nworkers="32" \
    --template="/home/stefano/CellContraction.jl/run/hpc/templates/db_mech_pop.slurm" \
    --idx="1-28"
```
Example 3: mechanism perturbation population simulation (9 mechanisms)
```
[stefano@server hpc]$ ./generate_slurm_from_template.sh \
    --julia="/home/stefano/.juliaup/bin/julia" \
    --project="/home/stefano/CellContraction.jl" \
    --script="/home/stefano/CellContraction.jl/run/hpc/scripts/run_mech_pop.jl" \
    --nworkers="32" \
    --template="/home/stefano/CellContraction.jl/run/hpc/templates/db_mech_pop.slurm" \
    --idx="1-9"
```
> Note: both drug block and mechanism perturbation population simulation batch scripts are generated from the same template `db_mech_pop.slurm`. This is because they both make use of job-arrays to run in parallel multiple ion channel block configurations (i.e., compounds) / multiple parameter scaling sets (i.e., mechanisms) (notice the *--idx* flag that controls which compounds / mechanisms to be run according to the row index they occupy in their respective input *\*.csv* files).

Finally, to run the simulations submit the generated batch script to Slurm via the `sbatch` command:
```
[stefano@server hpc]$ sbatch run_pop.slurm
```
Simulation results will be stored in the `run/hpc/output` folder and ready to be used for further analyses.

---
## Reproducing Publication
>Please be aware that your results may still vary slightly from those presented in the paper due to the inherent discrepancies in numerical rounding across different `CellContraction.jl` installations. Nevertheless, general trends and order of magnitudes of calculated quantities should be preserved.

### Negative inotropes: simulations + dose-response analysis
To reproduce the modelling performed on negative inotropic compounds from the paper, generate both the baseline population simulation and the drug block population simulation batch scripts and run them sequentially without modifying any of the input sources from the `run/hpc/data/` folder:
```
[stefano@server hpc]$ sbatch run_pop.slurm
```
await completion, then:
```
[stefano@server hpc]$ sbatch run_db_pop.slurm
```
> Note: you might need to manually adjust the slurm scripts according to the HPC computational resources available. We recommend using at least 32 cores per compound, given that we are simulating 15 concentrations to steady-state for 323 models. With this setup, each compound should run in ~20 minutes.

The output will have the following structure:
```
│run
│├── hpc
│    ├── output
│    │   ├── biomarkers.csv
│    │   ├── Astemizole
│    │   │   ├── biomarkers_0.1.csv
│    │   │   │   ...
│    │   │   ├── biomarkers_100000.0.csv
│    │   ├── Bepridil
│    │   │   ├── biomarkers_0.1.csv
│    │   │   │   ...
│    │   │   ├── biomarkers_100000.0.csv
│    │   │   .
│    │   │   .
│    │   │   .
│    │   ├── Verapamil
│    │   │   ├── biomarkers_0.1.csv
│    │   │   │   ...
│    │   │   ├── biomarkers_100000.0.csv
```

Each *\*.csv* file will contain biomarker values for each model in the population. The `biomarkers.csv` file outside compound-specific folders will be used to normalise drug block-associated values to control values.

#### Dose-response Analysis (Median / Bayesian)
In the `analyse/` folder, we have provided a script to analyse dose-response data and fit a Hill curve through them using either a classic non-linear least squares approach on median biomarker values from the population or a Bayesian approach that takes into account population variability.
1. `run_doseresponse_analysis.jl` - Runs a dose-response data analysis (to use one of the two above-mentioned approaches, you need to edit the script variable `method` to be equal to either `"Median"` or `"Bayesian"`)

This will generate in each compound-specific folder within the overarching output a folder named as the analysis which was run e.g., 'Bayesian' or 'Median', containing fitted IC50 parameter values and plots.

>Note: This can be done for any of the calculated biomarkers (deault: peak tension). For other biomarkers, please inspect the `biomarkers.csv` column headers and choose any other biomarker known to decay with increasing dose.

Additionally, we can compare the derived IC50 values to experimentally measured IC50 values (sarcomere shortening) from [Nguyen et al. 2017](https://doi.org/10.3389/fphys.2017.01073) by running:

2. `run_doseresponse_comparison_with_expdata.jl` - Performs a comparison with experimental IC50 values for the compounds that showed a negative inotropic effect.

Summary plots and csv files will be generated within the `analyse/output/` folder. A joint plot with all the simulated compounds can be obtained and saved in the same folder via:

3. `assemble_doseresponse_plots.jl` - Plots all the dose-response curves for all the simulated compounds (both neutral and negative inotropic ones) with experimental data highlighted in red.

### Positive inotropes: simulations + perturbation-response analysis
To reproduce the modelling performed on positive inotropic compounds from the paper, generate both the baseline population simulation and the mechanism perturbation population simulation batch scripts and run them sequentially without modifying any of the input sources from the `run/hpc/data/` folder:
```
[stefano@server hpc]$ sbatch run_pop.slurm
```
await completion, then:
```
[stefano@server hpc]$ sbatch run_mech_pop.slurm
```

The output will have the following structure:
```
│run
│├── hpc
│    ├── output
│    │   ├── biomarkers.csv
│    │   ├── Ca50_decrease
│    │   │   ├── biomarkers_m1_Ca50_0.1.csv
│    │   │   │   ...
│    │   │   ├── biomarkers_m9_Ca50_0.9.csv
│    │   ├── Cao_increase
│    │   │   ├── biomarkers_m1_Cao_1.0.csv
│    │   │   │   ...
│    │   │   ├── biomarkers_m9_Cao_3.0.csv
│    │   │   .
│    │   │   .
│    │   │   .
│    │   ├── SERCA_increase
│    │   │   ├── biomarkers_m1_GJup_1.0.csv
│    │   │   │   ...
│    │   │   ├── biomarkers_m9_GJup_3.0.csv
```

Each *\*.csv* file will contain biomarker values for each model in the population. The `biomarkers.csv` file outside mechanism-specific folders will be used to normalise mechanism perturbation-associated values to control values.

#### Perturbation-response Analysis (LR deg=2 fit and root finding)
In the `analyse/` folder, we have provided a script to analyse perturbation-response data and fit a parabola through them using non-linear least squares approach on median biomarker values from the population.
1. `run_perturbresponse_analysis.jl` - Runs a perturbation-response data analysis on mechanisms perturbation simulations.

Additionally, we can calculate what is the perturbation needed to achieve a given percent change from control observed experimentally (sarcomere shortening) from [Abi-Gerges et al. 2020](https://doi.org/10.1038/s41598-020-64657-2) by running these 2 scripts sequentially:

2. `run_perturbresponse_preprocess_expdata.jl` - Calculates the percent change from control in sarcomere shortening at the EC50 from digitized experimental dose-response curves.
3. `run_perturbresponse_comparison_with_expdata.jl` - Find the positive root of the perturbation-response parabola, corresponding to the scaling factor that can achieve the same observed tension perturbation from control.

Summary plots and csv files will be generated within the `analyse/output/` folder. A joint plot with all the simulated mechanisms (1-parameter perturbation only) can be obtained and saved in the same folder via:

4. `assemble_perturbresponse_plots.jl` - Plots all the perturbation-response curves for all the mechanisms simulated by altering one model parameter using multiple scaling factors.

For mechanisms simulated via perturbation of 2 parameters simultaneously (e.g., beta agonism), please use the equivalent 2D version `assemble_perturbresponse_plots_2D.jl` (only displays one 2D mechanisms at a time).

### Miscellanea
Other plots are available in the publication supplementary meterials. One series of plots display all the calculated biomarkers' perturbation-response curves at once (18 in total) for pure ion channel block simulations. These plots can be obtained by running `assemble_biomarker_plots.jl` for any given mechanism that is simulated via perturbation of a single parameter. Another plot compares the biomarker values distributions across the full population of models for 3 different implementations of the same cell electro-mechanics model: Virtual Assay, MATLAB, and Julia (this codebase). This shows that the 3 implementations are equivalent (apart from slight variations due to different numerical solvers adopted), and can be obtained by running `assemble_software_comp_plots.jl`. Again, plots will be automatically saved in the `analyse/output/` folder.

---
## License
This project is licensed under the Apache License 2.0. Please refer to the [LICENSE](LICENSE) file for more details.

---
*This README.md complies with [this project template](https://github.com/ShadowTemplate/project-template). Feel free to adopt it and reuse it.*
