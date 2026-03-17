# fractel-simulations

This repository contains the R-based simulation framework and analysis scripts used to evaluate the FRACTEL method.

> Note: This repository provides simulation and analysis code for the FRACTEL study. It is not the primary implementation of the FRACTEL test.

---

## Repository Structure

```

fractel-simulations/
├── config/        # Configuration files for simulation parameters
├── data/          # Input data (raw and processed)
├── scripts/       # Core simulation and analysis scripts
├── results/       # Simulation outputs and derived results
├── docs/          # Additional documentation

```

---

## Workflow Overview

The simulation pipeline proceeds in the following stages:

1. **Prepare input parameters**
```

Rscript scripts/01_prepare_input.R

```

2. **Compile simulation parameter grid**
```

Rscript scripts/02_compile_parameters.R

```

3. **Generate null distributions**
```

Rscript scripts/03_simulate_null.R

```

4. **Run alternative simulations**
```

Rscript scripts/04_simulate_alternative.R <task_id>

```

5. **Aggregate results and compute power**
```

Rscript scripts/05_aggregate_results.R

```

---

## Parallel Execution (Slurm)

Large-scale simulations were executed using Slurm array jobs.

The alternative simulation script is designed to run per-task:

```

Rscript scripts/04_simulate_alternative.R <task_id>

```

where `<task_id>` indexes rows of the parameter grid (`output/input_parameters.txt`).

Example Slurm usage:

```

sbatch --array=1-N run_alt_sim.sh

```

Users without access to Slurm can adapt this interface using local loops or job schedulers.

---

## Configuration

Simulation parameters and input file paths are defined in:

```

config/config.yaml

```

This includes:
- input data locations
- simulation sizes
- parameter ranges
- normalization and interpolation settings

---

## Requirements

- R (version 4.x recommended)
- Required packages:
  - yaml
  - MASS
  - dplyr
  - stats

---

## Reproducibility Notes

- All scripts are intended to be run from the repository root.
- Output files are written to the `results/` directory.
- Some simulations may be computationally intensive and require parallel execution.

---

## Data

Input data required for reproducing simulations should be placed in:

```

data/raw/

```

Processed intermediate files are written to:

```

data/processed/

```

(Details of expected formats are provided in `docs/workflow.md`.)

---

## License

This project is licensed under the MIT License.

---

## Citation

If you use this code, please cite the associated FRACTEL manuscript (preprint forthcoming).

