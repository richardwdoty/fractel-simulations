# fractel-simulations

This repository contains the R-based simulation framework and analysis scripts used to evaluate the FRACTEL method.

> Note: This repository provides simulation and analysis code for the FRACTEL study. It is not the primary implementation of the FRACTEL test.

## Repository Structure

```text
fractel-simulations/
├── config/                 # Configuration files and optional gene list
├── data/
│   ├── raw/                # Raw input files
│   └── processed/          # Derived intermediate files used by the pipeline
├── scripts/                # Core simulation and aggregation scripts
├── results/                # Final simulation outputs
└── docs/                   # Additional documentation
````

## Pipeline Overview

The simulation workflow proceeds in five stages:

1. **Prepare input parameters**

   ```bash
   Rscript scripts/01_prepare_input.R
   ```

2. **Build parameter grid**

   ```bash
   Rscript scripts/02_build_parameter_grid.R
   ```

3. **Simulate null reference distributions**

   ```bash
   Rscript scripts/03_simulate_null.R
   ```

4. **Run alternative simulations**

   ```bash
   Rscript scripts/04_simulate_alternative.R <task_id>
   ```

5. **Aggregate task-level results**

   ```bash
   Rscript scripts/05_aggregate_results.R
   ```

## Parallel Execution

Alternative simulations are designed to run as independent tasks, one parameter row per task. This makes the workflow suitable for Slurm array jobs or other parallel execution systems.

Example Slurm usage:

```bash
sbatch --array=1-N run_alt_sim.sh
```

where `N` is the number of rows in `data/processed/input_parameters.tsv`.

Each task writes a separate output file to:

```text
results/alternative/by_task/
```

These task-level outputs are then combined by `05_aggregate_results.R`.

## Configuration

The pipeline is controlled by:

```text
config/config.yaml
```

Key configuration sections include:

* `input_files`: raw beta0, dispersion, and counts input files
* `gene_list_file`: optional gene subset file used only when interpolation is disabled
* `interpolation`: settings for generating interpolated `(beta0, phi)` pairs
* `parameters`: simulation design parameters such as `n`, `q`, effect size, and `k_value`
* `simulation`: numbers of null and alternative simulations
* `alt_simulation`: choice of approximation or count-based simulation mode

## Input Files

Place the following raw input files in:

```text
data/raw/
```

Required files:

* beta0 table
* dispersion table
* counts table

The expected filenames are specified in `config/config.yaml`.

If using the non-interpolated observed-value workflow with gene subsetting, place the gene list file in:

```text
config/
```

### Negative Binomial Parameterization

Dispersion inputs are assumed to follow the parameterization:

$\text{E}[Y] = \mu$
$\text{Var}[Y] = \mu + \phi \mu^2$

i.e., the NB distribution is parameterized by mean $\mu$ and dispersion $\phi$.

Dispersion inputs may be provided as either $\phi$ or $\theta$ (size).
If $\theta$ is provided, it is converted internally via $\phi$ = 1 / $\theta$.

## Main Outputs

Final combined outputs are written to:

```text
results/alternative/
```

including:

* `power_by_k.tsv`: full power curve for each task across all values of `k = 1, ..., n`
* `power_at_configured_k.tsv`: one-row-per-task summary using the configured bound

Intermediate files written during preprocessing and null simulation are stored under `data/processed/`.

## Modeling Notes

The mean effect size is specified through `exp_beta1`, the fold change corresponding to `exp(beta1)`.

Guide-level effects are simulated on the log scale:

```text
beta1 ~ Normal(log(exp_beta1), beta1_sd^2)
```

so when `beta1_sd > 0`, the induced fold-change distribution is log-normal.

Interpolation is the default workflow. When interpolation is disabled, the pipeline can instead use exact fitted `(beta0, phi)` values, optionally restricted to a supplied gene list.

The bound parameter is configured through `k_value`:

* if `0 < k_value < 1`, then `k = ceiling(k_value * n)`
* if `k_value >= 1`, then `k = as.integer(k_value)`

## Requirements

* R (4.x recommended)
* Required packages:

  * `yaml`
  * `MASS`
  * `dplyr`
  * `stats`

## Reproducibility Notes

* All scripts are intended to be run from the repository root.
* Large simulation outputs and intermediate files are not intended to be tracked in Git.
* The count-based alternative simulation mode is substantially slower than the normal-approximation mode and is best treated as a validation path rather than the default workflow.

## License

This project is licensed under the MIT License.

## Citation

If you use this code, please cite the associated FRACTEL manuscript or preprint when available.