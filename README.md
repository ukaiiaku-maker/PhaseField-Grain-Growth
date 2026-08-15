# Event-resolved disconnection phase-field grain growth

This repository implements a two-dimensional, constrained multiphase-field
grain-growth model with persistent grain-boundary and triple-junction entities,
discrete disconnection modes, cumulative-hazard event clocks, shear-memory and
nonlocal elastic backends, free-volume/climb kinetics, and entity-attached
pinning controls. It independently reimplements and extends the phase-field and
shear-coupling physics studied in Qiu et al., *Why Grain Growth Is Not Curvature
Flow* (PNAS, DOI 10.1073/pnas.2500707122).

The mathematical model actually executed is in
[`docs/physics_equations.md`](docs/physics_equations.md). Source provenance,
including the externally retained Qiu code, is documented in
[`docs/source_manifest.md`](docs/source_manifest.md) and
[`docs/qiu_code_audit.md`](docs/qiu_code_audit.md). Current scientific gates and
preserved failed results are listed in
[`docs/validation_status.md`](docs/validation_status.md).

## Install and test

```bash
python -m pip install -e '.[analysis,test]'
pytest -q
PYTHONPATH=src python scripts/run_validations.py
PYTHONPATH=src python scripts/validate_stochastic_engine.py
PYTHONPATH=src python scripts/run_qiu_regressions.py
```

## Run and analyze campaigns

```bash
PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/mechanism_scaling_200.yaml --processes 10
PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/temperature_selected_200.yaml --processes 10
PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/temperature_fully_physical_200.yaml --processes 10
PYTHONPATH=src python scripts/analyze_campaign.py results/campaigns/<campaign>
PYTHONPATH=src python scripts/plot_campaign.py results/campaigns/<campaign> \
  --output results/plots/<campaign> \
  --summary results/production_summaries/<campaign>.csv
PYTHONPATH=src python scripts/pareto_campaign.py \
  results/production_summaries/<jerkiness-summary>.csv \
  --output results/production_summaries/<jerkiness-pareto>.csv
PYTHONPATH=src python scripts/aggregate_summaries.py \
  results/production_summaries/<mechanism>.csv \
  results/production_summaries/<temperature-isolation>.csv \
  results/production_summaries/<temperature-physical>.csv \
  results/production_summaries/<jerkiness>.csv \
  --output results/final_mechanism_summary.csv
```

The selected-temperature command is the mechanism-isolation experiment:
intrinsic mobility is temperature-independent except in the B1 control. The
fully physical command lets the selected barriers act together with a 0.45 eV
Arrhenius easy-mode mobility normalized to the validated 900 K baseline.

Completed runs can be extended without modifying their source trajectories:

```bash
PYTHONPATH=src python scripts/extend_campaign.py \
  results/campaigns/<source> --max-steps 3500 \
  --termination-grains 60 --processes 10
```

Every run manifest records the exact launch SHA and full configuration. Dense
checkpoints and raw trajectories remain outside Git under `results/campaigns/`
or `results/runs/`; compact validation reports, summaries, plots, and failure
records are versioned. Scientific analysis rejects incomplete campaigns by
default; `--allow-incomplete` exists only for explicitly labeled diagnostics.
