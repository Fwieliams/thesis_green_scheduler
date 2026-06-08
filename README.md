# Thesis Scheduler Simulator

This folder contains the working simulator for the thesis project.

The simulator does **not** run the real ChIP-seq workflow. It simulates scheduling
decisions for a reduced 15-task ChIP-seq DAG on a modeled heterogeneous cluster.
The trace files and carbon-intensity file needed for the trace-driven run are
included in this repository under `data/`, so the simulator can be reproduced
after cloning this repository only.

## What it compares

- `baseline`: performance-first scheduler that chooses the node with the earliest
  predicted finish time.
- `green_aware`: predictive green-aware scheduler that applies a performance
  guardrail and then chooses based on a weighted score for time, energy and carbon.

## Run

From this folder:

```bash
python3 simulator.py
```

To also run guardrail/weight sensitivity analysis:

```bash
python3 simulator.py --sensitivity
```

This reproduces the trace-driven CSV/text outputs in `results_trace_driven/`.
The simulator uses:

- `data/traces/workflows/chipseq/chipseq-1.csv`
- `data/traces/workflows/chipseq/chipseq-2.csv`
- `data/traces/workflows/chipseq/chipseq-3.csv`
- `data/intensity/out/de-15112023-08122023.csv`

To reproduce the older placeholder-only version:

```bash
python3 simulator.py --input-source placeholder --sensitivity
```

## Outputs

The script writes:

- `results_trace_driven/schedule_baseline.csv`
- `results_trace_driven/schedule_green_aware.csv`
- `results_trace_driven/summary_primary.csv`
- `results_trace_driven/comparison_primary.txt`
- `results_trace_driven/runtime_profile.csv`
- `results_trace_driven/input_profile.txt`
- `results_trace_driven/summary_sensitivity.csv` when `--sensitivity` is used

The older `results` folder contains preliminary placeholder results from the
first simulator version. Keep it as a record, but do not use it as the final
trace-driven result set.

## Reproduce from a fresh clone

```bash
git clone https://github.com/Fwieliams/thesis_green_scheduler.git
cd thesis_green_scheduler
python3 simulator.py --sensitivity
```

No external workflow repository is required for the default trace-driven run.

## Important thesis note

The simulator uses trace-derived predictions, not physical measurements on the
modeled Sherwood/Olympus/Atlantis/Camelot nodes. In the thesis, describe this as
a trace-driven simulation study using a reduced ChIP-seq DAG.
