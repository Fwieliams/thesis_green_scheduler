"""
Fardeen Wieliams - 13770160
Afstudeerproject Informatiekunde / Bachelor Thesis

Trace-driven scheduler simulator for the thesis project.

This script does not execute the real ChIP-seq workflow. Instead, it simulates
how a reduced DAG would be scheduled on a modeled heterogeneous cluster.

It compares:
1. A performance-first baseline scheduler.
2. A predictive green-aware scheduler with a performance guardrail.

"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable


# ---------------------------------------------------------------------------
# Configurable experiment assumptions
# ---------------------------------------------------------------------------


PRIMARY_GUARDRAIL = 0.10
PRIMARY_WEIGHTS = {
    "time": 0.40,
    "energy": 0.20,
    "carbon": 0.40,
}

SENSITIVITY_GUARDRAILS = [0.05, 0.10, 0.15]
SENSITIVITY_WEIGHT_SETS = [
    {"time": 0.60, "energy": 0.20, "carbon": 0.20},  # performance-focused
    {"time": 0.40, "energy": 0.20, "carbon": 0.40},  # primary setting
    {"time": 0.40, "energy": 0.40, "carbon": 0.20},  # energy-focused
    {"time": 0.20, "energy": 0.40, "carbon": 0.40},  # sustainability-focused
]

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
CHIPSEQ_TRACE_DIR = DATA_DIR / "traces" / "workflows" / "chipseq"
CHIPSEQ_TRACE_FILES = [
    CHIPSEQ_TRACE_DIR / "chipseq-1.csv",
    CHIPSEQ_TRACE_DIR / "chipseq-2.csv",
    CHIPSEQ_TRACE_DIR / "chipseq-3.csv",
]
CARBON_TRACE_FILE = (
    DATA_DIR
    / "intensity"
    / "out"
    / "de-15112023-08122023.csv"
)


# Hourly average carbon intensity in gCO2e/kWh.
# This is a deterministic placeholder profile for the first working simulator.
# Replace with region-specific data later if available.
PLACEHOLDER_CARBON_INTENSITY_G_PER_KWH = [
    420, 410, 405, 395, 385, 370,
    350, 330, 300, 270, 245, 230,
    220, 225, 240, 265, 300, 340,
    380, 410, 435, 450, 445, 430,
]

CARBON_INTENSITY_G_PER_KWH = list(PLACEHOLDER_CARBON_INTENSITY_G_PER_KWH)


# Reduced logical ChIP-seq DAG.
# base_runtime_min is the reference duration on a node with speed coefficient 1.0.
# These are modeling assumptions for the first simulator version and should be
# replaced by trace-derived medians when those values are available.
PLACEHOLDER_TASKS = {
    "RawQC": {
        "predecessors": [],
        "base_runtime_min": 12.0,
        "memory_gb": 8,
        "utilization": 0.35,
    },
    "Trim": {
        "predecessors": ["RawQC"],
        "base_runtime_min": 18.0,
        "memory_gb": 16,
        "utilization": 0.35,
    },
    "Align_1": {
        "predecessors": ["Trim"],
        "base_runtime_min": 95.0,
        "memory_gb": 32,
        "utilization": 0.75,
    },
    "Align_2": {
        "predecessors": ["Trim"],
        "base_runtime_min": 102.0,
        "memory_gb": 32,
        "utilization": 0.75,
    },
    "Align_3": {
        "predecessors": ["Trim"],
        "base_runtime_min": 88.0,
        "memory_gb": 32,
        "utilization": 0.75,
    },
    "Align_4": {
        "predecessors": ["Trim"],
        "base_runtime_min": 110.0,
        "memory_gb": 32,
        "utilization": 0.75,
    },
    "Filter_1": {
        "predecessors": ["Align_1"],
        "base_runtime_min": 24.0,
        "memory_gb": 16,
        "utilization": 0.55,
    },
    "Filter_2": {
        "predecessors": ["Align_2"],
        "base_runtime_min": 26.0,
        "memory_gb": 16,
        "utilization": 0.55,
    },
    "Filter_3": {
        "predecessors": ["Align_3"],
        "base_runtime_min": 22.0,
        "memory_gb": 16,
        "utilization": 0.55,
    },
    "Filter_4": {
        "predecessors": ["Align_4"],
        "base_runtime_min": 28.0,
        "memory_gb": 16,
        "utilization": 0.55,
    },
    "Merge": {
        "predecessors": ["Filter_1", "Filter_2", "Filter_3", "Filter_4"],
        "base_runtime_min": 35.0,
        "memory_gb": 64,
        "utilization": 0.55,
    },
    "MarkDuplicates": {
        "predecessors": ["Merge"],
        "base_runtime_min": 50.0,
        "memory_gb": 64,
        "utilization": 0.75,
    },
    "PeakCalling": {
        "predecessors": ["MarkDuplicates"],
        "base_runtime_min": 42.0,
        "memory_gb": 32,
        "utilization": 0.75,
    },
    "Annotation": {
        "predecessors": ["PeakCalling"],
        "base_runtime_min": 20.0,
        "memory_gb": 16,
        "utilization": 0.55,
    },
    "Report": {
        "predecessors": ["Annotation"],
        "base_runtime_min": 10.0,
        "memory_gb": 8,
        "utilization": 0.35,
    },
}


TASKS = {
    task: {**spec, "predecessors": list(spec["predecessors"])}
    for task, spec in PLACEHOLDER_TASKS.items()
}


NODE_TYPES = {
    "Sherwood": {
        "role": "low-power",
        "speed": 1.00,
        "memory_gb": 32,
        "idle_power_w": 25,
        "max_power_w": 55,
    },
    "Olympus": {
        "role": "legacy-balanced",
        "speed": 1.30,
        "memory_gb": 64,
        "idle_power_w": 55,
        "max_power_w": 110,
    },
    "Atlantis": {
        "role": "compute-oriented",
        "speed": 1.80,
        "memory_gb": 128,
        "idle_power_w": 85,
        "max_power_w": 170,
    },
    "Camelot": {
        "role": "memory-rich",
        "speed": 1.60,
        "memory_gb": 256,
        "idle_power_w": 95,
        "max_power_w": 160,
    },
}

NODES_PER_TYPE = 2


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Node:
    node_id: str
    node_type: str
    role: str
    speed: float
    memory_gb: int
    idle_power_w: float
    max_power_w: float


@dataclass(frozen=True)
class Candidate:
    task: str
    node: Node
    start_min: float
    finish_min: float
    runtime_min: float
    power_w: float
    energy_kwh: float
    carbon_g: float


@dataclass(frozen=True)
class ScheduleEntry:
    scheduler: str
    task: str
    node_id: str
    node_type: str
    node_role: str
    start_min: float
    finish_min: float
    runtime_min: float
    energy_kwh: float
    carbon_g: float
    memory_required_gb: int
    node_memory_gb: int


@dataclass(frozen=True)
class Summary:
    scheduler: str
    makespan_min: float
    total_energy_kwh: float
    total_carbon_g: float
    average_utilization: float
    guardrail: float | None = None
    weight_time: float | None = None
    weight_energy: float | None = None
    weight_carbon: float | None = None


@dataclass(frozen=True)
class InputProfile:
    runtime_source: str
    carbon_source: str
    trace_files: tuple[Path, ...] = ()
    carbon_file: Path | None = None
    stage_stats: dict[str, dict[str, float]] | None = None
    carbon_points: int = 0


# ---------------------------------------------------------------------------
# Trace-derived input preparation
# ---------------------------------------------------------------------------


def copy_task_table(source: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Create a mutable copy of the task table and predecessor lists."""
    return {
        task: {**spec, "predecessors": list(spec["predecessors"])}
        for task, spec in source.items()
    }


def reset_to_placeholder_inputs() -> InputProfile:
    """Reset the simulator to the built-in placeholder task and carbon inputs."""
    global TASKS, CARBON_INTENSITY_G_PER_KWH
    TASKS = copy_task_table(PLACEHOLDER_TASKS)
    CARBON_INTENSITY_G_PER_KWH = list(PLACEHOLDER_CARBON_INTENSITY_G_PER_KWH)
    return InputProfile(
        runtime_source="placeholder",
        carbon_source="placeholder",
        carbon_points=len(CARBON_INTENSITY_G_PER_KWH),
    )


def use_trace_inputs() -> InputProfile:
    """Load trace-derived task parameters and the external carbon profile."""
    global TASKS, CARBON_INTENSITY_G_PER_KWH
    tasks, stage_stats = load_trace_task_profile(CHIPSEQ_TRACE_FILES)
    carbon_profile = load_carbon_profile(CARBON_TRACE_FILE)
    TASKS = tasks
    CARBON_INTENSITY_G_PER_KWH = carbon_profile
    return InputProfile(
        runtime_source="trace-derived ChIP-seq medians/quantiles",
        carbon_source="German hourly carbon-intensity trace",
        trace_files=tuple(CHIPSEQ_TRACE_FILES),
        carbon_file=CARBON_TRACE_FILE,
        stage_stats=stage_stats,
        carbon_points=len(carbon_profile),
    )


def require_files(paths: list[Path]) -> None:
    """Raise an error if one or more required input files are missing."""
    missing = [path for path in paths if not path.exists()]
    if missing:
        joined = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(
            "Required input file(s) not found. The repository should include "
            "the trace files under data/traces/workflows/chipseq and the "
            "carbon-intensity file under data/intensity/out.\n"
            f"Missing:\n{joined}"
        )


def load_carbon_profile(path: Path) -> list[float]:
    """Read hourly carbon-intensity values from a CSV file."""
    require_files([path])
    values: list[float] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            raw_value = row.get("actual")
            if raw_value in (None, ""):
                continue
            values.append(float(raw_value))

    if not values:
        raise ValueError(f"No carbon-intensity values found in {path}")
    return values


def map_process_to_stage(process: str) -> str | None:
    """Map a trace-level process name to one reduced logical DAG stage."""
    if "FASTQC_TRIMGALORE:FASTQC" in process:
        return "RawQC"
    if "FASTQC_TRIMGALORE:TRIMGALORE" in process:
        return "Trim"
    if "ALIGN_BWA_MEM" in process:
        return "Align"
    if "FILTER_BAM_BAMTOOLS" in process or "BEDTOOLS_GENOMECOV" in process:
        return "Filter"
    if "PICARD_MERGESAMFILES" in process:
        return "Merge"
    if "MARK_DUPLICATES_PICARD" in process:
        return "MarkDuplicates"
    if "MACS2_CALLPEAK" in process or "MACS2_CONSENSUS" in process or "FRIP_SCORE" in process:
        return "PeakCalling"
    if (
        "HOMER_ANNOTATEPEAKS" in process
        or "ANNOTATE_BOOLEAN_PEAKS" in process
        or "SUBREAD_FEATURECOUNTS" in process
        or "DESEQ2_QC" in process
    ):
        return "Annotation"
    if (
        "MULTIQC" in process
        or "DEEPTOOLS" in process
        or process.endswith(":IGV")
        or "PLOT_" in process
        or "CUSTOM_DUMPSOFTWAREVERSIONS" in process
        or "PHANTOMPEAKQUALTOOLS" in process
        or "PICARD_COLLECTMULTIPLEMETRICS" in process
        or "UCSC_BEDGRAPHTOBIGWIG" in process
        or "PRESEQ_LCEXTRAP" in process
    ):
        return "Report"
    return None


def parse_runtime_minutes(row: dict[str, str]) -> float:
    """Extract a positive task runtime from a trace row and convert it to minutes."""
    for field in ("realtime", "duration"):
        raw_value = row.get(field)
        if raw_value not in (None, "", "-"):
            value = float(raw_value)
            if value > 0:
                return value / 60000.0
    return 0.0


def percentile(values: list[float], fraction: float) -> float:
    """Compute a percentile using linear interpolation over sorted values."""
    if not values:
        raise ValueError("Cannot compute percentile of an empty list")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    lower_value = ordered[lower]
    upper_value = ordered[upper]
    return lower_value + (upper_value - lower_value) * (position - lower)


def median_or_default(values: list[float], default: float) -> float:
    """Return the median of a list, or a default when the list is empty."""
    return percentile(values, 0.50) if values else default


def clamp(value: float, minimum: float, maximum: float) -> float:
    """Return the median of a list, or a default when the list is empty."""
    return max(minimum, min(maximum, value))


def load_trace_task_profile(
    trace_files: list[Path],
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, float]]]:
    """Build reduced-DAG task parameters from raw ChIP-seq trace files."""
    require_files(trace_files)

    stage_runtime_by_unit: dict[str, dict[str, float]] = {}
    stage_memory_values: dict[str, list[float]] = {}
    stage_utilization_values: dict[str, list[float]] = {}
    mapped_rows = 0
    skipped_rows = 0

    for trace_file in trace_files:
        with trace_file.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                stage = map_process_to_stage(row["process"])
                if stage is None:
                    skipped_rows += 1
                    continue

                runtime_min = parse_runtime_minutes(row)
                if runtime_min <= 0:
                    skipped_rows += 1
                    continue

                tag = row.get("tag") or row.get("task_id") or "unknown"
                unit_key = f"{trace_file.stem}:{tag}"
                stage_runtime_by_unit.setdefault(stage, {})
                stage_runtime_by_unit[stage][unit_key] = (
                    stage_runtime_by_unit[stage].get(unit_key, 0.0) + runtime_min
                )

                memory_bytes = float(row["memory"])
                stage_memory_values.setdefault(stage, []).append(memory_bytes / (1024**3))

                cpus = max(float(row["cpus"]), 1.0)
                cpu_percent = float(row["%cpu"])
                stage_utilization_values.setdefault(stage, []).append(
                    clamp(cpu_percent / (100.0 * cpus), 0.05, 1.0)
                )
                mapped_rows += 1

    stage_values = {
        stage: list(unit_runtimes.values())
        for stage, unit_runtimes in stage_runtime_by_unit.items()
    }
    required_stages = {
        "RawQC",
        "Trim",
        "Align",
        "Filter",
        "Merge",
        "MarkDuplicates",
        "PeakCalling",
        "Annotation",
        "Report",
    }
    missing_stages = sorted(stage for stage in required_stages if not stage_values.get(stage))
    if missing_stages:
        raise ValueError(f"No mapped trace values found for stages: {', '.join(missing_stages)}")

    stage_stats: dict[str, dict[str, float]] = {}
    for stage, values in stage_values.items():
        stage_stats[stage] = {
            "task_units": float(len(values)),
            "runtime_p25_min": percentile(values, 0.25),
            "runtime_p50_min": percentile(values, 0.50),
            "runtime_p75_min": percentile(values, 0.75),
            "runtime_p90_min": percentile(values, 0.90),
            "memory_p95_gb": percentile(stage_memory_values.get(stage, [1.0]), 0.95),
            "utilization_p50": median_or_default(stage_utilization_values.get(stage, []), 0.50),
        }
    stage_stats["_trace_rows"] = {
        "mapped_rows": float(mapped_rows),
        "skipped_rows": float(skipped_rows),
        "trace_files": float(len(trace_files)),
    }

    tasks = copy_task_table(PLACEHOLDER_TASKS)

    def apply_stage(
        task: str,
        stage: str,
        runtime_fraction: float,
        fallback_utilization: float,
    ) -> None:
        """Apply trace-derived runtime, memory, and utilization values to one task."""
        stats = stage_stats[stage]
        tasks[task]["base_runtime_min"] = max(0.05, percentile(stage_values[stage], runtime_fraction))
        tasks[task]["memory_gb"] = max(1, math.ceil(stats["memory_p95_gb"]))
        tasks[task]["utilization"] = clamp(stats["utilization_p50"], 0.10, 0.95) or fallback_utilization

    apply_stage("RawQC", "RawQC", 0.50, 0.35)
    apply_stage("Trim", "Trim", 0.50, 0.35)
    for task, fraction in zip(("Align_1", "Align_2", "Align_3", "Align_4"), (0.25, 0.50, 0.75, 0.90)):
        apply_stage(task, "Align", fraction, 0.75)
    for task, fraction in zip(("Filter_1", "Filter_2", "Filter_3", "Filter_4"), (0.25, 0.50, 0.75, 0.90)):
        apply_stage(task, "Filter", fraction, 0.55)
    apply_stage("Merge", "Merge", 0.50, 0.55)
    apply_stage("MarkDuplicates", "MarkDuplicates", 0.50, 0.75)
    apply_stage("PeakCalling", "PeakCalling", 0.50, 0.75)
    apply_stage("Annotation", "Annotation", 0.50, 0.55)
    apply_stage("Report", "Report", 0.50, 0.35)

    return tasks, stage_stats


# ---------------------------------------------------------------------------
# Model helpers
# ---------------------------------------------------------------------------


def build_nodes() -> list[Node]:
    """Create the modeled heterogeneous cluster from node-type assumptions."""
    nodes: list[Node] = []
    for node_type, spec in NODE_TYPES.items():
        for index in range(1, NODES_PER_TYPE + 1):
            nodes.append(
                Node(
                    node_id=f"{node_type}_{index}",
                    node_type=node_type,
                    role=str(spec["role"]),
                    speed=float(spec["speed"]),
                    memory_gb=int(spec["memory_gb"]),
                    idle_power_w=float(spec["idle_power_w"]),
                    max_power_w=float(spec["max_power_w"]),
                )
            )
    return nodes


def active_power_w(node: Node, utilization: float) -> float:
    """Estimate active node power with a linear idle-to-maximum power model."""
    utilization = max(0.0, min(1.0, utilization))
    return node.idle_power_w + utilization * (node.max_power_w - node.idle_power_w)


def runtime_on_node_min(task: str, node: Node) -> float:
    """Estimate task runtime on a node using the node speed coefficient."""
    return float(TASKS[task]["base_runtime_min"]) / node.speed


def interval_carbon_intensity(start_min: float, finish_min: float) -> float:
    """Return duration-weighted average carbon intensity over an interval."""
    if finish_min <= start_min:
        raise ValueError("finish_min must be greater than start_min")

    total_weighted_intensity = 0.0
    remaining_start = start_min
    while remaining_start < finish_min:
        hour_index = int(remaining_start // 60) % len(CARBON_INTENSITY_G_PER_KWH)
        next_hour_boundary = (int(remaining_start // 60) + 1) * 60
        segment_end = min(finish_min, next_hour_boundary)
        duration = segment_end - remaining_start
        total_weighted_intensity += CARBON_INTENSITY_G_PER_KWH[hour_index] * duration
        remaining_start = segment_end

    return total_weighted_intensity / (finish_min - start_min)


def estimate_candidate(
    task: str,
    node: Node,
    node_available_at: dict[str, float],
    finished_at: dict[str, float],
) -> Candidate:
    predecessors = TASKS[task]["predecessors"]
    predecessor_finish = max((finished_at[pred] for pred in predecessors), default=0.0)
    start_min = max(node_available_at[node.node_id], predecessor_finish)
    runtime_min = runtime_on_node_min(task, node)
    finish_min = start_min + runtime_min
    utilization = float(TASKS[task]["utilization"])
    power_w = active_power_w(node, utilization)
    energy_kwh = power_w * (runtime_min / 60.0) / 1000.0
    carbon_intensity = interval_carbon_intensity(start_min, finish_min)
    carbon_g = energy_kwh * carbon_intensity
    return Candidate(
        task=task,
        node=node,
        start_min=start_min,
        finish_min=finish_min,
        runtime_min=runtime_min,
        power_w=power_w,
        energy_kwh=energy_kwh,
        carbon_g=carbon_g,
    )


def feasible_nodes(task: str, nodes: list[Node]) -> list[Node]:
    """Return all nodes that satisfy the task memory requirement."""
    required_memory = int(TASKS[task]["memory_gb"])
    feasible = [node for node in nodes if node.memory_gb >= required_memory]
    if not feasible:
        raise ValueError(f"No feasible node found for task {task!r}")
    return feasible


def normalize(value: float, values: list[float]) -> float:
    """Normalize a value to the local 0-1 range of candidate values."""
    minimum = min(values)
    maximum = max(values)
    if abs(maximum - minimum) < 1e-12:
        return 0.0
    return (value - minimum) / (maximum - minimum)


# ---------------------------------------------------------------------------
# DAG ordering and validation
# ---------------------------------------------------------------------------


def successors_by_task() -> dict[str, list[str]]:
    """Build a mapping from each task to its direct successors."""
    successors = {task: [] for task in TASKS}
    for task, spec in TASKS.items():
        for predecessor in spec["predecessors"]:
            successors[predecessor].append(task)
    return successors


def topological_order() -> list[str]:
    """Return a dependency-respecting topological order of workflow tasks."""
    incoming_count = {
        task: len(spec["predecessors"])
        for task, spec in TASKS.items()
    }
    successors = successors_by_task()
    ready = [task for task, count in incoming_count.items() if count == 0]
    order: list[str] = []

    while ready:
        task = ready.pop(0)
        order.append(task)
        for successor in successors[task]:
            incoming_count[successor] -= 1
            if incoming_count[successor] == 0:
                ready.append(successor)

    if len(order) != len(TASKS):
        raise ValueError("The workflow graph contains a cycle or invalid dependency")
    return order


def upward_ranks() -> dict[str, float]:
    """Compute HEFT-style upward ranks from task runtimes and successors."""
    successors = successors_by_task()
    memo: dict[str, float] = {}

    def rank(task: str) -> float:
        """Recursively compute the upward rank for one task.
        I tried to do this with dynamic programming."""
        if task in memo:
            return memo[task]
        own_runtime = float(TASKS[task]["base_runtime_min"])
        if not successors[task]:
            memo[task] = own_runtime
        else:
            memo[task] = own_runtime + max(rank(successor) for successor in successors[task])
        return memo[task]

    return {task: rank(task) for task in TASKS}


def scheduling_order() -> list[str]:
    """HEFT-inspired static priority order using upward rank."""
    ranks = upward_ranks()
    topo_index = {task: index for index, task in enumerate(topological_order())}
    return sorted(TASKS, key=lambda task: (-ranks[task], topo_index[task]))


# ---------------------------------------------------------------------------
# Schedulers
# ---------------------------------------------------------------------------


def baseline_scheduler(candidates: list[Candidate]) -> Candidate:
    """Choose the feasible node with the earliest predicted finish time."""
    return min(
        candidates,
        key=lambda candidate: (
            candidate.finish_min,
            candidate.energy_kwh,
            candidate.node.node_id,
        ),
    )


def green_aware_scheduler(
    candidates: list[Candidate],
    guardrail: float,
    weights: dict[str, float],
) -> Candidate:
    """Choose a node using a performance guardrail and weighted green score."""
    validate_weights(weights)
    fastest_finish = min(candidate.finish_min for candidate in candidates)
    max_allowed_finish = fastest_finish * (1.0 + guardrail)
    guarded_candidates = [
        candidate
        for candidate in candidates
        if candidate.finish_min <= max_allowed_finish + 1e-12
    ]
    if not guarded_candidates:
        guarded_candidates = candidates

    finish_values = [candidate.finish_min for candidate in guarded_candidates]
    energy_values = [candidate.energy_kwh for candidate in guarded_candidates]
    carbon_values = [candidate.carbon_g for candidate in guarded_candidates]

    def score(candidate: Candidate) -> float:
        """Compute the weighted normalized score for one guarded candidate."""
        normalized_time = normalize(candidate.finish_min, finish_values)
        normalized_energy = normalize(candidate.energy_kwh, energy_values)
        normalized_carbon = normalize(candidate.carbon_g, carbon_values)
        return (
            weights["time"] * normalized_time
            + weights["energy"] * normalized_energy
            + weights["carbon"] * normalized_carbon
        )

    return min(
        guarded_candidates,
        key=lambda candidate: (
            score(candidate),
            candidate.finish_min,
            candidate.energy_kwh,
            candidate.node.node_id,
        ),
    )


def validate_weights(weights: dict[str, float]) -> None:
    """Validate that objective weights are non-negative and sum to one."""
    required = {"time", "energy", "carbon"}
    if set(weights) != required:
        raise ValueError(f"weights must contain exactly {sorted(required)}")
    total = sum(weights.values())
    if abs(total - 1.0) > 1e-9:
        raise ValueError(f"weights must sum to 1.0, got {total}")
    for name, value in weights.items():
        if value < 0:
            raise ValueError(f"weight {name!r} must be non-negative")


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------


def run_simulation(
    scheduler_name: str,
    scheduler: Callable[[list[Candidate]], Candidate],
) -> tuple[list[ScheduleEntry], Summary]:
    """Run one full workflow simulation with the provided scheduler policy."""
    nodes = build_nodes()
    order = scheduling_order()
    node_available_at = {node.node_id: 0.0 for node in nodes}
    finished_at: dict[str, float] = {}
    schedule: list[ScheduleEntry] = []

    for task in order:
        candidates = [
            estimate_candidate(task, node, node_available_at, finished_at)
            for node in feasible_nodes(task, nodes)
        ]
        selected = scheduler(candidates)
        node_available_at[selected.node.node_id] = selected.finish_min
        finished_at[task] = selected.finish_min
        schedule.append(
            ScheduleEntry(
                scheduler=scheduler_name,
                task=task,
                node_id=selected.node.node_id,
                node_type=selected.node.node_type,
                node_role=selected.node.role,
                start_min=selected.start_min,
                finish_min=selected.finish_min,
                runtime_min=selected.runtime_min,
                energy_kwh=selected.energy_kwh,
                carbon_g=selected.carbon_g,
                memory_required_gb=int(TASKS[task]["memory_gb"]),
                node_memory_gb=selected.node.memory_gb,
            )
        )

    validate_schedule(schedule)
    summary = summarize_schedule(scheduler_name, schedule, len(nodes))
    return schedule, summary


def summarize_schedule(
    scheduler_name: str,
    schedule: list[ScheduleEntry],
    node_count: int,
) -> Summary:
    """Aggregate task-level schedule entries into workflow-level summary metrics."""
    makespan = max(entry.finish_min for entry in schedule)
    total_runtime = sum(entry.runtime_min for entry in schedule)
    return Summary(
        scheduler=scheduler_name,
        makespan_min=makespan,
        total_energy_kwh=sum(entry.energy_kwh for entry in schedule),
        total_carbon_g=sum(entry.carbon_g for entry in schedule),
        average_utilization=total_runtime / (node_count * makespan) if makespan > 0 else 0.0,
    )


def validate_schedule(schedule: list[ScheduleEntry]) -> None:
    """Check dependencies, timing, memory feasibility, and node overlap constraints."""
    by_task = {entry.task: entry for entry in schedule}
    if set(by_task) != set(TASKS):
        raise ValueError("Schedule does not contain exactly all workflow tasks")

    for entry in schedule:
        if entry.finish_min <= entry.start_min:
            raise ValueError(f"Task {entry.task} has invalid timing")
        for predecessor in TASKS[entry.task]["predecessors"]:
            predecessor_finish = by_task[predecessor].finish_min
            if predecessor_finish > entry.start_min + 1e-9:
                raise ValueError(
                    f"Dependency violation: {entry.task} starts before {predecessor} finishes"
                )
        if entry.memory_required_gb > entry.node_memory_gb:
            raise ValueError(f"Memory violation for task {entry.task}")

    by_node: dict[str, list[ScheduleEntry]] = {}
    for entry in schedule:
        by_node.setdefault(entry.node_id, []).append(entry)
    for node_id, entries in by_node.items():
        entries.sort(key=lambda item: item.start_min)
        for left, right in zip(entries, entries[1:]):
            if left.finish_min > right.start_min + 1e-9:
                raise ValueError(f"Node overlap on {node_id}")


def run_primary_experiment() -> tuple[
    list[ScheduleEntry],
    Summary,
    list[ScheduleEntry],
    Summary,
]:
    """Run the baseline and primary green-aware scheduler comparison."""
    baseline_schedule, baseline_summary = run_simulation(
        "baseline",
        baseline_scheduler,
    )
    green_schedule, green_summary = run_simulation(
        "green_aware",
        lambda candidates: green_aware_scheduler(
            candidates,
            guardrail=PRIMARY_GUARDRAIL,
            weights=PRIMARY_WEIGHTS,
        ),
    )
    green_summary = replace(
        green_summary,
        guardrail=PRIMARY_GUARDRAIL,
        weight_time=PRIMARY_WEIGHTS["time"],
        weight_energy=PRIMARY_WEIGHTS["energy"],
        weight_carbon=PRIMARY_WEIGHTS["carbon"],
    )
    return baseline_schedule, baseline_summary, green_schedule, green_summary


def run_sensitivity_analysis() -> list[Summary]:
    """Run green-aware simulations for all guardrail and weight configurations."""
    summaries: list[Summary] = []
    for guardrail in SENSITIVITY_GUARDRAILS:
        for weights in SENSITIVITY_WEIGHT_SETS:
            _, summary = run_simulation(
                "green_aware",
                lambda candidates, guardrail=guardrail, weights=weights: green_aware_scheduler(
                    candidates,
                    guardrail=guardrail,
                    weights=weights,
                ),
            )
            summaries.append(
                replace(
                    summary,
                    guardrail=guardrail,
                    weight_time=weights["time"],
                    weight_energy=weights["energy"],
                    weight_carbon=weights["carbon"],
                )
            )
    return summaries


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def percent_change(new_value: float, baseline_value: float) -> float:
    """Calculate percentage change from a baseline value."""
    if abs(baseline_value) < 1e-12:
        return 0.0
    return (new_value - baseline_value) / baseline_value * 100.0


def percent_saving(new_value: float, baseline_value: float) -> float:
    """Calculate percentage saving relative to a baseline value."""
    if abs(baseline_value) < 1e-12:
        return 0.0
    return (baseline_value - new_value) / baseline_value * 100.0


def write_schedule_csv(path: Path, schedule: list[ScheduleEntry]) -> None:
    """Write a task-level schedule to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ScheduleEntry.__dataclass_fields__))
        writer.writeheader()
        for entry in schedule:
            writer.writerow(entry.__dict__)


def write_summary_csv(path: Path, summaries: list[Summary]) -> None:
    """Write workflow-level summary rows to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Summary.__dataclass_fields__))
        writer.writeheader()
        for summary in summaries:
            writer.writerow(summary.__dict__)


def source_stage_for_task(task: str) -> str:
    """Return the source logical stage for a reduced-DAG task."""
    if task.startswith("Align_"):
        return "Align"
    if task.startswith("Filter_"):
        return "Filter"
    return task


def write_runtime_profile_csv(path: Path) -> None:
    """Write the task runtime, memory, and utilization profile to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "task",
        "source_stage",
        "predecessors",
        "base_runtime_min",
        "memory_gb",
        "utilization",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for task, spec in TASKS.items():
            writer.writerow(
                {
                    "task": task,
                    "source_stage": source_stage_for_task(task),
                    "predecessors": ";".join(spec["predecessors"]),
                    "base_runtime_min": spec["base_runtime_min"],
                    "memory_gb": spec["memory_gb"],
                    "utilization": spec["utilization"],
                }
            )


def write_input_profile_txt(path: Path, input_profile: InputProfile) -> None:
    """Write a summary of trace and carbon input data in text."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "Simulator input profile",
        "=======================",
        "",
        f"Runtime source: {input_profile.runtime_source}",
        f"Carbon source: {input_profile.carbon_source}",
        f"Carbon intensity points: {input_profile.carbon_points}",
    ]
    if input_profile.trace_files:
        lines.append("")
        lines.append("Trace files:")
        for trace_file in input_profile.trace_files:
            lines.append(f"- {trace_file}")
    if input_profile.carbon_file is not None:
        lines.append("")
        lines.append(f"Carbon intensity file: {input_profile.carbon_file}")

    if input_profile.stage_stats:
        trace_rows = input_profile.stage_stats.get("_trace_rows", {})
        lines.extend(
            [
                "",
                f"Mapped trace rows: {trace_rows.get('mapped_rows', 0):.0f}",
                f"Skipped trace rows: {trace_rows.get('skipped_rows', 0):.0f}",
                "",
                "Stage runtime profile:",
            ]
        )
        for stage in (
            "RawQC",
            "Trim",
            "Align",
            "Filter",
            "Merge",
            "MarkDuplicates",
            "PeakCalling",
            "Annotation",
            "Report",
        ):
            stats = input_profile.stage_stats[stage]
            lines.append(
                "- "
                f"{stage}: p50={stats['runtime_p50_min']:.2f} min, "
                f"p75={stats['runtime_p75_min']:.2f} min, "
                f"p95 memory={stats['memory_p95_gb']:.2f} GB, "
                f"median utilization={stats['utilization_p50']:.2f}"
            )

    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_sensitivity_csv(
    path: Path,
    baseline: Summary,
    summaries: list[Summary],
) -> None:
    """Write sensitivity results with percentage changes against the baseline."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "scheduler",
        "guardrail",
        "weight_time",
        "weight_energy",
        "weight_carbon",
        "makespan_min",
        "total_energy_kwh",
        "total_carbon_g",
        "average_utilization",
        "makespan_change_percent",
        "energy_saving_percent",
        "carbon_saving_percent",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for summary in summaries:
            writer.writerow(
                {
                    "scheduler": summary.scheduler,
                    "guardrail": summary.guardrail,
                    "weight_time": summary.weight_time,
                    "weight_energy": summary.weight_energy,
                    "weight_carbon": summary.weight_carbon,
                    "makespan_min": summary.makespan_min,
                    "total_energy_kwh": summary.total_energy_kwh,
                    "total_carbon_g": summary.total_carbon_g,
                    "average_utilization": summary.average_utilization,
                    "makespan_change_percent": percent_change(
                        summary.makespan_min,
                        baseline.makespan_min,
                    ),
                    "energy_saving_percent": percent_saving(
                        summary.total_energy_kwh,
                        baseline.total_energy_kwh,
                    ),
                    "carbon_saving_percent": percent_saving(
                        summary.total_carbon_g,
                        baseline.total_carbon_g,
                    ),
                }
            )


def write_comparison_txt(
    path: Path,
    baseline: Summary,
    green: Summary,
    input_profile: InputProfile,
) -> None:
    """Write a readable text summary of the primary scheduler comparison."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "Primary scheduler comparison",
        "============================",
        "",
        f"Runtime source: {input_profile.runtime_source}",
        f"Carbon source: {input_profile.carbon_source}",
        "",
        f"Baseline makespan: {baseline.makespan_min:.2f} min",
        f"Green-aware makespan: {green.makespan_min:.2f} min",
        f"Makespan change: {percent_change(green.makespan_min, baseline.makespan_min):.2f}%",
        "",
        f"Baseline energy: {baseline.total_energy_kwh:.4f} kWh",
        f"Green-aware energy: {green.total_energy_kwh:.4f} kWh",
        f"Energy saving: {percent_saving(green.total_energy_kwh, baseline.total_energy_kwh):.2f}%",
        "",
        f"Baseline carbon: {baseline.total_carbon_g:.2f} gCO2e",
        f"Green-aware carbon: {green.total_carbon_g:.2f} gCO2e",
        f"Carbon saving: {percent_saving(green.total_carbon_g, baseline.total_carbon_g):.2f}%",
        "",
        f"Green-aware guardrail: {green.guardrail:.0%}",
        (
            "Green-aware weights: "
            f"time={green.weight_time:.2f}, "
            f"energy={green.weight_energy:.2f}, "
            f"carbon={green.weight_carbon:.2f}"
        ),
        "",
        "Interpretation note:",
        (
            "These values are model-based estimates from a reduced workflow DAG. "
            "When trace input is enabled, task runtimes are trace-derived predictions, "
            "not exact measurements for the modeled heterogeneous nodes."
        ),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def print_summary(baseline: Summary, green: Summary) -> None:
    """Print the primary comparison metrics to the terminal."""
    print("\nPrimary scheduler comparison")
    print("============================")
    print(f"Baseline makespan:    {baseline.makespan_min:8.2f} min")
    print(f"Green-aware makespan: {green.makespan_min:8.2f} min")
    print(f"Makespan change:      {percent_change(green.makespan_min, baseline.makespan_min):8.2f}%")
    print()
    print(f"Baseline energy:      {baseline.total_energy_kwh:8.4f} kWh")
    print(f"Green-aware energy:   {green.total_energy_kwh:8.4f} kWh")
    print(f"Energy saving:        {percent_saving(green.total_energy_kwh, baseline.total_energy_kwh):8.2f}%")
    print()
    print(f"Baseline carbon:      {baseline.total_carbon_g:8.2f} gCO2e")
    print(f"Green-aware carbon:   {green.total_carbon_g:8.2f} gCO2e")
    print(f"Carbon saving:        {percent_saving(green.total_carbon_g, baseline.total_carbon_g):8.2f}%")
    print()
    print(f"Green-aware guardrail: {green.guardrail:.0%}")
    print(
        "Green-aware weights:  "
        f"time={green.weight_time:.2f}, "
        f"energy={green.weight_energy:.2f}, "
        f"carbon={green.weight_carbon:.2f}"
    )


def run_and_write_outputs(
    out_dir: Path,
    include_sensitivity: bool,
    input_profile: InputProfile,
) -> None:
    """Run the experiments and write all requested output artifacts."""
    baseline_schedule, baseline_summary, green_schedule, green_summary = run_primary_experiment()

    write_schedule_csv(out_dir / "schedule_baseline.csv", baseline_schedule)
    write_schedule_csv(out_dir / "schedule_green_aware.csv", green_schedule)
    write_summary_csv(out_dir / "summary_primary.csv", [baseline_summary, green_summary])
    write_comparison_txt(
        out_dir / "comparison_primary.txt",
        baseline_summary,
        green_summary,
        input_profile,
    )
    write_runtime_profile_csv(out_dir / "runtime_profile.csv")
    write_input_profile_txt(out_dir / "input_profile.txt", input_profile)

    if include_sensitivity:
        sensitivity = run_sensitivity_analysis()
        write_sensitivity_csv(out_dir / "summary_sensitivity.csv", baseline_summary, sensitivity)

    print_summary(baseline_summary, green_summary)
    print(f"\nRuntime source: {input_profile.runtime_source}")
    print(f"Carbon source:  {input_profile.carbon_source}")
    print(f"\nWrote results to: {out_dir}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for input source and output location."""
    parser = argparse.ArgumentParser(
        description="Run the reduced ChIP-seq scheduling simulator.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Directory for CSV and text outputs. "
            "Default: results_trace_driven for trace input, results_placeholder for placeholder input."
        ),
    )
    parser.add_argument(
        "--input-source",
        choices=["trace", "placeholder"],
        default="trace",
        help="Use trace-derived inputs or the original placeholder assumptions. Default: trace",
    )
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Also run guardrail/weight sensitivity analysis.",
    )
    return parser.parse_args()


def main() -> None:
    """Speaks for itself really."""
    args = parse_args()
    if args.input_source == "trace":
        input_profile = use_trace_inputs()
        out_dir = args.out_dir or Path("results_trace_driven")
    else:
        input_profile = reset_to_placeholder_inputs()
        out_dir = args.out_dir or Path("results_placeholder")

    topological_order()
    validate_weights(PRIMARY_WEIGHTS)
    run_and_write_outputs(out_dir, args.sensitivity, input_profile)


if __name__ == "__main__":
    main()
