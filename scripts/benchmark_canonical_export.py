"""Bounded canonical-log benchmark; it never solves or rebuilds a network."""

from __future__ import annotations

import statistics
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

from vorpy.src.output import visualization
from vorpy.src.results import (
    INTERFACE_LOG_QUANTITY_UNITS,
    InterfaceLogRow,
    Provenance,
    Status,
)


ROWS = 106
SIMPLEX_IDS = 2_000
REPEATS = 7


def _synthetic_rows():
    provenances = [
        Provenance(
            "synthetic 1BNA",
            details={
                "source_role": "canonical-interface",
                "unresolved_simplex_ids": tuple(
                    f"1:{index},{index + 1}" for index in range(SIMPLEX_IDS)
                ),
            },
        ),
        Provenance("synthetic 1BNA", details={"source_role": "cached-dual"}),
        Provenance("synthetic 1BNA", details={"source_role": "cached-alpha"}),
    ]
    return [
        InterfaceLogRow(
            "1BNA-A-B",
            "physical_partition_interface",
            "shared",
            "area",
            250.0 + index,
            INTERFACE_LOG_QUANTITY_UNITS["area"],
            Status.PARTIAL,
            provenances[index % len(provenances)],
        )
        for index in range(ROWS)
    ]


def _rows(records, *, deduplicate_provenance):
    cache = {}
    rows = []
    for record in records:
        if deduplicate_provenance and record.provenance is not None:
            key = id(record.provenance)
            if key in cache:
                row = record.as_csv_row(provenance_json=cache[key])
            else:
                row = record.as_csv_row()
                cache[key] = row['provenance']
        else:
            row = record.as_csv_row()
        rows.append(row)
    return rows


def _legacy_once(records):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "logs.csv"
        original = visualization._compact_result_rows
        timing = {}

        def factory(iface):
            started = time.perf_counter()
            rows = _rows(records, deduplicate_provenance=False)
            timing['preparation_seconds'] = time.perf_counter() - started
            return rows

        visualization._compact_result_rows = factory
        try:
            started = time.perf_counter()
            cpu_started = time.process_time()
            visualization._compact_log(SimpleNamespace(), path)
            elapsed = time.perf_counter() - started
            cpu = time.process_time() - cpu_started
        finally:
            visualization._compact_result_rows = original
        return {
            "wall_seconds": elapsed,
            "cpu_seconds": cpu,
            "preparation_seconds": timing['preparation_seconds'],
            "bytes": path.stat().st_size,
            "files": 1,
        }


def _compact_once(records):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "logs.csv"
        original = visualization._compact_result_rows
        timing = {}

        def factory(iface, **kwargs):
            started = time.perf_counter()
            rows = _rows(records, deduplicate_provenance=True)
            timing['preparation_seconds'] = time.perf_counter() - started
            return rows

        visualization._compact_result_rows = factory
        try:
            started = time.perf_counter()
            cpu_started = time.process_time()
            profile = visualization._compact_log_bundle(SimpleNamespace(), path)
            elapsed = time.perf_counter() - started
            cpu = time.process_time() - cpu_started
        finally:
            visualization._compact_result_rows = original
        return {
            "wall_seconds": elapsed,
            "cpu_seconds": cpu,
            "preparation_seconds": timing['preparation_seconds'],
            "bytes": profile["bytes"]["total"],
            "compact_bytes": (
                profile["bytes"]["compact_csv"]
                + profile["bytes"]["provenance_sidecar"]
            ),
            "files": profile["counts"]["files_written"],
            "profile": profile,
        }


def _median(values):
    return statistics.median(values)


def main():
    records = _synthetic_rows()
    for _ in range(2):
        _legacy_once(records)
        _compact_once(records)
    legacy = [_legacy_once(records) for _ in range(REPEATS)]
    compact = [_compact_once(records) for _ in range(REPEATS)]
    profile = compact[-1]["profile"]
    result = {
        "rows": ROWS,
        "synthetic_unresolved_simplex_ids": SIMPLEX_IDS,
        "repeats": REPEATS,
        "legacy": {
            "wall_median_seconds": _median(item["wall_seconds"] for item in legacy),
            "cpu_median_seconds": _median(item["cpu_seconds"] for item in legacy),
            "preparation_median_seconds": _median(
                item["preparation_seconds"] for item in legacy
            ),
            "bytes": legacy[-1]["bytes"],
            "files": 1,
        },
        "compact_bundle": {
            "wall_median_seconds": _median(item["wall_seconds"] for item in compact),
            "cpu_median_seconds": _median(item["cpu_seconds"] for item in compact),
            "preparation_median_seconds": _median(
                item["preparation_seconds"] for item in compact
            ),
            "bytes_including_legacy_compatibility": compact[-1]["bytes"],
            "compact_csv_plus_sidecar_bytes": compact[-1]["compact_bytes"],
            "files": compact[-1]["files"],
            "profile": profile,
        },
        "memory": "not measured; optional psutil is not required by VorPy",
    }
    print(result)


if __name__ == "__main__":
    main()
