"""Qt-independent parser and typed model for VorPy ``info.txt`` metadata."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any

@dataclass(frozen=True)
class Measurement:
    value: float | None
    unit: str | None = None
    raw: str = ""
    warning: str | None = None

@dataclass(frozen=True)
class ChainRecord:
    chain: str
    atoms: int | None = None
    residues: int | None = None
    volume: Measurement = Measurement(None)
    boundary_area: Measurement = Measurement(None)
    inter_chain_area: Measurement = Measurement(None)
    solvent_interfacial_area: Measurement = Measurement(None)
    raw: str = ""

@dataclass(frozen=True)
class ResidueRecord:
    chain: str
    name: str
    identifier: str = ""
    atoms: int | None = None
    volume: Measurement = Measurement(None)
    boundary_area: Measurement = Measurement(None)
    inter_residue_area: Measurement = Measurement(None)
    solvent_interfacial_area: Measurement = Measurement(None)
    raw: str = ""

@dataclass
class NetworkSummary:
    source: Path | None = None
    section_source: dict[str, str] = field(default_factory=dict)
    system: str | None = None
    group: str | None = None
    composition: dict[str, Measurement] = field(default_factory=dict)
    build: dict[str, str | Measurement] = field(default_factory=dict)
    timing: dict[str, Measurement] = field(default_factory=dict)
    network: dict[str, Measurement] = field(default_factory=dict)
    geometry: dict[str, Measurement] = field(default_factory=dict)
    curvature: dict[str, Measurement] = field(default_factory=dict)
    energy: dict[str, Measurement] = field(default_factory=dict)
    classification: dict[str, Measurement] = field(default_factory=dict)
    chains: list[ChainRecord] = field(default_factory=list)
    residues: list[ResidueRecord] = field(default_factory=list)
    extras: dict[str, list[tuple[str, str]]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    missing_sections: list[str] = field(default_factory=list)
    legacy: bool = True

    @property
    def network_type(self) -> str | None:
        value = self.build.get("network_type")
        return value if isinstance(value, str) else None

NETWORK_NAMES = {"aw": "Atomic Voronoi", "pow": "Power", "prm": "Primitive"}
_SECTION_KEYS = {
    "COMPOSITION": "composition", "BUILD INFORMATION": "build", "BUILD TIMING": "timing",
    "VORONOI NETWORK": "network", "GROUP GEOMETRY": "geometry", "SURFACE CURVATURE": "curvature",
    "SURFACE ENERGY ESTIMATE": "energy", "SURFACE CLASSIFICATION": "classification",
}
_FIELD_KEYS = {
    "atoms": "atoms", "residues": "residues", "chains": "chains", "mass": "mass",
    "vertices": "vertices", "edges": "edges", "surfaces": "surfaces", "volume": "volume",
    "van der waals volume": "vdw_volume", "surface area": "surface_area", "density": "density",
    "surface area / volume": "surface_area_volume", "integrated mean curvature": "int_mean_curv",
    "integrated mean curvature squared": "int_mean_curv_sq", "integrated gaussian curvature": "int_gauss_curv",
    "area-weighted mean curvature": "area_mean_curv", "area-weighted gaussian curvature": "area_gauss_curv",
    "representative surface energy": "surf_energy", "representative energy / area": "energy_per_area",
    "mapped group atoms": "mapped_group_atoms", "mapped surrounding waters": "mapped_surrounding_waters",
    "vertex time": "vertex_time", "connection time": "connection_time", "surface building time": "surface_time",
    "analysis time": "analysis_time", "total time": "total_time",
}
_NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|nan|unavailable", re.I)
_NUM_TEXT = r"(?:[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|nan|unavailable)"

def _measurement(text: str) -> Measurement:
    text = text.strip()
    match = _NUMBER.search(text)
    if not match: return Measurement(None, None, text, "No numeric value")
    token = match.group(0); unit = text[match.end():].strip() or None
    try: value = float(token) if token.lower() not in {"nan", "unavailable"} else None
    except ValueError: value = None
    warning = None if value is not None else f"Unavailable numeric value: {token}"
    return Measurement(value, unit, text, warning)

def _rows(lines: list[str]) -> tuple[dict[str, list[tuple[str,str]]], list[str]]:
    sections: dict[str, list[tuple[str,str]]] = {}; current = "Overview"; sections[current] = []; warnings=[]
    for raw in lines:
        line=raw.strip()
        if not line or set(line) <= {"=", "-"}: continue
        if ":" not in line and line.isupper() and len(line)<100:
            current=line.title(); sections.setdefault(current, []); continue
        if ":" in line:
            if current in {"Chain Composition", "Residue Composition"}:
                sections.setdefault(current, []).append((line, ""))
            else:
                key,value=line.split(":",1); key=key.strip(); value=value.strip()
                if key: sections.setdefault(current, []).append((key,value))
    return sections,warnings

def _record(text: str, residue: bool = False):
    if residue:
        match=re.match(r"(?P<name>\S+)(?:\s+(?P<id>[^ ]+))?\s+(?P<rest>.*)", text.strip())
        if not match: return None
        name=match.group("name"); identifier=match.group("id") or ""; rest=match.group("rest")
        chain="Unassigned"; atoms=None
        am=re.search(r"(?P<atoms>\d+)\s+atoms", rest, re.I); atoms=int(am.group("atoms")) if am else None
        return name,identifier,chain,atoms,rest
    match=re.match(r"(?P<chain>.*?)\s*:\s*(?P<rest>.*)", text.strip())
    if not match: return None
    chain=match.group("chain").strip() or "Unassigned"; rest=match.group("rest")
    am=re.search(r"(\d+)\s+atoms",rest,re.I); rm=re.search(r"(\d+)\s+residues",rest,re.I)
    return chain,int(am.group(1)) if am else None,int(rm.group(1)) if rm else None,rest

def parse_network_summary(path: Path) -> NetworkSummary:
    summary=NetworkSummary(source=path); raw_sections, warnings = _rows(path.read_text(encoding="utf-8", errors="replace").splitlines())
    summary.warnings.extend(warnings); summary.missing_sections=[name.title() for name in _SECTION_KEYS if name.title() not in raw_sections]
    for section, rows in raw_sections.items():
        section_upper = section.upper()
        summary.section_source[section] = str(path)
        for key,value in rows:
            low=key.lower(); canonical=_FIELD_KEYS.get(low)
            if section == "Overview":
                if low == "system": summary.system=value
                elif low == "group": summary.group=value
                else: summary.extras.setdefault(section, []).append((key,value))
            elif section_upper == "BUILD INFORMATION":
                build_key={"network type":"network_type", "surface resolution":"surface_resolution", "box size":"box_size", "maximum allowable vertex":"max_vertex"}.get(low)
                if build_key == "network_type": summary.build[build_key]=value.lower()
                elif build_key: summary.build[build_key]=_measurement(value)
                else: summary.extras.setdefault(section, []).append((key,value))
            elif canonical:
                target=getattr(summary, _SECTION_KEYS.get(section_upper, "extras"), None)
                if isinstance(target, dict): target[canonical]=_measurement(value)
            elif section_upper not in {"CHAIN COMPOSITION", "RESIDUE COMPOSITION"}: summary.extras.setdefault(section, []).append((key,value))
    for key,value in raw_sections.get("Chain Composition", []):
        text = key if not value else f"{key}: {value}"
        if re.match(r"chain\s*:", text, re.I):
            text = re.sub(r"^chain\s*:", ":", text, flags=re.I)
        parsed=_record(text)
        if parsed:
            chain,atoms,residues,rest=parsed
            def m(label):
                x=re.search(label+r"\s*:\s*("+_NUM_TEXT+r")\s*([^ ]*)",rest,re.I); return _measurement(x.group(0).split(":",1)[1]) if x else Measurement(None)
            summary.chains.append(ChainRecord(chain,atoms,residues,m("Volume"),m("Total Boundary SA"),m("Inter-Chain SA"),m("Solvent-Interfacial SA"),f"{key}: {value}"))
    for key,value in raw_sections.get("Residue Composition", []):
        text = key if not value else f"{key} {value}"
        parsed=_record(text, True)
        if parsed:
            name,identifier,chain,atoms,rest=parsed
            def m(label):
                x=re.search(label+r"\s*:\s*("+_NUM_TEXT+r")\s*([^ ]*)",rest,re.I); return _measurement(x.group(0).split(":",1)[1]) if x else Measurement(None)
            summary.residues.append(ResidueRecord(chain,name,identifier,atoms,m("Volume"),m("Total Boundary SA"),m("Inter-Residue SA"),m("Solvent-Interfacial SA"),f"{key} {value}"))
    return summary

def parse_group_info(path: Path) -> dict[str, list[tuple[str, str]]]:
    """Backward-compatible raw-section view."""
    sections, _ = _rows(path.read_text(encoding="utf-8", errors="replace").splitlines())
    return {name: values for name, values in sections.items() if values}
