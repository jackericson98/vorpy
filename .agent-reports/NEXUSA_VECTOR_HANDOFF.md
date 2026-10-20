# VECTOR handoff to NEXUSA — 2026-10-19

- Classic interface logs are restored in the default interface bundle as
  `<net_type>_logs.csv`; canonical typed Results remain separate.
- Historical local samples verify the classic section order and older
  36/12/5/9 atom/surface/edge/vertex schema. Current committed output is the
  existing 47/17/17/11 header-compatible superset, not a fabricated format.
- Active incomplete cells retain `Complete Cell?=False`; missing curvature is
  blank rather than zero.
- Archive reload -> classic log export passes without `Network.build()`.
- Compact synthetic benchmark: legacy 101.2 ms and 1,097,266 bytes versus
  compact bundle 55.7 ms; compact CSV plus sidecar 42,425 bytes (96.1% less).
- Focused validation: 95 passed in 7.64 s; no real molecular solve.

Open decision for NEXUSA/ATLAS: certify whether consumers require the older
positional 36/12/5/9 layout or the current header-driven superset. Positional
compatibility is not claimed.
