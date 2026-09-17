"""Byte-indexed PDB frames: retain offsets, not a trajectory's atom tables."""
from pathlib import Path


def index_pdb_frames(path, progress=None):
    ranges = []
    start = 0
    has_atoms = False
    size = max(Path(path).stat().st_size, 1)
    with open(path, 'rb') as source:
        number = 0
        while True:
            offset = source.tell()
            line = source.readline()
            if not line:
                if has_atoms:
                    ranges.append((start, offset))
                break
            record = line[:6].strip().upper()
            if record == b'MODEL' and has_atoms:
                ranges.append((start, offset))
                start, has_atoms = offset, False
            has_atoms |= record in {b'ATOM', b'HETATM'}
            if record in {b'END', b'ENDMDL'} and has_atoms:
                ranges.append((start, source.tell()))
                start, has_atoms = source.tell(), False
            number += 1
            if progress and number % 4096 == 0:
                progress('Indexing PDB frames', min(20, source.tell() * 20 // size))
    return tuple(ranges)


def read_frame_lines(path, frame_range, progress=None):
    start, end = frame_range
    with open(path, 'rb') as source:
        source.seek(start)
        number = 0
        while source.tell() < end:
            line = source.readline()
            if not line:
                raise ValueError('The trajectory changed while loading. Reopen the structure.')
            yield line.decode('utf-8', errors='replace')
            number += 1
            if progress and number % 4096 == 0:
                progress('Reading selected frame', 20 + (source.tell() - start) * 50 // max(end - start, 1))
