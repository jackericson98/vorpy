"""Stream PDB trajectory frames through a reusable temporary input file."""

from pathlib import Path
from tempfile import TemporaryDirectory


def first_pdb_frame(lines, system=None):
    """Yield only frame 1 while counting atom-containing frames in a stream."""
    count = 0
    has_atoms = False
    for line in lines:
        record = line[:6].strip().upper()
        if record == 'MODEL' and has_atoms:
            count += 1
            has_atoms = False
        if count == 0:
            yield line
        has_atoms |= record in {'ATOM', 'HETATM'}
        if record in {'ENDMDL', 'END'} and has_atoms:
            count += 1
            has_atoms = False
    count += int(has_atoms)
    if system is not None:
        system.frame_count = count
        system.frame_index = 1
        system.loaded_frame_count = int(count > 0)


def count_pdb_frames(file):
    with open(file) as source:
        from types import SimpleNamespace
        metadata = SimpleNamespace()
        for _ in first_pdb_frame(source, metadata):
            pass
        return metadata.frame_count


def iter_pdb_frames(file):
    """Yield (ordinal, path) for each frame; paths are valid until advancing.

    MODEL/ENDMDL and concatenated END-delimited PDBs are supported. Only
    one frame is stored on disk at a time, and the source is read once.
    """
    with TemporaryDirectory(prefix="vorpy_frames_") as directory:
        frame_path = Path(directory) / (Path(file).stem + '.pdb')
        with open(file) as source:
            ordinal = 0
            output = None
            has_atoms = False
            try:
                for line in source:
                    record = line[:6].strip().upper()
                    if record == 'MODEL' and has_atoms:
                        output.close()
                        ordinal += 1
                        yield ordinal, str(frame_path)
                        output = None
                        has_atoms = False
                    if output is None:
                        output = frame_path.open('w')
                    output.write(line)
                    has_atoms |= record in {'ATOM', 'HETATM'}
                    if record in {'ENDMDL', 'END'} and has_atoms:
                        output.close()
                        ordinal += 1
                        yield ordinal, str(frame_path)
                        output = None
                        has_atoms = False
                if has_atoms:
                    output.close()
                    ordinal += 1
                    yield ordinal, str(frame_path)
                if not ordinal:
                    raise ValueError(f"No atom frames found in PDB file: {file}")
            finally:
                if output is not None:
                    output.close()
