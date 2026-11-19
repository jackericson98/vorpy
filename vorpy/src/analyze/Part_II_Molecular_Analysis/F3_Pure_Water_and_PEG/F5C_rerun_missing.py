"""Preview six missing/empty AW solvent exports; pass --run to rerun sequentially.

Uses original solute checkpoints, rather than reusing incomplete saved solvent
layers. New exports stay in this repository; existing scientific data are kept.
A successful process does not guarantee a complete layer: inspect layers.csv.
"""
import argparse
from pathlib import Path
import subprocess
import sys

REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents
                 if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(REPO_ROOT))
CHECKPOINT_ROOT = Path(
    r'C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack'
    r'\Ericson Voronoi DNA\P2\NewModels\PEG\Checkpoints\aw'
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import LEGACY_OUTPUT_ROOT

OUTPUT_ROOT = LEGACY_OUTPUT_ROOT / 'PEG/F5C/reruns'
CASES = (('PEG300', 3), ('PEG2000', 9), ('PEG4000', 2),
         ('PEG4000', 3), ('PEG4000', 7), ('PEG4000', 11))

from scripts.run_checkpoint_layers import build_command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true', help='Execute; default only validates and prints commands.')
    args = parser.parse_args()
    commands = []
    # Validate every input and destination before starting any expensive work.
    for model, number in CASES:
        frame = f'frame_{number:04d}'
        source = CHECKPOINT_ROOT / model / frame / f'{model}_checkpoint.vpy'
        output = OUTPUT_ROOT / model / frame
        if not source.is_file():
            raise FileNotFoundError(source)
        if output.exists():
            raise FileExistsError(f'Preserving existing rerun: {output}; choose a fresh OUTPUT_ROOT')
        command = build_command(source, output)
        commands.append(command)
        print(subprocess.list2cmdline(command), flush=True)
    if not args.run:
        print('Preview only. Pass --run to calculate four solvent layers for these six frames.')
        return
    for command in commands:
        subprocess.run(command, cwd=REPO_ROOT, check=True)


if __name__ == '__main__':
    main()
