"""Reproduce Phase 1 standard exports from an existing solved network.

By default this command reuses geometry without a network solve. The explicit
--build-interface option exercises the normal Interface.build() lifecycle.
"""
import argparse
import json
from pathlib import Path
import time

from vorpy.src.io import load_network
from vorpy.src.group import Group
from vorpy.src.interface.interface import Interface
from vorpy.src.interface.export import export_info
from vorpy.src.output.logs import write_interface_logs
from vorpy.src.inputs.logs import read_logs
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--network', required=True)
    first = parser.add_mutually_exclusive_group(required=True)
    first.add_argument('--group1-chains')
    first.add_argument('--group1-atoms', help='System indices, e.g. 0-31,35')
    second = parser.add_mutually_exclusive_group(required=True)
    second.add_argument('--group2-chains')
    second.add_argument('--group2-atoms')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--alpha', type=float, default=None)
    parser.add_argument('--build-interface', action='store_true')
    args = parser.parse_args()
    started = time.perf_counter()
    network = load_network(args.network)
    system = network.sys
    if system is None:
        raise ValueError('Saved network must retain its parent System metadata')
    groups = []
    for name, chains, atoms in (('Group_A', args.group1_chains, args.group1_atoms), ('Group_B', args.group2_chains, args.group2_atoms)):
        if chains is not None:
            chain_set = set(chains.split(','))
            indices = set(map(int, system.balls.index[system.balls['chain_name'].astype(str).str.strip().isin(chain_set)]))
        else:
            indices = set()
            for token in atoms.split(','):
                if '-' in token:
                    first_id, last_id = map(int, token.split('-', 1))
                    indices.update(range(first_id, last_id + 1))
                else:
                    indices.add(int(token))
            if not indices <= set(map(int, system.balls.index)):
                raise ValueError('Requested atoms are absent from the parent System')
        if not indices:
            raise ValueError(f'No System atoms match chains {chains}')
        group = Group(sys=system, name=name, settings=dict(network.settings), make_net=False, build_net=False, group_id=name)
        group.ball_ndxs = sorted(indices)
        groups.append(group)
    iface = Interface(system, *groups, name='phase1_interface_example')
    iface.net = network  # Reuse the complete solved physical geometry.
    iface.group1_indices, iface.group2_indices = map(set, (g.ball_ndxs for g in groups))
    iface.group2_name = groups[1].name
    iface.dir = str(Path(args.output_dir).resolve())
    Path(iface.dir).mkdir(parents=True, exist_ok=True)
    if args.build_interface:
        print('Building the normal Interface network from saved System geometry', flush=True)
        iface.make_net()
        try:
            iface.build()
        except (TypeError, ValueError, RuntimeError) as error:
            if iface.geometry_analysis is None:
                raise
            # Reproduction aid only: report a later water-stage failure while
            # exporting the already cached geometry; do not change the pipeline.
            warning = f'Post-geometry Interface.build failure: {type(error).__name__}: {error}'
            iface.geometry_analysis.metadata['build_warning'] = warning
            iface.geometry_analysis.unresolved.append(warning)
            print(warning, flush=True)
    print('Analyzing cached physical geometry; Phase 1 performs no network solves', flush=True)
    result = iface.analyze_geometry(alpha=args.alpha)
    export_info(iface)
    log = write_interface_logs(iface)
    for reader in (read_logs, read_logs2):
        parsed = reader(log, return_dict=True)
        assert parsed['interface geometry'] == json.loads(json.dumps(result.log_record()))
    Path(iface.dir, 'geometry_analysis.json').write_text(json.dumps(result.to_dict(), indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'output_dir': iface.dir, 'elapsed_seconds': time.perf_counter() - started,
                      'selection': result.alpha_selection, 'mean_certified': result.voronoi_side_1.get('mean_certified'),
                      'additional_network_solves': result.metadata['additional_network_solves']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
