"""Run with PyMOL's Python; check standalone PMLs and render their scenes."""
import argparse
import json
from pathlib import Path
import re
import shutil
import tempfile

from pymol import cmd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pml',nargs='+',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    results = []
    for pml in args.pml:
        name = pml.stem.removesuffix('_contacts')
        prefix = 'f3b_'+re.sub('[^A-Za-z0-9_]','_',name)
        cmd.reinitialize()
        cmd.set('max_threads',2)
        with tempfile.TemporaryDirectory() as directory:
            standalone = Path(directory)/pml.name
            shutil.copy2(pml,standalone)
            cmd.do('@'+str(standalone))
        assert set(cmd.get_scene_list()) == {'AW_Cell','Power_Cell','Contact_Environment'}
        assert cmd.count_atoms(prefix+'_central') == 1
        atoms = cmd.get_model(prefix+'_atoms').atom
        assert len({a.id for a in atoms}) == len(atoms)
        views = []
        for scene in ('AW_Cell','Power_Cell','Contact_Environment'):
            cmd.scene(scene,'recall',animate=0)
            enabled = set(cmd.get_names('objects',enabled_only=1))
            assert (prefix+'_aw_cell' in enabled) == (scene == 'AW_Cell')
            assert (prefix+'_pow_cell' in enabled) == (scene == 'Power_Cell')
            views.append(cmd.get_view())
            cmd.ray(1000,800)
            cmd.png(str((args.output/(name+'_'+scene+'.png')).resolve()))
        assert all(max(abs(a-b) for a,b in zip(views[0],v)) < 1e-6 for v in views)
        results.append(dict(example=name,pymol_version=cmd.get_version()[0],scenes=cmd.get_scene_list(),
                            unique_atoms=len(atoms),central_radius=cmd.get_model(prefix+'_central').atom[0].vdw,
                            shared_camera=True,standalone=True))
        print('REAL_PYMOL_SMOKE_OK',name,flush=True)
    (args.output/'smoke_results.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    main()
