"""Post-reaction XP scoring of the optimised covalent complexes.

This is the stage that produces the apparent-affinity score reported in the
manuscript, and it is separate from covalent_complex_optimization.py, which
only minimises. It reuses the installed CovDock scoreComplexes stage, which:

  1. splits the reacted ligand from the receptor and caps it with hydrogen,
  2. minimises only the added hydrogens (OPLS4),
  3. mutates the reactive residue to alanine in the receptor copy,
  4. runs Glide XP optandscore against the matching CovDock grid.

The apparent-affinity score returned by the vendor is

    apparent affinity = 0.5 * (pre-reaction sampling DockingScore
                               + post-reaction DockingScore)

and this identity is asserted for every complex below. No de novo docking,
reattachment, coordinate repair, new state selection or global pose sampling is
performed: this is a fixed-pose sensitivity analysis, not a CovDock-LO rerun.

Run with Schrodinger 2021-2 python3. Input files must be named
TARGET_SANCxxxxx_refined.maegz.
"""
from pathlib import Path
from argparse import ArgumentParser, Namespace
import csv
import hashlib
import json
import math
import os
import time
import traceback
from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.application.prime.packages.covalent_docking import CovalentDocking


def graph(st):
    """Atoms and bonds, for checking the molecular graph is unchanged."""
    atoms = [(a.index, a.atomic_number, a.formal_charge) for a in st.atom]
    bonds = sorted((min(b.atom1.index, b.atom2.index),
                    max(b.atom1.index, b.atom2.index), b.order) for b in st.bond)
    return atoms, bonds


def xyz(st):
    return [(a.x, a.y, a.z) for a in st.atom]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_options(jobname):
    # Use the vendor's pose-prediction scoring defaults, rather than duplicating
    # their Glide keyword list. No actual leadopt sampling is invoked.
    opts = Namespace(mode='leadopt', jobname=jobname, ncluster_init=None,
                     ncluster=None, rotamer_sampling=None, dist_constraint=None,
                     score_dock_options=None, sample_dock_options=None,
                     min_options={'OPLS_VERSION': 'OPLS4'}, calc_strain=False,
                     calc_mmgbsa=False, affinity=True)
    engine = object.__new__(CovalentDocking)
    engine.setDefaultsByMode(opts)
    return engine, opts


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True,
                        help='csv listing the verified covalent complexes to rescore')
    parser.add_argument('--native-dir', type=Path, required=True,
                        help='directory holding the covalent docking job output')
    parser.add_argument('--only', nargs='*',
                        help='optional exact TARGET_SANCxxxxx keys')
    parser.add_argument('--dry-run', action='store_true',
                        help='validate inputs and write plans without calculating')
    args = parser.parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.native_dir = args.native_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.manifest.open()))
    bykey = {r['target'] + '_' + r['compound']: r for r in rows}
    paths = sorted(args.input_dir.glob('*_refined.maegz'))
    if args.only:
        paths = [p for p in paths
                 if p.name.replace('_refined.maegz', '') in args.only]
    if not paths:
        raise RuntimeError('no matching refined structures')
    results=[]
    for path in paths:
        key=path.name.replace('_refined.maegz','')
        row=bykey[key]
        native=args.native_dir/(key+'-out.maegz')
        grid=args.native_dir/(key+'_00_glidegrid-out.zip')
        native_ct=list(structure.StructureReader(str(native)))[int(row['record'])-1]
        cts=list(structure.StructureReader(str(path)))
        assert len(cts)==1,'Expected one preselected native pose per refined file'
        ct=cts[0]
        assert graph(ct)==graph(native_ct),'Refined and original chemical graphs/atom indices differ'
        assert grid.is_file(),'Missing matching native CovDock grid'
        lig=analyze.evaluate_asl(ct,CovalentDocking.LIG_ASL)
        assert len(lig)>0 and len(lig)<ct.atom_total
        sg=ct.atom[int(row['SG_index'])]
        assert sg.pdbname.strip()=='SG' and sg.element=='S'
        rec_asl='chain.name '+sg.chain+' and res.num '+str(sg.resnum)
        sample=native_ct.property['r_i_sample_docking_score']
        ct.property['r_i_sample_docking_score']=sample
        ct.property['r_i_sample_glide_gscore']=native_ct.property['r_i_sample_glide_gscore']
        assert 'r_psp_Prime_Energy' in ct.property,'Input must have a newly calculated Prime energy'
        engine,opts=make_options(key+'_refinement_score')
        item={'target':row['target'],'compound':row['compound'],'native_record':int(row['record']),
            'refined_input':str(path),'refined_input_sha256':sha(path),'native_input':str(native),
            'native_input_sha256':sha(native),'matching_grid':str(grid),'grid_sha256':sha(grid),
            'sample_docking_score':sample,'prime_energy_within_complex_only':ct.property['r_psp_Prime_Energy'],
            'score_options':opts.score_dock_options,'hydrogen_capping_forcefield':'OPLS4',
            'scope':'Fixed selected native pose refinement sensitivity, not complete leadopt sampling'}
        work=args.output_dir/key;work.mkdir(exist_ok=True)
        (work/'plan.json').write_text(json.dumps(item,indent=2))
        if args.dry_run:
            item['status']='planned_only';results.append(item);continue
        start=time.time()
        oldcwd=Path.cwd()
        try:
            os.chdir(str(work))
            before_graph=graph(ct);before_xyz=xyz(ct)
            # This exact installed function splits/caps the reacted ligand,
            # minimizes only added H, prepares the Ala-mutated receptor and
            # invokes XP optandscore against the matching CovDock grid.
            scored=engine.scoreComplexes(opts,1,[ct],rec_asl,str(grid))
            assert len(scored)==1
            ct=scored[0]
            assert graph(ct)==before_graph,'Scoring changed the covalent complex graph'
            assert xyz(ct)==before_xyz,'Scoring unexpectedly moved covalent-complex coordinates'
            post=ct.property['r_i_score_docking_score']
            aff=ct.property['r_cdock_cdock_affinity']
            assert math.isfinite(post) and post != 999, 'No valid postreaction score returned'
            assert abs(aff-0.5*(sample+post))<1e-9,'Vendor apparent-affinity arithmetic mismatch'
            item.update({'status':'completed','postreaction_docking_score':post,
                'postreaction_glidescore':ct.property.get('r_i_score_glide_gscore'),
                'apparent_affinity_score':aff,'elapsed_seconds':time.time()-start})
            ct.property['s_user_refinement_scope']=item['scope']
            ct.property['s_user_score_provenance']='Installed CovDock 1.3 scoreComplexes, leadopt scoring defaults, native grid, retained native sample score'
            ct.property['s_user_native_pose_source']=str(native)
            ct.property['i_user_native_pose_record']=int(row['record'])
            ct.write(str(work/(key+'_refined_scored.maegz')))
        except Exception as exc:
            item.update({'status':'failed','error':str(exc),'traceback':traceback.format_exc(),
                         'elapsed_seconds':time.time()-start})
        finally:
            os.chdir(str(oldcwd))
        (work/'result.json').write_text(json.dumps(item,indent=2))
        results.append(item)
        (args.output_dir/'rescore_results.json').write_text(json.dumps(results,indent=2))
        print(json.dumps({k:item[k] for k in ['target','compound','status','elapsed_seconds']}),flush=True)
    (args.output_dir/'rescore_results.json').write_text(json.dumps(results,indent=2))
    if any(r['status']=='failed' for r in results):raise SystemExit(2)


if __name__=='__main__':main()
