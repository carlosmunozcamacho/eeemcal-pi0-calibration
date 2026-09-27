# SPDX-License-Identifier: LGPL-3.0-or-later
import copy
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from eeemcal_pi0.config import resolve,partition,cache_signature
from eeemcal_pi0.io import manifest,prepare,read_json,write_json
from eeemcal_pi0.candidates import load,CandidateSink
from eeemcal_pi0.extract import extract
from eeemcal_pi0.selection import compatible,choose
from eeemcal_pi0.normalization import calculate
from eeemcal_pi0 import batch
from fixture import reco

class ProducerTests(unittest.TestCase):
    def test_config_split_paths(self):
        c=resolve(overrides=['selection.energy_GeV=0.25','solver.min_usage=20','qa.halves=false'])
        self.assertEqual(c['selection']['energy_GeV'],.25);self.assertEqual(c['solver']['min_usage'],20);self.assertFalse(c['qa']['halves'])
        self.assertEqual(partition('immutable-file-id',9,c['split']),partition('immutable-file-id',9,c['split']))
        with self.assertRaises(ValueError):resolve(overrides=['solver.typo=2'])
        with self.assertRaises(ValueError):resolve(overrides=['split.training_fraction=.5'])
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'manifest').write_text('x.root\nx.root\n')
            with self.assertRaises(ValueError):manifest(p/'manifest')
            (p/'physical').mkdir();(p/'alias').symlink_to(p/'physical',target_is_directory=True)
            (p/'canonical').write_text(str(p/'alias/input.root')+'\n')
            self.assertEqual(manifest(p/'canonical'),[str(p/'alias/input.root')])
            with self.assertRaises(ValueError):prepare(p/'input',{},[p/'input'])
            prepare(p/'output',{},[])
            with self.assertRaises(FileExistsError):prepare(p/'output',{},[])
            with self.assertRaises(ValueError):prepare(p/'output',{'different':True},[],True)

    def test_truth_free_reco_and_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);expected=reco(p/'toy.root',40);c=resolve();out=p/'candidates.npz'
            extract([str(p/'toy.root')],out,c,CandidateSink(.15,cache_signature(c)))
            a=load(out);self.assertNotIn('truth_sim',a);np.testing.assert_array_equal(a['cell_id'],expected)
            self.assertEqual(a['cell_id'].dtype,np.uint64);compatible(a,c)
            selected,counts=choose(a,c,np.linspace(.04,.24,51));self.assertEqual(len(selected),40);self.assertEqual(counts['preselection_pairs'],40)
            c['selection']['energy_GeV']=.1
            with self.assertRaisesRegex(ValueError,'floor'):compatible(a,c)
            c=resolve();bad='_EcalEndcapNTrackClusterMatches_cluster/_EcalEndcapNTrackClusterMatches_cluster.index'
            reco(p/'bad.root',2,missing=bad)
            with self.assertRaisesRegex(ValueError,'track-match'):extract([str(p/'bad.root')],p/'bad.npz',c,CandidateSink(.15,cache_signature(c)))
            c['selection']['neutral']=False
            extract([str(p/'bad.root')],p/'no_veto.npz',c,CandidateSink(.15,cache_signature(c)))
            c['selection']['neutral']=True
            with self.assertRaisesRegex(ValueError,'track matching'):compatible(load(p/'no_veto.npz'),c)

    def test_normalization(self):
        ref={'source':'test bookkeeping','convention':'unweighted N/sigma','kind':'simulation','disjoint_generated_pieces':True,
             'components':[dict(piece_id='one',ensemble='same',weights='unit',unit_weight_convention_confirmed=True,
                 generated_events_represented=1000,count_convention='unskimmed',cross_section=dict(value=2,unit='nb'))]}
        x=calculate(ref,[1e33,1e34]);self.assertEqual(x['pb_inverse'],.5);self.assertEqual(x['nb_inverse'],500);self.assertAlmostEqual(x['equivalent_times'][0]['seconds'],500)
        self.assertEqual(calculate(None,[1e33])['status'],'unavailable')
        w=copy.deepcopy(ref);w['components'][0]['weights']='signed';self.assertEqual(calculate(w,[1e33])['status'],'unavailable')
        mix=copy.deepcopy(ref);mix['components'].append(dict(mix['components'][0],piece_id='two',ensemble='different'))
        self.assertNotIn('pb_inverse',calculate(mix,[1e33]))
        mix['components'][1]['ensemble']='same';mix['components'][1]['piece_id']='one';self.assertNotIn('pb_inverse',calculate(mix,[1e33]))
        d=dict(kind='data',source='run ledger',convention='recorded',recorded_luminosity=dict(value=2,unit='nb^-1'),delivered_luminosity=dict(value=3,unit='nb^-1'))
        self.assertEqual(calculate(d,[1e33])['pb_inverse'],.002)

    def test_scheduler_retry(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'batch').mkdir();(out/'logs').mkdir();site=batch.site_config(None);site.update(files_per_shard=1,files_per_task=1)
            p=batch.plan(out,['a','b','c'],None,out/'cache',out/'work',resolve(),site)
            write_json(out/'batch/submissions.json',[dict(kind='extract',job_id='100',tasks=[0,1,2]),dict(kind='analyze',job_id='101',tasks=[])])
            calls=[]
            def runner(cmd,**kw):
                calls.append(cmd)
                if cmd[0]=='squeue': return subprocess.CompletedProcess(cmd,0,'100_0|RUNNING\n100_1|PENDING\n101|PENDING\n','')
                return subprocess.CompletedProcess(cmd,0,'200\n' if '--array=2%16' in cmd else '201\n','')
            r=batch.submit(p,runner=runner)
            self.assertIn(['scancel','101'],calls);self.assertEqual(sum(x[0]=='squeue' for x in calls),1)
            submit=[x for x in calls if x[0]=='sbatch'];self.assertIn('--array=2%16',submit[0]);self.assertIn('--dependency=afterok:100_0:100_1:200',submit[1])
            self.assertIn('-r',calls[0]);self.assertNotIn(['scancel','100'],calls)
            with self.assertRaisesRegex(ValueError,'incomplete'):batch.worker(out/'batch/work_plan.json','analyze')
            # A compact-only workflow has exactly one submission and no extraction dependency.
            q=batch.plan(out,[],'/read-only/compact',None,out/'work',resolve(),site)
            (out/'batch/submissions.json').unlink();calls.clear();batch.submit(q,dry=True,runner=runner)
            self.assertFalse(calls);self.assertEqual(len(read_json(out/'batch/proposed_submissions.json')),1)

if __name__=='__main__':unittest.main()
