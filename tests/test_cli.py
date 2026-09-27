# SPDX-License-Identifier: LGPL-3.0-or-later
"""Installed CLI integration; no central dataset or scheduler required."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from fixture import reco

class CLITests(unittest.TestCase):
    def test_outside_tree_complete_outputs_and_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);reco(root/'toy.root');(root/'manifest').write_text(str(root/'toy.root')+'\n')
            command=[sys.executable,'-m','eeemcal_pi0','--input',str(root/'manifest'),'--output',str(root/'result'),'--cache-dir',str(root/'broad'),'--work-dir',str(root/'scratch'),'--set','qa.formats=[png]']
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            done=subprocess.run(command,cwd=root,env=env,capture_output=True,text=True)
            self.assertEqual(done.returncode,0,done.stdout+done.stderr)
            s=json.loads((root/'result/summary.json').read_text());self.assertEqual(s['normalization']['status'],'unavailable');self.assertEqual(s['active_parameters'],8)
            self.assertEqual(s['execution_status'],'PASS')
            for name in ['gains.csv','summary.json','provenance.json','config_resolved.yaml','inputs.txt','report.md','diagnostics/mass_spectra.csv','diagnostics/channel_status.csv','plots/gain_map.png','plots/pi0_mass_before_after_selection.png']:
                self.assertTrue((root/'result'/name).exists(),name)
            original_inputs=(root/'result/inputs.txt').read_text()
            refusal=subprocess.run(command,cwd=root,env=env,capture_output=True,text=True);self.assertEqual(refusal.returncode,2)
            again=subprocess.run(command+['--resume'],cwd=root,env=env,capture_output=True,text=True);self.assertEqual(again.returncode,0,again.stderr)
            self.assertIn('Already complete',again.stdout)
            self.assertEqual((root/'result/inputs.txt').read_text(),original_inputs)
            cache_before={str(f.relative_to(root/'broad')):f.read_bytes() for f in (root/'broad').rglob('*') if f.is_file()}
            # Retained cache supports a new selection without ROOT reading; deleting only the fixture proves it.
            (root/'toy.root').unlink()
            reuse=[sys.executable,'-m','eeemcal_pi0','--compact',str(root/'broad'),'--output',str(root/'reuse'),'--set','qa.halves=false','--set','qa.formats=[png]']
            r=subprocess.run(reuse,cwd=root,env=env,capture_output=True,text=True);self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual(cache_before,{str(f.relative_to(root/'broad')):f.read_bytes() for f in (root/'broad').rglob('*') if f.is_file()})
            st=json.loads((root/'reuse/summary.json').read_text())['stability'];self.assertEqual(st['status'],'NOT_EVALUATED');self.assertIsNone(st['below10'])
            regen=subprocess.run([sys.executable,'-m','eeemcal_pi0','--stage','plots','--output',str(root/'reuse')],cwd=root,env=env,capture_output=True,text=True);self.assertEqual(regen.returncode,0,regen.stderr)
            nopeak=subprocess.run([sys.executable,'-m','eeemcal_pi0','--compact',str(root/'broad'),'--output',str(root/'nopeak'),'--set','peak.criteria.entries=1000000','--set','qa.formats=[png]'],cwd=root,env=env,capture_output=True,text=True)
            self.assertEqual(nopeak.returncode,3,nopeak.stderr)
            failure=json.loads((root/'nopeak/summary.json').read_text());self.assertEqual(failure['numerical_status'],'FAIL');self.assertTrue((root/'nopeak/report.md').exists())

if __name__=='__main__':unittest.main()
