# SPDX-License-Identifier: LGPL-3.0-or-later
"""Portable CLI path propagation and compact-input protection; no dataset reads."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

from eeemcal_pi0.cli import main


class PathTests(unittest.TestCase):
    def test_resolved_paths_and_physics_reach_both_worker_scripts(self):
        with tempfile.TemporaryDirectory(prefix='pi0 paths ') as tmp:
            root = Path(tmp)
            (root / 'configs').mkdir()
            (root / 'inputs').mkdir()
            (root / 'inputs/manifest.txt').write_text('# example\n\nfile.root\n')
            (root / 'configs/physics.yaml').write_text('input:\n  inventory: cells.csv\n')
            plans = []
            for executor in ['local', 'slurm']:
                # Neither the missing ROOT file nor the inventory needs to be read in a dry run.
                command = [sys.executable, '-m', 'eeemcal_pi0',
                           '--input', 'inputs/manifest.txt', '--config', 'configs/physics.yaml',
                           '--output', executor + '/results', '--cache-dir', executor + '/cache',
                           '--work-dir', executor + '/scratch', '--executor', executor, '--dry-run']
                done = subprocess.run(command, cwd=root, capture_output=True, text=True)
                self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
                out = root / executor / 'results'
                planpath = out / 'batch/work_plan.json'
                plan = json.loads(planpath.read_text())
                plans.append(plan)
                self.assertEqual(plan['inputs'], [str(root / 'inputs/file.root')])
                self.assertEqual(plan['config']['input']['inventory'], str(root / 'configs/cells.csv'))
                for key, child in [('output', 'results'), ('cache', 'cache'), ('work', 'scratch')]:
                    self.assertEqual(plan[key], str(root / executor / child))
                self.assertEqual((out / 'batch/task_manifests/000000.txt').read_text(), str(root / 'inputs/file.root') + '\n')
                for kind in ['extract', 'analyze']:
                    argv = shlex.split((out / ('batch/' + kind + '.sh')).read_text().splitlines()[-1])
                    self.assertEqual(argv[0], 'exec')
                    self.assertEqual(argv[argv.index('--plan') + 1], str(planpath))
                    self.assertEqual(argv[argv.index('--worker-kind') + 1], kind)
                proposed = json.loads((out / 'batch/proposed_submissions.json').read_text())
                self.assertEqual(len(proposed), 2)
                self.assertTrue(all(any(arg.startswith('--output=' + str(out / 'logs/')) for arg in cmd) for cmd in proposed))
                self.assertFalse((out / 'batch/submissions.json').exists())
                self.assertFalse(list((root / executor / 'cache').rglob('*.npz')))
            self.assertEqual(plans[0]['config'], plans[1]['config'])

    def test_compact_input_overlap_refused_before_any_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            compact = root / 'compact'
            compact.mkdir()
            (compact / 'sentinel').write_bytes(b'unchanged compact input')
            alias = root / 'alias'
            alias.symlink_to(compact, target_is_directory=True)
            for output in [compact, compact / 'results', root, alias / 'results']:
                with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                    status = main(['--compact', str(compact), '--output', str(output), '--dry-run'])
                self.assertEqual(status, 2)
                self.assertEqual(list(compact.iterdir()), [compact / 'sentinel'])
                self.assertEqual((compact / 'sentinel').read_bytes(), b'unchanged compact input')
