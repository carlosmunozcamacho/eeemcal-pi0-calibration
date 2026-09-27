# SPDX-License-Identifier: LGPL-3.0-or-later
import argparse
from pathlib import Path
from full_regression import archive
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();archive(a.source,a.output)
