# License and attribution

The eeemcal-pi0 package is distributed under the GNU Lesser General Public
License, version 3 or (at your option) any later version
(**SPDX-License-Identifier: LGPL-3.0-or-later**). See [LICENSE](LICENSE) for
LGPLv3 and [COPYING](COPYING) for the GPLv3 terms incorporated by LGPLv3.
You may redistribute and modify this software under those terms. It is
provided without warranty, including implied warranties of merchantability
or fitness for a particular purpose; see the license texts for details.

The constrained kernel has documented ancestry in Axel Perez Ruiz's
`master_thesis/eeemcal_simulation/pi0_energy_calibration/pi0_calibration.ipynb`,
commit **c36e96ce463b3bfc29abac40243cf5825400c5d1**, through the EEEMCal project
solver. **The negative sign in D was already present in Axel's implementation**;
the project preserved it. [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
identifies the exact sources, adapted regions and substantive changes.

Scientific method: R. T. Jones et al., “A bootstrap method for gain calibration
and resolution determination of a lead-glass calorimeter,” Nuclear Instruments
and Methods A 566 (2006) 366–374,
[DOI 10.1016/j.nima.2006.07.061](https://doi.org/10.1016/j.nima.2006.07.061),
Section 3, equations (1)–(9). Section 4 resolution extraction is not implemented.
The article, upstream notebook/C++ files and third-party dependency sources
are not bundled. Dependencies retain their own licenses.
