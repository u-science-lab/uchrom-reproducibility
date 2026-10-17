"""``uchrom_recon`` -> the frozen engine build ``dipfrozen`` (uchrom-recon with its diploid extension), which predates
the rename and installs as ``uchrom_nucdyn``.  The runners in ``lib/`` differ from that build's only by the import name;
recon.sbatch puts this folder ahead of the build."""
import sys

import uchrom_nucdyn

sys.modules[__name__] = uchrom_nucdyn
