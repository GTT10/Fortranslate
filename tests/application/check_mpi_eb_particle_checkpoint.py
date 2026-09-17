"""Two-to-four-rank EB particle payload restart and rollback."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--exe', type=Path, required=True)
p.add_argument('--mpiexec', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
work = Path(tempfile.mkdtemp(prefix='particles_', dir=a.work.resolve()))
env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
path = work / 'saved'
for mode, ranks in [('write', 2), ('read', 4), ('reject', 4)]:
    target = path
    if mode == 'reject':
        target = work / 'truncated'
        target.write_bytes(path.read_bytes()[:-40])
    result = subprocess.run([str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), mode, str(target)],
                            env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=45)
    (work / f'{mode}.log').write_text(result.stdout)
    assert result.returncode == 0, result.stdout
print(f'PASS: immutable payload, MPI 2 to 4 restart, context/truncation rollback; {work}')
