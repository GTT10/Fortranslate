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
for film_mode in ['', 'films']:
    stem = film_mode or 'particles'
    path = work / f'{stem}_saved'
    modes = [('write', 2), ('read', 4), ('reject', 4)]
    if film_mode:
        modes.append(('normal', 4))
    for mode, ranks in modes:
        target = path
        operation = mode
        if mode == 'reject':
            target = work / f'{stem}_truncated'
            target.write_bytes(path.read_bytes()[:-40])
        elif mode == 'normal':
            target = work / 'invalid_normal'
            lines = path.read_text().splitlines()
            row = lines[-3].split()
            row[5] = '-1.0'
            lines[-3] = ' '.join(row)
            target.write_text('\n'.join(lines) + '\n')
            operation = 'reject'
        result = subprocess.run(
            [str(a.mpiexec.resolve()), '-n', str(ranks), str(a.exe.resolve()), operation, str(target), film_mode],
            env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=45)
        (work / f'{stem}_{mode}.log').write_text(result.stdout)
        assert result.returncode == 0, result.stdout
print(f'PASS: immutable particles/films, MPI 2 to 4 restart, context/truncation/geometry rollback; {work}')
