"""Active gas/held-energy repair followed by indexed changed-rank continuation.

The unit fixture supplies a conservative zero-divergence manufactured flux,
asserts actual repair against an equal-thermo nonexhausted control, and checks
accepted-step conservation and transactional failure rollback. Every invocation
also advances a complete serial authority and compares local MPI fields/fluxes.
This driver requires byte-identical checkpoint continuation through step two.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--preflag', action='append', default=[])
    parser.add_argument('--postflag', action='append', default=[])
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='active_repair_', dir=args.work.resolve()))
    exe = args.exe.resolve()
    records = []
    hashes = {}

    def run(clock, operation, ranks, incoming, outgoing):
        label = f'{clock}-{operation}-{ranks}'
        cmd = [args.mpiexec, args.numproc_flag, str(ranks), *args.preflag,
               str(exe), *args.postflag, 'repair_restart', clock, operation,
               str(incoming), str(outgoing)]
        result = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=60)
        log = result.stdout + result.stderr
        (work / f'{label}.log').write_text(log)
        record = dict(label=label, command=cmd, returncode=result.returncode)
        records.append(record)
        (work / 'commands.json').write_text(json.dumps(records, indent=2) + '\n')
        assert result.returncode == 0, (label, log)
        assert 'PASS_ACTIVE_HELD_REPAIR_RESTART' in log, (label, log)
        assert 'ACCEPTED_STEP' in log, (label, log)
        if operation != 'continue':
            assert 'ACTIVE_REPAIR' in log, (label, log)
        published = incoming if operation == 'write' else outgoing
        authority = Path(str(published) + '.serial')
        assert published.read_bytes() == authority.read_bytes(), label
        hashes[published.name] = digest(published)
        return published

    for clock in ('common', 'recursive'):
        first = work / f'{clock}-first.chk'
        run(clock, 'write', 2, first, work / 'unused')
        direct = run(clock, 'direct', 2, work / 'unused', work / f'{clock}-direct.chk')
        resumed = run(clock, 'continue', 3, first, work / f'{clock}-resumed3.chk')
        serial_resume = run(clock, 'continue', 1, Path(str(first) + '.serial'),
                            work / f'{clock}-resumed1.chk')
        assert direct.read_bytes() == resumed.read_bytes() == serial_resume.read_bytes(), clock
    result = dict(status='passed', work=str(work), executable=str(exe),
                  executable_sha256=digest(exe), checkpoint_sha256=hashes,
                  invocation_count=len(records),
                  contract='exact serial/MPI state, held energy, effective flux; byte-exact indexed restart 2->3 and serial->1')
    (work / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
