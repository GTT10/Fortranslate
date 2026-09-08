"""Compare generated build contracts before and after the one-time migration."""
import json
from pathlib import Path
import re
import subprocess
import sys

before, after, evidence = map(lambda s: Path(s).resolve(), sys.argv[1:4])
evidence.mkdir(parents=True, exist_ok=True)
configs = [('debug', 'Debug', False, True, False),
           ('release', 'Release', False, True, False),
           ('cantera', 'Debug', False, True, True),
           ('mpi-debug-selected', 'Debug', True, True, False),
           ('mpi-release-selected', 'Release', True, True, False),
           ('install-serial', 'Release', False, False, False),
           ('install-mpi', 'Release', True, False, False)]
report = []
for name, build_type, mpi, tests, cantera in configs:
    snapshots = []
    for label, source in [('before', before), ('after', after)]:
        build = evidence / label / name
        prefix = evidence / (label + '-install') / name
        args = ['cmake', '-S', str(source), '-B', str(build),
                '-G', 'Ninja', '-DCMAKE_Fortran_COMPILER=/usr/bin/gfortran',
                '-DCMAKE_BUILD_TYPE=' + build_type, '-DCMAKE_INSTALL_PREFIX=' + str(prefix),
                '-DPELEF_ENABLE_TESTS=' + ('ON' if tests else 'OFF'),
                '-DPELEF_ENABLE_MPI=' + ('ON' if mpi else 'OFF'),
                '-DPELEF_ENABLE_CANTERA_REFERENCE=' + ('ON' if cantera else 'OFF')]
        if mpi or not tests:
            args += ['-DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json',
                     '-DPELEF_MECHANISM_SOURCE=mechanisms/h2o2_cantera.yaml']
        if mpi:
            args += ['-DMPI_Fortran_COMPILER=/usr/bin/mpifort',
                     '-DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec']
        with (evidence / f'{label}-{name}.log').open('w') as log:
            subprocess.run(args, check=True, stdout=log, stderr=subprocess.STDOUT)
            subprocess.run(['cmake', '--graphviz=' + str(build / 'targets.dot'), str(build)],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
        def normalize(text):
            for path, token in [(build, '<BUILD>'), (prefix, '<PREFIX>'), (source, '<SOURCE>')]:
                text = text.replace(str(path), token)
            return re.sub(r' _BACKTRACE_TRIPLES "[^"]*"', '', text)
        snapshot = {str(p.relative_to(build)): normalize(p.read_text())
                    for pattern in ['CTestTestfile.cmake', 'cmake_install.cmake']
                    for p in build.rglob(pattern)}
        snapshot['targets.dot'] = normalize((build / 'targets.dot').read_text())
        tests_json = json.loads(subprocess.check_output(
            ['ctest', '--test-dir', str(build), '--show-only=json-v1'], text=True))
        snapshot['test_names'] = [test['name'] for test in tests_json['tests']]
        snapshots.append(snapshot)
    differences = [key for key in set(snapshots[0]) | set(snapshots[1])
                   if snapshots[0].get(key) != snapshots[1].get(key)]
    if differences:
        for label, snapshot in zip(['before', 'after'], snapshots):
            (evidence / f'{label}-{name}-manifest.json').write_text(json.dumps(snapshot, indent=2))
        raise RuntimeError(f'{name}: changed build contract: {differences}')
    record = {'configuration': name, 'tests': len(snapshots[0]['test_names']),
              'target_graph_identical': True, 'install_rules_identical': True,
              'ctest_commands_and_properties_identical': True}
    report.append(record)
    print(record, flush=True)
(evidence / 'layout-verification.json').write_text(json.dumps(report, indent=2) + '\n')
