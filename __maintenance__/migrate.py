"""One-time, hash-guarded, behavior-preserving CMake layout migration."""
from pathlib import Path
import hashlib
import json

ROOT = Path.cwd()
SPECS = {
    'CMakeLists.txt': ('6dd0a535c816b2ef1ce32b5ddbb99422a1f751c0a55b8c5411900e4cfcabd25f', [
        (0, 103, None), (103, 265, 'cmake/targets/Core.cmake'),
        (265, 878, 'cmake/targets/SelectedRuntime.cmake'),
        (878, 1002, 'cmake/targets/Applications.cmake'),
        (1002, 1007, None), (1007, 1909, 'cmake/targets/MPI.cmake'),
        (1909, 2613, 'cmake/targets/InstalledTests.cmake')]),
    'tests/CMakeLists.txt': ('d6d56f1fa092bab4a8a89fe0ccc68266ac19411b6acc3f6c1edca2420dd0970f', [
        (0, 599, 'tests/cmake/CoreUnitTests.cmake'),
        (599, 1242, 'tests/cmake/ReactiveUnitTests.cmake'),
        (1242, 1565, 'tests/cmake/MPIUnitTests.cmake'),
        (1565, 2383, 'tests/cmake/BaseRegressions.cmake'),
        (2383, 3036, 'tests/cmake/MechanismContracts.cmake'),
        (3036, 3147, 'tests/cmake/SelectedFixtures.cmake'),
        (3147, 3968, 'tests/cmake/SelectedMPI.cmake'),
        (3968, 4514, 'tests/cmake/SelectedRegularFlow.cmake'),
        (4514, 5126, 'tests/cmake/SelectedEB.cmake'),
        (5126, 5701, 'tests/cmake/SelectedEBAMR.cmake'),
        (5701, 6361, 'tests/cmake/SelectedEBMultilevel.cmake'),
        (6361, 6916, 'tests/cmake/SelectedAMR3D.cmake'),
        (6916, 7361, 'tests/cmake/SelectedMPIAMR3D.cmake'),
        (7361, 7737, 'tests/cmake/SelectedAMR1D.cmake'),
        (7737, 8033, 'tests/cmake/SelectedCVODE.cmake'),
        (8033, 8168, 'tests/cmake/SelectedZeroD.cmake'),
        (8168, 8530, 'tests/cmake/FixedReactive.cmake'),
        (8530, 8914, 'tests/cmake/FixedReactiveEB.cmake'),
        (8914, 9138, 'tests/cmake/EBAndMPIContext.cmake')]),
}

def migrate():
    pending = {}
    for filename, (expected, pieces) in SPECS.items():
        original = (ROOT / filename).read_bytes()
        if hashlib.sha256(original).hexdigest() != expected:
            raise RuntimeError(f'{filename}: source changed; refusing migration')
        text = original.decode('utf-8')
        if 'CMAKE_CURRENT_LIST_' in text or 'CMAKE_PARENT_LIST_' in text:
            raise RuntimeError('list-file-sensitive behavior needs manual review')
        lines = text.splitlines(keepends=True)
        assembled, wrapper, offset = [], [], 0
        for start, end, module in pieces:
            if start != offset or end <= start:
                raise RuntimeError('noncontiguous layout')
            chunk = ''.join(lines[start:end])
            assembled.append(chunk)
            if module is None:
                wrapper.append(chunk)
            else:
                if (ROOT / module).exists():
                    raise RuntimeError(f'refusing to overwrite {module}')
                pending[module] = chunk
                wrapper.append(f'include("${{PROJECT_SOURCE_DIR}}/{module}")\n')
            offset = end
        if offset != len(lines) or ''.join(assembled).encode() != original:
            raise RuntimeError('expanded source differs from the original')
        pending[filename] = ''.join(wrapper)
    for filename, content in pending.items():
        path = ROOT / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding='utf-8', newline='\n')
    print(json.dumps({'unchanged_expanded_cmake': True, 'files': sorted(pending)}, indent=2))

if __name__ == '__main__':
    migrate()
