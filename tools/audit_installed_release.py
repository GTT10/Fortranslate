#!/usr/bin/env python3
"""Audit a tests-disabled installed PeleF Release tree and its toolchain."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re
import subprocess

from check_no_executable_stack import stack_is_executable


ELF_MAGIC = b"\x7fELF"
EMPTY_ENVIRONMENT = {"PATH": "/usr/bin:/bin", "LC_ALL": "C"}


def run(command: list[str]) -> str:
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=EMPTY_ENVIRONMENT,
    )
    output = completed.stdout + completed.stderr
    if completed.returncode != 0:
        raise AssertionError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n{output}"
        )
    return output


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cache_value(cache_text: str, name: str) -> str:
    match = re.search(rf"^{re.escape(name)}:[^=]+=(.*)$", cache_text, re.MULTILINE)
    if not match:
        raise AssertionError(f"CMake cache lacks {name}")
    return match.group(1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--install-prefix", required=True, type=Path)
    parser.add_argument("--expected-elf-count", required=True, type=int)
    args = parser.parse_args()

    build_dir = args.build_dir.resolve()
    install_prefix = args.install_prefix.resolve()
    cache_path = build_dir / "CMakeCache.txt"
    manifest_path = build_dir / "install_manifest.txt"
    cache_text = cache_path.read_text(encoding="utf-8")
    manifest = manifest_path.read_text(encoding="utf-8").splitlines()
    if len(manifest) != args.expected_elf_count:
        raise AssertionError(
            f"install manifest has {len(manifest)} entries, "
            f"expected {args.expected_elf_count}"
        )

    expected_cache = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_Fortran_COMPILER": "/usr/bin/gfortran",
        "MPI_Fortran_COMPILER": "/usr/bin/mpifort",
        "MPIEXEC_EXECUTABLE": "/usr/bin/mpiexec",
        "PELEF_ENABLE_MPI": "ON",
        "PELEF_ENABLE_SUNDIALS": "OFF",
        "PELEF_ENABLE_TESTS": "OFF",
        "PELEF_ENABLE_CANTERA_REFERENCE": "OFF",
    }
    for name, expected in expected_cache.items():
        actual = cache_value(cache_text, name)
        if actual != expected:
            raise AssertionError(f"{name}={actual}, expected {expected}")
    if Path(cache_value(cache_text, "CMAKE_INSTALL_PREFIX")).resolve() != install_prefix:
        raise AssertionError("CMAKE_INSTALL_PREFIX does not match the audited tree")

    binaries = sorted((install_prefix / "bin").iterdir())
    binaries = [path for path in binaries if path.is_file()]
    if len(binaries) != args.expected_elf_count:
        raise AssertionError(
            f"installed bin has {len(binaries)} files, "
            f"expected {args.expected_elf_count}"
        )

    mpi_sonames: set[str] = set()
    for installed in binaries:
        if installed.read_bytes()[:4] != ELF_MAGIC:
            raise AssertionError(f"installed file is not ELF: {installed}")
        built = build_dir / installed.name
        if not built.is_file():
            raise AssertionError(f"build-tree counterpart is missing: {built}")
        if digest(installed) != digest(built):
            raise AssertionError(f"installed/build bytes differ: {installed.name}")
        if stack_is_executable(installed) is not False:
            raise AssertionError(f"invalid PT_GNU_STACK: {installed.name}")

        dynamic = run(["/usr/bin/readelf", "-d", str(installed)])
        if "(RPATH)" in dynamic or "(RUNPATH)" in dynamic:
            raise AssertionError(f"RPATH/RUNPATH found: {installed.name}")
        dependencies = run(["/usr/bin/ldd", str(installed)])
        if "not found" in dependencies:
            raise AssertionError(f"unresolved runtime dependency: {installed.name}")
        if installed.name.startswith("pelef_mpi_"):
            mpi_sonames.update(
                re.findall(r"\b(libmpi[^\s]*\.so\.\d+)\s+=>", dependencies)
            )
            if "libmpi.so.40" not in dependencies:
                raise AssertionError(f"OpenMPI libmpi.so.40 is missing: {installed.name}")
            if "libmpi.so.12" in dependencies:
                raise AssertionError(f"Intel MPI contamination: {installed.name}")

    compiler_files = sorted(
        (build_dir / "CMakeFiles").glob("*/CMakeFortranCompiler.cmake")
    )
    if len(compiler_files) != 1:
        raise AssertionError("expected one CMakeFortranCompiler.cmake")
    compiler_text = compiler_files[0].read_text(encoding="utf-8")
    audit_text = cache_text + "\n" + compiler_text
    if "oneapi" in audit_text.lower() or "libmpi.so.12" in audit_text:
        raise AssertionError("compiler/cache toolchain contamination")
    match = re.search(
        r'set\(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "([^"]*)"\)',
        compiler_text,
    )
    if not match:
        raise AssertionError("implicit Fortran link directories are missing")
    implicit_directories = match.group(1)
    expected_directories = (
        "/usr/lib/gcc/x86_64-linux-gnu/13;"
        "/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib"
    )
    if implicit_directories != expected_directories:
        raise AssertionError(
            f"unexpected implicit Fortran link directories: {implicit_directories}"
        )

    print(f"install_manifest_entries={len(manifest)}")
    print(f"installed_elf_count={len(binaries)}")
    print("installed_build_byte_identity=PASS")
    print(f"non_executable_stack_count={len(binaries)}")
    print("rpath_runpath=NONE")
    print("unresolved_runtime_dependencies=NONE")
    print(f"mpi_sonames={','.join(sorted(mpi_sonames))}")
    print(f"fortran_implicit_link_directories={implicit_directories}")
    print("oneapi_or_libmpi_so_12_contamination=NONE")
    print("installed_release_audit=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
