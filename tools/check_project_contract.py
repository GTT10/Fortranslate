#!/usr/bin/env python3
"""Check version and pinned-reference records for release-critical drift."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


SHA_PATTERN = re.compile(r"[0-9a-f]{40}")


def one_match(pattern: re.Pattern[str], text: str, label: str) -> str:
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise AssertionError(f"{label}: expected one match, found {matches}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()

    cmake_text = (source / "CMakeLists.txt").read_text(encoding="utf-8")
    constants_text = (source / "src/core/constants_mod.F90").read_text(
        encoding="utf-8"
    )
    readme_text = (source / "README.md").read_text(encoding="utf-8")

    cmake_version = one_match(
        re.compile(r"(?m)^\s*VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)\s*$"),
        cmake_text,
        "CMake project version",
    )
    runtime_version = one_match(
        re.compile(r'pelef_version\s*=\s*"([0-9]+\.[0-9]+\.[0-9]+)"'),
        constants_text,
        "runtime version",
    )
    readme_version = one_match(
        re.compile(r"The `([0-9]+\.[0-9]+\.[0-9]+)` milestone contains"),
        readme_text,
        "README current milestone",
    )
    versions = {cmake_version, runtime_version, readme_version}
    if len(versions) != 1:
        raise AssertionError(
            "version drift: "
            f"CMake={cmake_version}, runtime={runtime_version}, "
            f"README={readme_version}"
        )
    validation = source / "docs/validation" / f"{cmake_version}.md"
    if not validation.is_file():
        raise AssertionError(f"missing validation record: {validation}")

    baseline_path = source / "references/pelec_baseline.json"
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if baseline.get("schema_version") != 1:
        raise AssertionError("unsupported PeleC baseline schema")
    pelec = baseline.get("pelec", {})
    commit = pelec.get("commit", "")
    if not SHA_PATTERN.fullmatch(commit):
        raise AssertionError("PeleC baseline commit must be a 40-digit SHA")
    if pelec.get("repository") != "https://github.com/Pele-Suite/PeleC.git":
        raise AssertionError("unexpected PeleC baseline repository")
    if pelec.get("branch") != "development":
        raise AssertionError("unexpected PeleC baseline branch")

    submodules = baseline.get("submodules")
    if not isinstance(submodules, list) or not submodules:
        raise AssertionError("PeleC baseline must include recursive submodules")
    paths: set[str] = set()
    for entry in submodules:
        path = entry.get("path", "")
        sha = entry.get("commit", "")
        repository = entry.get("repository", "")
        if not path or path in paths:
            raise AssertionError(f"invalid or duplicate submodule path: {path!r}")
        if not SHA_PATTERN.fullmatch(sha):
            raise AssertionError(f"invalid submodule SHA for {path}")
        if not repository.startswith("https://github.com/"):
            raise AssertionError(f"invalid submodule repository for {path}")
        paths.add(path)

    roadmap = (source / "docs/completion_roadmap.md").read_text(encoding="utf-8")
    mapping = (source / "docs/pelec_mapping.md").read_text(encoding="utf-8")
    for label, document in (("completion roadmap", roadmap), ("mapping", mapping)):
        if commit not in document:
            raise AssertionError(f"{label} does not cite the pinned PeleC SHA")

    presets = json.loads((source / "CMakePresets.json").read_text(encoding="utf-8"))
    preset_names = {entry["name"] for entry in presets.get("configurePresets", [])}
    required_presets = {"debug", "release", "debug-cantera", "mpi-debug", "mpi-release"}
    missing_presets = sorted(required_presets - preset_names)
    if missing_presets:
        raise AssertionError(f"missing configure presets: {missing_presets}")

    print(
        f"project contract: PASS (version {cmake_version}, "
        f"PeleC {commit[:12]})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
