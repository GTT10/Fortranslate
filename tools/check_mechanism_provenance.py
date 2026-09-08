#!/usr/bin/env python3
"""Check source bytes and optionally regenerate a bundle without changing inputs."""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any


OPTIONS = (
    "module_name", "loader_name", "kernel_name", "jacobian_name",
    "thermo_loader_name", "transport_loader_name", "symbol_prefix",
    "chemistry_integrator", "description",
)


def reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON value: {value}")


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_bundle(path: Path) -> dict[str, Any]:
    bundle = json.loads(path.read_text(encoding="utf-8"),
                        parse_constant=reject_constant,
                        object_pairs_hook=unique_object)
    if not isinstance(bundle, dict) or not isinstance(bundle.get("source"), dict):
        raise ValueError("bundle must contain a source object")
    source = bundle["source"]
    for key in ("file", "sha256", "phase"):
        if not isinstance(source.get(key), str) or not source[key].strip():
            raise ValueError(f"source.{key} must be a nonempty string")
    if re.fullmatch(r"[0-9a-f]{64}", source["sha256"]) is None:
        raise ValueError("source.sha256 must contain 64 lowercase hex digits")
    return bundle


def source_path_for(bundle_path: Path, bundle: dict[str, Any],
                    source_path: Path | None) -> Path:
    if source_path is not None:
        return source_path.resolve()
    name = bundle["source"]["file"]
    if Path(name).name != name or name in {".", ".."} or "\\" in name:
        raise ValueError("non-basename source.file requires an explicit --source")
    path = (bundle_path.parent / name).resolve()
    if path.parent != bundle_path.parent.resolve():
        raise ValueError("source symlink leaves bundle directory; use --source")
    return path


def regenerate(bundle: dict[str, Any], source: Path) -> dict[str, Any]:
    # Optional dependency: byte-only preflight works without Cantera.
    from ingest_cantera_mechanism import build_bundle
    return build_bundle(source, phase=bundle["source"]["phase"],
                        source_file=bundle["source"]["file"],
                        **{key: bundle[key] for key in OPTIONS if key in bundle})


def stable_json(value: dict[str, Any]) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False,
                      allow_nan=False) + "\n"


def check(bundle_path: Path, source_path: Path | None = None,
          regenerate_output: bool = False, output: Path | None = None) -> bool:
    bundle_path = bundle_path.resolve()
    bundle = read_bundle(bundle_path)
    source = source_path_for(bundle_path, bundle, source_path)
    if output is not None:
        if not regenerate_output:
            raise ValueError("--output requires --regenerate")
        if output.resolve() in {bundle_path, source}:
            raise ValueError("--output must not overwrite the bundle or source")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    expected = bundle["source"]["sha256"]
    matches = digest == expected
    print(f"source:   {source}")
    print(f"expected: {expected}\nactual:   {digest}")
    if not matches:
        print("FAIL: source bytes do not match bundle provenance", file=sys.stderr)
    if regenerate_output:
        generated = regenerate(bundle, source)
        before, after = stable_json(bundle), stable_json(generated)
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(after.encode("utf-8"))
        if before != after:
            print("FAIL: regenerated bundle differs from committed bundle",
                  file=sys.stderr)
            sys.stdout.writelines(difflib.unified_diff(
                before.splitlines(keepends=True), after.splitlines(keepends=True),
                fromfile=str(bundle_path), tofile="regenerated"))
            matches = False
    if matches:
        print("PASS: mechanism provenance" +
              (" and regeneration" if regenerate_output else ""))
    return matches


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--output", type=Path,
                        help="write regenerated candidate for review, not in-place")
    args = parser.parse_args(argv)
    try:
        return 0 if check(args.bundle, args.source, args.regenerate, args.output) else 1
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        print(f"mechanism provenance: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
