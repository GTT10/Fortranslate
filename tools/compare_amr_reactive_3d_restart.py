#!/usr/bin/env python3
"""Require exact public-output parity for a 3D AMR checkpoint restart."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require_equal(reference: Path, restarted: Path, label: str) -> None:
    reference_data = reference.read_bytes()
    restarted_data = restarted.read_bytes()
    if reference_data != restarted_data:
        raise SystemExit(
            f"{label} restart mismatch: "
            f"reference_sha256={digest(reference_data)}, "
            f"restart_sha256={digest(restarted_data)}"
        )
    print(f"{label}_sha256={digest(reference_data)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-coarse", type=Path, required=True)
    parser.add_argument("--reference-fine", type=Path, required=True)
    parser.add_argument("--restart-coarse", type=Path, required=True)
    parser.add_argument("--restart-fine", type=Path, required=True)
    args = parser.parse_args()

    require_equal(args.reference_coarse, args.restart_coarse, "coarse")
    require_equal(args.reference_fine, args.restart_fine, "fine")
    print("amr_reactive_3d_restart_exact_match=PASS")


if __name__ == "__main__":
    main()
