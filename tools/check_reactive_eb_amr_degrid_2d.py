#!/usr/bin/env python3
"""Validate removal of an untagged reactive EB AMR fine patch."""

from __future__ import annotations

import argparse
from pathlib import Path

from check_reactive_eb_amr_2d import check_level, read_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    args = parser.parse_args()

    coarse_rows = read_rows(args.coarse, 12 * 12)
    check_level(
        coarse_rows, 12, 12, 0.0, 0.01, 0.0, 0.01, expected_time=2.0e-7
    )
    if args.fine.exists():
        raise AssertionError("inactive EB AMR fine patch wrote stale output")
    print("check_reactive_eb_amr_degrid_2d: PASS")


if __name__ == "__main__":
    main()
