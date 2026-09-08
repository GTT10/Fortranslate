#!/usr/bin/env python3
"""Reject ELF executables whose PT_GNU_STACK program header is executable."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path


ELF_MAGIC = b"\x7fELF"
PT_GNU_STACK = 0x6474E551
PF_X = 0x1


def stack_is_executable(path: Path) -> bool | None:
    with path.open("rb") as stream:
        identification = stream.read(16)
        if len(identification) < 16 or identification[:4] != ELF_MAGIC:
            return None
        elf_class = identification[4]
        data_encoding = identification[5]
        if data_encoding == 1:
            byte_order = "<"
        elif data_encoding == 2:
            byte_order = ">"
        else:
            raise AssertionError(f"{path}: unsupported ELF data encoding")

        if elf_class == 1:
            header_format = byte_order + "HHIIIIIHHHHHH"
            program_format = byte_order + "IIIIIIII"
            flags_index = 6
        elif elf_class == 2:
            header_format = byte_order + "HHIQQQIHHHHHH"
            program_format = byte_order + "IIQQQQQQ"
            flags_index = 1
        else:
            raise AssertionError(f"{path}: unsupported ELF class")

        header_size = struct.calcsize(header_format)
        header_data = stream.read(header_size)
        if len(header_data) != header_size:
            raise AssertionError(f"{path}: truncated ELF header")
        header = struct.unpack(header_format, header_data)
        program_offset = header[4]
        program_entry_size = header[8]
        program_count = header[9]
        expected_entry_size = struct.calcsize(program_format)
        if program_entry_size < expected_entry_size:
            raise AssertionError(f"{path}: invalid program-header size")

        for index in range(program_count):
            stream.seek(program_offset + index * program_entry_size)
            entry_data = stream.read(expected_entry_size)
            if len(entry_data) != expected_entry_size:
                raise AssertionError(f"{path}: truncated program header")
            entry = struct.unpack(program_format, entry_data)
            if entry[0] == PT_GNU_STACK:
                return bool(entry[flags_index] & PF_X)
    raise AssertionError(f"{path}: missing PT_GNU_STACK program header")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executables", nargs="+", type=Path)
    args = parser.parse_args()

    checked = 0
    failures: list[Path] = []
    for path in args.executables:
        if not path.is_file():
            raise AssertionError(f"missing executable: {path}")
        executable = stack_is_executable(path)
        if executable is None:
            continue
        checked += 1
        if executable:
            failures.append(path)
    if failures:
        names = ", ".join(str(path) for path in failures)
        raise AssertionError(f"executable stack requested by: {names}")
    print(f"non-executable stack: PASS ({checked} ELF executables)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
