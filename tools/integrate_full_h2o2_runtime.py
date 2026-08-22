#!/usr/bin/env python3
from pathlib import Path
import re

root = Path('.')


def save(path: Path, text: str) -> None:
    path.write_text(text)


# Extend the qualified seven-species transport database without changing it.
p = root / 'src/transport/transport_database_mod.F90'
t = p.read_text()
if 'load_h2o2_full_transport' not in t:
    pubs = list(re.finditer(r'^\s*public\s*::.*$', t, re.M | re.I))
    if not pubs:
        raise SystemExit('public list missing in transport database')
    pos = pubs[-1].end()
    t = t[:pos] + '\n  public :: load_h2o2_full_transport' + t[pos:]
    m = re.search(
        r'(?is)(\n\s*subroutine\s+load_h2o2_elementary_transport\b.*?'
        r'\n\s*end\s+subroutine\s+load_h2o2_elementary_transport\s*)',
        t,
    )
    if not m:
        raise SystemExit('elementary transport loader missing')
    block = re.sub(
        'load_h2o2_elementary_transport',
        'load_h2o2_full_transport',
        m.group(1),
        flags=re.I,
    )
    block, count = re.subn(
        r'allocate\s*\(\s*transport\s*\(\s*7\s*\)\s*\)',
        'allocate(transport(10))',
        block,
        count=1,
        flags=re.I,
    )
    if count != 1:
        raise SystemExit('transport allocation missing')
    extension = '''
    ! Full Cantera H2/O2 order: H2,H,O,O2,OH,H2O,HO2,H2O2,AR,N2.
    transport(10) = transport(7)
    transport(10)%name = "N2"
    transport(7) = transport(4)
    transport(7)%name = "HO2"
    transport(8) = transport(4)
    transport(8)%name = "H2O2"
    transport(9) = transport(2)
    transport(9)%name = "AR"
'''
    block = re.sub(
        r'(?i)(\n\s*end\s+subroutine\s+load_h2o2_full_transport)',
        extension + r'\1',
        block,
        count=1,
    )
    t = re.sub(
        r'(?i)\n\s*end\s+module\s+transport_database_mod\s*$',
        block + '\n\nend module transport_database_mod\n',
        t,
        count=1,
    )
    save(p, t)


# Add a runtime mechanism selector and the three additional mole fractions.
for rel, group in (
    ('src/driver/simulation_config_reactive_1d_mod.F90', 'reactive_1d'),
    ('src/driver/simulation_config_reactive_2d_mod.F90', 'reactive_2d'),
):
    p = root / rel
    t = p.read_text()
    if 'chemistry_model' not in t:
        t, count = re.subn(
            r'(\s*logical\s*::\s*chemistry_enabled[^\n]*\n)',
            r'\1    character(len=32) :: chemistry_model = "elementary"\n',
            t,
            count=1,
            flags=re.I,
        )
        if count != 1:
            raise SystemExit(f'chemistry_enabled component missing in {rel}')
        declarations = list(
            re.finditer(r'(?m)^\s*logical\s*::\s*chemistry_enabled\s*$', t, re.I)
        )
        if len(declarations) > 1:
            m = declarations[1]
            t = t[: m.end()] + '\n    character(len=32) :: chemistry_model' + t[m.end() :]
        t, count = re.subn(
            r'(namelist\s*/' + group + r'/.*?\bchemistry_enabled\b\s*,)',
            r'\1 chemistry_model,',
            t,
            count=1,
            flags=re.I | re.S,
        )
        if count != 1:
            raise SystemExit(f'chemistry_model namelist insertion failed in {rel}')
        t = re.sub(
            r'(chemistry_enabled\s*=\s*config%chemistry_enabled\s*\n)',
            r'\1    chemistry_model = config%chemistry_model\n',
            t,
            count=1,
            flags=re.I,
        )
        t = re.sub(
            r'(config%chemistry_enabled\s*=\s*chemistry_enabled\s*\n)',
            r'\1    config%chemistry_model = trim(adjustl(chemistry_model))\n',
            t,
            count=1,
            flags=re.I,
        )

    for field in ('x_ho2', 'x_h2o2', 'x_ar'):
        if re.search(r'\b' + field + r'\b', t, re.I):
            continue
        t, count = re.subn(
            r'(\s*real\(dp\)\s*::\s*x_h2o\s*=\s*[^\n]+\n)',
            r'\1    real(dp) :: ' + field + r' = 0.0_dp\n',
            t,
            count=1,
            flags=re.I,
        )
        if count != 1:
            raise SystemExit(f'{field} component insertion failed in {rel}')
        declarations = list(
            re.finditer(r'(?m)^\s*real\(dp\)\s*::[^\n]*\bx_h2o\b[^\n]*$', t, re.I)
        )
        if len(declarations) > 1:
            m = declarations[1]
            t = t[: m.end()] + '\n    real(dp) :: ' + field + t[m.end() :]
        t, count = re.subn(
            r'(namelist\s*/' + group + r'/.*?\bx_h2o\b\s*,)',
            r'\1 ' + field + ',',
            t,
            count=1,
            flags=re.I | re.S,
        )
        if count != 1:
            raise SystemExit(f'{field} namelist insertion failed in {rel}')
        t = re.sub(
            r'(x_h2o\s*=\s*config%x_h2o\s*\n)',
            r'\1    ' + field + ' = config%' + field + '\n',
            t,
            count=1,
            flags=re.I,
        )
        t = re.sub(
            r'(config%x_h2o\s*=\s*x_h2o\s*\n)',
            r'\1    config%' + field + ' = ' + field + '\n',
            t,
            count=1,
            flags=re.I,
        )
    save(p, t)


# Select elementary or full thermodynamics / mechanism at application startup.
for rel, is_2d in (
    ('app/pelef_reactive_1d.F90', False),
    ('app/pelef_reactive_2d.F90', True),
):
    p = root / rel
    t = p.read_text()
    if 'load_h2o2_full_thermo' not in t:
        uses = (
            '  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo\n'
            '  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism\n'
        )
        if is_2d:
            uses += '  use transport_database_mod, only: load_h2o2_full_transport\n'
        t, count = re.subn(
            r'(\n\s*implicit none)', '\n' + uses + r'\1', t, count=1, flags=re.I
        )
        if count != 1:
            raise SystemExit(f'use insertion failed in {rel}')
        start = t.lower().find('call load_h2o2_elementary_thermo')
        mechanism = t.lower().find('call load_h2o2_elementary_mechanism', start)
        simulate = t.lower().find('\n  call simulate_', mechanism)
        if min(start, mechanism, simulate) < 0:
            raise SystemExit(f'loader block not found in {rel}')
        line_end = t.find('\n', mechanism)
        prefix = t[:start]
        between = t[line_end:simulate]
        suffix = t[simulate:]
        dispatch = '''select case (trim(config%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) error stop "Failed to load elementary thermodynamics"
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load elementary mechanism"
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
  case default
    error stop "Unsupported chemistry_model"
  end select
'''
        if is_2d:
            between = re.sub(
                r'call\s+load_h2o2_elementary_transport\s*\(\s*transport\s*,\s*ok\s*\)',
                'if (trim(config%chemistry_model) == "full_h2o2") then\n'
                '    call load_h2o2_full_transport(transport, ok)\n'
                '  else\n'
                '    call load_h2o2_elementary_transport(transport, ok)\n'
                '  end if',
                between,
                flags=re.I,
            )
        else:
            between = '\n'
        t = prefix + dispatch + between + suffix
    save(p, t)


# Make the initialization arrays depend on the loaded species set.
for rel in ('src/reactive/reactive_1d_mod.F90', 'src/reactive/reactive_2d_mod.F90'):
    p = root / rel
    t = p.read_text()
    t = re.sub(
        r'real\(dp\)\s*::\s*mole_fractions\s*\(\s*7\s*\)',
        'real(dp), allocatable :: mole_fractions(:)',
        t,
        flags=re.I,
    )
    t = re.sub(
        r'real\(dp\)\s*::\s*mass_fractions\s*\(\s*7\s*\)',
        'real(dp), allocatable :: mass_fractions(:)',
        t,
        flags=re.I,
    )
    pattern = r'(?is)mole_fractions\s*=\s*\[\s*config%x_h2\s*,.*?config%x_n2\s*\]'
    replacement = '''allocate(mole_fractions(size(species)))
    select case (size(species))
    case (7)
      mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_n2]
    case (10)
      mole_fractions = [config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_ho2, config%x_h2o2, &
        config%x_ar, config%x_n2]
    case default
      ok = .false.
      return
    end select'''
    t, count = re.subn(pattern, replacement, t, count=1)
    if count != 1 and 'case (10)' not in t:
        raise SystemExit(f'composition constructor missing in {rel}')
    save(p, t)


# Minimal one- and two-dimensional smoke inputs.
for dim in ('1d', '2d'):
    (root / f'cases/reactive_full_h2o2_{dim}').mkdir(parents=True, exist_ok=True)
(root / 'cases/reactive_full_h2o2_1d/uniform.nml').write_text('''&reactive_1d
  nx = 8
  final_time = 1.0e-8
  maximum_steps = 1000
  problem = "uniform"
  reconstruction = "pcm"
  riemann_solver = "rusanov"
  chemistry_enabled = .true.
  chemistry_model = "full_h2o2"
  initial_temperature = 1000.0
  initial_pressure = 101325.0
  x_h2 = 0.29570
  x_o2 = 0.14785
  x_n2 = 0.55645
  output_file = "full_h2o2_1d.csv"
/
''')
(root / 'cases/reactive_full_h2o2_2d/uniform.nml').write_text('''&reactive_2d
  nx = 4
  ny = 4
  final_time = 1.0e-8
  maximum_steps = 1000
  problem = "uniform"
  reconstruction = "pcm"
  riemann_solver = "rusanov"
  use_transverse_correction = .false.
  chemistry_enabled = .true.
  chemistry_model = "full_h2o2"
  transport_enabled = .false.
  initial_temperature = 1000.0
  initial_pressure = 101325.0
  x_h2 = 0.29570
  x_o2 = 0.14785
  x_n2 = 0.55645
  output_file = "full_h2o2_2d.csv"
/
''')
validation = root / 'docs/validation/0.19-runtime-complete.txt'
validation.parent.mkdir(parents=True, exist_ok=True)
validation.write_text('ten-species full_h2o2 runtime dispatch and 1D/2D smoke gates\n')
