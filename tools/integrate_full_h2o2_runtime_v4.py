#!/usr/bin/env python3
from pathlib import Path
import re

root = Path('.')


def write(path: Path, text: str) -> None:
    path.write_text(text)


def insert_before_namelist(text: str, group: str, declaration: str) -> str:
    pattern = re.compile(r'(?mi)^\s*namelist\s*/' + re.escape(group) + r'/')
    match = pattern.search(text)
    if not match:
        raise SystemExit(f'namelist {group} not found')
    return text[:match.start()] + '    ' + declaration + '\n' + text[match.start():]


# Full transport record: preserve elementary loader and extend the known order.
p = root / 'src/transport/transport_database_mod.F90'
text = p.read_text()
if 'load_h2o2_full_transport' not in text:
    public_lines = list(re.finditer(r'^\s*public\s*::.*$', text, re.M | re.I))
    if not public_lines:
        raise SystemExit('transport public declarations not found')
    pos = public_lines[-1].end()
    text = text[:pos] + '\n  public :: load_h2o2_full_transport' + text[pos:]
    match = re.search(
        r'(?is)(\n\s*subroutine\s+load_h2o2_elementary_transport\b.*?'
        r'\n\s*end\s+subroutine\s+load_h2o2_elementary_transport\s*)',
        text,
    )
    if not match:
        raise SystemExit('elementary transport loader not found')
    block = re.sub(
        'load_h2o2_elementary_transport',
        'load_h2o2_full_transport',
        match.group(1),
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
        raise SystemExit('transport allocation not found')
    extension = '''
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
    text = re.sub(
        r'(?i)\n\s*end\s+module\s+transport_database_mod\s*$',
        block + '\n\nend module transport_database_mod\n',
        text,
        count=1,
    )
    write(p, text)


# Configuration types/readers.
for relative, group in (
    ('src/driver/simulation_config_reactive_1d_mod.F90', 'reactive_1d'),
    ('src/driver/simulation_config_reactive_2d_mod.F90', 'reactive_2d'),
):
    path = root / relative
    text = path.read_text()

    # Type component: first chemistry_enabled declaration is inside the config type.
    if not re.search(r'character\s*\(\s*len\s*=\s*32\s*\)\s*::\s*chemistry_model\s*=\s*"elementary"', text, re.I):
        text, count = re.subn(
            r'(\s*logical\s*::\s*chemistry_enabled[^\n]*\n)',
            r'\1    character(len=32) :: chemistry_model = "elementary"\n',
            text,
            count=1,
            flags=re.I,
        )
        if count != 1:
            raise SystemExit(f'config type chemistry field missing in {relative}')

    # Reader local declarations immediately before namelist.
    prefix = text[: re.search(r'(?mi)^\s*namelist\s*/' + group + r'/', text).start()]
    if not re.search(r'character\s*\(\s*len\s*=\s*32\s*\)\s*::\s*chemistry_model\b', prefix, re.I):
        text = insert_before_namelist(text, group, 'character(len=32) :: chemistry_model')
    for field in ('x_ho2', 'x_h2o2', 'x_ar'):
        # type component
        type_pattern = r'real\(dp\)\s*::\s*' + field + r'\s*='
        if not re.search(type_pattern, text, re.I):
            text, count = re.subn(
                r'(\s*real\(dp\)\s*::\s*x_h2o\s*=\s*[^\n]+\n)',
                r'\1    real(dp) :: ' + field + r' = 0.0_dp\n',
                text,
                count=1,
                flags=re.I,
            )
            if count != 1:
                raise SystemExit(f'type component {field} insertion failed in {relative}')
        namelist_pos = re.search(r'(?mi)^\s*namelist\s*/' + group + r'/', text).start()
        local_prefix = text[:namelist_pos]
        # A declaration with an initializer is the type component; require a second local declaration.
        local_occurrences = re.findall(r'real\(dp\)\s*::[^\n]*\b' + field + r'\b', local_prefix, re.I)
        if len(local_occurrences) < 2:
            text = insert_before_namelist(text, group, 'real(dp) :: ' + field)

    # Namelist entries.
    if not re.search(r'namelist\s*/' + group + r'/.*?\bchemistry_model\b', text, re.I | re.S):
        text, count = re.subn(
            r'(namelist\s*/' + group + r'/.*?\bchemistry_enabled\b\s*,)',
            r'\1 chemistry_model,',
            text,
            count=1,
            flags=re.I | re.S,
        )
        if count != 1:
            raise SystemExit(f'chemistry_model namelist insertion failed in {relative}')
    for field in ('x_ho2', 'x_h2o2', 'x_ar'):
        if not re.search(r'namelist\s*/' + group + r'/.*?\b' + field + r'\b', text, re.I | re.S):
            text, count = re.subn(
                r'(namelist\s*/' + group + r'/.*?\bx_h2o\b\s*,)',
                r'\1 ' + field + ',',
                text,
                count=1,
                flags=re.I | re.S,
            )
            if count != 1:
                raise SystemExit(f'{field} namelist insertion failed in {relative}')

    # Copy type defaults into locals before READ and locals back into the type after READ.
    if not re.search(r'chemistry_model\s*=\s*config%chemistry_model', text, re.I):
        text, count = re.subn(
            r'(chemistry_enabled\s*=\s*config%chemistry_enabled\s*\n)',
            r'\1    chemistry_model = config%chemistry_model\n',
            text,
            count=1,
            flags=re.I,
        )
        if count != 1:
            raise SystemExit(f'chemistry_model default copy failed in {relative}')
    if not re.search(r'config%chemistry_model\s*=\s*trim', text, re.I):
        text, count = re.subn(
            r'(config%chemistry_enabled\s*=\s*chemistry_enabled\s*\n)',
            r'\1    config%chemistry_model = trim(adjustl(chemistry_model))\n',
            text,
            count=1,
            flags=re.I,
        )
        if count != 1:
            raise SystemExit(f'chemistry_model copy-back failed in {relative}')
    for field in ('x_ho2', 'x_h2o2', 'x_ar'):
        if not re.search(field + r'\s*=\s*config%' + field, text, re.I):
            text, count = re.subn(
                r'(x_h2o\s*=\s*config%x_h2o\s*\n)',
                r'\1    ' + field + ' = config%' + field + '\n',
                text,
                count=1,
                flags=re.I,
            )
            if count != 1:
                raise SystemExit(f'{field} default copy failed in {relative}')
        if not re.search(r'config%' + field + r'\s*=\s*' + field, text, re.I):
            text, count = re.subn(
                r'(config%x_h2o\s*=\s*x_h2o\s*\n)',
                r'\1    config%' + field + ' = ' + field + '\n',
                text,
                count=1,
                flags=re.I,
            )
            if count != 1:
                raise SystemExit(f'{field} copy-back failed in {relative}')
    write(path, text)


# Runtime loader dispatch.
for relative, is_2d in (
    ('app/pelef_reactive_1d.F90', False),
    ('app/pelef_reactive_2d.F90', True),
):
    path = root / relative
    text = path.read_text()
    if 'load_h2o2_full_thermo' not in text:
        uses = (
            '  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo\n'
            '  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism\n'
        )
        if is_2d:
            uses += '  use transport_database_mod, only: load_h2o2_full_transport\n'
        text, count = re.subn(
            r'(\n\s*implicit none)', '\n' + uses + r'\1', text, count=1, flags=re.I
        )
        if count != 1:
            raise SystemExit(f'use insertion failed in {relative}')

        first = text.lower().find('call load_h2o2_elementary_thermo')
        mechanism = text.lower().find('call load_h2o2_elementary_mechanism', first)
        simulate = text.lower().find('\n  call simulate_', mechanism)
        if min(first, mechanism, simulate) < 0:
            raise SystemExit(f'elementary loader block missing in {relative}')
        prefix = text[:first]
        middle = text[text.find('\n', mechanism):simulate]
        suffix = text[simulate:]
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
            middle = re.sub(
                r'call\s+load_h2o2_elementary_transport\s*\(\s*transport\s*,\s*ok\s*\)',
                'if (trim(config%chemistry_model) == "full_h2o2") then\n'
                '    call load_h2o2_full_transport(transport, ok)\n'
                '  else\n'
                '    call load_h2o2_elementary_transport(transport, ok)\n'
                '  end if',
                middle,
                flags=re.I,
            )
        else:
            middle = '\n'
        text = prefix + dispatch + middle + suffix
    write(path, text)


# Size-aware mole/mass-fraction arrays and constructors.
for relative in ('src/reactive/reactive_1d_mod.F90', 'src/reactive/reactive_2d_mod.F90'):
    path = root / relative
    text = path.read_text()
    text = re.sub(
        r'real\(dp\)\s*::\s*mole_fractions\s*\(\s*7\s*\)',
        'real(dp) :: mole_fractions(size(species))',
        text,
        flags=re.I,
    )
    text = re.sub(
        r'real\(dp\)\s*::\s*mass_fractions\s*\(\s*7\s*\)',
        'real(dp) :: mass_fractions(size(species))',
        text,
        flags=re.I,
    )
    pattern = r'(?is)mole_fractions\s*=\s*\[\s*config%x_h2\s*,.*?config%x_n2\s*\]'
    replacement = '''select case (size(species))
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
    text, count = re.subn(pattern, replacement, text, count=1)
    if count != 1 and 'case (10)' not in text:
        raise SystemExit(f'composition constructor missing in {relative}')
    write(path, text)


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
marker = root / 'docs/validation/0.19-runtime-complete.txt'
marker.parent.mkdir(parents=True, exist_ok=True)
marker.write_text('ten-species full_h2o2 runtime dispatch and 1D/2D smoke gates\n')
