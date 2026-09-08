"""Independent tests for the elementary mechanism generator."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock


GENERATOR_PATH = (
    Path(__file__).resolve().parents[2] / "tools" / "generate_elementary_mechanism.py"
)
PROJECT_ROOT = GENERATOR_PATH.parents[1]
PROBE_TEMPLATE_PATH = PROJECT_ROOT / "app" / "pelef_mechanism_probe.F90.in"
REACTOR_TEMPLATE_PATH = PROJECT_ROOT / "app" / "pelef0d_selected.F90.in"
REACTIVE_1D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_reactive_1d_selected.F90.in"
)
AMR_REACTIVE_1D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_amr_reactive_1d_selected.F90.in"
)
MPI_REACTIVE_1D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_mpi_reactive_1d_selected.F90.in"
)
REACTIVE_2D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_reactive_2d_selected.F90.in"
)
REACTIVE_EB_2D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_reactive_eb_2d_selected.F90.in"
)
REACTIVE_EB_AMR_2D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_reactive_eb_amr_2d_selected.F90.in"
)
REACTIVE_3D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_reactive_3d_selected.F90.in"
)
REACTIVE_EB_3D_TEMPLATE_PATH = (
    PROJECT_ROOT / "app" / "pelef_reactive_eb_3d_selected.F90.in"
)
GENERATOR_SPEC = importlib.util.spec_from_file_location(
    "generate_elementary_mechanism", GENERATOR_PATH
)
if GENERATOR_SPEC is None or GENERATOR_SPEC.loader is None:
    raise ImportError(f"cannot load generator from {GENERATOR_PATH}")
generator = importlib.util.module_from_spec(GENERATOR_SPEC)
GENERATOR_SPEC.loader.exec_module(generator)


class GenerateMechanismBundleTest(unittest.TestCase):
    @staticmethod
    def bundle():
        return {
            "chemistry_integrator": "explicit",
            "module_name": "test_bundle_mechanism_mod",
            "loader_name": "load_test_bundle_mechanism",
            "kernel_name": "test_bundle_production_rates",
            "thermo_loader_name": "load_test_bundle_thermo",
            "transport_loader_name": "load_test_bundle_transport",
            "symbol_prefix": "test_bundle",
            "species": ["A", "B"],
            "thermo": [
                {
                    "name": "A",
                    "molecular_weight": 2.0,
                    "composition": {"X": 1.0},
                    "temperature_min": 200.0,
                    "temperature_mid": 1000.0,
                    "temperature_max": 3000.0,
                    "low_coefficients": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
                    "high_coefficients": [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0],
                },
                {
                    "name": "B",
                    "molecular_weight": 2.0,
                    "composition": {"X": 1.0},
                    "temperature_min": 250.0,
                    "temperature_mid": 1100.0,
                    "temperature_max": 3200.0,
                    "low_coefficients": [15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0],
                    "high_coefficients": [22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0],
                },
            ],
            "transport": [
                {
                    "name": "A",
                    "geometry": "atom",
                    "well_depth": 100.0,
                    "diameter": 2.0,
                    "dipole": 0.0,
                    "polarizability": 1.0,
                    "rotational_relaxation": 1.0,
                },
                {
                    "name": "B",
                    "geometry": "nonlinear",
                    "well_depth": 120.0,
                    "diameter": 2.5,
                    "dipole": 0.1,
                    "polarizability": 1.5,
                    "rotational_relaxation": 2.0,
                },
            ],
            "reactions": [
                {
                    "equation": "A <=> B",
                    "reactants": {"A": 1.0},
                    "products": {"B": 1.0},
                    "reversible": True,
                    "arrhenius": {"A": 1.0, "b": 0.0, "Ea": 0.0},
                }
            ],
        }

    def test_two_species_bundle_emits_thermo_and_transport_loaders(self):
        data = self.bundle()
        generator.validate(data)

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "test_bundle.F90"
            output.write_text(generator.generate(data), encoding="utf-8")
            source = output.read_text(encoding="utf-8")

        self.assertIn("public :: load_test_bundle_thermo", source)
        self.assertIn("public :: load_test_bundle_transport", source)
        self.assertIn("public :: test_bundle_chemistry_integrator", source)
        self.assertIn('"explicit"', source)
        self.assertIn(
            "subroutine load_test_bundle_thermo(species, ok)",
            source,
        )
        self.assertIn("subroutine load_test_bundle_transport( &", source)
        self.assertIn('species(1)%name = "A"', source)
        self.assertIn('names(2) = "B"', source)
        self.assertIn("geometries(1) = 0", source)
        self.assertIn("geometries(2) = 2", source)

    def test_complete_bundle_requires_explicit_integrator_metadata(self):
        data = self.bundle()
        data.pop("chemistry_integrator")
        with self.assertRaisesRegex(
            ValueError, "requires chemistry_integrator 'explicit' or 'implicit'"
        ):
            generator.validate(data)

        for value in ("explicit", "implicit"):
            data = self.bundle()
            data["chemistry_integrator"] = value
            generator.validate(data)
            self.assertIn(f'"{value}"', generator.generate(data))

        data = self.bundle()
        data["chemistry_integrator"] = "auto"
        with self.assertRaisesRegex(
            ValueError, "requires chemistry_integrator 'explicit' or 'implicit'"
        ):
            generator.validate(data)

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is not installed")
    def test_max_length_equation_emits_bounded_fortran_and_compiles(self):
        data = self.bundle()
        data["reactions"][0]["equation"] = "A" * generator.MAX_EQUATION_LENGTH
        generator.validate(data)
        generated = generator.generate(data)
        self.assertLessEqual(
            max(len(line) for line in generated.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )

        with tempfile.TemporaryDirectory() as directory:
            output_directory = Path(directory)
            generated_source = output_directory / "max_equation.F90"
            generated_source.write_text(generated, encoding="utf-8")
            sources = [
                PROJECT_ROOT / "src/core/precision_mod.F90",
                PROJECT_ROOT / "src/physics/nasa7_thermo_mod.F90",
                PROJECT_ROOT / "src/physics/mixture_thermo_mod.F90",
                PROJECT_ROOT / "src/chemistry/elementary_kinetics_mod.F90",
                generated_source,
            ]
            compiler = shutil.which("gfortran")
            assert compiler is not None
            for index, source in enumerate(sources):
                result = subprocess.run(
                    [
                        compiler,
                        "-std=f2018",
                        f"-J{output_directory}",
                        f"-I{output_directory}",
                        "-c",
                        str(source),
                        "-o",
                        str(output_directory / f"module_{index}.o"),
                    ],
                    cwd=output_directory,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(
                    result.returncode,
                    0,
                    f"gfortran failed for {source}:\n"
                    f"{result.stdout}{result.stderr}",
                )

    @staticmethod
    def render_probe(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_KERNEL_NAME@": data["kernel_name"],
            "@PELEF_BUNDLE_JACOBIAN_NAME@": data.get(
                "jacobian_name", f"{data['kernel_name']}_jacobian"
            ),
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_REQUIRE_ACTIVITY_LITERAL@": ".false.",
        }
        rendered = PROBE_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactor(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_KERNEL_NAME@": data["kernel_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
        }
        rendered = REACTOR_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactive_1d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = REACTIVE_1D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_amr_reactive_1d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = AMR_REACTIVE_1D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_mpi_reactive_1d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = MPI_REACTIVE_1D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactive_3d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = REACTIVE_3D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactive_eb_3d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = REACTIVE_EB_3D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactive_2d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = REACTIVE_2D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactive_eb_2d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = REACTIVE_EB_2D_TEMPLATE_PATH.read_text(encoding="utf-8")
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @staticmethod
    def render_reactive_eb_amr_2d(data):
        replacements = {
            "@PELEF_BUNDLE_MODULE_NAME@": data["module_name"],
            "@PELEF_BUNDLE_NSPECIES_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nspecies"
            ),
            "@PELEF_BUNDLE_NREACTIONS_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_nreactions"
            ),
            "@PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL@": (
                f"{data.get('symbol_prefix', 'h2o2')}_chemistry_integrator"
            ),
            "@PELEF_BUNDLE_LOADER_NAME@": data["loader_name"],
            "@PELEF_BUNDLE_THERMO_LOADER_NAME@": data[
                "thermo_loader_name"
            ],
            "@PELEF_BUNDLE_TRANSPORT_LOADER_NAME@": data[
                "transport_loader_name"
            ],
            "@PELEF_BUNDLE_BUNDLE_SHA256@": "0" * 64,
        }
        rendered = REACTIVE_EB_AMR_2D_TEMPLATE_PATH.read_text(
            encoding="utf-8"
        )
        for marker, value in replacements.items():
            rendered = rendered.replace(marker, value)
        return rendered

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is not installed")
    def test_long_interfaces_compile_with_probe_and_reactor(self):
        data = self.bundle()
        data.update({
            "module_name": "m" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH,
            "loader_name": "l" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH,
            "kernel_name": "k" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH,
            "jacobian_name": "j" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH,
            "thermo_loader_name": "t" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH,
            "transport_loader_name": "u" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH,
            "symbol_prefix": "p" * 30,
            "source": {
                "format": "Cantera YAML",
                "file": "fixture.yaml",
                "sha256": "0" * 64,
                "cantera_version": "version_" + "c" * 180,
                "runtime_cantera_version": "runtime_" + "r" * 180,
                "runtime_cantera_git_commit": "commit_" + "g" * 180,
                "phase": "phase_" + "p" * 180,
                "source_units": {"length": "cm"},
                "target_units": {
                    "length": "m",
                    "time": "s",
                    "quantity": "kmol",
                    "activation-energy": "J/kmol",
                },
            },
        })
        generator.validate(data)
        generated = generator.generate(data)
        probe = self.render_probe(data)
        reactor = self.render_reactor(data)
        reactive_1d = self.render_reactive_1d(data)
        amr_reactive_1d = self.render_amr_reactive_1d(data)
        mpi_reactive_1d = self.render_mpi_reactive_1d(data)
        reactive_2d = self.render_reactive_2d(data)
        reactive_eb_2d = self.render_reactive_eb_2d(data)
        reactive_eb_amr_2d = self.render_reactive_eb_amr_2d(data)
        reactive_3d = self.render_reactive_3d(data)
        reactive_eb_3d = self.render_reactive_eb_3d(data)
        self.assertLessEqual(
            max(len(line) for line in generated.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in probe.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactor.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactive_1d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in amr_reactive_1d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in mpi_reactive_1d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactive_2d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactive_eb_2d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactive_eb_amr_2d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactive_3d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        self.assertLessEqual(
            max(len(line) for line in reactive_eb_3d.splitlines()),
            generator.MAX_FORTRAN_LINE_LENGTH,
        )
        for rendered in (
            probe,
            reactor,
            reactive_1d,
            amr_reactive_1d,
            mpi_reactive_1d,
            reactive_2d,
            reactive_eb_2d,
            reactive_eb_amr_2d,
            reactive_3d,
            reactive_eb_3d,
        ):
            self.assertNotIn("@PELEF_BUNDLE_", rendered)

        with tempfile.TemporaryDirectory() as directory:
            output_directory = Path(directory)
            generated_source = output_directory / "long_bundle.F90"
            probe_source = output_directory / "long_probe.F90"
            reactor_source = output_directory / "long_reactor.F90"
            reactive_1d_source = output_directory / "long_reactive_1d.F90"
            amr_reactive_1d_source = (
                output_directory / "long_amr_reactive_1d.F90"
            )
            mpi_reactive_1d_source = (
                output_directory / "long_mpi_reactive_1d.F90"
            )
            reactive_2d_source = output_directory / "long_reactive_2d.F90"
            reactive_eb_2d_source = (
                output_directory / "long_reactive_eb_2d.F90"
            )
            reactive_eb_amr_2d_source = (
                output_directory / "long_reactive_eb_amr_2d.F90"
            )
            reactive_3d_source = output_directory / "long_reactive_3d.F90"
            reactive_eb_3d_source = (
                output_directory / "long_reactive_eb_3d.F90"
            )
            generated_source.write_text(generated, encoding="utf-8")
            probe_source.write_text(probe, encoding="utf-8")
            reactor_source.write_text(reactor, encoding="utf-8")
            reactive_1d_source.write_text(reactive_1d, encoding="utf-8")
            amr_reactive_1d_source.write_text(
                amr_reactive_1d, encoding="utf-8"
            )
            mpi_reactive_1d_source.write_text(
                mpi_reactive_1d, encoding="utf-8"
            )
            reactive_2d_source.write_text(reactive_2d, encoding="utf-8")
            reactive_eb_2d_source.write_text(
                reactive_eb_2d, encoding="utf-8"
            )
            reactive_eb_amr_2d_source.write_text(
                reactive_eb_amr_2d, encoding="utf-8"
            )
            reactive_3d_source.write_text(reactive_3d, encoding="utf-8")
            reactive_eb_3d_source.write_text(
                reactive_eb_3d, encoding="utf-8"
            )
            sources = [
                PROJECT_ROOT / "src/core/precision_mod.F90",
                PROJECT_ROOT / "src/core/constants_mod.F90",
                PROJECT_ROOT / "src/physics/nasa7_thermo_mod.F90",
                PROJECT_ROOT / "src/physics/mixture_thermo_mod.F90",
                PROJECT_ROOT / "src/chemistry/elementary_kinetics_mod.F90",
                PROJECT_ROOT / "src/chemistry/constant_volume_reactor_mod.F90",
                PROJECT_ROOT / "src/driver/selected_composition_mod.F90",
                PROJECT_ROOT / "src/transport/gas_transport_mod.F90",
                PROJECT_ROOT
                / "src/driver/selected_mechanism_runtime_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_selected_reactor_mod.F90",
                PROJECT_ROOT / "src/core/state_indices_mod.F90",
                PROJECT_ROOT / "src/transport/mixture_transport_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_reactive_1d_mod.F90",
                PROJECT_ROOT / "src/hydro/slope_limiter_mod.F90",
                PROJECT_ROOT / "src/hydro/reconstruction_weno_mod.F90",
                PROJECT_ROOT / "src/reactive/reactive_1d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_hierarchy_1d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_multipatch_1d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_regrid_1d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_reactive_1d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_multilevel_reactive_1d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_multipatch_reactive_1d_mod.F90",
                PROJECT_ROOT
                / "src/driver/amr_reactive_1d_application_mod.F90",
                PROJECT_ROOT / "src/core/mesh_mod.F90",
                PROJECT_ROOT / "src/core/mesh_3d_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_reactive_2d_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_reactive_3d_mod.F90",
                PROJECT_ROOT / "src/boundary/reactive_boundary_2d_mod.F90",
                PROJECT_ROOT / "src/transport/reactive_transport_2d_mod.F90",
                PROJECT_ROOT / "src/reactive/reactive_2d_mod.F90",
                PROJECT_ROOT / "src/driver/reactive_2d_application_mod.F90",
                PROJECT_ROOT / "src/eb/eb_geometry_2d_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_reactive_eb_2d_mod.F90",
                PROJECT_ROOT / "src/eb/reactive_eb_cfl_2d_mod.F90",
                PROJECT_ROOT / "src/eb/eb_reactive_wall_flux_2d_mod.F90",
                PROJECT_ROOT
                / "src/eb/eb_reactive_redistribution_2d_mod.F90",
                PROJECT_ROOT
                / "src/eb/eb_reactive_reconstruction_2d_mod.F90",
                PROJECT_ROOT / "src/eb/eb_reactive_hydro_2d_mod.F90",
                PROJECT_ROOT / "src/eb/eb_reactive_transport_2d_mod.F90",
                PROJECT_ROOT / "src/driver/reactive_eb_2d_driver_mod.F90",
                PROJECT_ROOT
                / "src/driver/reactive_eb_2d_application_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_reactive_eb_amr_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_hierarchy_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_patch_tree_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_multilevel_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_flux_register_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_reactive_2d_mod.F90",
                PROJECT_ROOT
                / "src/amr/amr_eb_multilevel_reactive_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_regrid_2d_mod.F90",
                PROJECT_ROOT / "src/amr/amr_eb_transport_2d_mod.F90",
                PROJECT_ROOT
                / "src/amr/amr_eb_multilevel_transport_2d_mod.F90",
                PROJECT_ROOT
                / "src/amr/amr_eb_multipatch_transport_2d_mod.F90",
                PROJECT_ROOT
                / "src/amr/amr_eb_patch_tree_reactive_2d_mod.F90",
                PROJECT_ROOT
                / "src/driver/reactive_eb_amr_2d_driver_mod.F90",
                PROJECT_ROOT
                / "src/driver/reactive_eb_amr_2d_application_mod.F90",
                PROJECT_ROOT
                / "src/reactive/reactive_directional_flux_3d_mod.F90",
                PROJECT_ROOT / "src/transport/reactive_transport_3d_mod.F90",
                PROJECT_ROOT / "src/reactive/reactive_3d_mod.F90",
                PROJECT_ROOT
                / "src/problems/reactive_entropy_wave_3d_problem_mod.F90",
                PROJECT_ROOT / "src/io/reactive_csv_io_3d_mod.F90",
                PROJECT_ROOT / "src/driver/reactive_3d_application_mod.F90",
                PROJECT_ROOT / "src/eb/eb_geometry_3d_mod.F90",
                PROJECT_ROOT
                / "src/driver/simulation_config_reactive_eb_3d_mod.F90",
                PROJECT_ROOT / "src/eb/reactive_eb_cfl_3d_mod.F90",
                PROJECT_ROOT / "src/eb/eb_reactive_wall_flux_3d_mod.F90",
                PROJECT_ROOT
                / "src/eb/eb_reactive_redistribution_3d_mod.F90",
                PROJECT_ROOT / "src/eb/eb_reactive_hydro_3d_mod.F90",
                PROJECT_ROOT / "src/eb/eb_reactive_transport_3d_mod.F90",
                PROJECT_ROOT / "src/io/reactive_eb_3d_checkpoint_mod.F90",
                PROJECT_ROOT / "src/driver/reactive_eb_3d_driver_mod.F90",
                PROJECT_ROOT
                / "src/driver/reactive_eb_3d_application_mod.F90",
                generated_source,
                probe_source,
                reactor_source,
                reactive_1d_source,
                amr_reactive_1d_source,
                reactive_2d_source,
                reactive_eb_2d_source,
                reactive_eb_amr_2d_source,
                reactive_3d_source,
                reactive_eb_3d_source,
            ]
            compiler = shutil.which("gfortran")
            assert compiler is not None
            for index, source in enumerate(sources):
                result = subprocess.run(
                    [
                        compiler,
                        "-std=f2018",
                        f"-J{output_directory}",
                        f"-I{output_directory}",
                        "-c",
                        str(source),
                        "-o",
                        str(output_directory / f"module_{index}.o"),
                    ],
                    cwd=output_directory,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(
                    result.returncode,
                    0,
                    f"gfortran failed for {source}:\n"
                    f"{result.stdout}{result.stderr}",
                )

            mpi_compiler = shutil.which("mpifort")
            if mpi_compiler is not None:
                mpi_sources = [
                    PROJECT_ROOT / "src/parallel/mpi_domain_1d_mod.F90",
                    PROJECT_ROOT
                    / "src/parallel/mpi_reactive_transport_1d_mod.F90",
                    PROJECT_ROOT / "src/parallel/mpi_reactive_1d_mod.F90",
                    PROJECT_ROOT
                    / "src/driver/mpi_reactive_1d_application_mod.F90",
                    mpi_reactive_1d_source,
                ]
                for index, source in enumerate(mpi_sources):
                    result = subprocess.run(
                        [
                            mpi_compiler,
                            "-std=f2018",
                            f"-J{output_directory}",
                            f"-I{output_directory}",
                            "-c",
                            str(source),
                            "-o",
                            str(output_directory / f"mpi_module_{index}.o"),
                        ],
                        cwd=output_directory,
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(
                        result.returncode,
                        0,
                        f"mpifort failed for {source}:\n"
                        f"{result.stdout}{result.stderr}",
                    )

    def test_fortran_identifier_namespace_validation(self):
        data = self.bundle()
        data["module_name"] = "_invalid_module_name"
        with self.assertRaisesRegex(ValueError, "valid Fortran identifier"):
            generator.validate(data)

        for dependency_module in generator.FORTRAN_DEPENDENCY_MODULES:
            data = self.bundle()
            data["module_name"] = dependency_module
            with self.subTest(dependency_module=dependency_module):
                with self.assertRaisesRegex(ValueError, "dependency module"):
                    generator.validate(data)

        for field, imported_identifier in (
            ("module_name", "nasa7_species"),
            ("loader_name", "dp"),
            ("kernel_name", "elementary_production_rates"),
            ("jacobian_name", "elementary_mass_fraction_jacobian"),
        ):
            data = self.bundle()
            data[field] = imported_identifier
            with self.subTest(field=field, imported_identifier=imported_identifier):
                with self.assertRaisesRegex(ValueError, "imported Fortran identifier"):
                    generator.validate(data)

        for field, probe_identifier in (
            ("module_name", "pelef_mechanism_probe"),
            ("module_name", "pelef0d_selected"),
            ("module_name", "pelef_reactive_1d_selected"),
            ("module_name", "pelef_amr_reactive_1d_selected"),
            ("module_name", "pelef_mpi_reactive_1d_selected"),
            ("module_name", "pelef_reactive_2d_selected"),
            ("module_name", "pelef_reactive_eb_2d_selected"),
            ("module_name", "pelef_reactive_eb_amr_2d_selected"),
            ("module_name", "pelef_reactive_3d_selected"),
            ("module_name", "pelef_reactive_eb_3d_selected"),
            ("loader_name", "ieee_is_finite"),
            ("kernel_name", "species"),
            ("jacobian_name", "jacobian"),
            ("thermo_loader_name", "species_index"),
            ("transport_loader_name", "ok"),
            ("kernel_name", "mass_fractions"),
        ):
            data = self.bundle()
            data[field] = probe_identifier
            with self.subTest(field=field, probe_identifier=probe_identifier):
                with self.assertRaisesRegex(ValueError, "probe namespace"):
                    generator.validate(data)

        data = self.bundle()
        data["symbol_prefix"] = "pelef_selected"
        with self.assertRaisesRegex(ValueError, "probe namespace"):
            generator.validate(data)

        for field, generated_local in (
            ("loader_name", "reactions"),
            ("kernel_name", "temperature"),
            ("jacobian_name", "molar_production_rates"),
            ("transport_loader_name", "names"),
            ("transport_loader_name", "geometries"),
            ("transport_loader_name", "values"),
            ("kernel_name", "density"),
        ):
            data = self.bundle()
            data[field] = generated_local
            with self.subTest(field=field, generated_local=generated_local):
                with self.assertRaisesRegex(ValueError, "probe namespace"):
                    generator.validate(data)

    def test_omitted_jacobian_fallback_is_length_validated(self):
        data = self.bundle()
        data["kernel_name"] = "k" * generator.MAX_FORTRAN_IDENTIFIER_LENGTH
        with self.assertRaisesRegex(
            ValueError,
            "derived jacobian_name is not a valid Fortran identifier",
        ):
            generator.validate(data)

    def test_omitted_jacobian_fallback_accepts_63_character_name(self):
        data = self.bundle()
        data["kernel_name"] = "k" * (
            generator.MAX_FORTRAN_IDENTIFIER_LENGTH - len("_jacobian")
        )
        generator.validate(data)
        self.assertIn(
            f"public :: {data['kernel_name']}_jacobian",
            generator.generate(data),
        )

    def test_case_insensitive_interface_names_are_validated(self):
        data = self.bundle()
        data.update({
            "module_name": "MixedCaseMechanism",
            "loader_name": "LoadMixedCaseMechanism",
            "kernel_name": "MixedCaseProductionRates",
            "jacobian_name": "MixedCaseJacobian",
            "thermo_loader_name": "LoadMixedCaseThermo",
            "transport_loader_name": "LoadMixedCaseTransport",
        })
        generator.validate(data)

        data = self.bundle()
        data["loader_name"] = "SpEcIeS"
        with self.assertRaisesRegex(ValueError, "probe namespace"):
            generator.validate(data)

    def test_species_bundle_validation_rejects_incompatible_records(self):
        cases = (
            ("thermo order", "thermo", "order mismatch"),
            ("transport order", "transport", "order mismatch"),
        )
        for label, records, message in cases:
            with self.subTest(label=label):
                data = self.bundle()
                data[records][0]["name"] = "B"
                with self.assertRaisesRegex(ValueError, message):
                    generator.validate(data)

        data = self.bundle()
        data["thermo"][0]["low_coefficients"][0] = "nan"
        with self.assertRaisesRegex(ValueError, "nonfinite value"):
            generator.validate(data)

        data = self.bundle()
        data["thermo"][1]["high_coefficients"] = [0.0] * 6
        with self.assertRaisesRegex(ValueError, "must contain seven values"):
            generator.validate(data)

        data = self.bundle()
        data["transport"][1]["geometry"] = "spherical"
        with self.assertRaisesRegex(ValueError, "unsupported transport geometry"):
            generator.validate(data)

        data = self.bundle()
        data["species"] = [f"S{index}" for index in range(33)]
        data["reactions"][0] = {
            "equation": "S0 <=> S1",
            "reactants": {"S0": 1.0},
            "products": {"S1": 1.0},
            "arrhenius": {"A": 1.0, "b": 0.0, "Ea": 0.0},
        }
        with self.assertRaisesRegex(ValueError, "mechanism exceeds 32 species"):
            generator.validate(data)

        data = self.bundle()
        data["species"][1] = "bad name"
        with self.assertRaisesRegex(ValueError, "Fortran-safe identifiers"):
            generator.validate(data)

        data = self.bundle()
        data["thermo"][0]["reference_pressure"] = 1.0e5
        with self.assertRaisesRegex(ValueError, "reference pressure"):
            generator.validate(data)

        data = self.bundle()
        data["thermo"][0].pop("composition")
        with self.assertRaisesRegex(ValueError, "element composition"):
            generator.validate(data)

        data = self.bundle()
        data["thermo"][1]["composition"] = {"X": 2.0}
        with self.assertRaisesRegex(ValueError, "element X is not balanced"):
            generator.validate(data)

        data = self.bundle()
        data["thermo"][1]["molecular_weight"] = 3.0
        with self.assertRaisesRegex(ValueError, "reaction mass is not balanced"):
            generator.validate(data)

        data = self.bundle()
        data["transport_loader_name"] = data["thermo_loader_name"].upper()
        with self.assertRaisesRegex(ValueError, "identifiers must be distinct"):
            generator.validate(data)

        data = self.bundle()
        data["thermo_loader_name"] = "module"
        with self.assertRaisesRegex(ValueError, "valid loader names"):
            generator.validate(data)

        data = self.bundle()
        data["thermo_loader_name"] = "x" * 64
        with self.assertRaisesRegex(ValueError, "valid loader names"):
            generator.validate(data)

        data = self.bundle()
        data["thermo_loader_name"] = "test_bundle_nspecies"
        with self.assertRaisesRegex(ValueError, "names collide"):
            generator.validate(data)

        data = self.bundle()
        data["reactions"][0]["equation"] = "A" * 129
        with self.assertRaisesRegex(ValueError, "equation exceeds"):
            generator.validate(data)

        data = self.bundle()
        data["reactions"][0]["reversible"] = "yes"
        with self.assertRaisesRegex(ValueError, "reversible flag"):
            generator.validate(data)

        data = self.bundle()
        data["module_name"] = "module"
        with self.assertRaisesRegex(ValueError, "reserved Fortran keyword"):
            generator.validate(data)

        data = self.bundle()
        data["symbol_prefix"] = "x" * 60
        with self.assertRaisesRegex(ValueError, "identifier length limit"):
            generator.validate(data)

        data = self.bundle()
        data["reactions"][0]["equation"] = "A \\ B"
        with self.assertRaisesRegex(ValueError, "Fortran-literal safe"):
            generator.validate(data)

    def test_source_provenance_requires_target_unit_contract(self):
        data = self.bundle()
        data["source"] = {
            "format": "Cantera YAML",
            "file": "test.yaml",
            "sha256": "0" * 64,
            "cantera_version": "2.5.0",
            "runtime_cantera_version": "3.2.0",
            "runtime_cantera_git_commit": "4a8358e",
            "phase": "gas",
            "source_units": {"length": "cm"},
            "target_units": {
                "length": "m",
                "time": "s",
                "quantity": "kmol",
                "activation-energy": "J/kmol",
            },
        }
        generator.validate(data)
        data["source"]["target_units"]["quantity"] = "mol"
        with self.assertRaisesRegex(ValueError, "target units"):
            generator.validate(data)

    def test_atomic_publication_preserves_existing_generated_source(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "generated.F90"
            output.write_text("sentinel\n", encoding="utf-8")
            with mock.patch.object(generator.os, "replace", side_effect=OSError("boom")):
                with self.assertRaisesRegex(OSError, "boom"):
                    generator.write_generated(output, "replacement\n")
            self.assertEqual(output.read_text(encoding="utf-8"), "sentinel\n")
            self.assertEqual([path.name for path in Path(directory).iterdir()], [
                "generated.F90"
            ])

    def test_legacy_kinetics_only_input_still_generates(self):
        data = self.bundle()
        for key in (
            "thermo",
            "transport",
            "thermo_loader_name",
            "transport_loader_name",
        ):
            data.pop(key)

        generator.validate(data)
        source = generator.generate(data)

        self.assertIn("subroutine load_test_bundle_mechanism", source)
        self.assertNotIn("load_test_bundle_thermo", source)
        self.assertNotIn("load_test_bundle_transport", source)
        self.assertNotIn("valid_nasa7_species", source)
        self.assertNotIn("test_bundle_chemistry_integrator", source)


if __name__ == "__main__":
    unittest.main()
