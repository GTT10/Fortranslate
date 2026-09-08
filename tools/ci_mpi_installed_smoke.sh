#!/usr/bin/env bash
# Run from any directory; the workflow uses a fresh build-mpi/install-mpi pair.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$ROOT"

cmake --install build-mpi
# 17 fixed serial + 10 selected serial + 11 fixed MPI + 4 selected MPI.
test "$(find install-mpi/bin -maxdepth 1 -type f -executable | wc -l)" -eq 42
for installed in install-mpi/bin/*; do
  name=${installed##*/}
  cmp "build-mpi/$name" "$installed"
done
python3 tools/check_no_executable_stack.py install-mpi/bin/*
if readelf -d install-mpi/bin/pelef_mpi_reactive_1d_selected | grep -Eq 'RPATH|RUNPATH'; then
  echo "selected MPI install contains RPATH/RUNPATH"
  exit 1
fi
cmake -DPELEF_MPI_EXECUTABLE=${ROOT}/install-mpi/bin/pelef_mpi_reactive_1d_selected -P cmake/CheckMpiRuntimeDependencies.cmake
mkdir installed-eb-amr-smoke
(
  cd installed-eb-amr-smoke
  ../install-mpi/bin/pelef_reactive_eb_amr_2d ../cases/selected_reactive_eb_amr_2d/h2o2_full_fixed.nml
  ../install-mpi/bin/pelef_reactive_eb_amr_2d_selected ../cases/selected_reactive_eb_amr_2d/h2o2_full_selected.nml
  cmp full_h2o2_reactive_eb_amr_2d_coarse.csv selected_full_h2o2_reactive_eb_amr_2d_coarse.csv
  cmp full_h2o2_reactive_eb_amr_2d_fine.csv selected_full_h2o2_reactive_eb_amr_2d_fine.csv
)
mkdir installed-eb3d-smoke
(
  cd installed-eb3d-smoke
  ../install-mpi/bin/pelef_reactive_eb_3d_selected \
    ../cases/selected_reactive_eb_3d/h2o2_full_transport_chemistry.nml
  python3 ../tools/check_selected_reactive_eb_3d.py \
    --input selected_full_reactive_eb_3d_transport_chemistry.csv \
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2 \
    --nx 10 --ny 8 --nz 6 --final-time 2.0e-9 \
    --x-upper 1.0e-4 --y-upper 1.0e-4 --z-upper 1.0e-4 \
    --plane-axis x --plane-position 3.95e-5 \
    --expected-regular 288 --expected-cut 48 \
    --expected-covered 144 --activity-species H2 \
    --initial-mass-fraction 0.028502723329253223 \
    --minimum-change 1.0e-9 \
    --minimum-active-density-span 2.0e-2
)
mkdir installed-mpi-eb-tree-smoke
(
  cd installed-mpi-eb-tree-smoke
  /usr/bin/mpiexec --oversubscribe -n 1 \
    ../install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d \
    ../cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_fixed.nml \
    installed_fixed.csv
  for ranks in 1 2 4; do
    /usr/bin/mpiexec --oversubscribe -n "$ranks" \
      ../install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d_selected \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml \
      "installed_selected_np${ranks}.csv"
  done
  python3 ../tools/check_selected_mpi_reactive_eb_patch_tree_2d.py \
    --reference installed_fixed.csv \
    --candidate installed_selected_np1.csv \
    --candidate installed_selected_np2.csv \
    --candidate installed_selected_np4.csv \
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2 \
    --expected-levels 0 1 2 3 --final-time 1.0e-12
  /usr/bin/mpiexec --oversubscribe -n 1 \
    ../install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d \
    ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_fixed_reference.nml \
    restart_fixed_np1.csv
  for ranks in 1 2 4; do
    /usr/bin/mpiexec --oversubscribe -n "$ranks" \
      ../install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d_selected \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_selected_reference.nml \
      "restart_selected_np${ranks}.csv"
  done
  /usr/bin/mpiexec --oversubscribe -n 1 \
    ../install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d_selected \
    ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_checkpoint_stop.nml
  for ranks in 2 4; do
    /usr/bin/mpiexec --oversubscribe -n "$ranks" \
      ../install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d_selected \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml \
      "restarted_np${ranks}.csv"
  done
  python3 ../tools/check_selected_mpi_reactive_eb_patch_tree_2d_restart.py \
    --checkpoint selected_mpi_eb_restart.chk \
    --reference restart_fixed_np1.csv \
    --candidate restart_selected_np1.csv \
    --candidate restart_selected_np2.csv \
    --candidate restart_selected_np4.csv \
    --candidate restarted_np2.csv \
    --candidate restarted_np4.csv \
    --stopped selected_mpi_eb_restart_stopped.csv \
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2 \
    --expected-levels 0 1 2 3 --final-time 1.0e-12 \
    --bundle-sha256 f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3 \
    --integrator implicit \
    --composition 0.29570 1.0e-5 1.0e-5 0.14784 1.0e-5 0.0 0.0 0.0 0.0 0.55643
  python3 ../tools/check_selected_mpi_reactive_eb_patch_tree_2d_restart_failures.py \
    --launcher /usr/bin/mpiexec --numproc-flag=-n --ranks 1 \
    --executable \
      "${ROOT}/install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d_selected" \
    --input \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml \
    --checkpoint selected_mpi_eb_restart.chk \
    --work-directory invalid-corruptions
  python3 ../tools/check_selected_mpi_reactive_eb_patch_tree_2d_aliases.py \
    --launcher /usr/bin/mpiexec --numproc-flag=-n --ranks 1 \
    --executable \
      "${ROOT}/install-mpi/bin/pelef_mpi_reactive_eb_patch_tree_2d_selected" \
    --reference-input \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_selected_reference.nml \
    --checkpoint-input \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_checkpoint_stop.nml \
    --restart-input \
      ../cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml \
    --work-directory invalid-aliases
)
for ranks in 1 2 4; do
  /usr/bin/mpiexec --oversubscribe -n "$ranks" install-mpi/bin/pelef_mpi_reactive_1d "installed_fixed_np${ranks}.csv"
  /usr/bin/mpiexec --oversubscribe -n "$ranks" install-mpi/bin/pelef_mpi_reactive_1d_selected "installed_selected_np${ranks}.csv"
done
python3 tools/compare_mpi_reactive_1d.py installed_fixed_np1.csv installed_fixed_np2.csv installed_fixed_np4.csv installed_selected_np1.csv installed_selected_np2.csv installed_selected_np4.csv
