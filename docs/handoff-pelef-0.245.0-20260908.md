# PeleF 0.245.0 completion qualification handoff (2026-09-08)

## 前提

- 対象: `/home/gtt13/projects/pelec`, branch `main`, HEAD `bc584fc65e7e945b98475cd1c0394b483c5cc685`。
- 目的: PeleC の Fortran 化である PeleF 0.245.0 を、実装、回帰試験、クリーンインストール監査、適用範囲の文書化まで進める。
- 作業木にはユーザーの既存変更が多数ある。reset、checkout、無関係な削除を行わない。2026-09-08のユーザー指示「ローカルに残さないでほしい」に基づくリモート退避を除き、commit、pushを行わない。

## 現状 (確定した事実)

- 今回の実装対象は `src/physics/mixture_thermo_mod.F90`、`src/physics/nasa7_thermo_mod.F90` と対応する2本の単体試験。IEEE非有限値、逆数境界、演算overflowを fail-closed にし、通常範囲の計算順序と出力hashを維持しつつ、重複検証による性能低下を解消した。
- このhandoff追加前の凍結worktree hashは `8c4fff6220d7f1d1ad8e2d51852aa41eab4e5df2656d227e2db80ff05e027565`。証跡は `/tmp/pelef-0.245.0-post-finiteness-perf-matrix.eFeuG8/source-freeze.txt`。この文書は試験停止後に追加したため、以後の全ファイルhashには本ファイル分だけ差が加わるが、Fortran、CMake、テスト入力は変更していない。
- 最終8構成CTest matrixは4570/4570 PASS。内訳は Debug 524/524、Release 524/524、Debug+Cantera 529/529、MPI-Debug 694/694、MPI-Release 694/694、CVODE-Debug 533/533、CVODE-Debug+Cantera 539/539、CVODE-Release 533/533。証跡rootは `/tmp/pelef-0.245.0-post-finiteness-perf-matrix.eFeuG8`。MPI-Release完走は `mpi-release-clean-rerun-summary.txt:3-9`。
- fresh MPI構成にはGNU/OpenMPIだけを使用した。最初のMPI-Releaseはユーザー指示により227/694で停止したため算入せず、最初から再実行した694/694だけを最終matrixへ採用した。
- fresh tests-disabled Release install監査は42/42 ELFでPASS。build/install byte identity、42本のnon-executable stack、RPATH/RUNPATHなし、未解決依存なし、OpenMPI ABI 40のみ、oneAPI/`libmpi.so.12`なし。証跡は `/tmp/pelef-0.245.0-post-finiteness-audit.GifJDS/tests-disabled-summary.txt:1-30` と `audit.log:1-10`。
- 別のtests-enabled Release prefixではinstalled 0.245 transport-restart testを23件列挙し、23/23 PASS。証跡は同audit rootの `installed-transport-0245-ctest-summary.txt:1-10`。tests-on/off prefix間のbyte identityは主張しない。
- tests-disabled installを使ったcheckpoint schema 2/3/4/5のserial/MPI np1/2/4/8 qualificationは、application 32/32、checker 8/8、checkpoint/coarse/fine byte comparison 12/12 PASS。全hash、schema 4/5 transport diagnostics、np8 zero-fine-plane ownershipも一致した。証跡は同audit rootの `manual-installed-smoke-report.md`。
- schema 2の探索中に誤ってliteral `${n}` prefixを使った出力は `manual-installed-smoke/operator-diagnostic-schema2-literal-prefix` に保存し、上記32/8/12件から明示的に除外した。
- 現在、`/home/gtt13/projects/pelec` または上記evidence rootを参照するPeleF/CTest/buildプロセスは0件。外部GTT/OpenFOAMジョブには触れていない。2名の補助エージェントも完了済み。
- 初回の既存 `build/mpi-debug` は古いCMake cacheの `libmpi.so.12` 汚染でbuild停止した。ソース不具合ではない。この木は再利用せず、fresh build `/tmp/pelef-0.245.0-post-finiteness-perf-matrix.eFeuG8/mpi-debug-clean-build` で694/694 PASSした。
- README/docsは最終成功値へまだ更新していない。資格試験は完了し、文書反映と最終静的監査だけが0.245の残作業。
- Releaseの `amr_eb_patch_tree_reactive_2d_mod.F90` に出る allocatable descriptor `-Wmaybe-uninitialized` は、leaf return後の全使用がallocationに支配され、旧版にも同じ警告があり、直接単体試験も通るためGCC 13.3のfalse positiveと判定。0.245 blockerではなく、警告隠しの変更も不要。

## 次の一手

1. 実装を変更せず、`README.md`、`docs/validation/0.245.0.md`、`docs/implementation_status.md`、`docs/completion_roadmap.md`、`docs/parity_strategy.md`、`docs/pelec_mapping.md`、`docs/porting_plan.md` の「0.245 full matrix/clean install未実施」という古い記載を最終値へ更新する。
2. `docs/validation/evidence/0.245.0-clean-release-install-audit.txt` を今回のaudit rootへ更新し、`docs/validation/evidence/0.245.0-final-matrix.txt` を追加する。後者には4570/4570、各構成件数・時間、凍結hash、fresh OpenMPI build、停止した227件を除外した事実を記録する。
3. 文書編集後に `git diff --check`、Fortran 132桁制限、`python3 -m py_compile tools/*.py`、`tools/check_project_contract.py`、古い未実施記述のgrepを行う。編集対象diffを限定して確認し、実装凍結後の差が文書のみであることを示す。
4. 0.245資格文書を閉じた後もプロジェクト全体は未完成。roadmapから次の小さな実装単位を選ぶ。dynamic 3D AMR topology、physical boundaries、3D EB AMR、LES、Lagrangian particles/spray、scalable I/O、production physical validation、external PeleC field parityが残る。licenseはowner判断待ちのrelease blockerとして残す。

## 検証方法

- 0.245文書が4570/4570、42/42、23/23、32/32、8/8、12/12を正確に示し、focused smoke、full matrix、install監査、数値/物理PeleC parityを混同していないこと。
- tests-on/offの別prefixを明記し、prefix間のbyte identityを主張していないこと。
- 文書上の合格範囲は static、strictly interior、periodic、two-level transport persistenceのみ。dynamic AMR/topology、physical coarse boundaries、boundary-touching refinement、3D EB AMR、scalable/distributed I/O、crash-atomic replacement、payload authentication、performance/scaling、detailed-fuel physical validation、ignition、external PeleC field parity、license owner choiceは未完として残す。
- `git diff --check`、Python構文検査、project contractがすべてPASSし、実装凍結後に文書以外の差が増えていないこと。

## 触ってはいけないもの

- ユーザーの既存dirty/untracked変更。特に `git reset --hard`、`git checkout --`、一括削除をしない。commit、pushは2026-09-08のユーザー指示に基づくリモート退避に限る。
- 外部GTT計算とそのPID、出力、tmux。PeleF再開時も停止・変更しない。
- `/tmp/pelef-0.245.0-post-finiteness-perf-matrix.eFeuG8` と `/tmp/pelef-0.245.0-post-finiteness-audit.GifJDS` の確定証跡。再利用・上書きせずread-onlyで扱う。
- 中断前のMPI-Release 227件や `manual-installed-smoke/operator-diagnostic-schema2-literal-prefix` を最終資格結果へ混ぜない。
- licenseをプロジェクトownerの明示判断なしに選ばない。
