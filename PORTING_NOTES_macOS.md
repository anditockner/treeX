# macOS port — what changed and why

treeX previously installed and ran only on Windows (and partly Linux). These
changes make it install and run on macOS (Intel and Apple Silicon) from
R/RStudio, without changing results on the other platforms.

## 1. Native code now compiles from source (`src/`)
The package relied on `useDynLib(edc)` with **precompiled binaries**
(`inst/libs/edc.so` for Linux x86-64, `edc.dll` for Windows). There was no
macOS binary, so the package could not load on a Mac at all.

* `src/oregMclust.cc` — the C++ implementation of `c_oregMcirc` (circular
  regression clustering), vendored from the CRAN package **edci 1.1-3**
  (GPL, Tim Garlipp), which is what the old `edc` binary was built from.
* `src/init.c` — routine registration under `R_init_treeX`.
* `NAMESPACE`: `useDynLib(edc)` → `useDynLib(treeX, .registration = TRUE)`
  (and the matching roxygen tag in `R/x1_treeLocalisationDensityClusters.R`).
* `inst/libs/` (prebuilt binaries) removed.

The code compiles cleanly with clang under C++17, so the standard Apple
command line tools (`xcode-select --install`) are all that is needed.

## 2. Installation no longer breaks on missing/Windows-only packages
Top-level `library()` calls in `R/treeExtraction_pls_at.R` ran at *install*
time and attached undeclared packages (`plyr`, `dplyr`, `sp`, `Morpho`,
`doMC`/`doParallel`). They were removed; all required packages are now
declared in `DESCRIPTION` (`plyr`, `dplyr`, `sp`, `Morpho`, `foreach` added
to Depends; optional ones moved to Suggests).

## 3. Parallel backend: doParallel everywhere
All seven `if (windows) doParallel else doMC` blocks now use
`doParallel::registerDoParallel(cores = n)` on Unix-likes, which forks
exactly like doMC did — `doMC` is no longer needed. The `mcprogress`
progress bar is now optional.

## 4. TreeLS is now optional
`TreeLS` is archived from CRAN and routinely fails to compile on recent
macOS toolchains. The seven `TreeLS::tlsSample(x, TreeLS::smp.voxelize(v))`
calls now go through the internal `voxelSample()` helper
(`R/zz_platform_compat.R`): it uses TreeLS when installed and otherwise
falls back to `lidR::decimate_points(random_per_voxel(res, 1))` — the same
one-random-point-per-voxel thinning.

The defensive package-loading block in `diameterBeast_i()` also no longer
hard-fails when archived packages (recexcavAAR, lmfor, TreeLS, …) are
missing; none of their functions are actually called.

## 5. CloudCompare integration is cross-platform
`stemSplit()` (connected components) and `openColoredLAZ()` hard-coded
`Cloudcompare.exe` under `C:/Program Files` and used the Windows-only
`system(..., show.output.on.console=)` argument. The internal helpers
`findCloudCompare()` and `runCommand()` now locate the binary per platform
(macOS: `/Applications/CloudCompare.app/Contents/MacOS/CloudCompare`) and
run it portably. A user-supplied path still wins.

## 6. Misc
* `pullFromFaroConnect()` overwrote its arguments with hard-coded `D:/…`
  paths (broken on every non-Windows machine) — the arguments are used now.
* `WORKFLOW_sample.R` / `WORKFLOW_treeX_cmd.R` use platform-neutral paths.
* `_install_for_treeX.R` rewritten: works on macOS, installs CRAN binaries,
  treats TreeLS as optional, skips doMC.

## Installing on macOS
```r
# once: xcode-select --install   (Terminal), XQuartz for rgl 3D plots
source("_install_for_treeX.R")            # dependencies
install.packages(".", repos = NULL, type = "source")   # from this folder
# or: pak::local_install(".")  /  remotes::install_local(".")
```
