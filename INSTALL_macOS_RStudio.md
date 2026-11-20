# Installing and Running treeX in RStudio on macOS

Step-by-step guide for the macOS-adapted version of **treeX**
(individual tree extraction from ground-based LiDAR point clouds).
Works on both Intel and Apple Silicon Macs.

---

## 1. One-time system setup (outside RStudio)

1. **Install R** (4.3 or newer) from <https://cran.r-project.org> — on
   Apple Silicon pick the **arm64** installer, on Intel the x86_64 one.
2. **Install RStudio Desktop** from <https://posit.co>.
3. **Apple command line tools** (needed to compile treeX's C++ code).
   In Terminal:

   ```sh
   xcode-select --install
   ```

   If it says "already installed", you're set.
4. **XQuartz** from <https://www.xquartz.org> — needed by the `rgl`
   package for 3D plots (log out and back in once after installing).
5. **CloudCompare** *(optional, but required for the `stemSplit()`
   connected-components step and the `openColoredLAZ()` viewer)*:
   download from <https://www.cloudcompare.org> or run
   `brew install --cask cloudcompare`, so it ends up at
   `/Applications/CloudCompare.app`.

---

## 2. Get the package source

Clone this repository (or unzip a downloaded copy) somewhere permanent,
e.g. so you end up with `~/treeX/` containing `DESCRIPTION`, `R/`, `src/`:

```sh
git clone https://github.com/antoniogrn/treeX.git ~/treeX
```

Alternatively, once published on GitHub you can skip steps 2 and 4 and
install directly with `pak::pak("antoniogrn/treeX")` (still run
step 3 first).

---

## 3. Install dependencies (in RStudio)

Open RStudio, then open the project via **File → Open Project…** →
`~/treeX/treeX.Rproj` (or set the working directory with
`setwd("~/treeX")`). Then run:

```r
source("_install_for_treeX.R")
```

This installs all CRAN dependencies (as macOS binaries, so no
compilation issues) and *tries* TreeLS from GitHub — if TreeLS fails to
compile, that is fine: treeX falls back to lidR automatically. Expect
this to take several minutes the first time.

---

## 4. Install treeX itself

Still from the `~/treeX` working directory:

```r
install.packages(".", repos = NULL, type = "source")
```

Equivalent alternatives: `pak::local_install(".")` or
`remotes::install_local(".", force = TRUE)`.

You will see clang compile `oregMclust.cc` and `init.c` in the output —
that is the native routine building for your Mac.

---

## 5. Verify the installation

```r
library(treeX)                 # should print "treeX version 1.0.27 loaded."
is.loaded("c_oregMcirc")       # should return TRUE (compiled circle-fitting routine)
```

If both work, the installation is complete.

---

## 6. Run it on a point cloud

Create a working folder, put your `.las`/`.laz` ground-based point cloud
in it, then follow `WORKFLOW_sample.R`:

```r
library(treeX)

dirPath <- path.expand("~/treeX_try/")        # all output lands here
if (!dir.exists(dirPath)) dir.create(dirPath)
setwd(dirPath)

LASfile    <- list.files(dirPath, pattern = "laz", full.names = TRUE)  # your input cloud
fileFinder <- "plot1"      # name used to tag all output files
voxelSize  <- 20           # cm; 20 is fast/rough, go down to 5 (or 2) for quality

extractVegetation(LASfile, fileFinder)        # 1. ground removal + DTM
clustSplit(fileFinder = fileFinder)           # 2. tree detection + DBH
crownFeel(fileFinder, voxelSize = voxelSize,  # 3. instance segmentation
          tileClipping = 0, retainPointClouds = TRUE)
computeTreeParams(fileFinder,                 # 4. height, crown metrics
                  voxelSize = voxelSize, crownParameters = TRUE)
```

Results (tree lists, measured parameters, segmented clouds) appear in
the folder structure treeX creates under `dirPath` — `trees_dbh.txt` is
the main tree list. During the parallel steps you should see
`using doParallel fork backend` — that is the macOS code path working.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| "C++ compiler not found" / clang errors in step 4 | Command line tools missing or license not accepted: rerun `xcode-select --install`, then `sudo xcodebuild -license accept`. |
| Linker fails with `tapi error: malformed file` / `unknown architecture` in `.tbd` files (often with `using SDK: 'NA'` earlier in the log) | The installed Command Line Tools are older than the macOS SDK (common after a macOS update). Reinstall them: `sudo rm -rf /Library/Developer/CommandLineTools` then `xcode-select --install`; verify with `xcrun --show-sdk-version`, restart RStudio and rerun the install. |
| A dependency tries to compile from source and fails (`sf`, `alphahull`, …) | Answer **No** when R asks "install from sources?", or run `options(pkgType = "binary")` before step 3; the CRAN binaries need no extra system libraries. |
| A parallel step crashes the RStudio session | Forked workers occasionally clash with the RStudio GUI; set the function's `nr_cores = 1` argument, or run the script with `Rscript` in Terminal instead. |
| `stemSplit()` stops with "CloudCompare executable not found" | Install CloudCompare (step 1.5); treeX looks in `/Applications` automatically, or pass the path via the `Path.CloudCompare` variable. |
| TreeLS fails to compile (`ld: library 'emutls_w' not found`, `/opt/gfortran` missing) | Expected and harmless — treeX uses its lidR fallback automatically. Only if you specifically want TreeLS: install the CRAN gfortran toolchain from <https://mac.r-project.org/tools/> and rerun the install script. |
| "package 'ADPclust' is not available for this version of R" | Harmless — ADPclust is archived from CRAN and treeX never calls it; the install script now tries GitHub (`ethanyxu/ADPclust`) and skips it cleanly if that fails. |

---

*See `PORTING_NOTES_macOS.md` inside the package for details on what was
changed to make treeX run on macOS.*
