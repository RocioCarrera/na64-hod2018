# NA64 Hod2018 TTree Producer

**Author:** R. Carrera  
**Date:** May 2026  
**Experiment:** NA64 @ CERN (H4 beamline)  
**Data:** 2018 test beam, 100 GeV electron beam  
**Runs of interest:** 3610, 3655, 3661, 3662, 3663, 3664, 3665

---

## Overview

This repository contains code to reconstruct NA64 raw data from 2018 and produce
ROOT TTrees, with focus on the **Hod2018** hodoscope (16 X strips + 16 Y strips,
MPPC-based scintillator detector built by the Chilean team, P. Ulloa et al.).

The hodoscope data required special treatment: it uses **SADC preprocessing
mode 0-33** (full waveform storage) instead of the standard mode 5 (peak only).
The standard p348reco reconstruction does not decode it correctly for 2018 runs.
The solution is to read the raw `ChipSADC::Digit` waveform samples directly.

---

## Repository Structure

```
na64-hod2018/
├── README.md
├── code/
│   ├── reco2018tree.cc       # Main TTree producer
│   └── printdigits.cc        # Diagnostic: prints SADC detector names per event
├── patches/
│   ├── conddb_run_range.patch  # Fix: missing geometry for runs 3574-3854
│   └── led_fallback.patch      # Fix: missing LED calibrations for run 3661
├── maps/
│   └── HODO_2018_fixed.xml   # Fixed HODO.xml with extended run range (3610+)
└── plots/
    ├── hod2_strip_hits_X.png  # Beam profile: hits per X strip
    └── hod2_2d_hitmap.png     # 2D hit map (X strip vs Y strip)
```

---

## Dependencies

All dependencies are available on CERN lxplus via CVMFS.

- **p348-daq** (p348reco framework): `gitlab.cern.ch/P348/p348-daq`
- **GenFit**, **tracking-tools**, **catsc**: built from source (see setup below)
- **ROOT** 6.38.04 (via CVMFS on lxplus)

---

## Setup on lxplus

### 1. Clone and build p348-daq

```bash
git clone https://gitlab.cern.ch/P348/p348-daq.git
cd p348-daq
./build.sh
cd p348reco
```

### 2. Build tracking libraries

```bash
unset CC; unset CXX

# GenFit
cd GenFit && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../../p348reco/GenFit/install ..
make -j4 && make install && cd ../../p348reco

# tracking-tools
cd tracking-tools && mkdir -p build && cd build
cmake \
  -DCMAKE_INSTALL_PREFIX=../../p348reco/tracking-tools/install \
  -DGENFIT_INCLUDE_DIR=../../p348reco/GenFit/install/include \
  -DGENFIT_LIBRARY=../../p348reco/GenFit/install/lib/libgenfit2.so \
  -DCMAKE_C_COMPILER=`which gcc` \
  -DCMAKE_CXX_COMPILER=`which g++` ..
make -j4 && make install && cd ../../p348reco

# catsc (check catsc/install/ already exists, otherwise build similarly)
```

### 3. Set environment variables

Add to `~/.bashrc` (stored on AFS — applies to all lxplus nodes):

```bash
export ENABLE_TRACKING=1
export GENFIT_DIR=/afs/cern.ch/user/r/rocarrer/na64/p348-daq/p348reco/GenFit/install
export TTOOLS_DIR=/afs/cern.ch/user/r/rocarrer/na64/p348-daq/p348reco/tracking-tools/install
export CATSC_DIR=/afs/cern.ch/user/r/rocarrer/na64/p348-daq/p348reco/catsc/install
source ~/.bashrc
```

> **Note:** The variable is `ENABLE_TRACKING=1`, not `GENFIT=1` as stated in
> some older documentation. The Makefile checks for `ENABLE_TRACKING`.

---

## Bugs Fixed

Three bugs in the official p348-daq code prevented processing of 2018
calibration runs. Patches are in `patches/`.

### Bug 1 — Missing geometry for early 2018 runs (`conddb.cc`)

**File:** `p348reco/conddb.cc`, line ~6634  
**Problem:** Runs 3574–3854 had no tracking geometry assigned. The condition
called `Init_Geometry_2018_invis100()` only for runs 3855–4207, leaving a gap.
`LoadFromJson("")` was called with an empty path, crashing the JSON parser.

**Fix:**
```cpp
// Before:
if (3855 <= run && run <= 4207) Init_Geometry_2018_invis100();
// After:
if (3574 <= run && run <= 4207) Init_Geometry_2018_invis100();
```

### Bug 2 — Missing LED calibrations for run 3661 (`led.cc`)

**File:** `p348reco/led.cc`  
**Problem:** Run 3661 has no LED calibration data in `ledspill-2018-v3.txt`.
The standard fallback (previous spills + run average) also failed, causing abort.

**Fix:** After the standard fallback loop, search for the nearest run with data:

```cpp
#include <climits>  // add at top of file

// Replace throw with nearest-run fallback:
int best_run = -1, best_dist = INT_MAX;
for (auto& entry : LEDcalib) {
  if (entry.first.spill != 0) continue;
  int dist = abs((int)entry.first.run - (int)eid.run);
  if (dist < best_dist) { best_dist = dist; best_run = entry.first.run; }
}
if (best_run > 0) {
  loadid = {best_run, 0};
  hasData = (LEDcalib.find(loadid) != LEDcalib.end());
  cout << "WARNING: using LED data from nearest run " << best_run << endl;
}
if (!hasData) { cout << "WARNING: No LED corrections, skipping." << endl; return; }
```

### Bug 3 — Hod2018 not decoded: run range missing in maps (`HODO.xml`)

**File:** `maps/2018.xml/HODO.xml`  
**Problem:** The ChipSADC map for Hod2018 (Chilean team, port 3, source 623)
covered only runs 3667–4223. Runs 3610, 3655, 3661 were excluded — the decoder
silently ignored all HOD2 digits.

**Diagnosis:** `printdigits.cc` listed all SADC detector names per event.
HOD2X/HOD2Y were absent. The HODO.xml run ranges confirmed the gap.

**Fix:**
```xml
<!-- Before: -->
runs = "3667-4223"
<!-- After: -->
runs = "3610-4223"
```

---

## HOD2 Special Treatment: SADC Preprocessing Mode 0-33

In standard **mode 5**, the SADC hardware computes and stores only the peak
amplitude → `Cell.amplitude` is filled by `RunCellReco`.

In **mode 0-33**, the SADC stores the full waveform (all time samples).
`RunCellReco` is not called with a valid 2018 HOD2 calibration, so
`Cell.amplitude = 0` for all strips.

**Solution in `reco2018tree.cc`:** After `RunP348Reco()`, iterate raw
`ChipSADC::Digit` objects for detector names `"HOD2X"` and `"HOD2Y"`:

```cpp
const vector<CS::uint16>& samples = sadc->GetSamples();
// pedestal = mean of first 5 samples
// amplitude = max(samples) - pedestal
```

This gives the `HOD2_xamp[16]` and `HOD2_yamp[16]` branches — the primary
source of hodoscope data. The standard reco output (`HOD2_xamp_reco`,
`HOD2_xene_reco`, etc.) is also saved for cross-checking.

---

## Compiling and Running

```bash
cd /path/to/p348reco
cp reco2018tree.cc .
make reco2018tree.exe

# Single file
./reco2018tree.exe /eos/experiment/na64/data/raw/2018/cdr01001-003661.dat

# Find all files for a set of runs
find /eos/experiment/na64/data/raw/2018/ \
  -regextype posix-extended -regex ".*00366[1-5]\.dat"
```

Output: `reco2018tree.root` containing TTree `tree`.

---

## TTree Contents

### Event identification and trigger

| Branch | Type | Description |
|--------|------|-------------|
| `run`, `spill`, `spillevent`, `runevent` | `UInt_t` | Event ID |
| `isPhysics`, `isCalibration`, `isOnSpill` | `bool` | Event type |
| `isTriggerPhysics`, `isTriggerBeam`, `isTriggerRandom` | `bool` | Trigger type |
| `isTriggerECAL`, `isTriggerMuonInv` | `bool` | Trigger type |
| `masterTime` | `double` | Master timing reference (ns) |

### Calorimetry and beam counters

| Branch | Type | Description |
|--------|------|-------------|
| `S0_amp`, `S1_amp`, `S2_amp` | `double` | Scintillator amplitudes |
| `S0_t0`, `S1_t0`, `S2_t0` | `double` | Scintillator times (ns) |
| `ecalTotal`, `ecal0`, `ecal1` | `double` | ECAL total energy (GeV) |
| `ECAL_ene[2][6][6]` | `double[]` | ECAL cell energies |
| `ECAL_t0[2][6][6]` | `double[]` | ECAL cell times (ns) |
| `hcal0`, `hcal1`, `hcal2` | `double` | HCAL module energies (NaN for calib runs) |
| `hcal3` | `double` | Zero-degree HCAL energy (has data) |
| `srdTotal` | `double` | SRD total energy |
| `SRD_ene[4]`, `SRD_t0[4]` | `double[]` | SRD per-cell energy and time |
| `VETO_ene[6]`, `VETO_t0[6]` | `double[]` | VETO energies and times |
| `showerChi2` | `double` | ECAL shower profile chi2 |

### Hod2018 (HOD2) — primary detector of interest

| Branch | Type | Description |
|--------|------|-------------|
| `HOD2_xamp[16]` | `double[]` | **X strip signal amplitude (raw SADC)** |
| `HOD2_yamp[16]` | `double[]` | **Y strip signal amplitude (raw SADC)** |
| `HOD2_xped[16]`, `HOD2_yped[16]` | `double[]` | Pedestal per strip |
| `HOD2_xmax[16]`, `HOD2_ymax[16]` | `double[]` | Raw max sample value |
| `HOD2_xnsamples[16]`, `HOD2_ynsamples[16]` | `int[]` | Waveform length |
| `HOD2_xamp_reco[16]`, `HOD2_yamp_reco[16]` | `double[]` | Amplitude from RunCellReco |
| `HOD2_xene_reco[16]`, `HOD2_yene_reco[16]` | `double[]` | Energy from RunCellReco |
| `HOD2_xt0_reco[16]`, `HOD2_yt0_reco[16]` | `double[]` | Time from RunCellReco |
| `HOD2_xhit_reco[16]`, `HOD2_yhit_reco[16]` | `bool[]` | hasDigit flag |
| `hod2_xhits_reco`, `hod2_yhits_reco` | `vector<int>` | Reconstructed hit positions |

---

## Quick Analysis in ROOT

```cpp
TFile f("reco2018tree.root");
TTree* tree = (TTree*)f.Get("tree");

// Events with HOD2 signal
tree->GetEntries("HOD2_xamp[0]>0")

// Amplitude distribution for strip X0
tree->Draw("HOD2_xamp[0]", "HOD2_xamp[0]>5")

// Beam profile: hits per X strip
TH1F* hx = new TH1F("hx","Beam profile X;Strip;Events",16,0,16);
for(int i=0; i<16; i++)
  hx->SetBinContent(i+1, tree->GetEntries(Form("HOD2_xamp[%d]>50",i)));
hx->Draw();

// 2D hit map (strip X vs strip Y)
{
  TH2F* h2 = new TH2F("h2","HOD2 hit map;X strip;Y strip",16,0,16,16,0,16);
  double xamp[16], yamp[16];
  tree->SetBranchAddress("HOD2_xamp", xamp);
  tree->SetBranchAddress("HOD2_yamp", yamp);
  for(Long64_t ev=0; ev<tree->GetEntries(); ev++) {
    tree->GetEntry(ev);
    int xstrip=-1; double xmax=50;
    for(int i=0;i<16;i++) if(xamp[i]>xmax){xmax=xamp[i]; xstrip=i;}
    int ystrip=-1; double ymax=50;
    for(int i=0;i<16;i++) if(yamp[i]>ymax){ymax=yamp[i]; ystrip=i;}
    if(xstrip>=0 && ystrip>=0) h2->Fill(xstrip,ystrip);
  }
  h2->Draw("COLZ");
}
```

---

## Known Issues / Pending

- **Runs 3662, 3665:** Segfault in `DoHitAccumulation` (Micromega reco).
  Likely a different MM strip configuration in the maps. Under investigation.
- **HOD2 32 strips:** The report (Ulloa, 2019) shows a 32×32 hit map.
  Current mapping covers only 16 strips per plane (port 3, channels 32–63).
  Port 5 (channels 0–31) may also have been active — needs confirmation with P. Ulloa.
- **HCAL NaN:** `hcalTotalEnergy(0/1/2)` returns NaN for calibration runs
  (B field = 0T, HCAL not calibrated for this mode). `hcal3` works normally.

---

## Bugs to Report Upstream (P348 GitLab)

1. **`conddb.cc`:** Gap in geometry coverage for runs 3574–3854.
2. **`maps/2018.xml/HODO.xml`:** Run range for Hod2018 (port 3) starts at 3667
   but detector was active from at least run 3610.

---

## References

- Ulloa, P. (2019). *Hodoscopes Results Report*. NA64 internal note.
- NA64 DAQ: `gitlab.cern.ch/P348/p348-daq`
- 2018 beam test: H4 beamline, CERN SPS, May 2018
